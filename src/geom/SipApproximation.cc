// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2017 LSST Corporation.
 *
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the LSST License Statement and
 * the GNU General Public License along with this program.  If not,
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */

#include <algorithm>
#include <vector>

#include "Eigen/SVD"
#include "Eigen/QR"
#include "lsst/afw/geom/SipApproximation.h"
#include "lsst/afw/math/polynomials/PolynomialFunction2d.h"

namespace lsst { namespace afw { namespace geom {

namespace poly = math::polynomials;

namespace {

std::pair<poly::PolynomialFunction2d, poly::PolynomialFunction2d> fitSipOneDirection(
    int order,
    Box2D const & box,
    double svdThreshold,
    std::vector<Point2D> const & input,
    std::vector<Point2D> const & output
) {
    // The scaled polynomial basis evaluates polynomials after mapping the
    // input coordinates from the given box to [-1, 1]x[-1, 1] (for numerical
    // stability).
    auto basis = poly::makeScaledPolynomialBasis2d(order, box);
    auto workspace = basis.makeWorkspace();
    Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(input.size(), basis.size());
    Eigen::VectorXd xRhs(input.size());
    Eigen::VectorXd yRhs(input.size());
    for (int i = 0; i < matrix.rows(); ++i) {
        basis.fill(input[i], matrix.row(i), workspace);
        auto rhs = output[i] - input[i];
        xRhs[i] = rhs.getX();
        yRhs[i] = rhs.getY();
    }
    // Since we're not trying to null the zeroth- and first-order terms, the
    // solution is just linear least squares, and we can do that with SVD.
    auto decomp = matrix.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);
    if (svdThreshold >= 0) {
        decomp.setThreshold(svdThreshold);
    }
    auto scaledX = makeFunction2d(basis, decomp.solve(xRhs));
    auto scaledY = makeFunction2d(basis, decomp.solve(yRhs));
    // On return, we simplify the polynomials by moving the remapping transform
    // into the coefficients themselves.
    return std::make_pair(simplify(scaledX), simplify(scaledY));
}

// Return a vector of points on a grid, covering the given bounding box.
std::vector<Point2D> makeGrid(Box2D const & bbox, Extent2I const & dimensions) {
    std::vector<Point2D> points;
    points.reserve(dimensions.getX()*dimensions.getY());
    double const dx = bbox.getWidth()/dimensions.getX();
    double const dy = bbox.getHeight()/dimensions.getY();
    for (int iy = 0; iy < dimensions.getY(); ++iy) {
        double const y = bbox.getMinY() + iy*dy;
        for (int ix = 0; ix < dimensions.getX(); ++ix) {
            points.emplace_back(bbox.getMinX() + ix*dx, y);
        }
    }
    return points;
}

poly::PolynomialFunction2d readFunctionFromCoeffMatrix(ndarray::Array<double const, 2> const & coeffs) {
    LSST_THROW_IF_NE(coeffs.getSize<0>(), coeffs.getSize<1>(), pex::exceptions::InvalidParameterError,
                     "Coefficient matrix must be square (%d != %d).");
    poly::PolynomialBasis2d basis(coeffs.getSize<0>() - 1);
    Eigen::VectorXd packed(basis.size());
    for (auto const & i : basis.getIndices()) {
        packed[i.flat] = coeffs[i.x][i.y];
    }
    return poly::PolynomialFunction2d(basis, packed);
}

} // anonymous

// Private implementation object for SipApproximation that manages the grid of points on which
// we evaluate the exact transform.
struct SipApproximation::Grid {

    // Set up the grid.
    Grid(Extent2I const & dims, SipApproximation const & parent);

    Extent2I const dimensions;  //  number of grid points in each dimension
    std::vector<Point2D> dpix1; //  [pixel coords] - CRPIX
    std::vector<Point2D> siwc;  //  CD^{-1}([intermediate world coords])
    std::vector<Point2D> dpix2; //  round-tripped version of dpix1 if useInverse, or exactly dpix1
};

// Private implementation object for SipApproximation that manages the solution
struct SipApproximation::Solution {

    static std::unique_ptr<Solution> fit(int order_, double svdThreshold, SipApproximation const & parent);

    Solution(poly::PolynomialFunction2d const & a_,
             poly::PolynomialFunction2d const & b_,
             poly::PolynomialFunction2d const & ap_,
             poly::PolynomialFunction2d const & bp_) :
        a(a_), b(b_), ap(ap_), bp(bp_)
    {
        LSST_THROW_IF_NE(a.getBasis().getOrder(), b.getBasis().getOrder(),
                         pex::exceptions::InvalidParameterError,
                         "A and B polynomials must have the same order (%d != %d).");
        LSST_THROW_IF_NE(a.getBasis().getOrder(), ap.getBasis().getOrder(),
                         pex::exceptions::InvalidParameterError,
                         "A and AP polynomials must have the same order (%d != %d).");
        LSST_THROW_IF_NE(a.getBasis().getOrder(), bp.getBasis().getOrder(),
                         pex::exceptions::InvalidParameterError,
                         "A and BP polynomials must have the same order (%d != %d).");
    }

    using Workspace = poly::PolynomialFunction2d::Workspace;

    Workspace makeWorkspace() const { return a.makeWorkspace(); }

    Point2D applyForward(Point2D const & dpix, Workspace & ws) const {
        return dpix + Extent2D(a(dpix, ws), b(dpix, ws));
    }

    Point2D applyInverse(Point2D const & siwc, Workspace & ws) const {
        return siwc + Extent2D(ap(siwc, ws), bp(siwc, ws));
    }

    poly::PolynomialFunction2d a;
    poly::PolynomialFunction2d b;
    poly::PolynomialFunction2d ap;
    poly::PolynomialFunction2d bp;
};

SipApproximation::Grid::Grid(Extent2I const & dims, SipApproximation const & parent) :
    dimensions(dims),
    dpix1(makeGrid(parent._bbox, dimensions)),
    siwc(parent._pixelToIwc->applyForward(dpix1))
{
    // Apply the CRPIX offset to make pix1 into dpix1 (in-place)
    std::for_each(dpix1.begin(), dpix1.end(), [&parent](Point2D & p){ p -= parent._crpix; });

    if (parent._useInverse) {
        // Set from the given inverse of the given pixels-to-iwc transform
        // Note that at this point, siwc is still just iwc, because the scaling by cdInv is later.
        dpix2 = parent._pixelToIwc->applyInverse(siwc);
        // Apply the CRPIX offset to make pix1 into dpix2 (in-place)
        std::for_each(dpix2.begin(), dpix2.end(), [&parent](Point2D & p){ p -= parent._crpix; });
    } else {
        // Just make dpix2 = dpix1, and hence fit to the true inverse of pixels-to-iwc.
        dpix2 = dpix1;
    }

    // Apply the CD^{-1} transform to siwc
    std::for_each(siwc.begin(), siwc.end(), [&parent](Point2D & p){ p = parent._cdInv(p); });
}

std::unique_ptr<SipApproximation::Solution> SipApproximation::Solution::fit(
    int order,
    double svdThreshold,
    SipApproximation const & parent
) {
    poly::PolynomialBasis2d basis(order);
    if (basis.size() > parent._grid->dpix1.size()) {
        throw LSST_EXCEPT(
            pex::exceptions::LogicError,
            (boost::format("Number of parameters (%d) is larger than number of data points (%d)")
             % (2*basis.size()) % (2*parent._grid->dpix1.size())).str()
        );
    }

    Box2D boxFwd(parent._bbox);
    boxFwd.shift(-parent._crpix);
    auto fwd = fitSipOneDirection(order, boxFwd, svdThreshold, parent._grid->dpix1, parent._grid->siwc);

    Box2D boxInv;
    for (auto const & point : parent._grid->siwc) {
        boxInv.include(point);
    }
    auto inv = fitSipOneDirection(order, boxInv, svdThreshold, parent._grid->siwc, parent._grid->dpix2);

    return std::make_unique<Solution>(fwd.first, fwd.second, inv.first, inv.second);
}

SipApproximation::SipApproximation(
    std::shared_ptr<TransformPoint2ToPoint2> pixelToIwc,
    Point2D const & crpix,
    LinearTransform const & cd,
    Box2D const & bbox,
    Extent2I const & gridDimensions,
    int order,
    bool useInverse,
    double svdThreshold
) :
    _useInverse(useInverse),
    _pixelToIwc(std::move(pixelToIwc)),
    _bbox(bbox),
    _crpix(crpix),
    _cdInv(cd.invert()),
    _grid(new Grid(gridDimensions, *this)),
    _solution(Solution::fit(order, svdThreshold, *this))
{}

SipApproximation::SipApproximation(
    std::shared_ptr<TransformPoint2ToPoint2> pixelToIwc,
    Point2D const & crpix,
    LinearTransform const & cd,
    Box2D const & bbox,
    Extent2I const & gridDimensions,
    ndarray::Array<double const, 2> const & a,
    ndarray::Array<double const, 2> const & b,
    ndarray::Array<double const, 2> const & ap,
    ndarray::Array<double const, 2> const & bp,
    bool useInverse
) :
    _useInverse(useInverse),
    _pixelToIwc(std::move(pixelToIwc)),
    _bbox(bbox),
    _crpix(crpix),
    _cdInv(cd.invert()),
    _grid(new Grid(gridDimensions, *this)),
    _solution(
        new Solution(
            readFunctionFromCoeffMatrix(a),
            readFunctionFromCoeffMatrix(b),
            readFunctionFromCoeffMatrix(ap),
            readFunctionFromCoeffMatrix(bp)
        )
    )
{}

SipApproximation::~SipApproximation() {}

int SipApproximation::getOrder() const {
    return _solution->a.getBasis().getOrder();
}

double SipApproximation::getA(int p, int q) const {
    return _solution->a.getCoefficient(_solution->a.getBasis().index(p, q));
}

double SipApproximation::getB(int p, int q) const {
    return _solution->b.getCoefficient(_solution->a.getBasis().index(p, q));
}

double SipApproximation::getAP(int p, int q) const {
    return _solution->ap.getCoefficient(_solution->a.getBasis().index(p, q));
}

double SipApproximation::getBP(int p, int q) const {
    return _solution->bp.getCoefficient(_solution->a.getBasis().index(p, q));
}

Point2D SipApproximation::applyForward(Point2D const & pix) const {
    auto cd = _cdInv.invert();
    auto ws = _solution->makeWorkspace();
    return cd(_solution->applyForward(pix - _crpix, ws));
}

std::vector<Point2D> SipApproximation::applyForward(std::vector<Point2D> const & pix) const {
    auto ws = _solution->makeWorkspace();
    std::vector<Point2D> iwc;
    iwc.reserve(pix.size());
    auto cd = _cdInv.invert();
    for (auto const & point : pix) {
        iwc.push_back(cd(_solution->applyForward(point - _crpix, ws)));
    }
    return iwc;
}

Point2D SipApproximation::applyInverse(Point2D const & iwc) const {
    auto ws = _solution->makeWorkspace();
    return _solution->applyInverse(_cdInv(iwc), ws) + _crpix;
}

std::vector<Point2D> SipApproximation::applyInverse(std::vector<Point2D> const & iwc) const {
    auto ws = _solution->makeWorkspace();
    std::vector<Point2D> pix;
    pix.reserve(iwc.size());
    for (auto const & point : iwc) {
        pix.push_back(_solution->applyInverse(_cdInv(point), ws) + _crpix);
    }
    return pix;
}

Extent2D SipApproximation::getGridStep() const {
    return Extent2D(_bbox.getWidth()/_grid->dimensions.getX(),
                          _bbox.getHeight()/_grid->dimensions.getY());
}

Extent2I SipApproximation::getGridDimensions() const {
    return _grid->dimensions;
}

void SipApproximation::updateGrid(Extent2I const & dimensions) {
    _grid = std::make_unique<Grid>(dimensions, *this);
}

void SipApproximation::refineGrid(int f) {
    // We shrink the grid spacing by the given factor, which is not the same
    // as increasing the number of grid points by that factor, because there
    // is one more grid point that step in each dimension.
    Extent2I unit(1);
    updateGrid((_grid->dimensions - unit)*f + unit);
}

void SipApproximation::fit(int order, double svdThreshold) {
    _solution = Solution::fit(order, svdThreshold, *this);
}

std::pair<double, double> SipApproximation::computeMaxDeviation() const {
    std::pair<double, double> maxDiff(0.0, 0.0);
    auto ws = _solution->makeWorkspace();
    for (std::size_t i = 0; i < _grid->dpix1.size(); ++i) {
        auto siwc2 = _solution->applyForward(_grid->dpix1[i], ws);
        auto dpix2 = _solution->applyInverse(_grid->siwc[i], ws);
        maxDiff.first = std::max(maxDiff.first, (_grid->siwc[i] - siwc2).computeNorm());
        maxDiff.second = std::max(maxDiff.second, (_grid->dpix2[i] - dpix2).computeNorm());
    }
    return maxDiff;
}

}}}  // namespace lsst::afw::geom
