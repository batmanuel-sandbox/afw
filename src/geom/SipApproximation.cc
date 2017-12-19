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

namespace lsst { namespace afw { namespace geom {

namespace {

//
// 2-D POLYNOMIAL PACKING
//
// 2-D polynomials in this file are packed into 1-d arrays with the following ordering:
// [(0,0), (0,1), (1,0), (0,2), (1,1), (2,0), ...]
//

// Index into the first element in a packed 2-d polynomial array with combined order n.
int index1(int n) {
    return n*(n+1)/2;
}

// Index into a packged 2-d polynomial array with powers p and q.
int index2(int p, int q) {
    return index1(p + q) + p;
}

// Class that evaluates 2-d polynomials efficiently.
class PolynomialEvaluator {
public:

    // Construct an evaluator that evaluates polynomials of the given combined order.
    // User must then call at() before any other methods.
    explicit PolynomialEvaluator(int order) :
        _order(order),
        _xPow(_order + 1),
        _yPow(_order + 1)
    {
        _xPow[0] = 1.0;
        _yPow[0] = 1.0;
    }

    // Set the evaluator to work at the given point in x and y.
    void at(Point2D const & point) {
        for (int n = 1; n <= _order; ++n) {
#ifdef USE_POW_FOR_POLYS
            _xPow[n] = std::pow(point.getX(), n);
            _yPow[n] = std::pow(point.getY(), n);
#else
            _xPow[n] = _xPow[n - 1]*point.getX();
            _yPow[n] = _yPow[n - 1]*point.getY();
#endif
        }
    }

    // Evaluate the polynomial with the given packed coefficient vector.
    template <typename Vector>
    double sumWith(Vector const & coefficients) const {
        double r = 0.0;
        for (int n = 0, j = 0; n <= _order; ++n) {
            for (int p = 0, q = n; p <= n; ++p, --q, ++j) {
                assert(index2(p, q) == j);
                r += coefficients[j]*_xPow[p]*_yPow[q];
            }
        }
        return r;
    }

    // Fill a vector such that it evaluates the polynomial when its inner
    // product with a packed coefficient vector is taken.
    template <typename Vector>
    void fillBasis(Vector && basis) const {
        // Need to use universal references because this might be a non-const
        // reference to a Matrix object or a temporary Block object passed as an
        // rvalue (that we can nevertheless write to because it's a view).
        for (int n = 0, j = 0; n <= _order; ++n) {
            for (int p = 0, q = n; p <= n; ++p, --q, ++j) {
                assert(index2(p, q) == j);
                std::forward<Vector>(basis)[j] = _xPow[p]*_yPow[q];
            }
        }
    }

private:
    int const _order;
    Eigen::VectorXd _xPow;
    Eigen::VectorXd _yPow;
};

// Construct a matrix that evaluates a polynomial at a list of points when multiplied by a packed
// coefficient vector.
Eigen::MatrixXd makePolynomialBasis(int const order, std::vector<Point2D> const & points) {
    Eigen::MatrixXd result = Eigen::MatrixXd::Zero(points.size(), index1(order + 1));
    PolynomialEvaluator poly(order);
    for (int i = 0; i < result.rows(); ++i) {
        poly.at(points[i]);
        poly.fillBasis(result.row(i));
    }
    return result;
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

template <typename Vector>
void subtractPointsIntoVectors(
    Vector && x,
    Vector && y,
    std::vector<Point2D> const & a,
    std::vector<Point2D> const & b
) {
    for (std::size_t i = 0; i < a.size(); ++i) {
        std::forward<Vector>(x)[i] = a[i].getX() - b[i].getX();
        std::forward<Vector>(y)[i] = a[i].getY() - b[i].getY();
    }
}

} // anonymous

// Private implementation object for SipApproximation that manages the grid of points on which
// we evaluate the exact transform.
struct SipApproximation::Grid {

    // Set up the grid.
    Grid(Extent2I const & dims, SipApproximation const & parent);

    Extent2I const dimensions; //  number of grid points in each dimension
    std::vector<Point2D> dpix; //  [pixel coords] - CRPIX
    std::vector<Point2D> siwc; //  CD^{-1}([intermediate world coords])
    std::vector<Point2D> dpixI; //  round-tripped version of dpix if useInverse, or exactly dpix
};

// Private implementation object for SipApproximation that manages the solution
struct SipApproximation::Solution {

    // Solve for the best-fit coefficients.
    Solution(int order_, double svdThreshold, SipApproximation const & parent);

    // Read in external coefficients.
    Solution(ndarray::Array<double const, 2> const & a,
             ndarray::Array<double const, 2> const & b,
             ndarray::Array<double const, 2> const & ap,
             ndarray::Array<double const, 2> const & bp);

    int const order;
    Eigen::VectorXd a;   // SIP coefficients, packed
    Eigen::VectorXd b;
    Eigen::VectorXd ap;
    Eigen::VectorXd bp;
};

SipApproximation::Grid::Grid(Extent2I const & dims, SipApproximation const & parent) :
    dimensions(dims),
    dpix(makeGrid(parent._bbox, dimensions)),
    siwc(parent._pixelToIwc->applyForward(dpix)),
    dpixI(parent._useInverse ? parent._pixelToIwc->applyInverse(siwc) : dpix)
{
    // Apply the CRPIX offset to dpix
    std::for_each(dpix.begin(), dpix.end(), [&parent](Point2D & p){ p -= parent._crpix; });
    // Apply the CRPIX offset to dpixI
    std::for_each(dpixI.begin(), dpixI.end(), [&parent](Point2D & p){ p -= parent._crpix; });
    // Apply the CD^{-1} transform to siwc
    std::for_each(siwc.begin(), siwc.end(), [&parent](Point2D & p){ p = parent._cdInv(p); });
}

SipApproximation::Solution::Solution(int order_, double svdThreshold, SipApproximation const & parent) :
    order(order_)
{
    if (index1(order + 1) > static_cast<int>(parent._grid->dpix.size())) {
        throw LSST_EXCEPT(
            pex::exceptions::LogicError,
            (boost::format("Number of parameters (%d) is larger than number of data points (%d)")
             % (2*index1(order + 1)) % (2*parent._grid->dpix.size())).str()
        );
    }

    Eigen::VectorXd xRhs(parent._grid->dpix.size());
    Eigen::VectorXd yRhs(parent._grid->dpix.size());

    // xRhs = siwc.x - dpix.x
    // yRhs = siwc.y - dpix.y
    subtractPointsIntoVectors(xRhs, yRhs, parent._grid->siwc, parent._grid->dpix);

    // Compute the basis matrices for the forward and inverse problems.
    auto fwdMatrix = makePolynomialBasis(order, parent._grid->dpix);

    // Since we're not trying to null the zeroth- and first-order terms, the
    // solution is just linear least squares, and we can do that with SVD.
    auto fwdDecomp = fwdMatrix.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);
    if (svdThreshold >= 0) {
        fwdDecomp.setThreshold(svdThreshold);
    }
    a = fwdDecomp.solve(xRhs);
    b = fwdDecomp.solve(yRhs);

    // xRhs = dpixI.x - siwc.x
    // yRhs = dpixI.y - siwc.y
    subtractPointsIntoVectors(xRhs, yRhs, parent._grid->dpixI, parent._grid->siwc);

    auto invMatrix = makePolynomialBasis(order, parent._grid->siwc);

    auto invDecomp = invMatrix.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);
    if (svdThreshold >= 0) {
        invDecomp.setThreshold(svdThreshold);
    }

    ap = invDecomp.solve(xRhs);
    bp = invDecomp.solve(yRhs);
}

SipApproximation::Solution::Solution(
    ndarray::Array<double const, 2> const & a_,
    ndarray::Array<double const, 2> const & b_,
    ndarray::Array<double const, 2> const & ap_,
    ndarray::Array<double const, 2> const & bp_
) :
    order(a_.getSize<0>() - 1),
    a(index1(order + 1)),
    b(a.size()),
    ap(a.size()),
    bp(a.size())
{
    LSST_THROW_IF_NE(a_.getSize<0>(), a_.getSize<1>(), pex::exceptions::InvalidParameterError,
                     "A matrix must be square (%d != %d).");
    LSST_THROW_IF_NE(b_.getSize<0>(), b_.getSize<1>(), pex::exceptions::InvalidParameterError,
                     "B matrix must be square (%d != %d).");
    LSST_THROW_IF_NE(ap_.getSize<0>(), ap_.getSize<1>(), pex::exceptions::InvalidParameterError,
                     "AP matrix must be square (%d != %d).");
    LSST_THROW_IF_NE(bp_.getSize<0>(), bp_.getSize<1>(), pex::exceptions::InvalidParameterError,
                     "BP matrix must be square (%d != %d).");
    LSST_THROW_IF_NE(a_.getSize<0>(), b_.getSize<0>(), pex::exceptions::InvalidParameterError,
                     "A and B matrices must have the same size (%d != %d).");
    LSST_THROW_IF_NE(a_.getSize<0>(), ap_.getSize<0>(), pex::exceptions::InvalidParameterError,
                     "A and AP matrices must have the same size (%d != %d).");
    LSST_THROW_IF_NE(a_.getSize<0>(), bp_.getSize<0>(), pex::exceptions::InvalidParameterError,
                     "A and BP matrices must have the same size (%d != %d).");
    for (int n = 0, j = 0; n <= order; ++n) {
        for (int p = 0, q = n; p <= n; ++p, --q, ++j) {
            assert(index2(p, q) == j);
            a[j] = a_[p][q];
            b[j] = b_[p][q];
            ap[j] = ap_[p][q];
            bp[j] = bp_[p][q];
        }
    }
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
    _cdInv(cd.invert())
{
    _grid = std::make_unique<Grid>(gridDimensions, *this);
    _solution = std::make_unique<Solution>(order, svdThreshold, *this);
}

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
    _cdInv(cd.invert())
{
    _grid = std::make_unique<Grid>(gridDimensions, *this);
    _solution = std::make_unique<Solution>(a, b, ap, bp);
}

SipApproximation::~SipApproximation() {}

int SipApproximation::getOrder() const {
    return _solution->order;
}

double SipApproximation::getA(int p, int q) const {
    return _solution->a[index2(p, q)];
}

double SipApproximation::getB(int p, int q) const {
    return _solution->b[index2(p, q)];
}

double SipApproximation::getAP(int p, int q) const {
    return _solution->ap[index2(p, q)];
}

double SipApproximation::getBP(int p, int q) const {
    return _solution->bp[index2(p, q)];
}


Point2D SipApproximation::applyForward(Point2D const & pix) const {
    std::vector<Point2D> vec(1, pix);
    return applyForward(vec).front();
}

std::vector<Point2D> SipApproximation::applyForward(std::vector<Point2D> const & pix) const {
    PolynomialEvaluator poly(_solution->order);
    std::vector<Point2D> iwc;
    iwc.reserve(pix.size());
    auto cd = _cdInv.invert();
    for (auto const & point : pix) {
        auto dpix = point - _crpix;
        poly.at(dpix);
        iwc.push_back(cd(dpix + Extent2D(poly.sumWith(_solution->a), poly.sumWith(_solution->b))));
    }
    return iwc;
}

Point2D SipApproximation::applyInverse(Point2D const & iwc) const {
    std::vector<Point2D> vec(1, iwc);
    return applyInverse(vec).front();
}

std::vector<Point2D> SipApproximation::applyInverse(std::vector<Point2D> const & iwc) const {
    PolynomialEvaluator poly(_solution->order);
    std::vector<Point2D> pix;
    pix.reserve(iwc.size());
    for (auto const & point : iwc) {
        auto siwc = _cdInv(point);
        poly.at(siwc);
        pix.push_back(siwc + Extent2D(poly.sumWith(_solution->ap), poly.sumWith(_solution->bp)) + _crpix);
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
    _solution = std::make_unique<Solution>(order, svdThreshold, *this);
}

namespace {

template <typename It1, typename It2, typename Func>
void forEachZip(It1 it1, It1 last1, It2 it2, Func func) {
    for (; it1 != last1; ++it1, ++it2) {
        func(*it1, *it2);
    }
}

} // anonymous

std::pair<double, double> SipApproximation::computeMaxDeviation() const {
    std::pair<double, double> maxDiff(0.0, 0.0);
    PolynomialEvaluator poly(_solution->order);
    for (std::size_t i = 0; i < _grid->dpix.size(); ++i) {
        poly.at(_grid->dpix[i]);
        auto siwc2 = _grid->dpix[i] +
            Extent2D(poly.sumWith(this->_solution->a), poly.sumWith(this->_solution->b));
        poly.at(_grid->siwc[i]);
        auto dpix2 = _grid->siwc[i] +
            Extent2D(poly.sumWith(this->_solution->ap), poly.sumWith(this->_solution->bp));
        maxDiff.first = std::max(maxDiff.first, (_grid->siwc[i] - siwc2).computeNorm());
        maxDiff.second = std::max(maxDiff.second, (_grid->dpixI[i] - dpix2).computeNorm());
    }
    return maxDiff;
}


}}}  // namespace lsst::afw::geom
