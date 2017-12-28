// -*- LSST-C++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2018 LSST/AURA.
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
#ifndef LSST_AFW_MATH_POLYNOMIALS_Basis2d_h_INCLUDED
#define LSST_AFW_MATH_POLYNOMIALS_Basis2d_h_INCLUDED
#ifdef DOXYGEN

namespace lsst { namespace afw { namespace math { namespace polynomials {

/**
 *  A basis concept for 2-d series expansions.
 *
 *  @note This class is only present in the documentation, as it represents an
 *        abstract concept.
 */
template <typename Basis1d>
class Basis2d {
public:

    /// A Function2d object that uses this basis.
    using Function = ...;

    /// The type returned by scale().
    using Scaled = ...;

    /// The type returned by makeWorkspace().
    using Workspace = ...;

    /// The type returned by getIndices().
    using IndexRange = ...;

    /// Return the maximum order of the basis.
    std::size_t getOrder() const;

    /// Return the number of basis functions.
    std::size_t size() const;

    /**
     *  Return a scaled basis that delegates to a copy of `this`.
     *
     *  The scaled basis will transform all points by the given scaling
     *  before evaluating the basis functions in the same way as `this`.
     */
    Scaled scale(Scaling2d const & first) const;

    /// Return the flattened index of the basis function with the given x and y orders.
    std::size_t index(std::size_t x, std::size_t y) const;

    /**
     *  Return a range of iterators that dereference to Index2d.
     *
     *  This provides the most natural way to interpret the packed coefficients and basis functions
     *  utilized by PackedBasis2d; for example,
     *  ```
     *      // Evaluate basis functions at a point.
     *      geom::Point2D point(x, y);
     *      PolynomialBasis2d basis(order);
     *      std::vector<double> values(basis.size());
     *      basis.fill(point, values);
     *      // Iterate over tuples of flattened indices and x and y orders.
     *      for (auto const & index : basis.getIndices()) {
     *          double a = values[index.flat];
     *          // standard polynomial basis functions are just powers
     *          double b = std::pow(point.getX(), index.x)*std::pow(point.getY(), index.y)];
     *          assert(std::pow(a - b) < std::numeric_limits<double>::epsilon());
     *      }
     *  ```
     */
    IndexRange getIndices() const;

    /// Allocate workspace that can be passed to sumWith() and fill() to avoid repeated memory allocations.
    Workspace makeWorkspace() const;

    /**
     *  Evaluate a basis expansion with the given coefficients.
     *
     *  If the 1-d basis elements are @f$B_n(x)@f$ and the given coefficients are
     *  a vector @f$a_{p, q}@f$, this computes
     *  @f[
     *      \sum_{p = 0, q = 0}^{p + q \le N} a_{p,q} B_{p}(x) B_{q}(y)
     *  @f]
     *
     *  @param[in] point         Point at which to evaluate the expansion.
     *  @param[in] coefficients  Flattened coefficients vector.  May be any
     *                           type for which `coefficients[n]` returns an
     *                           object convertible to `double` for all `n <=
     *                           getOrder()`.  This includes
     *                           `std::vector<double>`,
     *                           `ndarray::Array<double,1>`,
     *                           `Eigen::VectorXd`, and random access
     *                           iterators.  If a lazy expression template
     *                           object is passed, the elements of the
     *                           expression will be evaluated only once.
     */
    template <typename Vector>
    double sumWith(geom::Point2D const & point, Vector const & coefficients) const;

    /// Evaluate a basis expansion with the given coefficients (external workspace version).
    template <typename Vector>
    double sumWith(geom::Point2D const & point, Vector const & coefficients, Workspace & workspace) const;

    /**
     *  Evaluate the basis at a given point.
     *
     *  @param[in] point  Point at which to evaluate the basis functions.
     *  @param[in] basis  Flattened output vector.  May be any type for which
     *                    `coefficients[n]` returns a non-const reference to a
     *                    floating-point value.  This includes
     *                    `std::vector<double>`, `ndarray::Array<double,1>`,
     *                    `Eigen::VectorXd`, `Eigen` view expressions, and
     *                    mutable random access iterators.
     */
    template <typename Vector>
    void fill(geom::Point2D const & point, Vector && basis) const;

    /// Evaluate the basis at a given point (external workspace version).
    template <typename Vector>
    void fill(geom::Point2D const & point, Vector && basis, Workspace & workspace) const;

private:
    Basis1d _basis1d;
};

}}}} // namespace lsst::afw::math::polynomials

#endif // DOXYGEN
#endif // !LSST_AFW_MATH_POLYNOMIALS_Basis2d_h_INCLUDED
