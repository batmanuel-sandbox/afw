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
#ifndef LSST_AFW_MATH_POLYNOMIALS_Basis1d_h_INCLUDED
#define LSST_AFW_MATH_POLYNOMIALS_Basis1d_h_INCLUDED
#ifdef DOXYGEN

namespace lsst { namespace afw { namespace math { namespace polynomials {

/**
 *  A basis concept for 1-d series expansions.
 *
 *  @note This class is only present in the documentation, as it represents an
 *        abstract concept.
 */
class Basis1d {
public:

    /// A Function1d object that uses this basis.
    using Function = ...;

    /// The type returned by scale().
    using Scaled = ...;

    /// Return the order of the basis.
    std::size_t getOrder() const;

    /// Return the number of elements in the basis.
    std::size_t size() const;

    /**
     *  Return a scaled basis that delegates to a copy of `this`.
     *
     *  The scaled basis will transform all points by the given scaling
     *  before evaluating the basis functions in the same way as `this`.
     */
    Scaled scale(Scaling1d const & scaling) const;

    /**
     *  Evaluate a basis expansion with the given coefficients.
     *
     *  If the basis elements are @f$B_n(x)@f$ and the given coefficients are
     *  a vector @f$a_n@f$, this computes
     *  @f[
     *      \sum_{n = 0}^{n \le N} a_n B_n(x)
     *  @f]
     *
     *  @param[in] x             Point at which to evaluate the expansion.
     *  @param[in] coefficients  Coefficients vector.  May be any type for
     *                           which `coefficients[n]` returns an object
     *                           convertible to `double` for all `n <=
     *                           getOrder()`.  This includes
     *                           `std::vector<double>`,
     *                           `ndarray::Array<double,1>`,
     *                           `Eigen::VectorXd`, and random access
     *                           iterators.  If a lazy expression template
     *                           object is passed, the elements of the
     *                           expression will be evaluated only once.
     *
     *  @exceptsafe Does not throw unless `coefficients[n]` does, and provides
     *              the same exception safety as it if it does.
     */
    template <typename Vector>
    double sumWith(double x, Vector const & coefficients) const;

    /**
     *  Evaluate the basis at a given point.
     *
     *  @param[in] x      Point at which to evaluate the basis functions.
     *  @param[in] basis  Output vector.  May be any type for which
     *                    `coefficients[n]` returns a non-const reference to a
     *                    floating-point value.  This includes
     *                    `std::vector<double>`, `ndarray::Array<double,1>`,
     *                    `Eigen::VectorXd`, `Eigen` view expressions, and
     *                    mutable random access iterators.
     *
     *  @exceptsafe Does not throw unless `coefficients[n]` does, and provides
     *              basic exception safety if it does.
     */
    template <typename Vector>
    void fill(double x, Vector && basis) const;

};

}}}} // namespace lsst::afw::math::polynomials

#endif // DOXYGEN
#endif // !LSST_AFW_MATH_POLYNOMIALS_Basis1d_h_INCLUDED
