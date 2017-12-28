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
#ifndef LSST_AFW_MATH_POLYNOMIALS_ScaledBasis1d_h_INCLUDED
#define LSST_AFW_MATH_POLYNOMIALS_ScaledBasis1d_h_INCLUDED

#include "lsst/afw/math/polynomials/Scaling1d.h"

namespace lsst { namespace afw { namespace math { namespace polynomials {

template <typename Basis>
class Function1d;

/**
 *  A 1-d basis that transforms all input points before evaluating nested basis.
 *
 *  If the nested basis is defined by basis functions @f$B_n(x)@f$, the scaled
 *  basis functions are @f$B_n(S(x))@f$, where @f$S(x)@f$ is the scaling
 *  transform.
 *
 *  Both the nested basis and ScaledBasis1d itself are models of the Basis1d
 *  concept.
 */
template <typename Nested>
class ScaledBasis1d {
public:

    /// A Function1d object that uses this basis.
    using Function = Function1d<ScaledBasis1d>;

    /// The type returned by scale().
    using Scaled = ScaledBasis1d<Nested>;

    /// Construct a scaled basis from a nested basis and a scaling transform.
    explicit ScaledBasis1d(Nested const & nested, Scaling1d const & scaling) noexcept :
        _nested(nested),
        _scaling(scaling)
    {}

    /// Default copy constructor.
    ScaledBasis1d(ScaledBasis1d const &) = default;

    /// Default move constructor.
    ScaledBasis1d(ScaledBasis1d &&) = default;

    /// Default copy assignment.
    ScaledBasis1d & operator=(ScaledBasis1d const &) = default;

    /// Default move assignment.
    ScaledBasis1d & operator=(ScaledBasis1d &&) = default;

    /// Return the nested basis.
    Nested const & getNested() const noexcept { return _nested; }

    /// Return the scaling transform.
    Scaling1d const & getScaling() const noexcept { return _scaling; }

    /// Return the order of the basis.
    std::size_t getOrder() const noexcept { return getNested().getOrder(); }

    /// Return the number of elements in the basis.
    std::size_t size() const noexcept { return getNested().size(); }

    /**
     *  Return a further-scaled basis with the same order.
     *
     *  The scaled basis will transform all points by the given scaling
     *  before evaluating the basis functions in the same way as `this`.
     */
    Scaled scale(Scaling1d const & first) const {
        return getNested().scale(first.then(getScaling()));
    }

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
    double sumWith(double x, Vector const & coefficients) const {
        return getNested().sumWith(getScaling().applyForward(x), coefficients);
    }

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
    void fill(double x, Vector && basis) const {
        return getNested().fill(getScaling().applyForward(x), std::forward<Vector>(basis));
    }

private:
    Nested _nested;
    Scaling1d _scaling;
};

}}}} // namespace lsst::afw::math::polynomials

#endif // !LSST_AFW_MATH_POLYNOMIALS_ScaledBasis1d_h_INCLUDED
