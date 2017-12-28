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
#ifndef LSST_AFW_MATH_POLYNOMIALS_Function1d_h_INCLUDED
#define LSST_AFW_MATH_POLYNOMIALS_Function1d_h_INCLUDED

#include "Eigen/Core"
#include "lsst/afw/math/polynomials/ScaledBasis1d.h"

namespace lsst { namespace afw { namespace math { namespace polynomials {

/**
 *  A 1-d function defined by a series expansion and its coefficients.
 *
 *  A Function1d combines a Basis1d that defines basis functions @f$B_n(x)@f$
 *  with a vector of associated coefficients @f$a_n@f$.  Evaluating the
 *  function computes
 *  @f[
 *      \sum_{n=0}^{n \le N} a_n B_n(x)
 *  @f]
 */
template <typename Basis_>
class Function1d {
public:

    /// The basis type used by this function.
    using Basis = Basis_;

    /// Construct with zero-valued coefficients.
    explicit Function1d(Basis const & basis) :
        _basis(basis),
        _coefficients(Eigen::VectorXd::Zero(basis.size()))
    {}

    /// Construct with coefficients from an Eigen object.
    Function1d(Basis const & basis, Eigen::VectorXd const & coefficients) :
        _basis(basis),
        _coefficients(coefficients)
    {
        assert(basis.size() == static_cast<std::size_t>(_coefficients.size()));
    }

    /// Construct with coefficients from an STL iterator range.
    template <typename Iterator>
    Function1d(Basis const & basis, Iterator first, Iterator last) :
        _basis(basis),
        _coefficients(basis.size())
    {
        std::copy(first, last, &_coefficients[0]);
    }

    /// Default copy constructor.
    Function1d(Function1d const &) = default;

    /// Default move constructor.
    Function1d(Function1d &&) = default;

    /// Default copy assignment.
    Function1d & operator=(Function1d const &) = default;

    /// Default move assignment.
    Function1d & operator=(Function1d &&) = default;

    /// Return the associated Basis1d object.
    Basis const & getBasis() const { return _basis; }

    /// Evaluate the function at the given point.
    double operator()(double x) const {
        return _basis.sumWith(x, _coefficients);
    }

    /// Return the coefficient associated with the nth basis function.
    double getCoefficient(std::size_t n) const { return _coefficients[n]; }

    /// Set the coefficient associated with the nth basis function.
    void setCoefficient(std::size_t n, double v) { _coefficients[n] = v; }

    /// Return a new function that applies the given scaling to all points before evaluation.
    Function1d<typename Basis::Scaled> scale(Scaling1d const & scaling) const {
        return Function1d<typename Basis::Scaled>(getBasis().scale(scaling), _coefficients);
    }

private:
    Basis _basis;
    Eigen::VectorXd _coefficients;
};

/// Create a Function1d of the appropriate type from a Basis1d and an Eigen object containing coefficients.
template <typename Basis>
Function1d<Basis> makeFunction1d(Basis const & basis, Eigen::VectorXd const & coefficients) {
    return Function1d<Basis>(basis, coefficients);
}

/// Create a Function1d of the appropriate type from a Basis1d and an iterator range containing coefficients.
template <typename Basis, typename Iterator>
Function1d<Basis> makeFunction1d(Basis const & basis, Iterator first, Iterator last) {
    return Function1d<Basis>(basis, first, last);
}


}}}} // namespace lsst::afw::math::polynomials

#endif // !LSST_AFW_MATH_POLYNOMIALS_Function1d_h_INCLUDED
