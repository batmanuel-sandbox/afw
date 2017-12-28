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
#ifndef LSST_AFW_MATH_POLYNOMIALS_Function2d_h_INCLUDED
#define LSST_AFW_MATH_POLYNOMIALS_Function2d_h_INCLUDED

#include "Eigen/Core"
#include "lsst/afw/math/polynomials/ScaledBasis2d.h"

namespace lsst { namespace afw { namespace math { namespace polynomials {

/**
 *  A 2-d function defined by a series expansion and its coefficients.
 *
 *  A Function1d combines a Basis2d that defines basis functions @f$B_{n}(x, y)@f$
 *  with a flattened vector of associated coefficients @f$a_{n}@f$.  Evaluating the
 *  function computes
 *  @f[
 *      \sum_{n}^{n \le N} a_{n} B_{n}(x, y)
 *  @f]
 */
template <typename Basis_>
class Function2d {
public:

    /// The basis type used by this function.
    using Basis = Basis_;

    /// Type returned by makeWorkspace().
    using Workspace = typename Basis::Workspace;

    /// Construct with zero-valued coefficients.
    explicit Function2d(Basis const & basis) :
        _basis(basis),
        _coefficients(Eigen::VectorXd::Zero(basis.size()))
    {}

    /// Construct with coefficients from an Eigen object.
    explicit Function2d(Basis const & basis, Eigen::VectorXd const & coefficients) :
        _basis(basis),
        _coefficients(coefficients)
    {
        assert(basis.size() == static_cast<std::size_t>(_coefficients.size()));
    }

    /// Construct with coefficients from an STL iterator range.
    template <typename Iterator>
    explicit Function2d(Basis const & basis, Iterator first, Iterator last) :
        _basis(basis),
        _coefficients(basis.size())
    {
        std::copy(first, last, &_coefficients[0]);
    }

    /// Default copy constructor.
    Function2d(Function2d const &) = default;

    /// Default move constructor.
    Function2d(Function2d &&) = default;

    /// Default copy assignment.
    Function2d & operator=(Function2d const &) = default;

    /// Default move assignment.
    Function2d & operator=(Function2d &&) = default;

    /// Return the associated Basis2d object.
    Basis const & getBasis() const { return _basis; }

    /// Allocate workspace that can be passed to operator() to avoid repeated memory allocations.
    Workspace makeWorkspace() const { return _basis.makeWorkspace(); }

    /// Evaluate the function at the given point.
    double operator()(geom::Point2D const & point) const {
        return _basis.sumWith(point, _coefficients);
    }

    /// Evaluate the function at the given point.
    double operator()(geom::Point2D const & point, Workspace & workspace) const {
        return _basis.sumWith(point, _coefficients, workspace);
    }

    /// Return the coefficient with the given (possibly flattened) index.
    double getCoefficient(std::size_t i) const { return _coefficients[i]; }

    /// Return the coefficient with the given (possibly flattened) index.
    void setCoefficient(std::size_t i) const { return _coefficients[i]; }

    /// Return a new function that applies the given scaling to all points before evaluation.
    Function2d<typename Basis::Scaled> scale(Scaling2d const & scaling) const {
        return Function2d<typename Basis::Scaled>(getBasis().scale(scaling), _coefficients);
    }

private:
    Basis _basis;
    Eigen::VectorXd _coefficients;
};

/// Create a Function2d of the appropriate type from a Basis2d and an Eigen object containing coefficients.
template <typename Basis>
Function2d<Basis> makeFunction2d(Basis const & basis, Eigen::VectorXd const & coefficients) {
    return Function2d<Basis>(basis, coefficients);
}

/// Create a Function2d of the appropriate type from a Basis2d and an iterator range containing coefficients.
template <typename Basis, typename Iterator>
Function2d<Basis> makeFunction2d(Basis const & basis, Iterator first, Iterator last) {
    return Function2d<Basis>(basis, first, last);
}

}}}} // namespace lsst::afw::math::polynomials

#endif // !LSST_AFW_MATH_POLYNOMIALS_Function2d_h_INCLUDED
