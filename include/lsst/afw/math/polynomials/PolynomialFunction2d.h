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
#ifndef LSST_AFW_MATH_POLYNOMIALS_PolynomialFunction2d_h_INCLUDED
#define LSST_AFW_MATH_POLYNOMIALS_PolynomialFunction2d_h_INCLUDED

#include "lsst/afw/math/polynomials/Function2d.h"
#include "lsst/afw/math/polynomials/PolynomialBasis2d.h"

namespace lsst { namespace afw { namespace math { namespace polynomials {

/// A Function2d for standard polynomials.
using PolynomialFunction2d = Function2d<PolynomialBasis2d>;

/// A Function2d for scaled standard polynomials.
using ScaledPolynomialFunction2d = Function2d<ScaledPolynomialBasis2d>;

/**
 *  Calculate the standard polynomial function that is equivalent to a scaled
 *  standard polynomial function.
 *
 *  This operation is not numerically stable at high order, but when fitting
 *  standard polynomials, it is still much more stable to first fit in a scaled
 *  basis and then (if necessary) use `simplify` to compute the coefficients
 *  of an equivalent unscaled polynomial.
 */
PolynomialFunction2d simplify(ScaledPolynomialFunction2d const & f);

}}}} // namespace lsst::afw::math::polynomials

#endif // !LSST_AFW_MATH_POLYNOMIALS_PolynomialFunction2d_h_INCLUDED
