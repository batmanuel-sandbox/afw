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
#ifndef LSST_AFW_MATH_polynomials_h_INCLUDED
#define LSST_AFW_MATH_polynomials_h_INCLUDED

#include "lsst/afw/math/polynomials/BinomialMatrix.h"
#include "lsst/afw/math/polynomials/PolynomialFunction1d.h"
#include "lsst/afw/math/polynomials/Chebyshev1Function1d.h"
#include "lsst/afw/math/polynomials/PolynomialFunction2d.h"
#include "lsst/afw/math/polynomials/Chebyshev1Function2d.h"

namespace lsst { namespace afw { namespace math {

/**
 *  @namespace lsst::afw::math::polynomials Low-level polynomials (including special polynomials) in C++.
 *
 *  The afw::math::polynomials library provides low-level classes for
 *  efficiently evaluating polynomial basis functions and expansions in C++.
 *  The classes here:
 *   - are not available in Python;
 *   - are not polymorphic (no virtual functions);
 *   - provide workspace objects to minimize memory allocations when
 *     appropriate;
 *   - do not throw exceptions (users are responsible for providing valid
 *     inputs);
 *   - have almost no outside dependencies (just Eigen);
 *   - use templates to allow them to work with any array/vector objects (not
 *     just Eigen).
 *
 *  They are intended to be used as the building blocks of higher-level
 *  objects that are visible to Python users and persisted with our data
 *  products, such as ChebyshevBoundedField and the afw::math::Function
 *  hierarchy (DM-13156 and DM-13157, respectively).
 *
 *  At present, the library only includes support for 1-d and 2-d standard
 *  polynomials and Chebyshev polynomials of the first kind, but adding
 *  support for any other function defined by a recurrence relation (i.e. any
 *  other special polynomial) should be extremely easy (see
 *  RecurrenceBasis1d), and need not be done within the polynomials library
 *  itself.
 *
 *  For both 1-d and 2-d, the library contains the following kinds of objects:
 *
 *   - Basis1d and Basis2d: objects that evaluate basis functions at one
 *     point at a time.  The fundamental 1-d basis is the RecurrenceBasis1d
 *     template, while the fundamental 2-d basis is the PackedBasis2d template.
 *     @ref PolynomialBasis1d, @ref Chebyshev1Basis1d, @ref PolynomialBasis2d, and
 *     @ref Chebyshev1Basis2d provide typedefs to instantiations of these that
 *     should generally be preferred.
 *
 *   - Function1d and Function2d: templates that combine a Basis1d or Basis2d
 *     with a vector of coefficients.
 *
 *   - Scaling1d and Scaling2d: transformation objects that, with the
 *     ScaledBasis1d and ScaledBasis2d templates, allow a basis or function to
 *     be constructed that remaps points as part of evaluating the basis
 *     functions.  For example, a ScaledChebyshev1Basis1d combines both the a
 *     Chebyshev polynomial basis and the scaling from some domain to [-1, 1]
 *     that enables most special properties of Chebyshevs.  When necessary,
 *     the @ref simplify functions can be used to convert a scaled polynomial
 *     function to an equivalent unscaled function.
 *
 *  The library also includes a few utility classes and functions:
 *
 *   - BinomialMatrix, which provides an efficient and stable way to get
 *     binomial coefficients.
 *
 *   - PackedIndexIterator and PackedIndexRange, which handle the flattening
 *     of pairs of 1-d polynomials coefficients and basis functions into 1-d
 *     arrays.
 */
namespace polynomials {

}}}} // namespace lsst::afw::math::polynomials

#endif // !LSST_AFW_MATH_polynomials_h_INCLUDED