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
#ifndef LSST_AFW_MATH_POLYNOMIALS_Chebyshev1Basis1d_h_INCLUDED
#define LSST_AFW_MATH_POLYNOMIALS_Chebyshev1Basis1d_h_INCLUDED

#include "lsst/afw/math/polynomials/RecurrenceBasis1d.h"
#include "lsst/afw/math/polynomials/ScaledBasis1d.h"
#include "lsst/afw/math/polynomials/Scaling1d.h"

namespace lsst { namespace afw { namespace math { namespace polynomials {

/**
 *  A Recurrence for Chebyshev polynomials of the first kind.
 *
 *  @see RecurrenceBasis1d.
 */
class Chebyshev1Recurrence {
public:

    static double getB0(double x) { return 1; }

    static double getB1(double x) { return x; }

    static double next(double x, std::size_t n, double current, double previous) { return 2*x*current - previous; }

};

/// A Basis1d for Chebyshev polynomials of the first kind.
using Chebyshev1Basis1d = RecurrenceBasis1d<Chebyshev1Recurrence>;

/// A Basis1d for scaled Chebyshev polynomials of the first kind.
using ScaledChebyshev1Basis1d = ScaledBasis1d<Chebyshev1Basis1d>;

/**
 *  Construct a ScaledChebyshev1Basis1d that remaps the given interval to [-1, 1] before
 *  evaluating Chebyshev polynomials.
 *
 *  @param[in]  order    Maximum order of the basis (inclusive).
 *  @param[in]  min      Minimum point of the interval, mapped to -1.
 *  @param[in]  max      Maximum point of the interval, mapped to 1.
 */
inline ScaledChebyshev1Basis1d makeScaledChebyshev1Basis1d(std::size_t order, double min, double max) {
    return ScaledChebyshev1Basis1d(Chebyshev1Basis1d(order), makeUnitRangeScaling1d(min, max));
}

}}}} // namespace lsst::afw::math::polynomials

#endif // !LSST_AFW_MATH_POLYNOMIALS_Chebyshev1Basis1d_h_INCLUDED
