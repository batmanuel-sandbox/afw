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
#ifndef LSST_AFW_MATH_POLYNOMIALS_PolynomialBasis1d_h_INCLUDED
#define LSST_AFW_MATH_POLYNOMIALS_PolynomialBasis1d_h_INCLUDED

#include "lsst/afw/math/polynomials/Scaling1d.h"
#include "lsst/afw/math/polynomials/RecurrenceBasis1d.h"
#include "lsst/afw/math/polynomials/ScaledBasis1d.h"

namespace lsst { namespace afw { namespace math { namespace polynomials {

/**
 *  A Recurrence for standard polynomials.
 *
 *  @see RecurrenceBasis1d.
 */
class PolynomialRecurrence {
public:

    static double getB0(double x) noexcept { return 1; }

    static double getB1(double x) noexcept { return x; }

    static double next(double x, std::size_t n, double current, double previous) noexcept {
        return current*x;
    }

};

/// A Basis1d for standard polynomials.
using PolynomialBasis1d = RecurrenceBasis1d<PolynomialRecurrence>;

/// A ScaledBasis1d for standard polynomials.
using ScaledPolynomialBasis1d = ScaledBasis1d<PolynomialBasis1d>;

/**
 *  Construct a ScaledPolynomialBasis1d that remaps the given interval to [-1, 1] before
 *  evaluating polynomials.
 *
 *  @param[in]  order    Maximum order of the basis (inclusive).
 *  @param[in]  min      Minimum point of the interval, mapped to -1.
 *  @param[in]  max      Maximum point of the interval, mapped to 1.
 */
inline ScaledPolynomialBasis1d makeScaledPolynomialBasis1d(
    std::size_t order,
    double min, double max
) noexcept {
    return ScaledPolynomialBasis1d(PolynomialBasis1d(order), makeUnitRangeScaling1d(min, max));
}

}}}} // namespace lsst::afw::math::polynomials

#endif // !LSST_AFW_MATH_POLYNOMIALS_PolynomialBasis1d_h_INCLUDED
