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
#ifndef LSST_AFW_MATH_POLYNOMIALS_Chebyshev1Basis2d_h_INCLUDED
#define LSST_AFW_MATH_POLYNOMIALS_Chebyshev1Basis2d_h_INCLUDED

#include "lsst/afw/math/polynomials/Chebyshev1Basis1d.h"
#include "lsst/afw/math/polynomials/PackedBasis2d.h"
#include "lsst/afw/math/polynomials/ScaledBasis2d.h"
#include "lsst/afw/math/polynomials/Scaling2d.h"

namespace lsst { namespace afw { namespace math { namespace polynomials {

/// A Basis2d for Chebyshev polynomials of the first kind.
using Chebyshev1Basis2d = PackedBasis2d<Chebyshev1Basis1d>;

/// A Basis2d for scaled Chebyshev polynomials of the first kind.
using ScaledChebyshev1Basis2d = ScaledBasis2d<Chebyshev1Basis2d>;

/**
 *  Construct a ScaledChebyshev1Basis2d that remaps the given box to [-1, 1]x[-1, 1] before
 *  evaluating Chebyshev polynomials.
 *
 *  @param[in]  order    Maximum order of the basis (inclusive).
 *  @param[in]  box      Box to be mapped to [-1, 1]x[-1, 1]
 */
inline ScaledChebyshev1Basis2d makeScaledChebyshev1Basis2d(std::size_t order, geom::Box2D const & box) {
    return ScaledChebyshev1Basis2d(Chebyshev1Basis2d(order), makeUnitRangeScaling2d(box));
}

}}}} // namespace lsst::afw::math::polynomials

#endif // !LSST_AFW_MATH_POLYNOMIALS_Chebyshev1Basis2d_h_INCLUDED
