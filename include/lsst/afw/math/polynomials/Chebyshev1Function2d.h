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
#ifndef LSST_AFW_MATH_POLYNOMIALS_Chebyshev1Function2d_h_INCLUDED
#define LSST_AFW_MATH_POLYNOMIALS_Chebyshev1Function2d_h_INCLUDED

#include "lsst/afw/math/polynomials/Function2d.h"
#include "lsst/afw/math/polynomials/Chebyshev1Basis2d.h"

namespace lsst { namespace afw { namespace math { namespace polynomials {

/// A Function2d for Chebyshev polynomials of the first kind.
using Chebyshev1Function2d = Function2d<Chebyshev1Basis2d>;

/// A Function2d for scaled Chebyshev polynomials of the first kind.
using ScaledChebyshev1Function2d = Function2d<ScaledChebyshev1Basis2d>;

}}}} // namespace lsst::afw::math::polynomials

#endif // !LSST_AFW_MATH_POLYNOMIALS_Chebyshev1Function2d_h_INCLUDED
