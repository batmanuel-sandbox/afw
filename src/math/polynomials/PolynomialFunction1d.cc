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

#include "lsst/afw/math/polynomials/PolynomialFunction1d.h"
#include "lsst/afw/math/polynomials/BinomialMatrix.h"

namespace lsst { namespace afw { namespace math { namespace polynomials {

PolynomialFunction1d simplify(ScaledPolynomialFunction1d const & f) {
    Eigen::VectorXd result = Eigen::VectorXd::Zero(f.getBasis().size());
    double const s = f.getBasis().getScaling().getScale();
    double const v = f.getBasis().getScaling().getShift();
    double sn = 1; // s^n
    BinomialMatrix binomial(f.getBasis().getNested().getOrder());
    for (std::size_t n = 0; n < f.getBasis().size(); ++n, sn *= s) {
        double vk = 1; // v^k
        for (std::size_t k = 0; k <= n; ++k, vk *= v) {
            result[n - k] += sn*binomial(n, k)*f.getCoefficient(n)*vk;
        }
    }
    return makeFunction1d(f.getBasis().getNested(), result);
}

}}}} // namespace lsst::afw::math::polynomials
