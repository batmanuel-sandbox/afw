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

#include "lsst/afw/math/polynomials/PolynomialFunction2d.h"
#include "lsst/afw/math/polynomials/BinomialMatrix.h"

namespace lsst { namespace afw { namespace math { namespace polynomials {

namespace {

Eigen::VectorXd computePowers(double x, int n) {
    Eigen::VectorXd r(n + 1);
    r[0] = 1.0;
    for (int i = 1; i <= n; ++i) {
        r[i] = r[i - 1]*x;
    }
    return r;
}

} // anonymous


PolynomialFunction2d simplify(ScaledPolynomialFunction2d const & f) {
    Eigen::VectorXd result = Eigen::VectorXd::Zero(f.getBasis().size());
    std::size_t const n = f.getBasis().getOrder();
    auto rPow = computePowers(f.getBasis().getScaling().getX().getScale(), n);
    auto sPow = computePowers(f.getBasis().getScaling().getY().getScale(), n);
    auto uPow = computePowers(f.getBasis().getScaling().getX().getShift(), n);
    auto vPow = computePowers(f.getBasis().getScaling().getY().getShift(), n);
    BinomialMatrix binomial(f.getBasis().getNested().getOrder());
    for (auto const & i : f.getBasis().getIndices()) {
        for (std::size_t j = 0; j <= i.x; ++j) {
            double tmp = binomial(i.x, j)*uPow[j] *
                f.getCoefficient(i.flat)*rPow[i.x]*sPow[i.y];
            for (std::size_t k = 0; k <= i.y; ++k) {
                result[f.getBasis().index(i.x - j, i.y - k)] +=
                    binomial(i.y, k)*vPow[k]*tmp;
            }
        }
    }
    return makeFunction2d(f.getBasis().getNested(), result);
}

}}}} // namespace lsst::afw::math::polynomials
