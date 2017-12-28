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
#ifndef LSST_AFW_MATH_POLYNOMIALS_BinomialMatrix_h_INCLUDED
#define LSST_AFW_MATH_POLYNOMIALS_BinomialMatrix_h_INCLUDED

#include "Eigen/Core"

namespace lsst { namespace afw { namespace math { namespace polynomials {

/**
 *  A class that computes binomial coefficients up to a certain power.
 *
 *  The binomial coefficient is defined as:
 *  @f[
 *     \left(\begin{array}{ c }
 *       n \\
 *       k
 *     \end{array}\right)
 *     = \frac{n!}{k!(n-k)!}
 *  @f]
 *  with both @f$n@f$ and @f$k@f$ nonnegative integers and @f$k \le n@f$
 *
 *  This class uses recurrence relations to avoid computing factorials directly,
 *  making it both more efficient and numerically stable.
 */
class BinomialMatrix {
public:

    /**
     *  Construct an object that can compute binomial coefficients with @f$n@f$
     *  up to and including the given value.
     */
    explicit BinomialMatrix(int nMax) { extend(nMax); }

    /**
     *  Return the binomial coefficient.
     *
     *  No error checking is performed; the behavior of this method is is
     *  undefined if the given values do not satisfy
     *  @code
     *  n <= nMax && k <= n && n >=0 && k >= 0
     *  @endcode
     */
    double operator()(int n, int k) const {
        return getMatrix()(n, k);
    }

private:

    // Because binomial coefficients only depend on two ints, and we compute
    // them via a recurrence relation, we actually store all of the ones we
    // have every calculated in a static matrix.  We expand that matrix
    // when constructing a BinomialMatrix with nMax higher than we have seen
    // before, and wrap that expansion in a mutex for thread safety.
    // Reads from the matrix do not require going through the mutex.
    static Eigen::MatrixXd & getMatrix();

    // Ensure the static matrix supports n values up to and including the
    // given one, updating it if necessary.
    static void extend(int const n);
};

}}}} // namespace lsst::afw::math::polynomials

#endif // !LSST_AFW_MATH_POLYNOMIALS_BinomialMatrix_h_INCLUDED
