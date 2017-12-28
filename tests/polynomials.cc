/*
 * LSST Data Management System
 * Copyright 2008-2017 LSST/AURA.
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
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE polynomials

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-variable"
#include "boost/test/unit_test.hpp"
#pragma clang diagnostic pop

#include <vector>
#include <typeinfo>
#include "Eigen/Core"

#include "lsst/afw/math/polynomials.h"

using namespace lsst::afw::math::polynomials;
using lsst::afw::geom::Point2D;
using lsst::afw::geom::Box2D;

namespace {

constexpr double DEFAULT_RTOL = 2*std::numeric_limits<double>::epsilon();

// Wrapper around getType that returns a comparable object with operator<< support,
// so we can use it with BOOST_CHECK_EQUAL
template <typename T>
std::string getType(T const & v) {
    return std::string(typeid(T).name());
}

inline bool compare(double a, double b, double rtol) {
    if (!(std::abs(a - b) <= std::sqrt(a*b)*rtol)) {
        std::cerr << (boost::format("a=%0.16g, b=%0.16g, diff=%0.16g > %0.16g")
                      % a % b % std::abs(a - b) % (std::sqrt(a*b)*rtol)) << std::endl;
        return false;
    }
    return true;
}

inline bool compare(Point2D const & a, Point2D const & b, double rtol) {
    return compare(a.getX(), b.getX(), rtol) && compare(a.getY(), b.getY(), rtol);
}

inline bool compare(std::vector<double> const & a, std::vector<double> const & b, double rtol) {
    if (a.size() != b.size()) {
        return false;
    }
    for (std::size_t i = 0; i != a.size(); ++i) {
        if (!compare(a[i], b[i], rtol)) {
            return false;
        }
    }
    return true;
}

#define CUSTOM_CHECK_CLOSE(a, b, rtol) \
    do {                               \
        auto a1 = a;                   \
        auto b1 = b;                   \
        BOOST_CHECK( compare(a1, b1, rtol) );   \
    } while (false)

void testBinomialMatrix(int n) {
    auto factorial = [](int m) {
        double result = 1.0;
        for (int k = 1; k <= m; ++k) {
            result *= k;
        }
        return result;
    };
    BinomialMatrix binomial(n);
    for (int k = 0; k <= n; ++k) {
        CUSTOM_CHECK_CLOSE(binomial(n, k), factorial(n)/(factorial(k)*factorial(n - k)), DEFAULT_RTOL);
    }
}

template <typename Basis, typename Point>
void testBasis(Basis const & basis, Point const & point, std::vector<double> const & coefficients) {

    BOOST_CHECK_EQUAL(basis.size(), coefficients.size());

    // Test that we can call sumWith on:
    //  1) STL random-access containers
    std::vector<double> const coefficients1(coefficients);
    double z1 = basis.sumWith(point, coefficients1);
    //  2) STL random-access iterators
    std::vector<double> const coefficients2(coefficients);
    double z2 = basis.sumWith(point, coefficients2.begin());
    //  3) Eigen objects
    Eigen::VectorXd coefficients3(basis.size());
    std::copy(coefficients.begin(), coefficients.end(), &coefficients3[0]);
    double z3 = basis.sumWith(point, coefficients3);
    //  4) Eigen expressions
    Eigen::VectorXd coefficients4(basis.size() + 1);
    std::copy(coefficients1.begin(), coefficients1.end(), &coefficients4[0]);
    coefficients4[basis.size()] = std::numeric_limits<double>::quiet_NaN(); // should be ignored
    double z4 = basis.sumWith(point, coefficients4.head(basis.size())*2.0) / 2.0;
    // All of those evaluations should give the same results.
    BOOST_CHECK_EQUAL(z1, z2);
    BOOST_CHECK_EQUAL(z1, z3);
    BOOST_CHECK_EQUAL(z1, z4);

    // Test that we can call fill on the same:
    //  1) STL random-access containers
    std::vector<double> basis1(basis.size(), 0.0);
    basis.fill(point, basis1);
    //  2) STL random-access iterators
    std::vector<double> basis2(basis.size(), 0.0);
    basis.fill(point, basis2.begin());
    //  3) Eigen objects
    Eigen::VectorXd basis3(basis.size());
    basis.fill(point, basis3);
    //  4) Eigen *view* expressions (can't assign to rvalue expressions, of course)
    Eigen::VectorXd basis4(basis.size() + 1);
    basis4[basis.size()] = std::numeric_limits<double>::quiet_NaN(); // should be ignored
    basis.fill(point, basis4.head(basis.size()));
    // All of those basis vectors should have the same values
    BOOST_CHECK(std::equal(basis1.begin(), basis1.end(), basis2.begin()));
    BOOST_CHECK(std::equal(basis1.begin(), basis1.end(), &basis3[0]));
    BOOST_CHECK(std::equal(basis1.begin(), basis1.end(), &basis4[0]));

    // sum(basis * coefficients) should be equal to sumWith(), subject to round-off error
    CUSTOM_CHECK_CLOSE(basis3.dot(coefficients3), z1, 5*DEFAULT_RTOL);

    // Test using Function object to do the evaluation.
    typename Basis::Function func(basis, coefficients.begin(), coefficients.end());
    double z5 = func(point);
    BOOST_CHECK_EQUAL(z1, z5);

}

template <typename Basis, typename Point, typename Scaling>
void testScaledBasis(Basis const & basis, Point const & point, std::vector<double> const & coefficients,
                          Scaling const & scaling) {
    auto tBasis = basis.scale(scaling);
    auto tPoint = scaling.applyForward(point);

    // Run the regular basis tests on the scaled basis.
    testBasis(tBasis, tPoint, coefficients);

    double z1 = basis.sumWith(tPoint, coefficients);
    double z2 = tBasis.sumWith(point, coefficients);
    BOOST_CHECK_EQUAL(z1, z2);

    std::vector<double> basis1(basis.size(), 0.0);
    std::vector<double> basis2(tBasis.size(), 0.0);
    basis.fill(tPoint, basis1);
    tBasis.fill(point, basis2);

    CUSTOM_CHECK_CLOSE(basis1, basis2, 2*DEFAULT_RTOL);

    // Test using Function object to do the scaling and evaluation.
    typename Basis::Function func(basis, coefficients.begin(), coefficients.end());
    auto tFunc = func.scale(scaling);
    double z3 = tFunc(point);
    BOOST_CHECK_EQUAL(z1, z3);
}

} // anonymous

BOOST_AUTO_TEST_CASE(PackedIndex) {
    int const order = 6;
    PackedIndexRange const range(PackedIndexIterator(), PackedIndexIterator::makeEnd(order));
    std::size_t n = 0;
    std::size_t count = 0;
    for (auto const & index : range) {
        BOOST_CHECK_EQUAL(index.flat, PackedIndexRange::computeIndex(index.x, index.y));
        BOOST_CHECK(index.x + index.y >= n); // order is strictly increasing throughout iteration.
        n = index.x + index.y;
        BOOST_CHECK(n <= order);
        if (index.x == 0u) {
            BOOST_CHECK_EQUAL(index.flat, PackedIndexRange::computeOffset(n));
        }
        ++count;
    }
    BOOST_CHECK_EQUAL(count, range.size());
    BOOST_CHECK_EQUAL(count, PackedIndexRange::computeSize(order));
}

BOOST_AUTO_TEST_CASE(scalings1d) {
    Scaling1d affine(2, -0.5);
    auto inverse = affine.invert();
    auto identity = affine.then(inverse);
    for (double x = -0.5; x < 2; x += 0.3) {
        CUSTOM_CHECK_CLOSE(identity.applyForward(x), x, DEFAULT_RTOL);
        CUSTOM_CHECK_CLOSE(identity.applyInverse(x), x, DEFAULT_RTOL);
        double y = affine.applyForward(x);
        BOOST_CHECK_EQUAL(y, (x + affine.getShift())*affine.getScale());
        CUSTOM_CHECK_CLOSE(affine.applyInverse(y), x, DEFAULT_RTOL);
        CUSTOM_CHECK_CLOSE(inverse.applyForward(y), x, DEFAULT_RTOL);
        CUSTOM_CHECK_CLOSE(inverse.applyInverse(x), y, DEFAULT_RTOL);
    }

    double min=-0.5, max=2.0;
    auto toUnitRange = makeUnitRangeScaling1d(min, max);
    CUSTOM_CHECK_CLOSE(toUnitRange.applyForward(min), -1.0, DEFAULT_RTOL);
    CUSTOM_CHECK_CLOSE(toUnitRange.applyForward(max), 1.0, DEFAULT_RTOL);
    CUSTOM_CHECK_CLOSE(toUnitRange.applyInverse(-1.0), min, DEFAULT_RTOL);
    CUSTOM_CHECK_CLOSE(toUnitRange.applyInverse(1.0), max, DEFAULT_RTOL);
}

BOOST_AUTO_TEST_CASE(binomials) {
    testBinomialMatrix(3);
    testBinomialMatrix(5);
}

BOOST_AUTO_TEST_CASE(basis1d) {
    std::vector<double> coefficients = { 4.2, 1.6, -3.0, 0.2, -1.1, 0.8 };
    Scaling1d scaling(2.0, -1.0);
    double point = 1.5;
    double min = -0.5, max=2.0;
    int order = 5;

    // regular bases
    testBasis(PolynomialBasis1d(order), point, coefficients);
    testBasis(Chebyshev1Basis1d(order), point, coefficients);

    // scaled once
    testScaledBasis(PolynomialBasis1d(order), point, coefficients, scaling);
    testScaledBasis(Chebyshev1Basis1d(order), point, coefficients, scaling);

    // scaled twice
    testScaledBasis(makeScaledPolynomialBasis1d(order, min, max), point, coefficients, scaling);
    testScaledBasis(makeScaledChebyshev1Basis1d(order, min, max), point, coefficients, scaling);
}

BOOST_AUTO_TEST_CASE(basis2d) {
    std::vector<double> coefficients = { 4.2, 1.6, -3.0, 0.2, -1.1, 0.8 };
    Scaling2d scaling(
        Scaling1d(2.0, -1.0),
        Scaling1d(0.8, 0.6)
    );
    Point2D point(1.5, -0.3);
    Box2D box(Point2D(-4.0, -3.5), Point2D(2.2, 1.8));
    int order = 2;

    // regular bases
    testBasis(PolynomialBasis2d(order), point, coefficients);
    testBasis(Chebyshev1Basis2d(order), point, coefficients);

    // scaled once
    testScaledBasis(PolynomialBasis2d(order), point, coefficients, scaling);
    testScaledBasis(Chebyshev1Basis2d(order), point, coefficients, scaling);

    // scaled twice
    testScaledBasis(makeScaledPolynomialBasis2d(order, box), point, coefficients, scaling);
    testScaledBasis(makeScaledChebyshev1Basis2d(order, box), point, coefficients, scaling);
}

BOOST_AUTO_TEST_CASE(simplify1d) {
    std::vector<double> const coefficients({ 4.2, 1.6, -3.0, 0.2, -1.1, 0.8 });
    double min = -0.5, max=2.0;
    auto sfunc = ScaledPolynomialFunction1d(
        makeScaledPolynomialBasis1d(5, min, max),
        coefficients.begin(), coefficients.end()
    );
    auto func = simplify(sfunc);
    for (double x = min; x < max; x += 0.3) {
        CUSTOM_CHECK_CLOSE(sfunc(x), func(x), 2*DEFAULT_RTOL);
    }
}

BOOST_AUTO_TEST_CASE(simplify2d) {
    std::vector<double> const coefficients({ 4.2, 1.6, -3.0, 0.2, -1.1, 0.8, 1.2, 0.7, 1.9, -0.6 });
    Box2D box(Point2D(-4.0, -3.5), Point2D(2.2, 1.8));
    auto sfunc = ScaledPolynomialFunction2d(
        makeScaledPolynomialBasis2d(4, box),
        coefficients.begin(), coefficients.end()
    );
    auto func = simplify(sfunc);
    for (double x = box.getMinX(); x < box.getMaxX(); x += 0.3) {
        for (double y = box.getMinY(); y < box.getMaxY(); y += 0.3) {
            Point2D point(x, y);
            CUSTOM_CHECK_CLOSE(sfunc(point), func(point), 100*DEFAULT_RTOL);
        }
    }
}
