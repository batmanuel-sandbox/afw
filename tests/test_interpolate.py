#
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#

"""
Tests for Interpolate

Run with:
   ./Interpolate.py
or
   python
   >>> import Interpolate; Interpolate.run()
"""

import unittest

import numpy as np

import lsst.utils.tests
import lsst.afw.math as afwMath
import lsst.pex.exceptions as pexExcept


class InterpolateTestCase(lsst.utils.tests.TestCase):

    """A test case for Interpolate Linear"""

    def setUp(self):
        self.n = 10
        self.x = np.zeros(self.n, dtype=float)
        self.y1 = np.zeros(self.n, dtype=float)
        self.y2 = np.zeros(self.n, dtype=float)
        self.y0 = 1.0
        self.dydx = 1.0
        self.d2ydx2 = 0.5

        for i in range(0, self.n, 1):
            self.x[i] = i
            self.y1[i] = self.dydx*self.x[i] + self.y0
            self.y2[i] = self.d2ydx2*self.x[i] * \
                self.x[i] + self.dydx*self.x[i] + self.y0

        self.xtest = 4.5
        self.y1test = self.dydx*self.xtest + self.y0
        self.y2test = self.d2ydx2*self.xtest*self.xtest + self.dydx*self.xtest + self.y0

    def tearDown(self):
        del self.x
        del self.y1
        del self.y2

    def testLinearRamp(self):

        # === test the Linear Interpolator ============================
        # default is akima spline
        yinterpL = afwMath.makeInterpolate(self.x, self.y1)
        youtL = yinterpL.interpolate(self.xtest)

        self.assertEqual(youtL, self.y1test)

    def testNaturalSplineRamp(self):

        # === test the Spline interpolator =======================
        # specify interp type with the string interface
        yinterpS = afwMath.makeInterpolate(
            self.x, self.y1, afwMath.Interpolate.NATURAL_SPLINE)
        youtS = yinterpS.interpolate(self.xtest)

        self.assertEqual(youtS, self.y1test)

    def testAkimaSplineParabola(self):
        """test the Spline interpolator"""
        # specify interp type with the enum style interface
        yinterpS = afwMath.makeInterpolate(
            self.x, self.y2, afwMath.Interpolate.AKIMA_SPLINE)
        youtS = yinterpS.interpolate(self.xtest)

        self.assertEqual(youtS, self.y2test)

    def testConstant(self):
        """test the constant interpolator"""
        # [xy]vec:   point samples
        # [xy]vec_c: centered values
        xvec = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0])
        xvec_c = np.array(
            [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5])
        yvec = np.array([1.0, 2.4, 5.0, 8.4, 13.0,
                         18.4, 25.0, 32.6, 41.0, 50.6])
        yvec_c = np.array([1.0, 1.7, 3.7, 6.7, 10.7, 15.7,
                           21.7, 28.8, 36.8, 45.8, 50.6])

        interp = afwMath.makeInterpolate(
            xvec, yvec, afwMath.Interpolate.CONSTANT)

        for x, y in zip(xvec_c, yvec_c):
            self.assertAlmostEqual(interp.interpolate(x + 0.1), y)
            self.assertAlmostEqual(interp.interpolate(x), y)

        self.assertEqual(interp.interpolate(xvec[0] - 10), yvec[0])
        n = len(yvec)
        self.assertEqual(interp.interpolate(xvec[n - 1] + 10), yvec[n - 1])

        # test caching as we go backwards
        for x, y in reversed(list(zip(xvec_c, yvec_c))):
            self.assertAlmostEqual(interp.interpolate(x + 0.1), y)
            self.assertAlmostEqual(interp.interpolate(x), y)

        i = 2
        for x in np.arange(xvec_c[i], xvec_c[i + 1], 10):
            self.assertEqual(interp.interpolate(x), yvec_c[i])

    def testInvalidInputs(self):
        """Test that invalid inputs cause an abort"""

        with self.assertRaises(pexExcept.OutOfRangeError):
            afwMath.makeInterpolate(np.array([], dtype=float), np.array([], dtype=float),
                                    afwMath.Interpolate.CONSTANT)

        afwMath.makeInterpolate(np.array([0], dtype=float), np.array([1], dtype=float),
                                afwMath.Interpolate.CONSTANT)

        with self.assertRaises(pexExcept.OutOfRangeError):
            afwMath.makeInterpolate(np.array([0], dtype=float), np.array([1], dtype=float),
                                    afwMath.Interpolate.LINEAR)

    def testLookupMaxInterpStyle(self):
        for numPoints in range(1, 6):
            maxInterpStyle = afwMath.lookupMaxInterpStyle(numPoints)
            desiredMax = {
                1: afwMath.Interpolate.Style.CONSTANT,
                2: afwMath.Interpolate.Style.LINEAR,
                3: afwMath.Interpolate.Style.CUBIC_SPLINE,
                4: afwMath.Interpolate.Style.CUBIC_SPLINE,
            }.get(numPoints, afwMath.Interpolate.Style.AKIMA_SPLINE)
            self.assertEqual(maxInterpStyle, desiredMax)

        for badNumPoints in (-5, -1, 0):
            with self.assertRaises(pexExcept.InvalidParameterError):
                afwMath.lookupMaxInterpStyle(badNumPoints)

    def testLookupMinInterpPoints(self):
        for style in afwMath.Interpolate.Style.__members__.values():
            if style in (afwMath.Interpolate.Style.UNKNOWN, afwMath.Interpolate.Style.NUM_STYLES):
                with self.assertRaises(pexExcept.OutOfRangeError):
                    afwMath.lookupMinInterpPoints(style)
            else:
                minPoints = afwMath.lookupMinInterpPoints(style)
                desiredMin = {
                    afwMath.Interpolate.Style.CONSTANT: 1,
                    afwMath.Interpolate.Style.LINEAR: 2,
                    afwMath.Interpolate.Style.NATURAL_SPLINE: 3,
                    afwMath.Interpolate.Style.CUBIC_SPLINE: 3,
                    afwMath.Interpolate.Style.CUBIC_SPLINE_PERIODIC: 3,
                    afwMath.Interpolate.Style.AKIMA_SPLINE: 5,
                    afwMath.Interpolate.Style.AKIMA_SPLINE_PERIODIC: 5,
                }.get(style, None)
                if desiredMin is None:
                    self.fail("Unrecognized style: %s" % (style,))
                self.assertEqual(minPoints, desiredMin)

    def testStringToInterpStyle(self):
        for name, desiredStyle in afwMath.Interpolate.Style.__members__.items():
            if name in ("UNKNOWN", "NUM_STYLES"):
                with self.assertRaises(pexExcept.InvalidParameterError):
                    afwMath.stringToInterpStyle(name)
            else:
                style = afwMath.stringToInterpStyle(name)
                self.assertEqual(style, desiredStyle)

        for badName in ("BOGUS", ""):
            with self.assertRaises(pexExcept.InvalidParameterError):
                afwMath.stringToInterpStyle(badName)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
