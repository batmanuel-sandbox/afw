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
Tests for Angle

Run with:
   angle.py
or
   python
   >>> import angle; angle.run()
"""
import itertools
import math
import unittest

import numpy as np

import lsst.utils.tests
import lsst.afw.geom as afwGeom


class AngleTestCase(unittest.TestCase):
    """A test case for Angle"""

    def setUp(self):
        self.pi = afwGeom.Angle(math.pi, afwGeom.radians)
        self.d = 180*afwGeom.degrees

    def testCtor(self):
        self.assertEqual(self.pi, math.pi)
        self.assertEqual(self.pi, afwGeom.Angle(math.pi))
        self.assertEqual(self.pi, self.d)

        dd = afwGeom.Angle(180, afwGeom.degrees)
        self.assertEqual(self.d, dd)
        dd = afwGeom.Angle(60*180, afwGeom.arcminutes)
        self.assertEqual(self.d, dd)
        dd = afwGeom.Angle(60*60*180, afwGeom.arcseconds)
        self.assertEqual(self.d, dd)

    def testArithmetic(self):
        self.assertTrue(afwGeom.isAngle(self.pi))
        self.assertFalse(afwGeom.isAngle(self.pi.asRadians()))
        self.assertFalse(afwGeom.isAngle(math.pi))

        with self.assertRaises(TypeError):
            self.pi - math.pi           # subtracting a float from an Angle
        self.assertEqual(self.pi - math.pi*afwGeom.radians, 0)
        self.assertEqual(self.pi - self.d, 0)  # can subtract Angles

        with self.assertRaises(TypeError):
            self.pi + math.pi           # adding a float to an Angle

        with self.assertRaises(NotImplementedError):
            self.pi*afwGeom.degrees     # self.pi is already an Angle

        self.assertEqual((self.pi + self.d).asAngularUnits(afwGeom.degrees),
                         360)
        self.assertEqual((self.pi).asRadians(), math.pi)
        self.assertEqual((self.pi/2).asDegrees(), 90)
        self.assertEqual((self.pi*2).asArcminutes(), 360*60)
        self.assertEqual((self.pi*2).asArcseconds(), 360*60*60)
        self.assertEqual((-self.pi).asRadians(), -math.pi)

        with self.assertRaises(TypeError):
            2.0 / self.pi           # dividing a float by an Angle
        with self.assertRaises(TypeError):
            self.pi / self.pi

        # automatic conversion to double
        self.assertEqual(math.sin(self.pi/2), 1.0)

    def testAbs(self):
        self.assertEqual(abs(0.0*afwGeom.degrees - self.pi), self.pi)

    def testPi(self):
        self.assertEqual(afwGeom.PI, math.pi)

    def testComparison(self):
        a2 = 2.0 * afwGeom.arcseconds
        a1 = 0.5 * afwGeom.arcseconds
        a3 = 0.5 * afwGeom.arcseconds
        self.assertEqual(a1, a3)
        self.assertNotEqual(a1, a2)
        self.assertLessEqual(a1, a2)
        self.assertLess(a1, a2)
        self.assertGreater(a2, a1)
        self.assertGreaterEqual(a2, a1)

        self.assertFalse(a1 != a3)
        self.assertFalse(a1 == a2)
        self.assertFalse(a1 >= a2)
        self.assertFalse(a1 > a2)
        self.assertFalse(a2 < a1)
        self.assertFalse(a2 <= a1)

        self.assertTrue(a1 == float(a1))
        self.assertTrue(float(a1) == a1)

    def testTrig(self):
        self.assertEqual(math.cos(self.d), -1.0)
        self.assertAlmostEqual(math.sin(self.d), 0.0, places=15)
        thirty = 30.*afwGeom.degrees
        self.assertAlmostEqual(math.sin(thirty), 0.5, places=15)

    def testSeparation(self):
        """Tests whether angle differences are computed as expected.

        Wrapping accuracy is assumed tested by testWrap.
        """
        angleBase = 0.0*afwGeom.degrees
        angleWrap = 360.0*afwGeom.degrees
        angleHalf = -180.0*afwGeom.degrees
        angleOdd = 32.0*afwGeom.degrees

        self.checkWrappedAngle(angleBase.separation(angleWrap),
                               0.0*afwGeom.degrees)
        self.checkWrappedAngle(angleWrap.separation(angleBase),
                               0.0*afwGeom.degrees)

        self.checkWrappedAngle(angleBase.separation(angleHalf), angleHalf)
        self.checkWrappedAngle(angleHalf.separation(angleBase), angleHalf)

        self.checkWrappedAngle(angleWrap.separation(angleHalf), angleHalf)
        self.checkWrappedAngle(angleHalf.separation(angleWrap), angleHalf)

        self.checkWrappedAngle(angleOdd.separation(angleBase), angleOdd)
        self.checkWrappedAngle(angleBase.separation(angleOdd), -angleOdd)

        self.checkWrappedAngle(angleOdd.separation(angleWrap), angleOdd)
        self.checkWrappedAngle(angleWrap.separation(angleOdd), -angleOdd)

    def checkWrappedAngle(self, observed, expected):
        """Tests whether an angle wrapped to [-pi, pi) both matches its expected
        value and strictly satisfies its range restriction.
        """
        obs = observed.asRadians()
        exp = expected.asRadians()
        self.assertAlmostEqual(obs, exp, delta=np.finfo(float).eps)
        self.assertGreaterEqual(obs, -math.pi)
        self.assertLess(obs, math.pi)

    def testWrap(self):
        eps = np.finfo(float).eps
        self.assertNotEqual(1 + eps, eps)
        for wrap, offset, epsMult in itertools.product(
            (-1000, -10, -1, 0, 1, 10, 1000),
            (-2*math.pi, -math.pi, -math.pi*0.5, 0.0, math.pi*0.5, math.pi*0.75, math.pi, math.pi*2.0),
            (-3, -2, -1, 0, 1, 2, 3),
        ):
            angRad = (offset + (wrap*math.pi)) * (1 + (eps*epsMult))
            ang = angRad * afwGeom.radians
            angDeg = ang.asDegrees()
            sinAng = math.sin(angRad)
            cosAng = math.cos(angRad)

            posAng = (angRad * afwGeom.radians).wrap()
            posAngRad = posAng.asRadians()
            posAngDeg = posAng.asDegrees()
            posAngArcmin = posAng.asArcminutes()
            posAngArcsec = posAng.asArcseconds()
            # the code promises 0 <= posAng for all units
            self.assertGreaterEqual(posAngRad, 0)
            self.assertGreaterEqual(posAngDeg, 0)
            self.assertGreaterEqual(posAngArcmin, 0)
            self.assertGreaterEqual(posAngArcsec, 0)
            # wrap promises posAng < 2*pi only for radians,
            # but it seems to work for all units
            self.assertLess(posAngRad, 2*math.pi)
            self.assertLess(posAngDeg, 360)
            self.assertLess(posAngArcmin, 360 * 60)
            self.assertLess(posAngArcsec, 360 * 3600)
            # prove that posAngDeg and angDeg are the same angle
            posErrAng = ((posAngDeg - angDeg) * afwGeom.degrees) \
                .wrapCtr()
            self.assertAlmostEqual(posErrAng.asDegrees(), 0)
            # a sanity check in case wrapCtr gives the wrong answer
            self.assertAlmostEqual(math.sin(posAngRad), sinAng)
            self.assertAlmostEqual(math.cos(posAngRad), cosAng)

            ctrAng = (angRad * afwGeom.radians).wrapCtr()
            ctrAngRad = ctrAng.asRadians()
            ctrAngDeg = ctrAng.asDegrees()
            ctrAngArcmin = ctrAng.asArcminutes()
            ctrAngArcsec = ctrAng.asArcseconds()
            # wrapCtr promises -pi <= ctrAngRad < pi only for radians,
            # but it seems to work for all units
            self.assertGreaterEqual(ctrAngRad, -math.pi)
            self.assertGreaterEqual(ctrAngDeg, -180)
            self.assertGreaterEqual(ctrAngArcmin, -180 * 60)
            self.assertGreaterEqual(ctrAngArcsec, -180 * 3600)
            self.assertLess(ctrAngRad, math.pi)
            self.assertLess(ctrAngDeg, 180)
            self.assertLess(ctrAngArcmin, 180 * 60)
            self.assertLess(ctrAngArcsec, 180 * 3600)
            # prove that ctrAngDeg and ang are the same angle
            ctrErrAng = ((ctrAngDeg - angDeg) * afwGeom.degrees) \
                .wrapCtr()
            self.assertAlmostEqual(ctrErrAng.asDegrees(), 0)
            self.assertAlmostEqual(math.sin(ctrAngRad), sinAng)
            self.assertAlmostEqual(math.cos(ctrAngRad), cosAng)

            for refAngBase, refAngWrap, refEpsMult in itertools.product(
                (-math.pi, 0.0, math.pi, math.pi*2.0),
                range(-10, 11, 2),
                (-3, -2, -1, 0, 1, 2, 3),
            ):
                refAngRad = refAngBase * (1 + (eps * refEpsMult)) + refAngWrap * math.pi * 2.0
                refAng = refAngRad * afwGeom.radians
                refAngDeg = refAng.asDegrees()
                refAngArcmin = refAng.asArcminutes()
                refAngArcsec = refAng.asArcseconds()
                nearAng = (angRad * afwGeom.radians).wrapNear(refAng)
                nearAngRad = nearAng.asRadians()
                nearAngDeg = nearAng.asDegrees()
                nearAngArcmin = nearAng.asArcminutes()
                nearAngArcsec = nearAng.asArcseconds()

                fudgeFactor = 5  # 1 and 2 are too small for one case, 3 works; 5 has margin
                relEps = max(1, nearAngRad / math.pi, refAngRad / math.pi) * eps * fudgeFactor
                piWithSlop = math.pi * (1 + relEps)
                oneEightyWithSlop = 180 * (1 + relEps)

                # wrapNear promises nearAngRad - refAngRad >= -pi
                # for radians but has known failures due to
                # roundoff error for other units and for < pi
                self.assertGreaterEqual(nearAngRad - refAngRad, -math.pi)
                self.assertGreaterEqual(nearAngDeg - refAngDeg, -oneEightyWithSlop)
                self.assertGreaterEqual(nearAngArcmin - refAngArcmin, -oneEightyWithSlop * 60)
                self.assertGreaterEqual(nearAngArcsec - refAngArcsec, -oneEightyWithSlop * 3600)
                self.assertLessEqual(nearAngRad - refAngRad, piWithSlop)
                self.assertLessEqual(nearAngDeg - refAngDeg, oneEightyWithSlop)
                self.assertLessEqual(nearAngArcmin - refAngArcmin, oneEightyWithSlop * 60)
                self.assertLessEqual(nearAngArcsec - refAngArcsec, oneEightyWithSlop * 3600)
                # prove that nearAng and ang are the same angle
                nearErrAng = ((nearAngRad - angRad) * afwGeom.radians).wrapCtr()
                self.assertAlmostEqual(nearErrAng.asRadians(), 0)
                self.assertAlmostEqual(math.sin(nearAngRad), sinAng)
                self.assertAlmostEqual(math.cos(nearAngRad), cosAng)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
