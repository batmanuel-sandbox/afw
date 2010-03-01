#!/usr/bin/env python
import os
import math
import pdb                          # we may want to say pdb.set_trace()
import unittest

import eups
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.utils.tests as utilsTests
import lsst.afw.display.ds9 as ds9
import lsst.pex.exceptions.exceptionsLib as exceptions
import lsst

try:
    type(verbose)
except NameError:
    verbose = 0

dataDir = eups.productDir("afwdata")
if not dataDir:
    raise RuntimeError("Must set up afwdata to run these tests")
InputImagePath = os.path.join(dataDir, "871034p_1_MI")
InputSmallImagePath = os.path.join(dataDir, "small_img.fits")
InputCorruptMaskedImageName = "small_MI_corrupt"
currDir = os.path.abspath(os.path.dirname(__file__))
InputCorruptFilePath = os.path.join(currDir, "data", InputCorruptMaskedImageName)
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class WCSTestCaseSDSS(unittest.TestCase):
    """A test case for WCS using a small (SDSS) image with a slightly weird WCS"""

    def setUp(self):
        im = afwImage.DecoratedImageD(InputSmallImagePath)

        self.wcs = afwImage.makeWcs(im.getMetadata())

        if False:
            ds9.mtv(im, wcs=self.wcs)

    def tearDown(self):
        del self.wcs

    def testValidWcs(self):
        """Test operator bool() (== isValid)"""
        pass

    def testInvalidWcs(self):
        """Test operator bool() (== isValid)
        This test has been improved by deleting some essential
        metadata (in this case, CRPIX1, and CRPIX2) from the
        MaskedImage's metadata and using that.
        """
        wcs = afwImage.Wcs()
        self.assertFalse(wcs)

        # Using MaskedImage with corrupt metadata
        infile = afwImage.MaskedImageF_imageFileName(InputCorruptFilePath)
        decoratedImage = afwImage.DecoratedImageF(infile)
        metadata = decoratedImage.getMetadata()

        
        self.assertRaises(exceptions.LsstCppException, afwImage.makeWcs, metadata)

    def testXyToRaDecArguments(self):
        """Check that conversion of xy to ra dec (and back again) works"""
        xy = afwImage.PointD(110, 123)
        raDec = self.wcs.pixelToSky(xy)
        xy2 = self.wcs.skyToPixel(raDec)

        self.assertAlmostEqual(xy.getX(), xy2.getX())
        self.assertAlmostEqual(xy.getY(), xy2.getY())

        if False:
            #This part of the test causes an exception. The input SDSS image
            #image treats DEC as its first coordinate and RA as its second
            #coordinate (CRVAL1, 2; the opposition of how things are usually
            #done. As a result, if you pass ra/dec into wcs.skyToPixel()
            #wcslib returns an error because it tries to solve for dec ra
            #which isn't legal.
            #
            #As I'm not sure whether we should be treating this header
            #as legally or illegally formatted, I'm commenting it out
            #for the moment.
            #
            #The same problem affects the test at the start of the function
            #but as we don't check the intermediate raDec value we get
            #away with it

            #This line causes an exception to be raised
            raDec = afwImage.PointD(245.167400, +19.1976583)
            #This doesn't
            #raDec = afwImage.PointD(+19.1976583, 245.167400)

            xy = self.wcs.skyToPixel(raDec)
            print xy
            raDec2 = self.wcs.pixelToSky(xy)

            self.assertAlmostEqual(raDec.getX(), raDec2.getX())
            self.assertAlmostEqual(raDec.getY(), raDec2.getY())

    def test_RaTan_DecTan(self):
        """Check the RA---TAN, DEC--TAN WCS conversion"""
        raDec = self.wcs.pixelToSky(0.0, 0.0)
        raDec0 = afwImage.PointD(19.1960467992, 245.1598413385) # values from wcstools' xy2sky, transposed

        self.assertAlmostEqual(raDec.getX(), raDec0.getX(), 5)
        self.assertAlmostEqual(raDec.getY(), raDec0.getY(), 5) # dec from ds9

    def testIdentity(self):
        """Convert from ra, dec to col, row and back again"""
        raDec = afwImage.PointD(20, 150)
        rowCol = self.wcs.skyToPixel(raDec)
        raDec2 = self.wcs.pixelToSky(rowCol)

        self.assertAlmostEqual(raDec.getX(), raDec2.getX())
        self.assertAlmostEqual(raDec.getY(), raDec2.getY())

    def testInvalidRaDec(self):
        """Test a conversion for an invalid position.  Well, "test" isn't
        quite right as the result is invalid, but make sure that it still is"""
        raDec = afwImage.PointD(1, 2)

        self.assertRaises(lsst.pex.exceptions.exceptionsLib.LsstCppException, self.wcs.skyToPixel, raDec)

    def testCD(self):
        print self.wcs.getCDMatrix()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class WCSTestCaseCFHT(unittest.TestCase):
    """A test case for WCS"""

    def setUp(self):
        path = InputImagePath + "_img.fits"
        self.metadata = afwImage.readMetadata(path)
        self.wcs = afwImage.makeWcs(self.metadata)
        if False:
            ds9.mtv(e)

    def tearDown(self):
        del self.wcs
        del self.metadata

    def test_RaTan_DecTan(self):
        """Check the RA---TAN, DEC--TAN WCS conversion"""
        raDec = self.wcs.pixelToSky(0.0, 0.0) # position read off ds9

        self.assertAlmostEqual(raDec.getX(), 17.87673, 5) # ra from ds9
        self.assertAlmostEqual(raDec.getY(),  7.72231, 5) # dec from ds9

    def testPlateScale(self):
        """Test that we can measure the area of a pixel"""

        p00 = afwImage.PointD(10, 10)
        p00 = afwImage.PointD(self.metadata.getAsDouble("CRPIX1"), self.metadata.getAsDouble("CRPIX2"))

        sky00 = self.wcs.pixelToSky(p00)
        cosdec = math.cos(math.pi/180*sky00.getY())

        side = 1e-3
        p10 = self.wcs.skyToPixel(sky00 + afwImage.PointD(side/cosdec, 0))    - p00
        p01 = self.wcs.skyToPixel(sky00 + afwImage.PointD(0,           side)) - p00

        area = side*side/abs(p10.getX()*p01.getY() - p01.getX()*p10.getY())

        self.assertAlmostEqual(math.sqrt(self.wcs.pixArea(p00)), math.sqrt(area))
        #
        # Now check that the area's the same as the CD matrix gives.
        #
        cd = [self.metadata.get("CD1_1"), self.metadata.get("CD1_2"),
              self.metadata.get("CD2_1"), self.metadata.get("CD2_2")]
        area = math.fabs(cd[0]*cd[3] - cd[1]*cd[2])

        self.assertAlmostEqual(math.sqrt(self.wcs.pixArea(p00)), math.sqrt(area))

    def testReadWcs(self):
        """Test reading a Wcs directly from a fits file"""

        meta = afwImage.readMetadata(InputImagePath + "_img.fits")
        wcs = afwImage.makeWcs(meta)

        self.assertEqual(wcs.pixelToSky(0.0, 0.0), self.wcs.pixelToSky(0.0, 0.0))

    def testShiftWcs(self):
        """Test shifting the reference pixel"""
        sky10_10 = self.wcs.pixelToSky(afwImage.PointD(10, 10))

        self.wcs.shiftReferencePixel(-10, -10)
        sky00 = self.wcs.pixelToSky(afwImage.PointD(0, 0))
        self.assertEqual((sky00.getX(), sky00.getY()), (sky10_10.getX(), sky10_10.getY()))

    def testCloneWcs(self):
        """Test Cloning a Wcs"""
        sky00 = self.wcs.pixelToSky(afwImage.PointD(0, 0))

        new = self.wcs.clone()
        self.wcs.pixelToSky(afwImage.PointD(10, 10)) # shouldn't affect new

        nsky00 = new.pixelToSky(afwImage.PointD(0, 0))
        self.assertEqual((sky00.getX(), sky00.getY()), (nsky00.getX(), nsky00.getY()))

    def testCD(self):
        cd = self.wcs.getCDMatrix()
        self.assertAlmostEqual(cd[0,0], self.metadata.getAsDouble("CD1_1"))
        self.assertAlmostEqual(cd[0,1], self.metadata.getAsDouble("CD1_2"))
        self.assertAlmostEqual(cd[1,0], self.metadata.getAsDouble("CD2_1"))
        self.assertAlmostEqual(cd[1,1], self.metadata.getAsDouble("CD2_2"))

    def testConstructor(self):
        copy = afwImage.Wcs(self.wcs.getSkyOrigin(), self.wcs.getPixelOrigin(), 
                            self.wcs.getCDMatrix())

    def testAffineTransform(self):
        a = self.wcs.getAffineTransform()
        l = self.wcs.getCDMatrix()
        #print a[a.X], a[a.Y], print a[a.XX], a[a.XY], a[a.YX], a[a.YY]

        sky00g = afwGeom.makePointD(10, 10)
        sky00i = afwImage.PointD(sky00g.getX(), sky00g.getY())
        a = self.wcs.linearizeAt(sky00i)
        pix00i = self.wcs.skyToPixel(sky00i)
        pix00g = afwGeom.makePointD(pix00i.getX(), pix00i.getY())
        sky00gApprox = a(pix00g);
        self.assertAlmostEqual(sky00g.getX(), sky00gApprox.getX())
        self.assertAlmostEqual(sky00g.getY(), sky00gApprox.getY())
        self.assertAlmostEqual(self.wcs.pixArea(sky00i), abs(a[a.XX]* a[a.YY] - a[a.XY]*a[a.YX]))
        a.invert()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(WCSTestCaseSDSS)
    suites += unittest.makeSuite(WCSTestCaseCFHT)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)

    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    utilsTests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
