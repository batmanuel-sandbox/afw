import os.path
import lsst.afw.geom as afwGeom
from lsst.afw.table import AmpInfoCatalog
from lsst.afw.cameraGeom import FOCAL_PLANE, PUPIL, PIXELS, ACTUAL_PIXELS, CameraSys,\
                                Camera, Detector, Orientation, CameraTransformMap

__all__ = ["CameraFactoryTask"]

class CameraFactoryTask(object):
    """Make a camera

    Eventually we hope that camera data will be unpersisted using a butler,
    which is why this is written to look something like a Task.
    """
    cameraSysList = [PUPIL, FOCAL_PLANE, PIXELS, ACTUAL_PIXELS]
    cameraSysMap = dict((sys.getSysName(), sys) for sys in cameraSysList)

    def __init__(self):
        """Construct a CameraFactoryTask
        """
        pass

    def run(self, cameraConfig, ampInfoPath):
        """Construct a camera (lsst.afw.cameraGeom Camera)

        @param[in] cameraConfig: an instance of CameraConfig
        @param[in] ampInfoPath: path to find the persisted AmpInfoCatalogs
        @return camera (an lsst.afw.cameraGeom.Camera)
        """
        detectorList = []
        for detectorConfig in cameraConfig.detectorList.itervalues():
            #HACK until the mapper can get the short name
            nameEls = detectorConfig.name.split(" ")
            if len(nameEls[1]) == 7:
                nmap = {'A':'C0', 'B':'C1'}
                shortName = "R%s%s_S%s%s_%s"%(nameEls[0][2], nameEls[0][4], nameEls[1][2], nameEls[1][4], nmap[nameEls[1][6]])
            else:
                shortName = "R%s%s_S%s%s"%(nameEls[0][2], nameEls[0][4], nameEls[1][2], nameEls[1][4])
            ampCatPath = os.path.join(ampInfoPath, shortName + ".fits")
            ampInfoCat = AmpInfoCatalog.readFits(ampCatPath)
            detectorList.append(self.makeDetector(detectorConfig, ampInfoCat))
        nativeSys = self.cameraSysMap[cameraConfig.transformDict.nativeSys]
        transformDict = self.makeTransformDict(cameraConfig.transformDict.transforms)
        transformMap = CameraTransformMap(nativeSys, transformDict)
        return Camera(cameraConfig.name, detectorList, transformMap)

    def runCatDict(self, cameraConfig, ampInfoCatDict):
        """Construct a camera (lsst.afw.cameraGeom Camera)

        @param[in] cameraConfig: an instance of CameraConfig
        @param[in] ampInfoCatDict: a dictionary keyed on the detector name of AmpInfoCatalog objects
        @return camera (an lsst.afw.cameraGeom.Camera)
        """
        detectorList = []
        for detectorConfig in cameraConfig.detectorList.itervalues():
            detectorList.append(self.makeDetector(detectorConfig, ampInfoCatDict[detectorConfig.name]))
        nativeSys = self.cameraSysMap[cameraConfig.transformDict.nativeSys]
        transformDict = self.makeTransformDict(cameraConfig.transformDict.transforms)
        transformMap = CameraTransformMap(nativeSys, transformDict)
        return Camera(cameraConfig.name, detectorList, transformMap)

    def makeDetector(self, detectorConfig, ampInfoCatalog):
        """Make a detector object:

        @param detectorConfig -- config for this detector (an lsst.pex.config.Config)
        @param ampInfoCatalog -- amplifier information for this detector (an lsst.afw.table.AmpInfoCatalog)
        @return detector (an lsst.afw.cameraGeom.Detector)
        """
        orientation = self.makeOrientation(detectorConfig)
        pixelSize = afwGeom.Extent2D(detectorConfig.pixelSize_x, detectorConfig.pixelSize_y)
        transforms = self.makeTransformDict(detectorConfig.transformDict.transforms)
        transforms[FOCAL_PLANE] = orientation.makePixelFpTransform(pixelSize)
        llPoint = afwGeom.Point2I(detectorConfig.bbox_x0, detectorConfig.bbox_y0)
        urPoint = afwGeom.Point2I(detectorConfig.bbox_x1, detectorConfig.bbox_y1)
        bbox = afwGeom.Box2I(llPoint, urPoint)
        return Detector(detectorConfig.name, detectorConfig.detectorType, detectorConfig.serial, bbox, ampInfoCatalog, 
                             orientation, pixelSize, transforms)

    def makeOrientation(self, detectorConfig):
        """Make an instance of an Orientation class

        @param detectorConfig -- config for this detector (an lsst.pex.config.Config)
        @return orientation (an lsst.afw.cameraGeom.Orientation)
        """
        offset = afwGeom.Point2D(detectorConfig.offset_x, detectorConfig.offset_y)
        refPos = afwGeom.Point2D(detectorConfig.refpos_x, detectorConfig.refpos_y)
        yaw = afwGeom.Angle(detectorConfig.yawDeg, afwGeom.degrees)
        pitch = afwGeom.Angle(detectorConfig.pitchDeg, afwGeom.degrees)
        roll = afwGeom.Angle(detectorConfig.rollDeg, afwGeom.degrees)
        return Orientation(offset, refPos, yaw, pitch, roll)
        
    def makeTransformDict(self, transformConfigDict):
        """Make a dictionary of CameraSys: XYTransform.

        @param transformConfigDict -- an lsst.pex.config.ConfigDictField from an XYTransform registry;
            keys are camera system names.
        @return a dict of CameraSys or CameraSysPrefix: XYTransform
        """
        resMap = dict()
        if transformConfigDict is not None:
            for key in transformConfigDict:
                #TODO This needs to be handled by someone else.
                if key == "Pupil":
                    transform = afwGeom.InvertedXYTransform(transformConfigDict[key].transform.apply())
                else:
                    transform = transformConfigDict[key].transform.apply()
                resMap[CameraSys(key)] =  transform
        return resMap