#
# LSST Data Management System
# Copyright 2014 LSST Corporation.
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
from __future__ import absolute_import, division
from lsst.afw.geom import Box2D
from lsst.afw.cameraGeom import FOCAL_PLANE

class DetectorCollection(object):
    """An immutable collection of Detectors that can be accessed in various ways
    """
    def __init__(self, detectorList):
        """Construct a DetectorCollection
        
        @param[in] detectorList: a sequence of detectors in index order
        """
        self._detectorList = tuple(detectorList)
        self._nameDetectorDict = dict((d.getName(), d) for d in detectorList)
        self._fpBBox = Box2D()
        for detector in self._detectorList:
            for corner in detector.getCorners(FOCAL_PLANE):
                self._fpBBox.include(corner)

    def __iter__(self):
        """Return an iterator over all detectors in this collection"""
        return iter(self._detectorList)

    def __len__(self):
        """Return the number of detectors in this collection"""
        return len(self._detectorList)

    def __getitem__(self, key):
        """Return a detector given its name or index
        """
        if isinstance(key, basestring):
            return self._nameDetectorDict[key]
        else:
            return self._detectorList[key]
    
    def getFpBBox(self):
        """Return a focal plane bounding box that encompasses all detectors
        """
        return self._fpBBox