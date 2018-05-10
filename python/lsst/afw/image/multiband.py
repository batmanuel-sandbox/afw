#
# LSST Data Management System
# Copyright 2008-2017 LSST/AURA.
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

__all__ = ["MultibandBase", "MultibandImage"]

from collections import OrderedDict

import numpy as np

from .image.image import ImageF
from lsst.afw.geom import Point2I, Box2I


class MultibandBase(object):
    """Base class for multiband objects

    The LSST stack has a number of image-like classes that have
    data in multiple bands that are stored as separate objects.
    Analyzing the data can be easier using a Multiband object that
    wraps the underlying data as a single data cube that can be sliced and
    updated as a single object.

    `MultibandBase` is designed to contain the most important universal
    methods for initializing, slicing, and extracting common parameters
    (such as the bounding box or XY0 position) to all of the single band classes.
    """
    def __init__(self, singles, slices=None, **kwargs):
        """Initialize a `MultibandBase` object

        Must be overloaded in derived classes.
        At a minimum, the derived class must be initialized with
        a `_singles` attribute that is a list of single band objects,
        a `_filters` attribute that is a list of filter
        names, a `_bbox` attribute is a `Box2I`,
        and a `_xy0` attribute that is a `Point2I`.

        Below are the only required `__init__` parameters for all
        dervied classes.

        Parameters
        ----------
        singles: list, `OrderedDict`, or `MultibandBase`
            Either a list of single band objects or an `OrderedDict` with
            filter names as keys and single band objects as values, or
            an instance of an object that inherits from `MultibandBase`.
        slices: list
            Slices used to extract a subset of an image
            (used when `singles` inherits from `MultibandBase` to
            create a new object with a slice of the original).
        """
        err = "Inherited classes from `MultibandBase` require an `__init__` method"
        raise NotImplementedError(err)

    @property
    def filters(self):
        """List of filter names for the single band objects
        """
        return self._filters

    @property
    def singles(self):
        """List of afw single band objects
        """
        return self._singles

    @property
    def bbox(self):
        """Bounding box
        """
        return self._bbox

    def getBBox(self):
        """Get the bounding box
        """
        return self._bbox

    @property
    def XY0(self):
        """Minimum coordinate in the bounding box
        """
        return self._xy0

    def setXY0(self, xy0):
        """Update the XY0 position for each single band object
        """
        assert isinstance(xy0, Point2I)
        self._xy0 = xy0
        for n in range(len(self.singles)):
            self.singles[n].setXY0(xy0)

    @property
    def x0(self):
        """X0

        X component of XY0 `Point2I.getX()`
        """
        return self.XY0.getX()

    @property
    def y0(self):
        """Y0

        Y component of XY0 `Point2I.getY()`
        """
        return self.XY0.getY()

    @property
    def yx0(self):
        """Minimum (y,x) position
        """
        return (self.y0, self.x0)

    @property
    def xy0(self):
        """Minimum (x,y) position
        """
        return (self.x0, self.y0)

    @property
    def width(self):
        """Width of the images
        """
        return self.XY0.getWidth()

    @property
    def height(self):
        """Height of the images
        """
        return self.XY0.getHeight()

    def __getitem__(self, args):
        """Get a slice of the underlying array

        If there is only a single slice, which is a string,
        then the element of `singles` for the given filter is
        returned.
        """
        if np.issubdtype(type(args), np.integer) or isinstance(args, slice):
            slices = (args,)
        else:
            slices = args
        if isinstance(slices[0], str):
            idx = self.filters.index(slices[0])
            return self.singles[idx]
        return type(self)(self, slices=slices)

    def __repr__(self):
        result = "<{0}, filters={1}, yx0={2}, image shape={3}>".format(
            self.__class__.__name__, self.filters, self.yx0, self.array.shape[-2:])
        return result

    def __str__(self):
        return str(self.array)


class MultibandImage(MultibandBase):
    """Multiband Image class

    This class acts as a container for multiple `afw.Image<X>` objects,
    where `<X>` is the data type (for example `ImageF`).
    All images must be contained in the same bounding box,
    and eb the same data type.
    """
    def __init__(self, singles=None, filters=None, slices=None, deep=False, imageType=ImageF,
                 filterKwargs=None, **kwargs):
        """Initialize a `MultibandBase` object

        If `singles` is another `MultibandImage`
        then a duplicate `MultibandImage` is made,
        sliced appropriately.

        Parameters
        ----------
        singles: list, `OrderedDict`, or `MultibandBase`
            Either a list of single band objects or an `OrderedDict` with
            filter names as keys and single band objects as values, or
            an instance of an object that inherits from `MultibandBase`.
        filters: list
            List of filter names. If `singles` is an `OrderedDict` or
            inherits from `MultibandBase` then this arguement is ignored,
            otherwise it is required.
        slices: list
            Slices used to extract a subset of an image
            (used when `singles` inherits from `MultibandBase` to
            create a new object with a slice of the original).
        deep: bool
            Whether or not this is a deep copy.
            Only used when `singles` is a `MultibandBase`.
        kwargs: dict
            Keyword arguments to pass to `_fullInitialize` to
            initialize a new instance of an inherited class.
        """
        if isinstance(singles, MultibandBase):
            if not deep:
                self._filters = singles.filters
                self._singles = singles.singles
            else:
                self._filters = [f for f in singles.filters]
                self._singles = self._copySingles(singles)
            if slices is None:
                self._array = singles._array
            else:
                try:
                    # This will fail if the first slice is a string,
                    # or list of strings
                    self._array = singles._array[slices]
                    if len(singles.filters) > 1:
                        if isinstance(slices[0], slice):
                            self._filters = singles.filters[slices[0]]
                        elif np.issubdtype(type(slices[0]), np.integer):
                            self._filters = [singles.filters[slices[0]]]
                        else:
                            self._filters = [singles.filters[f] for f in slices[0]]
                    else:
                        self._filters = singles.filters
                except IndexError:
                    # Filter using the list of strings
                    _filters = slices[0]
                    slices = list(slices)
                    slices[0] = [self.filters.index(f) for f in slices[0]]
                    self._array = singles._array[slices]
                    self._filters = _filters
            # Set the bbox size based on the slices
            x0 = singles.x0
            y0 = singles.y0
            if len(slices) > len(singles.array.shape):
                err = "Too many indices, expected {0} but received {1}"
                raise IndexError(err.format(len(singles.array.shape), len(slices)))
            if len(slices) == len(singles.array.shape):
                xslice = slices[-1]
                if isinstance(xslice, slice):
                    if xslice.start is not None:
                        x0 += xslice.start
                else:
                    raise IndexError("x index must be a slice")
            if len(slices) >= len(singles.array.shape)-1:
                if len(slices) == len(singles.array.shape):
                    yslice = slices[-2]
                else:
                    yslice = slices[-1]
                if isinstance(yslice, slice):
                    if yslice.start is not None:
                        y0 += yslice.start
                else:
                    raise IndexError("y index must be a slice")
            xy0 = Point2I(x0, y0)
            xyF = Point2I(x0+self.array.shape[-1]-1, y0+self.array.shape[-2]-1)
            self._xy0 = xy0
            self._bbox = Box2I(xy0, xyF)
            for singleObj in self.singles:
                singleObj.setXY0(xy0)

            # Create instances of the single band objects
            # whose data points to the appropriate slice
            # of self.array
            self._updateSingles(singles.imageType)
        else:
            # Extract the single band objects and filters
            if isinstance(singles, OrderedDict):
                self._filters = list(singles.keys())
                _singles = list(singles.values())
            elif singles is not None:
                self._filters = filters
                _singles = singles
            # Set the required attributes
            self._bbox = _singles[0].getBBox()
            self._xy0 = _singles[0].getXY0()
            if _singles is not None and self.filters is not None:
                self.imageType = type(_singles[0])
                self._singles = _singles
            elif singles is None and self.filters is not None:
                # Attempt to load a set of images
                self.imageType = imageType
                self._singles = []
                for f in self.filters:
                    if filterKwargs is not None:
                        for key, value in filterKwargs:
                            kwargs[key] = value[f]
                    self._singles.append(imageType(**kwargs))
            elif singles is None:
                err = """Currently an OrderedDict
                         or a list of filters with a
                         list of lsst.afw.image.Image<X> objects or
                         kwargs to open a list of new Image<X> objects
                         is required"""
                raise NotImplementedError(err)

            assert all([img.getBBox() == self.bbox for img in self.singles])
            assert all([type(img) == self.imageType for img in self.singles])
            self._array = np.array([image.array for image in self.singles])
            self._updateSingles(type(self.singles[0]))

        # Make sure that all of the parameters have been setup appropriately
        assert isinstance(self._bbox, Box2I)
        assert isinstance(self.XY0, Point2I)

    @property
    def array(self):
        """Data cube array in multiple bands
        """
        return self._array

    @property
    def shape(self):
        """Shape of the Multiband Object
        """
        return self.array.shape

    def _updateSingles(self, imageType):
        """Update the Image<X> in each band

        This method is called when a `MultibandImage` is initialized.
        """
        self.imageType = imageType
        if len(self.array.shape) == 2:
            self._singles = [self.imageType(array=self.array, xy0=self.XY0)]
        else:
            self._singles = [self.imageType(array=self.array[n], xy0=self.XY0)
                             for n in range(len(self.filters))]

    def _copySingles(self, singles):
        """Perform a deep copy of the `Image<X>` objects
        """
        self.singles = [
            self.imageType(array=s.array.copy, deep=True, xy0=s.array.XY0)
            for s in singles
        ]
