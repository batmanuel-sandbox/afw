// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2017 LSST Corporation.
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

#ifndef LSST_AFW_GEOM_FRAMESETUTILS_H
#define LSST_AFW_GEOM_FRAMESETUTILS_H

#include <memory>
#include <vector>

#include "astshim.h"
#include "ndarray.h"

#include "lsst/afw/coord/Coord.h"
#include "lsst/afw/geom/Endpoint.h"
#include "lsst/afw/geom/Transform.h"
#include "lsst/daf/base/PropertyList.h"
#include "lsst/daf/base/PropertySet.h"

namespace lsst {
namespace afw {
namespace geom {
namespace detail {

/**
Read a FITS convention WCS FrameSet from FITS metadata

The resulting FrameSet may be any kind of WCS supported by FITS;
if it is a celestial WCS then 1,1 will be the lower left corner of the image
(the FITS convention, not the LSST convention).

This routine replaces RADECSYS with RADESYS if the former is present and the latter is not,
since that is a common misspelling in FITS headers.

The returned FrameSet will have an IWC (intermediate world coordinate system) frame.

@param[in,out] metadata  FITS header cards
@param[in] strip  If true: strip items from `metadata` used to create the WCS,
   such as RADESYS, EQUINOX, CTYPE12, CRPIX12, CRVAL12, etc.
   Always keep keywords that might be wanted for other purpposes, including NAXIS12
   and date-related keywords such as "DATE-OBS" and "TIMESYS" (but not "EQUINOX").

@throws pex::exceptions::TypeError if the metadata does not contain a FITS-WCS
*/
std::shared_ptr<ast::FrameSet> readFitsWcs(daf::base::PropertySet& metadata, bool strip = true);

/**
Read an LSST celestial WCS FrameDict from a FITS header.

Calls getImageXY0FromMetadata to determine image XY0.

Saves CRVAL by setting the output SkyFrame's SkyRef to CRVAL and SkyRefIs="Ignore"
(so SkyRef is not treated as an offset).

@warning Does not compensate for the offset between a subimage and its parent image;
callers must do that manually, e.g. by calling SkyWcs::copyAtShiftedPosition.

@warning the general SkyWcs generated by LSST software cannot be exactly represented using
standard WCS FITS cards, and so this function cannot read those. This function is intended
for two purposes:
- Read the standard FITS WCS found in raw data and other images from non-LSST observatories
    and convert it to the LSST pixel convention.
- Read the approximate FITS WCS that LSST writes to FITS images (for use by non-LSST code).

The frames of the returned WCS will be as follows:
- "PIXELS" (base frame): pixel frame with 0,0 the lower left corner of the image (LSST convention)
- "IWC": FITS WCS intermediate world coordinate system
- "SKY" (current frame): an ast::SkyFrame with domain "SKY": ICRS RA, Dec

All frames are instances of ast::Frame except the SKY frame. All have 2 axes.

@param[in,out] metadata  FITS header cards
@param[in] strip  If true: strip items from `metadata` used to create the WCS,
   such as RADESYS, EQUINOX, CTYPE12, CRPIX12, CRVAL12, etc.
   Always keep keywords that might be wanted for other purpposes, including NAXIS12
   and date-related keywords such as "DATE-OBS" and "TIMESYS" (but not "EQUINOX").

@throws lsst::pex::exceptions::TypeError if the metadata does not describe a celestial WCS.
*/
std::shared_ptr<ast::FrameDict> readLsstSkyWcs(daf::base::PropertySet& metadata, bool strip = true);

/**
Copy values from an AST FitsChan into a PropertyList

@warning COMMENT and HISTORY cards are treated as string values

@throws lsst::pex::exceptions::TypeError if `fitsChan` contains cards whose
    type is not supported by `PropertyList`: complex numbers, or cards with no value
*/
std::shared_ptr<daf::base::PropertyList> getPropertyListFromFitsChan(ast::FitsChan& fitsChan);

}  // namespace detail
}  // namespace geom
}  // namespace afw
}  // namespace lsst

#endif
