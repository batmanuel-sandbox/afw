/*
 * LSST Data Management System
 * Copyright 2014 LSST Corporation.
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

#if !defined(LSST_AFW_CAMERAGEOM_CAMERASYS_H)
#define LSST_AFW_CAMERAGEOM_CAMERASYS_H

#include <functional>
#include <string>
#include <ostream>
#include <sstream>

namespace lsst {
namespace afw {
namespace cameraGeom {

/**
 * Camera coordinate system prefix
 *
 * Used for coordinate systems that are detector-based before the detector name is known
 * (e.g. for constants such as PIXELS).
 *
 * This is Jim Bosch's clever idea for simplifying Detector.convert;
 * CameraSys is always complete and CameraSysPrefix is not.
 */
class CameraSysPrefix {
public:
    explicit CameraSysPrefix(std::string const &sysName  ///< coordinate system name
                             )
            : _sysName(sysName) {}

    ~CameraSysPrefix() = default;
    CameraSysPrefix(CameraSysPrefix const &) = default;
    CameraSysPrefix(CameraSysPrefix &&) = default;
    CameraSysPrefix &operator=(CameraSysPrefix const &) = default;
    CameraSysPrefix &operator=(CameraSysPrefix &&) = default;

    /**
     * Get coordinate system name
     */
    std::string getSysName() const { return _sysName; };

    bool operator==(CameraSysPrefix const &rhs) const { return _sysName == rhs.getSysName(); }

    bool operator!=(CameraSysPrefix const &rhs) const { return !(*this == rhs); }

    /**
     * Hash function for this object.
     *
     * @return a value that is guaranteed equal for any two equal
     *         CameraSysPrefix, and unlikely to be equal for any two unequal
     *         CameraSysPrefix.
     *
     * @note Workhorse for std::hash<CameraSysPrefix>.
     */
    size_t hash() const noexcept;

private:
    std::string _sysName;  ///< coordinate system name
};

/**
 * Camera coordinate system; used as a key in in TransformMap
 */
class CameraSys {
public:
    /**
     * Construct a CameraSys from a sysName and a detectorName
     */
    explicit CameraSys(std::string const &sysName,           ///< coordinate system name
                       std::string const &detectorName = ""  ///< detector name
                       )
            : _sysName(sysName), _detectorName(detectorName){};

    /**
     * Construct a CameraSys from a CameraSysPrefix and a detectorName
     */
    explicit CameraSys(CameraSysPrefix const &sysPrefix,     ///< coordinate system prefix
                       std::string const &detectorName = ""  ///< detector name
                       )
            : _sysName(sysPrefix.getSysName()), _detectorName(detectorName){};

    /// default constructor so SWIG can wrap a vector of pairs containing these
    CameraSys() : _sysName("?"), _detectorName(){};

    ~CameraSys() = default;
    CameraSys(CameraSys const &) = default;
    CameraSys(CameraSys &&) = default;
    CameraSys &operator=(CameraSys const &) = default;
    CameraSys &operator=(CameraSys &&) = default;

    /**
     * Get coordinate system name
     */
    std::string getSysName() const { return _sysName; };

    /**
     * Get detector name, or "" if not a detector-specific coordinate system
     */
    std::string getDetectorName() const { return _detectorName; };

    /**
     * Does this have a non-blank detector name?
     */
    bool hasDetectorName() const { return !_detectorName.empty(); }

    bool operator==(CameraSys const &rhs) const {
        return _sysName == rhs.getSysName() && _detectorName == rhs.getDetectorName();
    }

    bool operator!=(CameraSys const &rhs) const { return !(*this == rhs); }

    // less-than operator required for use in std::map
    bool operator<(CameraSys const &rhs) const {
        if (_sysName == rhs.getSysName()) {
            return _detectorName < rhs.getDetectorName();
        } else {
            return _sysName < rhs.getSysName();
        }
    }

    /**
     * Hash function for this object.
     *
     * @return a value that is guaranteed equal for any two equal CameraSys,
     *         and unlikely to be equal for any two unequal CameraSys.
     *
     * @note Workhorse for std::hash<CameraSys>.
     */
    size_t hash() const noexcept;

private:
    std::string _sysName;       ///< coordinate system name
    std::string _detectorName;  ///< detector name; "" if not a detector-specific coordinate system
};

// *** Standard camera coordinate systems ***

/**
 * Focal plane coordinates:
 * Rectilinear x, y (and z when talking about the location of a detector) on the camera focal plane (mm).
 * For z=0 choose a convenient point near the focus at x, y = 0.
 */
extern CameraSys const FOCAL_PLANE;

/**
 * Field angle coordinates:
 * Angular x,y offset from the optic axis (radians).
 */
extern CameraSys const FIELD_ANGLE;

/**
 * Nominal pixels on the detector (unbinned)
 * This ignores manufacturing imperfections, "tree ring" distortions and all other such effects.
 * It is a uniform grid of rectangular (usually square) pixels.
 *
 * This is a detector prefix; call Detector.makeCameraSys(PIXELS) to make a full coordsys.
 */
extern CameraSysPrefix const PIXELS;

/**
 * Tangent-plane pixels on the detector (unbinned)
 *
 * Converting from PIXELS to TAN_PIXELS has the effect of removing optical distortion
 * (and the distortion due to rectangular pixels)
 * with the point at the center of the detector being unaffected by the transformation.
 *
 * In detail, PIXELS->TAN_PIXELS is PIXELS->FIELD_ANGLE plus an affine transformation, such that:
 * * The x,y axes are parallel to the detector axes
 * * The dimensions are nominal pixels at the center of the focal plane
 *   (where nominal pixels size is mean of x, y pixel size).
 * * The point at the center of the detector has the same value in PIXELS and TAN_PIXELS
 *
 * This is a detector prefix; call Detector.makeCameraSys(TAN_PIXELS) to make a full coordsys.
 */
extern CameraSysPrefix const TAN_PIXELS;

/**
 * The actual pixels where the photon lands and electrons are generated (unbinned)
 * This takes into account manufacturing defects, "tree ring" distortions and other such effects.
 *
 * This is a detector prefix; call Detector.makeCameraSys(ACTUAL_PIXELS) to make a full coordsys.
 */
extern CameraSysPrefix const ACTUAL_PIXELS;

std::ostream &operator<<(std::ostream &os, CameraSysPrefix const &detSysPrefix);

std::ostream &operator<<(std::ostream &os, CameraSys const &cameraSys);
}
}
}

namespace std {
template <>
struct hash<lsst::afw::cameraGeom::CameraSysPrefix> {
    size_t operator()(lsst::afw::cameraGeom::CameraSysPrefix const &obj) const noexcept { return obj.hash(); }
};

template <>
struct hash<lsst::afw::cameraGeom::CameraSys> {
    size_t operator()(lsst::afw::cameraGeom::CameraSys const &obj) const noexcept { return obj.hash(); }
};
}

#endif
