/*
 * LSST Data Management System
 * Copyright 2008-2016  AURA/LSST.
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
 * see <https://www.lsstcorp.org/LegalNotices/>.
 */

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "numpy/arrayobject.h"
#include "ndarray/pybind11.h"

#include "lsst/afw/geom/SipApproximation.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst { namespace afw { namespace geom { namespace {

using PySipApproximation = py::class_<SipApproximation, std::shared_ptr<SipApproximation>>;

PYBIND11_PLUGIN(sipApproximation) {
    py::module mod("sipApproximation");

    py::module::import("lsst.afw.geom.coordinates");
    py::module::import("lsst.afw.geom.box");
    py::module::import("lsst.afw.geom.linearTransform");
    py::module::import("lsst.afw.geom.transform");

    // Need to import numpy for ndarray conversions
    if (_import_array() < 0) {
        PyErr_SetString(PyExc_ImportError, "numpy.core.multiarray failed to import");
        return nullptr;
    }


    PySipApproximation cls(mod, "SipApproximation");

    cls.def(
        py::init<
            std::shared_ptr<TransformPoint2ToPoint2>,
            Point2D const &,
            LinearTransform const &,
            Box2D const &,
            Extent2I const &,
            int,
            bool,
            double
        >(),
        "pixelToIwc"_a, "crpix"_a, "cd"_a, "bbox"_a, "gridDimensions"_a, "order"_a,
        "useInverse"_a=true, "svdThreshold"_a=-1
    );

    cls.def(
        py::init<
            std::shared_ptr<TransformPoint2ToPoint2>,
            Point2D const &,
            LinearTransform const &,
            Box2D const &,
            Extent2I const &,
            ndarray::Array<double const, 2> const &,
            ndarray::Array<double const, 2> const &,
            ndarray::Array<double const, 2> const &,
            ndarray::Array<double const, 2> const &,
            bool
        >(),
        "pixelToIwc"_a, "crpix"_a, "cd"_a, "bbox"_a, "gridDimensions"_a,
        "a"_a, "b"_a, "ap"_a, "bp"_a, "useInverse"_a=true
    );

    using ScalarTransform = Point2D (SipApproximation::*)(Point2D const &) const;
    using VectorTransform = std::vector<Point2D> (SipApproximation::*)(std::vector<Point2D> const &) const;

    cls.def("getOrder", &SipApproximation::getOrder);
    cls.def("getA", &SipApproximation::getA, "p"_a, "q"_a);
    cls.def("getB", &SipApproximation::getB, "p"_a, "q"_a);
    cls.def("getAP", &SipApproximation::getAP, "p"_a, "q"_a);
    cls.def("getBP", &SipApproximation::getBP, "p"_a, "q"_a);
    cls.def("applyForward", (ScalarTransform)&SipApproximation::applyForward);
    cls.def("applyForward", (VectorTransform)&SipApproximation::applyForward);
    cls.def("applyInverse", (ScalarTransform)&SipApproximation::applyInverse);
    cls.def("applyInverse", (VectorTransform)&SipApproximation::applyInverse);
    cls.def("getGridStep", &SipApproximation::getGridStep);
    cls.def("getGridDimensions", &SipApproximation::getGridDimensions);
    cls.def("updateGrid", &SipApproximation::updateGrid, "dimensions"_a);
    cls.def("refineGrid", &SipApproximation::refineGrid, "factor"_a=2);
    cls.def("fit", &SipApproximation::fit, "order"_a, "svdThreshold"_a=-1);
    cls.def("computeMaxDeviation", &SipApproximation::computeMaxDeviation);

    return mod.ptr();
}

}}}}  // namespace lsst::afw::<anonymous>
