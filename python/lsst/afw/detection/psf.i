// -*- lsst-C++ -*-

%{
#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/detection/DoubleGaussianPsf.h"
%}

%shared_ptr(lsst::afw::detection::Psf);
%shared_ptr(lsst::afw::detection::KernelPsf);
%shared_ptr(lsst::afw::detection::DoubleGaussianPsf);

%ignore PsfFactoryBase;
%ignore lsst::afw::detection::Psf::writeToRecords;
%ignore lsst::afw::detection::Psf::readFromRecords;
%include "lsst/afw/detection/Psf.h"
%include "lsst/afw/detection/DoubleGaussianPsf.h"

%lsst_persistable(lsst::afw::detection::Psf);
%lsst_persistable(lsst::afw::detection::DoubleGaussianPsf);
