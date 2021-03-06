// -*- lsst-c++ -*-

/*
 * LSST Data Management System
 * Copyright 2008, 2009, 2010, 2012 LSST Corporation.
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

%define ap_match_DOCSTRING
"
Access to association pipeline matching functionality.
"
%enddef

%feature("autodoc", "1");
%module(package="lsst.ap.match", docstring=ap_match_DOCSTRING) matchLib

// Suppress swig complaints
#pragma SWIG nowarn=314                 // print is a python keyword (--> _print)
#pragma SWIG nowarn=362                 // operator=  ignored

%{
#include "lsst/daf/base.h"
#include "lsst/pex/logging.h"
#include "lsst/afw/table.h"
#include "lsst/afw/image.h"
#include "lsst/afw/detection.h"
#include "lsst/afw/cameraGeom.h"
#include "lsst/ap/match/ExposureInfo.h"
#include "lsst/ap/match/ReferenceMatch.h"

#define PY_ARRAY_UNIQUE_SYMBOL LSST_AP_MATCH_NUMPY_ARRAY_API
#include "numpy/arrayobject.h"
#include "ndarray/swig.h"
#include "ndarray/swig/eigen.h"
%}

%init %{
    import_array();
%}

%declareNumPyConverters(Eigen::Matrix<double,2,1,Eigen::DontAlign>);
%declareNumPyConverters(Eigen::Matrix<double,3,1,Eigen::DontAlign>);

%include "lsst/p_lsstSwig.i"

%import "lsst/daf/base/baseLib.i"
%import "lsst/afw/geom/geomLib.i"
%import "lsst/afw/image/imageLib.i"
%import "lsst/ap/utils/utilsLib.i"

%include "ndarray.i"
%include "lsst/pex/config.h"

%lsst_exceptions()

%shared_ptr(lsst::ap::match::BBox);
%shared_ptr(lsst::ap::match::CatalogControl);
%shared_ptr(lsst::ap::match::ExposureInfo);
%shared_ptr(lsst::ap::match::ExposureInfoMap);

%include "lsst/ap/match/BBox.h"
%include "lsst/ap/match/CatalogControl.h"
%include "lsst/ap/match/ExposureInfo.h"
%include "lsst/ap/match/ReferenceMatch.h"

%template(ExposureInfoVector) std::vector<lsst::ap::match::ExposureInfo::Ptr>;

