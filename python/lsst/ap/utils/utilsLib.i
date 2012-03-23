// -*- lsst-c++ -*-

/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
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
 
%define ap_utils_DOCSTRING
"
Access to association pipeline utilities.
"
%enddef

%feature("autodoc", "1");
%module(package="lsst.ap.utils", docstring=ap_utils_DOCSTRING) utilsLib

// Suppress swig complaints
#pragma SWIG nowarn=314                 // print is a python keyword (--> _print)
#pragma SWIG nowarn=362                 // operator=  ignored

%{
#include "lsst/tr1/unordered_map.h"
#include "lsst/daf/base.h"
#include "lsst/pex/logging.h"
#include "lsst/pex/policy.h"
#include "lsst/afw/geom.h"
#include "lsst/afw/table.h"
#include "lsst/afw/image.h"
#include "lsst/afw/cameraGeom.h"
#include "lsst/ap/utils/CsvControl.h"
#include "lsst/ap/utils/ImageUtils.h"
#include "lsst/ap/utils/PT1SkyTile.h"
%}

%include "lsst/p_lsstSwig.i"

%import "lsst/daf/base/baseLib.i"
%import "lsst/afw/table/tableLib.i"
%import "lsst/afw/image/imageLib.i"

%include "lsst/pex/config.h"

%lsst_exceptions()

%shared_ptr(lsst::ap::utils::CsvControl);
%shared_ptr(lsst::ap::utils::PT1SkyTile);

%import "lsst/ap/Common.h"
%include "lsst/ap/utils/ImageUtils.h"
%include "lsst/ap/utils/PT1SkyTile.h"
%include "lsst/ap/utils/CsvControl.h"

