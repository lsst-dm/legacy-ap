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
// #include "lsst/tr1/unordered_map.h"
#include "lsst/daf/base.h"
#include "lsst/pex/policy.h"
#include "lsst/ap/match/ReferenceMatch.h"
%}

namespace boost {
#if defined(SWIGWORDSIZE64)
    typedef long int64_t;
#else
    typedef long long int64_t;
#endif
}

%include "lsst/p_lsstSwig.i"

%import "lsst/daf/base/baseLib.i"
%import "lsst/pex/policy/policyLib.i"

%lsst_exceptions()

%import "lsst/ap/Common.h"
%include "lsst/ap/match/ReferenceMatch.h"

