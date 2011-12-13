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
 
%define ap_DOCSTRING
"
Access to association pipeline persistable result objects and implementation methods.
"
%enddef

%feature("autodoc", "1");
%module(package="lsst.ap", docstring=ap_DOCSTRING) apLib

// Suppress swig complaints
#pragma SWIG nowarn=314                 // print is a python keyword (--> _print)
#pragma SWIG nowarn=362                 // operator=  ignored

%{
#include "lsst/daf/base.h"
#include "lsst/pex/policy.h"
#include "lsst/daf/data.h"
#include "lsst/afw.h"
#include "lsst/ap/Results.h"
#include "lsst/ap/io/ResultFormatters.h"
#include "lsst/ap/Stages.h"
#include "lsst/ap/Utils.h"
#include <sstream>
%}

%include "std_pair.i"
%include "lsst/p_lsstSwig.i"
%include "lsst/daf/base/persistenceMacros.i"

%import "lsst/daf/base/baseLib.i"
%import "lsst/pex/policy/policyLib.i"
%import "lsst/afw/detection/detectionLib.i"
%import "lsst/mops/mopsLib.i"

%lsst_exceptions()

%shared_ptr(lsst::ap::PersistableIdVector);
%shared_ptr(lsst::ap::PersistableIdPairVector);
%shared_ptr(lsst::ap::PersistableMatchPairVector);

%include "lsst/ap/Common.h"
%include "lsst/ap/Results.h"

%template(IdPair) std::pair<boost::int64_t, boost::int64_t>;
%extend std::pair<boost::int64_t, boost::int64_t> {
    std::string toString() {
        std::ostringstream os;
        os << "IdPair (" << $self->first << ", " << $self->second << ")";
        return os.str();
    }
};

%extend lsst::ap::MatchPair {
    std::string toString() {
        std::ostringstream os;
        os << "MatchPair (" << $self->getFirst() << ", " << $self->getSecond() << ", " <<
              $self->getDistance() << ")";
        return os.str();
    }
};

%pythoncode %{
IdPair.__str__ = IdPair.toString
MatchPair.__str__ = MatchPair.toString
%}

%template(IdPairVec) std::vector<std::pair<boost::int64_t, boost::int64_t> >;
%template(MatchPairVec) std::vector<lsst::ap::MatchPair>;

%lsst_persistable(lsst::ap::PersistableIdVector);
%lsst_persistable(lsst::ap::PersistableIdPairVector);
%lsst_persistable(lsst::ap::PersistableMatchPairVector);

%include "lsst/ap/CircularRegion.h"
%include "lsst/ap/RectangularRegion.h"
%include "lsst/ap/SpatialUtil.h"

%shared_ptr(VisitProcessingContext, lsst::ap::VisitProcessingContext);

%include "lsst/ap/Stages.h"

