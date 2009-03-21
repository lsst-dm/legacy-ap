// -*- lsst-c++ -*-
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

namespace boost {
#if defined(SWIGWORDSIZE64)
    typedef long int64_t;
#else
    typedef long long int64_t;
#endif
}

%include "std_pair.i"
%include "lsst/p_lsstSwig.i"
%include "lsst/daf/base/persistenceMacros.i"

%import "lsst/daf/base/baseLib.i"
%import "lsst/pex/policy/policyLib.i"
%import "lsst/afw/detection/detectionLib.i"
%import "lsst/mops/mopsLib.i"

%lsst_exceptions()

SWIG_SHARED_PTR(PersistableIdVec, lsst::ap::PersistableIdVector);
SWIG_SHARED_PTR(PersistableIdPairVec, lsst::ap::PersistableIdPairVector);
SWIG_SHARED_PTR(PersistableMatchPairVec, lsst::ap::PersistableMatchPairVector);

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

SWIG_SHARED_PTR(VisitProcessingContext, lsst::ap::VisitProcessingContext);

%include "lsst/ap/Stages.h"

