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
#include "lsst/daf/persistence.h"
#include "lsst/ap/Results.h"
//#include "lsst/ap/io/ResultFormatters.h"
#include "lsst/ap/Stages.h"
#include "lsst/ap/Utils.h"
#include <sstream>
%}

%inline %{
namespace boost { namespace filesystem {} }
%}

%include "lsst/p_lsstSwig.i"
%include "lsst/daf/base/persistenceMacros.i"

%import "lsst/daf/base/baseLib.i"
%import "lsst/pex/policy/policyLib.i"

%include <stdint.i>
%include <typemaps.i>
%include <std_pair.i>
%include <std_vector.i>

%import  "lsst/ap/Common.h"

%rename(IdVec)        lsst::ap::IdVector;
%rename(IdPairVec)    lsst::ap::IdPairVector;
%rename(MatchPairVec) lsst::ap::MatchPairVector;

%include "lsst/ap/Results.h"


%define %lsst_idpair(Type)
    %template(IdPair) std::pair<Type, Type>;

    %extend std::pair<Type, Type> {
        std::string toString() {
            std::ostringstream os;
            os << "IdPair (" << $self->first << ", " << $self->second << ")";
            return os.str();
        }
    };

    %pythoncode %{
    IdPair.__str__ = IdPair.toString
    %}
%enddef

#if defined(SWIGWORDSIZE64)
    %lsst_idpair(long);
#else
    %lsst_idpair(long long);
#endif

%extend lsst::ap::MatchPair {
    std::string toString() {
        std::ostringstream os;
        os << "MatchPair (" << $self->getFirst() << ", " << $self->getSecond() << ", " <<
              $self->getDistance() << ")";
        return os.str();
    }
};

%pythoncode %{
MatchPair.__str__ = MatchPair.toString
%}



// Export instantiations of boost::shared_ptr for persistable data vectors
%lsst_persistable(lsst::ap::IdVector);
%lsst_persistable(lsst::ap::IdPairVector);
%lsst_persistable(lsst::ap::MatchPairVector);

%include "lsst/ap/RectangularRegion.h"
%include "lsst/ap/CircularRegion.h"
%include "lsst/ap/SpatialUtil.h"
%include "lsst/ap/Stages.h"

