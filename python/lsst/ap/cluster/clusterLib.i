// -*- lsst-c++ -*-
%define ap_cluster_DOCSTRING
"
Access to association pipeline clustering functions.
"
%enddef

%feature("autodoc", "1");
%module(package="lsst.ap.cluster", docstring=ap_cluster_DOCSTRING) clusterLib

// Suppress swig complaints
#pragma SWIG nowarn=314                 // print is a python keyword (--> _print)
#pragma SWIG nowarn=362                 // operator=  ignored

%{
#include "lsst/pex/policy.h"
#include "lsst/afw/geom.h"
#include "lsst/afw/detection/Source.h"
#include "lsst/ap/cluster/PT1SkyTile.h"
#include "lsst/ap/cluster/SourceCluster.h"
%}

namespace boost {
#if defined(SWIGWORDSIZE64)
    typedef long int64_t;
#else
    typedef long long int64_t;
#endif
}

%include "lsst/p_lsstSwig.i"
// %include "lsst/daf/base/persistenceMacros.i" (will need this later)

%import "lsst/pex/policy/policyLib.i"
%import "lsst/afw/detection/detectionLib.i"

%lsst_exceptions()

SWIG_SHARED_PTR(PT1SkyTile, lsst::ap::cluster::PT1SkyTile);
SWIG_SHARED_PTR(SourceClusterAttributes, lsst::ap::cluster::SourceCluster);

%typemap(out) std::vector<lsst::afw::detection::SourceSet> {
    // $1 is a SwigValueWrapper, must dereference contents to call member functions
    int len = (&($1))->size();
    swig_type_info * info = SWIGTYPE_p_std__vectorT_boost__shared_ptrT_lsst__afw__detection__Source_t_std__allocatorT_boost__shared_ptrT_lsst__afw__detection__Source_t_t_t;
    // Something like the following should work but does not: 
    // swig_type_info * info =  SWIG_TypeQuery(
    //     "std::vector<boost::shared_ptr<lsst::afw::detection::Source>, ");
    //     "std::allocator< boost::shared_ptr< lsst::afw::detection::Source > > >");
    $result = PyList_New(len);
    for (int i = 0; i < len; i++) {
        lsst::afw::detection::SourceSet * sourceSet = 
            new lsst::afw::detection::SourceSet((&($1))->operator[](i));
        PyObject * obj = SWIG_NewPointerObj(SWIG_as_voidptr(sourceSet), info, SWIG_POINTER_OWN);
        PyList_SetItem($result, i, obj);
    }
}

%import "lsst/ap/Common.h"
%include "lsst/ap/cluster/PT1SkyTile.h"
%include "lsst/ap/cluster/SourceCluster.h"

%template(SourceClusterSet) std::vector<lsst::ap::cluster::SourceClusterAttributes::Ptr>;

