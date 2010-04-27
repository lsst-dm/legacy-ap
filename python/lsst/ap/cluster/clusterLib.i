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
#include "lsst/tr1/unordered_map.h"
#include "lsst/daf/base.h"
#include "lsst/pex/policy.h"
#include "lsst/afw/geom.h"
#include "lsst/afw/image.h"
#include "lsst/afw/detection.h"
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
%include "lsst/daf/base/persistenceMacros.i"

%import "lsst/daf/base/baseLib.i"
%import "lsst/pex/policy/policyLib.i"
%import "lsst/afw/detection/detectionLib.i"

%lsst_exceptions()

SWIG_SHARED_PTR(PT1SkyTile, lsst::ap::cluster::PT1SkyTile);
SWIG_SHARED_PTR(PerFilterSourceClusterAttributes, lsst::ap::cluster::PerFilterSourceClusterAttributes);
SWIG_SHARED_PTR(SourceClusterAttributes, lsst::ap::cluster::SourceClusterAttributes);

// Forward declaration for Nullable template
namespace lsst { namespace ap { namespace cluster {
    template <typename ScalarT> class Nullable;
}}}

%ignore lsst::ap::cluster::Nullable<float>;
%ignore lsst::ap::cluster::SourceClusterAttributes::PerFilterAttributeMap;

%typemap(out) lsst::ap::cluster::Nullable<float> const {
    if ((&$1)->isNull()) {
        $result = Py_None;
    } else {
        $result = PyFloat_FromDouble(static_cast<double>(*(&$1)));
    }
}

%typemap(in) lsst::ap::cluster::Nullable<float> const (lsst::ap::cluster::Nullable<float> temp) {
    if (PyFloat_CheckExact($input)) {
        temp = static_cast<float>(PyFloat_AsDouble($input));
    } else if ($input == Py_None) {
        temp.setNull();
    } else {
        SWIG_exception_fail(SWIG_TypeError, "failed to convert Python input to a lsst::ap::cluster::Nullable<float>");
    }
    $1 = &temp;
}

%typemap(typecheck, precedence=SWIG_TYPECHECK_FLOAT) lsst::ap::cluster::Nullable<float> {
    $1 = (PyFloat_CheckExact($input) || $input == Py_None) ? 1 : 0; 
}

%typemap(out) std::vector<int> {
    int len = ($1).size();
    $result = PyList_New(len);
    for (int i = 0; i < len; ++i) {
        PyList_SetItem($result, i, PyInt_FromLong(($1)[i]));
    }
}

%typemap(out) std::vector<lsst::afw::detection::SourceSet> const {
    // $1 is a SwigValueWrapper, must dereference contents to call member functions
    int len = (*(&$1)).size();
    swig_type_info * info = SWIGTYPE_p_std__vectorT_boost__shared_ptrT_lsst__afw__detection__Source_t_std__allocatorT_boost__shared_ptrT_lsst__afw__detection__Source_t_t_t;
    // Something like the following should work but does not: 
    // swig_type_info * info =  SWIG_TypeQuery(
    //     "std::vector<boost::shared_ptr<lsst::afw::detection::Source>, ");
    //     "std::allocator< boost::shared_ptr< lsst::afw::detection::Source > > >");
    $result = PyList_New(len);
    for (int i = 0; i < len; i++) {
        lsst::afw::detection::SourceSet * sourceSet = 
            new lsst::afw::detection::SourceSet((*(&$1))[i]);
        PyObject * obj = SWIG_NewPointerObj(SWIG_as_voidptr(sourceSet), info, SWIG_POINTER_OWN);
        PyList_SetItem($result, i, obj);
    }
}

%typemap(out) lsst::ap::cluster::SourceClusterAttributes::PerFilterAttributeMap const & {
    $result = PyDict_New();
    swig_type_info * info = SWIGTYPE_p_boost__shared_ptrT_lsst__ap__cluster__PerFilterSourceClusterAttributes_t;
    typedef lsst::ap::cluster::SourceClusterAttributes::PerFilterAttributeMap _PFAMap;
    typedef lsst::ap::cluster::PerFilterSourceClusterAttributes _PFA;
    for (_PFAMap::const_iterator i = ($1)->begin(), e = ($1)->end(); i != e; ++i) {
        _PFA::Ptr * pfa = new _PFA::Ptr(new _PFA(i->second));
        PyObject * obj = SWIG_NewPointerObj(SWIG_as_voidptr(pfa), info, SWIG_POINTER_OWN);
        PyDict_SetItem($result, PyInt_FromLong(i->first), obj);
    }
}

%import "lsst/ap/Common.h"
%include "lsst/ap/cluster/PT1SkyTile.h"
%include "lsst/ap/cluster/SourceCluster.h"

%lsst_persistable(lsst::ap::cluster::SourceClusterAttributes);
%template(SourceClusterAttributesSet) std::vector<lsst::ap::cluster::SourceClusterAttributes::Ptr>;

