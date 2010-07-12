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
#include "lsst/ap/cluster/Utils.h"
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
%import "lsst/afw/image/imageLib.i"

%lsst_exceptions()

SWIG_SHARED_PTR(PT1SkyTile, lsst::ap::cluster::PT1SkyTile);
SWIG_SHARED_PTR(PerFilterSourceClusterAttributes, lsst::ap::cluster::PerFilterSourceClusterAttributes);
SWIG_SHARED_PTR(SourceClusterAttributes, lsst::ap::cluster::SourceClusterAttributes);
SWIG_SHARED_PTR(PersistableSourceClusterVector, lsst::ap::cluster::PersistableSourceClusterVector);

%ignore lsst::ap::cluster::NullOr<float>;
%ignore lsst::ap::cluster::SourceClusterAttributes::PerFilterAttributeMap;

%typemap(out) lsst::ap::cluster::NullOr<float> const & {
    if (($1)->isNull()) {
        Py_INCREF(Py_None);
        $result = Py_None;
    } else {
        $result = PyFloat_FromDouble(static_cast<double>(*($1)));
    }
}

%typemap(in) lsst::ap::cluster::NullOr<float> const & (lsst::ap::cluster::NullOr<float> temp) {
    if (PyFloat_CheckExact($input)) {
        temp = static_cast<float>(PyFloat_AsDouble($input));
    } else if ($input == Py_None) {
        temp.setNull();
    } else {
        SWIG_exception_fail(SWIG_TypeError, "failed to convert Python input to a lsst::ap::cluster::NullOr<float>");
    }
    $1 = &temp;
}

%typemap(typecheck, precedence=SWIG_TYPECHECK_FLOAT) lsst::ap::cluster::NullOr<float> const & {
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

%typemap(out) lsst::ap::cluster::SourceClusterAttributes::PerFilterAttributesMap const & {
    $result = PyDict_New();
    swig_type_info * info = SWIGTYPE_p_boost__shared_ptrT_lsst__ap__cluster__PerFilterSourceClusterAttributes_t;
    typedef lsst::ap::cluster::SourceClusterAttributes::PerFilterAttributesMap _PFAMap;
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
%include "lsst/ap/cluster/Utils.h"

%lsst_persistable(lsst::ap::cluster::PersistableSourceClusterVector);
%template(SourceClusterVector) std::vector<lsst::ap::cluster::SourceClusterAttributes::Ptr>;

