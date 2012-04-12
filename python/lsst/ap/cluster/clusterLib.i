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
Access to association pipeline clustering functionality.
"
%enddef

%feature("autodoc", "1");
%module(package="lsst.ap.cluster", docstring=ap_cluster_DOCSTRING) clusterLib

// Suppress swig complaints
#pragma SWIG nowarn=314                 // print is a python keyword (--> _print)
#pragma SWIG nowarn=362                 // operator=  ignored

%{
#include "lsst/daf/base.h"
#include "lsst/pex/logging.h"
#include "lsst/afw/image.h"
#include "lsst/afw/cameraGeom.h"
#include "lsst/afw/table.h"
#include "lsst/ap/cluster/ClusteringControl.h"
#include "lsst/ap/cluster/SourceProcessingControl.h"
#include "lsst/ap/cluster/SourceCluster.h"
#include "lsst/ap/cluster/clustering.h"

#define PY_ARRAY_UNIQUE_SYMBOL LSST_AFW_TABLE_NUMPY_ARRAY_API
#include "numpy/arrayobject.h"
#include "ndarray/swig.h"
#include "ndarray/swig/eigen.h"
%}

%include "ndarray.i"
%init %{
    import_array();
%}

%declareNumPyConverters(ndarray::Array<int,1,1>);
%declareNumPyConverters(ndarray::Array<int const,1,1>);
%declareNumPyConverters(ndarray::Array<double,1>);
%declareNumPyConverters(ndarray::Array<double const,1>);
%declareNumPyConverters(Eigen::Matrix<double,2,2>);

%include "lsst/p_lsstSwig.i"

%import "lsst/daf/base/baseLib.i"
%import "lsst/afw/geom/geomLib.i"
%import "lsst/afw/coord/coordLib.i"
%import "lsst/afw/table/tableLib.i"

%include "lsst/pex/config.h"

%lsst_exceptions()

// -- Control objects --------

%shared_ptr(lsst::ap::cluster::ClusteringControl);
%shared_ptr(lsst::ap::cluster::SourceProcessingControl);
%include "lsst/ap/cluster/ClusteringControl.h"
%include "lsst/ap/cluster/SourceProcessingControl.h"


// -- SourceCluster table, record, and ID factory --------

%shared_ptr(lsst::ap::cluster::SourceClusterIdFactory);
%shared_ptr(lsst::ap::cluster::SourceClusterRecord);
%shared_ptr(lsst::ap::cluster::SourceClusterTable);

%include "lsst/ap/cluster/SourceCluster.h"

%addCastMethod(lsst::ap::cluster::SourceClusterTable, lsst::afw::table::BaseTable)
%addCastMethod(lsst::ap::cluster::SourceClusterRecord, lsst::afw::table::BaseRecord)

%template(SourceClusterColumnViewBase) lsst::afw::table::ColumnViewT<lsst::ap::cluster::SourceClusterRecord>;
%template(SourceClusterColumnView) lsst::ap::cluster::SourceClusterColumnViewT<lsst::ap::cluster::SourceClusterRecord>;

%pythondynamic;

%template(SourceClusterCatalogBase) lsst::afw::table::CatalogT<lsst::ap::cluster::SourceClusterRecord>;
%template(SourceClusterCatalog)  lsst::afw::table::SimpleCatalogT<lsst::ap::cluster::SourceClusterRecord>;
%extend lsst::afw::table::SimpleCatalogT<lsst::ap::cluster::SourceClusterRecord> {
    %pythoncode %{
        SourceClusterCatalog.Table = SourceClusterTable
        SourceClusterCatalog.Record = SourceClusterRecord
        SourceClusterCatalog.ColumnView = SourceClusterColumnView
    %}
}
%pythoncode %{
    SourceClusterRecord.Table = SourceClusterTable
    SourceClusterRecord.Catalog = SourceClusterCatalog
    SourceClusterRecord.ColumnView = SourceClusterColumnView
    SourceClusterTable.Record = SourceClusterRecord
    SourceClusterTable.Catalog = SourceClusterCatalog
    SourceClusterTable.ColumnView = SourceClusterColumnView
    SourceClusterColumnView.Record = SourceClusterRecord
    SourceClusterColumnView.Table = SourceClusterTable
    SourceClusterColumnView.Catalog = SourceClusterCatalog
%}

%pythonnondynamic;


// -- clustering and attribute computation --------

%import "lsst/ap/match/matchLib.i"
%import "lsst/ap/utils/utilsLib.i"

%typemap(out) std::vector<lsst::afw::table::SourceCatalog> const {
    // $1 is a SwigValueWrapper, must dereference contents to call member functions
    Py_ssize_t len = static_cast<Py_ssize_t>((*(&$1)).size());
    $result = PyList_New(len);
    for (Py_ssize_t i = 0; i < len; ++i) {
        lsst::afw::table::SourceCatalog * cat = 
            new lsst::afw::table::SourceCatalog((*(&$1))[i]);
        PyObject * obj = SWIG_NewPointerObj(SWIG_as_voidptr(cat),
            SWIGTYPE_p_lsst__afw__table__SimpleCatalogTT_lsst__afw__table__SourceRecord_t,
            SWIG_POINTER_OWN);
        PyList_SetItem($result, i, obj);
    }
}


%typemap(out) std::pair<boost::shared_ptr<lsst::afw::table::SourceTable>, lsst::afw::table::SchemaMapper> const {
    $result = PyTuple_New(2);
    boost::shared_ptr<lsst::afw::table::SourceTable> *table =
        new boost::shared_ptr<lsst::afw::table::SourceTable>((*(&$1)).first);
    PyObject * obj1 = SWIG_NewPointerObj(SWIG_as_voidptr(table),
        SWIGTYPE_p_boost__shared_ptrT_lsst__afw__table__SourceTable_t, SWIG_POINTER_OWN);
    PyTuple_SetItem($result, 0, obj1);
    lsst::afw::table::SchemaMapper * mapper =
        new lsst::afw::table::SchemaMapper((*(&$1)).second);
    PyObject * obj2 = SWIG_NewPointerObj(SWIG_as_voidptr(mapper),
        SWIGTYPE_p_lsst__afw__table__SchemaMapper, SWIG_POINTER_OWN);
    PyTuple_SetItem($result, 1, obj2);
}

%include "lsst/ap/cluster/clustering.h"

