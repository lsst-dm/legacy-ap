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
#include "lsst/afw/image/TanWcs.h"
#include "lsst/afw/table.h"
#include "lsst/ap/cluster/ClusteringControl.h"
#include "lsst/ap/cluster/SourceCluster.h"

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


// -- ClusteringControl --------

%shared_ptr(lsst::ap::cluster::ClusteringControl);
%include "lsst/ap/cluster/ClusteringControl.h"


// -- SourceCluster table, record, and ID factory --------

%shared_ptr(lsst::ap::cluster::SourceClusterIdFactory);
%shared_ptr(lsst::ap::cluster::SourceClusterRecord);
%shared_ptr(lsst::ap::cluster::SourceClusterTable);

%include "lsst/ap/cluster/SourceCluster.h"

%addCastMethod(lsst::ap::cluster::SourceClusterTable, lsst::afw::table::BaseTable)
%addCastMethod(lsst::ap::cluster::SourceClusterRecord, lsst::afw::table::BaseRecord)

%template(SourceClusterColumnViewBase) lsst::afw::table::ColumnViewT<lsst::ap::cluster::SourceClusterRecord>;
%template(SourceClusterColumnView) lsst::ap::cluster::SourceClusterColumnViewT<lsst::ap::cluster::SourceClusterRecord>;


%pythondynamic;  // We want to add attributes in Python for the classes wrapped here.

namespace lsst { namespace ap { namespace cluster {
    %template (SourceClusterCatalogBase) lsst::afw::table::CatalogT<SourceClusterRecord>;
    %declareCatalog(lsst::afw::table::SimpleCatalogT, SourceCluster);

    // The equivalent is performed by %declareCatalog(), but for some reason has no effect
    %pythoncode %{
        SourceClusterCatalog.Table = SourceClusterTable
        SourceClusterCatalog.Record = SourceClusterRecord
        SourceClusterCatalog.ColumnView = SourceClusterColumnView
    %}

}}}

%pythonnondynamic; // Re-enable attribute restriction

