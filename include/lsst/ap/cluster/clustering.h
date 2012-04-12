// -*- lsst-c++ -*-

/* 
 * LSST Data Management System
 * Copyright 2012 LSST Corporation.
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

/** @file
  * @brief  Functions for clustering sources.
  * @author Serge Monkewitz
  */
#ifndef LSST_AP_CLUSTER_CLUSTERING_H
#define LSST_AP_CLUSTER_CLUSTERING_H

#include <utility>
#include <vector>

#include "../match/ExposureInfo.h"
#include "../utils/PT1SkyTile.h"
#include "ClusteringControl.h"
#include "SourceProcessingControl.h"
#include "SourceCluster.h"


namespace lsst { namespace ap { namespace cluster {

/** Creates an output source table based on a prototypical input source table.
  * The schema and slot mappings from the input table are used as is, and exposure
  * related fields, required for database ingestion / clustering, are added.
  *
  * In particular, the following additions are performed:
  *
  *     - "<coord>.err"          (Cov<Point<F8>>) : Source sky-coordinate error
  *     - "<exposure>.id"        (I8)             : ID of exposure source was measured on.
  *     - "<exposure>.filter.id" (I4)             : ID of filter for exposure.
  *     - "<exposure>.time"      (F4)             : Exposure time (s).
  *     - "<exposure>.time.mid"  (F8)             : Middle of exposure time (MJD).
  *     - "<cluster>.id"         (I8)             : ID of cluster containing source.
  *     - "<cluster>.coord"      (Coord)          : ICRS coordinates of cluster containing source.
  *
  * "<exposure>" and "<cluster>" are field name prefixes obtained from SourceProcessingControl.
  * If either prefix is empty, the corresponding fields will not be added. Note however that
  * the exposure related fields are necessary for attribute computation and database ingest
  * into the LSST schema. The cluster related fields are required only for database ingest.
  *
  * The "<coord>" prefix is not configurable - it corresponds to the name of the Coord slot
  * in the prototypical input table.
  *
  * @return An output SourceTable and a SchemaMapper that can be used to copy common
  *         field values between input and output records.
  */
std::pair<boost::shared_ptr<lsst::afw::table::SourceTable>,
          lsst::afw::table::SchemaMapper> const
makeOutputSourceTable(
    ::lsst::afw::table::SourceTable const & prototype, /// @param[in] Provides prototypical schema and slot mappings.
    SourceProcessingControl const & control            /// @param[in] Source processing parameters.
);

/** Process input sourcees, distributing them to one of 3 output catalogs.
  *     - \a sources:         sources suitable for spatial clustering.
  *     - \a badSources:      sources ignored for the purposes of clustering;
  *                           identified by one or more flags.
  *     - \a invalidSources:  sources which do not have sky-coordinates or centroids.
  *
  * Sources not lying in the given sky-tile are discarded.
  *
  * Existing fields in the input sources are not touched in any way. However:
  *     - source footprints are discarded
  *     - source sky-coordinate errors are computed
  *     - exposure information (ID, filter ID, etc...) is added to each source.
  *
  * The output catalogs will typically have been constructed from tables obtained
  * from makeOutputSourceTable - their schemas and slot mappings must all be 
  * identical. The input table slots must match output table slots, and the
  * input schema must be fully contained in the output schema.
  */
void processSources(
    lsst::afw::table::SourceCatalog const & expSources, /// @param[in] Single exposure sources to process.
    lsst::ap::match::ExposureInfo const & expInfo,     /// @param[in] Exposure information.
    lsst::ap::utils::PT1SkyTile const & skyTile,       /// @param[in] Sky-tile being processed.
    SourceProcessingControl const & control,           /// @param[in] Source processing parameters.
    lsst::afw::table::SchemaMapper const & mapper,     /// @param[in] Maps between input and output source records.
    lsst::afw::table::SourceCatalog & sources,         /// @param[inout] Catalog for sources that will be clustered.
    lsst::afw::table::SourceCatalog & badSources,      /// @param[inout] Catalog for sources with bad measurement flags.
    lsst::afw::table::SourceCatalog & invalidSources   /// @param[inout] Catalog for sources with invalid measurements.
);

/** Spatially cluster sources using the OPTICS algorithm.
  */
std::vector<lsst::afw::table::SourceCatalog> const cluster(
    lsst::afw::table::SourceCatalog const & sources, ///< @param[in] Sources to cluster.
    ClusteringControl const & control                ///< @param[in] Clustering parameters.
);

}}} // namespace lsst::ap::cluster

#endif // LSST_AP_CLUSTER_CLUSTERING_H

