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

/** Create an output source table based on a prototypical input source table.
  * The schema and slot mappings from the input table are used as is, and exposure
  * related fields, required for database ingestion / clustering, are added.
  *
  * In particular, the following additions are performed:
  *
  *     - "coord.err"             (Cov<Point<F8>>) : Source sky-coordinate error
  *     - "<exposure>.id"         (I8)             : ID of exposure source was measured on.
  *     - "<exposure>.filter.id"  (I4)             : ID of filter for exposure.
  *     - "<exposure>.time"       (F4)             : Exposure time (s).
  *     - "<exposure>.time.mid"   (F8)             : Middle of exposure time (MJD).
  *     - "<cluster>.id"          (I8)             : ID of cluster containing source.
  *     - "<cluster>.coord"       (Coord)          : ICRS coordinates of cluster containing source.
  *
  * "<exposure>" and "<cluster>" are field name prefixes obtained from SourceProcessingControl.
  * If either prefix is empty, the corresponding fields will not be added. Note however that
  * the exposure related fields are necessary for attribute computation and database ingest
  * into the LSST schema. The cluster related fields are required only for database ingest.
  *
  * @param[in] prototype Prototypical schema/slot mappings.
  * @param[in] control   Source processing parameters.
  *
  * @return An output SourceTable and a SchemaMapper that can be used to copy common
  *         field values between input and output records.
  */
std::pair<boost::shared_ptr<lsst::afw::table::SourceTable>,
          lsst::afw::table::SchemaMapper> const
makeOutputSourceTable(
    lsst::afw::table::SourceTable const & prototype,
    SourceProcessingControl const & control
);

/** Create a source cluster table based on a prototypical input source table.
  * The names of flux/shape fields from SourceProcessingControl which correspond
  * to fields in the prototype are used to create flux/shape fields for every
  * filter that has been defined via the lsst::afw::image::Filter API. Slots
  * from the prototype are used to setup corresponding filter specific slots
  * in the source cluster table.
  *
  * Note that filters are typically defined by the data butler's mapper class
  * (see e.g. the obs_lsstSim package). The names of flux/shape fields in the
  * cluster schema are prefixed with filter names - for example, if the PSF
  * flux slot in the prototype is mapped to a field named "flux.naive" and
  * filters "u" and "g" are defined, then the cluster schema will contain fields
  * "u.flux.naive" and "g.flux.naive", with filter specific PSF flux slots
  * mapped accordingly.
  * 
  * The following additional fields will be created (names not configurable
  * at the moment):
  *
  *     - "id"                  (I8)             : cluster ID (from minimal schema)
  *     - "coord"               (Coord)          : mean sky-coordinates (from minimal schema)
  *     - "coord.err"           (Cov<Point<F8>>) : uncertainty of coord
  *     - "obs.count"           (I4)             : number of sources in cluster
  *     - "flag.noise"          (Flag)           : was cluster created from a single noise source?
  *     - "<filter>.obs.count"  (I4)             : numer of sources in a specific filter
  *
  * If the input table has valid centroid and centroid error slot keys,
  * then the following fields are added:
  *
  *     - "coord.weightedmean"       (Coord)          : inverse variance weighted mean sky-coordinates
  *     - "coord.weightedmean.err"   (Cov<Point<F8>>) : covariance matrix for coord.weighted
  *     - "coord.weightedmean.count" (I4)             : number of samples included in "coord.weighted"
  *
  * Finally, if the "<exposure>.time.mid" field exists in the input table, then:
  *
  *     - "obs.time.min"          (F8)             : earliest observation time of sources in cluster
  *     - "obs.time.mean"         (F8)             : mean observation time of sources in cluster
  *     - "obs.time.max"          (F8)             : latest observation time of sources in cluster
  *     - "<filter>.obs.time.min" (F8)             : earliest observation time in a specific filter
  *     - "<filter>.obs.time.max" (F8)             : latest observation time in a specific filter
  *
  * are added (where the "<exposure>" prefix is obtained from SourceProcessingControl).
  *
  * @param[in] prototype Prototypical schema/slot mappings.
  * @param[in] idFactory ID generator.
  * @param[in] control   Source processing parameters.
  *
  * @return A SourceClusterTable corresponding to clusters of sources having
  *         the given prototypical schema and slots.
  */
boost::shared_ptr<SourceClusterTable> const makeSourceClusterTable(
    lsst::afw::table::SourceTable const & prototype,
    boost::shared_ptr<lsst::afw::table::IdFactory> const & idFactory,
    SourceProcessingControl const & control
);

/** Process input sourcees, distributing them to one of 3 output catalogs.
  *     - \a sources:         sources suitable for spatial clustering.
  *     - \a badSources:      sources ignored for the purposes of clustering;
  *                           identified by one or more flags.
  *     - \a invalidSources:  sources which do not have sky-coordinates or centroids.
  *
  * Sources not lying in the given sky-tile are discarded.
  *
  * Existing fields in the input sources are modified as follows:
  *     - the MSBs of each source ID are set to the exposure ID
  *     - source footprints are discarded
  *     - source sky-coordinate errors are computed and added
  *     - exposure information (ID, filter ID, etc...) is added to each source.
  *     - the cluster coordinates for a source are set to the source coordinates.
  *
  * The output catalogs will typically have been constructed from tables obtained
  * via makeOutputSourceTable() - their schemas and slot mappings must all be 
  * identical. The input table slots must match output table slots, and the
  * input schema must be fully contained in the output schema.
  *
  * @param[in]  expSources     Single exposure sources to process.
  * @param[in]  expInfo        Exposure information.
  * @param[in]  skyTile        Sky-tile being processed.
  * @param[in]  control        Source processing parameters.
  * @param[in]  mapper         Maps between input and output source records.
  * @param[out] sources        Catalog for sources that will be clustered.
  * @param[out] badSources     Catalog for sources with bad measurement flags.
  * @param[out] invalidSources Catalog for sources with invalid measurements.
  */
void processSources(
    lsst::afw::table::SourceCatalog const & expSources,
    lsst::ap::match::ExposureInfo const & expInfo,
    lsst::ap::utils::PT1SkyTile const & skyTile,
    SourceProcessingControl const & control,
    lsst::afw::table::SchemaMapper const & mapper,
    lsst::afw::table::SourceCatalog & sources,
    lsst::afw::table::SourceCatalog & badSources,
    lsst::afw::table::SourceCatalog & invalidSources
);

/** Spatially cluster sources using the OPTICS algorithm. For details,
  * see the following paper:
  *
  * "OPTICS: Ordering Points To Identify the Clustering Structure".
  * Mihael Ankerst, Markus M. Breunig, Hans-Peter Kriegel, Jorg Sander (1999).
  * ACM SIGMOD international conference on Management of data.
  * ACM Press. pp. 49-60.
  *
  * http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.129.6542
  * http://en.wikipedia.org/wiki/OPTICS_algorithm
  * 
  * @param[in] sources  Sources to cluster.
  * @param[in] control  Clustering parameters.
  */
std::vector<lsst::afw::table::SourceCatalog> const cluster(
    lsst::afw::table::SourceCatalog const & sources,
    ClusteringControl const & control
);

/** Set the "<cluster>.id" and "<cluster>.coord" fields of each source
  * in the given catalog to the ID and sky-coordinates of the given cluster.
  * The "<cluster>" field name prefix is obtained from SourceProcessingControl.
  * If those fields are not setup in the given catalog, this function is a no-op.
  *
  * This is intended to support database ingest via qserv, where sources are
  * partitioned by their associated astrophysical object, and objects
  * (for which clusters are a temporary stand-in) are partitioned by position.
  * Denormalizing the output source schema by appending cluster position avoids
  * a potentially very expensive join during database ingest.
  *
  * @param[out] sources  Sources to update.
  * @param[in]  record   Cluster to obtain ID/sky-coordinates from.
  * @param[in]  control  Supplies cluster field name prefix.
  */  
void setClusterFields(
    lsst::afw::table::SourceCatalog & sources,
    SourceClusterRecord const & record,
    SourceProcessingControl const & control
);

}}} // namespace lsst::ap::cluster

#endif // LSST_AP_CLUSTER_CLUSTERING_H

