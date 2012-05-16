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
  * @brief  Functions for computing cluster attributes.
  * @author Serge Monkewitz
  */
#ifndef LSST_AP_CLUSTER_ATTRIBUTES_H
#define LSST_AP_CLUSTER_ATTRIBUTES_H 

#include <vector>

#include "../match/ExposureInfo.h"
#include "SourceProcessingControl.h"
#include "SourceCluster.h"


namespace lsst { namespace ap { namespace cluster {

/** Helper class that bundles a source, it's originating exposure, and an
  * affine transform from the pixel space of that exposure to the pixel space
  * of a tangent plane projection centered at the fiducial cluster position
  * and with the standard north, east basis.
  */
class SourceAndExposure {
public:
    SourceAndExposure();

    SourceAndExposure(
        boost::shared_ptr<lsst::afw::table::SourceRecord> const & source,
        boost::shared_ptr<lsst::ap::match::ExposureInfo> const & exposure,
        lsst::afw::geom::AffineTransform const & transform);

    ~SourceAndExposure();

    boost::shared_ptr<lsst::afw::table::SourceRecord> const getSource() const {
        return _source;
    }
    boost::shared_ptr<lsst::ap::match::ExposureInfo> const getExposureInfo() const {
        return _exposure;
    }
    lsst::afw::geom::AffineTransform const & getTransform() const {
        return _transform;
    }

private:
    boost::shared_ptr<lsst::afw::table::SourceRecord> _source;
    boost::shared_ptr<lsst::ap::match::ExposureInfo> _exposure;
    lsst::afw::geom::AffineTransform _transform;
};


/** Compute basic cluster attributes.
  *
  * The Coord slot of @c cluster is set to the unweighted mean of the
  * coordinates of the associated sources. This fiducial cluster position
  * is used as the center of a N,E tangent plane projection that serves
  * as a common coordinate system for averaging other measurements. Exposure
  * information for each input source is looked up, and used to create a
  * vector of tuples that bundle a source with it's originating exposure
  * and an affine transform from the pixel space coordinate system centered
  * on the source's centroid to the common coordinate system.
  *
  * This vector of tuples is sorted such that sources originating from
  * exposures in the same filter are all consecutive.
  *
  * Finally, several other cluster attributes are computed/set:
  *     - "obs.count" : number of sources in the cluster
  *     - "obs.time.min", "obs.time.max", "obs.time.mean" :
  *       observation time range
  *     - "coord.weighted", "coord.weighted.err", "coord.weighted.count" :
  *       inverse variance weighted mean coordinate, error, and sample count
  *       (computed only if input sources have a valid CoordErr slot mapping).
  *     - "\<filter\>.obs.count" : filter specific source count
  *     - "\<filter\>.time.min", "\<filter\>.time.max" : filter specific
  *       observation time range
  *
  * @param[out] cluster         SourceClusterRecord to store results in.
  * @param[in]  sources         Sources in cluster.
  * @param[in]  exposures       Map from exposure ID to exposure information.
  * @param[in]  exposurePrefix  Exposure column name prefix.
  *
  * @return    A list of (source, exposure information, transform) tuples,
  *            grouped by filter.
  */
boost::shared_ptr<std::vector<SourceAndExposure> > const computeBasicAttributes(
    SourceClusterRecord & cluster,
    lsst::afw::table::SourceCatalog const & sources,
    lsst::ap::match::ExposureInfoMap const & exposures,
    std::string const & exposurePrefix
);

/** Compute per-filter calibrated flux means for a cluster of sources.
  *
  * The input source vector is assumed to be grouped by filter. Flux
  * measurments and errors are obtained from fields named '\<fluxDef\>'
  * and '\<fluxDef\>.err'. Results are stored in fields named 
  * '\<filter\>.\<fluxDef\>' and '\<filter\>.\<fluxDef\>.err',
  * with a sample count stored in '\<filter\>.\<fluxDef\>.count' field.
  * For an input source to be included in the mean, the following must hold:
  *     - the source must have finite flux
  *     - the source must have positive flux error
  *     - none of the flag fields in skipFlags must be set.
  *
  * If, after applying these filtering critera, 2 or more flux samples
  * are available, then the inverse variance weighted mean of the samples
  * (after calibration) is computed and stored. If only one sample is
  * available, the input flux and error are calibrated and stored directly.
  * If no samples are available, the output flux, error, and sample count
  * fields are set to NaN, NaN, and 0.
  *
  * @param[out] cluster    Cluster to store results in.
  * @param[in]  sources    List of sources and their exposures.
  * @param[in]  fluxDef    Name of flux field to operate on.
  * @param[in]  skipFlags  Flags marking sources to skip.
  * @param[in]  fluxScale  Flux scaling factor.
  */
void computeFluxMean(
    SourceClusterRecord & cluster,
    std::vector<SourceAndExposure> const & sources,
    std::string const & fluxDef,
    std::vector<lsst::afw::table::Key<lsst::afw::table::Flag > > const & skipFlags,
    double fluxScale
);

/** Compute per-filter shape means for a cluster.
  *
  * The input source vector is assumed to be grouped by filter. Shape 
  * measurments and errors are obtained from fields named '\<shapeDef\>'
  * and '\<shapeDef\>.err'. Results are stored in fields named 
  * '\<filter\>.\<shapeDef\>' and '\<filter\>.\<shapeDef\>.err',
  * with a sample count stored in '\<filter\>.\<shapeDef\>.count' field.
  * For an input source to be included in the mean, the following must hold:
  *     - the source must have finite moments
  *     - the source must have positive entries in the covariance
  *       matrix diagonals
  *     - off diagonal elements of the covariance matrix must not be NaN
  *     - the covariance matrix must be symmetric.
  *     - none of the flag fields in skipFlags must be set.
  *
  * If, after applying these filtering critera, 2 or more shape samples
  * are available, then the inverse variance weighted mean of the samples
  * (after transformation to a coordinate system S centered on the cluster
  * position with the standard N,E basis) is computed and stored. If only
  * one sample is available, the input shape and error are transformed to S
  * and stored directly. If no samples are available, the output moments and
  * covariance matrix are set to NaNs, and the sample count is zeroed.
  *
  * @param[out] cluster    Cluster to store results in.
  * @param[in]  sources    List of sources and their exposures.
  * @param[in]  shapeDef   Name of shape field to operate on.
  * @param[in]  skipFlags  Flags marking sources to skip.
  */
void computeShapeMean(
    SourceClusterRecord & cluster,
    std::vector<SourceAndExposure> const & sources,
    std::string const & shapeDef,
    std::vector<lsst::afw::table::Key<lsst::afw::table::Flag > > const & skipFlags
);

}}} // namespace lsst::ap::cluster

#endif // LSST_AP_CLUSTER_ATTRIBUTES_H

