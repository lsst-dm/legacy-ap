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


/** Computes basic cluster attributes.
  *
  *
  *
  * @return    A list of source / exposure information pairs, grouped by filter.
  */
boost::shared_ptr<std::vector<SourceAndExposure> > const computeBasicAttributes(
    SourceClusterRecord & cluster,                      ///< @param[out] Cluster to store results in.
    lsst::afw::table::SourceCatalog const & sources,    ///< @param[in]  Sources in cluster.
    lsst::ap::match::ExposureInfoMap const & exposures, ///< @param[in]  Map from exposure ID to exposure information.
    std::string const & exposurePrefix                  ///< @param[in]  Exposure column name prefix.
);

/** Compute per-filter calibrated flux means for a cluster of sources.
  *
  * The input source vector is assumed to be grouped by filter. Flux
  * measurments and errors are obtained from fields named '<fluxDef>'
  * and '<fluxDef>.err'. Results are stored in fields named 
  * '<filter>.<fluxDef>' and '<filter>.<fluxDef>.err', with a sample count
  * stored in '<filter>.<fluxDef>.count' field. For an input source
  * to be included in the mean, the following must hold:
  *     - the source must have finite flux
  *     - the source must have positive flux error
  *     - none of the flag fields in skipFlags must be set.
  *
  * If, after applying these filtering critera, 2 or more flux samples
  * are available, then the inverse variance weighted mean of the samples
  * (after calibration) is computed and stored. If only one sample is
  * available, the input flux and error are calibrated and stored directly.
  * If no samples are available, the output flux, error, and sample count
  * fields are left untouched.
  */
void computeFluxMean(
    SourceClusterRecord & cluster,                  ///< @param[out] Cluster to store results in.
    std::vector<SourceAndExposure> const & sources, ///< @param[in]  List of sources and their exposures.
    std::string const & fluxDef,                    ///< @param[in]  Name of flux field to operate on.
    std::vector<lsst::afw::table::Key<lsst::afw::table::Flag > > const & skipFlags,
                                                    ///< @param[in]  Flags marking sources to skip.
    double fluxScale                                ///< @param[in]  Flux scaling factor.
);

/** Compute per-filter shape means for a cluster.
  *
  * The input source vector is assumed to be grouped by filter. Shape 
  * measurments and errors are obtained from fields named '<shapeDef>'
  * and '<shapeDef>.err'. Results are stored in fields named 
  * '<filter>.<shapeDef>' and '<filter>.<shapeDef>.err', with a sample
  * count stored in '<filter>.<shapeDef>.count' field. For an input
  * source to be included in the mean, the following must hold:
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
  * and stored directly. If no samples are available, the output moments,
  * covariance matrix, and sample count fields are left untouched.
  */
void computeShapeMean(
    SourceClusterRecord & cluster,                  ///< @param[out] Cluster to store results in.
    std::vector<SourceAndExposure> const & sources, ///< @param[in]  List of sources and their exposures.
    std::string const & shapeDef,                   ///< @param[in]  Name of shape field to operate on.
    std::vector<lsst::afw::table::Key<lsst::afw::table::Flag > > const & skipFlags
                                                    ///< @param[in]  Flags marking sources to skip.
);

}}} // namespace lsst::ap::cluster

#endif // LSST_AP_CLUSTER_ATTRIBUTES_H

