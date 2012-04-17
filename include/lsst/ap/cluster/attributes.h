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
  * TODO
  * an unweighted position mean for the given cluster,
  * and sets earliest/mean/latest observation times.
  */
std::vector<SourceAndExposure> const computeBasicAttributes(
    SourceClusterRecord & cluster,
    lsst::afw::table::SourceCatalog & sources,
    lsst::ap::match::ExposureInfoMap & exposures,
    SourceProcessingControl const & control
);

/** Compute per-filter calibrated flux means for a cluster.
  *  
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
  */
void computeShapeMean(
    SourceClusterRecord & cluster,
    SourceProcessingControl const & control,
    std::vector<SourceAndExposure> const & sources,
    std::string const & shapeDef,
    std::vector<lsst::afw::table::Key<lsst::afw::table::Flag > > const & skipFlags
);

}}} // namespace lsst::ap::cluster

#endif // LSST_AP_CLUSTER_ATTRIBUTES_H

