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
  * @brief  Implementation of cluster attribute computation.
  * @author Serge Monkewitz
  */
#include "lsst/ap/cluster/attributes.h"

#include <utility>

#include "lsst/utils/ieee.h"
#include "lsst/afw/table/Schema.h"


using lsst::pex::exceptions::InvalidParameterException;

using lsst::afw::geom::AffineTransform;
using lsst::afw::geom::LinearTransform;
using lsst::afw::geom::ellipses::Quadrupole;

using lsst::afw::table::Flag;
using lsst::afw::table::Flux;
using lsst::afw::table::Key;
using lsst::afw::table::SourceRecord;
using lsst::afw::table::Schema;
using lsst::afw::table::Shape;

using lsst::ap::match::ExposureInfo;
using lsst::ap::match::ExposureInfoMap;


namespace lsst { namespace ap { namespace cluster {

// -- SourceAndExposure implementation --------

SourceAndExposure::SourceAndExposure(
    boost::shared_ptr<SourceRecord> const & source,
    boost::shared_ptr<ExposureInfo> const & exposure,
    AffineTransform const & transform
) : _source(source), _exposure(exposure), _transform(transform) { }

SourceAndExposure::~SourceAndExposure() { }


// -- Flux --------

namespace {

    typedef std::pair<double, double> FluxAndErr;

    FluxAndErr const inverseVarianceWeightedMean(std::vector<FluxAndErr> const & samples) {
        typedef std::vector<FluxAndErr>::const_iterator Iter;
        if (samples.empty()) {
            throw LSST_EXCEPT(InvalidParameterException, "No samples supplied");
        } else if (samples.size() == 1) {
            return std::make_pair(samples[0].first, std::sqrt(samples[0].second));
        }
        double wsum = 0.0;
        double wmean = 0.0;
        for (Iter i = samples.begin(), e = samples.end(); i != e; ++i) {
            double w = 1.0/i->second;
            wmean += w*i->first;
            wsum += w;
        }
        wmean /= wsum;
        double vwm = 0.0;
        for (Iter i = samples.begin(), e = samples.end(); i != e; ++i) {
            double d = i->first - wmean;
            vwm += (d*d)/i->second;
        }
        return FluxAndErr(wmean, std::sqrt(vwm/(wsum*(samples.size() - 1))));
    }

} // namespace <anonymous>


void computeFluxMean(
    SourceClusterRecord & cluster,
    std::vector<SourceAndExposure> const & sources,
    std::string const & fluxDef,
    std::vector<Key<Flag > > const & skipFlags,
    double fluxScale)
{
    typedef std::vector<SourceAndExposure>::const_iterator Iter;
    typedef std::vector<Key<Flag > >::const_iterator FlagIter;

    if (sources.empty()) {
        return;
    }
    Schema const sourceSchema = sources[0].getSource()->getSchema();
    Schema const clusterSchema = cluster.getSchema();
    Flux::MeasKey const sourceFluxKey = sourceSchema[fluxDef];
    Flux::ErrKey const sourceFluxErrKey = sourceSchema[fluxDef + ".err"];

    Iter i = sources.begin();
    Iter const e = sources.end();
    while (i != e) {
        // Determine range [i, j) containing sources from the same filter
        int const filterId = i->getExposureInfo()->getFilter().getId();
        Iter j = i;
        for (++j; j != e && filterId == j->getExposureInfo()->getFilter().getId(); ++j) { }
        // iterate over [i, j), storing all valid flux/error pairs
        // after calibrating them.
        std::vector<FluxAndErr> samples;
        samples.reserve(j - i);
        for (; i < j; ++i) {
            SourceRecord const * r = i->getSource().get();
            double flux = r->get(sourceFluxKey);
            double fluxErr = r->get(sourceFluxErrKey);
            if (!lsst::utils::isfinite(flux) || lsst::utils::isnan(fluxErr) || fluxErr <= 0.0) {
                continue;
            }
            bool skip = false;
            for (FlagIter f = skipFlags.begin(), fe = skipFlags.end(); f != fe; ++f) {
                if (r->get(*f)) {
                    skip = true;
                    break;
                }
            }
            if (skip) {
                continue;
            }
            samples.push_back(i->getExposureInfo()->calibrateFlux(flux, fluxErr, fluxScale));
        }
        if (!samples.empty()) {
            FluxAndErr mean = inverseVarianceWeightedMean(samples);
            std::string const filter = i->getExposureInfo()->getFilter().getName();
            Flux::MeasKey const measKey = clusterSchema[filter + "." + fluxDef];
            cluster.set(measKey, mean.first);
            Flux::ErrKey const errKey = clusterSchema[filter + "." + fluxDef + ".err"];
            cluster.set(errKey, mean.second);
            Key<int> const countKey = clusterSchema[filter + "." + fluxDef + ".count"];
            cluster.set(countKey, static_cast<int>(samples.size()));
        }
    }
}


// -- Shape --------

namespace {

    typedef std::pair<Quadrupole, Eigen::Matrix3d> ShapeAndErr;

    ShapeAndErr const inverseVarianceWeightedMean(std::vector<ShapeAndErr> const & samples) {
        typedef std::vector<ShapeAndErr>::const_iterator Iter;
        if (samples.empty()) {
            throw LSST_EXCEPT(InvalidParameterException, "No samples supplied");
        } else if (samples.size() == 1) {
            return samples[0];
        }

        Eigen::Matrix3d invCovSum = Eigen::Matrix3d::Zero();
        Eigen::Vector3d wmean = Eigen::Vector3d::Zero();
        for (Iter i = samples.begin(), e = samples.end(); i != e; ++i) {
            Eigen::Matrix3d m = i->second.inverse();
            invCovSum += m;
            wmean += m * (i->first.getParameterVector());
        }
        Eigen::Matrix3d cov = invCovSum.inverse();
        wmean = cov * wmean;
        return ShapeAndErr(Quadrupole(wmean), cov);
    }

} // namespace <anonymous>


void computeShapeMean(
    SourceClusterRecord & cluster,
    std::vector<SourceAndExposure> const & sources,
    std::string const & shapeDef,
    std::vector<lsst::afw::table::Key<lsst::afw::table::Flag > > const & skipFlags)
{
    typedef std::vector<SourceAndExposure>::const_iterator Iter;
    typedef std::vector<Key<Flag > >::const_iterator FlagIter;

    if (sources.empty()) {
        return;
    }
    Schema const sourceSchema = sources[0].getSource()->getSchema();
    Schema const clusterSchema = cluster.getSchema();
    Shape::MeasKey const sourceShapeKey = sourceSchema[shapeDef];
    Shape::ErrKey const sourceShapeErrKey = sourceSchema[shapeDef + ".err"];

    Iter i = sources.begin();
    Iter const e = sources.end();
    while (i != e) {
        // Determine range [i, j) containing sources from the same filter
        int const filterId = i->getExposureInfo()->getFilter().getId();
        Iter j = i;
        for (++j; j != e && filterId == j->getExposureInfo()->getFilter().getId(); ++j) { }
        // iterate over [i, j) storing all valid shape/error pairs in samples
        std::vector<ShapeAndErr> samples;
        samples.reserve(j - i); 
        for (; i < j; ++i) {
            SourceRecord const * r = i->getSource().get();
            bool skip = false;
            for (FlagIter f = skipFlags.begin(), fe = skipFlags.end(); f != fe; ++f) {
                if (r->get(*f)) {
                    skip = true;
                    break;
                }
            }
            if (skip) {
                continue;
            }
            Quadrupole q = r->get(sourceShapeKey);
            if (!lsst::utils::isfinite(q.getIxx()) ||
                !lsst::utils::isfinite(q.getIyy()) ||
                !lsst::utils::isfinite(q.getIxy())) {
                continue; // shape contains NaNs
            }
            Eigen::Matrix3d cov = r->get(sourceShapeErrKey);
            if (lsst::utils::isnan(cov(0,0)) ||
                lsst::utils::isnan(cov(0,1)) ||
                lsst::utils::isnan(cov(0,2)) ||
                lsst::utils::isnan(cov(1,1)) ||
                lsst::utils::isnan(cov(1,2)) ||
                lsst::utils::isnan(cov(2,2))) {
                continue; // covariance matrix contains NaNs
            } else if (cov(0,0) <= 0.0 ||
                       cov(1,1) <= 0.0 ||
                       cov(2,2) <= 0.0) {
                continue; // negative variance ?!
            }
            if (cov(0,1) != cov(1,0) ||
                cov(0,2) != cov(2,0) ||
                cov(1,2) != cov(2,1)) {
                continue; // not symmetric!
            }
            // transform quadrupole moments and covariance matrix to N,E basis
            LinearTransform const * xform = &(i->getTransform().getLinear());
            Eigen::Matrix3d j = q.transform(*xform).d();
            q.transform(*xform).inPlace();
            cov = j * cov * j.transpose();
            samples.push_back(ShapeAndErr(q, cov));
        }
        if (!samples.empty()) {
            ShapeAndErr mean = inverseVarianceWeightedMean(samples);
            std::string const filter = i->getExposureInfo()->getFilter().getName();
            Shape::MeasKey const measKey = clusterSchema[filter + "." + shapeDef];
            cluster.set(measKey, mean.first);
            Shape::ErrKey const errKey = clusterSchema[filter + "." + shapeDef + ".err"];
            cluster.set(errKey, mean.second);
            Key<int> const countKey = clusterSchema[filter + "." + shapeDef + ".count"];
            cluster.set(countKey, static_cast<int>(samples.size()));
        }
    }
}

}}} // namespace lsst::ap::cluster

