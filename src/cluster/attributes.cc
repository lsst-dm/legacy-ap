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
#include "lsst/ap/utils/SpatialUtils.h"


using lsst::pex::exceptions::InvalidParameterException;
using lsst::pex::exceptions::NotFoundException;

using lsst::afw::coord::IcrsCoord;

using lsst::afw::geom::HALFPI;
using lsst::afw::geom::PI;
using lsst::afw::geom::TWOPI;
using lsst::afw::geom::AffineTransform;
using lsst::afw::geom::LinearTransform;
using lsst::afw::geom::Point2D;
using lsst::afw::geom::radians;
using lsst::afw::geom::ellipses::Quadrupole;

using lsst::afw::image::Wcs;

using lsst::afw::table::Covariance;
using lsst::afw::table::Flag;
using lsst::afw::table::Flux;
using lsst::afw::table::Key;
using lsst::afw::table::SourceCatalog;
using lsst::afw::table::SourceRecord;
using lsst::afw::table::Schema;
using lsst::afw::table::Shape;

using lsst::ap::match::ExposureInfo;
using lsst::ap::match::ExposureInfoMap;


namespace lsst { namespace ap { namespace cluster {

// -- SourceAndExposure implementation --------

SourceAndExposure::SourceAndExposure() : _source(), _exposure(), _transform() { }

SourceAndExposure::SourceAndExposure(
    boost::shared_ptr<SourceRecord> const & source,
    boost::shared_ptr<ExposureInfo> const & exposure,
    AffineTransform const & transform
) : _source(source), _exposure(exposure), _transform(transform) { }

SourceAndExposure::~SourceAndExposure() { }


// -- Basic attributes ---------

namespace {

    // Compute unweighted position mean for the given sources.
    void meanCoord(SourceClusterRecord & cluster, SourceCatalog const & sources) {
        typedef SourceCatalog::const_iterator Iter;
        // compute point v such that the sum of the squared angular separations
        // between v and each source position is minimized.
        Eigen::Vector3d v = Eigen::Vector3d::Zero();
        for (Iter i = sources.begin(), e = sources.end(); i != e; ++i) {
            v += i->getCoord().getVector().asEigen();
        }
        IcrsCoord c = lsst::ap::utils::cartesianToIcrs(v);
        cluster.setCoord(c);
        Eigen::Matrix2d cov;
        if (sources.size() == 1) {
            // copy position error from source (if available)
            try {
                Key<Covariance<lsst::afw::table::Point<float> > > coordErrKey =
                    sources[0].getSchema()[std::string("coord.err")];
                cov = sources[0].get(coordErrKey).cast<double>();
            } catch (NotFoundException &) {
                cov = Eigen::Matrix2d::Constant(std::numeric_limits<double>::quiet_NaN());
            }
        } else {
            // compute sample covariance matrix
            cov = Eigen::Matrix2d::Zero();
            for (Iter i = sources.begin(), e = sources.end(); i != e; ++i) {
                double dra = i->getRa().asRadians() - c.getRa().asRadians();
                if (dra > PI) {
                    dra = dra - TWOPI;
                } else if (dra < -PI) {
                    dra = dra + TWOPI;
                }
                double dde = i->getDec().asRadians() - c.getDec().asRadians();
                cov(0,0) += dra * dra;
                cov(0,1) += dra * dde;
                cov(1,1) += dde * dde;
            }
            cov(1,0) = cov(0,1);
            // ... and map to unbiased estimate of the covariance
            // of the unweighted coordinate mean
            cov /= static_cast<double>(sources.size() - 1);
        }
        cluster.setCoordErr(cov);
    }


    // A tangent plane projection with the standard N,E basis.
    class NeTanProj {
    public:
        NeTanProj(IcrsCoord const & center);

        // Map N,E coordinates to sky coordinates.
        Eigen::Vector2d const neToSky(Eigen::Vector2d const & ne) const {
            return lsst::ap::utils::cartesianToSpherical(
                _origin + ne.x() * _east + ne.y() * _north);
        }

        // Map covariance matrix in the N,E tangent plane projection to sky coordinates.
        Eigen::Matrix2d const neToSky(Eigen::Matrix2d const & cov) const {
            return _jDiag.asDiagonal() * cov * _jDiag.asDiagonal();
        }

        // Compute a mapping from the given pixel space to this tangent plane.
        AffineTransform const pixelToNeTransform(
            Eigen::Vector2d const &point, Wcs const &wcs) const;

    private:
        IcrsCoord _center;
        Eigen::Vector3d _origin;
        Eigen::Vector3d _north;
        Eigen::Vector3d _east;
        Eigen::Vector2d _jDiag;
    };

    NeTanProj::NeTanProj(IcrsCoord const & center) : _center(center) {
        double sinLon = std::sin(center.getLongitude().asRadians());
        double cosLon = std::cos(center.getLongitude().asRadians());
        double sinLat = std::sin(center.getLatitude().asRadians());
        double cosLat = std::cos(center.getLatitude().asRadians());
        _origin = Eigen::Vector3d(cosLat * cosLon, cosLat * sinLon, sinLat);
        if (std::fabs(center.getLatitude().asRadians()) == HALFPI) {
            _north = Eigen::Vector3d(-1.0, 0.0, 0.0);
            _east = Eigen::Vector3d(0.0, 1.0, 0.0);
        } else {
            _north = Eigen::Vector3d(-sinLat * cosLon, -sinLat * sinLon, cosLat);
            _east = Eigen::Vector3d(-sinLon, cosLon, 0.0);
        }
        // The Jacobian of the N,E to sky transform evaluated at (x,y) = (0,0) is
        // [ 1/cos(lat), 0
        //   0,          1 ]
        // which is undefined when the tangent plane is centered at a pole.
        // Maybe the LSST schema should specify position errors in the N,E basis,
        // e.g. the usual semi major/minor axis lengths + position angle?
        _jDiag = Eigen::Vector2d(1.0 / cosLat, 1.0);
    }

    AffineTransform const NeTanProj::pixelToNeTransform(
        Eigen::Vector2d const &point, Wcs const &wcs) const
    {
        double const pix = 8.0;
        double const invPix = 1.0 / pix; // exact

        Eigen::Vector3d p = wcs.pixelToSky(
            point.x(), point.y())->toIcrs().getVector().asEigen();
        Eigen::Vector3d px = wcs.pixelToSky(
            point.x() + pix, point.y())->toIcrs().getVector().asEigen();
        Eigen::Vector3d py = wcs.pixelToSky(
            point.x(), point.y() + pix)->toIcrs().getVector().asEigen();

        Eigen::Vector2d pne(p.dot(_east), p.dot(_north));
        Eigen::Vector2d pxne(px.dot(_east), px.dot(_north));
        Eigen::Vector2d pyne(py.dot(_east), py.dot(_north));
        pne /= p.dot(_origin);
        pxne /= px.dot(_origin);
        pyne /= py.dot(_origin);

        Eigen::Matrix2d m;
        m.col(0) = (pxne - pne)*invPix;
        m.col(1) = (pyne - pne)*invPix;

        return AffineTransform(m, pne - m*point);
    }


    // Min, mean, and max value of a source field.
    struct StatTuple {
        StatTuple();
        void compute(
            std::vector<SourceAndExposure>::const_iterator first,
            std::vector<SourceAndExposure>::const_iterator last,
            Key<double> const & key);

        double min;
        double mean;
        double max;
    };

    StatTuple::StatTuple() :
        min(std::numeric_limits<double>::quiet_NaN()),
        mean(std::numeric_limits<double>::quiet_NaN()),
        max(std::numeric_limits<double>::quiet_NaN())
    { }

    void StatTuple::compute(
        std::vector<SourceAndExposure>::const_iterator first,
        std::vector<SourceAndExposure>::const_iterator last,
        Key<double> const & key)
    {
        if (first == last) {
            min = std::numeric_limits<double>::quiet_NaN();
            mean = std::numeric_limits<double>::quiet_NaN();
            max = std::numeric_limits<double>::quiet_NaN();
        } else {
            min = std::numeric_limits<double>::infinity();
            max = -min;
            mean = 0.0;
            double ns = static_cast<double>(last - first);
            for (; first != last; ++first) {
                double val = first->getSource()->get(key);
                if (val < min) {
                    min = val;
                }
                if (val > max) {
                    max = val;
                }
                mean += val;
            }
            mean /= ns;
            if (mean < min) {
                mean = min;
            } else if (mean > max) {
                mean = max;
            }
        }
    }


    // Compute weighted position mean for the given sources.
    void weightedMeanCoord(
        SourceClusterRecord & cluster,
        std::vector<SourceAndExposure> const & sources,
        NeTanProj const & proj)
    {
        typedef std::vector<SourceAndExposure>::const_iterator Iter;

        int n = 0;
        Eigen::Vector2d wmean = Eigen::Vector2d::Zero();
        Eigen::Matrix2d invCovSum = Eigen::Matrix2d::Zero();
        for (Iter i = sources.begin(), e = sources.end(); i != e; ++i) {
            SourceRecord const & r = *(i->getSource());
            Eigen::Matrix2d cov = r.getCentroidErr().cast<double>();
            if (lsst::utils::isnan(cov(0,0)) ||
                lsst::utils::isnan(cov(1,1))) {
                continue; // variance is NaN
            } else if (cov(0,0) <= 0.0 || cov(1,1) <= 0.0) {
                continue; // negative variance
            }
            // FIXME: for now, no measurement algorithms actually
            //        compute full covariance matrixes, and the
            //        off diagonal matrix elements are always NaN.
            //        Longer term, we will need some algorithm metadata
            //        to decide whether an off-diagonal NaN means a
            //        sample should be ignored because a computation failed,
            //        or whether it should be zeroed because the algorithm
            //        never computes it.
            if (cov(0,1) != cov(1,0)) {
                if (lsst::utils::isnan(cov(0,1)) && lsst::utils::isnan(cov(1,0))) {
                    cov(0,1) = 0.0; cov(1,0) = 0.0;
                } else {
                    continue; // covariance matrix not symmetric
                }
            }
            Point2D p = r.getCentroid();
            Eigen::Matrix2d m = i->getTransform().getLinear().getMatrix();
            Eigen::Matrix2d invCov = (m * cov * m.transpose()).inverse();
            invCovSum += invCov;
            p = i->getTransform()(p);
            wmean += invCov * p.asEigen();
            ++n;
        }
        if (n > 0) {
            Eigen::Matrix2d cov = invCovSum.inverse();
            wmean = cov * wmean;
            wmean = proj.neToSky(wmean);
            cluster.setWeightedMeanCoord(IcrsCoord(
                wmean.x() * radians, wmean.y() * radians));
            cluster.setWeightedMeanCoordErr(proj.neToSky(cov));
            // chi-squared correction for cov?
        } else {
            cluster.setWeightedMeanCoord(IcrsCoord(
                std::numeric_limits<double>::quiet_NaN() * radians,
                std::numeric_limits<double>::quiet_NaN() * radians));
            cluster.setWeightedMeanCoordErr(Eigen::Matrix2d::Constant(
                std::numeric_limits<double>::quiet_NaN()));
        }
        cluster.setWeightedMeanCoordCount(n);
    }


    // Comparison functor for sorting std::vector<SourceAndExposure> by filter ID.
    struct CmpSourceAndExposure {
        bool operator()(SourceAndExposure const & a, SourceAndExposure const & b) const {
            return a.getExposureInfo()->getFilter().getId() <
                   b.getExposureInfo()->getFilter().getId();
        }
    };

} // namespace <anonymous>


boost::shared_ptr<std::vector<SourceAndExposure> > const computeBasicAttributes(
    SourceClusterRecord & cluster,
    SourceCatalog const & sources,
    ExposureInfoMap const & exposures,
    std::string const & exposurePrefix)
{
    typedef SourceCatalog::const_iterator SourceIter;
    typedef std::vector<SourceAndExposure>::const_iterator SeIter;

    if (sources.empty()) {
        throw LSST_EXCEPT(InvalidParameterException, "No sources in cluster");
    }
    // unweighted mean coordinates
    meanCoord(cluster, sources);

    boost::shared_ptr<std::vector<SourceAndExposure> > se(
        new std::vector<SourceAndExposure>);
    NeTanProj const proj(cluster.getCoord());
    Schema const schema = sources.getSchema();
    Key<int64_t> const expIdKey = schema.find<int64_t>(exposurePrefix + ".id").key;
    Key<double> const expTimeMidKey = schema.find<double>(exposurePrefix + ".time.mid").key;
    // fill in vector of source,exposure pairs
    se->reserve(sources.size());
    for (SourceIter i = sources.begin(), e = sources.end(); i != e; ++i) {
        boost::shared_ptr<ExposureInfo> exp =
            const_cast<ExposureInfoMap &>(exposures).get(i->get(expIdKey));
        if (!exp) {
            throw LSST_EXCEPT(NotFoundException, "No ExposureInfo for source");
        }
        se->push_back(SourceAndExposure(
            i, exp, proj.pixelToNeTransform(i->getCentroid().asEigen(), *exp->getWcs())));
    }
    // source counts and observation time range
    StatTuple t;
    t.compute(se->begin(), se->end(), expTimeMidKey);
    cluster.setNumSources(static_cast<int>(sources.size()));
    cluster.setTimeMin(t.min);
    cluster.setTimeMean(t.mean);
    cluster.setTimeMax(t.max);
    // inverse variance weighted mean coordinates
    if (sources.getTable()->getCentroidKey().isValid() &&
        sources.getTable()->getCentroidErrKey().isValid()) {
        weightedMeanCoord(cluster, *se, proj);
    }
    // sort source,exposure vector by filter
    std::sort(se->begin(), se->end(), CmpSourceAndExposure());
    // per-filter source counts and observation time range
    SeIter i = se->begin();
    SeIter const e = se->end();
    while (i != e) {
        int const filterId = i->getExposureInfo()->getFilter().getId();
        SeIter j = i;
        for (++j; j != e && filterId == j->getExposureInfo()->getFilter().getId(); ++j) { }
        std::string const filterName = i->getExposureInfo()->getFilter().getName();
        cluster.setNumSources(filterName, static_cast<int>(j - i));
        t.compute(i, j, expTimeMidKey);
        cluster.setTimeMin(filterName, t.min);
        cluster.setTimeMax(filterName, t.max);
        i = j;
    }
    return se;
}


// -- Flux --------

namespace {

    typedef std::pair<double, double> FluxAndErr;
    typedef std::pair<double, double> FluxAndVariance;

    FluxAndErr const weightedMeanFlux(std::vector<FluxAndVariance> const & samples)
    {
        typedef std::vector<FluxAndVariance>::const_iterator Iter;
        if (samples.empty()) {
            return FluxAndErr(std::numeric_limits<double>::quiet_NaN(),
                              std::numeric_limits<double>::quiet_NaN());
        } else if (samples.size() == 1) {
            return FluxAndErr(samples[0].first, std::sqrt(samples[0].second));
        }
        // compute weighted mean x_w = sum(w_i * x_i) / sum(w_i), where
        // w_i = 1/Var(x_i)
        double wsum = 0.0;
        double wmean = 0.0;
        for (Iter i = samples.begin(), e = samples.end(); i != e; ++i) {
            double w = 1.0/i->second;
            wmean += w*i->first;
            wsum += w;
        }
        wmean /= wsum;
        // Var(sum(w_i * x_i) / sum(w_i)) = sum(w_i)^-1
        // Multiply by chi-squared over the number of degrees of freedom to
        // account for errors in the variances of the individual samples,
        // yielding Var(x_w) = sum(w_i)^-1 * sum(w_i * (x_i - x_w)) / (n - 1)
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
        throw LSST_EXCEPT(InvalidParameterException, "No sources in cluster");
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
        std::vector<FluxAndVariance> samples;
        samples.reserve(j - i);
        std::string const filter = i->getExposureInfo()->getFilter().getName();
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
        FluxAndErr mean = weightedMeanFlux(samples);
        Flux::MeasKey const measKey = clusterSchema[filter + "." + fluxDef];
        cluster.set(measKey, mean.first);
        Flux::ErrKey const errKey = clusterSchema[filter + "." + fluxDef + ".err"];
        cluster.set(errKey, mean.second);
        Key<int> const countKey = clusterSchema[filter + "." + fluxDef + ".count"];
        cluster.set(countKey, static_cast<int>(samples.size()));
    }
}


// -- Shape --------

namespace {

    typedef std::pair<Quadrupole, Eigen::Matrix<double, 3, 3, Eigen::DontAlign> > ShapeAndErr;

    void weightedMeanShape(ShapeAndErr & result, std::vector<ShapeAndErr> const & samples)
    {
        typedef std::vector<ShapeAndErr>::const_iterator Iter;
        if (samples.empty()) {
            result.first.setParameterVector(
                Eigen::Vector3d::Constant(std::numeric_limits<double>::quiet_NaN()));
            result.second = Eigen::Matrix3d::Constant(std::numeric_limits<double>::quiet_NaN());
            return;
        } else if (samples.size() == 1) {
            result = samples[0];
            return;
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
        result.first = Quadrupole(wmean);
        result.second = cov;
        // chi-squared correction for cov?
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
        throw LSST_EXCEPT(InvalidParameterException, "No sources in cluster");
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
        std::string const filter = i->getExposureInfo()->getFilter().getName();
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
            Eigen::Matrix3d cov = r->get(sourceShapeErrKey).cast<double>();
            if (lsst::utils::isnan(cov(0,0)) ||
                lsst::utils::isnan(cov(1,1)) ||
                lsst::utils::isnan(cov(2,2))) {
                continue; // covariance matrix contains NaNs
            } else if (cov(0,0) <= 0.0 ||
                       cov(1,1) <= 0.0 ||
                       cov(2,2) <= 0.0) {
                continue; // negative variance
            }
            // FIXME: for now, no measurement algorithms actually
            //        compute full covariance matrixes, and the
            //        off diagonal matrix elements are always NaN.
            //        Longer term, we will need some algorithm metadata
            //        to decide whether an off-diagonal NaN means a
            //        sample should be ignored because a computation failed,
            //        or whether it should be zeroed because the algorithm
            //        never computes it.
            if (cov(0,1) != cov(1,0) || cov(0,2) != cov(2,0) || cov(1,2) != cov(2,1)) {
                if (lsst::utils::isnan(cov(0,1)) && lsst::utils::isnan(cov(1,0)) &&
                    lsst::utils::isnan(cov(0,2)) && lsst::utils::isnan(cov(2,0)) &&
                    lsst::utils::isnan(cov(1,2)) && lsst::utils::isnan(cov(2,1))) {
                    cov(0,1) = 0.0; cov(1,0) = 0.0;
                    cov(0,2) = 0.0; cov(2,0) = 0.0;
                    cov(1,2) = 0.0; cov(2,1) = 0.0;
                } else {
                    continue; // covariance matrix not symmetric
                }
            }
            // transform moments and covariance matrix to N,E basis
            LinearTransform const * xform = &(i->getTransform().getLinear());
            Eigen::Matrix3d j = q.transform(*xform).d();
            q.transform(*xform).inPlace();
            cov = j * cov * j.transpose();
            samples.push_back(ShapeAndErr(q, cov));
        }
        ShapeAndErr mean;
        weightedMeanShape(mean, samples);
        Shape::MeasKey const measKey = clusterSchema[filter + "." + shapeDef];
        cluster.set(measKey, mean.first);
        Shape::ErrKey const errKey = clusterSchema[filter + "." + shapeDef + ".err"];
        cluster.set(errKey, mean.second);
        Key<int> const countKey = clusterSchema[filter + "." + shapeDef + ".count"];
        cluster.set(countKey, static_cast<int>(samples.size()));
    }
}

}}} // namespace lsst::ap::cluster

