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
 
/** @file
  * @brief Implementation of high-level source clustering/attributes API. 
  *
  * @ingroup ap
  * @author Serge Monkewitz
  */
#include "lsst/ap/cluster/SourceCluster.h"

#include <cmath>
#include <limits>
#include <utility>

#include "Eigen/LU"

#include "boost/scoped_array.hpp"

#include "lsst/pex/exceptions.h"
#include "lsst/pex/policy.h"
#include "lsst/afw/coord/Coord.h"
#include "lsst/afw/geom/ellipses/Distortion.h"
#include "lsst/afw/image/ImageUtils.h"
#include "lsst/afw/image/Wcs.h"

#include "lsst/ap/Common.h"
#include "lsst/ap/cluster/optics/Metrics.h"
#include "lsst/ap/cluster/optics/Optics.cc"
#include "lsst/ap/utils/SpatialUtils.h"


namespace coord = lsst::afw::coord;
namespace detection = lsst::afw::detection;
namespace except = lsst::pex::exceptions;
namespace geom = lsst::afw::geom;
namespace image = lsst::afw::image;
namespace logging = lsst::pex::logging;
namespace policy = lsst::pex::policy;
namespace match = lsst::ap::match;

using std::sqrt;

using lsst::ap::utils::cartesianToSpherical;
using lsst::ap::utils::degrees;
using lsst::ap::utils::sphericalToCartesian;


namespace lsst { namespace ap { namespace cluster {

namespace {

typedef optics::Point<3, detection::Source::Ptr> Point;
typedef optics::Optics<3, detection::Source::Ptr> Optics;

/** @internal
  * Initializes a Point from a Source::Ptr.
  */
inline void initPoint(Point & entry,
                      lsst::afw::detection::Source::Ptr const & ptr)
{
    // Note: this will change once the new Source class heirarchy is in.
    double ra = ptr->getRa();
    double dec = ptr->getDec();
    double cosDec = std::cos(dec);
    entry.coords.coeffRef(0) = std::sin(ra) * cosDec;
    entry.coords.coeffRef(1) = std::cos(ra) * cosDec;
    entry.coords.coeffRef(2) = std::sin(dec);
    entry.data = &ptr;
}

/** @internal
  * Maximum number of sources that can be processed at once.
  */
unsigned int const MAX_SOURCES =
    static_cast<unsigned int>(std::numeric_limits<int>::max());

/** @internal
  * A tangent plane projection centered at a point with the standard N,E basis.
  */
class LSST_AP_LOCAL NeTanProj {
public:
    NeTanProj(lsst::afw::coord::Coord::Ptr center);

    /** Maps N,E coordinates to sky coordinates
      */
    Eigen::Vector2d const neToSky(Eigen::Vector2d const & ne) const {
        return cartesianToSpherical(
            _origin + ne.x()*_east + ne.y()*_north);
    }

    /** Maps a covariance matrix in the N,E tangent plane projection
      * to sky coordinates.
      */
    Eigen::Matrix2d const neToSky(Eigen::Matrix2d const & cov) const {
        return _jDiag.asDiagonal()*cov*_jDiag.asDiagonal();
    }

    lsst::afw::geom::AffineTransform const pixelToNeTransform(
        Eigen::Vector2d const &point,
        lsst::afw::image::Wcs const &wcs) const;

private:
    lsst::afw::coord::Coord::Ptr _center;
    Eigen::Vector3d _origin;
    Eigen::Vector3d _north;
    Eigen::Vector3d _east;
    Eigen::Vector2d _jDiag;
};

/** Creates a N,E tangent plane projection with the given center.
  */
NeTanProj::NeTanProj(lsst::afw::coord::Coord::Ptr center) :
    _center(center)
{
    double sinLon = std::sin(center->getLongitude(coord::RADIANS));
    double cosLon = std::cos(center->getLongitude(coord::RADIANS));
    double sinLat = std::sin(center->getLatitude(coord::RADIANS));
    double cosLat = std::cos(center->getLatitude(coord::RADIANS));
    _origin = Eigen::Vector3d(cosLat*cosLon, cosLat*sinLon, sinLat);
    if (std::fabs(center->getLatitude(coord::RADIANS)) == M_PI_2) {
        _north = Eigen::Vector3d(-1.0, 0.0, 0.0);
        _east = Eigen::Vector3d(0.0, 1.0, 0.0);
    } else {
        _north = Eigen::Vector3d(-sinLat*cosLon, -sinLat*sinLon, cosLat);
        _east = Eigen::Vector3d(-sinLon, cosLon, 0.0);
    }
    // The Jacobian of the N,E to sky transform evaluated at (x,y) = (0,0) is
    // [ 1/cos(lat), 0 
    //   0,          1 ]
    // which is undefined when the tangent plane is centered at a pole.
    // Maybe the LSST schema should specify position errors in the N,E basis,
    // or via the usual semi major/minor axis lengths and position angle?
    _jDiag = Eigen::Vector2d(1.0/cosLat, 1.0);
}

/** Computes an affine transform from the given pixel space to this tangent plane.
  */
lsst::afw::geom::AffineTransform const NeTanProj::pixelToNeTransform(
    Eigen::Vector2d const &point,    ///< Pixel coordinates to linearize around
    lsst::afw::image::Wcs const &wcs ///< WCS describing pixel/sky transformation
) const {
    double const pix = 8.0;
    double const invPix = 1.0/pix; // exact

    Eigen::Vector3d p = wcs.pixelToSky(
        point.x(), point.y())->getVector().asVector();
    Eigen::Vector3d px = wcs.pixelToSky(
        point.x() + pix, point.y())->getVector().asVector();
    Eigen::Vector3d py = wcs.pixelToSky(
        point.x(), point.y() + pix)->getVector().asVector();

    Eigen::Vector2d pne(p.dot(_east), p.dot(_north));
    Eigen::Vector2d pxne(px.dot(_east), px.dot(_north));
    Eigen::Vector2d pyne(py.dot(_east), py.dot(_north));
    pne /= p.dot(_origin);
    pxne /= px.dot(_origin);
    pyne /= py.dot(_origin);

    Eigen::Matrix2d m;
    m.col(0) = (pxne - pne)*invPix;
    m.col(1) = (pyne - pne)*invPix;

    return geom::AffineTransform(m, pne - m*point);
}

/** @internal
  * Computes the inverse variance weighted mean and the associated error
  * from a set of (sample, variance) pairs.
  */
std::pair<double, double> const inverseVarianceWeightedMean(
    std::vector<std::pair<double, double> > const & samples
) {
    typedef std::vector<std::pair<double, double> >::const_iterator Iter;
    if (samples.empty()) {
        throw LSST_EXCEPT(except::InvalidParameterException,
                          "No samples supplied");
    } else if (samples.size() == 1) {
        return std::make_pair(samples[0].first, sqrt(samples[0].second));
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
    return std::make_pair(wmean, sqrt(vwm/(wsum*(samples.size() - 1))));
}

/** @internal
  * Computes the position of a set of sources by combining the positions
  * using inverse variance weighting. The method is similar to the one
  * described by John E. Davis here:
  * ftp://space.mit.edu/pub/davis/misc/ellipse.pdf
  *
  * The difference is that the input is a set of pixel space position
  * covariance matrices, rather than position error ellipses in sky coordinates.
  * The covariance matrices are transformed directly to a common N,E tangent
  * plane projection, considerably simplifying the computation.
  */
bool combinePositions(std::vector<SourceAndExposure> const & sources,
                      NeTanProj const & proj,
                      Eigen::Vector2d & position,
                      Eigen::Matrix2d & covariance)
{
    typedef std::vector<SourceAndExposure>::const_iterator Iter;
    if (sources.empty()) {
        return false;
    }
    int ns = 0;
    Eigen::Vector2d wmean = Eigen::Vector2d::Zero();
    Eigen::Matrix2d invCovSum = Eigen::Matrix2d::Zero();
    for (Iter i = sources.begin(), e = sources.end(); i != e; ++i) {
        detection::Source const & s = *(i->getSource());
        if (s.isNull(detection::X_ASTROM_ERR) ||
            s.isNull(detection::Y_ASTROM_ERR)) {
            continue;
        }
        geom::PointD p = geom::makePointD(s.getXAstrom(), s.getYAstrom());
        Eigen::Vector2d var(s.getXAstromErr(), s.getYAstromErr());
        var = var.cwise().square();
        geom::AffineTransform const & transform = i->getTransform();
        Eigen::Matrix2d const & m = transform.getLinear().getMatrix();
        Eigen::Matrix2d invCov = (m*var.asDiagonal()*m.transpose()).inverse();
        invCovSum += invCov;
        p = transform(p);
        wmean += invCov*p.asVector();
        ++ns;
    }
    if (ns > 0) {
        Eigen::Matrix2d cov = invCovSum.inverse();
        wmean = cov*wmean;
        position = proj.neToSky(wmean);
        covariance = proj.neToSky(cov);
        return true;
    }
    return false;
}

} // namespace


/** Clusters a set of sources using the OPTICS algorithm. The following
  * parameters are read from @c policy:
  *
  * @li @c "epsilon" (double) : generating distance for clusters (arcsec).
  * @li @c "minNeighbors" (int) : minimum number of points that must be in an
  *     epsilon neighborhood of a point P for P to be assigned to a cluster.
  * @li @c "pointsPerLeaf" (int) : a performance tuning parameter that
  *     specifies the target number of points per leaf of the k-d tree used
  *     by OPTICS internally.
  *
  * @param[in] sources      The sources to cluster.
  * @param[in] policy       Policy containing clustering parameters.
  *
  * @return     A vector of SourceCluster objects.
  */
LSST_AP_API std::vector<lsst::afw::detection::SourceSet> const cluster(
    lsst::afw::detection::SourceSet const & sources,
    lsst::pex::policy::Policy::Ptr policy)
{
    typedef detection::SourceSet::const_iterator Iter;
    if (sources.size() > MAX_SOURCES) {
        throw LSST_EXCEPT(except::InvalidParameterException,
                          "too many sources to cluster");
    }
    double epsilon = policy->getDouble("epsilonArcsec");
    double leafExtentThreshold = policy->getDouble("leafExtentThresholdArcsec");
    if (epsilon < 0.0) {
        throw LSST_EXCEPT(except::InvalidParameterException,
                          "OPTICS epsilon (clustering distance) policy "
                          "parameter value is negative");
    }

    boost::scoped_array<Point> entries(new Point[sources.size()]);
    std::vector<detection::SourceSet> clusters;
    // Transform sources into a form the OPTICS implementation understands
    int i = 0;
    for (Iter s = sources.begin(), e = sources.end(); s != e; ++s, ++i) {
        initPoint(entries[i], *s);
    }
    if (i > 0) {
        // Convert epsilon and leafExtentThreshold to radians, and account
        // for the fact that our metric is the squared euclidian distance,
        // not angular separation.
        epsilon = std::sin(0.5 * RADIANS_PER_ARCSEC * epsilon);
        epsilon = 4.0 * epsilon * epsilon;
        if (leafExtentThreshold > 0.0) {
            leafExtentThreshold = std::sin(
                0.5 * RADIANS_PER_ARCSEC * leafExtentThreshold);
            leafExtentThreshold = 4.0 * leafExtentThreshold * leafExtentThreshold;
        }
        Optics optics(entries.get(),
                      i,
                      policy->getInt("minNeighbors"),
                      epsilon,
                      leafExtentThreshold,
                      policy->getInt("pointsPerLeaf"));
        optics.run(clusters, optics::SquaredEuclidianDistanceOverSphere());
    }
    return clusters;
}

/** Sets the object id and object position of each of the given sources
  * to the id of the given cluster. This is necessary because when exported
  * to the database, sources will be partitioned by their associated object
  * (cluster), and objects are partitioned by position.
  */
LSST_AP_API void updateSources(SourceClusterAttributes const & cluster,
                               lsst::afw::detection::SourceSet & sources)
{
    typedef detection::SourceSet::iterator Iter;
    for (Iter i = sources.begin(), e = sources.end(); i != e; ++i) {
        (*i)->setObjectId(cluster.getClusterId());
        (*i)->setRaObject(cluster.getRaPs());
        (*i)->setDecObject(cluster.getDecPs());
    }
}

/** Computes sky coordinates and errors for sources.  Sources without positions
  * (those with astrometric x/y equal to NULL/NaN), with pixel coordinates
  * outside the originating exposure or with negative flux or shape errors are
  * removed from @a sources and appended to @a invalidSources.
  */
LSST_AP_API void locateAndFilterSources(
    lsst::afw::detection::SourceSet & sources,
    lsst::afw::detection::SourceSet & invalidSources,
    lsst::ap::match::ExposureInfoMap const & exposures
) {
    typedef detection::SourceSet::iterator Iter;

    if (sources.empty()) {
        return;
    }

    logging::Log log(logging::Log::getDefaultLog(), "lsst.ap.cluster");
    Iter j = sources.begin();;
    Iter end = sources.end();
    double const xyMin = image::indexToPosition(0) - 0.5;
    int64_t id = (*j)->getAmpExposureId() + 1;
    double xMax = 0.0;
    double yMax = 0.0;
    image::Wcs::ConstPtr wcs;
 
    for (Iter i = j; i != end; ++i) {
        if ((*i)->getAmpExposureId() != id) {
            id = (*i)->getAmpExposureId();
            match::ExposureInfo::ConstPtr einfo = exposures.get(id);
            if (!einfo) {
                throw LSST_EXCEPT(except::NotFoundException,
                                  "No ExposureInfo found for source");
            }
            wcs = einfo->getWcs();
            xMax = image::indexToPosition(einfo->getWidth() - 1) + 0.5;
            yMax = image::indexToPosition(einfo->getHeight() - 1) + 0.5;
        }
        // Map peak x/y to sky coords (no errors)
        bool invalid = false;
        double x = (*i)->getXPeak();
        double y = (*i)->getYPeak();
        if (!lsst::utils::isnan(x) && !lsst::utils::isnan(y) &&
            !(*i)->isNull(detection::X_PEAK) &&
            !(*i)->isNull(detection::Y_PEAK)) {
            if (x < xyMin || x > xMax || y < xyMin || y > yMax) {
                log.format(logging::Log::WARN, "Peak (x, y) = (%.3f, %.3f) of "
                           "source %lld does not lie on exposure %lld", x, y,
                           static_cast<long long>((*i)->getSourceId()),
                           static_cast<long long>(id));
                invalid = true;
            } else {
                coord::Coord::Ptr sky = wcs->pixelToSky(x, y);
                (*i)->setRaPeak(sky->getLongitude(coord::RADIANS));
                (*i)->setDecPeak(sky->getLatitude(coord::RADIANS));
            }
        }
        // Map astrom x/y to sky coords
        x = (*i)->getXAstrom();
        y = (*i)->getYAstrom();
        if (lsst::utils::isnan(x) || lsst::utils::isnan(y) ||
            (*i)->isNull(detection::X_ASTROM) ||
            (*i)->isNull(detection::Y_ASTROM) ||
            x < xyMin || x > xMax || y < xyMin || y > yMax) {
            log.format(logging::Log::WARN, "Astrom (x, y) = (%.3f, %.3f) of "
                       "source %lld is NULL or does not lie on exposure %lld",
                       x, y,
                       static_cast<long long>((*i)->getSourceId()),
                       static_cast<long long>(id));
            invalid = true;
        } else {
            coord::Coord::Ptr sky = wcs->pixelToSky(x, y);
            (*i)->setRaAstrom(sky->getLongitude(coord::RADIANS));
            (*i)->setDecAstrom(sky->getLatitude(coord::RADIANS));
            Eigen::Vector2d v((*i)->getXAstromErr(), (*i)->getYAstromErr());
            if (!lsst::utils::isnan(v.x()) && !lsst::utils::isnan(v.y()) &&
                !(*i)->isNull(detection::X_ASTROM_ERR) &&
                !(*i)->isNull(detection::Y_ASTROM_ERR)) {
                if (v.x() < 0.0 || v.y() < 0.0) {
                    log.format(logging::Log::WARN, "Source %lld in exposure "
                               "%lld has negative astrom x/y error",
                               static_cast<long long>((*i)->getSourceId()),
                               static_cast<long long>(id));
                    invalid = true;
                } else {
                    geom::AffineTransform xform =
                        wcs->linearizePixelToSky(sky, coord::RADIANS);
                    v = (xform.getLinear().getMatrix().cwise() *
                         xform.getLinear().getMatrix()) * v.cwise().square();
                    (*i)->setRaAstromErr(sqrt(v.x()));
                    (*i)->setDecAstromErr(sqrt(v.y()));
                }
            }
        }
        // Map flux x/y to sky coords.
        x = (*i)->getXFlux();
        y = (*i)->getYFlux();
        if (!lsst::utils::isnan(x) && !lsst::utils::isnan(y) &&
            !(*i)->isNull(detection::X_FLUX) &&
            !(*i)->isNull(detection::Y_FLUX)) {
            if (x < xyMin || x > xMax || y < xyMin || y > yMax) {
                log.format(logging::Log::WARN, "Flux (x, y) = (%.3f, %.3f) of "
                           "source %lld does not lie on exposure %lld", x, y,
                           static_cast<long long>((*i)->getSourceId()),
                           static_cast<long long>(id));
                invalid = true;
            } else {
                coord::Coord::Ptr sky = wcs->pixelToSky(x, y);
                (*i)->setRaFlux(sky->getLongitude(coord::RADIANS));
                (*i)->setDecFlux(sky->getLatitude(coord::RADIANS));
                Eigen::Vector2d v((*i)->getXFluxErr(), (*i)->getYFluxErr());
                if (!lsst::utils::isnan(v.x()) && !lsst::utils::isnan(v.y()) &&
                    !(*i)->isNull(detection::X_FLUX_ERR) &&
                    !(*i)->isNull(detection::Y_FLUX_ERR)) {
                    if (v.x() < 0.0 || v.y() < 0.0) {
                        log.format(logging::Log::WARN, "Source %lld in "
                                   "exposure %lld has negative flux x/y error",
                                   static_cast<long long>((*i)->getSourceId()),
                                   static_cast<long long>(id));
                        invalid = true;
                    } else {
                        geom::AffineTransform xform =
                            wcs->linearizePixelToSky(sky, coord::RADIANS);
                        v = (xform.getLinear().getMatrix().cwise() *
                             xform.getLinear().getMatrix()) * v.cwise().square();
                        (*i)->setRaFluxErr(sqrt(v.x()));
                        (*i)->setDecFluxErr(sqrt(v.y()));
                    }
                }
            }
        }
        // Check that PSF and aperture flux errors are positive
        if ((!lsst::utils::isnan((*i)->getPsfFluxErr()) &&
             (*i)->getPsfFluxErr() < 0.0) ||
            (!lsst::utils::isnan((*i)->getApFluxErr()) &&
             (*i)->getApFluxErr() < 0.0)) {
            log.format(logging::Log::WARN, "Source %lld in exposure %lld has "
                       "negative PSF and/or aperture flux error",
                       static_cast<long long>((*i)->getSourceId()),
                       static_cast<long long>(id));
            invalid = true;
        }
        // Check that errors on adaptive second moments are positive
        if ((!lsst::utils::isnan((*i)->getIxxErr()) &&
             !(*i)->isNull(detection::IXX_ERR) && (*i)->getIxxErr() < 0.0) ||
            (!lsst::utils::isnan((*i)->getIyyErr()) &&
             !(*i)->isNull(detection::IYY_ERR) && (*i)->getIyyErr() < 0.0) ||
            (!lsst::utils::isnan((*i)->getIxyErr()) &&
             !(*i)->isNull(detection::IXY_ERR) && (*i)->getIxyErr() == 0.0)) {
            log.format(logging::Log::WARN, "Source %lld in exposure %lld has "
                       "negative IXX and/or IYY error or 0 IXY error",
                       static_cast<long long>((*i)->getSourceId()),
                       static_cast<long long>(id));
            invalid = true;
        }
        // Store source in appropriate vector 
        if (invalid) {
            invalidSources.push_back(*i);
        } else {
            if (i != j) {
                *j = *i;
            }
            ++j;
        }
    }
    sources.erase(j, end);
}

/** Removes invalid sources without positions from @a sources and appends them
  * to @a invalidSources. A source is invalid if it has NaN or out-of-bounds
  * right ascension and/or declination, or negative errors.
  */
LSST_AP_API void segregateInvalidSources(
    lsst::afw::detection::SourceSet & sources,
    lsst::afw::detection::SourceSet & invalidSources)
{
    typedef detection::SourceSet::iterator Iter;

    if (sources.empty()) {
        return;
    }

    logging::Log log(logging::Log::getDefaultLog(), "lsst.ap.cluster");
    Iter j = sources.begin();
    Iter end = sources.end();
    for (Iter i = j; i != end; ++i) {
        bool invalid = false;
        double ra = (*i)->getRa();
        double dec = (*i)->getDec();
        if (lsst::utils::isnan(ra) ||
            lsst::utils::isnan(dec) ||
            ra < 0.0 || ra >= 2.0 * M_PI ||
            dec < -M_PI_2 || dec > M_PI_2) {
            log.format(logging::Log::WARN, "Source %lld in exposure %lld has "
                       "invalid position (NaN or out-of-bounds coordinate value(s))",
                       static_cast<long long>((*i)->getSourceId()),
                       static_cast<long long>((*i)->getAmpExposureId()));
            invalid = true;
        }
        // Check that X/Y astrom errors are non-negative
        double xerr = (*i)->getXAstromErr();
        double yerr = (*i)->getYAstromErr();
        if ((!lsst::utils::isnan(xerr) &&
             !(*i)->isNull(detection::X_ASTROM_ERR) && xerr < 0.0) ||
            (!lsst::utils::isnan(yerr) &&
             !(*i)->isNull(detection::Y_ASTROM_ERR) && yerr < 0.0)) {
            log.format(logging::Log::WARN, "Source %lld in exposure %lld has "
                       "negative X and/or Y astrom error",
                       static_cast<long long>((*i)->getSourceId()),
                       static_cast<long long>((*i)->getAmpExposureId()));
            invalid = true;
        }
        // Check that X/Y flux errors are non-negative
        xerr = (*i)->getXFluxErr();
        yerr = (*i)->getYFluxErr();
        if ((!lsst::utils::isnan(xerr) &&
             !(*i)->isNull(detection::X_FLUX_ERR) && xerr < 0.0) ||
            (!lsst::utils::isnan(yerr) &&
             !(*i)->isNull(detection::Y_FLUX_ERR) && yerr < 0.0)) {
            log.format(logging::Log::WARN, "Source %lld in exposure %lld has "
                       "negative X and/or Y flux error",
                       static_cast<long long>((*i)->getSourceId()),
                       static_cast<long long>((*i)->getAmpExposureId()));
            invalid = true;
        }
        // Check that PSF and aperture flux errors are non-negative
        double pfe = (*i)->getPsfFluxErr();
        double afe = (*i)->getApFluxErr();
        if ((!lsst::utils::isnan(pfe) && pfe < 0.0) ||
            (!lsst::utils::isnan(afe) && afe < 0.0)) {
            log.format(logging::Log::WARN, "Source %lld in exposure %lld has "
                       "negative PSF and/or aperture flux error",
                       static_cast<long long>((*i)->getSourceId()),
                       static_cast<long long>((*i)->getAmpExposureId()));
            invalid = true;
        }
        // Check that errors on adaptive second moments are non-negative
        double ixxe = (*i)->getIxxErr();
        double iyye = (*i)->getIyyErr();
        double ixye = (*i)->getIxyErr();
        if ((!lsst::utils::isnan(ixxe) &&
             !(*i)->isNull(detection::IXX_ERR) && ixxe < 0.0) ||
            (!lsst::utils::isnan(iyye) &&
             !(*i)->isNull(detection::IYY_ERR) && iyye < 0.0) ||
            (!lsst::utils::isnan(ixye) &&
             !(*i)->isNull(detection::IXY_ERR) && ixye == 0.0)) {
            log.format(logging::Log::WARN, "Source %lld in exposure %lld has "
                       "negative IXX and/or IYY error, or 0 IXY error",
                       static_cast<long long>((*i)->getSourceId()),
                       static_cast<long long>((*i)->getAmpExposureId()));
            invalid = true;
        }
        if (invalid) {
            invalidSources.push_back(*i);
        } else {
            if (i != j) {
                *j = *i;
            }
            ++j;
        }
    }
    sources.erase(j, end);
}

/** Removes "bad" sources from @a sources and appends them to @a badSources.
  * A source is considered bad when one or more of the bits set in its detection
  * flags matches a bit set in @a badSourceMask.
  */
LSST_AP_API void segregateBadSources(
    lsst::afw::detection::SourceSet & sources,
    lsst::afw::detection::SourceSet & badSources,
    int const badSourceMask
) {
    typedef detection::SourceSet::iterator Iter;

    Iter j = sources.begin();
    Iter end = sources.end();
    for (Iter i = j; i != end; ++i) {
        if (((*i)->getFlagForDetection() & badSourceMask) != 0) {
            badSources.push_back(*i);
        } else {
            if (i != j) {
                *j = *i;
            }
            ++j;
        }
    }
    sources.erase(j, end);
}

/** Sets the object id of each of the sources to NULL. Also sets
  * the object position of each of bad source to the source position.
  * This allows bad sources to be stored in the same table as good
  * sources (which are partitioned by the position of their associated
  * objects/clusters).
  */
LSST_AP_API void updateUnclusteredSources(
    lsst::afw::detection::SourceSet & badSources
) {
    typedef detection::SourceSet::iterator Iter;
    for (Iter i = badSources.begin(), e = badSources.end(); i != e; ++i) {
        (*i)->setNull(detection::OBJECT_ID, true);
        (*i)->setRaObject((*i)->getRa());
        (*i)->setDecObject((*i)->getDec());
    }
}


// -- PerFilterSourceClusterAttributes ----

PerFilterSourceClusterAttributes::PerFilterSourceClusterAttributes() :
    _filterId(0),
    _numObs(0),
    _flags(0),
    _earliestObsTime(0.0), _latestObsTime(0.0),
    _flux(), _fluxSigma(),
    _e1(), _e2(), _radius(),
    _e1Sigma(), _e2Sigma(), _radiusSigma()
{ }

/** Sets the flux of a source cluster in a particular filter to the
  * calibrated flux of the given source.
  */
void PerFilterSourceClusterAttributes::computeFlux(
    SourceAndExposure const & source, ///< Source to compute flux from.
    double fluxScale,         ///< Flux scaling factor, must be \> 0.0
    int fluxIgnoreMask        ///< Detection flag bitmask identifying sources
                              ///  that should be ignored when determining
                              ///  cluster fluxes.
) {
    detection::Source const & s = *source.getSource();
    match::ExposureInfo const & e = *source.getExposureInfo();
    if (!lsst::utils::isnan(s.getPsfFlux()) &&
        s.getPsfFluxErr() > 0.0 &&
        (s.getFlagForDetection() & fluxIgnoreMask) == 0) {
        std::pair<double, double> f =
            e.calibrateFlux(s.getPsfFlux(), s.getPsfFluxErr(), fluxScale);
        setFlux(static_cast<float>(f.first),
                static_cast<float>(sqrt(f.second)));
        setNumFluxSamples(1);
    }
}

/** Sets the flux of a source cluster in a particular filter to the inverse
  * variance weighted sample mean of the calibrated fluxes of the given
  * sources. The flux error is set to the error of the weighted flux mean.
  */
void PerFilterSourceClusterAttributes::computeFlux(
    std::vector<SourceAndExposure> const & sources, ///< Sources to compute
                                                    ///  flux from.
    double fluxScale,         ///< Flux scaling factor, must be \> 0.0
    int fluxIgnoreMask        ///< Detection flag bitmask identifying sources
                              ///  that should be ignored when determining
                              ///  cluster fluxes.
) {
    typedef std::vector<SourceAndExposure>::const_iterator Iter;
    if (sources.empty()) {
        return;
    }
    std::vector<std::pair<double, double> > samples;
    samples.reserve(sources.size());

    for (Iter i = sources.begin(), e = sources.end(); i != e; ++i) {
        detection::Source const & s = *(i->getSource());
        if (!lsst::utils::isnan(s.getPsfFlux()) &&
            (s.getFlagForDetection() & fluxIgnoreMask) == 0 &&
            s.getPsfFluxErr() > 0.0) {
            match::ExposureInfo const & exp = *(i->getExposureInfo());
            samples.push_back(exp.calibrateFlux(s.getPsfFlux(),
                                                s.getPsfFluxErr(),
                                                fluxScale));
        }
    }
    if (!samples.empty()) {
        setNumFluxSamples(static_cast<int>(samples.size()));
        std::pair<double, double> f = inverseVarianceWeightedMean(samples);
        setFlux(static_cast<float>(f.first), static_cast<float>(f.second));
    }
}

/** Sets the ellipticities of a source cluster in a particular filter to the
  * ellipticities derived from the adaptive moments of the given source. These
  * are expressed in a tangent plane projection centered at the fiducial
  * cluster position with the standard N,E basis. The scaling of the projection
  * is such that the mean ellipse radius R is expressed in radians.
  */
void PerFilterSourceClusterAttributes::computeEllipticity(
    SourceAndExposure const & source, ///< Source to compute ellipticities from
    int ellipticityIgnoreMask ///< Detection flag bitmask identifying sources
                              ///  that should be ignored when determining
                              ///  ellipticities.
) {
    using geom::ellipses::Distortion;

    detection::Source const & s = *source.getSource();
    if (!s.isNull(detection::IXX) &&
        !s.isNull(detection::IYY) &&
        !s.isNull(detection::IXY) &&
        !lsst::utils::isnan(s.getIxx()) &&
        !lsst::utils::isnan(s.getIyy()) &&
        !lsst::utils::isnan(s.getIxy()) &&
        (s.getFlagForDetection() & ellipticityIgnoreMask) == 0) {

        double ixx = s.getIxx();
        double iyy = s.getIyy();
        double ixy = s.getIxy();
        double t = ixx + iyy;
        if (t == 0.0) {
            return;
        }
        Distortion distortion((ixx - iyy)/t, 2.0*ixy/t, sqrt(t));
        distortion.transform(source.getTransform()).inPlace();
        // store results
        setNumEllipticitySamples(1);
        setEllipticity(static_cast<float>(distortion[Distortion::E1]),
                       static_cast<float>(distortion[Distortion::E2]),
                       static_cast<float>(distortion[Distortion::R]),
                       NullOr<float>(), NullOr<float>(), NullOr<float>());
    }
#if 0
    if (!s.isNull(detection::IXX) &&
        !s.isNull(detection::IYY) &&
        !s.isNull(detection::IXY) &&
        !lsst::utils::isnan(s.getIxx()) &&
        !lsst::utils::isnan(s.getIyy()) &&
        !lsst::utils::isnan(s.getIxy()) &&
        !s.isNull(detection::IXX_ERR) &&
        !s.isNull(detection::IYY_ERR) &&
        !s.isNull(detection::IXY_ERR) &&
        s.getIxxErr() > 0.0 &&
        s.getIyyErr() > 0.0 &&
        s.getIxyErr() != 0.0 &&
        (s.getFlagForDetection() & ellipticityIgnoreMask) == 0) {

        double ixx = s.getIxx();
        double iyy = s.getIyy();
        double ixy = s.getIxy();
        double t = ixx + iyy;
        if (t == 0.0) {
            return;
        }
        double tinv = 1.0/t;
        double rt = sqrt(t);
        Distortion distortion((ixx - iyy)*tinv, 2.0*ixy*tinv, rt);
        // v.asDiagonal() gives the covariance matrix of (IXX, IYY, IXY) -
        // this is incomplete, but these are the only bits of the actual
        // covariance matrix that survive persistence.
        Eigen::Vector3d v(s.getIxxErr(), s.getIyyErr(), s.getIxyErr());
        v = v.cwise().square();
        // Compute Jacobian J of mapping from (IXX, IYY, IXY) to (E1, E2, R)
        Eigen::Matrix3d j;
        j(0,2) = 0.0;
        j(1,2) = 2.0*tinv;
        j(2,2) = 0.0;
        double d = 2.0*tinv*tinv;
        j(0,0) = d*iyy;
        j(0,1) = - d*ixx;
        j(1,0) = - d*ixy;
        j(1,1) = - d*ixy;
        d = 0.5/rt;
        j(2,0) = d;
        j(2,1) = d;
        // Get covariance matrix of (E1, E2, R) by linearizing
        // moment to ellipticity mapping using J
        Eigen::Matrix3d cov = j*v.asDiagonal()*j.transpose();
        // transform to N,E coordinate system
        j = distortion.transform(source.getTransform()).d();
        distortion.transform(source.getTransform()).inPlace();
        v = (j*cov*j.transpose()).diagonal();
        // store results
        setNumEllipticitySamples(1);
        setEllipticity(static_cast<float>(distortion[Distortion::E1]),
                       static_cast<float>(distortion[Distortion::E2]),
                       static_cast<float>(distortion[Distortion::R]),
                       static_cast<float>(sqrt(v.x())),
                       static_cast<float>(sqrt(v.y())),
                       static_cast<float>(sqrt(v.z())));
    }
#endif
}

/** Sets the ellipticities of a source cluster in a particular filter to
  * the inverse variance weighted mean of the ellipticities derived from
  * the adaptive moments of the given sources. These are all expressed 
  * in a tangent plane projection centered at the fiducial cluster position
  * with the standard N,E basis. The scaling of the projection is such that
  * that the mean ellipse radius R is expressed in radians.
  */
void PerFilterSourceClusterAttributes::computeEllipticity(
    std::vector<SourceAndExposure> const & sources, ///< Sources to compute
                                                    ///  ellipticities from
    int ellipticityIgnoreMask ///< Detection flag bitmask identifying sources
                              ///  that should be ignored when determining
                              ///  ellipticities.
) {
    using geom::ellipses::Distortion;
    typedef std::vector<SourceAndExposure>::const_iterator SeIter;
    typedef std::vector<Distortion>::const_iterator DistortionIter;

    if (sources.empty()) {
        return;
    }
    std::vector<Distortion> samples;
    samples.reserve(sources.size());
    Eigen::Vector3d mean = Eigen::Vector3d::Zero();
    for (SeIter i = sources.begin(), e = sources.end(); i != e; ++i) {
        detection::Source const & s = *(i->getSource());
        if (!s.isNull(detection::IXX) &&
            !s.isNull(detection::IYY) &&
            !s.isNull(detection::IXY) &&
            !lsst::utils::isnan(s.getIxx()) &&
            !lsst::utils::isnan(s.getIyy()) &&
            !lsst::utils::isnan(s.getIxy()) &&
            (s.getFlagForDetection() & ellipticityIgnoreMask) == 0) {

            double ixx = s.getIxx();
            double iyy = s.getIyy();
            double ixy = s.getIxy();
            double t = ixx + iyy;
            if (t == 0.0) {
                continue;
            }
            Distortion distortion((ixx - iyy)/t, 2.0*ixy/t, sqrt(t));
            distortion.transform(i->getTransform()).inPlace();
            mean += distortion.getVector();
            samples.push_back(distortion);
        }
    }
    if (!samples.empty()) {
        setNumEllipticitySamples(static_cast<int>(samples.size()));
        if (samples.size() == 1) {
            setEllipticity(static_cast<float>(mean[Distortion::E1]),
                           static_cast<float>(mean[Distortion::E2]),
                           static_cast<float>(mean[Distortion::R]),
                           NullOr<float>(),
                           NullOr<float>(),
                           NullOr<float>());
        } else {
            double n = samples.size();
            mean /= n;
            Eigen::Vector3d var = Eigen::Vector3d::Zero();
            for (DistortionIter i = samples.begin(), e = samples.end();
                 i != e; ++i) {
                var += (i->getVector() - mean).cwise().square();
            }
            setEllipticity(
                static_cast<float>(mean[Distortion::E1]),
                static_cast<float>(mean[Distortion::E2]),
                static_cast<float>(mean[Distortion::R]),
                static_cast<float>(sqrt(var[Distortion::E1]/(n*(n - 1.0)))),
                static_cast<float>(sqrt(var[Distortion::E2]/(n*(n - 1.0)))),
                static_cast<float>(sqrt(var[Distortion::R]/(n*(n - 1.0)))));
        }
    }

#if 0
    int ns = 0;
    Eigen::Matrix3d invCovSum = Eigen::Matrix3d::Zero();
    Eigen::Vector3d wmean = Eigen::Vector3d::Zero();
    for (SeIter i = sources.begin(), e = sources.end(); i != e; ++i) {
        detection::Source const & s = *(i->getSource());
        if (!s.isNull(detection::IXX) &&
            !s.isNull(detection::IYY) &&
            !s.isNull(detection::IXY) &&
            !lsst::utils::isnan(s.getIxx()) &&
            !lsst::utils::isnan(s.getIyy()) &&
            !lsst::utils::isnan(s.getIxy()) &&
            !s.isNull(detection::IXX_ERR) &&
            !s.isNull(detection::IYY_ERR) &&
            !s.isNull(detection::IXY_ERR) &&
            s.getIxxErr() > 0.0 &&
            s.getIyyErr() > 0.0 &&
            s.getIxyErr() != 0.0 &&
            (s.getFlagForDetection() & ellipticityIgnoreMask) == 0) {

            double ixx = s.getIxx();
            double iyy = s.getIyy();
            double ixy = s.getIxy();
            double t = ixx + iyy;
            if (t == 0.0) {
                continue;
            }
            double tinv = 1.0/t;
            double rt = sqrt(t);
            Distortion distortion((ixx - iyy)*tinv, 2.0*ixy*tinv, rt);
            // v.asDiagonal() gives the covariance matrix of (IXX, IYY, IXY) -
            // this is incomplete, but these are the only bits of the actual
            // covariance matrix that survive persistence.
            Eigen::Vector3d v(s.getIxxErr(), s.getIyyErr(), s.getIxyErr());
            v = v.cwise().square();
            // Compute Jacobian of mapping from (IXX, IYY, IXY) to (E1, E2, R)
            Eigen::Matrix3d j;
            j(0,2) = 0.0;
            j(1,2) = 2.0*tinv;
            j(2,2) = 0.0;
            double d = 2.0*tinv*tinv;
            j(0,0) = d*iyy;
            j(0,1) = - d*ixx;
            j(1,0) = - d*ixy;
            j(1,1) = - d*ixy;
            d = 0.5/rt;
            j(2,0) = d;
            j(2,1) = d;
            // Get covariance matrix of (E1, E2, R)
            Eigen::Matrix3d cov = j*v.asDiagonal()*j.transpose();
            // transform to N,E coordinate system
            j = distortion.transform(i->getTransform()).d();
            distortion.transform(i->getTransform()).inPlace();
            Eigen::Matrix3d invCov = (j*cov*j.transpose()).inverse();
            // add inverse covariance matrix to running sum, multiply sample
            // by inverse cov and add to weighted mean.
            invCovSum += invCov;
            wmean += invCov*distortion.getVector();
            ++ns;
        }
    }
    if (ns > 0) {
        Eigen::Matrix3d cov = invCovSum.inverse();
        wmean = cov*wmean;
        setNumEllipticitySamples(ns);
        setEllipticity(static_cast<float>(wmean.x()),
                       static_cast<float>(wmean.y()),
                       static_cast<float>(wmean.z()),
                       static_cast<float>(sqrt(cov(0,0))),
                       static_cast<float>(sqrt(cov(1,1))),
                       static_cast<float>(sqrt(cov(2,2))));
    }
#endif
}

PerFilterSourceClusterAttributes::~PerFilterSourceClusterAttributes() { }

/** Sets the earliest and latest cluster observation times in this filter.
  */
void PerFilterSourceClusterAttributes::setObsTimeRange(double earliest,
                                                       double latest)
{
    if (earliest > latest) {
        throw LSST_EXCEPT(except::InvalidParameterException, "earliest "
                          "observation time is after latest observation time");
    }
    _earliestObsTime = earliest;
    _latestObsTime = latest;
}

/** Returns the number of smaples (sources) used to determine the PSF
  * flux sample mean.
  *
  * @li If this number is zero, then none of the sources satisified the
  *     the criteria for being included in the PSF flux sample mean, and
  *     both the flux and its uncertainty are invalid (NULL/NaN).
  * @li If this number is one, then the flux uncertainty is set to the
  *     uncertainty of the PSF flux for that single source, rather than
  *     to an estimate of the standard deviation of the sample mean.
  */
int PerFilterSourceClusterAttributes::getNumFluxSamples() const {
    return (_flags >> FLUX_NSAMPLE_OFF) & NSAMPLE_MASK;
}

/** Sets the number of samples (sources) used to determine the PSF flux 
  * sample mean.
  */
void PerFilterSourceClusterAttributes::setNumFluxSamples(int samples)
{
    if (samples < 0 || samples > NSAMPLE_MASK) {
        throw LSST_EXCEPT(except::InvalidParameterException, "number of "
                          "flux samples (sources) is negative or too large");
    }
    _flags = (_flags & ~(NSAMPLE_MASK << FLUX_NSAMPLE_OFF)) |
             (samples << FLUX_NSAMPLE_OFF);
}

/** Sets the PSF flux and uncertainty.
  */
void PerFilterSourceClusterAttributes::setFlux(
    NullOr<float> const & flux,
    NullOr<float> const & fluxSigma)
{
    if (flux.isNull() && !fluxSigma.isNull()) {
        throw LSST_EXCEPT(except::InvalidParameterException,
                          "flux is null but uncertainty is not");
    }
    if (fluxSigma < 0.0) {
        throw LSST_EXCEPT(except::InvalidParameterException,
                          "negative flux uncertainty");
    }
    _flux = flux;
    _fluxSigma = fluxSigma;
}

/** Returns the number of smaples (sources) used to determine the ellipticity
  * parameter sample means.
  *
  * @li If this number is zero, then none of the sources satisified the
  *     the criteria for being included in the ellipticity parameter sample
  *     means, and the ellipticity parameters and uncertainties are invalid
  *     (NULL/NaN).
  */
int PerFilterSourceClusterAttributes::getNumEllipticitySamples() const {
    return (_flags >> ELLIPTICITY_NSAMPLE_OFF) & NSAMPLE_MASK;
}

/** Sets the number of samples (sources) used to determine the
  * ellipticity parameter sample means.
  */
void PerFilterSourceClusterAttributes::setNumEllipticitySamples(int samples)
{
    if (samples < 0 || samples > NSAMPLE_MASK) {
        throw LSST_EXCEPT(except::InvalidParameterException,
                          "number of ellipticity parameter samples (sources) "
                          "is negative or too large");
    }
    _flags = (_flags & ~(NSAMPLE_MASK << ELLIPTICITY_NSAMPLE_OFF)) |
             (samples << ELLIPTICITY_NSAMPLE_OFF);
}

/** Sets all ellipticity parameters and uncertainties to null.
  */
void PerFilterSourceClusterAttributes::setEllipticity()
{
    _e1.setNull();
    _e2.setNull();
    _radius.setNull();
    _e1Sigma.setNull();
    _e2Sigma.setNull();
    _radiusSigma.setNull();
}

/** Sets the ellipticity parameters and uncertainties. 
  */
void PerFilterSourceClusterAttributes::setEllipticity(
    NullOr<float> const & e1,
    NullOr<float> const & e2,
    NullOr<float> const & radius,
    NullOr<float> const & e1Sigma,
    NullOr<float> const & e2Sigma,
    NullOr<float> const & radiusSigma)
{
    if (e1.isNull() != e2.isNull() || e2.isNull() != radius.isNull()) {
        throw LSST_EXCEPT(except::InvalidParameterException, "ellipticity "
                          "parameters e1, e2, radius must either all be null "
                          "or all be valid");
    }
    if (e1Sigma.isNull() != e2Sigma.isNull() ||
        e2Sigma.isNull() != radiusSigma.isNull()) {
        throw LSST_EXCEPT(except::InvalidParameterException, "ellipticity "
                          "uncertainties for e1, e2, radius must either all "
                          "be null or all be valid");
    }
    if (e1.isNull() && !e1Sigma.isNull()) {
        throw LSST_EXCEPT(except::InvalidParameterException, "ellipticity "
                          "parameters are null, but their uncertainties are "
                          "not");
    }
    if (!e1Sigma.isNull() &&
        (e1Sigma < 0.0 || e2Sigma < 0.0 || radiusSigma < 0.0)) {
        throw LSST_EXCEPT(except::InvalidParameterException,
                          "negative ellipticity parameter uncertainty");
    }
    _e1 = e1;
    _e2 = e2;
    _radius = radius;
    _e1Sigma = e1Sigma;
    _e2Sigma = e2Sigma;
    _radiusSigma = radiusSigma;
}

bool PerFilterSourceClusterAttributes::operator==(
    PerFilterSourceClusterAttributes const & attributes) const
{
    return (_filterId == attributes._filterId &&
            _numObs == attributes._numObs &&
            _flags == attributes._flags &&
            _earliestObsTime == attributes._earliestObsTime &&
            _latestObsTime == attributes._latestObsTime &&
            _flux == attributes._flux &&
            _fluxSigma == attributes._fluxSigma &&
            _e1 == attributes._e1 &&
            _e2 == attributes._e2 &&
            _radius == attributes._radius &&
            _e1Sigma == attributes._e1Sigma &&
            _e2Sigma == attributes._e2Sigma &&
            _radiusSigma == attributes._radiusSigma);
}


// -- SourceClusterAttributes ----

SourceClusterAttributes::SourceClusterAttributes() :
    _clusterId(0),
    _numObs(0),
    _flags(0),
    _earliestObsTime(0.0), _latestObsTime(0.0),
    _raPs(0.0), _decPs(0.0),
    _raPsSigma(), _decPsSigma(), _raDecPsCov(),
    _raSg(0.0), _decSg(0.0),
    _raSgSigma(), _decSgSigma(), _raDecSgCov(),
    _perFilterAttributes()
{ }

SourceClusterAttributes::SourceClusterAttributes(int64_t id) :
    _clusterId(id),
    _numObs(0),
    _flags(0),
    _earliestObsTime(0.0), _latestObsTime(0.0),
    _raPs(0.0), _decPs(0.0),
    _raPsSigma(), _decPsSigma(), _raDecPsCov(),
    _raSg(0.0), _decSg(0.0),
    _raSgSigma(), _decSgSigma(), _raDecSgCov(),
    _perFilterAttributes()
{ }

SourceClusterAttributes::~SourceClusterAttributes() { }

bool SourceClusterAttributes::operator==(
    SourceClusterAttributes const & attributes) const
{
    if (_clusterId != attributes._clusterId ||
        _numObs != attributes._numObs ||
        _flags != attributes._flags ||
        _earliestObsTime != attributes._earliestObsTime ||
        _latestObsTime != attributes._latestObsTime ||
        _raPs != attributes._raPs ||
        _decPs != attributes._decPs ||
        _raPsSigma != attributes._raPsSigma ||
        _decPsSigma != attributes._decPsSigma ||
        _raDecPsCov != attributes._raDecPsCov ||
        _raSg != attributes._raSg ||
        _decSg != attributes._decSg ||
        _raSgSigma != attributes._raSgSigma ||
        _decSgSigma != attributes._decSgSigma ||
        _raDecSgCov != attributes._raDecSgCov) {
        return false;
    }
    return _perFilterAttributes == attributes._perFilterAttributes;
}

/** Returns @c true if the cluster has attributes specific to the given filter.
  */
bool SourceClusterAttributes::hasFilter(int filterId) const {
    return _perFilterAttributes.find(filterId) != _perFilterAttributes.end();
}

/** Returns a vector of the filter ids the cluster has attributes for.
  */
std::vector<int> const SourceClusterAttributes::getFilterIds() const {
    typedef PerFilterAttributesMap::const_iterator Iter;
    std::vector<int> filterIds;
    for (Iter i = _perFilterAttributes.begin(), e = _perFilterAttributes.end();
         i != e; ++i) {
        filterIds.push_back(i->first);
    }
    return filterIds;
}

/** Returns the filter specific cluster attributes for the given filter,
  * or throws an exception if none are available.
  */
PerFilterSourceClusterAttributes const &
SourceClusterAttributes::getPerFilterAttributes(int filterId) const
{
    typedef PerFilterAttributesMap::const_iterator Iter;
    Iter i = _perFilterAttributes.find(filterId);
    if (i == _perFilterAttributes.end()) {
        throw LSST_EXCEPT(except::NotFoundException, "source cluster "
                          "contains no attributes for specified filter");
    }
    return i->second;
}

/** Sets the earliest and latest cluster observation times.
  */
void SourceClusterAttributes::setObsTimeRange(double earliest, double latest)
{
    if (earliest > latest) {
        throw LSST_EXCEPT(except::InvalidParameterException, "earliest "
                          "observation time is after latest observation time");
    }
    _earliestObsTime = earliest;
    _latestObsTime = latest;
}

/** Sets the point source model cluster position and error.
  */
void SourceClusterAttributes::setPsPosition(double ra,
                                            double dec,
                                            NullOr<float> const & raSigma,
                                            NullOr<float> const & decSigma,
                                            NullOr<float> const & raDecCov)
{
    if (lsst::utils::isnan(ra) || lsst::utils::isnan(dec)) {
        throw LSST_EXCEPT(except::InvalidParameterException,
                          "Right ascension or declination is NaN");
    }
    if (raSigma.isNull() != decSigma.isNull()) {
        throw LSST_EXCEPT(except::InvalidParameterException,
                          "Right ascension and declination errors"
                          "must both be null or both be valid");
    }
    if (!raSigma.isNull() && (raSigma < 0.0 || decSigma < 0.0)) {
        throw LSST_EXCEPT(except::InvalidParameterException,
                          "negative position error");
    }
    if (raSigma.isNull() && !raDecCov.isNull()) {
        throw LSST_EXCEPT(except::InvalidParameterException,
                          "Position errors are null "
                          "but covariance is not.");
    }
    _raPs = ra;
    _decPs = dec;
    _raPsSigma = raSigma;
    _decPsSigma = decSigma;
    _raDecPsCov = raDecCov;
}

/** Sets the small galaxy model cluster position and error.
  */
void SourceClusterAttributes::setSgPosition(NullOr<double> const & ra,
                                            NullOr<double> const & dec,
                                            NullOr<float> const & raSigma,
                                            NullOr<float> const & decSigma,
                                            NullOr<float> const & raDecCov)
{
    if (ra.isNull() != dec.isNull()) {
        throw LSST_EXCEPT(except::InvalidParameterException,
                          "Right ascension and declination must both "
                          "be null or both be valid");
    }
    if (raSigma.isNull() != decSigma.isNull()) {
        throw LSST_EXCEPT(except::InvalidParameterException,
                          "Right ascension and declination errors "
                          "must both be null or both be valid");
    }
    if (!raSigma.isNull() && (raSigma < 0.0 || decSigma < 0.0)) {
        throw LSST_EXCEPT(except::InvalidParameterException,
                          "negative position error");
    }
    if (ra.isNull() && (!raSigma.isNull() || !raDecCov.isNull())) {
        throw LSST_EXCEPT(except::InvalidParameterException,
                          "Position is NULL but errors are not");
    }
    if (raSigma.isNull() && !raDecCov.isNull()) {
        throw LSST_EXCEPT(except::InvalidParameterException,
                          "Position errors are null "
                          "but covariance is not.");
    }
    _raSg = ra;
    _decSg = dec;
    _raSgSigma = raSigma;
    _decSgSigma = decSigma;
    _raDecSgCov = raDecCov;
}

/** Removes all filter specific attributes from the source cluster.
  */
void SourceClusterAttributes::clearPerFilterAttributes() {
    _perFilterAttributes.clear();
}

/** Adds or replaces filter specific attributes for the source cluster.
  *
  * @return @c true if the attributes were added, @c false if existing
  *         attributes were replaced.
  */
bool SourceClusterAttributes::setPerFilterAttributes(
    PerFilterSourceClusterAttributes const & attributes)
{
    typedef PerFilterAttributesMap::iterator Iter;
    int filterId = attributes.getFilterId();
    std::pair<Iter, bool> v = _perFilterAttributes.insert(
        std::make_pair(filterId, attributes));
    if (!v.second) {
        v.first->second = attributes;
    }
    return v.second;
}

/** Removes cluster attributes specific to the given filter.
  *
  * @return  @c false if no attributes for the given filter were found.
  */
bool SourceClusterAttributes::removePerFilterAttributes(int filterId) {
    return _perFilterAttributes.erase(filterId) != 0;
}

/** Computes cluster attributes from the given source.
  */
void SourceClusterAttributes::computeAttributes(
    boost::shared_ptr<lsst::afw::detection::Source const> source,
    lsst::ap::match::ExposureInfoMap const & exposures,
    double fluxScale,         ///< Flux scaling factor
    int fluxIgnoreMask,       ///< Detection flag bitmask identifying sources
                              ///  to ignore when computing fluxes
    int ellipticityIgnoreMask ///< Detection flag bitmask identifying sources
                              ///  to ignore when computing ellipticities
) {
    setNumObs(1);
    setObsTimeRange(source->getTaiMidPoint(), source->getTaiMidPoint());
    match::ExposureInfo::ConstPtr exposure =
        exposures.get(source->getAmpExposureId());
    if (!exposure) {
        throw LSST_EXCEPT(except::NotFoundException,
                          "No ExposureInfo found for source");
    }
    image::Wcs::ConstPtr wcs = exposure->getWcs();
    // get a Coord in the same coordinate system as the input WCSes
    coord::Coord::Ptr center = wcs->getSkyOrigin()->clone();
    double ra = source->getRa();
    double dec = source->getDec();
    center->reset(degrees(ra), degrees(dec), center->getEpoch());
    NeTanProj proj(center);

    if (source->isNull(detection::RA_ASTROM_ERR) ||
        source->isNull(detection::DEC_ASTROM_ERR) ||
        lsst::utils::isnan(ra) ||
        source->getRaAstromErr() < 0.0f ||
        lsst::utils::isnan(dec) ||
        source->getDecAstromErr() < 0.0f) {
        setPsPosition(ra, dec, NullOr<float>(), NullOr<float>(), NullOr<float>());
        setSgPosition(ra, dec, NullOr<float>(), NullOr<float>(), NullOr<float>());
    } else {
        setPsPosition(ra, dec, 
                      source->getRaAstromErr(), source->getDecAstromErr(),
                      NullOr<float>());
        setSgPosition(ra, dec, 
                      source->getRaAstromErr(), source->getDecAstromErr(),
                      NullOr<float>());
    }
    SourceAndExposure se(source, exposure, proj.pixelToNeTransform(
        Eigen::Vector2d(source->getXAstrom(), source->getYAstrom()), *wcs));
    PerFilterSourceClusterAttributes pfa;
    pfa.setFilterId(source->getFilterId());
    pfa.setNumObs(1);
    pfa.setObsTimeRange(source->getTaiMidPoint(), source->getTaiMidPoint());
    pfa.computeFlux(se, fluxScale, fluxIgnoreMask);
    pfa.computeEllipticity(se, ellipticityIgnoreMask);
    _perFilterAttributes.insert(std::make_pair(source->getFilterId(), pfa));
}

/** Computes cluster attributes from a set of sources.
  */
void SourceClusterAttributes::computeAttributes(
    lsst::afw::detection::SourceSet const & sources,
    lsst::ap::match::ExposureInfoMap const & exposures,
    double fluxScale,
    int fluxIgnoreMask,
    int ellipticityIgnoreMask
) {
    typedef detection::SourceSet::const_iterator SourceIter;
    typedef std::vector<SourceAndExposure>::const_iterator SeIter;
    typedef std::tr1::unordered_map<int, std::vector<SourceAndExposure> >
        FilterMap;

    if (sources.empty()) {
        throw LSST_EXCEPT(except::InvalidParameterException, "Cannot compute "
                          "cluster attributes from an empty SourceSet");
    } else if (sources.size() == 1) {
        computeAttributes(sources.front(), exposures, fluxScale,
                          fluxIgnoreMask, ellipticityIgnoreMask);
        return;
    }
 
    setNumObs(static_cast<int>(sources.size()));
    {
        double tbeg = sources.front()->getTaiMidPoint();
        double tend = tbeg;
        for (SourceIter i = sources.begin() + 1, e = sources.end(); i != e; ++i) {
            double t = (*i)->getTaiMidPoint();
            if (t < tbeg) {
                tbeg = t;
            } else if (t > tend) {
                tend = t;
            }
        }
        setObsTimeRange(tbeg, tend);
    }
    // compute fiducial cluster position (unweighted mean,
    // the point source model position).
    _computePsPosition(sources);
    // get a Coord in the same coordinate system as the input WCSes
    // and equal to the fiducial cluster position
    coord::Coord::Ptr center;
    {
        match::ExposureInfo::ConstPtr exposure = exposures.get(
            sources.front()->getAmpExposureId());
        if (!exposure) {
            throw LSST_EXCEPT(except::NotFoundException,
                              "No ExposureInfo found for source");
        }
        center = exposure->getWcs()->getSkyOrigin()->clone();
        center->reset(degrees(_raPs), degrees(_decPs), center->getEpoch());
    }
    // construct N,E tan projection at fiducial cluster position
    NeTanProj proj(center);
    // build a vector of SourceAndExposure objects
    std::vector<SourceAndExposure> se;
    se.reserve(sources.size());
    for (SourceIter i = sources.begin(), e = sources.end(); i != e; ++i) {
        match::ExposureInfo::ConstPtr exposure = exposures.get(
            (*i)->getAmpExposureId());
        if (!exposure) {
            throw LSST_EXCEPT(except::NotFoundException,
                              "No ExposureInfo found for source");
        }
        se.push_back(SourceAndExposure(*i, exposure, proj.pixelToNeTransform(
            Eigen::Vector2d((*i)->getXAstrom(), (*i)->getYAstrom()),
            *exposure->getWcs())));
    }
    // compute small galaxy model position
    Eigen::Vector2d sgPos;
    Eigen::Matrix2d sgCov;
    if (combinePositions(se, proj, sgPos, sgCov)) {
        setSgPosition(sgPos.x(), sgPos.y(),
                      sqrt(sgCov(0,0)),
                      sqrt(sgCov(1,1)),
                      sgCov(0,1));
    }
    // bin sources by filter - does not rely on filterId being in range [0, N)
    FilterMap filters;
    for (SeIter i = se.begin(), e = se.end(); i != e; ++i) {
        int filterId = i->getSource()->getFilterId();
        FilterMap::iterator j = filters.find(filterId);
        if (j == filters.end()) {
            j = filters.insert(std::make_pair(
                filterId, std::vector<SourceAndExposure>())).first;
        }
        j->second.push_back(*i);
    }
    // compute per-filter properties
    for (FilterMap::iterator i = filters.begin(), e = filters.end();
         i != e; ++i) {
        PerFilterSourceClusterAttributes pfa;
        pfa.setFilterId(i->first);
        pfa.setNumObs(static_cast<int>(i->second.size()));
        double tbeg = i->second.front().getSource()->getTaiMidPoint();
        double tend = tbeg;
        for (SeIter ii = i->second.begin() + 1, ee = i->second.end();
             ii != ee; ++ii) {
            double t = ii->getSource()->getTaiMidPoint();
            if (t < tbeg) {
                tbeg = t;
            } else if (t > tend) {
                tend = t;
            }
        }
        pfa.setObsTimeRange(tbeg, tend);
        pfa.computeFlux(i->second, fluxScale, fluxIgnoreMask);
        pfa.computeEllipticity(i->second, ellipticityIgnoreMask);
        _perFilterAttributes.insert(std::make_pair(i->first, pfa));
    }
}

/** Computes the point source model position of the cluster by finding the
  * position P minimizing the sum of the squared angular separations between
  * P and each source position. This an unweighted mean position.
  */
void SourceClusterAttributes::_computePsPosition(
    lsst::afw::detection::SourceSet const & sources)
{
    typedef detection::SourceSet::const_iterator Iter;
    size_t const n = sources.size();
    if (n == 1) {
        detection::Source const & source = *sources.front();
        if (source.isNull(detection::RA_ASTROM_ERR) ||
            source.isNull(detection::DEC_ASTROM_ERR) ||
            lsst::utils::isnan(source.getRaAstromErr()) ||
            source.getRaAstromErr() <= 0.0f ||
            lsst::utils::isnan(source.getDecAstromErr()) ||
            source.getDecAstromErr() <= 0.0f) {
            setPsPosition(source.getRa(),
                          source.getDec(),
                          NullOr<float>(),
                          NullOr<float>(),
                          NullOr<float>());
        } else {
            setPsPosition(source.getRa(),
                          source.getDec(),
                          source.getRaAstromErr(),
                          source.getDecAstromErr(),
                          NullOr<float>());
        }
        return;
    }
    // compute point P such that the sum of the squared angular separations
    // between P and each source position is minimized.
    Eigen::Vector3d v = Eigen::Vector3d::Zero();
    for (Iter i = sources.begin(), e = sources.end(); i != e; ++i) {
        v += sphericalToCartesian((*i)->getRa(), (*i)->getDec());
    }
    Eigen::Vector2d sc = cartesianToSpherical(v);
    // compute covariance matrix
    double covRaRa = 0.0;
    double covDecDec = 0.0;
    double covRaDec = 0.0;
    for (Iter i = sources.begin(), e = sources.end(); i != e; ++i) {
       // use minimum delta ra to avoid huge errors for samples spanning
       // the 0/2*M_PI discontinuity
       double dr = std::fabs((*i)->getRa() - sc.x());
       double dd = (*i)->getDec() - sc.y();
       if (dr > M_PI) {
          dr = 2*M_PI - dr;
       } 
       covRaRa += dr * dr;
       covDecDec += dd * dd;
       covRaDec += dr * dd;
    }
    // set the ra/dec uncertainties to an estimate of the
    // standard deviation of the sample mean
    setPsPosition(sc.x(), sc.y(),
                  static_cast<float>(sqrt(covRaRa / (n * (n - 1)))),
                  static_cast<float>(sqrt(covDecDec / (n * (n - 1)))),
                  static_cast<float>(covRaDec / n));
}

}}} // namespace lsst:ap::cluster

