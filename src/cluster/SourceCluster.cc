// -*- lsst-c++ -*-
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

#include "boost/scoped_array.hpp"

#include "lsst/pex/exceptions.h"
#include "lsst/pex/policy.h"
#include "lsst/ap/Common.h"
#include "lsst/ap/cluster/optics/Metrics.h"
#include "lsst/ap/cluster/optics/Optics.cc"


namespace detection = lsst::afw::detection;
namespace except = lsst::pex::exceptions;
namespace policy = lsst::pex::policy;

using std::sqrt;


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

} // namespace


/** Clusters a set of sources using the OPTICS algorithm. The following
  * parameters are read from @c policy:
  *
  * @li @c "epsilon" (double) : generating distance for clusters (arcsec).
  * @li @c "minPoints" (int) : minimum number of points that must be in an
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
                      policy->getInt("minPoints"),
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
        (*i)->setRaObject(cluster.getRa());
        (*i)->setDecObject(cluster.getDec());
    }
}

/** Removes sources without positions (those with longitude and/or latitude
  * angles that are NaN or out of bounds) from @c sources and appends them
  * to @c badSources. Valid longitude angles are in range [0, 2*pi) and valid
  * latitude angles are in range [-pi/2, pi/2].
  */
LSST_AP_API void segregateInvalidSources(
    lsst::afw::detection::SourceSet & sources,
    lsst::afw::detection::SourceSet & badSources)
{
    size_t const n = sources.size();
    size_t j = 0;
    for (size_t i = 0; i < n; ++i) {
        double ra = sources[i]->getRa();
        double dec = sources[i]->getDec();
        if (lsst::utils::isnan(ra) || lsst::utils::isnan(dec) ||
            ra < 0.0 || ra >= 2.0 * M_PI ||
            dec < -0.5 * M_PI || dec > 0.5 * M_PI) {
            badSources.push_back(sources[i]);
        } else {
            if (j != i) {
                sources[j] = sources[i];
            }
            ++j;
        }
    }
    if (j < n) {
        sources.erase(sources.begin() + j, sources.end());
    }
}

/** Removes "bad" sources from @c sources and appends them to @c badSources.
  * A source is considered bad when one or more of the bits set in its detection
  * flags matches a bit set in @c badSourceMask.
  */
LSST_AP_API void segregateBadSources(
    lsst::afw::detection::SourceSet & sources,
    lsst::afw::detection::SourceSet & badSources,
    int const badSourceMask)
{
    size_t const n = sources.size();
    size_t j = 0;
    for (size_t i = 0; i < n; ++i) {
        if ((sources[i]->getFlagForDetection() & badSourceMask) != 0) {
            badSources.push_back(sources[i]);
        } else {
            if (j != i) {
                sources[j] = sources[i];
            }
            ++j;
        }
    }
    if (j < n) {
        sources.erase(sources.begin() + j, sources.end());
    }
}

/** Sets the object id of each of the bad sources to NULL. Also sets
  * the object position of each of bad source to the source position.
  * This allows bad sources to be stored in the same table as good
  * sources (which are partitioned by the position of their associated
  * objects/clusters).
  */
LSST_AP_API void updateBadSources(lsst::afw::detection::SourceSet & badSources)
{
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

/** Creates a new PerFilterSourceClusterAttributes, computing attributes
  * from the given source.
  *
  * @param[in] source                   Source to compute attributes from.
  * @param[in] fluxIgnoreMask           Detection flag bitmask identifying
  *                                     sources that should be ignored when
  *                                     determining cluster fluxes.
  * @param[in] ellipticityIgnoreMask    Detection flag bitmask identifying
  *                                     sources that should be ignore when
  *                                     determining cluster ellipticities 

  */
PerFilterSourceClusterAttributes::PerFilterSourceClusterAttributes(
    lsst::afw::detection::Source const & source,
    int fluxIgnoreMask,
    int ellipticityIgnoreMask
) :
    _filterId(source.getFilterId()),
    _numObs(1),
    _flags(0),
    _earliestObsTime(source.getTaiMidPoint()),
    _latestObsTime(source.getTaiMidPoint()),
    _flux(), _fluxSigma(),
    _e1(), _e2(), _radius(),
    _e1Sigma(), _e2Sigma(), _radiusSigma()
{
    if (!lsst::utils::isnan(source.getPsfFlux()) &&
        (source.getFlagForDetection() & fluxIgnoreMask) == 0) {
        setFlux(source.getPsfFlux(), source.getPsfFluxErr());
        setNumFluxSamples(1);
    }
    if (source.isNull(detection::IXX) ||
        source.isNull(detection::IYY) ||
        source.isNull(detection::IXY) ||
        lsst::utils::isnan(source.getIxx()) ||
        lsst::utils::isnan(source.getIyy()) ||
        lsst::utils::isnan(source.getIxy()) ||
        (source.getFlagForDetection() & ellipticityIgnoreMask) != 0) {
        return;
    }
    double mxx = source.getIxx();
    double myy = source.getIyy();
    double mxy = source.getIxy();
    double t = mxx + myy;
    if (t != 0.0) {
        setNumEllipticitySamples(1);
        setEllipticity(static_cast<float>((mxx - myy) / t),
                       static_cast<float>(2.0 * mxy / t),
                       static_cast<float>(sqrt(t)),
                       NullOr<float>(),
                       NullOr<float>(),
                       NullOr<float>());
    }
}

/** Creates a new PerFilterSourceClusterAttributes, computing attributes
  * from the given set of sources.
  *
  * @param[in] sources                  Sources to compute attributes from.
  * @param[in] fluxIgnoreMask           Detection flag bitmask identifying
  *                                     sources that should be ignored when
  *                                     determining cluster fluxes.
  * @param[in] ellipticityIgnoreMask    Detection flag bitmask identifying
  *                                     sources that should be ignore when
  *                                     determining cluster ellipticities 

  */
PerFilterSourceClusterAttributes::PerFilterSourceClusterAttributes(
    lsst::afw::detection::SourceSet const & sources,
    int fluxIgnoreMask,
    int ellipticityIgnoreMask
) :
    _numObs(static_cast<int>(sources.size())),
    _flags(0),
    _flux(), _fluxSigma(),
    _e1(), _e2(), _radius(),
    _e1Sigma(), _e2Sigma(), _radiusSigma()
{
    typedef detection::SourceSet::const_iterator Iter;
    if (sources.size() == 0) {
        throw LSST_EXCEPT(except::InvalidParameterException,
                          "empty source set");
    }
    _filterId = sources.front()->getFilterId();
    double tbeg = std::numeric_limits<double>::infinity();
    double tend = -tbeg;
    for (Iter i = sources.begin(), e = sources.end(); i != e; ++i) {
        if ((*i)->getFilterId() != _filterId) {
            throw LSST_EXCEPT(except::InvalidParameterException,
                              "sources originate from multiple filters");
        }
        double t = (*i)->getTaiMidPoint();
        tbeg = std::min(tbeg, t);
        tend = std::max(tend, t);
    }
    _earliestObsTime = tbeg;
    _latestObsTime = tend;
    computeFlux(sources, fluxIgnoreMask);
    computeEllipticity(sources, ellipticityIgnoreMask);
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
  * @li If this number is one, then the ellipticity parameter uncertainties
  *     are invalid (NULL/NaN). In principle, uncertainties could be derived
  *     from the moment covariance matrix of the source, but sources do not
  *     carry enough of this matrix for this to be possible.
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

/** Computes the sample mean and standard error of the fluxes of each
  * source from this filter.
  */
void PerFilterSourceClusterAttributes::computeFlux(
    lsst::afw::detection::SourceSet const & sources,
    int const fluxIgnoreMask)
{
    typedef detection::SourceSet::const_iterator Iter;
    if (sources.empty()) {
        throw LSST_EXCEPT(except::InvalidParameterException,
                          "empty source set");
    }
    int ns = 0;
    double flux = 0.0;
    for (Iter i = sources.begin(), e = sources.end(); i != e; ++i) {
        double f = (*i)->getPsfFlux();
        if (lsst::utils::isnan(f) ||
            ((*i)->getFlagForDetection() & fluxIgnoreMask) != 0) {
            continue;
        }
        flux += f;
        ++ns;
    }
    setNumFluxSamples(ns);
    if (ns == 1) {
        // set uncertainty to PSF flux uncertainty of the only source
        // available
        for (Iter i = sources.begin(), e = sources.end(); i != e; ++i) {
            float f = (*i)->getPsfFlux();
            if (!lsst::utils::isnan(f) && 
                ((*i)->getFlagForDetection() & fluxIgnoreMask) == 0) {
                setFlux(f, (*i)->getPsfFluxErr());
                break;
            }
        }
    } else if (ns > 1) {
        // set flux to the sample mean
        flux /= ns;
        // set flux uncertainty to an estimate of the standard
        // deviation of the sample mean 
        double ff = 0.0;
        for (Iter i = sources.begin(), e = sources.end(); i != e; ++i) {
            double f = (*i)->getPsfFlux();
            if (lsst::utils::isnan(f) ||
                ((*i)->getFlagForDetection() & fluxIgnoreMask) != 0) {
                continue;
            }
            double df = f - flux;
            ff += df * df;
        }
        setFlux(static_cast<float>(flux),
                static_cast<float>(sqrt(ff / (ns * (ns - 1)))));
    }
}

/** Computes the sample means and standard errors of the ellipticity
  * parameters of each source from this filter. Ellipticity parameters
  * are derived from the sources adaptive moments - see 
  *
  * Image Ellipticity From Atmospheric Aberrations
  * W. H. de Vries, S. S. Olivier, S. J. Asztalos, L. J. Rosenberg, and K. L. Baker
  * The Astrophysical Journal, 622:744-749, 2007 June 10
  *
  * and the implementation of @c lsst::meas::algorithms::Shape .
  */
void PerFilterSourceClusterAttributes::computeEllipticity(
    lsst::afw::detection::SourceSet const & sources,
    int const ellipticityIgnoreMask)
{
    typedef detection::SourceSet::const_iterator Iter;
    if (sources.empty()) {
        throw LSST_EXCEPT(except::InvalidParameterException,
                          "empty source set");
    }
    boost::scoped_array<Eigen::Vector3d> eparams(
        new Eigen::Vector3d[sources.size()]);
    Eigen::Vector3d m(0.0, 0.0, 0.0);
    int ns = 0;
    for (Iter i = sources.begin(), e = sources.end(); i != e; ++i) {
       if ((*i)->isNull(detection::IXX) ||
           (*i)->isNull(detection::IYY) ||
           (*i)->isNull(detection::IXY) ||
           ((*i)->getFlagForDetection() & ellipticityIgnoreMask) != 0) {
           continue;
       }
       double mxx = (*i)->getIxx();
       double myy = (*i)->getIyy();
       double mxy = (*i)->getIxy();
       double t = mxx + myy;
       // make sure the moments aren't NaN
       if (lsst::utils::isnan(mxx) || lsst::utils::isnan(myy) ||
           lsst::utils::isnan(mxy) || (t == 0.0)) {
           continue;
       }
       // compute ellipticity parameters from moments
       Eigen::Vector3d ep((mxx - myy) / t, 2.0 * mxy / t, sqrt(t));
       // store for later
       eparams[ns++] = ep;
       // and add to running sum
       m += ep;
    }
    setNumEllipticitySamples(ns);
    if (ns > 0) {
        // set cluster ellipticity parameters to sample means
        m /= static_cast<double>(ns);
        if (ns == 1) {
            // Ellipticity uncertainties could be estimated from the covariance
            // matrix of the moments. Unfortunately, only a subset of this
            // matrix is stored in Source, so unless ns > 1, the ellipticity
            // uncertainties are NULL/NaN.
            setEllipticity(static_cast<float>(m(0)),
                           static_cast<float>(m(1)),
                           static_cast<float>(m(2)),
                           NullOr<float>(),
                           NullOr<float>(),
                           NullOr<float>());
        } else {
            // compute standard errors
            Eigen::Vector3d v(0.0, 0.0, 0.0);
            for (int i = 0; i < ns; ++i) {
                v += (eparams[i] - m).cwise() * (eparams[i] - m);
            }
            setEllipticity(static_cast<float>(m(0)), 
                           static_cast<float>(m(1)), 
                           static_cast<float>(m(2)),
                           static_cast<float>(sqrt(v[0] / (ns * (ns - 1)))),
                           static_cast<float>(sqrt(v[1] / (ns * (ns - 1)))),
                           static_cast<float>(sqrt(v[2] / (ns * (ns - 1)))));
        }
    }
}


// -- SourceClusterAttributes ----

SourceClusterAttributes::SourceClusterAttributes() :
    _clusterId(0),
    _numObs(0),
    _flags(0),
    _earliestObsTime(0.0), _latestObsTime(0.0),
    _ra(0.0), _dec(0.0),
    _raSigma(0.0), _decSigma(0.0),
    _raDecCov(),
    _perFilterAttributes()
{ }

/** Creates a new SourceClusterAttributes with the given id, computing
  * attributes from the given source.
  *
  * @param[in] source   Source to compute attributes from.
  * @param[in] id       Source cluster id
  * @param[in] fluxIgnoreMask           Detection flag bitmask identifying
  *                                     sources that should be ignored when
  *                                     determining cluster fluxes.
  * @param[in] ellipticityIgnoreMask    Detection flag bitmask identifying
  *                                     sources that should be ignore when
  *                                     determining cluster ellipticity
  *                                     parameters. 
  */
SourceClusterAttributes::SourceClusterAttributes(
    lsst::afw::detection::Source const & source,
    int64_t id,
    int fluxIgnoreMask,
    int ellipticityIgnoreMask
) :
    _clusterId(id),
    _numObs(1),
    _flags(0),
    _earliestObsTime(source.getTaiMidPoint()),
    _latestObsTime(source.getTaiMidPoint()),
    _raDecCov(),
    _perFilterAttributes()
{
    if (source.isNull(detection::RA_ASTROM_ERR) ||
        source.isNull(detection::DEC_ASTROM_ERR) ||
        lsst::utils::isnan(source.getRaAstromErr()) ||
        source.getRaAstromErr() < 0.0f ||
        lsst::utils::isnan(source.getDecAstromErr()) ||
        source.getDecAstromErr() < 0.0f) {
        setPosition(source.getRa(),
                    source.getDec(),
                    NullOr<float>(),
                    NullOr<float>(),
                    NullOr<float>());
    } else {
        setPosition(source.getRa(),
                    source.getDec(), 
                    source.getRaAstromErr(),
                    source.getDecAstromErr(),
                    NullOr<float>());
    }
    _perFilterAttributes.insert(std::make_pair(source.getFilterId(),
            PerFilterSourceClusterAttributes(
                source, fluxIgnoreMask, ellipticityIgnoreMask)));
}

/** Creates a new SourceClusterAttributes with the given id, computing
  * attributes from the given set of sources.
  *
  * @param[in] sources  Sources to compute attributes from.
  * @param[in] id       Source cluster id
  * @param[in] fluxIgnoreMask           Detection flag bitmask identifying
  *                                     sources that should be ignored when
  *                                     determining cluster fluxes.
  * @param[in] ellipticityIgnoreMask    Detection flag bitmask identifying
  *                                     sources that should be ignore when
  *                                     determining cluster ellipticity
  *                                     parameters. 
  */
SourceClusterAttributes::SourceClusterAttributes(
    lsst::afw::detection::SourceSet const & sources,
    int64_t id,
    int fluxIgnoreMask,
    int ellipticityIgnoreMask
) :
    _clusterId(id),
    _numObs(static_cast<int>(sources.size())),
    _flags(0),
    _raDecCov(),
    _perFilterAttributes()
{
    typedef detection::SourceSet::const_iterator SourceIter;
    typedef std::tr1::unordered_map<int, detection::SourceSet>::iterator
        HashIter;
    if (sources.empty()) {
        throw LSST_EXCEPT(except::InvalidParameterException,
                          "empty source set");
    }
    double earliest = std::numeric_limits<double>::infinity();
    double latest = -earliest;
    for (SourceIter i = sources.begin(), e = sources.end(); i != e; ++i) {
        double t = (*i)->getTaiMidPoint();
        earliest = std::min(earliest, t);
        latest = std::max(latest, t);
    }
    _earliestObsTime = earliest;
    _latestObsTime = latest;
    computePosition(sources);
    // bin sources by filter
    std::tr1::unordered_map<int, detection::SourceSet> filters;
    for (SourceIter i = sources.begin(), e = sources.end(); i != e; ++i) {
        int filterId = (*i)->getFilterId();
        HashIter j = filters.find(filterId);
        if (j == filters.end()) {
            j = filters.insert(std::make_pair(
                filterId, detection::SourceSet())).first;
        }
        j->second.push_back(*i);
    }
    // compute per-filter properties
    for (HashIter i = filters.begin(), e = filters.end(); i != e; ++i) {
        _perFilterAttributes.insert(std::make_pair(i->first,
            PerFilterSourceClusterAttributes(
                i->second, fluxIgnoreMask, ellipticityIgnoreMask)));
    }
}

SourceClusterAttributes::~SourceClusterAttributes() { }

bool SourceClusterAttributes::operator==(
    SourceClusterAttributes const & attributes) const
{
    if (_clusterId != attributes._clusterId ||
        _numObs != attributes._numObs ||
        _flags != attributes._flags ||
        _earliestObsTime != attributes._earliestObsTime ||
        _latestObsTime != attributes._latestObsTime ||
        _ra != attributes._ra ||
        _dec != attributes._dec ||
        _raSigma != attributes._raSigma ||
        _decSigma != attributes._decSigma ||
        _raDecCov != attributes._raDecCov) {
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

/** Sets the cluster position and position uncertainties.
  */
void SourceClusterAttributes::setPosition(double ra,
                                          double dec,
                                          NullOr<float> const & raSigma,
                                          NullOr<float> const & decSigma,
                                          NullOr<float> const & raDecCov)
{
    if (lsst::utils::isnan(ra) || lsst::utils::isnan(dec)) {
        throw LSST_EXCEPT(except::InvalidParameterException,
                          "Longitude and/or latitude angle is NaN");
    }
    if (raSigma.isNull() != decSigma.isNull()) {
        throw LSST_EXCEPT(except::InvalidParameterException,
                          "Longitude/latitude angle uncertainties must both "
                          "be null, or both be valid");
    }
    if (!raSigma.isNull() && (raSigma < 0.0 || decSigma < 0.0)) {
        throw LSST_EXCEPT(except::InvalidParameterException,
                          "negative position uncertainty");
    }
    _ra = ra;
    _dec = dec;
    _raSigma = raSigma;
    _decSigma = decSigma;
    _raDecCov = raDecCov;
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

/** Computes the position of the cluster by finding the position P minimizing
  * the sum of the squared angular separations between P and each source
  * position.
  */
void SourceClusterAttributes::computePosition(
    lsst::afw::detection::SourceSet const & sources)
{
    typedef detection::SourceSet::const_iterator Iter;
    size_t const n = sources.size();
    if (n == 0) {
        throw LSST_EXCEPT(except::InvalidParameterException,
                          "empty source set");
    }
    if (n == 1) {
        detection::Source const & source = *sources.front();
        if (source.isNull(detection::RA_ASTROM_ERR) ||
            source.isNull(detection::DEC_ASTROM_ERR) ||
            lsst::utils::isnan(source.getRaAstromErr()) ||
            source.getRaAstromErr() < 0.0f ||
            lsst::utils::isnan(source.getDecAstromErr()) ||
            source.getDecAstromErr() < 0.0f) {
            setPosition(source.getRa(),
                        source.getDec(),
                        NullOr<float>(),
                        NullOr<float>(),
                        NullOr<float>());
        } else {
            setPosition(source.getRa(),
                        source.getDec(),
                        source.getRaAstromErr(),
                        source.getDecAstromErr(),
                        NullOr<float>());
        }
        return;
    }
    // Compute point P such that the sum of the squared angular separations
    // between P and each source position is minimized.
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    for (Iter i = sources.begin(), e = sources.end(); i != e; ++i) {
        double ra = (*i)->getRa();
        double dec = (*i)->getDec();
        double cosDec = std::cos(dec);
        x += std::cos(ra) * cosDec;
        y += std::sin(ra) * cosDec;
        z += std::sin(dec); 
    }
    double d2 = x * x + y * y; 
    double ra = (d2 == 0.0) ? 0.0 : atan2(y, x);
    double dec = (z == 0.0) ? 0.0 : atan2(z, std::sqrt(d2)); 
    // Compute covariance matrix
    double covRaRa = 0.0;
    double covDecDec = 0.0;
    double covRaDec = 0.0;
    for (Iter i = sources.begin(), e = sources.end(); i != e; ++i) {
       double dr = (*i)->getRa() - ra;
       double dd = (*i)->getDec() - dec;
       covRaRa += dr * dr;
       covDecDec += dd * dd;
       covRaDec += dr * dd;
    }
    // Set the ra/dec uncertainties to an estimate of the
    // standard deviation of the sample mean
    setPosition(ra,
                dec,
                static_cast<float>(std::sqrt(covRaRa / (n * (n - 1)))),
                static_cast<float>(std::sqrt(covDecDec / (n * (n - 1)))),
                static_cast<float>(covRaDec / n));
}

}}} // namespace lsst:ap::cluster

