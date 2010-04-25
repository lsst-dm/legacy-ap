// -*- lsst-c++ -*-
/** @file
  * @brief Implementation of the OPTICS algorithm.
  *
  * For details of the algorithm, see the following paper:
  *
  * "OPTICS: Ordering Points To Identify the Clustering Structure".
  * Mihael Ankerst, Markus M. Breunig, Hans-Peter Kriegel, Jörg Sander (1999).
  * ACM SIGMOD international conference on Management of data.
  * ACM Press. pp. 49–60.
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
#include "lsst/meas/algorithms/Measure.h"
#include "lsst/ap/Common.h"
#include "lsst/ap/cluster/optics/Metrics.h"
#include "lsst/ap/cluster/optics/Optics.cc"


namespace detection = lsst::afw::detection;
namespace except = lsst::pex::exceptions;
namespace policy = lsst::pex::policy;

using lsst::meas::algorithms::Flags;


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

/// Maximum number of sources that can be processed at once.
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


// -- PerFilterSourceClusterAttributes ----

PerFilterSourceClusterAttributes::PerFilterSourceClusterAttributes() :
    _filterId(0),
    _numObs(0),
    _earliestObsTime(0.0), _latestObsTime(0.0),
    _flux(0.0), _fluxSigma(0.0),
    _e1(), _e1Sigma(),
    _e2(), _e2Sigma(),
    _radius(), _radiusSigma()
{ }

/** Creates a new PerFilterSourceClusterAttributes, computing attributes
  * from the given set of sources.
  */
PerFilterSourceClusterAttributes::PerFilterSourceClusterAttributes(
    lsst::afw::detection::SourceSet const & sources) :
    _numObs(static_cast<int>(sources.size())),
    _e1(), _e1Sigma(),
    _e2(), _e2Sigma(),
    _radius(), _radiusSigma()
{
    typedef detection::SourceSet::const_iterator Iter;
    if (sources.size() == 0) {
        throw LSST_EXCEPT(except::InvalidParameterException,
                          "empty source set");
    }
    _filterId = sources.front()->getFilterId();
    _numObs = static_cast<int>(sources.size());
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
    computeFlux(sources);
    computeEllipticities(sources);
}

PerFilterSourceClusterAttributes::~PerFilterSourceClusterAttributes() { }

/** Sets the earliest and latest cluster observation times in this filter.
  */
void PerFilterSourceClusterAttributes::setObsTimeRange(
    double earliest, double latest)
{
    if (earliest > latest) {
        throw LSST_EXCEPT(except::InvalidParameterException, "earliest "
                          "observation time is after latest observation time");
    }
    _earliestObsTime = earliest;
    _latestObsTime = latest;
}

/** Sets the flux and flux uncertainty for the cluster in this filter.
  */
void PerFilterSourceClusterAttributes::setFlux(double flux, double fluxSigma)
{
    if (fluxSigma < 0.0) {
        throw LSST_EXCEPT(except::InvalidParameterException,
                          "negative flux sigma");
    }
    _flux = flux;
    _fluxSigma = fluxSigma;
}

/** Sets all ellipticity parameters and uncertainties to null.
  */
void PerFilterSourceClusterAttributes::setEllipticities()
{
    _e1.setNull();
    _e1Sigma.setNull();
    _e2.setNull();
    _e2Sigma.setNull();
    _radius.setNull();
    _radiusSigma.setNull();
}

/** Sets the ellipticity parameters and uncertainties for the
  * cluster in this filter.
  */
void PerFilterSourceClusterAttributes::setEllipticities(
    Nullable<double> const & e1,
    Nullable<double> const & e1Sigma,
    Nullable<double> const & e2,
    Nullable<double> const & e2Sigma,
    Nullable<double> const & radius,
    Nullable<double> const & radiusSigma)
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
        (e1Sigma() < 0.0 || e2Sigma() < 0.0 || radiusSigma() < 0.0)) {
        throw LSST_EXCEPT(except::InvalidParameterException,
                          "negative ellipticity parameter uncertainty");
    }
    _e1 = e1;
    _e1Sigma = e1Sigma;
    _e2 = e2;
    _e2Sigma = e2Sigma;
    _radius = radius;
    _radiusSigma = radiusSigma;
}

/** Computes the sample mean and standard error of the fluxes of each
  * source from this filter.
  */
void PerFilterSourceClusterAttributes::computeFlux(
    lsst::afw::detection::SourceSet const & sources)
{
    typedef detection::SourceSet::const_iterator Iter;
    size_t const n = sources.size();
    if (n == 0) {
        throw LSST_EXCEPT(except::InvalidParameterException,
                          "empty source set");
    } else if (n == 1) {
        detection::Source::Ptr src = sources.front();
        _flux = src->getPsfFlux();
        _fluxSigma = src->getPsfFluxErr();
        return;
    }
    // set filter flux to the sample mean
    double flux = 0.0;
    for (Iter i = sources.begin(), e = sources.end(); i != e; ++i) {
        flux += (*i)->getPsfFlux();
    }
    _flux = flux / n;
    // set filter flux uncertainty to an estimate of the standard
    // deviation of the sample mean 
    double ff = 0.0;
    for (Iter i = sources.begin(), e = sources.end(); i != e; ++i) {
        double df = (*i)->getPsfFlux() - _flux;
        ff += df * df;
    }
    _fluxSigma = std::sqrt(ff / (n * (n - 1)));
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
void PerFilterSourceClusterAttributes::computeEllipticities(
    lsst::afw::detection::SourceSet const & sources)
{
    typedef detection::SourceSet::const_iterator Iter;
    if (sources.empty()) {
        throw LSST_EXCEPT(except::InvalidParameterException,
                          "empty source set");
    }
    boost::scoped_array<Eigen::Vector3d> ellipticities(
        new Eigen::Vector3d[sources.size()]);
    Eigen::Vector3d m(0.0, 0.0, 0.0);
    int ns = 0;
    // TODO: should sources flagged as SHAPE_UNWEGHTED also be skipped?
    int skipMask = Flags::SHAPE_UNWEIGHTED_BAD;
    for (Iter i = sources.begin(), e = sources.end(); i != e; ++i) {
       if ((*i)->isNull(detection::IXX) ||
           (*i)->isNull(detection::IYY) ||
           (*i)->isNull(detection::IXY) ||
           ((*i)->getFlagForDetection() & skipMask) != 0) {
           continue;
       }
       double mxx = (*i)->getIxx();
       double myy = (*i)->getIyy();
       double mxy = (*i)->getIxy();
       // compute ellipticities from moments
       double t = mxx + myy;
       Eigen::Vector3d ep((mxx - myy) / t, 2.0 * mxy / t, std::sqrt(t));
       // store for later
       ellipticities[ns++] = ep;
       // and add to running sum
       m += ep;
    }
    if (ns > 0) {
        // set cluster ellipticities to sample means
        m /= static_cast<double>(ns);
        _e1 = m(0);
        _e2 = m(1);
        _radius = m(2);
        // If ns = 1, we could estimate ellipticity uncertainties from the
        // covariance matrix of the moments. Unfortunately, only a subset of
        // this matrix is stored in Source, so unless ns > 1, the ellipticity
        // uncertainties are null.
        if (ns > 1) {
            // compute standard errors
            Eigen::Vector3d v(0.0, 0.0, 0.0);
            for (int i = 0; i < ns; ++i) {
                v += (ellipticities[i] - m).cwise() * (ellipticities[i] - m);
            }
            _e1Sigma = std::sqrt(v[0] / (ns * (ns - 1)));
            _e2Sigma = std::sqrt(v[1] / (ns * (ns - 1)));
            _radiusSigma = std::sqrt(v[2] / (ns * (ns - 1)));
        }
    }
}


// -- SourceClusterAttributes ----

SourceClusterAttributes::SourceClusterAttributes() :
    _clusterId(0),
    _numObs(0),
    _earliestObsTime(0.0), _latestObsTime(0.0),
    _ra(0.0), _dec(0.0),
    _raSigma(0.0), _decSigma(0.0),
    _raDecCov(),
    _perFilterAttributes()
{ }

/** Creates a new SourceClusterAttributes with the given id, computing
  * attributes from the given set of sources.
  */
SourceClusterAttributes::SourceClusterAttributes(
    lsst::afw::detection::SourceSet const & sources,
    int64_t id) :
    _clusterId(id),
    _numObs(static_cast<int>(sources.size())),
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
            j = filters.insert(std::make_pair(filterId, detection::SourceSet())).first;
        }
        j->second.push_back(*i);
    }
    // compute per-filter properties
    for (HashIter i = filters.begin(), e = filters.end(); i != e; ++i) {
        _perFilterAttributes.insert(std::make_pair(
            i->first, PerFilterSourceClusterAttributes(i->second)));
    }
}

SourceClusterAttributes::~SourceClusterAttributes() { }

/** Returns @c true if the cluster has attributes specific to the given filter.
  */
bool SourceClusterAttributes::hasFilter(int filterId) const {
    return _perFilterAttributes.find(filterId) != _perFilterAttributes.end();
}

/** Returns a vector of the filter ids the cluster has attributes for. 
  */
std::vector<int> const SourceClusterAttributes::getFilterIds() const {
    typedef PerFilterAttributeMap::const_iterator Iter;
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
    typedef PerFilterAttributeMap::const_iterator Iter;
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
                                          double raSigma,
                                          double decSigma,
                                          Nullable<double> const & raDecCov)
{
    if (raSigma < 0.0 || decSigma < 0.0) {
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
    typedef PerFilterAttributeMap::iterator Iter;
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
        detection::Source::Ptr src = sources.front();
        _ra = src->getRa();
        _dec = src->getDec();
        _raSigma = src->getRaAstromErr();
        _decSigma = src->getDecAstromErr();
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
    _ra = (d2 == 0.0) ? 0.0 : atan2(y, x);
    _dec = (z == 0.0) ? 0.0 : atan2(z, std::sqrt(d2)); 
    // Compute covariance matrix and store mean position
    // as the object position of each source.
    double covRaRa = 0.0;
    double covDecDec = 0.0;
    double covRaDec = 0.0;
    for (Iter i = sources.begin(), e = sources.end(); i != e; ++i) {
       double dr = (*i)->getRa() - _ra;
       double dd = (*i)->getDec() - _dec;
       covRaRa += dr * dr;
       covDecDec += dd * dd;
       covRaDec += dr * dd;
    }
    // Set the ra/dec uncertainties to an estimate of the
    // of the standard deviation of the sample mean
    _raSigma = std::sqrt(covRaRa / (n * (n - 1)));
    _decSigma = std::sqrt(covDecDec / (n * (n - 1)));
    _raDecCov = covRaDec / n;
}

}}} // namespace lsst:ap::cluster

