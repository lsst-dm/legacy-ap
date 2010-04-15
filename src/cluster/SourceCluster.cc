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

#include "boost/scoped_array.hpp"

#include "lsst/pex/exceptions.h"
#include "lsst/pex/policy.h"

#include "lsst/ap/Common.h"
#include "lsst/ap/cluster/optics/Metrics.h"
#include "lsst/ap/cluster/optics/Optics.cc"


namespace except = lsst::pex::exceptions;
namespace detection = lsst::afw::detection;
namespace policy = lsst::pex::policy;


namespace lsst { namespace ap { namespace cluster {

SourceClusterAttributes::SourceClusterAttributes() { }

SourceClusterAttributes::~SourceClusterAttributes() { }


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
    double ra = ptr->getRa() * RADIANS_PER_DEGREE;
    double dec = ptr->getDec() * RADIANS_PER_DEGREE;
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


/// @name OPTICS clustering functions.
//@{
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
LSST_AP_API std::vector<lsst::afw::detection::SourceSet> cluster(
    lsst::afw::detection::SourceSet const & sources,
    lsst::pex::policy::Policy::Ptr policy)
{
    typedef detection::SourceSet::const_iterator Iter;
    if (sources.size() > MAX_SOURCES) {
        throw LSST_EXCEPT(except::InvalidParameterException,
                          "too many sources to cluster");
    }
    boost::scoped_array<Point> entries(new Point[sources.size()]);
    std::vector<lsst::afw::detection::SourceSet> clusters;
    // Transform sources into a form the OPTICS implementation understands
    int i = 0;
    for (Iter s = sources.begin(), e = sources.end(); s != e; ++s, ++i) {
        initPoint(entries[i], *s);
    }
    if (i > 0) {
        double epsilon = policy->getDouble("epsilon");
        double leafExtentThreshold = policy->getDouble("leafExtentThreshold");
        // Convert epsilon and leafExtentThreshold to radians, and account
        // for the fact that our metric is the squared euclidian distance,
        // not angular separation.
        epsilon = std::sin(0.5 * RADIANS_PER_ARCSEC * epsilon);
        leafExtentThreshold = std::sin(
            0.5 * RADIANS_PER_ARCSEC * leafExtentThreshold);
        epsilon = 4.0 * epsilon * epsilon;
        leafExtentThreshold = 4.0 * leafExtentThreshold * leafExtentThreshold;
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

LSST_AP_API std::vector<lsst::afw::detection::SourceSet> cluster(
    std::vector<lsst::afw::detection::SourceSet> & sources,
    lsst::pex::policy::Policy::Ptr policy)
{
    typedef std::vector<detection::SourceSet>::const_iterator SetIter;
    typedef detection::SourceSet::const_iterator Iter;
    // find total number of sources
    size_t n = 0;
    for (SetIter ss = sources.begin(), es = sources.end(); ss != es; ++ss) {
        n += ss->size();
    }
    if (n > MAX_SOURCES) {
        throw LSST_EXCEPT(except::InvalidParameterException,
                          "too many sources to cluster");
    }
    boost::scoped_array<Point> entries(new Point[n]);
    std::vector<lsst::afw::detection::SourceSet> clusters;
    // Transform sources into a form the OPTICS implementation understands
    int i = 0;
    for (SetIter ss = sources.begin(), es = sources.end(); ss != es; ++ss) {
        for (Iter s = ss->begin(), e = ss->end(); s != e; ++s, ++i) {
            initPoint(entries[i], *s);
        }
    }
    if (i > 0) {
        double epsilon = policy->getDouble("epsilon");
        double leafExtentThreshold = policy->getDouble("leafExtentThreshold");
        // Convert epsilon and leafExtentThreshold to radians, and account
        // for the fact that our metric is the squared euclidian distance,
        // not angular separation.
        epsilon = std::sin(0.5 * RADIANS_PER_ARCSEC * epsilon);
        leafExtentThreshold = std::sin(
            0.5 * RADIANS_PER_ARCSEC * leafExtentThreshold);
        epsilon = 4.0 * epsilon * epsilon;
        leafExtentThreshold = 4.0 * leafExtentThreshold * leafExtentThreshold;
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
//@}

}}} // namespace lsst:ap::cluster

