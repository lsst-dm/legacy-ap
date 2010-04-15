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
#ifndef LSST_AP_OPTICS_OPTICS_H
#define LSST_AP_OPTICS_OPTICS_H

#include "lsst/ap/optics/Optics.h"

#include <cmath>
#include <algorithm>
#include <limits>

#include "boost/scoped_ptr.hpp"

#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Log.h"
#include "lsst/pex/policy.h"

#include "lsst/ap/Common.h"
#include "lsst/ap/optics/SourceCluster.h"
#include "lsst/ap/optics/detail/KDTree.cc"
#include "lsst/ap/optics/detail/Metrics.h"
#include "lsst/ap/optics/detail/SeedList.cc"


namespace except = lsst::pex::exceptions;
namespace detection = lsst::afw::detection;
namespace policy = lsst::pex::policy;

using lsst::pex::logging::Log;


// explicit instantiations for implementation templates


namespace lsst { namespace ap { namespace optics {

namespace {

typedef detail::Point<3, detection::Source::Ptr> Point;
typedef detail::KDTree<3, detection::Source::Ptr> KDTree;
typedef detail::SeedList<3, detection::Source::Ptr> SeedList;

/** @internal
  * Class that encapsulates parameters and state operated on by the
  * OPTICS algorithm.
  */
class LSST_AP_LOCAL Optics {
public:
    Optics(Point * points,
           int numPoints,
           lsst::pex::policy::Policy::ConstPtr policy);
    ~Optics();

    void run(SourceClusterVector & clusters);

private:
    Point * _points;
    boost::scoped_ptr<KDTree> _tree;
    boost::scoped_ptr<SeedList> _seeds;
    boost::scoped_array<double> _distances;
    double _epsilon;
    int _numPoints;
    int _minPoints;
    bool _ran;
    Log _log;

    void expandClusterOrder(int i);
};

/** Initializes data structures required by OPTICS for a set of input
  * points. The following parameters are read from @c policy:
  *
  * @li @c "epsilon" (double) : generating distance for clusters (arcsec).
  * @li @c "minPoints" (int) : minimum number of points that must be in an
  *     epsilon neighborhood of a point P for P to be assigned to a cluster.
  * @li @c "pointsPerLeaf" (int) : a performance tuning parameter that
  *     specifies the target number of points per leaf of the k-d tree used
  *     by OPTICS internally.
  */
Optics::Optics(Point * points,
               int numPoints,
               lsst::pex::policy::Policy::ConstPtr policy) :
    _points(points),
    _tree(),
    _seeds(),
    _distances(),
    _epsilon(RADIANS_PER_DEGREE * (policy->getDouble("epsilon") / 3600.0)),
    _numPoints(numPoints),
    _minPoints(policy->getInt("minPoints")),
    _ran(false),
    _log(Log::getDefaultLog(), "lsst.ap.optics")
{
    if (_points == 0) {
        throw LSST_EXCEPT(except::InvalidParameterException,
                          "Input point array is null");
    }
    if (_numPoints <= 0) {
        throw LSST_EXCEPT(except::InvalidParameterException,
                          "Number of input points must be at least 1");
    }
    if (_minPoints < 0) {
        throw LSST_EXCEPT(except::InvalidParameterException, "OPTICS "
                          "\"minPoints\" policy parameter value is negative");
    }
    if (_epsilon < 0.0) {
        throw LSST_EXCEPT(except::InvalidParameterException, "OPTICS "
                          "\"epsilon\" (clustering distance) policy "
                          "parameter value is negative");
    }
    int pointsPerLeaf = policy->getInt("pointsPerLeaf");
    if (pointsPerLeaf <= 0) {
        throw LSST_EXCEPT(except::InvalidParameterException, "OPTICS "
                          "\"pointsPerLeaf\" policy parameter must be positive");
    }
    // Set leaf extent threshold to twice the clustering distance. Also
    // account for the fact that our metric is the squared euclidian distance.
    double leafExtentThreshold = std::sin(_epsilon);
    leafExtentThreshold = 4.0 * leafExtentThreshold * leafExtentThreshold;
    double e = std::sin(0.5 * _epsilon);
    _epsilon = 4.0 * e * e;

    _log.log(Log::INFO, "Building k-d tree for sources");
    boost::scoped_ptr<KDTree> tree(
        new KDTree(points, numPoints, pointsPerLeaf, leafExtentThreshold));
    _log.format(Log::INFO, "Created k-d tree for %d sources", numPoints);

    boost::scoped_ptr<SeedList> seeds(new SeedList(points, numPoints));
    boost::scoped_array<double> distances(new double[_minPoints]);
    using std::swap;
    swap(_tree, tree);
    swap(_seeds, seeds);
    swap(_distances, distances);
}

Optics::~Optics() { }

/** Runs the OPTICS algorithm, appending clusters to @c clusters.
  * This method may only be called once for a given Optics instance.
  */
void Optics::run(SourceClusterVector & clusters)
{
    if (_ran) {
        throw LSST_EXCEPT(except::LogicErrorException,
                          "OPTICS has already been run");
    }
    SourceCluster cluster;
    size_t s = clusters.size();
    int scanFrom = 0;

    _log.log(Log::INFO, "Clustering sources using OPTICS");
    _ran = true;

    while (true) {
        int i;
        if (_seeds->empty()) {
            // add next unprocessed point to seed list
            for (i = scanFrom; i < _numPoints; ++i) {
                if (_points[i].state == Point::UNPROCESSED) {
                    scanFrom = i + 1;
                    break;
                }
            }
            if (i == _numPoints) {
                break;
            }
        } else {
            i = _seeds->pop();
        }
        _points[i].state = Point::PROCESSED;
        expandClusterOrder(i);
        if (_points[i].reach == std::numeric_limits<double>::infinity()) {
            if (clusters.size() != 0) {
                clusters.push_back(cluster);
            }
            cluster.sources.clear();
        }
        cluster.sources.push_back(*(_points[i].data));
    }
    if (cluster.sources.size() > 0) {
        clusters.push_back(cluster);
    }
    _log.format(Log::INFO, "Produced %d clusters",
                static_cast<int>(clusters.size() - s));
}

void Optics::expandClusterOrder(int i)
{
    // find epsilon neighborhood of point i
    int range = _tree->inRange(_points[i].coords, _epsilon,
                               detail::SquaredEuclidianDistanceOverSphere());
    // compute core-distance
    int n = 0;
    int j = range;
    while (j != -1) {
        Point * p = _points + j;
        if (j != i) {
            double d = p->dist;
            if (n < _minPoints) {
                _distances[n++] = d;
                std::push_heap(_distances.get(), _distances.get() + n,
                               std::greater<double>());
            } else if (_distances[0] > d) {
                std::pop_heap(_distances.get(), _distances.get() + n,
                              std::greater<double>());
                _distances[n - 1] = d;
                std::push_heap(_distances.get(), _distances.get() + n,
                               std::greater<double>());
            }
        }
        j = p->next;
    }
    if (n == _minPoints) {
        // point i is a core-object. Update reachability-distance of all
        // points in the epsilon-neighborhood of point i.
        double coreDist = _distances[0];
        j = range;
        while (j != -1) {
            Point * p = _points + j;
            if (j != i) {
                _seeds->update(i, std::max(coreDist, p->dist));
            }
            j = p->next;
        }
    }
}

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
/** @brief Clusters a set of sources using the OPTICS algorithm.
  *
  * @param[in] sources      The sources to cluster.
  * @param[in] epsilon      The clustering distance to use (radians).
  * @param[in] minPoints    The minimum number of points that must be in an
  *                         epsilon-neighborhood of a source S for S to be
  *                         assigned to a cluster.
  *
  * @return     A vector of SourceCluster objects.
  */
LSST_AP_API SourceClusterVector cluster(
    lsst::afw::detection::SourceSet const & sources,
    lsst::pex::policy::Policy::Ptr policy)
{
    typedef detection::SourceSet::const_iterator Iter;
    if (sources.size() > MAX_SOURCES) {
        throw LSST_EXCEPT(except::InvalidParameterException,
                          "too many sources to cluster");
    }
    boost::scoped_array<Point> entries(new Point[sources.size()]);
    SourceClusterVector clusters;
    // Transform sources into a form the OPTICS implementation understands
    int i = 0;
    for (Iter s = sources.begin(), e = sources.end(); s != e; ++s, ++i) {
        initPoint(entries[i], *s);
    }
    if (i > 0) {
        Optics optics(entries.get(), i, policy);
        optics.run(clusters);
    }
    return clusters;
}

LSST_AP_API SourceClusterVector cluster(
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
    SourceClusterVector clusters;
    // Transform sources into a form the OPTICS implementation understands
    int i = 0;
    for (SetIter ss = sources.begin(), es = sources.end(); ss != es; ++ss) {
        for (Iter s = ss->begin(), e = ss->end(); s != e; ++s, ++i) {
            initPoint(entries[i], *s);
        }
    }
    if (i > 0) {
        Optics optics(entries.get(), i, policy);
        optics.run(clusters);
    }
    return clusters;
}
//@}

}}} // namespace lsst:ap::optics

#endif // LSST_AP_OPTICS_OPTICS_H
