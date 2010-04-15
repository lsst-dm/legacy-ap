// -*- lsst-c++ -*-
/** @file
  * @brief Implementation of the Optics class.
  *
  * @ingroup ap
  * @author Serge Monkewitz
  */
#ifndef LSST_AP_OPTICS_DETAIL_OPTICS_CC
#define LSST_AP_OPTICS_DETAIL_OPTICS_CC

#include "Optics.h"

#include <algorithm>
#include <limits>

#include "lsst/pex/exceptions.h"


namespace lsst { namespace ap { namespace optics { namespace detail {

namespace except = ::lsst::pex::exceptions;

/** Initializes data structures required by the OPTICS to run over the given
  * set of points.
  */
template <int K, typename DataT>
Optics<K, DataT>::Optics(Point<K, DataT> * points,
                         int numPoints,
                         int minPoints,
                         double epsilon,
                         double leafExtentThreshold,
                         int pointsPerLeaf) :
    _points(points),
    _tree(),
    _seeds(),
    _distances(),
    _epsilon(epsilon),
    _numPoints(numPoints),
    _minPoints(minPoints),
    _ran(false),
    _log(lsst::pex::logging::Log::getDefaultLog(), "lsst.ap.optics")
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
                          "minPoints parameter value is negative");
    }
    if (_epsilon < 0.0) {
        throw LSST_EXCEPT(except::InvalidParameterException, "OPTICS "
                          "epsilon (clustering distance) parameter value "
                          "is negative");
    }
    if (pointsPerLeaf <= 0) {
        throw LSST_EXCEPT(except::InvalidParameterException, "K-D tree "
                          "pointsPerLeaf parameter must be positive");
    }

    _log.log(Log::INFO, "Building k-d tree for sources");
    boost::scoped_ptr<KDTree<K, DataT> > tree(new KDTree<K, DataT>(
        points, numPoints, pointsPerLeaf, leafExtentThreshold));
    _log.format(Log::INFO, "Created k-d tree for %d sources", numPoints);

    boost::scoped_ptr<SeedList<K, DataT> > seeds(new SeedList<K, DataT>(
        points, numPoints));
    boost::scoped_array<double> distances(new double[_minPoints]);
    using std::swap;
    swap(_tree, tree);
    swap(_seeds, seeds);
    swap(_distances, distances);
}

template <int K, typename DataT>
Optics<K, DataT>::~Optics() { }

/** Runs the OPTICS algorithm, appending clusters to @c clusters.
  * This method may only be called once for a given Optics instance.
  */
template <int K, typename DataT>
    template <typename MetricT>
void Optics<K, DataT>::run(std::vector<std::vector<DataT> > & clusters,
                           MetricT const & metric)
{
    if (_ran) {
        throw LSST_EXCEPT(except::LogicErrorException,
                          "OPTICS has already been run");
    }
    std::vector<DataT> cluster;
    size_t s = clusters.size();
    int scanFrom = 0;

    _log.log(Log::INFO, "Clustering sources using OPTICS");
    _ran = true;

    while (true) {
        int i;
        if (_seeds->empty()) {
            // find next unprocessed point
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
        expandClusterOrder(i, metric);
        if (_points[i].reach == std::numeric_limits<double>::infinity()) {
            if (clusters.size() != 0) {
                clusters.push_back(cluster);
            }
            cluster.clear();
        }
        cluster.push_back(*(_points[i].data));
    }
    if (cluster.size() > 0) {
        clusters.push_back(cluster);
    }
    _log.format(Log::INFO, "Produced %d clusters",
                static_cast<int>(clusters.size() - s));
}

template <int K, typename DataT>
    template <typename MetricT>
void Optics<K, DataT>::expandClusterOrder(int i, MetricT const & metric)
{
    // find epsilon neighborhood of point i
    int range = _tree->inRange(_points[i].coords, _epsilon, metric);
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

}}}} // namespace lsst:ap::optics::detail

#endif // LSST_AP_OPTICS_DETAIL_OPTICS_CC