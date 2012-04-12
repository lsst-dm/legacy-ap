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
  * @brief Implementation of the Optics class.
  *
  * @ingroup ap
  * @author Serge Monkewitz
  */
#ifndef LSST_AP_CLUSTER_OPTICS_OPTICS_CC
#define LSST_AP_CLUSTER_OPTICS_OPTICS_CC

#include "Optics.h"

#include <algorithm>
#include <limits>

#include "lsst/pex/exceptions.h"

#include "KDTree.cc"
#include "SeedList.cc"


namespace lsst { namespace ap { namespace cluster { namespace optics {


/** Initializes data structures required by the OPTICS to run over the given
  * set of points.
  */
template <int K, typename RecordT>
Optics<K, RecordT>::Optics(Point<K, boost::shared_ptr<RecordT> > * points,
                         int numPoints,
                         int minNeighbors,
                         double epsilon,
                         double leafExtentThreshold,
                         int pointsPerLeaf) :
    _points(points),
    _tree(),
    _seeds(),
    _distances(),
    _epsilon(epsilon),
    _numPoints(numPoints),
    _minNeighbors(minNeighbors),
    _ran(false),
    _log(lsst::pex::logging::Log::getDefaultLog(), "lsst.ap.cluster.optics")
{
    if (_points == 0) {
        throw LSST_EXCEPT(lsst::pex::exceptions::InvalidParameterException,
                          "Input point array is null");
    }
    if (_numPoints <= 0) {
        throw LSST_EXCEPT(lsst::pex::exceptions::InvalidParameterException,
                          "Number of input points must be at least 1");
    }
    if (_minNeighbors < 0) {
        throw LSST_EXCEPT(lsst::pex::exceptions::InvalidParameterException,
                          "OPTICS minNeighbors parameter value is negative");
    }
    if (_epsilon < 0.0) {
        throw LSST_EXCEPT(lsst::pex::exceptions::InvalidParameterException,
                          "OPTICS epsilon parameter value is negative");
    }
    if (pointsPerLeaf <= 0) {
        throw LSST_EXCEPT(lsst::pex::exceptions::InvalidParameterException, 
                          "K-D tree pointsPerLeaf parameter must be positive");
    }

    _log.log(lsst::pex::logging::Log::INFO, "Building k-d tree for sources");
    boost::scoped_ptr<KDTree<K, DataT> > tree(new KDTree<K, DataT>(
        points, numPoints, pointsPerLeaf, leafExtentThreshold));
    _log.format(lsst::pex::logging::Log::INFO,
                "Created k-d tree for %d sources", numPoints);

    boost::scoped_ptr<SeedList<K, DataT> > seeds(new SeedList<K, DataT>(
        points, numPoints));
    boost::scoped_array<double> distances(new double[_minNeighbors]);
    using std::swap;
    swap(_tree, tree);
    swap(_seeds, seeds);
    swap(_distances, distances);
}

template <int K, typename RecordT>
Optics<K, RecordT>::~Optics() { }

/** Runs the OPTICS algorithm, appending clusters to @c clusters.
  * This method may only be called once for a given Optics instance.
  */
template <int K, typename RecordT>
    template <typename MetricT>
void Optics<K, RecordT>::run(boost::shared_ptr<typename RecordT::Table> table,
                             std::vector<typename RecordT::Catalog> & clusters,
                             MetricT const & metric)
{
    if (_ran) {
        throw LSST_EXCEPT(lsst::pex::exceptions::LogicErrorException,
                          "OPTICS has already been run");
    }
    typename RecordT::Catalog cluster(table);
    size_t const s = clusters.size();
    int scanFrom = 0;

    _log.log(lsst::pex::logging::Log::INFO, "Clustering sources using OPTICS");
    _ran = true;

    while (true) {
        int i;
        if (_seeds->empty()) {
            // find next unprocessed point
            for (i = scanFrom; i < _numPoints; ++i) {
                if (_points[i].state == Point<K, DataT>::UNPROCESSED) {
                    scanFrom = i + 1;
                    break;
                }
            }
            if (i == _numPoints) {
                break;
            }
            _points[i].state = Point<K, DataT>::PROCESSED;
            expandClusterOrder(i, metric);
            if (cluster.size() > 0) {
                // clusters of size 1 are generated for noise sources
                clusters.push_back(cluster);
                cluster.clear();
            }
            cluster.push_back(_points[i].data);
        } else {
            // expand cluster around seed with smallest reachability-distance
            i = _seeds->pop();
            expandClusterOrder(i, metric);
            assert(_points[i].reach != std::numeric_limits<double>::infinity());
            cluster.push_back(_points[i].data);
        }
    }
    clusters.push_back(cluster);
    _log.format(lsst::pex::logging::Log::INFO, "Produced %d clusters",
                static_cast<int>(clusters.size() - s));
}

template <int K, typename RecordT>
    template <typename MetricT>
void Optics<K, RecordT>::expandClusterOrder(int i, MetricT const & metric)
{
    // find epsilon neighborhood of point i
    int range = _tree->inRange(_points[i].coords, _epsilon, metric);
    // compute core-distance
    int n = 0;
    int j = range;
    while (j != -1) {
        Point<K, DataT> * p = _points + j;
        if (j != i) {
            double d = p->dist;
            if (n < _minNeighbors) {
                _distances[n++] = d;
                std::push_heap(_distances.get(), _distances.get() + n);
            } else if (_distances[0] > d) {
                std::pop_heap(_distances.get(), _distances.get() + n);
                _distances[n - 1] = d;
                std::push_heap(_distances.get(), _distances.get() + n);
            }
        }
        j = p->next;
    }
    if (n == _minNeighbors) {
        // point i is a core-object. Update reachability-distance of all
        // points in the epsilon-neighborhood of point i.
        double coreDist = _distances[0];
        j = range;
        while (j != -1) {
            Point<K, DataT> * p = _points + j;
            if (p->state != Point<K, DataT>::PROCESSED) {
                _seeds->update(j, std::max(coreDist, p->dist));
            }
            j = p->next;
        }
    }
}

}}}} // namespace lsst:ap::cluster::optics

#endif // LSST_AP_CLUSTER_OPTICS_OPTICS_CC
