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
  * @brief Main OPTICS algorithm implementation class.
  *
  * For details of the algorithm, see the following paper:
  *
  * "OPTICS: Ordering Points To Identify the Clustering Structure".
  * Mihael Ankerst, Markus M. Breunig, Hans-Peter Kriegel, Jorg Sander (1999).
  * ACM SIGMOD international conference on Management of data.
  * ACM Press. pp. 49-60.
  *
  * @ingroup ap
  * @author Serge Monkewitz
  */
#ifndef LSST_AP_CLUSTER_DETAIL_OPTICS_H
#define LSST_AP_CLUSTER_DETAIL_OPTICS_H

#include <vector>

#include "boost/scoped_array.hpp"
#include "boost/scoped_ptr.hpp"
#include "boost/shared_ptr.hpp"

#include "lsst/pex/logging/Log.h"

#include "KDTree.h"
#include "SeedList.h"


namespace lsst { namespace ap { namespace cluster { namespace detail {

/** @internal
  * Class that encapsulates parameters and state operated on by the
  * OPTICS algorithm.
  */
template <int K, typename RecordT>
class Optics {
public:
    Optics(Point<K, boost::shared_ptr<RecordT> > * points,
           int numPoints,
           int minNeighbors,
           double epsilon,
           double leafExtentThreshold,
           int pointsPerLeaf);
    ~Optics();

    template <typename MetricT>
    void run(boost::shared_ptr<typename RecordT::Table> table,
             std::vector<typename RecordT::Catalog> & clusters,
             MetricT const & metric);

private:
    typedef boost::shared_ptr<RecordT> DataT;

    Point<K, DataT> * _points;
    boost::scoped_ptr<KDTree<K, DataT> > _tree;
    boost::scoped_ptr<SeedList<K, DataT> > _seeds;
    boost::scoped_array<double> _distances;
    double _epsilon;
    int _numPoints;
    int _minNeighbors;
    bool _ran;
    lsst::pex::logging::Log _log;

    template <typename MetricT>
    void expandClusterOrder(int i, MetricT const & metric);
};

}}}} // namespace lsst:ap::cluster::detail

#endif // LSST_AP_CLUSTER_DETAIL_OPTICS_H
