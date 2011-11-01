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
  * @brief Class that maintains the OPTICS seed list.
  *
  * @ingroup ap
  * @author Serge Monkewitz
  */
#ifndef LSST_AP_CLUSTER_OPTICS_SEEDLIST_H
#define LSST_AP_CLUSTER_OPTICS_SEEDLIST_H

#include "../../Common.h"
#include "KDTree.h"


namespace lsst { namespace ap { namespace cluster { namespace optics {

/** Class for managing the OPTICS ordered seeds. Methods to add a seed,
  * remove the seed with smallest reachability distance and to decrease the 
  * reachability of a seed are provided.
  */
template <int K, typename DataT>
class SeedList {
public:
    SeedList(Point<K, DataT> * points, int numPoints);
    ~SeedList();

    inline bool empty() const;
    inline int size() const;
    inline int capacity() const;
    inline int pop();
    inline void add(int i);
    inline void update(int i, double reach);

    bool checkInvariants() const;

private:
    boost::scoped_array<int> _heap;
    Point<K, DataT> * _points;
    int _size;
    int _numPoints;

    inline void siftUp(int heapIndex, int pointIndex);
    inline void siftDown(int pointIndex);
};

}}}} // namespace lsst:ap::cluster::optics

#endif // LSST_AP_CLUSTER_OPTICS_SEEDLIST_H
