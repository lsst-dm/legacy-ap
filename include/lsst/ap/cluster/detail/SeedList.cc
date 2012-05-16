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
  * @brief Implementation of the SeedList class template.
  *
  * @ingroup ap
  * @author Serge Monkewitz
  */
#ifndef LSST_AP_CLUSTER_DETAIL_SEEDLIST_CC
#define LSST_AP_CLUSTER_DETAIL_SEEDLIST_CC

#include "SeedList.h"


namespace lsst { namespace ap { namespace cluster { namespace detail {

template <int K, typename DataT>
SeedList<K, DataT>::SeedList(Point<K, DataT> * points, int numPoints) :
    _heap(new int[numPoints]),
    _points(points),
    _size(0),
    _numPoints(numPoints)
{ }

template <int K, typename DataT>
SeedList<K, DataT>::~SeedList() { }

/// Returns @c true if this seed list contains no entries.
template <int K, typename DataT>
inline bool SeedList<K, DataT>::empty() const {
    return _size == 0;
}

/// Returns the number of entries in this seed list.
template <int K, typename DataT>
inline int SeedList<K, DataT>::size() const {
    return _size;
}

/// Returns the capacity of this seed list.
template <int K, typename DataT>
inline int SeedList<K, DataT>::capacity() const {
    return _numPoints;
}

/** Removes and returns the point with smallest reachability-distance.
  *
  * @return The index of the point with smallest reachability-distance,
  *         or -1 if the seed list is empty.
  */
template <int K, typename DataT>
inline int SeedList<K, DataT>::pop() {
    int s = _size;
    if (s == 0) {
        return -1;
    }
    int smallest = _heap[0];
    _points[smallest].state = Point<K, DataT>::PROCESSED;
    _size = --s;
    if (s > 1) {
        siftDown(_heap[s]);
    } else if (s == 1) {
        int i = _heap[1];
        _heap[0] = i;
        _points[i].state = 0;
    }
    return smallest;
}

/** Adds the i-th point to this SeedList.
  *
  * @pre i >= 0 && i < capacity()
  * @pre size() < capacity()
  */
template <int K, typename DataT>
inline void SeedList<K, DataT>::add(int i) {
    assert(i >= 0 && i < _numPoints);
    assert(_size < _numPoints);
    int s = _size;
    _size = s + 1;
    if (s == 0) {
        _heap[0] = i;
        _points[i].state = 0;
    } else {
        siftUp(s, i);
    }
}

/** Updates the reachability-distance of the i-th point. If it is not
  * already in the seed list, then it is added. Otherwise, if the
  * new reachability-distance is smaller than the current one, the
  * i-th points reachability-distance is updated.
  *
  * @pre i >= 0 && i < capacity()
  */
template <int K, typename DataT>
inline void SeedList<K, DataT>::update(int i, double reach) {
    assert(i >= 0 && i < _numPoints);
    int heapIndex = _points[i].state;
    if (heapIndex >= 0) {
        assert(_heap[heapIndex] == i);
        // the i-th point is already in the seed list
        if (reach < _points[i].reach) {
            _points[i].reach = reach;
            siftUp(heapIndex, i);
        }
    } else {
        // add i-th point to the seed list
        _points[i].reach = reach;
        add(i);
    }
}

template <int K, typename DataT>
inline void SeedList<K, DataT>::siftUp(int heapIndex, int pointIndex) {
    assert(heapIndex < _size);
    assert(pointIndex >= 0 && pointIndex < _numPoints);
    double reach = _points[pointIndex].reach;
    while (heapIndex > 0) {
        int parentHeapIndex = (heapIndex - 1) >> 1;
        int parentPointIndex = _heap[parentHeapIndex];
        if (_points[parentPointIndex].reach <= reach) {
            break;
        }
        _heap[heapIndex] = parentPointIndex;
        _points[parentPointIndex].state = heapIndex;
        heapIndex = parentHeapIndex;
    }
    _heap[heapIndex] = pointIndex;
    _points[pointIndex].state = heapIndex;
}

template <int K, typename DataT>
inline void SeedList<K, DataT>::siftDown(int pointIndex) {
    assert(pointIndex >= 0 && pointIndex < _numPoints);
    double reach = _points[pointIndex].reach;
    int halfSize = _size >> 1;
    int heapIndex = 0;
    while (heapIndex < halfSize) {
        int childHeapIndex = (heapIndex << 1) + 1;
        int siblingHeapIndex = childHeapIndex + 1;
        int childPointIndex = _heap[childHeapIndex];
        double childReach = _points[childPointIndex].reach;
        if (siblingHeapIndex < _size) {
            int siblingPointIndex = _heap[siblingHeapIndex];
            double siblingReach = _points[siblingPointIndex].reach;
            if (siblingReach < childReach) {
                childReach = siblingReach;
                childPointIndex = siblingPointIndex;
                childHeapIndex = siblingHeapIndex;
            }
        }
        if (reach <= childReach) {
            break;
        }
        _heap[heapIndex] = childPointIndex;
        _points[childPointIndex].state = heapIndex;
        heapIndex = childHeapIndex;
    }
    _heap[heapIndex] = pointIndex;
    _points[pointIndex].state = heapIndex;
}

/** Returns @c true if implementation defined invariants over internal state
  * hold. To be used by unit tests.
  */
template <int K, typename DataT>
bool SeedList<K, DataT>::checkInvariants() const {
    // check that each point knows its location in the seed list
    for (int i = 0; i < _numPoints; ++i) {
        int h = _points[i].state;
        if (h >= 0) {
            if (h >= _size) {
                // point has an invalid index into the seed list
                return false;
            }
            if (_heap[h] != i) {
                // point has an incorrect index into the seed list
                return false;
            }
        }
    }
    for (int i = 0; i < _size; ++i) {
        int p = _heap[i];
        if (p < 0 || p >= _numPoints) {
            // heap contains an invalid point index
            return false;
        }
        if (_points[p].state != i) {
            // point has an incorrect index into the seed list
            return false;
        }
    }
    // check the heap invariant
    for (int i = 0; i < _size >> 1; ++i) {
        double reach = _points[_heap[i]].reach;
        int h = (i << 1) + 1;
        int p = _heap[h];
        if (_points[p].reach < reach) {
            return false;
        }
        h += 1;
        if (h < _size) {
            int p = _heap[h];
            if (_points[p].reach < reach) {
                return false;
            }
        }
    }
    return true;
}

}}}} // namespace lsst:ap::cluster::detail

#endif // LSST_AP_CLUSTER_DETAIL_SEEDLIST_CC
