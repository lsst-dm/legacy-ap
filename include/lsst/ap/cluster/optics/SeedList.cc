// -*- lsst-c++ -*-
/** @file
  * @brief Implementation of the SeedList class template.
  *
  * @ingroup ap
  * @author Serge Monkewitz
  */
#ifndef LSST_AP_CLUSTER_OPTICS_SEEDLIST_CC
#define LSST_AP_CLUSTER_OPTICS_SEEDLIST_CC

#include "SeedList.h"


namespace lsst { namespace ap { namespace cluster { namespace optics {

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
    int i = _heap[--s];
    _size = s;
    if (s > 0) {
        siftDown(i);
    } else {
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
        if (_points[parentPointIndex].reach < reach) {
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

}}}} // namespace lsst:ap::cluster::optics

#endif // LSST_AP_CLUSTER_OPTICS_SEEDLIST_CC
