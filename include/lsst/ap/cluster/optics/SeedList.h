// -*- lsst-c++ -*-
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
class LSST_AP_LOCAL SeedList {
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
