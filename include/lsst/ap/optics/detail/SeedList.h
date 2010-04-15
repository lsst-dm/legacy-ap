// -*- lsst-c++ -*-
/** @file
  * @brief Class that maintains the OPTICS seed list.
  *
  * @ingroup ap
  * @author Serge Monkewitz
  */
#ifndef LSST_AP_OPTICS_DETAIL_SEEDLIST_H
#define LSST_AP_OPTICS_DETAIL_SEEDLIST_H

#include "../../Common.h"
#include "KDTree.h"


namespace lsst { namespace ap { namespace optics { namespace detail {

/** TODO
  *
  */
template <int K, typename DataT>
class LSST_AP_LOCAL SeedList {
public:
    SeedList(Point<K, DataT> * points, int numPoints);
    ~SeedList();

    bool empty() const;
    int size() const;
    inline int capacity() const;
    inline int pop();
    inline void add(int i);
    inline void update(int i, double reach);

private:
    boost::scoped_array<int> _heap;
    Point<K, DataT> * _points;
    int _size;
    int _numPoints;

    inline void siftUp(int heapIndex, int pointIndex);
    inline void siftDown(int pointIndex);
};

}}}} // namespace lsst:ap::optics::detail

#endif // LSST_AP_OPTICS_DETAIL_SEEDLIST_H
