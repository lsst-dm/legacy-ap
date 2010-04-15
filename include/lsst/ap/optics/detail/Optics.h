// -*- lsst-c++ -*-
/** @file
  * @brief Class that encapsulates the state required by the OPTICS algorithm.
  *
  * @ingroup ap
  * @author Serge Monkewitz
  */
#ifndef LSST_AP_OPTICS_DETAIL_OPTICS_H
#define LSST_AP_OPTICS_DETAIL_OPTICS_H

#include <vector>

#include "boost/scoped_array.hpp"
#include "boost/scoped_ptr.hpp"

#include "lsst/pex/logging/Log.h"

#include "KDTree.h"
#include "SeedList.h"


namespace lsst { namespace ap { namespace optics { namespace detail {

/** @internal
  * Class that encapsulates parameters and state operated on by the
  * OPTICS algorithm.
  */
template <int K, typename DataT>
class LSST_AP_LOCAL Optics {
public:
    Optics(Point<K, DataT> * points,
           int numPoints,
           int minPoints,
           double epsilon,
           double leafExtentThreshold,
           int pointsPerLeaf);
    ~Optics();

    template <typename MetricT>
    void run(std::vector<std::vector<DataT> > & clusters,
             MetricT const & metric);

private:
    Point<K, DataT> * _points;
    boost::scoped_ptr<KDTreePoint<K, DataT> > _tree;
    boost::scoped_ptr<SeedListPoint<K, DataT> > _seeds;
    boost::scoped_array<double> _distances;
    double _epsilon;
    int _numPoints;
    int _minPoints;
    bool _ran;
    lsst::pex::logging::Log _log;

    template <typename MetricT>
    void expandClusterOrder(int i, MetricT const & metric);
};

}}}} // namespace lsst:ap::optics::detail

#endif // LSST_AP_OPTICS_DETAIL_OPTICS_H