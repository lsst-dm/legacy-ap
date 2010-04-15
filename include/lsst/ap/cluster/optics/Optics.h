// -*- lsst-c++ -*-
/** @file
  * @brief Main OPTICS algorithm implementation class.
  *
  * For details of the algorithm, see the following paper:
  *
  * "OPTICS: Ordering Points To Identify the Clustering Structure".
  * Mihael Ankerst, Markus M. Breunig, Hans-Peter Kriegel, JÃÂ¶rg Sander (1999).
  * ACM SIGMOD international conference on Management of data.
  * ACM Press. pp. 49Ã¢<93>60.
  *
  * @ingroup ap
  * @author Serge Monkewitz
  */
#ifndef LSST_AP_CLUSTER_OPTICS_OPTICS_H
#define LSST_AP_CLUSTER_OPTICS_OPTICS_H

#include <vector>

#include "boost/scoped_array.hpp"
#include "boost/scoped_ptr.hpp"

#include "lsst/pex/logging/Log.h"

#include "KDTree.h"
#include "SeedList.h"


namespace lsst { namespace ap { namespace cluster { namespace optics {

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
    boost::scoped_ptr<KDTree<K, DataT> > _tree;
    boost::scoped_ptr<SeedList<K, DataT> > _seeds;
    boost::scoped_array<double> _distances;
    double _epsilon;
    int _numPoints;
    int _minPoints;
    bool _ran;
    lsst::pex::logging::Log _log;

    template <typename MetricT>
    void expandClusterOrder(int i, MetricT const & metric);
};

}}}} // namespace lsst:ap::cluster::optics

#endif // LSST_AP_CLUSTER_OPTICS_OPTICS_H
