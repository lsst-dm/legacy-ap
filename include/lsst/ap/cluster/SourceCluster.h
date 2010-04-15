// -*- lsst-c++ -*-
/** @file
  * @brief Functions for clustering points (using the OPTICS algorithm)
  *        and for computing cluster attributes.
  *
  * @ingroup ap
  * @author Serge Monkewitz
  */
#ifndef LSST_AP_CLUSTER_SOURCECLUSTER_H
#define LSST_AP_CLUSTER_SOURCECLUSTER_H

#include <vector>

#include "boost/shared_ptr.hpp"

#include "lsst/pex/policy.h"
#include "lsst/afw/detection/Source.h"

#include "../Common.h"


namespace lsst { namespace ap { namespace cluster {

/** Attributes derived from a cluster of Sources; i.e. Sources that
  * have been determined to be observations of the same physical object
  * according to some clustering algorithm.
  */
class LSST_AP_API SourceClusterAttributes {
public:
    typedef boost::shared_ptr<SourceClusterAttributes> Ptr;

    SourceClusterAttributes();
    ~SourceClusterAttributes();

private:
    int64_t _clusterId;

    // average position, proper motion estimate and errors
    double _raAvg;
    double _declAvg;
    double _raAvgSigma;
    double _declAvgSigma;
    double _raDeclAvgCov;
    double _muRa;
    double _muDecl;
    double _muRaSigma;
    double _muDeclSigma;
    double _muRaDeclCov;

    // per filter properties
    // TODO
};

LSST_AP_API std::vector<lsst::afw::detection::SourceSet> cluster(
    lsst::afw::detection::SourceSet const & sources,
    lsst::pex::policy::Policy::Ptr policy);

LSST_AP_API std::vector<lsst::afw::detection::SourceSet> cluster(
    std::vector<lsst::afw::detection::SourceSet> const & sources,
    lsst::pex::policy::Policy::Ptr policy);

}}} // namespace lsst:ap::cluster

#endif // LSST_AP_CLUSTER_SOURCECLUSTER_H
