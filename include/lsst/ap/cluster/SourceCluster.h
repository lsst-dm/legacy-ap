// -*- lsst-c++ -*-
/** @file
  * @brief Classes that represent clusters of sources and their attributes.
  *
  * @ingroup ap
  * @author Serge Monkewitz
  */
#ifndef LSST_AP_OPTICS_SOURCECLUSTER_H
#define LSST_AP_OPTICS_SOURCECLUSTER_H

#include <vector>

#include "boost/shared_ptr.hpp"

#include "lsst/afw/detection/Source.h"

#include "../Common.h"


namespace lsst { namespace ap { namespace optics {

/** Attributes derived from a cluster of Sources; i.e. Sources that
  * have been determined to be observations of the same physical object
  * according to some algorithm.
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


/** A class that packages up the sources in a cluster along with cluster attributes.
  */
struct LSST_AP_API SourceCluster {
    SourceClusterAttributes::Ptr attributes;
    lsst::afw::detection::SourceSet sources;

    SourceCluster() : attributes(), sources() { }
    ~SourceCluster();
};

typedef std::vector<SourceCluster> SourceClusterVector;

}}} // namespace lsst:ap::optics

#endif // LSST_AP_OPTICS_SOURCECLUSTER_H
