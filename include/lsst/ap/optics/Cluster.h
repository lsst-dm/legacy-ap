// -*- lsst-c++ -*-
/** @file
  * @brief Functions for clustering points using the OPTICS algorithm.
  *
  * For details of the algorithm, see the following paper:
  *
  * "OPTICS: Ordering Points To Identify the Clustering Structure".
  * Mihael Ankerst, Markus M. Breunig, Hans-Peter Kriegel, Jörg Sander (1999).
  * ACM SIGMOD international conference on Management of data.
  * ACM Press. pp. 49–60.
  *
  * @ingroup ap
  * @author Serge Monkewitz
  */
#ifndef LSST_AP_OPTICS_CLUSTER_H
#define LSST_AP_OPTICS_CLUSTER_H

#include <vector>

#include "lsst/pex/policy/Policy.h"
#include "lsst/afw/detection/Source.h"


namespace lsst { namespace ap { namespace optics {

LSST_AP_API std::vector<lsst::afw::detection::SourceSet> cluster(
    lsst::afw::detection::SourceSet const & sources,
    lsst::pex::policy::Policy::Ptr policy);

LSST_AP_API std::vector<lsst::afw::detection::SourceSet> cluster(
    std::vector<lsst::afw::detection::SourceSet> const & sources,
    lsst::pex::policy::Policy::Ptr policy);

}}} // namespace lsst:ap::optics

#endif // LSST_AP_OPTICS_CLUSTER_H
