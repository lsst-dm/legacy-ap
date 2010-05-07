// -*- lsst-c++ -*-
/** @file
  * @brief Utility methods.
  *
  * @ingroup ap
  * @author Serge Monkewitz
  */
#ifndef LSST_AP_CLUSTER_UTILS_H
#define LSST_AP_CLUSTER_UTILS_H

#include "lsst/afw/detection/Source.h"
#include "lsst/afw/image/Image.h"
#include "lsst/afw/image/Wcs.h"

#include "../Common.h"


namespace lsst { namespace ap { namespace cluster {

LSST_AP_API void makeSourceHistogram(
    lsst::afw::detection::SourceSet const & sources,
    lsst::afw::image::Image<unsigned short>::Ptr histogram,
    lsst::afw::image::Wcs::Ptr wcs,
    bool ignoreOffImage);

}}} // namespace lsst:ap::cluster

#endif // LSST_AP_CLUSTER_UTILS_H
