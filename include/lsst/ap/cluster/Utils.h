// -*- lsst-c++ -*-
/** @file
  * @brief Utility methods.
  *
  * @ingroup ap
  * @author Serge Monkewitz
  */
#ifndef LSST_AP_CLUSTER_UTILS_H
#define LSST_AP_CLUSTER_UTILS_H

#include <vector>

#include "lsst/afw/geom/Point.h"
#include "lsst/afw/detection/Source.h"
#include "lsst/afw/image/Image.h"
#include "lsst/afw/image/Wcs.h"

#include "../Common.h"


namespace lsst { namespace ap { namespace cluster {

LSST_AP_API void makeSourceHistogram(
    lsst::afw::image::Image<unsigned short>::Ptr histogram,
    lsst::afw::detection::SourceSet const & sources,
    lsst::afw::image::Wcs::Ptr wcs,
    bool ignoreOffImage);

LSST_AP_API void rasterizePolygon(
    std::vector<lsst::afw::geom::Point2D> const &verts,
    lsst::afw::image::Image<float>::Ptr img);

LSST_AP_API void updateCoverageMap(
    lsst::afw::image::Image<float>::Ptr covMap,
    lsst::afw::image::Wcs::Ptr covMapWcs,
    lsst::afw::image::Wcs::Ptr wcs,
    int width,
    int height,
    bool cornersOnly);

}}} // namespace lsst:ap::cluster

#endif // LSST_AP_CLUSTER_UTILS_H
