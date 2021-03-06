// -*- lsst-c++ -*-

/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
 * 
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the LSST License Statement and 
 * the GNU General Public License along with this program.  If not, 
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */
 
/** @file
  * @brief Image utility methods.
  *
  * @ingroup ap
  * @author Serge Monkewitz
  */
#ifndef LSST_AP_UTILS_IMAGEUTILS_H
#define LSST_AP_UTILS_IMAGEUTILS_H

#include <vector>

#include "lsst/afw/geom/Point.h"
#include "lsst/afw/table/Source.h"
#include "lsst/afw/image/Image.h"
#include "lsst/afw/image/Wcs.h"


namespace lsst { namespace ap { namespace utils {

void makeSourceHistogram(
    lsst::afw::image::Image<unsigned short>::Ptr histogram,
    lsst::afw::table::SourceCatalog const & sources,
    lsst::afw::image::Wcs::Ptr wcs,
    bool ignoreOffImage);

void rasterizePolygon(
    std::vector<lsst::afw::geom::Point2D> const &verts,
    lsst::afw::image::Image<float>::Ptr img);

void updateCoverageMap(
    lsst::afw::image::Image<float>::Ptr covMap,
    lsst::afw::image::Wcs::Ptr covMapWcs,
    lsst::afw::image::Wcs::Ptr wcs,
    int width,
    int height,
    int step);

}}} // namespace lsst:ap::utils

#endif // LSST_AP_UTILS_IMAGEUTILS_H
