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
  * @brief  Class that bundles together the WCS, extents, time, and calibration
  *         information from an image (typically a CCD).
  * @author Serge Monkewitz
  */
#ifndef LSST_AP_MATCH_IMAGEINFO_H
#define LSST_AP_MATCH_IMAGEINFO_H

#include "lsst/afw/geom/Extent.h"
#include "lsst/afw/image/Calib.h"
#include "lsst/afw/image/Wcs.h"


namespace lsst { namespace ap { namespace match {

/** Class that bundles together the WCS, extents, time and flux calibration
  * information from an image (typically a CCD). No pixel access is provided.
  * The lsst::daf::base::PropertySet containing the FITS header cards from
  * which instances are created is not stored.
  *
  * This provides a relatively memory efficient representation of the CCD
  * metadata AP cares about, which is important since metadata for tens of
  * thousands of CCDs may need to be kept in memory simultaneously.
  */
class ImageInfo {
public:
     typedef boost::shared_ptr<ImageInfo> Ptr;
     typedef boost::shared_ptr<ImageInfo const> ConstPtr;

     ImageInfo(int64_t id, lsst::daf::base::PropertySet::Ptr props);
     ~ImageInfo();

     // image id
     inline int64_t getId() const { return _id; }

     // image extent
     inline int getWidth() const  { return _extent.getX(); }
     inline int getHeight() const { return _extent.getY(); }
#ifndef SWIG
     lsst::afw::geom::Extent2I const getExtent() const { return _extent; }
#endif
     // image time and flux calibration metadata
     inline bool canCalibrateFlux() const { return _canCalibrateFlux; }
#ifndef SWIG
     inline lsst::afw::image::Calib::ConstPtr getCalib() const {
         return _calib;
     }
#endif
     inline lsst::afw::image::Calib::Ptr getCalib() {
         return _calib;
     }

     // image WCS
#ifndef SWIG
     inline lsst::afw::image::Wcs::ConstPtr getWcs() const {
         return _wcs;
     }
#endif
     inline lsst::afw::image::Wcs::Ptr getWcs() {
         return _wcs;
     }

private:
     int64_t _id;
     lsst::afw::geom::Extent2I _extent;
     lsst::afw::image::Calib::Ptr _calib;
     lsst::afw::image::Wcs::Ptr _wcs;
     bool _canCalibrateFlux;
};

}}} // namespace lsst::ap::match

#endif // LSST_AP_MATCH_IMAGEINFO_H

