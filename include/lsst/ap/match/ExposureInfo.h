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
#ifndef LSST_AP_MATCH_EXPOSUREINFO_H
#define LSST_AP_MATCH_EXPOSUREINFO_H

#include <utility>

#include "Eigen/Core"

#include "lsst/afw/geom/Extent.h"
#include "lsst/afw/image/Calib.h"
#include "lsst/afw/image/Wcs.h"

#include "../Common.h"
#include "../utils/EarthPosition.h"
#include "BBox.h"


namespace lsst { namespace ap { namespace match {

/** Class that bundles together the WCS, extents, time and flux calibration
  * information from an exposure (typically a CCD). No pixel access is provided.
  * The lsst::daf::base::PropertySet containing the FITS header cards from
  * which instances are created is not stored.
  *
  * This provides a relatively memory efficient representation of the CCD
  * metadata AP cares about, which is important since metadata for tens of
  * thousands of CCDs may need to be kept in memory simultaneously.
  */
class LSST_AP_API ExposureInfo : public BBox {
public:
     typedef boost::shared_ptr<ExposureInfo> Ptr;
     typedef boost::shared_ptr<ExposureInfo const> ConstPtr;

     ExposureInfo(lsst::daf::base::PropertySet::Ptr props,
                  std::string const &idKey=std::string("scienceCcdExposureId"));
     ~ExposureInfo();

     /** Returns a unique integer identifier for the exposure.
       */
     inline int64_t getId() const { return _id; }

     /** Returns the filter-id of the exposure (0-5).
       */
     inline int getFilterId() const { return _filterId; }

     /** Returns the exposure mid-point, MJD TAI.
       */
     inline double getEpoch() const { return _epoch; }

     /** Returns the exposure time, s.
       */
     inline double getExposureTime() const { return _expTime; }

     /** Returns the ICRS coordinates of the image center (rad).
       */
     inline Eigen::Vector2d const & getCenter() const {
         return _center;
     }

     /** Returns the SSB coordinates of the earth at t = getEpoch().
       */
     inline Eigen::Vector3d const & getEarthPosition() const {
         if (!_epValid) {
             _earthPos = lsst::ap::utils::earthPosition(_epoch);
             _epValid = true;
         }
         return _earthPos;
     }

     /** Gets exposure width and/or height.
       */
     ///@{
     inline int getWidth() const  { return _extent.getX(); }
     inline int getHeight() const { return _extent.getY(); }
     lsst::afw::geom::Extent2I const getExtent() const { return _extent; }
     ///@}

     /** Is there enough information to calibrate fluxes?
       */
     inline bool canCalibrateFlux() const { return _canCalibrateFlux; }

     double calibrateFlux(double flux, double fluxScale) const;

     std::pair<double, double> const calibrateFlux(double flux,
                                                   double fluxSigma,
                                                   double fluxScale) const;

     /** Returns the exposure WCS.
       */
     ///@{
     inline lsst::afw::image::Wcs::ConstPtr getWcs() const {
         return _wcs;
     }
     inline lsst::afw::image::Wcs::Ptr getWcs() {
         return _wcs;
     }
     ///@}

     // BBox API
     virtual double getMinCoord0() const;
     virtual double getMaxCoord0() const;
     virtual double getMinCoord1() const;
     virtual double getMaxCoord1() const;

#ifndef SWIG
     EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#endif

private:
     Eigen::Vector3d const _pixToSky(double x, double y) const;

     Eigen::Vector2d _center;
     double _radius;
     double _alpha;
     mutable Eigen::Vector3d _earthPos;
     int64_t _id;
     double _epoch;
     double _expTime;
     double _fluxMag0;
     double _fluxMag0Sigma;
     lsst::afw::geom::Extent2I _extent;
     lsst::afw::image::Wcs::Ptr _wcs;
     int _filterId;
     bool _canCalibrateFlux;
     mutable bool _epValid;
};


/** A map from exposure ids to ExposureInfo objects.
  */
class LSST_AP_API ExposureInfoMap {
public:
     typedef boost::shared_ptr<ExposureInfoMap> Ptr;
     typedef boost::shared_ptr<ExposureInfoMap const> ConstPtr;

     ExposureInfoMap();
     ~ExposureInfoMap();

     size_t size() const { return _map.size(); }
     bool empty() const { return _map.empty(); }
     bool contains(int64_t id) const { return _map.find(id) != _map.end(); }

     ExposureInfo::Ptr get(int64_t id) {
         Map::const_iterator i = _map.find(id);
         return (i == _map.end()) ? ExposureInfo::Ptr() : i->second;
     }
     ExposureInfo::ConstPtr get(int64_t id) const {
         Map::const_iterator i = _map.find(id);
         return (i == _map.end()) ? ExposureInfo::ConstPtr() : i->second;
     }

     void insert(ExposureInfo::Ptr info);
     void clear();
     bool erase(int64_t id);

private:
     typedef std::tr1::unordered_map<int64_t, ExposureInfo::Ptr> Map;
     Map _map;
};


}}} // namespace lsst::ap::match

#endif // LSST_AP_MATCH_EXPOSUREINFO_H

