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

#include <cstddef>
#include <utility>

#include "Eigen/Core"

#include "lsst/afw/geom/Angle.h"
#include "lsst/afw/geom/Extent.h"
#include "lsst/afw/coord/Coord.h"
//#include "lsst/afw/image/Calib.h"
#include "lsst/afw/image/Filter.h"
#include "lsst/afw/image/Wcs.h"

#include "../utils/Csv.h"
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
class ExposureInfo : public BBox {
public:
    typedef boost::shared_ptr<ExposureInfo> Ptr;
    typedef boost::shared_ptr<ExposureInfo const> ConstPtr;

    static std::string const DEF_ID_KEY;

    ExposureInfo(lsst::daf::base::PropertySet::Ptr props,
                 std::string const &idKey=DEF_ID_KEY);
    ExposureInfo(lsst::daf::base::PropertySet::Ptr props,
                 int64_t id);
    ~ExposureInfo();

    /// @brief Return the unique integer identifier for the exposure.
    int64_t getId() const { return _id; }

    /// @brief Return the filter of the exposure. UNKNOWN for multi-band
    //         exposures (e.g. chi-squared coadds).
    lsst::afw::image::Filter const & getFilter() const { return _filter; }

    /// Return the exposure mid-point, MJD TAI. NaN for coadds.
    double getEpoch() const { return _epoch; }

    /// @brief Return the exposure time, s. NaN for coadds.
    double getExposureTime() const { return _expTime; }

    /// @brief Return the ICRS coordinates of the image center (rad).
    lsst::afw::coord::IcrsCoord const & getCenter() const {
        return _center;
    }

    /// @brief Return the SSB coordinates of the earth at t = getEpoch().
    Eigen::Vector3d const & getEarthPosition() const {
        if (!_epValid) {
            _earthPos = lsst::ap::utils::earthPosition(_epoch);
            _epValid = true;
        }
        return _earthPos;
    }

    /// @brief Get exposure width and/or height.
    ///@{
    int getWidth() const  { return _extent.getX(); }
    int getHeight() const { return _extent.getY(); }
    lsst::afw::geom::Extent2I const getExtent() const { return _extent; }
    ///@}

    /// @brief Get offset that, when added to parent exposure pixel coordinates,
    ///        transforms to local exposure pixel coordinates (LTV[12]).
    ///@{
    int getOffsetX() const { return _offset.getX(); }
    int getOffsetY() const { return _offset.getY(); }
    lsst::afw::geom::Point2I const getOffset() const { return _offset; }
    ///@}

    /// @brief Is there enough information to calibrate fluxes?
    bool canCalibrateFlux() const { return _canCalibrateFlux; }

    /** Returns a flux value calibrated using the flux of a zero magnitude object
      * associated with this exposure.
      */
    double calibrateFlux(
        double flux,     ///< flux to calibrate, DN
        double fluxScale ///< flux scaling factor, must be \> 0.0
    ) const;

    /** Returns a calibrated flux and its variance. The flux of a zero magnitude
      * object associated with this exposure is used to perform the calibration.
      */
    std::pair<double, double> const calibrateFlux(
        double flux,       ///< flux to calibrate, DN
        double fluxSigma,  ///< standard deviation of @a flux
        double fluxScale   ///< flux scaling factor, must be \> 0.0
    ) const;

    /// @brief Return the exposure WCS.
    ///@{
    lsst::afw::image::Wcs::ConstPtr getWcs() const {
        return _wcs;
    }
    lsst::afw::image::Wcs::Ptr getWcs() {
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
    void _init(lsst::daf::base::PropertySet::Ptr props);

    lsst::afw::coord::IcrsCoord _center;
    lsst::afw::geom::Angle _radius;
    lsst::afw::geom::Angle _alpha;
    mutable Eigen::Vector3d _earthPos;
    int64_t _id;
    double _epoch;
    double _expTime;
    double _fluxMag0;
    double _fluxMag0Sigma;
    lsst::afw::geom::Extent2I _extent;
    lsst::afw::geom::Point2I _offset;
    lsst::afw::image::Wcs::Ptr _wcs;
    lsst::afw::image::Filter _filter;
    int _filterId;
    bool _canCalibrateFlux;
    mutable bool _epValid;
};


/** A map from exposure ids to ExposureInfo objects.
  */
class ExposureInfoMap {
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


/** Reads an exposure metadata key-value CSV file (where metadata keys
  * must have been grouped by exposure id). An ExposureInfo object
  * is created for each input exposure and appended to @a exposures.
  *
  * @param[in,out] exposures  ExposureInfo objects are appended to this vector.
  * @param[in]     csvFile    Metadata table path.
  * @param[in]     control    Metadata table CSV format.
  * @param[in]     idColumn   Name of ID column, e.g. "scienceCcdExposureId".
  */
void readExposureInfos(
    std::vector<ExposureInfo::Ptr> & exposures,
    std::string const & csvFile,
    lsst::ap::utils::CsvControl const &control,
    std::string const & idColumn);

}}} // namespace lsst::ap::match

#endif // LSST_AP_MATCH_EXPOSUREINFO_H

