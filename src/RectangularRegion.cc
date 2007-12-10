// -*- lsst-c++ -*-
//
//##====----------------                                ----------------====##/
//
//! \file   RectangularRegion.cc
//
//##====----------------                                ----------------====##/

#include <stdexcept>

#include <lsst/ap/Exceptions.h>
#include <lsst/ap/CircularRegion.h>
#include <lsst/ap/RectangularRegion.h>
#include <lsst/ap/SpatialUtil.h>


namespace lsst {
namespace ap {


RectangularRegion::RectangularRegion(
    double const minRa,
    double const maxRa,
    double const minDec,
    double const maxDec
) :
    _minRa(minRa),
    _maxRa(maxRa),
    _minDec(minDec),
    _maxDec(maxDec)
{
    if (minRa < 0.0 || minRa >= 360.0 || maxRa < 0.0 || maxRa >= 360.0) {
        LSST_AP_THROW(OutOfRange, "right ascension must be in range [0, 360) degrees");
    }
    if (minDec < -90.0 || minDec > 90.0 || maxDec < -90.0 || maxDec > 90.0) {
        LSST_AP_THROW(OutOfRange, "declination must be in range [-90, 90] degrees");
    }
    if (maxDec < minDec) {
        LSST_AP_THROW(InvalidParameter, "minimum declination greater than maximum declination");
    }
}


RectangularRegion::RectangularRegion(
    double const centerRa,
    double const centerDec,
    double const radius
) {
    fromCircle(centerRa, centerDec, radius);
}


RectangularRegion::RectangularRegion(CircularRegion const & region) {
    fromCircle(region.getCenterRa(), region.getCenterDec(), region.getRadius());
}


void RectangularRegion::fromCircle(double const ra, double const dec, double const radius) {
    if (ra < 0.0 || ra >= 360.0) {
        LSST_AP_THROW(OutOfRange, "right ascension must be in range [0, 360) degrees");
    }
    if (dec < -90.0 || dec > 90.0) {
        LSST_AP_THROW(OutOfRange, "declination must be in range  [-90, 90] degrees");
    }
    if (radius < 0.0 || radius > 90.0) {
        LSST_AP_THROW(OutOfRange, "circle radius must be in range  [0, 90] degrees");
    }
    double alpha = maxAlpha(radius, dec);
    _minRa = ra - alpha;
    if (_minRa < 0.0) {
        _minRa += 360.0;
    }
    _maxRa = ra + alpha;
    if (_maxRa >= 360.0) {
        _maxRa -= 360.0;
    }
    _minDec = dec - radius;
    if (_minDec < -90.0) {
        _minDec = -90.0;
    }
    _maxDec = dec + radius;
    if (_maxDec > 90.0) {
        _maxDec = 90.0;
    }
}


}}  // end of namespace lsst::ap

