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
 

/**
 * @file
 * @brief   Implementation of Point class.
 *
 * @ingroup ap
 */

#include <cassert>
#include <cmath>

#include <iostream>

#include "lsst/ap/Common.h"
#include "lsst/ap/Point.h"
#include "lsst/ap/SpatialUtil.h"
#include "lsst/ap/Time.h"
#include "lsst/afw/geom/Angle.h"

using lsst::afw::math::Random;
namespace afwGeom = lsst::afw::geom;

namespace lsst { namespace ap { namespace {

double randomDec(Random & rng, double const decMin, double const decMax) {

    assert(decMin < decMax && decMin < 90.0 && decMax > -90.0);

    double min  = (decMin < -90.0) ? -90.0 : decMin;
    double max  = (decMax >  90.0) ?  90.0 : decMax;
    double z    = rng.flat(std::sin(afwGeom::degToRad(min)), std::sin(afwGeom::degToRad(max)));
    double res  = afwGeom::radToDeg(std::asin(z));
    if (res < decMin) {
        return decMin;
    } else if (res > decMax) {
        return decMax;
    }
    return res;
}

}}} // end of namespace lsst::ap::<anonymous>


/**
 * Randomly perturbs the point such that the results are distributed according to a normal
 * distribution centered on the original point and having a standard deviation of @a sigma
 * degrees.
 */
lsst::ap::Point & lsst::ap::Point::perturb(Random & rng, double const sigma) {
    return perturb(rng, sigma, rng.uniform()*360.0);
}


/**
 * Randomly perturbs the point in the direction given by the specified position angle so that the
 * distance to the original point is normally distributed with a standard deviation of @a sigma degrees.
 */
lsst::ap::Point & lsst::ap::Point::perturb(Random & rng, double const sigma, double const pa) {

    double sra  = std::sin(afwGeom::degToRad(_ra));
    double cra  = std::cos(afwGeom::degToRad(_ra));
    double sdec = std::sin(afwGeom::degToRad(_dec));
    double cdec = std::cos(afwGeom::degToRad(_dec));

    // original position p
    double x = cra*cdec;
    double y = sra*cdec;
    double z = sdec;

    double spa = std::sin(afwGeom::degToRad(pa));
    double cpa = std::cos(afwGeom::degToRad(pa));

    // north vector tangential to p
    double nx = - cra*sdec;
    double ny = - sra*sdec;
    double nz =   cdec;

    // east vector tangential to p
    double ex = - sra;
    double ey =   cra;
    double ez =   0.0;

    // rotate north vector at V by minus position angle
    double tx = spa*ex + cpa*nx;
    double ty = spa*ey + cpa*ny;
    double tz = spa*ez + cpa*nz;

    // perturb in this direction by a random angle that is normally
    // distributed with a standard deviation of sigma degrees
    double mag  = afwGeom::degToRad(rng.gaussian()*sigma);
    double smag = std::sin(mag);
    double cmag = std::cos(mag);

    // obtain the perturbed position
    x = x*cmag + tx*smag;
    y = y*cmag + ty*smag;
    z = z*cmag + tz*smag;
    // finally, convert back to spherical coordinates (in degrees)
    _ra = afwGeom::radToDeg(std::atan2(y, x));
    if (_ra < 0.0) {
        _ra += 360.0;
    }
    _dec = afwGeom::radToDeg(std::asin(z));
    if (_dec <= -90.0) {
        _dec = -90.0;
    } else if (_dec >= 90.0) {
        _dec = 90.0;
    }
    return *this;
}


/** Returns the angular distance to the given point (in degrees). */
double lsst::ap::Point::distance(Point const & p) const {

    double sra  = std::sin(afwGeom::degToRad(_ra));
    double cra  = std::cos(afwGeom::degToRad(_ra));
    double sdec = std::sin(afwGeom::degToRad(_dec));
    double cdec = std::cos(afwGeom::degToRad(_dec));

    double x = cra*cdec;
    double y = sra*cdec;
    double z = sdec;

    sra  = std::sin(afwGeom::degToRad(p._ra));
    cra  = std::cos(afwGeom::degToRad(p._ra));
    sdec = std::sin(afwGeom::degToRad(p._dec));
    cdec = std::cos(afwGeom::degToRad(p._dec));

    x *= cra*cdec;
    y *= sra*cdec;
    z *= sdec;

    return afwGeom::radToDeg(std::acos(x + y + z));
}


/** Picks a point uniformly at random on the unit sphere. */
lsst::ap::Point const lsst::ap::Point::random(Random & rng) {
    double z = rng.flat(-1.0, 1.0);
    return Point(rng.flat(0.0, 360.0), afwGeom::radToDeg(std::asin(z)));
}


/** Picks a point uniformly at random in the specified dec band. */
lsst::ap::Point const lsst::ap::Point::random(
    Random & rng,
    double const decMin,
    double const decMax
) {
    return Point(rng.flat(0.0, 360.0), randomDec(rng, decMin, decMax));
}


/** Picks a point uniformly at random in the specified box. */
lsst::ap::Point const lsst::ap::Point::random(
    Random & rng, 
    double const raMin,
    double const raMax,
    double const decMin,
    double const decMax
) {
    assert(raMin >= 0.0 && raMin <= 360.0);
    assert(raMax >= 0.0 && raMax <= 360.0);

    double ra;
    if (raMin < raMax) {
        ra = rng.flat(raMin, raMin);
        if (ra > raMax) {
            ra = raMax;
        }
    } else {
        // wrap-around
        double m = raMin - 360.0;
        ra = rng.flat(m, raMax);
        if (ra < 0) {
            ra += 360.0;
        }
    }
    return Point(ra, randomDec(rng, decMin, decMax));
}

