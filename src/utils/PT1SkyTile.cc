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
  * @brief Implementation of the PT1 sky tile class.
  *
  * @ingroup ap
  * @author Serge Monkewitz
  */
#include "lsst/ap/utils/PT1SkyTile.h"

#include <cmath>

#include "lsst/pex/exceptions.h"
#include "lsst/afw/geom/Angle.h"


namespace except = lsst::pex::exceptions;

using lsst::afw::coord::Coord;
using lsst::afw::geom::ONE_OVER_PI;
using lsst::afw::geom::PI;

namespace lsst { namespace ap { namespace utils {

/** Creates a new sky-tile.
  *
  * @param[in] resolution   Quad-sphere resolution - must be at least 3.
  * @param[in] root         Root pixel number - must be in range [0, 6).
  * @param[in] x            X coordinate within root pixel -
  *                         must be in range [0, resolution).
  * @param[in] y            Y coordinate within root pixel -
  *                         must be in range [0, resolution).
  * @param[in] id           A unique integer identifier for the sky-tile.
  */
PT1SkyTile::PT1SkyTile(int resolution, int root, int x, int y, int id) :
    _resolution(resolution),
    _root(root),
    _x(x),
    _y(y),
    _id(id)
{
    if (resolution < 3) {
        throw LSST_EXCEPT(except::InvalidParameterException,
                          "invalid quad-sphere resolution");
    }
    if (root < 0 || root > 5) {
        throw LSST_EXCEPT(except::InvalidParameterException,
                          "quad-sphere root pixel number not in range [0, 6)");
    }
    if (x < 0 || x >= resolution || y < 0 || y >= resolution) {
        throw LSST_EXCEPT(except::InvalidParameterException,
                          "quad sphere sky-pixel coordinates out of bounds");
    }
}

PT1SkyTile::~PT1SkyTile() { }

/** Tests whether a point given in spherical coordinates is inside the sky tile.
  *
  * @param[in] theta    Longitude angle (radians).
  * @param[in] phi      Latitude angle (radians).
  *
  * @return     @c true iff @c (theta,phi) is inside this sky-tile.
  */
bool PT1SkyTile::contains(Coord const & coord) const {
    double theta = coord.getLongitude().asRadians();
    double phi = coord.getLatitude().asRadians();
    int root = static_cast<int>(std::fmod(0.5 + 2.0 * ONE_OVER_PI * theta, 4.0));
    double theta1 = theta - 0.5 * PI * root;
    double tanPhi = std::tan(phi);
    double x, y = tanPhi / std::cos(theta1);
    if (y > 1.0) {
        if (_root != 0) {
            return false;
        }
        x = -std::sin(theta) / tanPhi;
        y = std::cos(theta) / tanPhi;
    } else if (y < -1.0) {
        if (_root != 5) {
            return false;
        }
        x = std::sin(theta) / tanPhi;
        y = std::cos(theta) / tanPhi;
    } else {
        if (_root != root + 1) {
            return false;
        }
        x = std::tan(theta1);
    }
    int ix = static_cast<int>(std::floor(_resolution * 0.5 * (x + 1.0)));
    int iy = static_cast<int>(std::floor(_resolution * 0.5 * (y + 1.0)));
    if (ix >= _resolution) {
        ix = _resolution - 1;
    }
    if (iy >= _resolution) {
        iy = _resolution - 1;
    }
    return _x == ix && _y == iy;
}

}}} // namespace lsst:ap::utils
