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
  * @brief Class representing a PT1 sky tile.
  *
  * @ingroup ap
  * @author Serge Monkewitz
  */
#ifndef LSST_AP_UTILS_PT1SKYTILE_H
#define LSST_AP_UTILS_PT1SKYTILE_H

#include "lsst/afw/geom/Angle.h"
#include "lsst/afw/detection/Source.h"

#include "../Common.h"


namespace lsst { namespace ap { namespace utils {

/** A PT1 sky-tile is a single pixel obtained from the quadrilateralized
  * spherical cube. Post-PT1, this code will be generalized and moved
  * to lsst::skypix. For now, this small class provides the minimal
  * functionality required by the PT1 SourceAssoc pipeline.
  */
class PT1SkyTile {
public:
    PT1SkyTile(int resolution, int root, int x, int y, int id);
    ~PT1SkyTile();

    int getId() const { return _id; }

    bool contains(lsst::afw::geom::Angle theta, lsst::afw::geom::Angle phi) const;

    void prune(lsst::afw::detection::SourceSet & sources) const;

private:
    int _resolution;
    int _root;
    int _x;
    int _y;
    int _id;
};

}}} // namespace lsst:ap::cluster

#endif // LSST_AP_UTILS_PT1SKYTILE_H
