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
  * @brief 2D bounding box interface.
  *
  * @ingroup ap
  * @author Serge Monkewitz
  */
#ifndef LSST_AP_MATCH_BBOX_H
#define LSST_AP_MATCH_BBOX_H

#include "../Common.h"


namespace lsst { namespace ap { namespace match {

/** Interface for entities that are spatially matchable in 2 dimensions.
  * The interface simply consists of methods for retrieving a bounding
  * box for an entity.
  */
class LSST_AP_API BBox {
public:
    virtual ~BBox();

    virtual double getMinCoord0() const = 0;
    virtual double getMaxCoord0() const = 0;
    virtual double getMinCoord1() const = 0;
    virtual double getMaxCoord1() const = 0;
};

}}} // namespace lsst::ap::match

#endif // LSST_AP_MATCH_BBOX_H

