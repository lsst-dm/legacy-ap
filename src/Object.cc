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
 * @brief   Implementation of Object class.
 *
 * @ingroup ap
 */

#include "lsst/ap/Object.h"


bool lsst::ap::operator==(lsst::ap::Object const & o1, lsst::ap::Object const & o2) {
    if (o1._objectId == o2._objectId && o1._ra == o2._ra && o1._decl == o2._decl) {
        for (int i = 0; i < lsst::ap::Object::NUM_FILTERS; ++i) {
            if (o1._varProb[i] != o2._varProb[i]) {
                return false;
            }
        }
        return true;
    }
    return false;
}

