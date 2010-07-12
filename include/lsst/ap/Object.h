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
 * @brief   C++ representations of an LSST Object.
 *
 * @ingroup ap
 */

#ifndef LSST_AP_OBJECT_H
#define LSST_AP_OBJECT_H

#include "lsst/daf/base/DateTime.h"

#include "lsst/afw/image/Filter.h"

#include "Common.h"
#include "Bitset.h"


namespace lsst { namespace ap {

/**
 * @brief   A partial representation of a full LSST Object containing only id,
 *          position, proper motions, and per-filter variability probabilities.
 *
 * This is sufficient for performing spatial crosss matches. None of these fields
 * may ever be NULL.
 */
struct LSST_AP_API Object {
    static int const NUM_FILTERS = 6;
    boost::int64_t _objectId;
    double _ra;
    double _decl;
    double _muRa;
    double _muDecl;
    double _parallax;
    double _radialVelocity;
    boost::int16_t _varProb[NUM_FILTERS];

    boost::int64_t getId() const {
        return _objectId;
    }
    /** @return	Right ascension (degrees) */
    double getRa() const {
        return _ra;
    }
    /** @return Declination (degrees) */
    double getDec() const {
        return _decl;
    }
    /** @return Proper motion (right ascension scaled by cos(dec),  mas/year) */
    double getMuRa() const {
        return _muRa;
    }
    /** @return Proper motion (declination, mas/year) */
    double getMuDecl() const {
        return _muDecl;
    }
    /** @return Parallax (mas) */
    double getParallax() const {
        return _parallax;
    }
    /** @return Radial velocity (km/s) */
    double getRadialVelocity() const {
        return _radialVelocity;
    }
    /**
     * Returns the epoch of the position (getRa(), getDec()). New objects
     * created from difference sources are created with zero proper motion
     * so returning an incorrect epoch for those positions has no effect on
     * association. Therefore, this method can always returns J2000 in MJD(TAI)
     * and there is no need to store epochs in chunk files: object positions
     * measured by deep detection are assumed to be in equatorial
     * J2000.0 coordinates at epoch = J2000.
     */
    double getEpoch() const {
        // 51544.5 - 32.184/86400; this is JD 2451545.0 (TT) converted to MJD(TAI)
        return 51544.4996275;
    }

    boost::int16_t getVarProb(lsst::afw::image::Filter const & f) const {
        return _varProb[f.getId()];
    }
};

LSST_AP_API bool operator==(Object const & o1, Object const & o2);

inline bool operator!=(Object const & o1, Object const & o2) {
    return !(o1 == o2);
}

}}  // end of namespace lsst::ap

#endif // LSST_AP_OBJECT_H
