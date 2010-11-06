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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE EarthPostion

#include "boost/test/unit_test.hpp"
#include <cmath>
#include <vector>

#include "lsst/afw/math/Random.h"
#include "lsst/ap/Time.h"
#include "lsst/ap/utils/EarthPosition.h"

// Include IAU SOFA reference implementation for comparison
#include "iausofa/sofam.h"
extern "C" {
    #include "iausofa/epv00.c"
}

using std::sqrt;
using std::vector;

using lsst::afw::math::Random;
using lsst::ap::Stopwatch;
using lsst::ap::utils::earthPosition;


BOOST_AUTO_TEST_CASE(compareToReference) {
    double pvh[2][3];
    double pvb[2][3];
    Random rng;
    vector<double> distKm;
    for (double days= -100.0*DJY; days <= 100.0*DJY; days += 10.0*rng.uniform()) {
        Eigen::Vector3d p = earthPosition(DJM00 + days);
        iauEpv00(DJ00, days, pvh, pvb);
        BOOST_CHECK_CLOSE(pvb[0][0], p.x(), 1.0e-7);
        BOOST_CHECK_CLOSE(pvb[0][1], p.y(), 1.0e-7);
        BOOST_CHECK_CLOSE(pvb[0][2], p.z(), 1.0e-7);
        double d0 = p[0] - pvb[0][0];
        double d1 = p[1] - pvb[0][1];
        double d2 = p[2] - pvb[0][2];
        distKm.push_back(sqrt(d0*d0 + d1*d1 + d2*d2)*DAU/1000.0);
    }
    double meanDistKm = 0.0;
    double maxDistKm = 0.0;
    double distKmRms = 0.0;
    for (vector<double>::const_iterator i = distKm.begin(), e = distKm.end();
         i != e; ++i) {
        meanDistKm += *i;
        if (*i > maxDistKm) {
            maxDistKm = *i;
        }
    }
    meanDistKm /= distKm.size();
    for (vector<double>::const_iterator i = distKm.begin(), e = distKm.end();
         i != e; ++i) {
        distKmRms += (*i - meanDistKm)*(*i - meanDistKm);
    }
    distKmRms = sqrt(distKmRms/distKm.size());
    BOOST_TEST_MESSAGE("Over " << distKm.size() << " epochs from 1900 to 2100:");
    BOOST_TEST_MESSAGE("maximum distance between IAU and LSST barycentric "
                       "earth coordinates: " << maxDistKm << " km");
    BOOST_TEST_MESSAGE("mean distance between IAU and LSST barycentric "
                       "earth coordinates: " << meanDistKm << " km");
    BOOST_TEST_MESSAGE("RMS of distance between IAU and LSST barycentric "
                       "earth coordinates: " << distKmRms << " km");
}

BOOST_AUTO_TEST_CASE(speed) {
    {
        Stopwatch watch(true);
        volatile double x = 0.0;
        volatile double y = 0.0;
        volatile double z = 0.0;
        int n = 0;
        for (double days= -100.0*DJY; days <= 100.0*DJY; days += 10.0, ++n) {
            Eigen::Vector3d p = earthPosition(DJM00 + days);
            x += p[0]; y += p[1]; z += p[2];
        }
        BOOST_TEST_MESSAGE("Computed " << n << " earth positions in " << watch);
    }
    {
        Stopwatch watch(true);
        double pvh[2][3];
        double pvb[2][3];
        volatile double x = 0.0;
        volatile double y = 0.0;
        volatile double z = 0.0;
        int n = 0;
        for (double days= -100.0*DJY; days <= 100.0*DJY; days += 10.0, ++n) {
            iauEpv00(DJ00, days, pvh, pvb);
            x += pvb[0][0]; y += pvb[0][1]; z += pvb[0][2];
        }
        BOOST_TEST_MESSAGE("Computed " << n << " earth positions in " << watch << " [IAU]");
    }
}
