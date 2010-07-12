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
 * @brief   Convenience wrapper for the C library timespec struct and a simple profiling class.
 *
 * @ingroup ap
 */

#ifndef LSST_AP_TIME_H
#define LSST_AP_TIME_H

#include <time.h>

#include <iosfwd>

#include "Common.h"


namespace lsst {
namespace ap {


/** @brief  Wraps the C library timespec struct. */
class LSST_AP_API TimeSpec : public ::timespec {

public :

    TimeSpec() {
        tv_sec  = 0;
        tv_nsec = 0;
    }

    TimeSpec(double const seconds);

    TimeSpec & operator+=(TimeSpec const & ts) {
        tv_sec  += ts.tv_sec;
        tv_nsec += ts.tv_nsec;
        if (tv_nsec >= 1000000000l) {
            tv_nsec -= 1000000000l;
            ++tv_sec;
        }
        return *this;
    }

    TimeSpec & operator-=(TimeSpec const & ts) {
        tv_sec  -= ts.tv_sec;
        tv_nsec -= ts.tv_nsec;
        if (tv_nsec < 0l) {
            tv_nsec += 1000000000l;
            --tv_sec;
        }
        return *this;
    }

    TimeSpec & operator=(double const seconds);

    TimeSpec & operator+=(double const seconds);

    TimeSpec & operator-=(double const seconds) {
        return operator+=(-seconds);
    }

    double seconds() const {
        return static_cast<double>(tv_sec) + static_cast<double>(tv_nsec)/1e9;
    }

    TimeSpec & now();

    TimeSpec & systemTime();
};


/** @brief  Utility class for profiling. */
class LSST_AP_API Stopwatch {

public :

    explicit Stopwatch(bool const go = true);

    void start();
    void stop();

    std::string const toString() const;
    double seconds() const;

    friend std::ostream & operator<<(std::ostream &, Stopwatch const &);

private :

    TimeSpec _ts;
    bool     _stopped;
};


}} // end of namespace lsst::ap

#endif // LSST_AP_TIME_H
