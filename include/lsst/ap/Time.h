// -*- lsst-c++ -*-
//
//##====----------------                                ----------------====##/
//
//! \file   Time.h
//! \brief  Convenience wrapper for the C library timespec struct and
//!         a simple profiling class.
//
//##====----------------                                ----------------====##/

#ifndef LSST_AP_TIME_H
#define LSST_AP_TIME_H

#include <time.h>

#include <iosfwd>

#include <lsst/ap/Common.h>


namespace lsst {
namespace ap {


/*! \brief  Wraps the C library timespec struct. */
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


/*! \brief  Utility class for profiling. */
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
