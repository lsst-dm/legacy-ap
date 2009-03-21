// -*- lsst-c++ -*-

/**
 * @file
 * @brief   TimeSpec and Stopwatch utility class implementations.
 *
 * @ingroup ap
 */

#include <sys/time.h> // for ::gettimeofday()
#if LSST_AP_HAVE_CLOCK_GETTIME
#   include <time.h> // for ::clock_gettime()
#endif

#include <ostream>
#include <sstream>

#include <boost/numeric/conversion/converter.hpp>

#include <lsst/ap/Time.h>


namespace lsst {
namespace ap {


// -- TimeSpec ----------------

/// @cond
typedef boost::numeric::converter<
    time_t,
    double,
    boost::numeric::conversion_traits<time_t, double>,
    boost::numeric::def_overflow_handler,
    boost::numeric::Floor<double>
> DoubleToTime;
/// @endcond

TimeSpec::TimeSpec(double const seconds) {
    time_t const sec = DoubleToTime::convert(seconds);
    tv_sec  = sec;
    tv_nsec = static_cast<long>((seconds - static_cast<double>(sec))*1e9);
}


TimeSpec & TimeSpec::operator=(double const seconds) {
    time_t const sec = DoubleToTime::convert(seconds);
    tv_sec  = sec;
    tv_nsec = static_cast<long>((seconds - static_cast<double>(sec))*1e9);
    return *this;
}


TimeSpec & TimeSpec::operator+=(double const seconds) {
    time_t const sec = DoubleToTime::convert(seconds);
    tv_sec  += sec;
    tv_nsec += static_cast<long>((seconds - static_cast<double>(sec))*1e9);
    if (tv_nsec >=  1000000000l) {
        tv_nsec -= 1000000000l;
        ++tv_sec;
    }
    return *this;
}


/**
 * Sets the TimeSpec to the number of seconds and nanoseconds that have elapsed
 * since some arbitrary point in the past.
 */
TimeSpec & TimeSpec::now() {
#if LSST_AP_HAVE_CLOCK_GETTIME
#   if defined(CLOCK_MONOTONIC)
    clock_gettime(CLOCK_MONOTONIC, this);
#   elif defined(CLOCK_HIGHRES)
    clock_gettime(CLOCK_HIGHRES, this);
#   else
    clock_gettime(CLOCK_REALTIME, this); // Use system clock if nothing else is available
#   endif
    return *this;
#else
    return systemTime();
#endif
}


/**
 * Sets the TimeSpec to the number of seconds and nanoseconds that have elapsed
 * since midnight (0 hour), January 1, 1970.
 */
TimeSpec & TimeSpec::systemTime() {
    ::timeval tv;
    ::gettimeofday(&tv, 0);
    tv_sec  = tv.tv_sec;
    tv_nsec = tv.tv_usec*1000l;
    return *this;
}


// -- Stopwatch ----------------

Stopwatch::Stopwatch(bool const go) : _ts(), _stopped(true) {
    if (go) {
        start();
    }
}


void Stopwatch::start() {
    if (_stopped) {
        _ts.now();
        _stopped = false;
    }
}


void Stopwatch::stop() {
    if (!_stopped) {
        TimeSpec t;
        t.now() -= _ts;
        _ts      = t;
        _stopped = true;
    }
}


std::string const Stopwatch::toString() const {
    std::ostringstream os;
    os << *this;
    return os.str();
}


double Stopwatch::seconds() const {
    TimeSpec t;
    if (_stopped) {
        t = _ts;
    } else {
        t.now();
        t -= _ts;
    }
    return t.seconds();
}


std::ostream & operator<<(std::ostream & os, Stopwatch const & watch) {
    TimeSpec t;
    if (watch._stopped) {
        t = watch._ts;
    } else {
        t.now();
        t -= watch._ts;
    }
    long nsec = t.tv_nsec % 1000000000l;
    long sec  = t.tv_sec + t.tv_nsec/1000000000l;
    long min  = (sec/60) % 60;
    long hour = sec/3600;

    double s = static_cast<double>(sec % 60) + 1e-9 * static_cast<double>(nsec);

    if (hour > 0) {
        os << hour << "hr ";
    }
    if (min > 0) {
        os << min << "min ";
    }
    os << s << "sec";
    return os;
}


}} // end of namespace lsst::ap
