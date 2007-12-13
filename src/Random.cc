// -*- lsst-c++ -*-
//
//##====----------------                                ----------------====##/
//
//! \file   Random.cc
//! \brief  Implementation of random number generation and point
//!         generation/perturbation.
//!
//! This code is nearly identical to code by Makoto Matsumoto and Takuji Nishimura, which
//! carries the following license:
//!
//! <pre>
//! Copyright (C) 2004, Makoto Matsumoto and Takuji Nishimura,
//! All rights reserved.
//!
//! Redistribution and use in source and binary forms, with or without
//! modification, are permitted provided that the following conditions
//! are met:
//!
//!   1. Redistributions of source code must retain the above copyright
//!      notice, this list of conditions and the following disclaimer.
//!
//!   2. Redistributions in binary form must reproduce the above copyright
//!      notice, this list of conditions and the following disclaimer in the
//!      documentation and/or other materials provided with the distribution.
//!
//!   3. The names of its contributors may not be used to endorse or promote
//!      products derived from this software without specific prior written
//!      permission.
//!
//! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
//! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
//! A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
//! CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
//! EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
//! PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
//! PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
//! LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
//! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//!
//! References:
//! T. Nishimura, ``Tables of 64-bit Mersenne Twisters''
//!   ACM Transactions on Modeling and
//!   Computer Simulation 10. (2000) 348--357.
//! M. Matsumoto and T. Nishimura,
//!   ``Mersenne Twister: a 623-dimensionally equidistributed
//!     uniform pseudorandom number generator''
//!   ACM Transactions on Modeling and
//!   Computer Simulation 8. (Jan. 1998) 3--30.
//!
//! Any feedback is very welcome.
//! http://www.math.hiroshima-u.ac.jp/~m-mat/MT/emt.html
//! email: m-mat @ math.sci.hiroshima-u.ac.jp (remove spaces)
//! </pre>
//
//##====----------------                                ----------------====##/

#include <cassert>
#include <cmath>

#include <iostream>

#include <lsst/ap/Common.h>
#include <lsst/ap/Random.h>
#include <lsst/ap/SpatialUtil.h>
#include <lsst/ap/Time.h>


namespace lsst {
namespace ap {

namespace {

#define NN       312
#define MM       156
#define MATRIX_A UINT64_C(0xB5026F5AA96619E9)
#define UM       UINT64_C(0xFFFFFFFF80000000) // Most significant 33 bits
#define LM       UINT64_C(0x7FFFFFFF)         // Least significant 31 bits

// The array for the state vector
static uint64_t mt[NN];

// mti==NN+1 means mt[NN] is not initialized
static int mti = NN + 1;


// Initializes mt[NN] with a seed
static void initRand64(uint64_t seed) {
    mt[0] = seed;
    for (mti = 1; mti < NN; mti++) {
        mt[mti] = UINT64_C(6364136223846793005) * (mt[mti-1] ^ (mt[mti-1] >> 62)) + mti;
    }
}


// Initializes with an array of seeds
static void initRand64(uint64_t * keys, uint64_t nkeys) {
    uint64_t i, j, k;
    initRand64(UINT64_C(19650218));
    i = 1;
    j = 0;
    k = (NN > nkeys ? NN : nkeys);
    for (; k; k--) {
        // non linear
        mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 62)) * UINT64_C(3935559000370003845))) + keys[j] + j;
        i++;
        j++;
        if (i >= NN) {
            mt[0] = mt[NN - 1];
            i = 1;
        }
        if (j >= nkeys) {
            j = 0;
        }
    }
    for (k = NN - 1; k; k--) {
        // non linear
        mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 62)) * UINT64_C(2862933555777941757))) - i;
        i++;
        if (i >= NN){
            mt[0] = mt[NN - 1];
            i = 1;
        }
    }

    mt[0] = UINT64_C(1) << 63; // MSB is 1; assuring non-zero initial array
}


// Generates a random number in the interval [0, 2^64 - 1]
static uint64_t genRand64() {
    static uint64_t const mag01[2] = { 0, MATRIX_A };

    int i;
    uint64_t x;

    if (mti >= NN) {
        // generate NN words at a time

        // if initRand64() has not been called
        // a default initial seed is used
        if (mti == NN + 1) {
            initRand64(UINT64_C(5489));
        }
        for (i = 0; i < NN - MM; i++) {
            x = (mt[i] & UM) | (mt[i + 1] & LM);
            mt[i] = mt[i + MM] ^ (x >> 1) ^ mag01[static_cast<int>(x & UINT64_C(1))];
        }
        for (;i<NN-1;i++) {
            x = (mt[i] & UM) | (mt[i + 1] & LM);
            mt[i] = mt[i + (MM - NN)] ^ (x >> 1) ^ mag01[static_cast<int>(x & UINT64_C(1))];
        }
        x = (mt[NN - 1] & UM) | (mt[0] & LM);
        mt[NN - 1] = mt[MM - 1] ^ (x >> 1) ^ mag01[static_cast<int>(x & UINT64_C(1))];

        mti = 0;
    }

    x = mt[mti++];

    x ^= (x >> 29) & UINT64_C(0x5555555555555555);
    x ^= (x << 17) & UINT64_C(0x71D67FFFEDA60000);
    x ^= (x << 37) & UINT64_C(0xFFF7EEE000000000);
    x ^= (x >> 43);

    return x;
}


static void initRand(TimeSpec const & ts) {
    uint64_t init[4] = { UINT64_C(0x12345), UINT64_C(0x23456), UINT64_C(0x34567), UINT64_C(0x45678) };
    init[0] += static_cast<uint64_t>(ts.tv_sec);
    init[2] += static_cast<uint64_t>(ts.tv_nsec);
    initRand64(init, 4);
}


} // end of anonymous namespace


/*!
    Initializes the random number generation functions with the current system time. Calls subsequent
    to the first initRandom() call will have no effect. N.b. random number generation routines are 
    \b not thread safe!
 */
LSST_AP_API void initRandom() {
    static bool called = false;
    if (!called) {
        TimeSpec ts;
#if !defined(LSST_AP_INIT_RANDOM_SEC) || !defined(LSST_AP_INIT_RANDOM_NSEC)
        ts.systemTime();
        std::clog << "\n"
            << "     /\n"
            << "    | Note: Random number generator initialized. Pass\n"
            << "    |    -DLSST_AP_INIT_RANDOM_SEC=" << ts.tv_sec  << "\n"
            << "    |    -DLSST_AP_INIT_RANDOM_NSEC=" << ts.tv_nsec << "\n"
            << "    | when compiling the association pipeline\n"
            << "    | (see SConstruct) to reproduce results\n"
            << "     \\\n" << std::endl;
#else
        ts.tv_sec  = LSST_AP_INIT_RANDOM_SEC;
        ts.tv_nsec = LSST_AP_INIT_RANDOM_NSEC;
#endif
        initRand(ts);
        called = true;
    }
}


/*! Returns uniformly distributed random numbers in the range [0,1). */
LSST_AP_API double uniformRandom() {
    return (genRand64() >> 11) * (1.0/9007199254740992.0);
}


/*! Returns uniformly distributed random numbers in the range [0,1]. */
LSST_AP_API double uniformRandom1() {
    return (genRand64() >> 11) * (1.0/9007199254740991.0);
}


/*! Returns normally distributed random numbers (with a standard deviation of 1). */
LSST_AP_API double normalRandom() {
    // Use Box-Muller transform to get a normal distribution
    // from the uniform distribution
    double x, y, s;

    do {
        x = 2.0*uniformRandom1() - 1.0;
        y = 2.0*uniformRandom1() - 1.0;
        s = x*x + y*y;
    } while (s > 1.0 || s == 0.0);
    // Even though s == 0.0 is extremely unlikely, check for it and
    // discard (since it will produce garbage below)
    return x * std::sqrt((- 2.0*std::log(s)) / s);
}


/*! Performs a weighted coin toss: returns \c true with probability \a p and \c false
    with probability 1 - \a p. */
LSST_AP_API bool coinToss(double const p) {
    return (uniformRandom() <= p);
}


// -- Point ----------------

namespace {

static double randomDec(double const decMin, double const decMax) {

    assert(decMin < decMax && decMin < 90.0 && decMax > -90.0);

    double min  = (decMin < -90.0) ? -90.0 : decMin;
    double max  = (decMax >  90.0) ?  90.0 : decMax;
    double smin = std::sin(radians(min));
    double smax = std::sin(radians(max));
    double z    = smin + uniformRandom1() * (smax - smin);

    double res  = degrees(std::asin(z));
    if (res < decMin) {
        return decMin;
    } else if (res > decMax) {
        return decMax;
    }
    return res;
}

}


/*!
    Randomly perturbs the point such that the results are distributed according to a normal
    distribution centered on the original point and having a standard deviation of \a sigma
    degrees.
 */
Point & Point::perturb(double const sigma) {
    return perturb(sigma, uniformRandom()*360.0);
}


/*!
    Randomly perturbs the point in the direction given by the specified position angle so that the
    distance to the original point is normally distributed with a standard deviation of \a sigma degrees.
 */
Point & Point::perturb(double const sigma, double const pa) {

    double sra  = std::sin(radians(_ra));
    double cra  = std::cos(radians(_ra));
    double sdec = std::sin(radians(_dec));
    double cdec = std::cos(radians(_dec));

    // original position p
    double x = cra*cdec;
    double y = sra*cdec;
    double z = sdec;

    double spa = std::sin(radians(pa));
    double cpa = std::cos(radians(pa));

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
    double mag  = radians(sigma * normalRandom());
    double smag = std::sin(mag);
    double cmag = std::cos(mag);

    // obtain the perturbed position
    x = x*cmag + tx*smag;
    y = y*cmag + ty*smag;
    z = z*cmag + tz*smag;
    // finally, convert back to spherical coordinates (in degrees)
    _ra = degrees(std::atan2(y, x));
    if (_ra < 0.0) {
        _ra += 360.0;
    }
    _dec = degrees(std::asin(z));
    if (_dec <= -90.0) {
        _dec = -90.0;
    } else if (_dec >= 90.0) {
        _dec = 90.0;
    }
    return *this;
}


/*! Returns the angular distance to the given point (in degrees). */
double Point::distance(Point const & p) const {

    double sra  = std::sin(radians(_ra));
    double cra  = std::cos(radians(_ra));
    double sdec = std::sin(radians(_dec));
    double cdec = std::cos(radians(_dec));

    double x = cra*cdec;
    double y = sra*cdec;
    double z = sdec;

    sra  = std::sin(radians(p._ra));
    cra  = std::cos(radians(p._ra));
    sdec = std::sin(radians(p._dec));
    cdec = std::cos(radians(p._dec));

    x *= cra*cdec;
    y *= sra*cdec;
    z *= sdec;

    return degrees(std::acos(x + y + z));
}


/*! Picks a point uniformly at random on the unit sphere.*/
Point const Point::random() {
    double z   = -1.0 + 2.0*uniformRandom1();
    double dec = degrees(std::asin(z));
    double ra  = uniformRandom1()*360.0;
    return Point(ra, dec);
}


/*! Picks a point uniformly at random in the specified dec band. */
Point const Point::random(double const decMin, double const decMax) {
    return Point(uniformRandom1()*360.0, randomDec(decMin, decMax));
}


/*! Picks a point uniformly at random in the specified box. */
Point const Point::random(
    double const raMin,
    double const raMax,
    double const decMin,
    double const decMax
) {
    assert(raMin >= 0.0 && raMin <= 360.0);
    assert(raMax >= 0.0 && raMax <= 360.0);

    double ra;
    if (raMin < raMax) {
        ra = raMin + uniformRandom1()*(raMax - raMin);
        if (ra > raMax) {
            ra = raMax;
        }
    } else {
        // wrap-around
        double m = raMin - 360.0;
        ra = m + uniformRandom1()*(raMax - m);
        if (ra < 0) {
            ra += 360.0;
        }
    }
    return Point(ra, randomDec(decMin, decMax));
}


}} // end of namespace lsst::ap
