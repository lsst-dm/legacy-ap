// -*- lsst-c++ -*-

/**
 * @file
 * @brief   Class for representing points on the sky, with support
 *          for random perturbations.
 *
 * @ingroup associate
 */

#ifndef LSST_AP_POINT_H
#define LSST_AP_POINT_H

#include "lsst/afw/math/Random.h"

namespace lsst { namespace ap {

/**
 * @brief   A point on the unit sphere (sky), specified in spherical polar coordinates.
 *
 * The units of all angles (stored or passed to member functions) are degrees.
 */
struct LSST_AP_API Point {

    double _ra;
    double _dec;

    Point() : _ra(0.0), _dec(0.0) {}
    Point(double const ra, double const dec) : _ra(ra), _dec(dec) {}

    Point & perturb(lsst::afw::math::Random & rng, double const sigma);
    Point & perturb(lsst::afw::math::Random & rng, double const sigma, double pa);

    double distance(Point const & p) const;

    static Point const random(lsst::afw::math::Random & rng);

    static Point const random(
        lsst::afw::math::Random & rng,
        double const decMin,
        double const decMax
    );

    static Point const random(
        lsst::afw::math::Random & rng, 
        double const raMin,
        double const raMax,
        double const decMin,
        double const decMax
    );
};

}} // end of namespace lsst::ap

#endif // LSST_AP_POINT_H
