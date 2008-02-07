// -*- lsst-c++ -*-

/**
 * @file
 * @brief   Methods for random number generation and randomly perturbing points,
 *          based on Mersenne-Twister code by T. Nishimura and M. Matsumoto.
 *
 * @ingroup associate
 */

#ifndef LSST_AP_RANDOM_H
#define LSST_AP_RANDOM_H


namespace lsst {
namespace ap {


LSST_AP_API void initRandom();

LSST_AP_API double uniformRandom();

LSST_AP_API double uniformRandom1();

LSST_AP_API double normalRandom();

LSST_AP_API bool coinToss(double const p);


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

    Point & perturb(double const sigma);

    Point & perturb(double const sigma, double pa);

    double distance(Point const & p) const;

    static Point const random();

    static Point const random(double const decMin, double const decMax);

    static Point const random(
        double const raMin,
        double const raMax,
        double const decMin,
        double const decMax
    );
};


}} // end of namespace lsst::ap

#endif // LSST_AP_RANDOM_H
