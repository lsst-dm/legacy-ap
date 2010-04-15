// -*- lsst-c++ -*-
/** @file
  * @brief Implementation of the PT1 sky tile class.
  *
  * @ingroup ap
  * @author Serge Monkewitz
  */
#include "lsst/ap/cluster/PT1SkyTile.h"

#include <cmath>

#include "lsst/pex/exceptions.h"
#include "lsst/ap/Common.h"


namespace except = lsst::pex::exceptions;

namespace lsst { namespace ap { namespace cluster {

/** Creates a new sky-tile.
  *
  * @param[in] resolution   Quad-sphere resolution - must be at least 3.
  * @param[in] root         Root pixel number - must be in range [0, 6).
  * @param[in] x            X coordinate within root pixel -
  *                         must be in range [0, resolution).
  * @param[in] y            Y coordinate within root pixel -
  *                         must be in range [0, resolution).
  * @param[in] id           A unique integer identifier for the sky-tile.
  */
PT1SkyTile::PT1SkyTile(int resolution, int root, int x, int y, int id) :
    _resolution(resolution),
    _root(root),
    _x(x),
    _y(y),
    _id(id)
{
    if (resolution < 3) {
        throw LSST_EXCEPT(except::InvalidParameterException,
                          "invalid quad-sphere resolution");
    }
    if (root < 0 || root > 5) {
        throw LSST_EXCEPT(except::InvalidParameterException,
                          "quad-sphere root pixel number not in range [0, 6)");
    }
    if (x < 0 || x >= resolution || y < 0 || y >= resolution) {
        throw LSST_EXCEPT(except::InvalidParameterException,
                          "quad sphere sky-pixel coordinates out of bounds");
    }
}

PT1SkyTile::~PT1SkyTile() { }

/** Tests whether a point given in spherical coordinates is inside the sky tile.
  *
  * @param[in] theta    Longitude angle (radians).
  * @param[in] phi      Latitude angle (radians).
  *
  * @return     @c true iff @c (theta,phi) is inside this sky-tile.
  */
bool PT1SkyTile::contains(double theta, double phi) const {
    int root = static_cast<int>(std::fmod(0.5 + 2.0 * ONE_OVER_PI * theta, 4.0));
    double theta1 = theta - 0.5 * M_PI * root;
    double tanPhi = std::tan(phi);
    double x, y = tanPhi / std::cos(theta1);
    if (y > 1.0) {
        if (_root != 0) {
            return false;
        }
        x = -std::sin(theta) / tanPhi;
        y = std::cos(theta) / tanPhi;
    } else if (y < -1.0) {
        if (_root != 5) {
            return false;
        }
        x = std::sin(theta) / tanPhi;
        y = std::cos(theta) / tanPhi;
    } else {
        if (_root != root + 1) {
            return false;
        }
        x = std::tan(theta1);
    }
    int ix = static_cast<int>(std::floor(_resolution * 0.5 * (x + 1.0)));
    int iy = static_cast<int>(std::floor(_resolution * 0.5 * (y + 1.0)));
    if (ix >= _resolution) {
        ix = _resolution - 1;
    }
    if (iy >= _resolution) {
        iy = _resolution - 1;
    }
    return _x == ix && _y == iy;
}

///@name Removing sources falling outside of a sky-tile
//@{
/** Removes sources falling outside of this sky-tile from the input
  * source list.
  *
  * @param[inout] sources   The sources to prune.
  *
  * @return     A @c std::pair<size_t,size_t>, where the first element
  *             of the pair contains the number of input sources falling
  *             inside the sky-tile, and the second contains the total
  *             number of input sources.
  */
std::pair<size_t, size_t> PT1SkyTile::prune(
    lsst::afw::detection::SourceSet & sources) const
{
    size_t n = sources.size();
    size_t j = 0;
    for (size_t i = 0; i < n; ++i) {
        double ra = sources[i]->getRa() * RADIANS_PER_DEGREE;
        double dec = sources[i]->getDec() * RADIANS_PER_DEGREE;
        if (contains(ra, dec)) {
            if (j != i) {
                sources[j] = sources[i];
                sources[i].reset();
            }
            ++j;
        }
    }
    if (j < n) {
        sources.erase(sources.begin() + j, sources.end());
    }
    return std::make_pair(j, n);
}

std::pair<size_t, size_t> PT1SkyTile::prune(
    std::vector<lsst::afw::detection::SourceSet> & sources) const
{
    typedef std::vector<lsst::afw::detection::SourceSet>::iterator SetIter;
    std::pair<size_t, size_t> p(0, 0);
    // prune each SourceSet
    for (SetIter i = sources.begin(), e = sources.end(); i != e; ++i) {
        std::pair<size_t, size_t> t = prune(*i);
        p.first += t.first;
        p.second += t.second;
    }
    // remove empty SourceSets
    size_t n = sources.size();
    size_t j = 0;
    for (size_t i = 0; i < n; ++i) {
        if (sources[i].size() == 0) {
            if (j != i) {
                sources[j] = sources[i];
                sources[i].clear();
            }
            ++j;
        }
    }
    if (j < n) {
        sources.erase(sources.begin() + j, sources.end());
    }
    return p;
}
//@}

}}} // namespace lsst:ap::cluster
