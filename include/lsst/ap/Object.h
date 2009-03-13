// -*- lsst-c++ -*-

/**
 * @file
 * @brief   C++ representations of an LSST Object.
 *
 * @ingroup associate
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
    boost::int64_t _objectId;
    double _ra;
    double _decl;
    double _pmRa;
    double _pmDecl;
    boost::int16_t _varProb[lsst::afw::image::Filter::NUM_FILTERS];

    boost::int64_t getId() const {
        return _objectId;
    }
    double getRa() const {
        return _ra;
    }
    double getDec() const {
        return _decl;
    }
    boost::int16_t getVarProb(lsst::afw::image::Filter const f) const {
        return _varProb[f];
    }
};

LSST_AP_API bool operator==(Object const & o1, Object const & o2);

inline bool operator!=(Object const & o1, Object const & o2) {
    return !(o1 == o2);
}

}}  // end of namespace lsst::ap

#endif // LSST_AP_OBJECT_H
