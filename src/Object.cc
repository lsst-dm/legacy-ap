// -*- lsst-c++ -*-
//
//##====----------------                                ----------------====##/
//
//! \file   Object.cc
//
//##====----------------                                ----------------====##/

#include <lsst/ap/Object.h>


namespace lsst {
namespace ap {

// -- SimpleObject ----------------

LSST_AP_API bool operator==(SimpleObject const & o1, SimpleObject const & o2) {
    if (o1._objectId == o2._objectId && o1._ra == o2._ra && o1._decl == o2._decl) {
        for (int i = 0; i < lsst::fw::Filter::NUM_FILTERS; ++i) {
            if (o1._varProb[i] != o2._varProb[i]) {
                return false;
            }
        }
        return true;
    }
    return false;
}


}} // end of namespace lsst::ap

