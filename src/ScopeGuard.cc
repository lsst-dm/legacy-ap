// -*- lsst-c++ -*-
//
//##====----------------                                ----------------====##/
//
//! \file   ScopeGuard.cc
//
//##====----------------                                ----------------====##/

#include <lsst/ap/ScopeGuard.h>

namespace lsst {
namespace ap {

int volatile ScopeGuard::_numGuards = 0;

}} // end of namespace lsst::ap

