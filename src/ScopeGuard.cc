// -*- lsst-c++ -*-

/**
 * @file
 * @brief   Static members for ScopeGuard.
 *
 * @ingroup associate
 */

#include <lsst/ap/ScopeGuard.h>

namespace lsst {
namespace ap {

int volatile ScopeGuard::_numGuards = 0;

}} // end of namespace lsst::ap

