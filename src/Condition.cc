// -*- lsst-c++ -*-

/**
 * @file
 * @brief   Implementation of POSIX condition variable wrappers.
 *
 * @ingroup associate
 */

#include <pthread.h>

#include "boost/bind.hpp"
#include "boost/format.hpp"

#include "lsst/pex/exceptions.h"

#include "lsst/ap/Condition.h"
#include "lsst/ap/Mutex.h"
#include "lsst/ap/ScopeGuard.h"

namespace ex = lsst::pex::exceptions;

namespace lsst { namespace ap {


template <>
Condition<Mutex>::Condition() {
    int err = ::pthread_cond_init(&_condition, 0);
    if (err != 0) {
        throw LSST_EXCEPT(ex::RuntimeErrorException,
            (boost::format("pthread_cond_init() failed, return code: %1%") % err).str());
    }
}


template <>
Condition<SharedMutex>::Condition() {
    ::pthread_condattr_t attr;
    int err = ::pthread_condattr_init(&attr);
    if (err != 0) {
        throw LSST_EXCEPT(ex::RuntimeErrorException,
            (boost::format("pthread_condattr_init() failed, return code: %1%") % err).str());
    }
    ScopeGuard attrGuard(boost::bind(::pthread_condattr_destroy, &attr));
    ::pthread_condattr_setpshared(&attr, PTHREAD_PROCESS_SHARED);
    err = ::pthread_cond_init(&_condition, &attr);
    if (err != 0) {
        throw LSST_EXCEPT(ex::RuntimeErrorException,
            (boost::format("pthread_cond_init() failed, return code: %1%") % err).str());
    }
}

}} // end of namespace lsst::ap

