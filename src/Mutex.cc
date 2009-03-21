// -*- lsst-c++ -*-

/**
 * @file
 * @brief   Implementation of the Mutex and SharedMutex classes
 *
 * @ingroup ap
 */

#include <pthread.h>

#include "boost/bind.hpp"
#include "boost/format.hpp"

#include "lsst/pex/exceptions.h"

#include "lsst/ap/Mutex.h"
#include "lsst/ap/ScopeGuard.h"

namespace ex = lsst::pex::exceptions;


lsst::ap::SharedMutex::SharedMutex() {
    ::pthread_mutexattr_t attr;

    int err = ::pthread_mutexattr_init(&attr);
    if (err != 0) {
        throw LSST_EXCEPT(ex::RuntimeErrorException,
            (boost::format("pthread_mutexattr_init() failed, return code: %1%") % err).str());
    }
    ScopeGuard attrGuard(boost::bind(::pthread_mutexattr_destroy, &attr));
    ::pthread_mutexattr_setpshared(&attr, PTHREAD_PROCESS_SHARED);
    err = ::pthread_mutex_init(&_mutex, &attr);
    if (err != 0) {
        throw LSST_EXCEPT(ex::RuntimeErrorException,
            (boost::format("pthread_mutex_init() failed, return code: %1%") % err).str());
    }
}

