// -*- lsst-c++ -*-

/**
 * @file
 * @brief   Implementation of POSIX condition variable wrappers.
 *
 * @ingroup associate
 */

#include <pthread.h>

#include <boost/bind.hpp>

#include <lsst/ap/Condition.h>
#include <lsst/ap/Exceptions.h>
#include <lsst/ap/Mutex.h>
#include <lsst/ap/ScopeGuard.h>


namespace lsst {
namespace ap {


template <>
Condition<Mutex>::Condition() {
    int err = ::pthread_cond_init(&_condition, 0);
    if (err != 0) {
        LSST_AP_THROW_ERR(Runtime, "pthread_cond_init() failed", err);
    }
}


template <>
Condition<SharedMutex>::Condition() {
    ::pthread_condattr_t attr;
    int err = ::pthread_condattr_init(&attr);
    if (err != 0) {
        LSST_AP_THROW_ERR(Runtime, "pthread_condattr_init() failed", err);
    }
    ScopeGuard attrGuard(boost::bind(::pthread_condattr_destroy, &attr));
    ::pthread_condattr_setpshared(&attr, PTHREAD_PROCESS_SHARED);
    err = ::pthread_cond_init(&_condition, &attr);
    if (err != 0) {
        LSST_AP_THROW_ERR(Runtime, "pthread_cond_init() failed", err);
    }
}


}} // end of namespace lsst::ap

