// -*- lsst-c++ -*-
//
//##====----------------                                ----------------====##/
//
//! \file   Mutex.cc
//
//##====----------------                                ----------------====##/

#include <pthread.h>

#include <stdexcept>

#include <boost/bind.hpp>

#include <lsst/ap/Exceptions.h>
#include <lsst/ap/Mutex.h>
#include <lsst/ap/ScopeGuard.h>


namespace lsst {
namespace ap {


SharedMutex::SharedMutex() {
    ::pthread_mutexattr_t attr;

    int err = ::pthread_mutexattr_init(&attr);
    if (err != 0) {
        LSST_AP_THROW_ERR(Runtime, "pthread_mutexattr_init() failed", err);
    }
    ScopeGuard attrGuard(boost::bind(::pthread_mutexattr_destroy, &attr));
    ::pthread_mutexattr_setpshared(&attr, PTHREAD_PROCESS_SHARED);
    err = ::pthread_mutex_init(&_mutex, &attr);
    if (err != 0) {
        LSST_AP_THROW_ERR(Runtime, "pthread_mutex_init() failed", err);
    }
}


}} // end of namespace lsst::ap

