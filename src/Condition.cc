// -*- lsst-c++ -*-

/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
 * 
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the LSST License Statement and 
 * the GNU General Public License along with this program.  If not, 
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */
 

/**
 * @file
 * @brief   Implementation of POSIX condition variable wrappers.
 *
 * @ingroup ap
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

