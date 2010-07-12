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
 * @brief   Simple wrapper class for POSIX condition variables.
 *
 * @ingroup ap
 */

#ifndef LSST_AP_CONDITION_H
#define LSST_AP_CONDITION_H

#include <cassert>
#include <errno.h>
#include <pthread.h>

#include "boost/noncopyable.hpp"

#include "Mutex.h"
#include "Time.h"


namespace lsst { namespace ap {

/** @brief  Encapsulates a POSIX condition variable. */
template <typename MutexT>
class Condition : private boost::noncopyable {

public :

    typedef ScopedLock<MutexT> Lock;

    Condition();

    ~Condition() {
        int result = ::pthread_cond_destroy(&_condition);
        assert(result == 0);
    }

    /**
     * Waits on the condition variable until the calling thread is woken as a result
     * of another thread calling notify() or notifyAll(). Spurious wakeup can occur.
     *
     * @pre     @a lock has been successfully acquired
     */
    void wait(Lock & lock) {
        assert(lock.isAcquired());
        int result = ::pthread_cond_wait(&_condition, lock.getPosixMutex());
        assert(result == 0);
    }

    /**
     * Waits on the condition variable until the given predicate evaluates to @c true.
     *
     * @pre     @a lock has been successfully acquired
     */
    template <typename P>
    void wait(Lock & lock, P predicate)
    {
        assert(lock.isAcquired());
        while (!predicate()) {
            int result = ::pthread_cond_wait(&_condition, lock.getPosixMutex());
            assert(result == 0);
        }
    }

    /**
     * Waits on this condition variable until either the given deadline expires or the calling
     * thread is woken as a result of another thread calling notify() or notifyAll(). Spurious
     * wakeup can occur.
     *
     * @pre     @a lock has been successfully acquired
     * @return  @c false if the deadline was missed, and @c true otherwise.
     */
    bool wait(Lock & lock, TimeSpec const & ts) {
        assert(lock.isAcquired());
        int result = ::pthread_cond_timedwait(&_condition, lock.getPosixMutex(), &ts);
        if (result == ETIMEDOUT) {
            return false;
        }
        assert(result == 0);
        return true;
    }

    /**
     * Waits on this condition variable until the given predicate evaluates to @c true
     * or the given deadline is missed.
     *
     * @pre     @a lock has been successfully acquired
     * @return  @c true if the predicate became @c true before the deadline expired,
     *          and @c false if the deadline was missed.
     */
    template <typename P>
    bool wait(Lock & lock, P predicate, TimeSpec const & deadline) {
        assert(lock.isAcquired());
        while (!predicate()) {
            int result = ::pthread_cond_timedwait(&_condition, lock.getPosixMutex(), &deadline);
            if (result == ETIMEDOUT) {
                return false;
            }
            assert(result == 0);
        }
        return true;
    }

    /**
     * Wakes up at least one thread waiting on the condition. For predictable scheduling, the
     * mutex associated with the condition should be acquired prior to calling this method.
     */
    void notify() {
        int result = ::pthread_cond_signal(&_condition);
        assert(result == 0);
    }

    /**
     * Wakes up all threads waiting on the condition. For predictable scheduling, the
     * mutex associated with the condition should be acquired prior to calling this method.
     */
    void notifyAll() {
        int result = ::pthread_cond_broadcast(&_condition);
        assert(result == 0);
    }

private :

    ::pthread_cond_t _condition;
};


}} // end of namespace lsst::ap

#endif // LSST_AP_CONDITION_H
