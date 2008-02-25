// -*- lsst-c++ -*-

/**
 * @file
 * @brief   Simple wrappers for POSIX mutual exclusion locks.
 *
 * @ingroup associate
 */

#ifndef LSST_AP_MUTEX_H
#define LSST_AP_MUTEX_H

#include <pthread.h>
#include <errno.h>

#include <cassert>

#include <boost/noncopyable.hpp>

#include "Common.h"
#include "Exceptions.h"


namespace lsst {
namespace ap {


template <typename MutexT> class ScopedLock;
template <typename MutexT> class Condition;


/** @brief A wrapper for a process private POSIX mutual exclusion lock. */
class LSST_AP_LOCAL Mutex : private boost::noncopyable {

public :

    Mutex() {
        int err = ::pthread_mutex_init(&_mutex, 0);
        if (err != 0) {
            LSST_AP_THROW_ERR(Runtime, "pthread_mutex_init() failed", err);
        }
    }

    ~Mutex() {
        int result = ::pthread_mutex_destroy(&_mutex);
        assert(result == 0);
    }

private :

    ::pthread_mutex_t _mutex;

    void acquire() {
        int result = ::pthread_mutex_lock(&_mutex);
        assert(result == 0);
    }

    bool tryAcquire() {
        return ::pthread_mutex_trylock(&_mutex) == 0;
    }

    void release() {
        int result = ::pthread_mutex_unlock(&_mutex);
        assert(result == 0);
    }

    friend class ScopedLock<Mutex>;
    friend class Condition<Mutex>;
};


/** @brief A wrapper for a POSIX process shared mutual exclusion lock. */
class LSST_AP_LOCAL SharedMutex : private boost::noncopyable {

public :

    SharedMutex();

    ~SharedMutex() {
        int result = ::pthread_mutex_destroy(&_mutex);
        assert(result == 0);
    }

private :

    ::pthread_mutex_t  _mutex;

    void acquire() {
        int result = ::pthread_mutex_lock(&_mutex);
        assert(result == 0);
    }

    bool tryAcquire() {
        return ::pthread_mutex_trylock(&_mutex) == 0;
    }

    void release() {
        int result = ::pthread_mutex_unlock(&_mutex);
        assert(result == 0);
    }

    friend class ScopedLock<SharedMutex>;
    friend class Condition<SharedMutex>;
};


/** @brief Grants access to a mutex, enforcing the RAII principle. */
template <typename MutexT>
class ScopedLock : private boost::noncopyable {

public :

    ScopedLock() : _mutex(0) {}

    explicit ScopedLock(MutexT & m) : _mutex(&m) {
        m.acquire();
    }

    ~ScopedLock() {
        if (_mutex != 0) {
            _mutex->release();
            _mutex = 0;
        }
    }

    /// Acquires the given Mutex.
    void acquire(MutexT & m) {
        assert(_mutex == 0);
        _mutex = &m;
        m.acquire();
    }

    /**
     * Attempts to acquire the given Mutex, returning immediately if this is not possible.
     *
     * @pre     Any previously acquired Mutex was released
     * @return  @c true if the mutual exclusion lock was acquired, @c false otherwise.
     */
    bool tryAcquire(MutexT & m) {
        assert(_mutex == 0);
        if (m.tryAcquire()) {
            _mutex = &m;
            return true;
        }
        return false;
    }

    /**
     * Releases a previously acquired Mutex.
     *
     * @pre     A Mutex was previously acquired via acquire(Mutex &) or tryAcquire(Mutex &).
     */
    void release() {
        assert(_mutex != 0);
        _mutex->release();
        _mutex = 0;
    }

    /**
     * Returns @c true if and only if a Mutex is currently held; that is, if a Mutex was
     * obtained via acquire(Mutex &) or tryAcquire(Mutex &) but not yet released via release().
     */
    bool isAcquired() const {
        return _mutex != 0;
    }

    // implicit conversion to "bool"
    typedef MutexT * ScopedLock::* UnspecifiedBool;

    operator UnspecifiedBool() const {
        return _mutex == 0 ? 0 : &ScopedLock::_mutex;
    }

    bool operator!() const {
        return _mutex == 0;
    }

private :

    MutexT * _mutex;

    pthread_mutex_t * getPosixMutex() {
        return &(_mutex->_mutex);
    }

    friend class Condition<MutexT>;
};


}} // end of namespace lsst::ap

#endif  // LSST_AP_MUTEX_H
