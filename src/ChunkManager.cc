// -*- lsst-c++ -*-

/**
 * @file
 * @brief   Implementation of shared memory chunk manager for Object instances.
 *
 * @ingroup ap
 */

#include <sys/mman.h>   // for mmap, munmap, shm_open, shm_unlink
#include <time.h>       // for nanosleep
#include <unistd.h>     // for ftruncate
#include <fcntl.h>
#include <errno.h>

#include <iostream>

#include "boost/bind.hpp"
#include "boost/format.hpp"

#include "lsst/pex/exceptions.h"

#include "lsst/ap/Chunk.cc"
#include "lsst/ap/ChunkManagerImpl.cc"
#include "lsst/ap/ChunkManager.h"
#include "lsst/ap/ScopeGuard.h"
#include "lsst/ap/SpatialUtil.h"

namespace ex = lsst::pex::exceptions;


namespace lsst { namespace ap { namespace detail {

// -- Explicit instantiations ----------------

typedef BlockAllocator<SharedMutex, Object> ObjAllocator;

/// @cond
template class BlockAllocator<SharedMutex, Object>;
template class SubManager<SharedMutex, Object>;
template class ChunkManagerImpl<SharedMutex, Object>;
/// @endcond

typedef ChunkManagerImpl<SharedMutex, Object> ObjChunkMgr;

} // end of namespace detail

///@ cond
template class ChunkRef<detail::ObjAllocator, Object>;
///@ endcond


namespace detail {

// -- Shared memory implementation details ----------------

/**
 * @brief   Mutual exclusion lock suitable for initializing shared memory objects.
 *
 * Before standard POSIX shared memory mutexes can be used for inter-process synchronization, shared
 * memory must be allocated and initialized. Therefore, such mutexes cannot be used to protect allocation
 * and initialization of shared memory blocks themselves. BootstrapLock instances get around this by
 * spinning on exclusive creation of a zero-size shared memory object.
 */
struct LSST_AP_LOCAL BootstrapLock {

    char const * const _name;
    int                _fd;

    BootstrapLock(char const * const name);
    ~BootstrapLock();
};


BootstrapLock::BootstrapLock(char const * const name) : _name(name), _fd(-1) {

    assert(name != 0 && name[0] == '/' && "BootstrapLock: illegal name");

    TimeSpec ts;
    ts.tv_nsec = 50000; // start off with 50usec sleep time
    int tries  = 22;
    int fd     = 0;

    while ((fd = ::shm_open(name, O_RDONLY | O_CREAT | O_EXCL, S_IRUSR | S_IRGRP | S_IROTH)) == -1) {
        if (tries == 0) {
            throw LSST_EXCEPT(ex::TimeoutException,
                (boost::format("Unable to obtain shared memory bootstrap lock %1%") % name).str());
        }
        ::nanosleep(&ts, 0);
        ts += ts; // exponential backoff
        --tries;
    }
    _fd = fd;
}


BootstrapLock::~BootstrapLock() {
    if (_fd != -1) {
        ::close(_fd);
        ::shm_unlink(_name);
    }
}


template <typename ManagerT>
LSST_AP_LOCAL ManagerT * getSingleton(char const * const shmObjName, char const * const shmLockName) {

    static Mutex mutex;
    static ManagerT * singleton = 0;

    ScopedLock<Mutex> lck(mutex);
    if (singleton != 0) {
        return singleton;
    }

    // Block interference from other processes
    BootstrapLock blck(shmLockName);

    std::size_t const pageSize  = static_cast<std::size_t>(::getpagesize());
    std::size_t const numBytes  = (ManagerT::size() + pageSize - 1) & ~(pageSize - 1);

    // 1. try to create shared memory object
    int fd = ::shm_open(shmObjName, O_RDWR | O_CREAT | O_EXCL, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);

    if (fd == -1) {

        // failed ...
        if (errno != EEXIST) {
            throw LSST_EXCEPT(ex::RuntimeErrorException,
                (boost::format("shm_open(): failed to create shared memory object %1%, errno: %2%") % shmObjName % errno).str());
        }
        // because the shared memory object already exists -- so try to open existing object instead
        fd = ::shm_open(shmObjName, O_RDWR, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
        if (fd == -1) {
            throw LSST_EXCEPT(ex::RuntimeErrorException,
                (boost::format("shm_open(): failed to open existing shared memory object %1%, errno: %2%") % shmObjName % errno).str());
        }
        ScopeGuard g(boost::bind(::close, fd));

        void * mem = ::mmap(0, numBytes, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
        if (mem == 0) {
            throw LSST_EXCEPT(ex::RuntimeErrorException,
                (boost::format("mmap(): failed to map %1% bytes of shared memory object %2%, errno: %3%") % numBytes % shmObjName % errno).str());
        }
        singleton = static_cast<ManagerT *>(mem);
        g.dismiss();

    } else {

        // shared memory object was created (initially of zero size)
        ScopeGuard g1(boost::bind(::shm_unlink, shmObjName));
        ScopeGuard g2(boost::bind(::close, fd));

        // set size of shared memory object
        if (::ftruncate(fd, numBytes) != 0) {
            throw LSST_EXCEPT(ex::RuntimeErrorException,
                (boost::format("ftruncate(): failed to set size of shared memory object %1% to %2% bytes, errno: %3%") % shmObjName % numBytes % errno).str());
        }

        // map shared memory object
        void * mem = ::mmap(0, numBytes, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
        if (mem == 0) {
            throw LSST_EXCEPT(ex::RuntimeErrorException,
                (boost::format("mmap(): failed to map %1% bytes of shared memory object %2%, errno: %3%") % numBytes % shmObjName % errno).str());
        }
        ScopeGuard g3(boost::bind(::munmap, mem, numBytes));

        new (mem) ManagerT();
        singleton = static_cast<ManagerT *>(mem);
        g3.dismiss();
        g2.dismiss();
        g1.dismiss();
    }

    // Try to lock chunk pages into memory to avoid swapping, but ignore failure to do so
    // - Solaris      : process must be run as root
    // - Linux/Darwin : process must be run as root or RLIMIT_MEMLOCK should be set to a large value
    ::mlock(static_cast<void *>(singleton), numBytes);

    // Note: we rely on the system to munmap and close the file descriptor above at process exit.

    return singleton;
}


// Note: GCC does not support setting template instantiation visibility via
// __attribute__((visibility)) : use pragma instead
#if defined(__GNUC__) && __GNUC__ > 3
#   pragma GCC visibility push(hidden)
#endif
/// @cond
template ObjChunkMgr * getSingleton<ObjChunkMgr>(char const * const, char const * const);
/// @endcond
#if defined(__GNUC__) && __GNUC__ > 3
#   pragma GCC visibility pop
#endif

} // end of namespace detail


// Warning: for portability reasons, these names should be no more than 14 characters long
static char const * const sSharedObjLock = "/ap_shlck";
static char const * const sSharedPrefix  = "/ap_";


// -- SharedObjectChunkManager ----------------

SharedObjectChunkManager::SharedObjectChunkManager(std::string const & name) : _manager(instance(name)) {}


LSST_AP_LOCAL detail::ObjChunkMgr * SharedObjectChunkManager::instance(std::string const & name) {
    std::string actualName(sSharedPrefix);
    actualName += name;
    return detail::getSingleton<detail::ObjChunkMgr>(actualName.c_str(), sSharedObjLock);
}


/**
 * Unlinks the shared memory object underlying all manager instances. The associated memory
 * is not returned to the system until all client processes have relinquished references to it.
 */
void SharedObjectChunkManager::destroyInstance(std::string const & name) {
    std::string actualName(sSharedPrefix);
    actualName += name;
    int res = ::shm_unlink(actualName.c_str());
    // Note: shm_unlink is broken on Mac OSX 10.4 - in violation of the documentation and standard,
    // EINVAL (rather than ENOENT) is returned when trying to shm_unlink a non-existant shared
    // memory object.
    if (res != 0 && errno != ENOENT && errno != EINVAL) {
        throw LSST_EXCEPT(ex::RuntimeErrorException,
            (boost::format("shm_unlink(): failed to unlink shared memory object %1%, errno: %2%") % name % errno).str());
    }
}


/// Returns the size in bytes of the underlying chunk manager and pool of memory blocks.
std::size_t SharedObjectChunkManager::size() { return detail::ObjChunkMgr::size(); }


}} // end of namespace lsst::ap
