// -*- lsst-c++ -*-

/**
 * @file
 * @brief   Class for managing chunks of Object instances in shared memory.
 *
 * @ingroup associate
 */

#ifndef LSST_AP_CHUNK_MANAGER_H
#define LSST_AP_CHUNK_MANAGER_H

#include <iosfwd>

#include "ChunkManagerImpl.h"
#include "Object.h"


namespace lsst { namespace ap {

/** @brief  A manager for Object chunks that exist in shared memory. */
class LSST_AP_API SharedObjectChunkManager {

private :

    typedef detail::ChunkManagerImpl<SharedMutex, Object> Manager;

    Manager * _manager;

    static LSST_AP_LOCAL Manager * instance(std::string const & name);

public :

    typedef Manager::Chunk ObjectChunk;

    SharedObjectChunkManager(std::string const & name);

    bool isVisitInFlight(int const visitId) {
        return _manager->isVisitInFlight(visitId);
    }
    void registerVisit(int const visitId) {
        _manager->registerVisit(visitId);
    }
    void failVisit(int const visitId) {
        _manager->failVisit(visitId);
    }

    void startVisit(
        std::vector<ObjectChunk> & toRead,
        std::vector<ObjectChunk> & toWaitFor,
        int const visitId,
        std::vector<int> const & chunkIds
    ) {
        _manager->startVisit(toRead, toWaitFor, visitId, chunkIds);
    }

    void waitForOwnership(
        std::vector<ObjectChunk> & toRead,
        std::vector<ObjectChunk> & toWaitFor,
        int const visitId,
        TimeSpec const & deadline
    ) {
        _manager->waitForOwnership(toRead, toWaitFor, visitId, deadline);
    }

    void getChunks(
        std::vector<ObjectChunk> & chunks,
        std::vector<int> const & chunkIds
    ) {
        _manager->getChunks(chunks, chunkIds);
    }

    bool endVisit(int const visitId, bool const rollback) {
        return _manager->endVisit(visitId, rollback);
    }

    void printVisits(std::ostream & os) const {
        _manager->printVisits(os);
    }
    void printChunks(std::ostream & os) const {
        _manager->printChunks(os);
    }
    void printVisit(int const visitId, std::ostream & os) const {
        _manager->printVisit(visitId, os);
    }
    void printChunk(int const chunkId, std::ostream & os) const {
        _manager->printChunk(chunkId, os);
    }

    static void destroyInstance(std::string const & name);

    static std::size_t size();
};


}} // end of namespace lsst::ap

#endif // LSST_AP_CHUNK_MANAGER_H
