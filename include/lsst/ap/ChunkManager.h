// -*- lsst-c++ -*-
//
//##====----------------                                ----------------====##/
//
//! \file   ChunkManager.h
//
//##====----------------                                ----------------====##/

#ifndef LSST_AP_CHUNK_MANAGER_H
#define LSST_AP_CHUNK_MANAGER_H

#include <iosfwd>

#include "ChunkManagerImpl.h"
#include "Object.h"


namespace lsst {
namespace ap {


/*! \brief  A manager for SimpleObject chunks that exist in shared memory. */
class LSST_AP_API SharedSimpleObjectChunkManager {

private :

    typedef detail::ChunkManagerSingleImpl<SharedMutex, SimpleObject> ManagerType;

    ManagerType * _manager;

    static LSST_AP_LOCAL ManagerType * instance(std::string const & name);

public :

    typedef ManagerType::ChunkType SimpleObjectChunkType;

    SharedSimpleObjectChunkManager(std::string const & name);

    bool isVisitInFlight(int64_t const visitId) { return _manager->isVisitInFlight(visitId); }
    void registerVisit  (int64_t const visitId) { _manager->registerVisit(visitId);          }
    void failVisit      (int64_t const visitId) { _manager->failVisit(visitId);              }

    void startVisit(
        std::vector<SimpleObjectChunkType>       & toRead,
        std::vector<SimpleObjectChunkType>       & toWaitFor,
        int64_t                            const   visitId,
        std::vector<int64_t>               const & chunkIds
    ) {
        _manager->startVisit(toRead, toWaitFor, visitId, chunkIds);
    }

    void waitForOwnership(
        std::vector<SimpleObjectChunkType> & toRead,
        std::vector<SimpleObjectChunkType> & toWaitFor,
        int64_t                      const   visitId,
        TimeSpec                     const & deadline
    ) {
        _manager->waitForOwnership(toRead, toWaitFor, visitId, deadline);
    }

    void getChunks(
        std::vector<SimpleObjectChunkType> & chunks,
        std::vector<int64_t>         const & chunkIds
    ) {
        _manager->getChunks(chunks, chunkIds);
    }

    bool endVisit(int64_t const visitId, bool const rollback) {
        return _manager->endVisit(visitId, rollback);
    }

    void printVisits(std::ostream & os) const { _manager->printVisits(os); }
    void printChunks(std::ostream & os) const { _manager->printChunks(os); }
    void printVisit(int64_t const visitId, std::ostream & os) const { _manager->printVisit(visitId, os); }
    void printChunk(int64_t const chunkId, std::ostream & os) const { _manager->printChunk(chunkId, os); }

    static void destroyInstance(std::string const & name);

    static size_t size();
};


}} // end of namespace lsst::ap

#endif // LSST_AP_CHUNK_MANAGER_H
