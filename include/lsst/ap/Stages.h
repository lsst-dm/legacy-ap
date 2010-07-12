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
 * @brief   C++ pipeline stage implementation methods.
 *
 * @ingroup ap
 */

#ifndef LSST_AP_STAGES_H
#define LSST_AP_STAGES_H

#include <string>
#include <vector>

#include "boost/noncopyable.hpp"
#include "boost/shared_ptr.hpp"

#include "lsst/daf/base/Citizen.h"
#include "lsst/daf/base/PropertySet.h"
#include "lsst/pex/policy/Policy.h"

#include "lsst/afw/detection/DiaSource.h"
#include "lsst/afw/image/Filter.h"
#include "lsst/mops/MovingObjectPrediction.h"

#include "Common.h"
#include "ChunkManager.h"
#include "CircularRegion.h"
#include "Results.h"
#include "SpatialUtil.h"
#include "Time.h"
#include "ZoneTypes.h"


namespace lsst { namespace ap {

#ifndef SWIG
struct DiaSourceChunk {
    typedef lsst::afw::detection::DiaSource Entry;
};
#endif


/** @brief  Container for inter-stage association pipeline state. */
class LSST_AP_API VisitProcessingContext :
    public  lsst::daf::base::Citizen,
    private boost::noncopyable
{
public :
    VisitProcessingContext(
        lsst::pex::policy::Policy::Ptr const policy,
        lsst::daf::base::PropertySet::Ptr const event,
        std::string const & runId,
        int const workerId,
        int const numWorkers
    );

    ~VisitProcessingContext();

    void setDiaSources(boost::shared_ptr<lsst::afw::detection::PersistableDiaSourceVector> diaSources);

#ifndef SWIG

    typedef SharedObjectChunkManager::ObjectChunk ObjectChunk;

    typedef ZoneEntry<ObjectChunk> ObjectEntry;
    typedef ZoneEntry<DiaSourceChunk> DiaSourceEntry;
    typedef ZoneIndex<ObjectEntry> ObjectIndex;
    typedef ZoneIndex<DiaSourceEntry> DiaSourceIndex;

    std::vector<int> const & getChunkIds() const {
        return _chunkIds;
    }
    std::vector<ObjectChunk> const & getChunks() const {
        return _chunks;
    }
    std::vector<int> & getChunkIds() {
        return _chunkIds;
    }
    std::vector<ObjectChunk> & getChunks() {
        return _chunks;
    }

    ObjectIndex & getObjectIndex() {
        return _objectIndex;
    }
    DiaSourceIndex & getDiaSourceIndex() {
        return _diaSourceIndex;
    }

    void buildObjectIndex();

    ZoneStripeChunkDecomposition const & getDecomposition() const {
        return _objectIndex.getDecomposition();
    }
    CircularRegion const & getFov() const {
        return _fov;
    }
    TimeSpec const & getDeadline() const {
        return _deadline;
    }
    lsst::afw::image::Filter getFilter() const {
        return _filter;
    }

#endif

    boost::shared_ptr<lsst::pex::policy::Policy> getPipelinePolicy() {
        return _policy;
    }
    std::string const & getRunId() const {
        return _runId;
    }
    double getMatchRadius() const {
        return _matchRadius;
    }
    double getEllipseScalingFactor() const {
        return _ellipseScalingFactor;
    }
    int getVisitId() const {
        return _visitId;
    }
    int getFilterId() const {
        return _filter.getId();
    }
    int getWorkerId() const {
        return _workerId;
    }
    int getNumWorkers() const {
        return _numWorkers;
    }
    bool debugSharedMemory() const {
        return _debugSharedMemory;
    }

private :

    lsst::pex::policy::Policy::Ptr _policy;
    std::vector<int> _chunkIds;
    std::vector<ObjectChunk> _chunks;
    ObjectIndex _objectIndex;
    DiaSourceIndex _diaSourceIndex;
    std::vector<lsst::afw::detection::DiaSource::Ptr> _diaSources;

    TimeSpec _deadline;
    CircularRegion _fov;
    std::string _runId;
    int _visitId;
    double _matchRadius;
    double _ellipseScalingFactor;
    double _visitTime;
    lsst::afw::image::Filter _filter;
    int _workerId;
    int _numWorkers;
    bool _debugSharedMemory;
};


LSST_AP_API void initialize(std::string const & runId);

LSST_AP_API void registerVisit(VisitProcessingContext & context);

LSST_AP_API void loadSliceObjects(VisitProcessingContext & context);

LSST_AP_API void buildObjectIndex(VisitProcessingContext & context);

LSST_AP_API void matchDiaSources(
    MatchPairVector & matches,
    VisitProcessingContext & context
);

LSST_AP_API void matchMops(
    MatchPairVector & matches,
    IdPairVector & newObjects,
    VisitProcessingContext & context,
    lsst::mops::MovingObjectPredictionVector & predictions
);

LSST_AP_API void storeSliceObjects(VisitProcessingContext & context);

LSST_AP_API void failVisit(VisitProcessingContext & context);

LSST_AP_API bool endVisit(VisitProcessingContext & context, bool const rollback);


}} // end of namespace lsst::ap

#endif // LSST_AP_STAGES_H

