// -*- lsst-c++ -*-

/**
 * @file
 * @brief   C++ pipeline stage implementation methods.
 *
 * @ingroup associate
 */

#ifndef LSST_AP_STAGES_H
#define LSST_AP_STAGES_H

#include <string>
#include <vector>

#include <boost/shared_ptr.hpp>

#include <lsst/daf/base/Citizen.h>
#include <lsst/daf/base/DataProperty.h>
#include <lsst/pex/policy/Policy.h>

#include <lsst/afw/detection/Source.h>
#include <lsst/afw/image/Filter.h>
#include <lsst/mops/MovingObjectPrediction.h>

#include "Common.h"
#include "ChunkManager.h"
#include "CircularRegion.h"
#include "Results.h"
#include "SpatialUtil.h"
#include "Time.h"
#include "ZoneTypes.h"


namespace lsst {
namespace ap {


#ifndef SWIG
struct DiaSourceChunk {
    typedef lsst::afw::detection::Source Entry;
};
#endif


/** @brief  Container for inter-stage association pipeline state. */
class LSST_AP_API VisitProcessingContext : public lsst::daf::base::Citizen {

public :

    VisitProcessingContext(
        lsst::daf::base::DataProperty::PtrType const & event,
        std::string const & runId,
        int         const   workerId,
        int         const   numWorkers
    );

    ~VisitProcessingContext();

    void setDiaSources(lsst::afw::detection::SourceVector & vec);

#ifndef SWIG

    typedef SharedSimpleObjectChunkManager::SimpleObjectChunk SimpleObjectChunk;

    typedef ZoneEntry<SimpleObjectChunk> SimpleObjectEntry;
    typedef ZoneEntry<DiaSourceChunk>    DiaSourceEntry;
    typedef ZoneIndex<SimpleObjectEntry> SimpleObjectIndex;
    typedef ZoneIndex<DiaSourceEntry>    DiaSourceIndex;

    std::vector<int64_t>           const & getChunkIds() const { return _chunkIds; }
    std::vector<SimpleObjectChunk> const & getChunks()   const { return _chunks;   }
    std::vector<int64_t>                 & getChunkIds()       { return _chunkIds; }
    std::vector<SimpleObjectChunk>       & getChunks()         { return _chunks;   }

    SimpleObjectIndex & getObjectIndex()    { return _objectIndex;    }
    DiaSourceIndex    & getDiaSourceIndex() { return _diaSourceIndex; }

    void buildObjectIndex();

    ZoneStripeChunkDecomposition const & getDecomposition() const {
        return _objectIndex.getDecomposition();
    }

    CircularRegion const & getFov()      const { return _fov;      }
    TimeSpec       const & getDeadline() const { return _deadline; }

    lsst::afw::image::Filter getFilter() const { return _filter;      }

#endif

    std::string const & getRunId() const { return _runId; }

    int64_t getVisitId()     const { return _visitId;        }
    double  getMatchRadius() const { return _matchRadius;    }
    int     getFilterId()    const { return _filter.getId(); }
    int     getWorkerId()    const { return _workerId;       }
    int     getNumWorkers()  const { return _numWorkers;     }

private :

    // inhibit copy construction and assignment
    VisitProcessingContext(VisitProcessingContext const &);
    VisitProcessingContext & operator=(VisitProcessingContext const &);

    std::vector<int64_t>           _chunkIds;
    std::vector<SimpleObjectChunk> _chunks;
    SimpleObjectIndex              _objectIndex;
    DiaSourceIndex                 _diaSourceIndex;
    lsst::afw::detection::SourceVector::Ptr _diaSources;

    TimeSpec         _deadline;
    CircularRegion   _fov;
    std::string      _runId;
    int64_t          _visitId;
    double           _matchRadius;
    lsst::afw::image::Filter _filter;
    int              _workerId;
    int              _numWorkers;
};


LSST_AP_API void initialize(lsst::pex::policy::Policy const * policy, std::string const & runId);

LSST_AP_API void registerVisit(VisitProcessingContext & context);

LSST_AP_API void loadSliceObjects(VisitProcessingContext & context);

LSST_AP_API void buildObjectIndex(VisitProcessingContext & context);

LSST_AP_API void matchDiaSources(
    boost::shared_ptr<MatchPairVector> & matches,
    VisitProcessingContext             & context
);

LSST_AP_API void matchMops(
    boost::shared_ptr<MatchPairVector>     & matches,
    boost::shared_ptr<IdPairVector>        & newObjects,
    VisitProcessingContext                 & context,
    lsst::mops::MovingObjectPredictionVector & predictions
);

LSST_AP_API void storeSliceObjects(VisitProcessingContext & context);

LSST_AP_API void failVisit(VisitProcessingContext & context);

LSST_AP_API bool endVisit(VisitProcessingContext & context, bool const rollback);


}} // end of namespace lsst::ap

#endif // LSST_AP_STAGES_H

