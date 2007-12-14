// -*- lsst-c++ -*-
//
//##====----------------                                ----------------====##/
//
//! \file   Stages.cc
//! \brief  Implementation of association pipeline stages.
//!
//! \todo   How is the time of observation for a visit specified? And how is
//!         the association pipeline deadline derived from this? Possibly post-DC2.
//!
//! \todo   [Post DC2] How are variability probabilities set for new objects? Currently,
//!         the probability is set to 100% for any new Object.
//!
//! \todo   [Post DC2] StoreStage can issue multiple database requests in parallel
//!
//! \todo   The end of the pipeline needs to send out a 'triggerAlertGeneration' event.
//!         Possibly post-DC2, since there is no alert generation pipeline.
//!
//! \todo   [Post DC2] Don't create new objects from difference sources which convincingly
//!         match a moving object prediction.
//!
//! \todo   [Post DC2] Better exception handling on the Python side
//!
//! \todo   [Post DC2] Add a version number to chunk delta file names. Use a database
//!         table to track the version numbers of valid delta files in transactional
//!         fashion. Note that the transaction which updates these version numbers
//!         should include the insert to the object table as well as the append to the
//!         historical DIASource table. Use MySQL API directly for this?
//
//##====----------------                                ----------------====##/

#if LSST_HAVE_OPEN_MP
#   include <omp.h>
#endif

#include <algorithm>
#include <memory>

#include <boost/any.hpp>
#include <boost/format.hpp>
#include <boost/scoped_array.hpp>

#include <lsst/mwi/logging/Log.h>

#include <lsst/fw/DiaSource.h>
#include <lsst/fw/Filter.h>

#include <lsst/ap/ChunkManager.h>
#include <lsst/ap/ChunkToNameMappings.h>
#include <lsst/ap/Match.h>
#include <lsst/ap/Stages.h>
#include <lsst/ap/Time.h>
#include <lsst/ap/Utils.h>


namespace lsst {
namespace ap {

using lsst::mwi::data::DataProperty;
using lsst::mwi::logging::Log;
using lsst::mwi::logging::Rec;


// -- Constants ----------------

static uint32_t const HAS_MATCH                    = 1;
static uint32_t const HAS_KNOWN_VARIABLE_MATCH     = 2;
static int8_t   const DEF_VP_THRESH                = 90;
static int32_t  const DEF_ZONES_PER_DEGREE         = 180;   // 20 arc-second zone height
static int32_t  const DEF_ZONES_PER_STRIPE         = 63;    // 0.35 degree stripes
static int32_t  const DEF_MAX_ENTRIES_PER_ZONE_EST = 4096;
static double   const DEF_SMAA_THRESH              = 300;   // arc-seconds == 5 arc-minutes
static double   const DEF_MATCH_RADIUS             = 0.05;  // arc-seconds

static char const * const DEF_OBJ_PATTERN       = "/tmp/%1%/objchunk/stripe_%2%/objref_chunk%3%";
static char const * const DEF_OBJ_DELTA_PATTERN = "/tmp/%1%/objdelta/stripe_%2%/objdelta_chunk%3%";
static char const * const DEF_DB_LOCATION       = "mysql://lsst10.ncsa.uiuc.edu:3306/test";

// -- Adjustable (by policy) configuration parameters ----------------

static int32_t sZonesPerDegree       = DEF_ZONES_PER_DEGREE;
static int32_t sZonesPerStripe       = DEF_ZONES_PER_STRIPE;
static int32_t sMaxEntriesPerZoneEst = DEF_MAX_ENTRIES_PER_ZONE_EST;


/*! The radius of a visit FOV. */
static double sFovRadius = FOV_RADIUS;

/*!
    Any error ellipses for predicted moving object positions with a semi-major axis length greater
    than this threshold (in arc-seconds) are ignored by the association pipeline.
 */
static double sSemiMajorAxisThreshold = DEF_SMAA_THRESH;

/*!
    Unless overridden by the detection pipeline/policy, this is the default match radius
    (in arc-seconds) for matching difference sources to objects.
 */
static double sMatchRadius = DEF_MATCH_RADIUS;

/*!
    The pattern (\c boost::format compatible) of simple object chunk file names. The pattern
    is supplied with 2 parameters - an integer stripe id and an integer chunk id - which may be used to
    construct a unique file name.
 */
static std::string sObjFilePattern(DEF_OBJ_PATTERN);

/*!
    The pattern (\c boost::format compatible) of simple object delta chunk file names. The pattern is
    supplied with 3 parameters - an integer stripe id, an integer chunk id, and an integer version id -
    which may be used to construct a file name.
 */
static std::string sObjDeltaFilePattern(DEF_OBJ_DELTA_PATTERN);

/*!
    Per-filter variable-object probability threshold above which an object is considered to be a
    "known variable". DiaSources matching "known variables" are not considered when matching against
    predicted positions of moving objects.
 */
static int8_t sVarProbThreshold[lsst::fw::Filter::NUM_FILTERS] = {
    DEF_VP_THRESH,
    DEF_VP_THRESH,
    DEF_VP_THRESH,
    DEF_VP_THRESH,
    DEF_VP_THRESH,
    DEF_VP_THRESH
};

/*! The location of the \b Filter table for mapping filter identifiers to names and vice versa. */
static std::string sFilterTableLocation(DEF_DB_LOCATION);


namespace detail {

// -- Typedefs and templates for chunks and spatial indexes ----------------

typedef SharedSimpleObjectChunkManager::SimpleObjectChunkType SimpleObjectChunk;

typedef std::vector<SimpleObjectChunk>  SimpleObjectChunkVector;

typedef ZoneEntry<SimpleObjectChunk>    SimpleObjectEntry;
typedef ZoneEntry<DiaSourceChunk>       DiaSourceEntry;

typedef Ellipse<lsst::fw::MovingObjectPrediction> MovingObjectEllipse;


#if defined(__GNUC__) && __GNUC__ > 3
#   pragma GCC visibility push(hidden)
#endif
//! \cond
template class ZoneEntry<SimpleObjectChunk>;
template class ZoneEntry<DiaSourceChunk>;

template class Zone<SimpleObjectEntry>;
template class Zone<DiaSourceEntry>;

template class ZoneIndex<SimpleObjectEntry>;
template class ZoneIndex<DiaSourceEntry>;

template class Ellipse<lsst::fw::MovingObjectPrediction>;
template class EllipseList<lsst::fw::MovingObjectPrediction>;
//! \endcond
#if defined(__GNUC__) && __GNUC__ > 3
#   pragma GCC visibility pop
#endif


// -- Match processors ----------------

/*! \brief  Processor for lists of objects matching a difference source */
template <typename ZoneEntryType>
class ObjectMatchProcessor {

public :

    typedef MatchWithDistance<ZoneEntryType>          MatchType;
    typedef typename std::vector<MatchType>::iterator MatchIteratorType;

    MatchPairVector        & _matches;
    lsst::fw::Filter const   _filter;
    int8_t           const   _threshold;

    ObjectMatchProcessor(
        MatchPairVector        & matches,
        lsst::fw::Filter const   filter
    ) :
        _matches(matches),
        _filter(filter),
        _threshold(sVarProbThreshold[filter])
    {}

    void operator()(DiaSourceEntry & ds, MatchIteratorType begin, MatchIteratorType end) {
        uint32_t flags = HAS_MATCH;
        do {
            ZoneEntryType * const obj = begin->_match;
            ++begin;
            // record match results (to be persisted later)
            _matches.push_back(MatchPair(ds._data->getId(), obj->_data->getId(), degrees(obj->_distance)));
            if (obj->_data->getVarProb(_filter) >= _threshold) {
                // flag ds as matching a known variable
                flags |= HAS_KNOWN_VARIABLE_MATCH;
            }
        } while (begin != end);
        ds._flags = flags;
    }
};


/*! \brief  Processor for matches between moving object predictions and difference sources. */
struct LSST_AP_LOCAL MovingObjectPredictionMatchProcessor {

    typedef DiaSourceEntry *                 MatchType;
    typedef std::vector<MatchType>::iterator MatchIteratorType;

    MatchPairVector & _results;

    explicit MovingObjectPredictionMatchProcessor(MatchPairVector & results) : _results(results) {}

    void operator()(MovingObjectEllipse & ell, DiaSourceEntry & ds) {
        double dx   = ell._cosRa * ell._cosDec - ds._x;
        double dy   = ell._sinRa * ell._cosDec - ds._y;
        double dz   = ell._sinDec              - ds._z;
        double dist = 2.0*std::asin(0.5*std::sqrt(dx*dx + dy*dy + dz*dz));
        // record match results (to be persisted later)
        _results.push_back(MatchPair(ell._data->getId(), ds._data->getId(), degrees(dist)));
    }
};


// -- Zone index filters and functors ----------------

/*! \brief  Filter which discards difference sources matching known variable objects. */
struct LSST_AP_LOCAL DiscardKnownVariableFilter {
    bool operator()(DiaSourceEntry const & ds) {
        return (ds._flags & HAS_KNOWN_VARIABLE_MATCH) == 0;
    }
};


/*! \brief  Filter which discards predicted moving objects with large position error ellipses. */
struct LSST_AP_LOCAL DiscardLargeEllipseFilter {
    bool operator()(MovingObjectEllipse const & ell) {
        return ell._data->getSemiMajorAxisLength() < sSemiMajorAxisThreshold;
    }
};


/*! \brief  Records ids of difference sources with no matches. */
struct LSST_AP_LOCAL NewObjectCreator {

    typedef std::map<int64_t, SimpleObjectChunk> ChunkMapType;
    typedef ChunkMapType::value_type             ChunkMapValueType;
    typedef ChunkMapType::iterator               ChunkMapIteratorType;

    IdPairVector                       & _results;
    ZoneStripeChunkDecomposition const & _zsc;
    std::map<int64_t, SimpleObjectChunk> _chunks;
    lsst::fw::Filter const               _filter;
    int64_t          const               _idNamespace;

    explicit NewObjectCreator(
        IdPairVector                       & results,
        SimpleObjectChunkVector            & chunks,
        ZoneStripeChunkDecomposition const & zsc,
        lsst::fw::Filter             const   filter
    ) :
        _results(results),
        _zsc(zsc),
        _chunks(),
        _filter(filter),
        _idNamespace(static_cast<int64_t>(filter + 1) << 56)
    {
        // build a map of ids to chunks
        SimpleObjectChunkVector::iterator const end = chunks.end();
        for (SimpleObjectChunkVector::iterator i = chunks.begin(); i != end; ++i) {
            _chunks.insert(ChunkMapValueType(i->getId(), *i));
        }
    }

    void operator()(DiaSourceEntry const & entry) {
        static int64_t const idLimit = INT64_C(1) << 56;

        if ((entry._flags & (HAS_MATCH | HAS_KNOWN_VARIABLE_MATCH)) == 0) {
            // difference source had no matches - record it as the source of a new object
            int64_t id = entry._data->getId();
            if (id >= idLimit) {
                LSST_AP_THROW(OutOfRange, "DiaSource id doesn't fit in 56 bits");
            }
            // generate a new simplified object (id, position, variability probabilities only)
            // and assign it to the appropriate chunk
            SimpleObject obj;
            
            obj._objectId           = id | _idNamespace;
            obj._ra                 = entry._data->getRa();
            obj._decl               = entry._data->getDec();
            obj._varProb[lsst::fw::Filter::U] = 0;
            obj._varProb[lsst::fw::Filter::G] = 0;
            obj._varProb[lsst::fw::Filter::R] = 0;
            obj._varProb[lsst::fw::Filter::I] = 0;
            obj._varProb[lsst::fw::Filter::Z] = 0;
            obj._varProb[lsst::fw::Filter::Y] = 0;
            obj._varProb[_filter]   = 100;

            _results.push_back(IdPair(id, obj._objectId));

            // find the chunk the new object belongs to and insert the new object into it
            ChunkMapIteratorType c = _chunks.find(_zsc.radecToChunk(obj._ra, obj._decl));
            if (c == _chunks.end()) {
                LSST_AP_THROW(Runtime, "new object falls outside of object chunks covering the FOV");
            }
            c->second.insert(obj);
        }
    }
};


// -- Index creation ----------------

template <typename EntryType>
void buildZoneIndex(
    ZoneIndex<EntryType>                             & index,
    std::vector<typename EntryType::ChunkType> const & chunks
) {
    typedef typename EntryType::DataType             DataType;
    typedef typename EntryType::ChunkType            ChunkType;
    typedef          std::vector<ChunkType>          ChunkVectorType;
    typedef typename ChunkVectorType::const_iterator ChunkIteratorType;
    typedef typename ChunkVectorType::size_type      SizeType;

    index.clear();
    if (chunks.empty()) {
        return;
    }

    Stopwatch watch(true);

    // determine stripe bounds for the input chunks
    ZoneStripeChunkDecomposition const & zsc = index.getDecomposition();
    int32_t minStripe = 0x7FFFFFFF;
    int32_t maxStripe = -1 - minStripe;
    ChunkIteratorType const end(chunks.end());
    for (ChunkIteratorType c(chunks.begin()); c != end; ++c) {
        int32_t const stripeId = ZoneStripeChunkDecomposition::chunkToStripe(c->getId());
        if (stripeId > maxStripe) {
            maxStripe = stripeId;
        }
        if (stripeId < minStripe) {
            minStripe = stripeId;
        }
    }
    assert(maxStripe >= minStripe && "invalid stripe bounds for chunk list");

    // Partition input chunks into stripes
    int32_t const numStripes = maxStripe - minStripe + 1;
    boost::scoped_array<ChunkVectorType> stripes(new ChunkVectorType[numStripes]);
    for (ChunkIteratorType c(chunks.begin()); c != end; ++c) {
        int32_t const stripeId = ZoneStripeChunkDecomposition::chunkToStripe(c->getId());
        assert(stripeId >= minStripe && stripeId <= maxStripe && "stripe id out of bounds");
        stripes[stripeId - minStripe].push_back(*c);
    }

    double minDec = zsc.getStripeDecMin(minStripe) - 0.001;
    double maxDec = zsc.getStripeDecMax(maxStripe) + 0.001;
    index.setDecBounds(std::max(minDec, -90.0), std::min(maxDec, 90.0));

    // Loop over stripes: since each stripe maps to a distinct set of adjacent zones,
    // multiple stripes can be added to the index in parallel without any locking.
    volatile bool failed = false;
#if LSST_HAVE_OPEN_MP
#   pragma omp parallel for shared(numStripes, stripes, failed) \
               schedule(dynamic,1)
#endif
    for (int32_t s = 0; s < numStripes; ++s) {
        try {
            ChunkVectorType & vec       = stripes[s];
            SizeType const    numChunks = vec.size();

            // Loop over chunks in stripe
            for (SizeType c = 0; c < numChunks; ++c) {

                ChunkType * const ch         = &vec[c];
                uint32_t    const numBlocks  = ch->blocks();
                uint32_t          i          = 0;

                // loop over blocks in chunk
                for (uint32_t b = 0; b < numBlocks; ++b) {
                    uint32_t const numEntries = ch->entries(b);
                    DataType                 * const block = ch->getBlock(b);
                    ChunkEntryFlagType const * const flags = ch->getFlagBlock(b);

                    // loop over entries in block
                    for (uint32_t e = 0; e < numEntries; ++e, ++i) {
                        if ((flags[e] & ChunkType::DELETED) == 0) {
                            index.insert(&block[e], ch, i);
                        }
                    }
                }
            }
        } catch(...) {
            // Don't throw in threaded section of code!
            failed = true;
        }
    } // end of omp parallel for

    // If any worker failed, throw an exception
    if (failed) {
        index.clear();
        LSST_AP_THROW(Runtime, "Failed to build zone index");
    }
    watch.stop();
    Log log(Log::getDefaultLog(), "associate");
    size_t numElements = index.size();
    Rec(log, Log::INFO) << "inserted elements into zone index" <<
        DataProperty("numElements", static_cast<long>(numElements)) <<
        DataProperty("time", watch.seconds()) << Rec::endr;

    // zone structure is filled, sort individual zones (on right ascension)
    watch.start();
    index.sort();
    watch.stop();
    Rec(log, Log::INFO) << "sorted zone index" <<
        DataProperty("numElements", static_cast<long>(numElements)) <<
        DataProperty("time", watch.seconds()) << Rec::endr;
}


// -- Template instantiations ----------------

#if defined(__GNUC__) && __GNUC__ > 3
#   pragma GCC visibility push(hidden)
#endif
//! \cond
template class ObjectMatchProcessor<SimpleObjectEntry>;

template size_t distanceMatch<
    DiaSourceEntry,
    SimpleObjectEntry,
    PassthroughFilter<DiaSourceEntry>,
    PassthroughFilter<SimpleObjectEntry>,
    ObjectMatchProcessor<SimpleObjectEntry>
>(
    ZoneIndex<DiaSourceEntry> &,
    ZoneIndex<SimpleObjectEntry> &,
    double const,
    PassthroughFilter<DiaSourceEntry> &,
    PassthroughFilter<SimpleObjectEntry> &,
    ObjectMatchProcessor<SimpleObjectEntry> &
);

template size_t ellipseMatch<
    lsst::fw::MovingObjectPrediction,
    DiaSourceEntry,
    DiscardLargeEllipseFilter,
    DiscardKnownVariableFilter,
    MovingObjectPredictionMatchProcessor
>(
    EllipseList<lsst::fw::MovingObjectPrediction> &,
    ZoneIndex<DiaSourceEntry> &,
    DiscardLargeEllipseFilter &,
    DiscardKnownVariableFilter &,
    MovingObjectPredictionMatchProcessor &
);
//! \endcond
#if defined(__GNUC__) && __GNUC__ > 3
#   pragma GCC visibility pop
#endif


} // end of namespace detail


// -- VisitProcessingContext ----------------

VisitProcessingContext::VisitProcessingContext(
    lsst::mwi::data::DataProperty::PtrType const & event,
    std::string                            const & runId,
    int const workerId,
    int const numWorkers
) :
    _chunkIds(),
    _chunks(),
    _objectIndex(sZonesPerDegree, sZonesPerStripe, sMaxEntriesPerZoneEst),
    _diaSourceIndex(sZonesPerDegree, sZonesPerStripe, sMaxEntriesPerZoneEst),
    _deadline(),
    _fov(),
    _runId(runId),
    _visitId(-1),
    _matchRadius(sMatchRadius),
    _filter(),
    _workerId(workerId),
    _numWorkers(numWorkers)
{
    lsst::mwi::data::DataProperty::PtrType dp = extractRequired(event, "FOVRA");
    double ra = boost::any_cast<double>(dp->getValue());
    dp = extractRequired(event, "FOVDec");
    double dec = boost::any_cast<double>(dp->getValue());
    _fov = CircularRegion(ra, dec, sFovRadius);
    dp = extractRequired(event, "visitId");
    _visitId = anyToInteger<int64_t>(dp->getValue());
    dp = event->findUnique("matchRadius");
    if (dp) {
        _matchRadius = boost::any_cast<double>(dp->getValue());
    }
    // DC2: set association pipeline deadline to 10 minutes
    // after creation of a visit processing context.
    _deadline.systemTime();
    _deadline.tv_sec += 600;
    dp = extractRequired(event, "filterName");
    std::string filterName = boost::any_cast<std::string>(dp->getValue());
    lsst::mwi::persistence::LogicalLocation location(sFilterTableLocation);
    _filter = lsst::fw::Filter(location, filterName);
}


VisitProcessingContext::~VisitProcessingContext() {}


void VisitProcessingContext::setDiaSources(lsst::fw::DiaSourceVector & vec) {
    _diaSourceIndex.clear();

    Stopwatch watch(true);
    lsst::fw::DiaSourceVector::size_type const sz = vec.size();
    double minDec =  90.0;
    double maxDec = -90.0;
    for (lsst::fw::DiaSourceVector::size_type i = 0; i < sz; ++i) {
        double dec = vec[i].getDec();
        if (dec < minDec) {
            minDec = dec;
        }
        if (dec > maxDec) {
            maxDec = dec;
        }
    }
    assert(maxDec >= minDec && "invalid dec bounds for DiaSource list");
    _diaSourceIndex.setDecBounds(minDec, maxDec);
    watch.stop();
    Log log(Log::getDefaultLog(), "associate");
    Rec(log, Log::INFO) <<  "set dec bounds for difference source index" <<
        DataProperty("numElements", static_cast<long>(sz)) <<
        DataProperty("time", watch.seconds()) << Rec::endr;
    watch.start();
    try {
        for (lsst::fw::DiaSourceVector::size_type i = 0; i < sz; ++i) {
            _diaSourceIndex.insert(&vec[i], 0, 0);
        }
    } catch (...) {
        _diaSourceIndex.clear();
        throw;
    }
    watch.stop();
    Rec(log, Log::INFO) << "inserted difference sources into zone index" <<
        DataProperty("numElements", static_cast<long>(sz)) <<
        DataProperty("time", watch.seconds()) << Rec::endr;
    watch.start();
    _diaSourceIndex.sort();
    watch.stop();
    Rec(log, Log::INFO) << "sorted difference source zone index" <<
        DataProperty("numElements", static_cast<long>(sz)) <<
        DataProperty("time", watch.seconds()) << Rec::endr;
}


void VisitProcessingContext::buildObjectIndex() {
    detail::buildZoneIndex(_objectIndex, _chunks);
}


// -- Load stage ----------------

/*!
    Sets up all fundamental visit processing parameters using a policy and ensure
    that a reference to the shared memory object used for chunk storage exists.
 */
LSST_AP_API void initialize(lsst::mwi::policy::Policy const * policy, std::string const & runId) {

    volatile SharedSimpleObjectChunkManager manager(runId);

    sZonesPerDegree         = DEF_ZONES_PER_DEGREE;
    sZonesPerStripe         = DEF_ZONES_PER_STRIPE;
    sMaxEntriesPerZoneEst   = DEF_MAX_ENTRIES_PER_ZONE_EST;
    sSemiMajorAxisThreshold = DEF_SMAA_THRESH;
    sObjFilePattern         = DEF_OBJ_PATTERN;
    sObjDeltaFilePattern    = DEF_OBJ_DELTA_PATTERN;

    if (policy) {
        sFovRadius              = policy->getDouble("fovRadius",              FOV_RADIUS);
        sSemiMajorAxisThreshold = policy->getDouble("semiMajorAxisThreshold", DEF_SMAA_THRESH);
        sMatchRadius            = policy->getDouble("matchRadius",            DEF_MATCH_RADIUS);

        sZonesPerDegree       = policy->getInt("zonesPerDegree",            DEF_ZONES_PER_DEGREE);
        sZonesPerStripe       = policy->getInt("zonesPerStripe",            DEF_ZONES_PER_STRIPE);
        sMaxEntriesPerZoneEst = policy->getInt("maxEntriesPerZoneEstimate", DEF_MAX_ENTRIES_PER_ZONE_EST);

        sObjFilePattern      = policy->getString("objectChunkFileNamePattern",      DEF_OBJ_PATTERN);
        sObjDeltaFilePattern = policy->getString("objectDeltaChunkFileNamePattern", DEF_OBJ_DELTA_PATTERN);
        sFilterTableLocation = policy->getString("filterTableLocation",             DEF_DB_LOCATION);

        sVarProbThreshold[lsst::fw::Filter::U] = policy->getInt("uVarProbThreshold", DEF_VP_THRESH);
        sVarProbThreshold[lsst::fw::Filter::G] = policy->getInt("gVarProbThreshold", DEF_VP_THRESH);
        sVarProbThreshold[lsst::fw::Filter::R] = policy->getInt("rVarProbThreshold", DEF_VP_THRESH);
        sVarProbThreshold[lsst::fw::Filter::I] = policy->getInt("iVarProbThreshold", DEF_VP_THRESH);
        sVarProbThreshold[lsst::fw::Filter::Z] = policy->getInt("zVarProbThreshold", DEF_VP_THRESH);
        sVarProbThreshold[lsst::fw::Filter::Y] = policy->getInt("yVarProbThreshold", DEF_VP_THRESH);
    }
}


/*!
    Computes ids for all object chunks covering the visit FOV and
    registers the visit with the shared memory chunk manager.
 */
LSST_AP_API void registerVisit(VisitProcessingContext & context) {
    context.getChunkIds().clear();
    computeChunkIds(context.getChunkIds(), context.getFov(), context.getDecomposition(), 0, 1);
    SharedSimpleObjectChunkManager manager(context.getRunId());
    manager.registerVisit(context.getVisitId());
}


/*!
    Ensures that object data for the chunks assigned to the calling slice has been read in or is
    owned by the given visit.
 */
LSST_AP_API void loadSliceObjects(VisitProcessingContext & context) {

    typedef VisitProcessingContext::SimpleObjectChunk ChunkType;
    typedef std::vector<ChunkType>                    ChunkVectorType;
    typedef std::vector<ChunkType>::iterator          ChunkIteratorType;

    SharedSimpleObjectChunkManager manager(context.getRunId());
    Log log(Log::getDefaultLog(), "associate");

    try {

        // From FOV, get chunk ids to be handled by this worker
        Stopwatch watch(true);
        computeChunkIds(
            context.getChunkIds(),
            context.getFov(),
            context.getDecomposition(),
            context.getWorkerId(),
            context.getNumWorkers()
        );
        watch.stop();
        Rec(log, Log::INFO) << "computed chunk ids in FOV for worker" <<
            DataProperty("numChunks", static_cast<long>(context.getChunkIds().size())) <<
            DataProperty("time", watch.seconds()) << Rec::endr;

        // Register interest in or create chunks via the chunk manager
        std::vector<ChunkType> toRead;
        std::vector<ChunkType> toWaitFor;
        watch.start();
        manager.startVisit(toRead, toWaitFor, context.getVisitId(), context.getChunkIds());
        watch.stop();
        Rec(log, Log::INFO) << "started processing visit" <<
            DataProperty("time", watch.seconds()) << Rec::endr;

        // record pointers to all chunks being handled by the slice
        ChunkVectorType & chunks = context.getChunks();
        chunks.clear();
        chunks.insert(chunks.end(), toRead.begin(),    toRead.end());
        chunks.insert(chunks.end(), toWaitFor.begin(), toWaitFor.end());

        // Read data files
        watch.start();
        ChunkToFileNameMapping     refNames(sObjFilePattern);
        ChunkToFileNameMapping     deltaNames(sObjDeltaFilePattern);
        ChunkIteratorType          end(toRead.end());
        ChunkVectorType::size_type numToRead(toRead.size());
        for (ChunkIteratorType i(toRead.begin()); i != end; ++i) {
            ChunkType & c = *i;
            c.read(
                refNames.getName(context.getRunId(), context.getDecomposition(), c.getId()),
                false
            );
            c.readDelta(
                deltaNames.getName(context.getRunId(), context.getDecomposition(), c.getId()),
                false
            );
            c.setUsable();
        }
        watch.stop();
        Rec(log, Log::INFO) << "read chunk files" <<
            DataProperty("numChunks", static_cast<long>(numToRead)) <<
            DataProperty("time", watch.seconds()) << Rec::endr;

        toRead.clear();
        ChunkVectorType::size_type numToWaitFor = toWaitFor.size();
        if (numToWaitFor > 0) {
            watch.start();
            // Wait for chunks that are owned by another visit
            manager.waitForOwnership(toRead, toWaitFor, context.getVisitId(), context.getDeadline());
            watch.stop();
            Rec(log, Log::INFO) << "acquired ownership of pre-existing chunks" <<
                DataProperty("numChunks", static_cast<long>(numToWaitFor)) <<
                DataProperty("time", watch.seconds()) << Rec::endr;

            // Read in chunks that were not successfully read by the previous owner
            watch.start();
            numToRead = toRead.size();
            end = toRead.end();
            for (ChunkIteratorType i(toRead.begin()); i != end; ++i) {
                ChunkType & c = *i;
                c.read(
                    refNames.getName(context.getRunId(), context.getDecomposition(), c.getId()),
                    false
                );
                c.readDelta(
                    deltaNames.getName(context.getRunId(), context.getDecomposition(), c.getId()),
                    false
                );
                c.setUsable();
            }
            watch.stop();
            Rec(log, Log::INFO) << "read straggling chunks" <<
                DataProperty("numChunks", static_cast<long>(numToRead)) <<
                DataProperty("time", watch.seconds()) << Rec::endr;
        }

    } catch (lsst::mwi::exceptions::ExceptionStack & ex) {
        Rec(log, Log::FATAL) << ex.what() << *(ex.getStack()) << Rec::endr;
        manager.failVisit(context.getVisitId());
    } catch (std::exception & ex) {
        log.log(Log::FATAL, ex.what());
        manager.failVisit(context.getVisitId());
    } catch (...) {
        log.log(Log::FATAL, "caught unknown exception");
        manager.failVisit(context.getVisitId());
    }
}


/*!
    Checks to see if all parallel workers succeeded in loading their chunks and, if so, builds
    a zone index for the objects just loaded.

    \param[in, out] context     State involved in processing a single visit.
 */
LSST_AP_API void buildObjectIndex(VisitProcessingContext & context) {
    // if the shared memory object used for chunk storage hasn't yet been unlinked, do so now
    SharedSimpleObjectChunkManager::destroyInstance(context.getRunId());
    SharedSimpleObjectChunkManager manager(context.getRunId());
    if (manager.isVisitInFlight(context.getVisitId())) {
        // Build zone index on objects
        manager.getChunks(context.getChunks(), context.getChunkIds());
        context.buildObjectIndex();
    } else {
        // One or more workers failed in the load phase - rollback the visit
        manager.endVisit(context.getVisitId(), true);
        LSST_AP_THROW(Runtime, "Association pipeline failed to read Object data for FOV");
    }
}


// -- MatchDiaSource stage ----------------

/*!
    Matches difference sources for a visit (obtained from the detection pipeline) against the objects
    in the visit FOV.

    \param[out]     matches     Set to a list of difference source to object match pairs.
    \param[in, out] context     State involved in processing a single visit.
 */
LSST_AP_API void matchDiaSources(
    boost::shared_ptr<MatchPairVector> & matches,
    VisitProcessingContext             & context
) {
    SharedSimpleObjectChunkManager manager(context.getRunId());

    try {

        matches.reset(new MatchPairVector());
        matches->reserve(65536);

        detail::ObjectMatchProcessor<detail::SimpleObjectEntry> mlp(*matches, context.getFilter());
        PassthroughFilter<detail::DiaSourceEntry>    pdf;
        PassthroughFilter<detail::SimpleObjectEntry> pof;

        Stopwatch watch(true);
        size_t nm = distanceMatch<
            detail::DiaSourceEntry,
            detail::SimpleObjectEntry,
            PassthroughFilter<detail::DiaSourceEntry>,
            PassthroughFilter<detail::SimpleObjectEntry>,
            detail::ObjectMatchProcessor<detail::SimpleObjectEntry>
        >(
            context.getDiaSourceIndex(),
            context.getObjectIndex(),
            context.getMatchRadius()/3600.0,    // match routine expects degrees, not arc-seconds
            pdf,
            pof,
            mlp
        );
        watch.stop();
        Log log(Log::getDefaultLog(), "associate");
        Rec(log, Log::INFO) << "matched difference sources to objects" <<
            DataProperty("numDiaSources", context.getDiaSourceIndex().size()) <<
            DataProperty("numObjects", context.getObjectIndex().size()) <<
            DataProperty("numMatches", static_cast<long>(nm)) <<
            DataProperty("time", watch.seconds()) << Rec::endr;

    } catch (...) {
        manager.endVisit(context.getVisitId(), true);
        throw;
    }
}


// -- MatchMop stage ----------------

/*!
    Matches moving object predictions falling within the FOV of a visit against the difference sources
    for that visit.

    \param[out]     matches     Set to a list of moving object prediction to difference source match pairs.
    \param[out]     newObjects  Set to a list of (difference source id, object id) pairs that specifies
                                which difference sources should be used to create new object and what
                                the id of each new object should be set to.
    \param[in, out] context     State involved in processing a single visit.
    \param[in]      predictions The list of moving object predictions to match against difference sources.
 */
LSST_AP_API void matchMops(
    boost::shared_ptr<MatchPairVector>     & matches,
    boost::shared_ptr<IdPairVector>        & newObjects,
    VisitProcessingContext                 & context,
    lsst::fw::MovingObjectPredictionVector & predictions
) {
    SharedSimpleObjectChunkManager manager(context.getRunId());

    try {

        matches.reset(new MatchPairVector);
        matches->reserve(8192);
        newObjects.reset(new IdPairVector);
        newObjects->reserve(8192);

        detail::MovingObjectPredictionMatchProcessor mpp(*matches);

        // discard difference sources with known variable matches
        Stopwatch watch(true);
        detail::DiscardKnownVariableFilter dvf;
        size_t nr = context.getDiaSourceIndex().pack(dvf);
        watch.stop();
        Log log(Log::getDefaultLog(), "associate");
        Rec(log, Log::INFO) << "removed difference sources matching known variables from index" <<
            DataProperty("numRemoved", static_cast<long>(nr)) <<
            DataProperty("time", watch.seconds()) << Rec::endr;

        // build ellipses required for matching from predictions
        watch.start();
        EllipseList<lsst::fw::MovingObjectPrediction> ellipses;
        ellipses.reserve(predictions.size());
        lsst::fw::MovingObjectPredictionVector::iterator const end = predictions.end();
        for (lsst::fw::MovingObjectPredictionVector::iterator i = predictions.begin(); i != end; ++i) {
            ellipses.push_back(*i);
        }
        watch.stop();
        Rec(log, Log::INFO) << "built list of match parameters for moving object predictions" <<
            DataProperty("numPredictions", static_cast<long>(ellipses.size())) <<
            DataProperty("time", watch.seconds()) << Rec::endr;

        // match them against difference sources
        detail::DiscardLargeEllipseFilter dlef;
        PassthroughFilter<detail::DiaSourceEntry> pf;
        watch.start();
        size_t nm = ellipseMatch<
            lsst::fw::MovingObjectPrediction,
            detail::DiaSourceEntry,
            detail::DiscardLargeEllipseFilter,
            PassthroughFilter<detail::DiaSourceEntry>,
            detail::MovingObjectPredictionMatchProcessor
        >(
            ellipses,
            context.getDiaSourceIndex(),
            dlef,
            pf,
            mpp
        );
        watch.stop();
        Rec(log, Log::INFO) << "matched moving object predictions to difference sources" <<
            DataProperty("numPredictions", static_cast<long>(ellipses.size())) <<
            DataProperty("numDiaSources", context.getDiaSourceIndex().size()) <<
            DataProperty("numMatches", static_cast<long>(nm)) <<
            DataProperty("time", watch.seconds()) << Rec::endr;

        // Create new objects from difference sources with no matches
        watch.start();
        detail::NewObjectCreator createObjects(*newObjects, context.getChunks(), context.getDecomposition(), context.getFilter());
        context.getDiaSourceIndex().apply(createObjects);
        watch.stop();
        Rec(log, Log::INFO) << "created new objects" <<
            DataProperty("numObjects", static_cast<long>(newObjects->size())) <<
            DataProperty("time", watch.seconds()) << Rec::endr;
    } catch (...) {
        manager.endVisit(context.getVisitId(), true);
        throw;
    }
}


// -- Store stage ----------------

/*!
    Stores any new objects that have been added to the FOV of the visit.

    \param[in, out] context     State involved in processing a single visit.
 */
LSST_AP_API void storeSliceObjects(VisitProcessingContext & context) {

    typedef VisitProcessingContext::SimpleObjectChunk ChunkType;
    typedef std::vector<ChunkType>                    ChunkVectorType;
    typedef std::vector<ChunkType>::iterator          ChunkIteratorType;

    SharedSimpleObjectChunkManager manager(context.getRunId());
    Log log(Log::getDefaultLog(), "associate");

    try {
        Stopwatch watch(true);
        ChunkToFileNameMapping deltaNames(sObjDeltaFilePattern);
        ChunkVectorType &      chunks = context.getChunks();
        ChunkIteratorType      end(chunks.end());
        for (ChunkIteratorType i(chunks.begin()); i != end; ++i) {
            ChunkType & c = *i;
            std::string file(deltaNames.getName(context.getRunId(), context.getDecomposition(), c.getId()));
            verifyPathName(file);
            c.writeDelta(file, true, false);
        }
        watch.stop();
        Rec(log, Log::INFO) << "wrote chunk delta files" <<
            DataProperty("numChunks", static_cast<long>(chunks.size())) <<
            DataProperty("time", watch.seconds()) << Rec::endr;

    } catch (lsst::mwi::exceptions::ExceptionStack & ex) {
        Rec(log, Log::FATAL) << ex.what() << *(ex.getStack()) << Rec::endr;
        manager.failVisit(context.getVisitId());
    } catch (std::exception & ex) {
        log.log(Log::FATAL, ex.what());
        manager.failVisit(context.getVisitId());
    } catch (...) {
        log.log(Log::FATAL, "caught unknown exception");
        manager.failVisit(context.getVisitId());
    }
}


/*!
    Marks processing for the given visit as a failure.

    \param[in, out] context     State involved in processing a single visit.
 */
LSST_AP_API void failVisit(VisitProcessingContext & context) {
    SharedSimpleObjectChunkManager manager(context.getRunId());
    manager.failVisit(context.getVisitId());
}


/*!
    Ends visit processing - in memory changes are rolled back if the given visit failed in any way, or
    if \a rollback is \c true .

    \param[in, out] context     State involved in processing a single visit.
    \param[in]      rollback    Indicates whether or not results for the visit should be rolled back

    \return     \c true if the visit existed, was not marked as a failure and was committed,
                \c false otherwise.
 */
LSST_AP_API bool endVisit(VisitProcessingContext & context, bool const rollback) {
    SharedSimpleObjectChunkManager manager(context.getRunId());
    bool committed = manager.endVisit(context.getVisitId(), rollback);
    Log log(Log::getDefaultLog(), "associate");
    if (committed) {
        log.log(Log::INFO, "Committed visit");
    } else {
        log.log(Log::FATAL, "Rolled back visit");
    }
    return committed;
}


}} // end of namespace lsst::ap

