// -*- lsst-c++ -*-

/**
 * @file
 * @brief   Implementation of association pipeline stages.
 *
 * @ingroup ap
 *
 * @todo    How is the time of observation for a visit specified? And how is
 *          the association pipeline deadline derived from this? Possibly post-DC2.
 *
 * @todo    [Post DC3a] How are variability probabilities set for new objects? Currently,
 *          the probability is set to 100% for any new Object.
 *
 * @todo    [Post DC3a] StoreStage can issue multiple database requests in parallel
 *
 * @todo    The end of the pipeline needs to send out a 'triggerAlertGeneration' event.
 *          Post-DC3a, since there is no alert generation pipeline.
 *
 * @todo    [Post DC3a] Add a version number to chunk delta file names. Use a database
 *          table to track the version numbers of valid delta files in transactional
 *          fashion. Note that the transaction which updates these version numbers
 *          should include the insert to the object table as well as the append to the
 *          historical DIASource table. Use MySQL API directly for this?
 */

#if LSST_HAVE_OPEN_MP
#   include <omp.h>
#endif

#include <algorithm>
#include <memory>
#include <utility>

#include "boost/format.hpp"
#include "boost/scoped_array.hpp"

#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Log.h"

#include "lsst/afw/detection/DiaSource.h"
#include "lsst/afw/image/Filter.h"
#include "lsst/mops/MovingObjectPrediction.h"

#include "lsst/ap/ChunkManager.h"
#include "lsst/ap/Match.h"
#include "lsst/ap/Point.h"
#include "lsst/ap/Stages.h"
#include "lsst/ap/Time.h"
#include "lsst/ap/Utils.h"

using lsst::daf::base::PropertySet;
using lsst::pex::logging::Log;
using lsst::pex::logging::Rec;
using lsst::pex::logging::Prop;
using lsst::pex::policy::Policy;
using lsst::daf::persistence::LogicalLocation;
using lsst::afw::detection::DiaSource;
using lsst::afw::detection::DiaSourceSet;
using lsst::afw::detection::PersistableDiaSourceVector;
using lsst::afw::image::Filter;
using lsst::mops::MovingObjectPrediction;

namespace ex = lsst::pex::exceptions;


namespace lsst { namespace ap {

// -- Constants ----------------

static boost::uint32_t const HAS_MATCH  = 1;
static boost::uint32_t const HAS_KNOWN_VARIABLE_MATCH = 2;
static boost::uint32_t const HAS_MOVING_OBJECT_MATCH = 4;
static char const * const VAR_PROB_THRESH_KEY[6] = {
    "uVarProbThreshold",
    "gVarProbThreshold",
    "rVarProbThreshold",
    "iVarProbThreshold",
    "zVarProbThreshold",
    "yVarProbThreshold"
};

// -- typedefs and templates for chunks and spatial indexes ----------------

namespace detail {

typedef SharedObjectChunkManager::ObjectChunk ObjectChunk;

typedef std::vector<ObjectChunk> ObjectChunkVector;

typedef ZoneEntry<ObjectChunk> ObjectEntry;
typedef ZoneEntry<DiaSourceChunk> DiaSourceEntry;

typedef Ellipse<MovingObjectPrediction> MovingObjectEllipse;

} // end of namespace detail

#if defined(__GNUC__) && __GNUC__ > 3
#   pragma GCC visibility push(hidden)
#endif
/// @cond

template class ZoneEntry<detail::ObjectChunk>;
template class ZoneEntry<DiaSourceChunk>;

template class ZoneEntryArray<detail::ObjectEntry>;
template class ZoneEntryArray<detail::DiaSourceEntry>;

template class ZoneIndex<detail::ObjectEntry>;
template class ZoneIndex<detail::DiaSourceEntry>;

template class Ellipse<lsst::mops::MovingObjectPrediction>;
template class EllipseList<lsst::mops::MovingObjectPrediction>;
/// @endcond
#if defined(__GNUC__) && __GNUC__ > 3
#   pragma GCC visibility pop
#endif


namespace detail {

// -- Proper motion correction for objects ----------------

/** @return the proper motion corrected position of the given object */
LSST_AP_LOCAL std::pair<double, double> correctProperMotion(Object const& obj, double const epoch) {
    static double const RAD_PER_MAS = (RADIANS_PER_DEGREE/360000.0);
    // (rad/mas)*(sec/year)/(km/AU)
    static double const SCALE = RAD_PER_MAS*(365.25*86400/149597870.691);

    double ra   = radians(obj.getRa());   // rad
    double decl = radians(obj.getDec()); // rad

    // Convert ra, dec to unit vector in cartesian coordinate system
    double const sinDecl = std::sin(decl);
    double const cosDecl = std::cos(decl);
    double const sinRa   = std::sin(ra);
    double const cosRa   = std::cos(ra);
    double x = cosRa * cosDecl;
    double y = sinRa * cosDecl;
    double z = sinDecl;

    // compute space motion vector (radians per year)
    double const pmRa   = obj.getMuRa()*RAD_PER_MAS/cosDecl;
    double const pmDecl = obj.getMuDecl()*RAD_PER_MAS;
    // divide radial velociy by distance to source
    double const w = obj.getParallax()*obj.getRadialVelocity()*SCALE;

    double const mx = - pmRa*y - pmDecl*cosRa   + x*w;
    double const my =   pmRa*x - pmDecl*sinRa   + y*w;
    double const mz =            pmDecl*cosDecl + z*w;

    // Linear interpolation of position
    double const dt = (epoch - obj.getEpoch()) * (1/365.25); // julian years
    x += mx*dt;
    y += my*dt;
    z += mz*dt;

    // Store unit vector for corrected position, convert back to
    // spherical coords and store scaled integer ra/dec
    double d2 = x*x + y*y;
    ra   = (d2 == 0.0) ? 0 : degrees(std::atan2(y, x));
    decl = (z  == 0.0) ? 0 : degrees(std::atan2(z, std::sqrt(d2)));
    if (ra < 0.0) {
        ra += 360.0;
    }
    return std::make_pair(ra, decl);
}


// -- Match processors ----------------

/** @brief  Processor for lists of objects matching a difference source */
template <typename ZoneEntryT>
class LSST_AP_LOCAL ObjectMatchProcessor {

public :

    typedef MatchWithDistance<ZoneEntryT> Match;
    typedef typename std::vector<Match>::iterator MatchIterator;

    MatchPairVector & _matches;
    Filter const _filter;
    int const _threshold;

    ObjectMatchProcessor(
        VisitProcessingContext & context,
        MatchPairVector & matches,
        Filter const filter
    ) :
        _matches(matches),
        _filter(filter),
        _threshold(context.getPipelinePolicy()->getInt(VAR_PROB_THRESH_KEY[filter]))
    {}

    void operator()(DiaSourceEntry & ds, MatchIterator begin, MatchIterator end) {
        boost::uint32_t flags = HAS_MATCH;
        do {
            ZoneEntryT * const obj  = begin->_match;
            double const dist = degrees(begin->_distance);
            ++begin;
            // record match results (to be persisted later)
            _matches.push_back(MatchPair(ds._data->getId(), obj->_data->getId(), dist));
            if (obj->_data->getVarProb(_filter) >= _threshold) {
                // flag ds as matching a known variable
                flags |= HAS_KNOWN_VARIABLE_MATCH;
            }
        } while (begin != end);
        ds._flags = flags;
    }
};


/** @brief  Processor for matches between moving object predictions and difference sources. */
struct LSST_AP_LOCAL MovingObjectPredictionMatchProcessor {

    typedef DiaSourceEntry * Match;
    typedef std::vector<Match>::iterator MatchIterator;

    MatchPairVector & _results;

    explicit MovingObjectPredictionMatchProcessor(MatchPairVector & results) : _results(results) {}

    void operator()(MovingObjectEllipse & ell, DiaSourceEntry & ds) {
        ds._flags  |= HAS_MOVING_OBJECT_MATCH;
        double dx   = ell._cosRa * ell._cosDec - ds._x;
        double dy   = ell._sinRa * ell._cosDec - ds._y;
        double dz   = ell._sinDec              - ds._z;
        double dist = 2.0*std::asin(0.5*std::sqrt(dx*dx + dy*dy + dz*dz));
        // record match results (to be persisted later)
        _results.push_back(MatchPair(ell._data->getId(), ds._data->getId(), degrees(dist)));
    }
};


// -- Zone index filters and functors ----------------

/** @brief  Filter which discards difference sources matching known variable objects. */
struct LSST_AP_LOCAL DiscardKnownVariableFilter {
    bool operator()(DiaSourceEntry const & ds) {
        return (ds._flags & HAS_KNOWN_VARIABLE_MATCH) == 0;
    }
};


/** @brief  Filter which discards predicted moving objects with large position error ellipses. */
struct LSST_AP_LOCAL DiscardLargeEllipseFilter {
    double semiMajorAxisThreshold;

    DiscardLargeEllipseFilter(Policy::Ptr const policy) :
        semiMajorAxisThreshold(policy->getDouble("semiMajorAxisThreshold")) {}

    bool operator()(lsst::mops::MovingObjectPrediction const & p) {
        return p.getSemiMajorAxisLength() < semiMajorAxisThreshold;
    }
};


/** @brief  Records ids of difference sources with no matches. */
struct LSST_AP_LOCAL NewObjectCreator {

    typedef std::map<int, ObjectChunk> ChunkMap;
    typedef ChunkMap::value_type ChunkMapValue;
    typedef ChunkMap::iterator ChunkMapIterator;

    IdPairVector & _results;
    ZoneStripeChunkDecomposition const & _zsc;
    Point const _fovCen;
    double const _fovRad;
    ChunkMap _chunks;
    lsst::afw::image::Filter const _filter;
    boost::int64_t const _idNamespace;

    NewObjectCreator(
        IdPairVector & results,
        VisitProcessingContext & context
    ) :
        _results(results),
        _zsc(context.getDecomposition()),
        _fovCen(context.getFov().getCenterRa(), context.getFov().getCenterDec()),
        _fovRad(context.getFov().getRadius()),
        _chunks(),
        _filter(context.getFilter()),
        _idNamespace(static_cast<boost::int64_t>(context.getFilter() + 1) << 56)
    {
        ObjectChunkVector & chunks = context.getChunks();
        // build a map of ids to chunks
        for (ObjectChunkVector::iterator i(chunks.begin()), end(chunks.end()); i != end; ++i) {
            _chunks.insert(ChunkMapValue(i->getId(), *i));
        }
    }

    void operator()(DiaSourceEntry const & entry) {
        // TODO - this logic should be moved into Python to make it easier to change and more configureable
        static boost::int64_t const SHAPE_DIFFERS_IN_BOTH_EXPOSURES_MASK = 1 << 2;
        static boost::int64_t const POSITIVE_FLUX_EXCURSION_MASK = (1 << 3) | (1<< 4);
        static boost::int64_t const idLimit = INT64_C(1) << 56;
        // Generate at most 1 object for a pair of difference sources
        if ((entry._data->getDiaSourceToId() & (1 << 30)) != 0) {
            return;
        }
        boost::int64_t classFlags = entry._data->getFlagClassification();
        // Don't generate objects for cosmic rays
        if (((classFlags & SHAPE_DIFFERS_IN_BOTH_EXPOSURES_MASK) != 0) &&
            ((classFlags & POSITIVE_FLUX_EXCURSION_MASK) != 0)) {
            return;
        }
        // TODO: Don't generate objects for fast movers. Requires knowledge of
        // ellipticity of difference source after PSF deconvolution, which is not
        // available for DC3a.

        if ((entry._flags & (HAS_MATCH | HAS_KNOWN_VARIABLE_MATCH)) == 0) {
            // difference source had no matches - record it as the source of a new object
            boost::int64_t id = entry._data->getId();
            if (id >= idLimit) {
                throw LSST_EXCEPT(ex::RangeErrorException, "DiaSource id doesn't fit in 56 bits");
            }
            // generate a new simplified object (id, position, proper motions, variability probabilities)
            // and assign it to the appropriate chunk
            Object obj;
            
            obj._objectId           = id | _idNamespace;
            obj._ra                 = entry._data->getRa();
            obj._decl               = entry._data->getDec();
            obj._muRa               = 0.0;
            obj._muDecl             = 0.0;
            obj._parallax           = 0.0;
            obj._radialVelocity     = 0.0; 
            obj._varProb[Filter::U] = 0;
            obj._varProb[Filter::G] = 0;
            obj._varProb[Filter::R] = 0;
            obj._varProb[Filter::I] = 0;
            obj._varProb[Filter::Z] = 0;
            obj._varProb[Filter::Y] = 0;
            obj._varProb[_filter]   = 100;

            _results.push_back(IdPair(id, obj._objectId));

            // find the chunk the new object belongs to and insert the new object into it
            int const chunkId = _zsc.radecToChunk(obj._ra, obj._decl);
            ChunkMapIterator c = _chunks.find(chunkId);
            if (c == _chunks.end()) {
                throw LSST_EXCEPT(ex::RuntimeErrorException,
                    (boost::format("new object from DIASource %1% ra,dec=(%2%, %3%) x,y=(%4%, %5%) "
                                   "not in any chunk overlapping FOV with center (%6%, %7%), "
                                   "radius=%8%: new object would go to chunk %9%; distance to FOV "
                                   "center is %10% deg") % id % obj._ra % obj._decl %
                     entry._data->getXAstrom() % entry._data->getYAstrom() % _fovCen._ra %
                     _fovCen._dec % _fovRad % chunkId % _fovCen.distance(Point(obj._ra, obj._decl))).str());
            }
            c->second.insert(obj);
        }
    }
};


// -- Index creation ----------------

template <typename EntryT>
void LSST_AP_LOCAL buildZoneIndex(
    ZoneIndex<EntryT> & index,
    std::vector<typename EntryT::Chunk> const & chunks,
    double const epoch
) {
    typedef typename EntryT::Data Data;
    typedef typename EntryT::Chunk Chunk;
    typedef std::vector<Chunk> ChunkVector;
    typedef typename ChunkVector::const_iterator ChunkIterator;
    typedef typename ChunkVector::size_type Size;

    index.clear();
    if (chunks.empty()) {
        return;
    }

    Stopwatch watch(true);

    // determine stripe bounds for the input chunks
    ZoneStripeChunkDecomposition const & zsc = index.getDecomposition();
    int minStripe = 0x7FFFFFFF;
    int maxStripe = -1 - minStripe;
    ChunkIterator const end(chunks.end());
    for (ChunkIterator c(chunks.begin()); c != end; ++c) {
        int const stripeId = ZoneStripeChunkDecomposition::chunkToStripe(c->getId());
        if (stripeId > maxStripe) {
            maxStripe = stripeId;
        }
        if (stripeId < minStripe) {
            minStripe = stripeId;
        }
    }
    assert(maxStripe >= minStripe && "invalid stripe bounds for chunk list");

    // Partition input chunks into stripes
    int const numStripes = maxStripe - minStripe + 1;
    boost::scoped_array<ChunkVector> stripes(new ChunkVector[numStripes]);
    for (ChunkIterator c(chunks.begin()); c != end; ++c) {
        int const stripeId = ZoneStripeChunkDecomposition::chunkToStripe(c->getId());
        assert(stripeId >= minStripe && stripeId <= maxStripe && "stripe id out of bounds");
        stripes[stripeId - minStripe].push_back(*c);
    }

    double minDec = zsc.getStripeDecMin(minStripe) - 0.001;
    double maxDec = zsc.getStripeDecMax(maxStripe) + 0.001;
    index.setDecBounds(std::max(minDec, -90.0), std::min(maxDec, 90.0));

    // Loop over stripes. Note that due to proper motion, objects from
    // different stripes can end up in the same zone. Objects are never
    // migrated across chunks/stripes; rather, their J2000 coordinates
    // determine the chunks that they will belong to. As time goes on,
    // the spatial extent of a chunk will grow, where the rate of growth
    // is bounded by Barnard's star with a proper motion of ~10.5 arcsec/year.
    // When intersecting the bounding circle of a FOV with chunk boundaries
    // to determine which chunks must be loaded for association, the bounding
    // circle must be padded to ensure all relevant objects are loaded.
    try {
        for (int s = 0; s < numStripes; ++s) {
            ChunkVector &  vec = stripes[s];
            Size const numChunks = vec.size();

            // Loop over chunks in stripe
            for (Size c = 0; c < numChunks; ++c) {
                Chunk * const ch = &vec[c];
                int const numBlocks  = ch->blocks();
                int i = 0;

                // loop over blocks in chunk
                for (int b = 0; b < numBlocks; ++b) {
                    int const numEntries = ch->entries(b);
                    Data * const block = ch->getBlock(b);
                    ChunkEntryFlag const * const flags = ch->getFlagBlock(b);

                    // loop over entries in block
                    for (int e = 0; e < numEntries; ++e, ++i) {
                        if ((flags[e] & Chunk::DELETED) == 0) {
                            std::pair<double, double> pos = correctProperMotion(block[e], epoch);
                            index.insert(pos.first, pos.second, &block[e], ch, i);
                        }
                    }
                }
            }
        }
    } catch(...) {
       index.clear();
       throw LSST_EXCEPT(ex::RuntimeErrorException, "Failed to build zone index"); 
    }

    watch.stop();
    Log log(Log::getDefaultLog(), "lsst.ap");
    int numElements = static_cast<int>(index.size());
    Rec(log, Log::INFO) << "inserted elements into zone index" <<
        Prop<int>("numElements", numElements) <<
        Prop<double>("time", watch.seconds()) << Rec::endr;

    // zone structure is filled, sort individual zones (on right ascension)
    watch.start();
    index.sort();
    watch.stop();
    Rec(log, Log::INFO) << "sorted zone index" <<
        Prop<int>("numElements", numElements) <<
        Prop<double>("time", watch.seconds()) << Rec::endr;
}

} // end of namespace detail


// -- Template instantiations ----------------

#if defined(__GNUC__) && __GNUC__ > 3
#   pragma GCC visibility push(hidden)
#endif
/// @cond
template class detail::ObjectMatchProcessor<detail::ObjectEntry>;

template std::size_t distanceMatch<
    detail::DiaSourceEntry,
    detail::ObjectEntry,
    PassthroughFilter<detail::DiaSourceEntry>,
    PassthroughFilter<detail::ObjectEntry>,
    detail::ObjectMatchProcessor<detail::ObjectEntry>
>(
    ZoneIndex<detail::DiaSourceEntry> &,
    ZoneIndex<detail::ObjectEntry> &,
    double const,
    PassthroughFilter<detail::DiaSourceEntry> &,
    PassthroughFilter<detail::ObjectEntry> &,
    detail::ObjectMatchProcessor<detail::ObjectEntry> &
);

template std::size_t ellipseMatch<
    MovingObjectPrediction,
    detail::DiaSourceEntry,
    PassthroughFilter<detail::MovingObjectEllipse>,
    PassthroughFilter<detail::DiaSourceEntry>,
    detail::MovingObjectPredictionMatchProcessor
>(
    EllipseList<MovingObjectPrediction> &,
    ZoneIndex<detail::DiaSourceEntry> &,
    PassthroughFilter<detail::MovingObjectEllipse> &,
    PassthroughFilter<detail::DiaSourceEntry> &,
    detail::MovingObjectPredictionMatchProcessor &
);
/// @endcond
#if defined(__GNUC__) && __GNUC__ > 3
#   pragma GCC visibility pop
#endif


// -- VisitProcessingContext ----------------

VisitProcessingContext::VisitProcessingContext(
    Policy::Ptr const policy,
    PropertySet::Ptr const event,
    std::string const & runId,
    int const workerId,
    int const numWorkers
) :
    lsst::daf::base::Citizen(typeid(*this)),
    _policy(policy),
    _chunkIds(),
    _chunks(),
    _objectIndex(policy->getInt("zonesPerDegree"),
                 policy->getInt("zonesPerStripe"),
                 policy->getInt("maxEntriesPerZoneEstimate")),
    _diaSourceIndex(policy->getInt("zonesPerDegree"),
                    policy->getInt("zonesPerStripe"),
                    policy->getInt("maxEntriesPerZoneEstimate")),
    _deadline(),
    _fov(),
    _runId(runId),
    _visitId(-1),
    _matchRadius(policy->getDouble("matchRadius")),
    _ellipseScalingFactor(policy->getDouble("ellipseScalingFactor")),
    _filter(),
    _workerId(workerId),
    _numWorkers(numWorkers),
    _debugSharedMemory(policy->getBool("debugSharedMemory"))
{
    double ra = event->getAsDouble("ra");
    double dec = event->getAsDouble("decl");
    _fov = CircularRegion(ra, dec, policy->getDouble("fovRadius"));
    _visitId = event->getAsInt("visitId");
    if (event->exists("matchRadius")) {
        _matchRadius = event->getAsDouble("matchRadius");
    }
    _visitTime = event->getAsDouble("dateObs");

    // DC3a: set association pipeline deadline to 10 minutes
    // after creation of a visit processing context.
    _deadline.systemTime();
    _deadline.tv_sec += 600;
    std::string filterName = event->getAsString("filter");
    LogicalLocation location(policy->getString("filterTableLocation"));
    _filter = Filter(location, filterName);
}


VisitProcessingContext::~VisitProcessingContext() {}


void VisitProcessingContext::setDiaSources(boost::shared_ptr<PersistableDiaSourceVector> diaSources) {
    _diaSources = diaSources->getSources();
    _diaSourceIndex.clear();
    int const sz = static_cast<int>(_diaSources.size());
    if (sz == 0) {
        return;
    }
    Stopwatch watch(true);
    double minDec =  90.0;
    double maxDec = -90.0;
    for (int i = 0; i < sz; ++i) {
        double dec = _diaSources[i]->getDec();
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
    Log log(Log::getDefaultLog(), "lsst.ap");
    Rec(log, Log::INFO) << "set dec bounds for difference source index" <<
        Prop<int>("numElements", sz) <<
        Prop<double>("time", watch.seconds()) << Rec::endr;
    watch.start();
    try {
        for (int i = 0; i < sz; ++i) {
            _diaSourceIndex.insert(_diaSources[i]->getRa(), _diaSources[i]->getDec(),
                                   _diaSources[i].get(), 0, 0);
        }
    } catch (...) {
        _diaSourceIndex.clear();
        throw;
    }
    watch.stop();
    Rec(log, Log::INFO) << "inserted difference sources into zone index" <<
        Prop<int>("numElements", sz) <<
        Prop<double>("time", watch.seconds()) << Rec::endr;
    watch.start();
    _diaSourceIndex.sort();
    watch.stop();
    Rec(log, Log::INFO) << "sorted difference source zone index" <<
        Prop<int>("numElements", sz) <<
        Prop<double>("time", watch.seconds()) << Rec::endr;
}


void VisitProcessingContext::buildObjectIndex() {
    detail::buildZoneIndex(_objectIndex, _chunks, _visitTime);
}


// -- Load stage ----------------

/**
 * Sets up all fundamental visit processing parameters using a policy and ensure
 * that a reference to the shared memory object used for chunk storage exists.
 */
LSST_AP_API void initialize(std::string const & runId) {
    // create shared memory object if it doesn't already exist
    volatile SharedObjectChunkManager manager(runId);
}


/**
 * Computes ids for all object chunks covering the visit FOV and
 * registers the visit with the shared memory chunk manager.
 */
LSST_AP_API void registerVisit(VisitProcessingContext & context) {
    context.getChunkIds().clear();
    computeChunkIds(context.getChunkIds(), context.getFov(), context.getDecomposition(), 0, 1);
    SharedObjectChunkManager manager(context.getRunId());
    manager.registerVisit(context.getVisitId());
}


/**
 * Ensures that object data for the chunks assigned to the calling slice has been read in or is
 * owned by the given visit.
 */
LSST_AP_API void loadSliceObjects(VisitProcessingContext & context) {

    typedef VisitProcessingContext::ObjectChunk Chunk;
    typedef std::vector<Chunk>                        ChunkVector;
    typedef std::vector<Chunk>::iterator              ChunkIterator;

    SharedObjectChunkManager manager(context.getRunId());
    Log log(Log::getDefaultLog(), "lsst.ap");

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
        Rec(log, Log::INFO) << "computed chunk ids in FOV for worker " <<
            Prop<int>("numChunks", static_cast<int>(context.getChunkIds().size())) <<
            Prop<double>("time", watch.seconds()) << Rec::endr;

        // Register interest in or create chunks via the chunk manager
        std::vector<Chunk> toRead;
        std::vector<Chunk> toWaitFor;
        watch.start();
        manager.startVisit(toRead, toWaitFor, context.getVisitId(), context.getChunkIds());
        watch.stop();
        Rec(log, Log::INFO) << "started processing visit" <<
            Prop<double>("time", watch.seconds()) << Rec::endr;

        // record pointers to all chunks being handled by the slice
        ChunkVector & chunks = context.getChunks();
        chunks.clear();
        chunks.insert(chunks.end(), toRead.begin(),    toRead.end());
        chunks.insert(chunks.end(), toWaitFor.begin(), toWaitFor.end());

        // Read data files
        watch.start();
        std::string refNamePattern(context.getPipelinePolicy()->getString("objectChunkFileNamePattern"));
        std::string deltaNamePattern(context.getPipelinePolicy()->getString("objectDeltaChunkFileNamePattern"));
        PropertySet::Ptr ps(new PropertySet);
        ps->set<std::string>("runId", context.getRunId());
        ChunkVector::size_type numToRead(toRead.size());
        for (ChunkIterator i(toRead.begin()), end(toRead.end()); i != end; ++i) {
            Chunk & c = *i;
            ps->set<int>("chunkId", c.getId());
            ps->set<int>("stripeId", ZoneStripeChunkDecomposition::chunkToStripe(c.getId()));
            ps->set<int>("chunkSeqNum", ZoneStripeChunkDecomposition::chunkToSequence(c.getId()));
            c.read(LogicalLocation(refNamePattern, ps).locString(), false);
            c.readDelta(LogicalLocation(deltaNamePattern, ps).locString(), false);
            c.setUsable();
        }
        watch.stop();
        Rec(log, Log::INFO) << "read chunk files" <<
            Prop<int>("numChunks", static_cast<int>(numToRead)) <<
            Prop<double>("time", watch.seconds()) << Rec::endr;

        toRead.clear();
        ChunkVector::size_type numToWaitFor = toWaitFor.size();
        if (numToWaitFor > 0) {
            watch.start();
            // Wait for chunks that are owned by another visit
            manager.waitForOwnership(toRead, toWaitFor, context.getVisitId(), context.getDeadline());
            watch.stop();
            Rec(log, Log::INFO) << "acquired ownership of pre-existing chunks" <<
                Prop<int>("numChunks", static_cast<int>(numToWaitFor)) <<
                Prop<double>("time", watch.seconds()) << Rec::endr;

            // Read in chunks that were not successfully read by the previous owner
            watch.start();
            numToRead = toRead.size();
            for (ChunkIterator i(toRead.begin()), end(toRead.end()); i != end; ++i) {
                Chunk & c = *i;
                ps->set<int>("chunkId", c.getId());
                ps->set<int>("stripeId", ZoneStripeChunkDecomposition::chunkToStripe(c.getId()));
                ps->set<int>("chunkSeqNum", ZoneStripeChunkDecomposition::chunkToSequence(c.getId()));
                c.read(LogicalLocation(refNamePattern, ps).locString(), false);
                c.readDelta(LogicalLocation(deltaNamePattern, ps).locString(), false);
                c.setUsable();
            }
            watch.stop();
            Rec(log, Log::INFO) << "read straggling chunks" <<
                Prop<int>("numChunks", static_cast<int>(numToRead)) <<
                Prop<double>("time", watch.seconds()) << Rec::endr;
        }

    } catch (ex::Exception & except) {
        Rec(log, Log::FATAL) << except.what() << Rec::endr;
        manager.failVisit(context.getVisitId());
    } catch (std::exception & except) {
        log.log(Log::FATAL, except.what());
        manager.failVisit(context.getVisitId());
    } catch (...) {
        log.log(Log::FATAL, "caught unknown exception");
        manager.failVisit(context.getVisitId());
    }
}


/**
 * Checks to see if all parallel workers succeeded in loading their chunks and, if so, builds
 * a zone index for the objects just loaded.
 *
 * @param[in, out] context  State involved in processing a single visit.
 */
LSST_AP_API void buildObjectIndex(VisitProcessingContext & context) {
    if (!context.debugSharedMemory()) {
        // if the shared memory object used for chunk storage hasn't yet been unlinked, do so now
        SharedObjectChunkManager::destroyInstance(context.getRunId());
    }
    SharedObjectChunkManager manager(context.getRunId());
    if (manager.isVisitInFlight(context.getVisitId())) {
        try {
            // Build zone index on objects
            manager.getChunks(context.getChunks(), context.getChunkIds());
            context.buildObjectIndex();
        } catch(...) {
            manager.endVisit(context.getVisitId(), true);
            throw;
        }
    } else {
        // One or more workers failed in the load phase - rollback the visit
        manager.endVisit(context.getVisitId(), true);
        throw LSST_EXCEPT(ex::RuntimeErrorException,
                          "Association pipeline failed to read Object data for FOV");
    }
}


// -- MatchDiaSource stage ----------------

/**
 * Matches difference sources for a visit (obtained from the detection pipeline) against the objects
 * in the visit FOV.
 *
 * @param[out]     matches  Set to a list of difference source to object match pairs.
 * @param[in, out] context  State involved in processing a single visit.
 */
LSST_AP_API void matchDiaSources(
    MatchPairVector & matches,
    VisitProcessingContext & context
) {
    SharedObjectChunkManager manager(context.getRunId());

    try {

        matches.clear();
        matches.reserve(65536);

        detail::ObjectMatchProcessor<detail::ObjectEntry> mlp(context, matches, context.getFilter());
        PassthroughFilter<detail::DiaSourceEntry> pdf;
        PassthroughFilter<detail::ObjectEntry> pof;

        Stopwatch watch(true);
        std::size_t nm = distanceMatch<
            detail::DiaSourceEntry,
            detail::ObjectEntry,
            PassthroughFilter<detail::DiaSourceEntry>,
            PassthroughFilter<detail::ObjectEntry>,
            detail::ObjectMatchProcessor<detail::ObjectEntry>
        >(
            context.getDiaSourceIndex(),
            context.getObjectIndex(),
            context.getMatchRadius()/3600.0,    // match routine expects degrees, not arc-seconds
            pdf,
            pof,
            mlp
        );
        watch.stop();
        Log log(Log::getDefaultLog(), "lsst.ap");
        Rec(log, Log::INFO) << "matched difference sources to objects" <<
            Prop<int>("numDiaSources", context.getDiaSourceIndex().size()) <<
            Prop<int>("numObjects", context.getObjectIndex().size()) <<
            Prop<int>("numMatches", static_cast<int>(nm)) <<
            Prop<double>("time", watch.seconds()) << Rec::endr;

    } catch (...) {
        manager.endVisit(context.getVisitId(), true);
        throw;
    }
}


// -- MatchMop stage ----------------

/**
 * Matches moving object predictions falling within the FOV of a visit against the difference sources
 * for that visit.
 *
 * @param[out]     matches     Set to a list of moving object prediction to difference source match pairs.
 * @param[out]     newObjects  Set to a list of (difference source id, object id) pairs that specifies
 *                             which difference sources should be used to create new object and what
 *                             the id of each new object should be set to.
 * @param[in, out] context     State involved in processing a single visit.
 * @param[in]      predictions The list of moving object predictions to match against difference sources.
 */
LSST_AP_API void matchMops(
    MatchPairVector & matches,
    IdPairVector & newObjects,
    VisitProcessingContext & context,
    lsst::mops::MovingObjectPredictionVector & predictions
) {
    SharedObjectChunkManager manager(context.getRunId());

    try {
        double const ellScale = context.getPipelinePolicy()->getDouble("ellipseScalingFactor");
        double const smaaClamp = context.getPipelinePolicy()->getDouble("semiMajorAxisClamp");
        double const smiaClamp = context.getPipelinePolicy()->getDouble("semiMinorAxisClamp");
        
        matches.clear();
        matches.reserve(8192);
        newObjects.clear();
        newObjects.reserve(8192);

        detail::MovingObjectPredictionMatchProcessor mpp(matches);

        // discard difference sources with known variable matches
        Stopwatch watch(true);
        detail::DiscardKnownVariableFilter dvf;
        int nr = context.getDiaSourceIndex().pack(dvf);
        watch.stop();
        Log log(Log::getDefaultLog(), "lsst.ap");
        Rec(log, Log::INFO) << "removed difference sources matching known variables from index" <<
            Prop<int>("numRemoved", nr) <<
            Prop<double>("time", watch.seconds()) << Rec::endr;

        // build ellipses required for matching from predictions
        watch.start();
        EllipseList<lsst::mops::MovingObjectPrediction::MovingObjectPrediction> ellipses;
        ellipses.reserve(predictions.size());
        detail::DiscardLargeEllipseFilter elf(context.getPipelinePolicy());
        lsst::mops::MovingObjectPredictionVector::iterator const end = predictions.end();
        for (lsst::mops::MovingObjectPredictionVector::iterator i = predictions.begin(); i != end; ++i) {
            i->setSemiMajorAxisLength(i->getSemiMajorAxisLength() * ellScale);
            i->setSemiMinorAxisLength(i->getSemiMinorAxisLength() * ellScale);
            if (elf(*i)) {
                // clamp error ellipses if necessary
                if (smaaClamp > 0.0 && i->getSemiMajorAxisLength() > smaaClamp) {
                    i->setSemiMajorAxisLength(smaaClamp);
                }
                if (smiaClamp > 0.0 && i->getSemiMinorAxisLength() > smiaClamp) {
                    i->setSemiMinorAxisLength(smiaClamp);
                }
                ellipses.push_back(Ellipse<lsst::mops::MovingObjectPrediction>(*i));
            }
        }
        watch.stop();
        Rec(log, Log::INFO) << "built list of match parameters for moving object predictions" <<
            Prop<int>("numPredictions", static_cast<int>(ellipses.size())) <<
            Prop<double>("time", watch.seconds()) << Rec::endr;

        // match them against difference sources
        PassthroughFilter<detail::MovingObjectEllipse> pef;
        PassthroughFilter<detail::DiaSourceEntry>      pdf;
        watch.start();
        std::size_t nm = ellipseMatch<
            MovingObjectPrediction,
            detail::DiaSourceEntry,
            PassthroughFilter<detail::MovingObjectEllipse>,
            PassthroughFilter<detail::DiaSourceEntry>,
            detail::MovingObjectPredictionMatchProcessor
        >(
            ellipses,
            context.getDiaSourceIndex(),
            pef,
            pdf,
            mpp
        );
        watch.stop();
        Rec(log, Log::INFO) << "matched moving object predictions to difference sources" <<
            Prop<int>("numPredictions", static_cast<int>(ellipses.size())) <<
            Prop<int>("numDiaSources", context.getDiaSourceIndex().size()) <<
            Prop<int>("numMatches", static_cast<int>(nm)) <<
            Prop<double>("time", watch.seconds()) << Rec::endr;

        // Create new objects from difference sources with no matches
        watch.start();
        detail::NewObjectCreator createObjects(newObjects, context);
        context.getDiaSourceIndex().apply(createObjects);
        watch.stop();
        Rec(log, Log::INFO) << "created new objects" <<
            Prop<int>("numObjects", static_cast<int>(newObjects.size())) <<
            Prop<double>("time", watch.seconds()) << Rec::endr;
    } catch (...) {
        manager.endVisit(context.getVisitId(), true);
        throw;
    }
}


// -- Store stage ----------------

/**
 * Stores any new objects that have been added to the FOV of the visit.
 *
 * @param[in, out] context  State involved in processing a single visit.
 */
LSST_AP_API void storeSliceObjects(VisitProcessingContext & context) {

    typedef VisitProcessingContext::ObjectChunk Chunk;
    typedef std::vector<Chunk> ChunkVector;
    typedef std::vector<Chunk>::iterator ChunkIterator;

    SharedObjectChunkManager manager(context.getRunId());
    Log log(Log::getDefaultLog(), "lsst.ap");
    try {
        Stopwatch watch(true);
        std::string deltaNamePattern = context.getPipelinePolicy()->getString("objectDeltaChunkFileNamePattern");
        PropertySet::Ptr ps(new PropertySet);
        ps->set<std::string>("runId", context.getRunId());
        ChunkVector & chunks = context.getChunks();
        for (ChunkIterator i(chunks.begin()), end(chunks.end()); i != end; ++i) {
            Chunk & c = *i;
            ps->set<int>("chunkId", c.getId());
            ps->set<int>("stripeId", ZoneStripeChunkDecomposition::chunkToStripe(c.getId()));
            ps->set<int>("chunkSeqNum", ZoneStripeChunkDecomposition::chunkToSequence(c.getId()));
            std::string file = LogicalLocation(deltaNamePattern, ps).locString();
            verifyPathName(file);
            c.writeDelta(file, true, false);
        }
        watch.stop();
        Rec(log, Log::INFO) << "wrote chunk delta files" <<
            Prop<int>("numChunks", static_cast<int>(chunks.size())) <<
            Prop<double>("time", watch.seconds()) << Rec::endr;
    } catch (ex::Exception & except) {
        Rec(log, Log::FATAL) << except.what() << Rec::endr;
        manager.failVisit(context.getVisitId());
    } catch (std::exception & except) {
        log.log(Log::FATAL, except.what());
        manager.failVisit(context.getVisitId());
    } catch (...) {
        log.log(Log::FATAL, "caught unknown exception");
        manager.failVisit(context.getVisitId());
    }
}


/**
 * Marks processing for the given visit as a failure.
 *
 * @param[in, out] context  State involved in processing a single visit.
 */
LSST_AP_API void failVisit(VisitProcessingContext & context) {
    SharedObjectChunkManager manager(context.getRunId());
    manager.failVisit(context.getVisitId());
}


/**
 * Ends visit processing - in memory changes are rolled back if the given visit failed in any way, or
 * if @a rollback is @c true .
 *
 * @param[in, out] context  State involved in processing a single visit.
 * @param[in]      rollback Indicates whether or not results for the visit should be rolled back
 *
 * @return  @c true if the visit existed, was not marked as a failure and was committed,
 *          @c false otherwise.
 */
LSST_AP_API bool endVisit(VisitProcessingContext & context, bool const rollback) {
    SharedObjectChunkManager manager(context.getRunId());
    bool committed = manager.endVisit(context.getVisitId(), rollback);
    Log log(Log::getDefaultLog(), "lsst.ap");
    if (committed) {
        log.log(Log::INFO, "Committed visit");
    } else {
        log.log(Log::FATAL, "Rolled back visit");
    }
    return committed;
}


}} // end of namespace lsst::ap

