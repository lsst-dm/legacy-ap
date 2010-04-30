// -*- lsst-c++ -*-
/** @file
  * @brief Formatter implementations for persistable classes in lsst::ap::cluster.
  *
  * @ingroup ap
  * @author Serge Monkewitz
  */
#include "lsst/ap/cluster/Formatters.h"

#include "boost/make_shared.hpp"
#include "boost/serialization/utility.hpp"
#include "boost/serialization/collections_save_imp.hpp"
#include "boost/serialization/collections_load_imp.hpp"
#include "boost/serialization/shared_ptr.hpp"
#include "boost/serialization/split_free.hpp"
#include "boost/serialization/vector.hpp"

#include "lsst/pex/exceptions.h"
#include "lsst/daf/persistence/BoostStorage.h"
#include "lsst/daf/persistence/DbStorage.h"
#include "lsst/daf/persistence/DbTsvStorage.h"
#include "lsst/daf/persistence/FormatterImpl.h"
#include "lsst/daf/persistence/LogicalLocation.h"
#include "lsst/afw/image/Filter.h"

#include "lsst/ap/Common.h"
#include "lsst/ap/cluster/SourceCluster.h"


namespace except = lsst::pex::exceptions;

using lsst::daf::base::Persistable;
using lsst::daf::base::PropertySet;
using lsst::daf::persistence::BoostStorage;
using lsst::daf::persistence::DbStorage;
using lsst::daf::persistence::DbTsvStorage;
using lsst::daf::persistence::Formatter;
using lsst::daf::persistence::FormatterRegistration;
using lsst::daf::persistence::LogicalLocation;
using lsst::daf::persistence::Storage;
using lsst::pex::policy::Policy;
using lsst::afw::image::Filter;


namespace lsst { namespace ap { namespace cluster {

// -- Helper functions ----

namespace {

/** @internal
  * Returns the name of the database table to read/write to.
  */
LSST_AP_LOCAL std::string const getTableName(Policy::Ptr policy,
                                             PropertySet::Ptr additionalData)
{
    std::string itemName = additionalData->getAsString("itemName");
    std::string pattern = policy->getString(itemName + ".tableNamePattern");
    PropertySet::Ptr props = additionalData->deepCopy();
    return LogicalLocation(pattern, props).locString();
}

/** @internal
  * Returns the name of the template database table to copy.
  */
LSST_AP_LOCAL std::string const getTemplateTableName(Policy::Ptr policy,
                                                     PropertySet::Ptr additionalData)
{
    std::string itemName = additionalData->getAsString("itemName");
    return policy->getString(itemName + ".templateTableName");
}

template <typename StorageT, typename T>
inline void insertFloat(StorageT & db, char const * const col, T const & val) {
    if (isNaN(val)) {
        db.setColumnToNull(col);
    } else {
        db.template setColumn<T>(col, val);
    }
}

template <typename StorageT, typename T>
inline void insertFloat(StorageT & db, char const * const col, NullOr<T> const & val) {
    insertFloat(db, col, static_cast<T>(val));
}

template <typename FloatT> inline FloatT radians(FloatT deg) {
    return static_cast<FloatT>(deg * (M_PI/180.0));
}

template <typename FloatT> inline FloatT degrees(FloatT rad) {
    return static_cast<FloatT>(rad * (180.0/M_PI));
}

template <typename FloatT> inline FloatT rangeReducedDegrees(FloatT rad) {
    FloatT deg = std::fmod(degrees(rad), static_cast<FloatT>(360));
    return deg < static_cast<FloatT>(0) ? deg + static_cast<FloatT>(360) : deg;
}

/** @internal
  * Stores a mapping between 0, 1, 2, 3, 4, 5 and the ids of the filters
  * named "u", "g", "r", "i", "z", "y" | "i2" in @c filterIds. (If no filter
  * named "y" is found, "i2" is used instead).
  */
LSST_AP_LOCAL void getFilterIds(int * const filterIds) {
    filterIds[0] = Filter("u").getId();
    filterIds[1] = Filter("g").getId();
    filterIds[2] = Filter("r").getId();
    filterIds[3] = Filter("i").getId();
    filterIds[4] = Filter("z").getId();
    filterIds[5] = Filter("y").getId();
    if (filterIds[5] == Filter::UNKNOWN) {
        filterIds[5] = Filter("i2").getId();
    }
}

} // namespace


// -- serialize() functions required by boost::serialization ----
//
// Re-use STL machinery in boost::serialization to persist
// std::tr1::unordered_map. This can't be in an anonymous
// namespace (template name lookups fail).

template <typename Archive, typename FloatT>
inline void serializeFloat(Archive & ar, FloatT value) {
    bool null = isNaN(value);
    ar & null;
    if (null) {
        ar & value;
    }
}

template<class Archive>
inline void save(
    Archive & ar,
    std::tr1::unordered_map<int, PerFilterSourceClusterAttributes> const & map,
    unsigned int const)
{
    boost::serialization::stl::save_collection<
        Archive,
        std::tr1::unordered_map<int, PerFilterSourceClusterAttributes>
    >(ar, map);                                  
}   

template<class Archive>
inline void load(
    Archive & ar,
    std::tr1::unordered_map<int, PerFilterSourceClusterAttributes> & map,
    unsigned int const)
{
    boost::serialization::stl::load_collection<
        Archive,
        std::tr1::unordered_map<int, PerFilterSourceClusterAttributes>,
        boost::serialization::stl::archive_input_map<
            Archive, std::tr1::unordered_map<int, PerFilterSourceClusterAttributes>
        >,
        boost::serialization::stl::no_reserve_imp<
            std::tr1::unordered_map<int, PerFilterSourceClusterAttributes>
        >
    >(ar, map);
}

template <typename Archive>
inline void serialize(
    Archive & ar,
    std::tr1::unordered_map<int, PerFilterSourceClusterAttributes> & map,
    unsigned int const version)
{
    boost::serialization::split_free(ar, map, version);
}

template <typename FloatT>
    template <typename Archive>
void NullOr<FloatT>::serialize(Archive & ar, unsigned int const) {
    serializeFloat(ar, _value);
}

template <typename Archive>
void PerFilterSourceClusterAttributes::serialize(Archive & ar,
                                                 unsigned int const) {
    ar & _filterId;
    ar & _numObs;
    serializeFloat(ar, _earliestObsTime);
    serializeFloat(ar, _latestObsTime);
    ar & _flags;
    ar & _flux;
    ar & _fluxSigma;
    ar & _e1;
    ar & _e2;
    ar & _radius;
    ar & _e1Sigma;
    ar & _e2Sigma;
    ar & _radiusSigma;
}

template <typename Archive>
void SourceClusterAttributes::serialize(Archive & ar, unsigned int const) {
    ar & _clusterId;
    ar & _numObs;
    ar & _flags;
    serializeFloat(ar, _earliestObsTime);
    serializeFloat(ar, _latestObsTime);
    serializeFloat(ar, _ra);
    serializeFloat(ar, _dec);
    ar & _raSigma;
    ar & _decSigma;
    ar & _raDecCov;
    ar & _perFilterAttributes;
}


// -- Helpers to deal with the Object database table ----

namespace {

template <typename DbStorageT>
void insertObjectRow(DbStorageT & db,
                     SourceClusterAttributes const & attributes,
                     int const * const filterIds)
{
    // Set uncomputed columns to null
    db.setColumnToNull("iauId");
    db.setColumnToNull("muRa_PS");
    db.setColumnToNull("muRa_PS_Sigma");
    db.setColumnToNull("muDecl_PS");
    db.setColumnToNull("muDecl_PS_Sigma");
    db.setColumnToNull("muRaDecl_PS_Cov");
    db.setColumnToNull("ra_SG");
    db.setColumnToNull("ra_SG_Sigma");
    db.setColumnToNull("decl_SG");
    db.setColumnToNull("decl_SG_Sigma");
    db.setColumnToNull("radecl_SG_Cov");
    db.setColumnToNull("raRange");
    db.setColumnToNull("declRange");
    db.setColumnToNull("parallax_PS");
    db.setColumnToNull("parallax_PS_Sigma");
    db.setColumnToNull("cannonicalFilterId");
    db.setColumnToNull("extendedness");
    db.setColumnToNull("varProb");

#define LSST_AP_PF_ALWAYS_NULL(filter) \
    do { \
        db.setColumnToNull(#filter "Extendedness"); \
        db.setColumnToNull(#filter "VarProb"); \
        db.setColumnToNull(#filter "RaOffset_PS"); \
        db.setColumnToNull(#filter "RaOffset_PS_Sigma"); \
        db.setColumnToNull(#filter "DeclOffset_PS"); \
        db.setColumnToNull(#filter "DeclOffset_PS_Sigma"); \
        db.setColumnToNull(#filter "RaDeclOffset_PS_Cov"); \
        db.setColumnToNull(#filter "RaOffset_SG"); \
        db.setColumnToNull(#filter "RaOffset_SG_Sigma"); \
        db.setColumnToNull(#filter "DeclOffset_SG"); \
        db.setColumnToNull(#filter "DeclOffset_SG_Sigma"); \
        db.setColumnToNull(#filter "RaDeclOffset_SG_Cov"); \
        db.setColumnToNull(#filter "LnL_PS"); \
        db.setColumnToNull(#filter "LnL_SG"); \
        db.setColumnToNull(#filter "Flux_SG"); \
        db.setColumnToNull(#filter "Flux_SG_Sigma"); \
        db.setColumnToNull(#filter "Flux_CSG"); \
        db.setColumnToNull(#filter "Flux_CSG_Sigma"); \
        db.setColumnToNull(#filter "Timescale"); \
        db.setColumnToNull(#filter "SersicN_SG"); \
        db.setColumnToNull(#filter "SersicN_SG_Sigma"); \
    } while (false)

    LSST_AP_PF_ALWAYS_NULL(u);
    LSST_AP_PF_ALWAYS_NULL(g);
    LSST_AP_PF_ALWAYS_NULL(r);
    LSST_AP_PF_ALWAYS_NULL(i);
    LSST_AP_PF_ALWAYS_NULL(z);
    LSST_AP_PF_ALWAYS_NULL(y);

#undef LSST_AP_PF_ALWAYS_NULL

    // set filter-agnostic columns
    db.template setColumn<int64_t>("objectId", attributes.getClusterId());
    insertFloat(db, "ra_PS", rangeReducedDegrees(attributes.getRa()));
    insertFloat(db, "ra_PS_Sigma", degrees(attributes.getRaSigma()));
    insertFloat(db, "decl_PS", degrees(attributes.getDec()));
    insertFloat(db, "decl_PS_Sigma", degrees(attributes.getDecSigma()));
    insertFloat(db, "radecl_PS_Cov", degrees(attributes.getRaDecCov()));
    insertFloat(db, "earliestObsTime", attributes.getEarliestObsTime());
    insertFloat(db, "latestObsTime", attributes.getLatestObsTime());
    db.template setColumn<int>("flags", attributes.getFlags());

    // set filter-specific columns
#define LSST_AP_INSERT_FILTER(filter, index) \
    do { \
        int const id = filterIds[index]; \
        if (attributes.hasFilter(id)) { \
            PerFilterSourceClusterAttributes const & pfa = \
                attributes.getPerFilterAttributes(id); \
            db.template setColumn<int>(#filter "NumObs", pfa.getNumObs()); \
            db.template setColumn<int>(#filter "Flags", pfa.getFlags()); \
            insertFloat(db, #filter "EarliestObsTime", pfa.getEarliestObsTime()); \
            insertFloat(db, #filter "LatestObsTime", pfa.getLatestObsTime()); \
            insertFloat(db, #filter "Flux_PS", pfa.getFlux()); \
            insertFloat(db, #filter "Flux_PS_Sigma", pfa.getFluxSigma()); \
            insertFloat(db, #filter "E1_SG", pfa.getE1()); \
            insertFloat(db, #filter "E1_SG_Sigma", pfa.getE1Sigma()); \
            insertFloat(db, #filter "E2_SG", pfa.getE2()); \
            insertFloat(db, #filter "E2_SG_Sigma", pfa.getE2Sigma()); \
            insertFloat(db, #filter "Radius_SG", pfa.getRadius()); \
            insertFloat(db, #filter "Radius_SG_Sigma", pfa.getRadiusSigma()); \
        } else { \
            db.setColumnToNull(#filter "NumObs"); \
            db.setColumnToNull(#filter "Flags"); \
            db.setColumnToNull(#filter "EarliestObsTime"); \
            db.setColumnToNull(#filter "LatestObsTime"); \
            db.setColumnToNull(#filter "Flux_PS"); \
            db.setColumnToNull(#filter "Flux_PS_Sigma"); \
            db.setColumnToNull(#filter "E1_SG"); \
            db.setColumnToNull(#filter "E1_SG_Sigma"); \
            db.setColumnToNull(#filter "E2_SG"); \
            db.setColumnToNull(#filter "E2_SG_Sigma"); \
            db.setColumnToNull(#filter "Radius_SG"); \
            db.setColumnToNull(#filter "Radius_SG_Sigma"); \
        } \
    } while (false)

    LSST_AP_INSERT_FILTER(u, 0);
    LSST_AP_INSERT_FILTER(g, 1);
    LSST_AP_INSERT_FILTER(r, 2);
    LSST_AP_INSERT_FILTER(i, 3);
    LSST_AP_INSERT_FILTER(z, 4);
    LSST_AP_INSERT_FILTER(y, 5);

#undef LSST_AP_INSERT_FILTER
    // Note - do not set _chunkId or _subChunkId. These are added
    // by the partitioner.
}

/** @internal
  * A simple C++ realization of an Object database table (row). The formatter
  * uses this type to read Object rows from the database, then translates to a
  * SourceClusterAttributes.
  */
struct LSST_AP_LOCAL ObjectRow
{
    int64_t objectId;
    double  earliestObsTime;
    double  latestObsTime;
    double  ra_PS;
    double  decl_PS;
    float   ra_PS_Sigma;
    float   decl_PS_Sigma;
    float   radecl_PS_Cov;
    int     flags;
    double  uEarliestObsTime;
    double  uLatestObsTime;
    int     uNumObs;
    int     uFlags;
    float   uFlux_PS;
    float   uFlux_PS_Sigma;
    float   uE1_SG;
    float   uE1_SG_Sigma;
    float   uE2_SG;
    float   uE2_SG_Sigma;
    float   uRadius_SG;
    float   uRadius_SG_Sigma;
    double  gEarliestObsTime;
    double  gLatestObsTime;
    int     gNumObs;
    int     gFlags;
    float   gFlux_PS;
    float   gFlux_PS_Sigma;
    float   gE1_SG;
    float   gE1_SG_Sigma;
    float   gE2_SG;
    float   gE2_SG_Sigma;
    float   gRadius_SG;
    float   gRadius_SG_Sigma;
    double  rEarliestObsTime;
    double  rLatestObsTime;
    int     rNumObs;
    int     rFlags;
    float   rFlux_PS;
    float   rFlux_PS_Sigma;
    float   rE1_SG;
    float   rE1_SG_Sigma;
    float   rE2_SG;
    float   rE2_SG_Sigma;
    float   rRadius_SG;
    float   rRadius_SG_Sigma;
    double  iEarliestObsTime;
    double  iLatestObsTime;
    int     iNumObs;
    int     iFlags;
    float   iFlux_PS;
    float   iFlux_PS_Sigma;
    float   iE1_SG;
    float   iE1_SG_Sigma;
    float   iE2_SG;
    float   iE2_SG_Sigma;
    float   iRadius_SG;
    float   iRadius_SG_Sigma;
    double  zEarliestObsTime;
    double  zLatestObsTime;
    int     zNumObs;
    int     zFlags;
    float   zFlux_PS;
    float   zFlux_PS_Sigma;
    float   zE1_SG;
    float   zE1_SG_Sigma;
    float   zE2_SG;
    float   zE2_SG_Sigma;
    float   zRadius_SG;
    float   zRadius_SG_Sigma;
    double  yEarliestObsTime;
    double  yLatestObsTime;
    int     yNumObs;
    int     yFlags;
    float   yFlux_PS;
    float   yFlux_PS_Sigma;
    float   yE1_SG;
    float   yE1_SG_Sigma;
    float   yE2_SG;
    float   yE2_SG_Sigma;
    float   yRadius_SG;
    float   yRadius_SG_Sigma;

    ObjectRow() { }

    template <typename DbStorageT> void setupFetch(DbStorageT * db);

    void to(DbStorage * db,
            SourceClusterAttributes & attributes,
            int const * const filterIds);
};

template <typename DbStorageT>
void ObjectRow::setupFetch(DbStorageT * db)
{
    db->outParam("objectId", &objectId);
    db->outParam("earliestObsTime", &earliestObsTime);
    db->outParam("latestObsTime", &latestObsTime);
    db->outParam("ra_PS", &ra_PS);
    db->outParam("decl_PS", &decl_PS);
    db->outParam("ra_PS_Sigma", &ra_PS_Sigma);
    db->outParam("decl_PS_Sigma", &decl_PS_Sigma);
    db->outParam("radecl_PS_Cov", &radecl_PS_Cov);
    db->outParam("flags", &flags);
    db->outParam("uEarliestObsTime", &uEarliestObsTime);
    db->outParam("uLatestObsTime", &uLatestObsTime);
    db->outParam("uNumObs", &uNumObs);
    db->outParam("uFlags", &uFlags);
    db->outParam("uFlux_PS", &uFlux_PS);
    db->outParam("uFlux_PS_Sigma", &uFlux_PS_Sigma);
    db->outParam("uE1_SG", &uE1_SG);
    db->outParam("uE1_SG_Sigma", &uE1_SG_Sigma);
    db->outParam("uE2_SG", &uE2_SG);
    db->outParam("uE2_SG_Sigma", &uE2_SG_Sigma);
    db->outParam("uRadius_SG", &uRadius_SG);
    db->outParam("uRadius_SG_Sigma", &uRadius_SG_Sigma);
    db->outParam("gEarliestObsTime", &gEarliestObsTime);
    db->outParam("gLatestObsTime", &gLatestObsTime);
    db->outParam("gNumObs", &gNumObs);
    db->outParam("gFlags", &gFlags);
    db->outParam("gFlux_PS", &gFlux_PS);
    db->outParam("gFlux_PS_Sigma", &gFlux_PS_Sigma);
    db->outParam("gE1_SG", &gE1_SG);
    db->outParam("gE1_SG_Sigma", &gE1_SG_Sigma);
    db->outParam("gE2_SG", &gE2_SG);
    db->outParam("gE2_SG_Sigma", &gE2_SG_Sigma);
    db->outParam("gRadius_SG", &gRadius_SG);
    db->outParam("gRadius_SG_Sigma", &gRadius_SG_Sigma);
    db->outParam("rEarliestObsTime", &rEarliestObsTime);
    db->outParam("rLatestObsTime", &rLatestObsTime);
    db->outParam("rNumObs", &rNumObs);
    db->outParam("rFlags", &rFlags);
    db->outParam("rFlux_PS", &rFlux_PS);
    db->outParam("rFlux_PS_Sigma", &rFlux_PS_Sigma);
    db->outParam("rE1_SG", &rE1_SG);
    db->outParam("rE1_SG_Sigma", &rE1_SG_Sigma);
    db->outParam("rE2_SG", &rE2_SG);
    db->outParam("rE2_SG_Sigma", &rE2_SG_Sigma);
    db->outParam("rRadius_SG", &rRadius_SG);
    db->outParam("rRadius_SG_Sigma", &rRadius_SG_Sigma);
    db->outParam("iEarliestObsTime", &iEarliestObsTime);
    db->outParam("iLatestObsTime", &iLatestObsTime);
    db->outParam("iNumObs", &iNumObs);
    db->outParam("iFlags", &iFlags);
    db->outParam("iFlux_PS", &iFlux_PS);
    db->outParam("iFlux_PS_Sigma", &iFlux_PS_Sigma);
    db->outParam("iE1_SG", &iE1_SG);
    db->outParam("iE1_SG_Sigma", &iE1_SG_Sigma);
    db->outParam("iE2_SG", &iE2_SG);
    db->outParam("iE2_SG_Sigma", &iE2_SG_Sigma);
    db->outParam("iRadius_SG", &iRadius_SG);
    db->outParam("iRadius_SG_Sigma", &iRadius_SG_Sigma);
    db->outParam("zEarliestObsTime", &zEarliestObsTime);
    db->outParam("zLatestObsTime", &zLatestObsTime);
    db->outParam("zNumObs", &zNumObs);
    db->outParam("zFlags", &zFlags);
    db->outParam("zFlux_PS", &zFlux_PS);
    db->outParam("zFlux_PS_Sigma", &zFlux_PS_Sigma);
    db->outParam("zE1_SG", &zE1_SG);
    db->outParam("zE1_SG_Sigma", &zE1_SG_Sigma);
    db->outParam("zE2_SG", &zE2_SG);
    db->outParam("zE2_SG_Sigma", &zE2_SG_Sigma);
    db->outParam("zRadius_SG", &zRadius_SG);
    db->outParam("zRadius_SG_Sigma", &zRadius_SG_Sigma);
    db->outParam("yEarliestObsTime", &yEarliestObsTime);
    db->outParam("yLatestObsTime", &yLatestObsTime);
    db->outParam("yNumObs", &yNumObs);
    db->outParam("yFlags", &yFlags);
    db->outParam("yFlux_PS", &yFlux_PS);
    db->outParam("yFlux_PS_Sigma", &yFlux_PS_Sigma);
    db->outParam("yE1_SG", &yE1_SG);
    db->outParam("yE1_SG_Sigma", &yE1_SG_Sigma);
    db->outParam("yE2_SG", &yE2_SG);
    db->outParam("yE2_SG_Sigma", &yE2_SG_Sigma);
    db->outParam("yRadius_SG", &yRadius_SG);
    db->outParam("yRadius_SG_Sigma", &yRadius_SG_Sigma);
}

template <typename FloatT> void nullToNaN(DbStorage * db, int col, FloatT & field) {
    if (db->columnIsNull(col)) {
        field = std::numeric_limits<FloatT>::quiet_NaN();
    }
}

void ObjectRow::to(DbStorage * db,
                   SourceClusterAttributes & attributes,
                   int const * const filterIds)
{
    if (db->columnIsNull(0)) {
        objectId = -1; 
    }
    nullToNaN(db, 1, earliestObsTime);
    nullToNaN(db, 2, latestObsTime);
    nullToNaN(db, 3, ra_PS);
    nullToNaN(db, 4, decl_PS);
    nullToNaN(db, 5, ra_PS_Sigma);
    nullToNaN(db, 6, decl_PS_Sigma);
    nullToNaN(db, 7, radecl_PS_Cov);
    if (db->columnIsNull(8)) {
        flags = 0;
    }
    attributes.setClusterId(objectId);
    attributes.setFlags(flags);
    attributes.setObsTimeRange(earliestObsTime, latestObsTime);
    attributes.setPosition(radians(ra_PS), radians(decl_PS),
                           radians(ra_PS_Sigma), radians(decl_PS_Sigma),
                           radians(radecl_PS_Cov));

    // per-filter columns
#define LSST_AP_HANDLE_PF_NULLS(i, filter) \
    do { \
        nullToNaN(db, i + 0, filter ## EarliestObsTime); \
        nullToNaN(db, i + 1, filter ## LatestObsTime); \
        if (db->columnIsNull(i + 2)) { \
            filter ## NumObs = 0; \
        } \
        if (db->columnIsNull(i + 3)) { \
            filter ## Flags = 0; \
        } \
        nullToNaN(db, i + 4, filter ## Flux_PS); \
        nullToNaN(db, i + 5, filter ## Flux_PS_Sigma); \
        nullToNaN(db, i + 6, filter ## E1_SG); \
        nullToNaN(db, i + 7, filter ## E1_SG_Sigma); \
        nullToNaN(db, i + 8, filter ## E2_SG); \
        nullToNaN(db, i + 9, filter ## E2_SG_Sigma); \
        nullToNaN(db, i + 10, filter ## Radius_SG); \
        nullToNaN(db, i + 11, filter ## Radius_SG_Sigma); \
    } while(false)

    LSST_AP_HANDLE_PF_NULLS(9 + 0*12, u);
    LSST_AP_HANDLE_PF_NULLS(9 + 1*12, g);
    LSST_AP_HANDLE_PF_NULLS(9 + 2*12, r);
    LSST_AP_HANDLE_PF_NULLS(9 + 3*12, i);
    LSST_AP_HANDLE_PF_NULLS(9 + 4*12, z);
    LSST_AP_HANDLE_PF_NULLS(9 + 5*12, y);

#undef LSST_AP_HANDLE_PF_NULLS

    attributes.setNumObs(uNumObs + gNumObs + rNumObs + iNumObs + zNumObs + yNumObs);

#define LSST_AP_TO_PFA(i, filter) \
    do { \
        if (filter ## NumObs > 0) { \
            PerFilterSourceClusterAttributes pfa; \
            pfa.setFilterId(filterIds[i]); \
            pfa.setNumObs(filter ## NumObs); \
            pfa.setFlags(filter ## Flags); \
            pfa.setObsTimeRange(filter ## EarliestObsTime, \
                                filter ## LatestObsTime); \
            pfa.setFlux(filter ## Flux_PS, filter ## Flux_PS_Sigma); \
            pfa.setEllipticity(filter ## E1_SG, \
                               filter ## E2_SG, \
                               filter ## Radius_SG, \
                               filter ## E1_SG_Sigma, \
                               filter ## E2_SG_Sigma, \
                               filter ## Radius_SG_Sigma); \
            attributes.setPerFilterAttributes(pfa); \
        } \
    } while (false)

    LSST_AP_TO_PFA(0, u);
    LSST_AP_TO_PFA(1, g);
    LSST_AP_TO_PFA(2, r);
    LSST_AP_TO_PFA(3, i);
    LSST_AP_TO_PFA(4, z);
    LSST_AP_TO_PFA(5, y);

#undef LSST_AP_TO_PFA
}

} // namespace


// -- SourceClusterVectorFormatter ----

SourceClusterVectorFormatter::SourceClusterVectorFormatter(Policy::Ptr policy) :
    Formatter(typeid(this)),
    _policy(policy)
{ }

SourceClusterVectorFormatter::~SourceClusterVectorFormatter() { }

Formatter::Ptr SourceClusterVectorFormatter::createInstance(Policy::Ptr policy) {
    return Formatter::Ptr(new SourceClusterVectorFormatter(policy));
}

FormatterRegistration SourceClusterVectorFormatter::registration(
    "PersistableSourceClusterVector",
    typeid(PersistableSourceClusterVector),
    createInstance);

template <class Archive>
void SourceClusterVectorFormatter::delegateSerialize(Archive & archive,
                                                     unsigned int const,
                                                     Persistable * persistable)
{
    PersistableSourceClusterVector * p =
        dynamic_cast<PersistableSourceClusterVector*>(persistable);
    archive & boost::serialization::base_object<Persistable>(*p);
    archive & p->_clusters;
}

/// @cond
template void SourceClusterVectorFormatter::delegateSerialize<boost::archive::text_oarchive>(
    boost::archive::text_oarchive &, unsigned int const, lsst::daf::base::Persistable *);
template void SourceClusterVectorFormatter::delegateSerialize<boost::archive::text_iarchive>(
    boost::archive::text_iarchive &, unsigned int const, lsst::daf::base::Persistable *);
/// @endcond

lsst::daf::base::Persistable * SourceClusterVectorFormatter::read(
    Storage::Ptr storage,
    PropertySet::Ptr additionalData)
{
    if (!storage) {
        throw LSST_EXCEPT(except::InvalidParameterException, "No Storage provided");
    }
    std::auto_ptr<PersistableSourceClusterVector> p(new PersistableSourceClusterVector);
    if (typeid(*storage) == typeid(BoostStorage)) {
        BoostStorage * b = dynamic_cast<BoostStorage *>(storage.get());
        b->getIArchive() & p->_clusters;
    } else if (typeid(*storage) == typeid(DbStorage) ||
               typeid(*storage) == typeid(DbTsvStorage)) {
        if (!additionalData) {
            throw LSST_EXCEPT(except::InvalidParameterException, "No PropertySet provided");
        }
        int filterIds[6];
        getFilterIds(filterIds);
        DbStorage * db = dynamic_cast<DbStorage *>(storage.get());
        std::string tableName = getTableName(_policy, additionalData);
        SourceClusterVector clusters;
        ObjectRow row;
        db->setTableForQuery(tableName);
        if (typeid(*storage) == typeid(DbStorage)) {
            row.setupFetch(db);
        } else {
            row.setupFetch(dynamic_cast<DbTsvStorage *>(storage.get()));
        }
        db->query();
        while (db->next()) {
            SourceClusterAttributes::Ptr sca(new SourceClusterAttributes);
            row.to(db, *sca, filterIds);
            clusters.push_back(sca);
        }
        p->setClusters(clusters);
    } else {
        throw LSST_EXCEPT(except::InvalidParameterException, "Unsupported Storage type");
    }
    return p.release();
}


void SourceClusterVectorFormatter::write(Persistable const * persistable,
                                         Storage::Ptr storage,
                                         PropertySet::Ptr additionalData)
{
    typedef SourceClusterVector::const_iterator Iter;

    if (!persistable) {
        throw LSST_EXCEPT(except::InvalidParameterException, "No Persistable provided");
    }
    if (!storage) {
        throw LSST_EXCEPT(except::InvalidParameterException, "No Storage provided");
    }
    PersistableSourceClusterVector const * p =
        dynamic_cast<PersistableSourceClusterVector const *>(persistable);
    if (p == 0) {
        throw LSST_EXCEPT(except::RuntimeErrorException, "Persistable was not "
                          "of concrete type PersistableSourceClusterVector");
    }

    if (typeid(*storage) == typeid(BoostStorage)) {
        BoostStorage * b = dynamic_cast<BoostStorage *>(storage.get());
        b->getOArchive() & p->_clusters;
    } else if (typeid(*storage) == typeid(DbStorage) || typeid(*storage) == typeid(DbTsvStorage)) {
        if (!additionalData) {
            throw LSST_EXCEPT(except::InvalidParameterException, "No PropertySet provided");
        }
        int filterIds[6];
        getFilterIds(filterIds);
        std::string tableName = getTableName(_policy, additionalData);
        std::string templateTableName = getTemplateTableName(_policy, additionalData);
        if (typeid(*storage) == typeid(DbStorage)) {
            DbStorage * db = dynamic_cast<DbStorage *>(storage.get());
            db->createTableFromTemplate(tableName, templateTableName, true);
            db->setTableForInsert(tableName);
            for (Iter i = p->_clusters.begin(), e = p->_clusters.end(); i != e; ++i) {
                insertObjectRow(*db, **i, filterIds);
            }
        } else {
            DbTsvStorage * db = dynamic_cast<DbTsvStorage *>(storage.get());
            db->createTableFromTemplate(tableName, templateTableName, true);
            db->setTableForInsert(tableName);
            for (Iter i = p->_clusters.begin(), e = p->_clusters.end(); i != e; ++i) {
                insertObjectRow(*db, **i, filterIds);
            }
        }
    } else {
        throw LSST_EXCEPT(except::InvalidParameterException, "Unsupported Storage type");
    }
}


void SourceClusterVectorFormatter::update(Persistable *, Storage::Ptr, PropertySet::Ptr)
{
    throw LSST_EXCEPT(except::RuntimeErrorException, "Updates not supported");
}

}}} // namespace lsst::ap::cluster
