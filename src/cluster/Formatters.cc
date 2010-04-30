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

// -- serialize() functions required by boost::serialization ----

template <typename Archive, typename FloatT>
inline void serializeFloat(Archive & ar, FloatT value) {
    bool null = isnan(value);
    ar & null;
    if (null) {
        ar & value;
    }
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

// re-use STL machinery in boost::serialization to persist std::tr1::unordered_map
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
void serialize(
    Archive & ar,
    std::tr1::unordered_map<int, PerFilterSourceClusterAttributes> & map,
    unsigned int const version)
{
    boost::serialization::split_free(ar, map, version);
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


// -- Serialization helper functions and classes ----

namespace {

/** @internal
  * Returns the name of the database table to read/write to.
  */
std::string const getTableName(Policy::Ptr policy,
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
std::string const getTemplateTableName(Policy::Ptr policy,
                                       PropertySet::Ptr additionalData)
{
    std::string itemName = additionalData->getAsString("itemName");
    return policy->getString(itemName + ".templateTableName");
}


template <typename StorageT, typename T>
inline void insertFloat(StorageT & db, char const * const col, T const & val) {
    if (isnan(val)) {
        db.setColumnToNull(col);
    } else {
        db.template setColumn<T>(col, val);
    }
}

template <typename FloatT> FloatT radians(FloatT deg) {
    return static_cast<FloatT>(deg * (M_PI/180.0));
}

template <typename FloatT> FloatT degrees(FloatT rad) {
    return static_cast<FloatT>(rad * (180.0/M_PI));
}

template <typename FloatT> FloatT rangeReducedDegrees(FloatT rad) {
    FloatT deg = std::fmod(degrees(rad), static_cast<FloatT>(360));
    return deg < static_cast<FloatT>(0) ? deg + static_cast<FloatT>(360) : deg;
}

template <typename StorageT>
void insertRow(StorageT & db,
               SourceClusterAttributes const & attributes,
               int filterIds[6])
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

#define LSST_AP_ALWAYS_NULL(filter) \
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

    LSST_AP_ALWAYS_NULL(u);
    LSST_AP_ALWAYS_NULL(g);
    LSST_AP_ALWAYS_NULL(r);
    LSST_AP_ALWAYS_NULL(i);
    LSST_AP_ALWAYS_NULL(z);
    LSST_AP_ALWAYS_NULL(y);

#undef LSST_AP_ALWAYS_NULL

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
#define LSST_AP_SET_FILTER(filter, index) \
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

    LSST_AP_SET_FILTER(u, 0);
    LSST_AP_SET_FILTER(g, 1);
    LSST_AP_SET_FILTER(r, 2);
    LSST_AP_SET_FILTER(i, 3);
    LSST_AP_SET_FILTER(z, 4);
    LSST_AP_SET_FILTER(y, 5);

#undef LSST_AP_SET_FILTER
}

/** @internal
  * A simple C++ realization of a row in the Object database table.
  */
struct Object {

    Object() { }

    Object(SourceClusterAttributes const & attributes) {
    }

    void fill(SourceClusterAttributes & attributes) {
    }
};

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
    } else if (typeid(*storage) == typeid(DbStorage) || typeid(*storage) == typeid(DbTsvStorage)) {
        if (!additionalData) {
            throw LSST_EXCEPT(except::InvalidParameterException, "No PropertySet provided");
        }
        DbStorage * db = dynamic_cast<DbStorage *>(storage.get());
        std::string tableName = getTableName(_policy, additionalData);
         
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
        int filterIds[6] = { Filter("u").getId(),
                             Filter("g").getId(),
                             Filter("r").getId(),
                             Filter("i").getId(),
                             Filter("z").getId(),
                             Filter("y").getId() };
        if (filterIds[5] == Filter::UNKNOWN) {
            filterIds[5] = Filter("i2").getId();
        }
        std::string tableName = getTableName(_policy, additionalData);
        std::string templateTableName = getTemplateTableName(_policy, additionalData);
        if (typeid(*storage) == typeid(DbStorage)) {
            DbStorage * db = dynamic_cast<DbStorage *>(storage.get());
            db->createTableFromTemplate(tableName, templateTableName, true);
            db->setTableForInsert(tableName);
            for (Iter i = p->_clusters.begin(), e = p->_clusters.end(); i != e; ++i) {
                insertRow(*db, **i, filterIds);
            }
        } else {
            DbTsvStorage * db = dynamic_cast<DbTsvStorage *>(storage.get());
            db->createTableFromTemplate(tableName, templateTableName, true);
            db->setTableForInsert(tableName);
            for (Iter i = p->_clusters.begin(), e = p->_clusters.end(); i != e; ++i) {
                insertRow(*db, **i, filterIds);
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
