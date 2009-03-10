// -*- lsst-c++ -*-

/**
 * @file
 * @brief   Persistence implementation for association pipeline result vectors.
 *
 * @ingroup associate
 */

#include <memory>

#include "boost/format.hpp"
//#include "boost/serialization/export.hpp"

#include "lsst/pex/exceptions.h"
#include "lsst/daf/base/Persistable.h"
#include "lsst/daf/base/PropertySet.h"
#include "lsst/daf/persistence/BoostStorage.h"
#include "lsst/daf/persistence/DbStorage.h"
#include "lsst/daf/persistence/DbTsvStorage.h"
#include "lsst/daf/persistence/FormatterImpl.h"

#include "lsst/afw/formatters/Utils.h"

#include "lsst/ap/Common.h"
#include "lsst/ap/Results.h"
#include "lsst/ap/Utils.h"
#include "lsst/ap/io/ResultFormatters.h"

//BOOST_CLASS_EXPORT(lsst::ap::MatchPairVector)
//BOOST_CLASS_EXPORT(lsst::ap::IdVector)

using lsst::daf::base::Persistable;
using lsst::daf::base::PropertySet;
using lsst::pex::policy::Policy;
using lsst::daf::persistence::BoostStorage;
using lsst::daf::persistence::DbStorage;
using lsst::daf::persistence::DbTsvStorage;
using lsst::daf::persistence::Formatter;
using lsst::daf::persistence::FormatterRegistration;
using lsst::daf::persistence::Storage;

namespace ex = lsst::pex::exceptions;
namespace fmt = lsst::afw::formatters;

// -- MatchPairVectorFormatter ----------------

lsst::ap::io::MatchPairVectorFormatter::MatchPairVectorFormatter(Policy::Ptr policy) :
    Formatter(typeid(*this)),
    _policy(policy)
{}


lsst::ap::io::MatchPairVectorFormatter::~MatchPairVectorFormatter() {}


FormatterRegistration lsst::ap::io::MatchPairVectorFormatter::registration(
    "PersistableMatchPairVector",
    typeid(PersistableMatchPairVector),
    createInstance
);


Formatter::Ptr lsst::ap::io::MatchPairVectorFormatter::createInstance(Policy::Ptr policy) {
    return Formatter::Ptr(new MatchPairVectorFormatter(policy));
}


template <class Archive>
void lsst::ap::io::MatchPairVectorFormatter::delegateSerialize(
    Archive &          archive,
    unsigned int const version,
    Persistable *      persistable
) {
    PersistableMatchPairVector * p = dynamic_cast<PersistableMatchPairVector *>(persistable);
    MatchPairVector::size_type sz;
    MatchPairVector & vec = p->getMatchPairs();

    if (Archive::is_loading::value) {
        MatchPair data;
        archive & sz;
        vec.reserve(sz);
        for (; sz > 0; --sz) {
            archive & data;
            vec.push_back(data);
        }
    } else {
        sz = vec.size();
        archive & sz;
        for (MatchPairVector::iterator i(vec.begin()), end(vec.end()); i != end; ++i) {
            archive & *i;
        }
    }
}

/// @cond
template void lsst::ap::io::MatchPairVectorFormatter::delegateSerialize<boost::archive::text_oarchive>(
    boost::archive::text_oarchive &, unsigned int const, Persistable *
);
template void lsst::ap::io::MatchPairVectorFormatter::delegateSerialize<boost::archive::text_iarchive>(
    boost::archive::text_iarchive &, unsigned int const, Persistable *
);
//template void lsst::ap::io::MatchPairVector::delegateSerialize<boost::archive::binary_oarchive>(
//    boost::archive::binary_oarchive &, unsigned int const, Persistable *
//);
//template void lsst::ap::io::MatchPairVector::delegateSerialize<boost::archive::binary_iarchive>(
//    boost::archive::binary_iarchive &, unsigned int const, Persistable *
//);
/// @endcond


void lsst::ap::io::MatchPairVectorFormatter::write(
    Persistable const * persistable,
    Storage::Ptr        storage,
    PropertySet::Ptr    additionalData
) {
    if (persistable == 0) {
        throw LSST_EXCEPT(ex::InvalidParameterException, "No Persistable provided");
    }
    if (!storage) {
        throw LSST_EXCEPT(ex::InvalidParameterException, "No Storage provided");
    }

    PersistableMatchPairVector const * p = dynamic_cast<PersistableMatchPairVector const *>(persistable);
    if (p == 0) {
        throw LSST_EXCEPT(ex::RuntimeErrorException,
                          "Persistable was not of concrete type PersistableMatchPairVector");
    }

    if (typeid(*storage) == typeid(BoostStorage)) {
        BoostStorage * bs = dynamic_cast<BoostStorage *>(storage.get());
        if (bs == 0) {
            throw LSST_EXCEPT(ex::RuntimeErrorException, "Didn't get BoostStorage");
        }
        bs->getOArchive() & *p;
    } else if (typeid(*storage) == typeid(DbStorage) || typeid(*storage) == typeid(DbTsvStorage)) {
        std::string name = fmt::getVisitSliceTableName(_policy, additionalData);
        std::string model = _policy->getString(_policy->getString("itemName") + ".templateTableName");
        MatchPairVector const & v = p->getMatchPairs();
        if (typeid(*storage) == typeid(DbStorage)) {
            DbStorage * db = dynamic_cast<DbStorage *>(storage.get());
            if (db == 0) {
                throw LSST_EXCEPT(ex::RuntimeErrorException, "Didn't get DbStorage");
            }
            db->createTableFromTemplate(name, model, false);
            db->setTableForInsert(name);
            for (MatchPairVector::const_iterator i(v.begin()), end(v.end()); i != end; ++i) {
                db->setColumn<boost::int64_t>("first", i->_first);
                db->setColumn<boost::int64_t>("second", i->_second);
                db->setColumn<double> ("distance", i->_distance);
                db->insertRow();
            }
        } else {
            DbTsvStorage * db = dynamic_cast<DbTsvStorage *>(storage.get());
            if (db == 0) {
                throw LSST_EXCEPT(ex::RuntimeErrorException, "Didn't get DbTsvStorage");
            }
            db->createTableFromTemplate(name, model, false);
            db->setTableForInsert(name);
            for (MatchPairVector::const_iterator i(v.begin()), end(v.end()); i != end; ++i) {
                db->setColumn<boost::int64_t>("first", i->_first);
                db->setColumn<boost::int64_t>("second", i->_second);
                db->setColumn<double> ("distance", i->_distance);
                db->insertRow();
            }
        }
    } else {
        throw LSST_EXCEPT(ex::InvalidParameterException, "Storage type is not supported");
    }
}


Persistable * lsst::ap::io::MatchPairVectorFormatter::read(
    Storage::Ptr     storage,
    PropertySet::Ptr additionalData
) {
    std::auto_ptr<PersistableMatchPairVector> p(new PersistableMatchPairVector);

    if (typeid(*storage) == typeid(BoostStorage)) {
        BoostStorage * bs = dynamic_cast<BoostStorage *>(storage.get());
        if (bs == 0) {
            throw LSST_EXCEPT(ex::RuntimeErrorException, "Didn't get BoostStorage");
        }
        bs->getIArchive() & *p;
    } else if (typeid(*storage) == typeid(DbStorage) || typeid(*storage) == typeid(DbTsvStorage)) {
        DbStorage * db = dynamic_cast<DbStorage *>(storage.get());
        if (db == 0) {
            throw LSST_EXCEPT(ex::RuntimeErrorException, "Didn't get DbStorage");
        }
        db->setTableForQuery(fmt::getVisitSliceTableName(_policy, additionalData));
        MatchPair data;
        db->outParam("first",    &(data._first));
        db->outParam("second",   &(data._second));
        db->outParam("distance", &(data._distance));
        db->query();
        while (db->next()) {
            if (db->columnIsNull(0)) {
                throw LSST_EXCEPT(ex::RuntimeErrorException, "null column \"first\"");
            }
            if (db->columnIsNull(1)) {
                throw LSST_EXCEPT(ex::RuntimeErrorException, "null column \"second\"");
            }
            if (db->columnIsNull(2)) {
                throw LSST_EXCEPT(ex::RuntimeErrorException,  "null column \"distance\"");
            }
            p->getMatchPairs().push_back(data);
        }
        db->finishQuery();
    } else {
        throw LSST_EXCEPT(ex::InvalidParameterException, "Storage type is not supported");
    }
    return p.release();
}


void lsst::ap::io::MatchPairVectorFormatter::update(Persistable *, Storage::Ptr, PropertySet::Ptr) {
    throw LSST_EXCEPT(ex::RuntimeErrorException, "MatchPairVectorFormatter: updates not supported");
}


// -- IdPairVectorFormatter ----------------

lsst::ap::io::IdPairVectorFormatter::IdPairVectorFormatter(Policy::Ptr policy) :
    Formatter(typeid(*this)),
    _policy(policy)
{}


lsst::ap::io::IdPairVectorFormatter::~IdPairVectorFormatter() {}


FormatterRegistration lsst::ap::io::IdPairVectorFormatter::registration(
    "IdPairVector",
    typeid(IdPairVector),
    createInstance
);


Formatter::Ptr lsst::ap::io::IdPairVectorFormatter::createInstance(Policy::Ptr policy) {
    return Formatter::Ptr(new IdPairVectorFormatter(policy));
}


template <class Archive>
void lsst::ap::io::IdPairVectorFormatter::delegateSerialize(
    Archive &          archive,
    unsigned int const version,
    Persistable *      persistable
) {
    PersistableIdPairVector * p = dynamic_cast<PersistableIdPairVector *>(persistable);
    IdPairVector::size_type sz;
    IdPairVector & vec = p->getIdPairs();

    if (Archive::is_loading::value) {
        IdPair data;
        archive & sz;
        vec.reserve(sz);
        for (; sz > 0; --sz) {
            archive & data.first;
            archive & data.second;
            vec.push_back(data);
        }
    } else {
        sz = vec.size();
        archive & sz;
        for (IdPairVector::iterator i(vec.begin()), end(vec.end()); i != end; ++i) {
            archive & i->first;
            archive & i->second;
        }
    }
}

/// @cond
template void lsst::ap::io::IdPairVectorFormatter::delegateSerialize<boost::archive::text_oarchive>(
    boost::archive::text_oarchive &, unsigned int const, Persistable *
);
template void lsst::ap::io::IdPairVectorFormatter::delegateSerialize<boost::archive::text_iarchive>(
    boost::archive::text_iarchive &, unsigned int const, Persistable *
);
//template void lsst::ap::io::IdPairVectorFormatter::delegateSerialize<boost::archive::binary_oarchive>(
//    boost::archive::binary_oarchive &, unsigned int const, Persistable *
//);
//template void lsst::ap::io::IdPairVectorFormatter::delegateSerialize<boost::archive::binary_iarchive>(
//    boost::archive::binary_iarchive &, unsigned int const, Persistable *
//);
/// @endcond


void lsst::ap::io::IdPairVectorFormatter::write(
    Persistable const * persistable,
    Storage::Ptr        storage,
    PropertySet::Ptr    additionalData
) {
    if (persistable == 0) {
        throw LSST_EXCEPT(ex::InvalidParameterException, "No Persistable provided");
    }
    if (!storage) {
        throw LSST_EXCEPT(ex::InvalidParameterException, "No Storage provided");
    }

    PersistableIdPairVector const * p = dynamic_cast<PersistableIdPairVector const *>(persistable);
    if (p == 0) {
        throw LSST_EXCEPT(ex::RuntimeErrorException, "Persistable was not of concrete type IdPairVector");
    }

    if (typeid(*storage) == typeid(BoostStorage)) {
        BoostStorage * bs = dynamic_cast<BoostStorage *>(storage.get());
        if (bs == 0) {
            throw LSST_EXCEPT(ex::RuntimeErrorException, "Didn't get BoostStorage");
        }
        bs->getOArchive() & *p;
    } else if (typeid(*storage) == typeid(DbStorage) || typeid(*storage) == typeid(DbTsvStorage)) {
        std::string name = fmt::getVisitSliceTableName(_policy, additionalData);
        std::string model = _policy->getString(_policy->getString("itemName") + ".templateTableName");
        IdPairVector const & v = p->getIdPairs();
        if (typeid(*storage) == typeid(DbStorage)) { 
            DbStorage * db = dynamic_cast<DbStorage *>(storage.get());
            if (db == 0) {
                throw LSST_EXCEPT(ex::RuntimeErrorException, "Didn't get DbStorage");
            }
            db->createTableFromTemplate(name, model, false);
            db->setTableForInsert(name);
            for (IdPairVector::const_iterator i(v.begin()), end(v.end()); i != end; ++i) {
                db->setColumn<boost::int64_t>("first",  i->first);
                db->setColumn<boost::int64_t>("second", i->second);
                db->insertRow();
            }
        } else {
            DbTsvStorage * db = dynamic_cast<DbTsvStorage *>(storage.get());
            if (db == 0) {
                throw LSST_EXCEPT(ex::RuntimeErrorException, "Didn't get DbTsvStorage");
            }
            db->createTableFromTemplate(name, model, false);
            db->setTableForInsert(name);
            for (IdPairVector::const_iterator i(v.begin()), end(v.end()); i != end; ++i) {
                db->setColumn<boost::int64_t>("first",  i->first);
                db->setColumn<boost::int64_t>("second", i->second);
                db->insertRow();
            }
        }
    } else {
        throw LSST_EXCEPT(ex::InvalidParameterException, "Storage type is not supported");
    }
}


Persistable * lsst::ap::io::IdPairVectorFormatter::read(
    Storage::Ptr     storage,
    PropertySet::Ptr additionalData
) {
    std::auto_ptr<PersistableIdPairVector> p(new PersistableIdPairVector);

    if (typeid(*storage) == typeid(BoostStorage)) {
        BoostStorage* bs = dynamic_cast<BoostStorage *>(storage.get());
        if (bs == 0) {
            throw LSST_EXCEPT(ex::RuntimeErrorException, "Didn't get BoostStorage");
        }
        bs->getIArchive() & *p;
    } else if (typeid(*storage) == typeid(DbStorage) || typeid(*storage) == typeid(DbTsvStorage)) {
        DbStorage * db = dynamic_cast<DbStorage *>(storage.get());
        if (db == 0) {
            throw LSST_EXCEPT(ex::RuntimeErrorException, "Didn't get DbStorage");
        }

        db->setTableForQuery(fmt::getVisitSliceTableName(_policy, additionalData));
        IdPair data;
        db->outParam("first",  &data.first);
        db->outParam("second", &data.second);
        db->query();
        while (db->next()) {
            if (db->columnIsNull(0)) {
                throw LSST_EXCEPT(ex::RuntimeErrorException, "null column \"first\"");
            }
            if (db->columnIsNull(1)) {
                throw LSST_EXCEPT(ex::RuntimeErrorException, "null column \"second\"");
            }
            p->getIdPairs().push_back(data);
        }
        db->finishQuery();
    } else {
        throw LSST_EXCEPT(ex::InvalidParameterException, "Storage type is not supported");
    }
    return p.release();
}


void lsst::ap::io::IdPairVectorFormatter::update(Persistable *, Storage::Ptr, PropertySet::Ptr) {
    throw LSST_EXCEPT(ex::RuntimeErrorException, "IdPairVectorFormatter: updates not supported");
}


// -- IdVectorFormatter ----------------

lsst::ap::io::IdVectorFormatter::IdVectorFormatter(Policy::Ptr policy) :
    Formatter(typeid(*this)),
    _policy(policy)
{}


lsst::ap::io::IdVectorFormatter::~IdVectorFormatter() {}


FormatterRegistration lsst::ap::io::IdVectorFormatter::registration(
    "IdVector",
    typeid(IdVector),
    createInstance
);


Formatter::Ptr lsst::ap::io::IdVectorFormatter::createInstance(Policy::Ptr policy) {
    return Formatter::Ptr(new IdVectorFormatter(policy));
}


template <class Archive>
void lsst::ap::io::IdVectorFormatter::delegateSerialize(
    Archive &          archive,
    unsigned int const version,
    Persistable *      persistable
) {
    PersistableIdVector * p = dynamic_cast<PersistableIdVector *>(persistable);
    IdVector::size_type sz;
    IdVector & vec = p->getIds();

    if (Archive::is_loading::value) {
        boost::int64_t data;
        archive & sz;
        vec.reserve(sz);
        for (; sz > 0; --sz) {
            archive & data;
            vec.push_back(data);
        }
    } else {
        sz = vec.size();
        archive & sz;
        for (IdVector::iterator i(vec.begin()), end(vec.end()); i != end; ++i) {
            archive & *i;
        }
    }
}

/// @cond
template void lsst::ap::io::IdVectorFormatter::delegateSerialize<boost::archive::text_oarchive>(
    boost::archive::text_oarchive &, unsigned int const, Persistable *
);
template void lsst::ap::io::IdVectorFormatter::delegateSerialize<boost::archive::text_iarchive>(
    boost::archive::text_iarchive &, unsigned int const, Persistable *
);
//template void lsst::ap::io::IdVectorFormatter::delegateSerialize<boost::archive::binary_oarchive>(
//    boost::archive::binary_oarchive &, unsigned int const, Persistable *
//);
//template void lsst::ap::io::IdVectorFormatter::delegateSerialize<boost::archive::binary_iarchive>(
//    boost::archive::binary_iarchive &, unsigned int const, Persistable *
//);
/// @endcond


void lsst::ap::io::IdVectorFormatter::write(
    Persistable const * persistable,
    Storage::Ptr        storage,
    PropertySet::Ptr    additionalData
) {
    if (persistable == 0) {
        throw LSST_EXCEPT(ex::InvalidParameterException, "No Persistable provided");
    }
    if (!storage) {
        throw LSST_EXCEPT(ex::InvalidParameterException, "No Storage provided");
    }

    PersistableIdVector const * p = dynamic_cast<PersistableIdVector const *>(persistable);
    if (p == 0) {
        throw LSST_EXCEPT(ex::RuntimeErrorException, "Persistable was not of concrete type IdVector");
    }

    if (typeid(*storage) == typeid(BoostStorage)) {
        BoostStorage * bs = dynamic_cast<BoostStorage *>(storage.get());
        if (bs == 0) {
            throw LSST_EXCEPT(ex::RuntimeErrorException, "Didn't get BoostStorage");
        }
        bs->getOArchive() & *p;
    } else if (typeid(*storage) == typeid(DbStorage) || typeid(*storage) == typeid(DbTsvStorage)) {
        std::string name = fmt::getVisitSliceTableName(_policy, additionalData);
        std::string model = _policy->getString(_policy->getString("itemName") + ".templateTableName");
        IdVector const & v = p->getIds();
        if (typeid(*storage) == typeid(DbStorage)) {
            DbStorage * db = dynamic_cast<DbStorage *>(storage.get());
            if (db == 0) {
                throw LSST_EXCEPT(ex::RuntimeErrorException, "Didn't get DbStorage");
            }
            db->createTableFromTemplate(name, model, false);
            db->setTableForInsert(name);
            for (IdVector::const_iterator i(v.begin()), end(v.end()); i != end; ++i) {
                db->setColumn<boost::int64_t>("id", *i);
                db->insertRow();
            }
        } else {
            DbTsvStorage * db = dynamic_cast<DbTsvStorage *>(storage.get());
            if (db == 0) {
                throw LSST_EXCEPT(ex::RuntimeErrorException, "Didn't get DbTsvStorage");
            }
            db->createTableFromTemplate(name, model, false);
            db->setTableForInsert(name);
            for (IdVector::const_iterator i(v.begin()), end(v.end()); i != end; ++i) {
                db->setColumn<long long>("id", *i);
                db->insertRow();
            }
        }
    } else {
        throw LSST_EXCEPT(ex::InvalidParameterException, "Storage type is not supported");
    }
}


Persistable * lsst::ap::io::IdVectorFormatter::read(
    Storage::Ptr     storage,
    PropertySet::Ptr additionalData
) {
    std::auto_ptr<PersistableIdVector> p(new PersistableIdVector);

    if (typeid(*storage) == typeid(BoostStorage)) {
        BoostStorage* bs = dynamic_cast<BoostStorage *>(storage.get());
        if (bs == 0) {
            throw LSST_EXCEPT(ex::RuntimeErrorException, "Didn't get BoostStorage");
        }
        bs->getIArchive() & *p;
    } else if (typeid(*storage) == typeid(DbStorage) || typeid(*storage) == typeid(DbTsvStorage)) {
        DbStorage * db = dynamic_cast<DbStorage *>(storage.get());
        if (db == 0) {
            throw LSST_EXCEPT(ex::RuntimeErrorException, "Didn't get DbStorage");
        }

        db->setTableForQuery(fmt::getVisitSliceTableName(_policy, additionalData));
        boost::int64_t data;
        db->outParam("id", &data);
        db->query();
        while (db->next()) {
            if (db->columnIsNull(0)) {
                throw LSST_EXCEPT(ex::RuntimeErrorException, "null column \"id\"");
            }
            p->getIds().push_back(data);
        }
        db->finishQuery();
    } else {
        throw LSST_EXCEPT(ex::InvalidParameterException, "Storage type is not supported");
    }
    return p.release();
}


void lsst::ap::io::IdVectorFormatter::update(Persistable *, Storage::Ptr, PropertySet::Ptr) {
    throw LSST_EXCEPT(ex::RuntimeErrorException, "IdVectorFormatter: updates not supported");
}

