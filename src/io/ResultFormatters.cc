// -*- lsst-c++ -*-
//
//##====----------------                                ----------------====##/
//
//! \file   ResultFormatters.cc
//! \brief  Implementation of persistence for association pipeline results.
//
//##====----------------                                ----------------====##/

#include <memory>

#include <lsst/mwi/persistence/BoostStorage.h>
#include <lsst/mwi/persistence/DbStorage.h>
#include <lsst/mwi/persistence/DbTsvStorage.h>
#include <lsst/mwi/persistence/FormatterImpl.h>

#include <boost/any.hpp>
#include <boost/format.hpp>
//#include <boost/serialization/export.hpp>

#include <lsst/ap/Common.h>
#include <lsst/ap/Exceptions.h>
#include <lsst/ap/Results.h>
#include <lsst/ap/Utils.h>
#include <lsst/ap/io/ResultFormatters.h>

//BOOST_CLASS_EXPORT(lsst::ap::MatchPairVector)
//BOOST_CLASS_EXPORT(lsst::ap::IdVector)


namespace lsst {
namespace ap {
namespace io {

using lsst::mwi::persistence::BoostStorage;
using lsst::mwi::persistence::DbTsvStorage;

// -- MatchPairVectorFormatter ----------------

MatchPairVectorFormatter::MatchPairVectorFormatter(Policy::Ptr const & policy) :
    Formatter(typeid(*this)),
    _policy(policy)
{}


MatchPairVectorFormatter::~MatchPairVectorFormatter() {}


FormatterRegistration MatchPairVectorFormatter::registration(
    "MatchPairVector",
    typeid(MatchPairVector),
    createInstance
);


Formatter::Ptr MatchPairVectorFormatter::createInstance(Policy::Ptr policy) {
    return Formatter::Ptr(new MatchPairVectorFormatter(policy));
}


template <class Archive>
void MatchPairVectorFormatter::delegateSerialize(
    Archive &          archive,
    unsigned int const version,
    Persistable *      persistable
) {
    MatchPairVector * p = dynamic_cast<MatchPairVector *>(persistable);
    archive & boost::serialization::base_object<Persistable>(*p);
    MatchPairVector::size_type sz;

    if (Archive::is_loading::value) {
        MatchPair data;
        archive & sz;
        p->reserve(sz);
        for (; sz > 0; --sz) {
            archive & data;
            p->push_back(data);
        }
    } else {
        sz = p->size();
        archive & sz;
        MatchPairVector::iterator const end(p->end());
        for (MatchPairVector::iterator i = p->begin(); i != end; ++i) {
            archive & *i;
        }
    }
}

//! \cond
template void MatchPairVectorFormatter::delegateSerialize<boost::archive::text_oarchive>(
    boost::archive::text_oarchive &, unsigned int const, Persistable *
);
template void MatchPairVectorFormatter::delegateSerialize<boost::archive::text_iarchive>(
    boost::archive::text_iarchive &, unsigned int const, Persistable *
);
//template void MatchPairVector::delegateSerialize<boost::archive::binary_oarchive>(
//    boost::archive::binary_oarchive &, unsigned int const, Persistable *
//);
//template void MatchPairVector::delegateSerialize<boost::archive::binary_iarchive>(
//    boost::archive::binary_iarchive &, unsigned int const, Persistable *
//);
//! \endcond


void MatchPairVectorFormatter::write(
    Persistable const *   persistable,
    Storage::Ptr          storage,
    DataProperty::PtrType additionalData
) {
    if (persistable == 0) {
        LSST_AP_THROW(InvalidParameter, "No Persistable provided");
    }
    if (!storage) {
        LSST_AP_THROW(InvalidParameter, "No Storage provided");
    }

    MatchPairVector const * p = dynamic_cast<MatchPairVector const *>(persistable);
    if (p == 0) {
        LSST_AP_THROW(Runtime, "Persistable was not of concrete type MatchPairVector");
    }

    if (typeid(*storage) == typeid(BoostStorage)) {
        BoostStorage * bs = dynamic_cast<BoostStorage *>(storage.get());
        if (bs == 0) {
            LSST_AP_THROW(Runtime, "Didn't get BoostStorage");
        }
        bs->getOArchive() & *p;
    } else if (typeid(*storage) == typeid(DbStorage) || typeid(*storage) == typeid(DbTsvStorage)) {
        std::string name(getTableName(_policy, additionalData));
        std::string model(getTableTemplateName(_policy, additionalData));
        if (typeid(*storage) == typeid(DbStorage)) {
            DbStorage * db = dynamic_cast<DbStorage *>(storage.get());
            if (db == 0) {
                LSST_AP_THROW(Runtime, "Didn't get DbStorage");
            }
            db->createTableFromTemplate(name, model, false);
            db->setTableForInsert(name);
            MatchPairVector::const_iterator const end(p->end());
            for (MatchPairVector::const_iterator i = p->begin(); i != end; ++i) {
                db->setColumn<int64_t>("first",    i->_first);
                db->setColumn<int64_t>("second",   i->_second);
                db->setColumn<double> ("distance", i->_distance);
                db->insertRow();
            }
        } else {
            DbTsvStorage * db = dynamic_cast<DbTsvStorage *>(storage.get());
            if (db == 0) {
                LSST_AP_THROW(Runtime, "Didn't get DbTsvStorage");
            }
            db->createTableFromTemplate(name, model, false);
            db->setTableForInsert(name);
            MatchPairVector::const_iterator const end(p->end());
            for (MatchPairVector::const_iterator i = p->begin(); i != end; ++i) {
                db->setColumn<int64_t>("first",    i->_first);
                db->setColumn<int64_t>("second",   i->_second);
                db->setColumn<double> ("distance", i->_distance);
                db->insertRow();
            }
        }
    } else {
        LSST_AP_THROW(InvalidParameter, "Storage type is not supported");
    }
}


Persistable * MatchPairVectorFormatter::read(
    Storage::Ptr          storage,
    DataProperty::PtrType additionalData
) {
    std::auto_ptr<MatchPairVector> p(new MatchPairVector);

    if (typeid(*storage) == typeid(BoostStorage)) {
        BoostStorage * bs = dynamic_cast<BoostStorage *>(storage.get());
        if (bs == 0) {
            LSST_AP_THROW(Runtime, "Didn't get BoostStorage");
        }
        bs->getIArchive() & *p;
    } else if (typeid(*storage) == typeid(DbStorage) || typeid(*storage) == typeid(DbTsvStorage)) {
        DbStorage * db = dynamic_cast<DbStorage *>(storage.get());
        if (db == 0) {
            LSST_AP_THROW(Runtime, "Didn't get DbStorage");
        }
        db->setTableForQuery(getTableName(_policy, additionalData));
        MatchPair data;
        db->outParam("first",    &(data._first));
        db->outParam("second",   &(data._second));
        db->outParam("distance", &(data._distance));
        db->query();
        while (db->next()) {
            if (db->columnIsNull(0)) { LSST_AP_THROW(Runtime, "null column \"first\"");    }
            if (db->columnIsNull(1)) { LSST_AP_THROW(Runtime, "null column \"second\"");   }
            if (db->columnIsNull(2)) { LSST_AP_THROW(Runtime, "null column \"distance\""); }
            p->push_back(data);
        }
        db->finishQuery();
    } else {
        LSST_AP_THROW(InvalidParameter, "Storage type is not supported");
    }
    return p.release();
}


void MatchPairVectorFormatter::update(Persistable *, Storage::Ptr, DataProperty::PtrType) {
    LSST_AP_THROW(Runtime, "MatchPairVectorFormatter: updates not supported");
}


// -- IdPairVectorFormatter ----------------

IdPairVectorFormatter::IdPairVectorFormatter(Policy::Ptr const & policy) :
    Formatter(typeid(*this)),
    _policy(policy)
{}


IdPairVectorFormatter::~IdPairVectorFormatter() {}


FormatterRegistration IdPairVectorFormatter::registration("IdPairVector", typeid(IdPairVector), createInstance);


Formatter::Ptr IdPairVectorFormatter::createInstance(Policy::Ptr policy) {
    return Formatter::Ptr(new IdPairVectorFormatter(policy));
}


template <class Archive>
void IdPairVectorFormatter::delegateSerialize(
    Archive &          archive,
    unsigned int const version,
    Persistable *      persistable
) {
    IdPairVector * p = dynamic_cast<IdPairVector *>(persistable);
    archive & boost::serialization::base_object<Persistable>(*p);
    IdPairVector::size_type sz;

    if (Archive::is_loading::value) {
        IdPair data;
        archive & sz;
        p->reserve(sz);
        for (; sz > 0; --sz) {
            archive & data.first;
            archive & data.second;
            p->push_back(data);
        }
    } else {
        sz = p->size();
        archive & sz;
        IdPairVector::iterator const end(p->end());
        for (IdPairVector::iterator i = p->begin(); i != end; ++i) {
            archive & i->first;
            archive & i->second;
        }
    }
}

//! \cond
template void IdPairVectorFormatter::delegateSerialize<boost::archive::text_oarchive>(
    boost::archive::text_oarchive &, unsigned int const, Persistable *
);
template void IdPairVectorFormatter::delegateSerialize<boost::archive::text_iarchive>(
    boost::archive::text_iarchive &, unsigned int const, Persistable *
);
//template void IdPairVectorFormatter::delegateSerialize<boost::archive::binary_oarchive>(
//    boost::archive::binary_oarchive &, unsigned int const, Persistable *
//);
//template void IdPairVectorFormatter::delegateSerialize<boost::archive::binary_iarchive>(
//    boost::archive::binary_iarchive &, unsigned int const, Persistable *
//);
//! \endcond


void IdPairVectorFormatter::write(
    Persistable const *   persistable,
    Storage::Ptr          storage,
    DataProperty::PtrType additionalData
) {
    if (persistable == 0) {
        LSST_AP_THROW(InvalidParameter, "No Persistable provided");
    }
    if (!storage) {
        LSST_AP_THROW(InvalidParameter, "No Storage provided");
    }

    IdPairVector const * p = dynamic_cast<IdPairVector const *>(persistable);
    if (p == 0) {
        LSST_AP_THROW(Runtime, "Persistable was not of concrete type IdPairVector");
    }

    if (typeid(*storage) == typeid(BoostStorage)) {
        BoostStorage * bs = dynamic_cast<BoostStorage *>(storage.get());
        if (bs == 0) {
            LSST_AP_THROW(Runtime, "Didn't get BoostStorage");
        }
        bs->getOArchive() & *p;
    } else if (typeid(*storage) == typeid(DbStorage) || typeid(*storage) == typeid(DbTsvStorage)) {
        std::string name(getTableName(_policy, additionalData));
        std::string model(getTableTemplateName(_policy, additionalData));
        if (typeid(*storage) == typeid(DbStorage)) {
            DbStorage * db = dynamic_cast<DbStorage *>(storage.get());
            if (db == 0) {
                LSST_AP_THROW(Runtime, "Didn't get DbStorage");
            }
            db->createTableFromTemplate(name, model, false);
            db->setTableForInsert(name);
            IdPairVector::const_iterator const end(p->end());
            for (IdPairVector::const_iterator i = p->begin(); i != end; ++i) {
                db->setColumn<int64_t>("first",  i->first);
                db->setColumn<int64_t>("second", i->second);
                db->insertRow();
            }
        } else {
            DbTsvStorage * db = dynamic_cast<DbTsvStorage *>(storage.get());
            if (db == 0) {
                LSST_AP_THROW(Runtime, "Didn't get DbTsvStorage");
            }
            db->createTableFromTemplate(name, model, false);
            db->setTableForInsert(name);
            IdPairVector::const_iterator const end(p->end());
            for (IdPairVector::const_iterator i = p->begin(); i != end; ++i) {
                db->setColumn<int64_t>("first",  i->first);
                db->setColumn<int64_t>("second", i->second);
                db->insertRow();
            }
        }
    } else {
        LSST_AP_THROW(InvalidParameter, "Storage type is not supported");
    }
}


Persistable * IdPairVectorFormatter::read(
    Storage::Ptr          storage,
    DataProperty::PtrType additionalData
) {
    std::auto_ptr<IdPairVector> p(new IdPairVector);

    if (typeid(*storage) == typeid(BoostStorage)) {
        BoostStorage* bs = dynamic_cast<BoostStorage *>(storage.get());
        if (bs == 0) {
            LSST_AP_THROW(Runtime, "Didn't get BoostStorage");
        }
        bs->getIArchive() & *p;
    } else if (typeid(*storage) == typeid(DbStorage) || typeid(*storage) == typeid(DbTsvStorage)) {
        DbStorage * db = dynamic_cast<DbStorage *>(storage.get());
        if (db == 0) {
            LSST_AP_THROW(Runtime, "Didn't get DbStorage");
        }

        db->setTableForQuery(getTableName(_policy, additionalData));
        IdPair data;
        db->outParam("first",  &data.first);
        db->outParam("second", &data.second);
        db->query();
        while (db->next()) {
            if (db->columnIsNull(0)) {
                LSST_AP_THROW(Runtime, "null column \"first\"");
            }
            if (db->columnIsNull(1)) {
                LSST_AP_THROW(Runtime, "null column \"second\"");
            }
            p->push_back(data);
        }
        db->finishQuery();
    } else {
        LSST_AP_THROW(InvalidParameter, "Storage type is not supported");
    }
    return p.release();
}


void IdPairVectorFormatter::update(Persistable *, Storage::Ptr, DataProperty::PtrType) {
    LSST_AP_THROW(Runtime, "IdPairVectorFormatter: updates not supported");
}


// -- IdVectorFormatter ----------------

IdVectorFormatter::IdVectorFormatter(Policy::Ptr const & policy) :
    Formatter(typeid(*this)),
    _policy(policy)
{}


IdVectorFormatter::~IdVectorFormatter() {}


FormatterRegistration IdVectorFormatter::registration("IdVector", typeid(IdVector), createInstance);


Formatter::Ptr IdVectorFormatter::createInstance(Policy::Ptr policy) {
    return Formatter::Ptr(new IdVectorFormatter(policy));
}


template <class Archive>
void IdVectorFormatter::delegateSerialize(
    Archive &          archive,
    unsigned int const version,
    Persistable *      persistable
) {
    IdVector * p = dynamic_cast<IdVector *>(persistable);
    archive & boost::serialization::base_object<Persistable>(*p);
    IdVector::size_type sz;

    if (Archive::is_loading::value) {
        int64_t data;
        archive & sz;
        p->reserve(sz);
        for (; sz > 0; --sz) {
            archive & data;
            p->push_back(data);
        }
    } else {
        sz = p->size();
        archive & sz;
        IdVector::iterator const end(p->end());
        for (IdVector::iterator i = p->begin(); i != end; ++i) {
            archive & *i;
        }
    }
}

//! \cond
template void IdVectorFormatter::delegateSerialize<boost::archive::text_oarchive>(
    boost::archive::text_oarchive &, unsigned int const, Persistable *
);
template void IdVectorFormatter::delegateSerialize<boost::archive::text_iarchive>(
    boost::archive::text_iarchive &, unsigned int const, Persistable *
);
//template void IdVectorFormatter::delegateSerialize<boost::archive::binary_oarchive>(
//    boost::archive::binary_oarchive &, unsigned int const, Persistable *
//);
//template void IdVectorFormatter::delegateSerialize<boost::archive::binary_iarchive>(
//    boost::archive::binary_iarchive &, unsigned int const, Persistable *
//);
//! \endcond


void IdVectorFormatter::write(
    Persistable const *   persistable,
    Storage::Ptr          storage,
    DataProperty::PtrType additionalData
) {
    if (persistable == 0) {
        LSST_AP_THROW(InvalidParameter, "No Persistable provided");
    }
    if (!storage) {
        LSST_AP_THROW(InvalidParameter, "No Storage provided");
    }

    IdVector const * p = dynamic_cast<IdVector const *>(persistable);
    if (p == 0) {
        LSST_AP_THROW(Runtime, "Persistable was not of concrete type IdVector");
    }

    if (typeid(*storage) == typeid(BoostStorage)) {
        BoostStorage * bs = dynamic_cast<BoostStorage *>(storage.get());
        if (bs == 0) {
            LSST_AP_THROW(Runtime, "Didn't get BoostStorage");
        }
        bs->getOArchive() & *p;
    } else if (typeid(*storage) == typeid(DbStorage) || typeid(*storage) == typeid(DbTsvStorage)) {
        std::string name(getTableName(_policy, additionalData));
        std::string model(getTableTemplateName(_policy, additionalData));
        if (typeid(*storage) == typeid(DbStorage)) {
            DbStorage * db = dynamic_cast<DbStorage *>(storage.get());
            if (db == 0) {
                LSST_AP_THROW(Runtime, "Didn't get DbStorage");
            }
            db->createTableFromTemplate(name, model, false);
            db->setTableForInsert(name);
            IdVector::const_iterator const end(p->end());
            for (IdVector::const_iterator i = p->begin(); i != end; ++i) {
                db->setColumn<int64_t>("id", *i);
                db->insertRow();
            }
        } else {
            DbTsvStorage * db = dynamic_cast<DbTsvStorage *>(storage.get());
            if (db == 0) {
                LSST_AP_THROW(Runtime, "Didn't get DbTsvStorage");
            }
            db->createTableFromTemplate(name, model, false);
            db->setTableForInsert(name);
            IdVector::const_iterator const end(p->end());
            for (IdVector::const_iterator i = p->begin(); i != end; ++i) {
                db->setColumn<int64_t>("id", *i);
                db->insertRow();
            }
        }
    } else {
        LSST_AP_THROW(InvalidParameter, "Storage type is not supported");
    }
}


Persistable * IdVectorFormatter::read(
    Storage::Ptr          storage,
    DataProperty::PtrType additionalData
) {
    std::auto_ptr<IdVector> p(new IdVector);

    if (typeid(*storage) == typeid(BoostStorage)) {
        BoostStorage* bs = dynamic_cast<BoostStorage *>(storage.get());
        if (bs == 0) {
            LSST_AP_THROW(Runtime, "Didn't get BoostStorage");
        }
        bs->getIArchive() & *p;
    } else if (typeid(*storage) == typeid(DbStorage) || typeid(*storage) == typeid(DbTsvStorage)) {
        DbStorage * db = dynamic_cast<DbStorage *>(storage.get());
        if (db == 0) {
            LSST_AP_THROW(Runtime, "Didn't get DbStorage");
        }

        db->setTableForQuery(getTableName(_policy, additionalData));
        int64_t data;
        db->outParam("id", &data);
        db->query();
        while (db->next()) {
            if (db->columnIsNull(0)) {
                LSST_AP_THROW(Runtime, "null column \"id\"");
            }
            p->push_back(data);
        }
        db->finishQuery();
    } else {
        LSST_AP_THROW(InvalidParameter, "Storage type is not supported");
    }
    return p.release();
}


void IdVectorFormatter::update(Persistable *, Storage::Ptr, DataProperty::PtrType) {
    LSST_AP_THROW(Runtime, "IdVectorFormatter: updates not supported");
}


}}} // end of namespace lsst::ap::io

