// -*- lsst-c++ -*-
//
//##====----------------                                ----------------====##/
//
//! \file   Filter.cc
//! \brief  Implements looking up a filter identifier by name.
//
//##====----------------                                ----------------====##/

#include <boost/bind.hpp>

#include <lsst/mwi/persistence/Persistence.h>
#include <lsst/mwi/persistence/Storage.h>
#include <lsst/mwi/policy/Policy.h>

#include <lsst/ap/Exceptions.h>
#include <lsst/ap/Filter.h>
#include <lsst/ap/ScopeGuard.h>


namespace lsst {
namespace ap {

using lsst::mwi::persistence::Persistence;
using lsst::mwi::persistence::Storage;
using lsst::mwi::policy::Policy;


/*!
    Creates a Filter with the given name, using the \c Filter table in the database given by
    \a location to map the filter name to an integer identifier.
 */
Filter::Filter(lsst::mwi::persistence::LogicalLocation const & location, std::string const & name) {
    Policy::Ptr      noPolicy;
    Persistence::Ptr persistence = Persistence::getPersistence(noPolicy);
    Storage::Ptr     storage     = persistence->getRetrieveStorage("DbStorage", location);
    lsst::mwi::persistence::DbStorage * db = dynamic_cast<lsst::mwi::persistence::DbStorage *>(storage.get());
    if (db == 0) {
        LSST_AP_THROW(Runtime, "Didn't get DbStorage");
    }
    db->startTransaction();
    ScopeGuard g(boost::bind(&lsst::mwi::persistence::DbStorage::endTransaction, db));
    _id = nameToId(*db, name);
}


/*!
    Returns the name of the filter, using the \b Filter table in the database given by
    \a location to map the filter identifier to a name.
 */
std::string const Filter::toString(lsst::mwi::persistence::LogicalLocation const & location) {
    Policy::Ptr      noPolicy;
    Persistence::Ptr persistence = Persistence::getPersistence(noPolicy);
    Storage::Ptr     storage     = persistence->getRetrieveStorage("DbStorage", location);
    lsst::mwi::persistence::DbStorage * db = dynamic_cast<lsst::mwi::persistence::DbStorage *>(storage.get());
    if (db == 0) {
        LSST_AP_THROW(Runtime, "Didn't get DbStorage");
    }
    db->startTransaction();
    ScopeGuard g(boost::bind(&lsst::mwi::persistence::DbStorage::endTransaction, db));
    return toString(*db);
}


/*!
    Returns the name of the filter, using the \b Filter table in the database currently
    set on the given DbStorage to map the filter identifier to a name.
 */
std::string const Filter::toString(lsst::mwi::persistence::DbStorage & db) {
    db.setTableForQuery("Filter");
    db.outColumn("filtName");
    // CORAL always maps MYSQL_TYPE_LONG (MySQL internal type specifier for INTEGER columns) to long
    db.condParam<long>("id", static_cast<long>(_id));
    db.setQueryWhere("filterId = :id");
    db.query();
    ScopeGuard g(boost::bind(&lsst::mwi::persistence::DbStorage::finishQuery, &db));
    if (!db.next() || db.columnIsNull(0)) {
        LSST_AP_THROW(IoError, "Failed to get name for filter " + _id);
    }
    std::string filterName = db.getColumnByPos<std::string>(0);
    if (db.next()) {
        LSST_AP_THROW(IoError, "Multiple names for filter " + _id);
    }
    return filterName;
}


int Filter::nameToId(lsst::mwi::persistence::DbStorage & db, std::string const & name) {
    db.setTableForQuery("Filter");
    db.outColumn("filterId");
    db.condParam<std::string>("name", name);
    db.setQueryWhere("filtName = :name");
    db.query();
    ScopeGuard g(boost::bind(&lsst::mwi::persistence::DbStorage::finishQuery, &db));
    if (!db.next() || db.columnIsNull(0)) {
        LSST_AP_THROW(IoError, "Failed to get id for filter named " + name);
    }
    int filterId = static_cast<int>(db.getColumnByPos<long>(0));
    if (db.next()) {
        LSST_AP_THROW(IoError, "Multiple ids for filter named " + name);
    }
    if (filterId < U || filterId >= NUM_FILTERS) {
        LSST_AP_THROW(OutOfRange, "Invalid filter id for filter named " + name);
    }
    return filterId;
}


}}  // end of namespace lsst::ap

