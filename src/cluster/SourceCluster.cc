// -*- lsst-c++ -*-
/* 
 * LSST Data Management System
 * Copyright 2012 LSST Corporation.
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

/** @file
  * @brief Source cluster table and record class implementation.
  *
  * @ingroup ap
  * @author Serge Monkewitz
  */

#include "lsst/ap/cluster/SourceCluster.h"


namespace lsst { namespace ap { namespace cluster {

namespace except = lsst::pex::exceptions;
namespace table = lsst::afw::table;


// -- SourceClusterRecordImpl and SourceClusterTableImpl --------
//
// These are a private table/record pair -- they're what you actually get when
// you do SourceClusterTable::make(), but are hidden to avoid the friending that
// would be necessary if they had to make their constructors private or protected.

namespace {

    class SourceClusterTableImpl;

    class SourceClusterRecordImpl : public SourceClusterRecord {
    public:
        explicit SourceClusterRecordImpl(PTR(SourceClusterTable) const & table) :
            SourceClusterRecord(table) { }
    };

    class SourceClusterTableImpl : public SourceClusterTable {
    public:
        explicit SourceClusterTableImpl(
            lsst::afw::table::Schema const & schema,
            PTR(lsst::afw::table::IdFactory) const & idFactory
        ) : SourceClusterTable(schema, idFactory) { }

        SourceClusterTableImpl(SourceClusterTableImpl const & other) :
            SourceClusterTable(other) { }

    private:
        virtual PTR(lsst::afw::table::BaseTable) _clone() const {
            return boost::make_shared<SourceClusterTableImpl>(*this);
        }

        virtual PTR(lsst::afw::table::BaseRecord) _makeRecord() {
            PTR(SourceClusterRecord) record =
                boost::make_shared<SourceClusterRecordImpl>(
                    getSelf<SourceClusterTableImpl>());
            if (getIdFactory()) {
                record->setId((*getIdFactory())());
            }
            return record;
        }
    };

} // namespace <anonymous>


// -- SourceClusterRecord implementation --------

SourceClusterRecord::~SourceClusterRecord() { }

SourceClusterRecord::SourceClusterRecord(PTR(SourceClusterTable) const & table) :
    lsst::afw::table::SimpleRecord(table) { }


// -- SourceClusterTable implementation --------

PTR(SourceClusterTable) SourceClusterTable::make(
    lsst::afw::table::Schema const & schema,
    PTR(lsst::afw::table::IdFactory) const & idFactory)
{
    if (!checkSchema(schema)) {
        throw LSST_EXCEPT(except::InvalidParameterException,
            "Schema for SourceClusterTable must contain at least the keys "
            "defined by getMinimalSchema().");
    }
    return boost::make_shared<SourceClusterTableImpl>(schema, idFactory);
}

SourceClusterTable::SourceClusterTable(
    lsst::afw::table::Schema const & schema,
    PTR(lsst::afw::table::IdFactory) const & idFactory
) : lsst::afw::table::SimpleTable(schema, idFactory),
    _keyCoordErr(),
    _keyCoord2(),
    _keyCoord2Err(),
    _keyNumSources(),
    _keyEarliestTime(),
    _keyLatestTime(),
    _filterSlots()
{ }

SourceClusterTable::SourceClusterTable(SourceClusterTable const & other) :
    lsst::afw::table::SimpleTable(other),
    _keyCoordErr(other._keyCoordErr),
    _keyCoord2(other._keyCoord2),
    _keyCoord2Err(other._keyCoord2Err),
    _keyNumSources(other._keyNumSources),
    _keyEarliestTime(other._keyEarliestTime),
    _keyLatestTime(other._keyLatestTime),
    _filterSlots(other._filterSlots)
{ }

SourceClusterTable::~SourceClusterTable() { }

PTR(lsst::afw::table::io::FitsWriter) SourceClusterTable::makeFitsWriter(
    lsst::afw::table::io::FitsWriter::Fits * fits) const
{
    //TODO return boost::make_shared<SourceClusterFitsWriter>(fits);
    return PTR(lsst::afw::table::io::FitsWriter)();
}

SourceClusterTable::FilterSlots const & SourceClusterTable::getFilterSlots(std::string const & filter) const {
    FilterSlotsMap::const_iterator i = _filterSlots.find(filter);
    if (i == _filterSlots.end()) {
        throw LSST_EXCEPT(except::NotFoundException, "SourceClusterTable "
            "contains no slot mappings for the filter named " + filter);
    }
    return i->second;
}

SourceClusterTable::FilterSlots::FilterSlots() :
    keyEarliestTime(),
    keyLatestTime(),
    keyNumSources(),
    keyPsfFlux(),
    keyModelFlux(),
    keyApFlux(),
    keyInstFlux(),
    keyShape()
{ }

SourceClusterTable::FilterSlots::~FilterSlots() { }

// -- SourceClusterIdFactory implementation --------

SourceClusterIdFactory::SourceClusterIdFactory(int skyTileId) :
    _id(static_cast<lsst::afw::table::RecordId>(skyTileId) << 32),
    _skyTileId(skyTileId)
{ }

SourceClusterIdFactory::~SourceClusterIdFactory() { }

lsst::afw::table::RecordId SourceClusterIdFactory::operator()() {
    lsst::afw::table::RecordId id = _id + 1;
    if (static_cast<int>(id >> 32) != _skyTileId) {
        throw LSST_EXCEPT(except::OverflowErrorException,
            "Source cluster ID space exhausted! Note that SourceClusterIdFactory "
            "can hand out a maximum of 2^32 - 1 IDs for a given sky-tile. If there "
            "are more than that many clusters, the sky-tile size must be reduced.");
    }
    _id = id;
    return id;
}

void SourceClusterIdFactory::notify(lsst::afw::table::RecordId id) {
    throw LSST_EXCEPT(except::LogicErrorException,
        "SourceClusterIdFactory does not support the notify() method "
        "of lsst::afw::table::IdFactory");
}


}}} // namespace lsst::ap::cluster

// -- Explicit instantiations

template class lsst::afw::table::CatalogT<lsst::ap::cluster::SourceClusterRecord>;
template class lsst::afw::table::CatalogT<lsst::ap::cluster::SourceClusterRecord const>;
template class lsst::afw::table::SimpleCatalogT<lsst::ap::cluster::SourceClusterRecord>;
template class lsst::afw::table::SimpleCatalogT<lsst::ap::cluster::SourceClusterRecord const>;

