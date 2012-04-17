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

#include <boost/algorithm/string/case_conv.hpp>


// boilerplate macros for saving/loading slot definitionsa

#define SAVE_SLOT(NAME, Name) \
    if (table->get ## Name ## Key().isValid()) { \
        std::string s = table->getSchema().find(table->get ## Name ## Key()).field.getName(); \
        std::replace(s.begin(), s.end(), '.', '_'); \
        _fits->writeKey(#NAME "_SLOT", s.c_str(), "Defines the " #Name " slot"); \
    }

#define SAVE_FILTER_SLOT(FILTER, filter, NAME, Name) \
    if (table->get ## Name ## Key(filter).isValid()) { \
        std::string s = table->getSchema().find(table->get ## Name ## Key(filter)).field.getName(); \
        std::replace(s.begin(), s.end(), '.', '_'); \
        _fits->writeKey(FILTER + "_" #NAME "_SLOT", s.c_str(), \
                        "Defines the " #Name " slot in the " + filter + " filter"); \
    }

#define SAVE_COMPOUND_SLOT(FILTER, filter, NAME, Name) \
    SAVE_FILTER_SLOT(FILTER, filter, NAME, Name) \
    SAVE_FILTER_SLOT(FILTER, filter, NAME ## _ERR, Name ## Err) \
    SAVE_FILTER_SLOT(FILTER, filter, NAME ## _COUNT, Name ## Count)

#define LOAD_SLOT(NAME, Name) \
    { \
        _fits->behavior &= ~lsst::afw::fits::Fits::AUTO_CHECK; \
        std::string s; \
        _fits->readKey(#NAME "_SLOT", s); \
        if (_fits->status == 0) { \
            metadata->remove(#NAME "_SLOT"); \
            std::replace(s.begin(), s.end(), '_', '.'); \
            table->define ## Name (schema[s]); \
        } else { \
            _fits->status = 0; \
        } \
        _fits->behavior |= lsst::afw::fits::Fits::AUTO_CHECK; \
    }

#define LOAD_FILTER_SLOT(FILTER, filter, NAME, Name) \
    { \
        _fits->behavior &= ~lsst::afw::fits::Fits::AUTO_CHECK; \
        std::string s; \
        _fits->readKey(FILTER + "_" #NAME "_SLOT", s); \
        if (_fits->status == 0) { \
            metadata->remove(FILTER + "_" #NAME "_SLOT"); \
            std::replace(s.begin(), s.end(), '_', '.'); \
            table->define ## Name (filter, schema[s]); \
        } else { \
            _fits->status = 0; \
        } \
        _fits->behavior |= lsst::afw::fits::Fits::AUTO_CHECK; \
    }

#define LOAD_COMPOUND_SLOT(FILTER, filter, NAME, Name) \
    { \
        _fits->behavior &= ~lsst::afw::fits::Fits::AUTO_CHECK; \
        std::string s, sErr, sCount; \
        _fits->readKey(FILTER + "_" #NAME "_SLOT", s); \
        _fits->readKey(FILTER + "_" #NAME "_ERR_SLOT", sErr); \
        _fits->readKey(FILTER + "_" #NAME "_COUNT_SLOT", sCount); \
        if (_fits->status == 0) { \
            metadata->remove(FILTER + "_" #NAME "_SLOT"); \
            metadata->remove(FILTER + "_" #NAME "_ERR_SLOT"); \
            metadata->remove(FILTER + "_" #NAME "_COUNT_SLOT"); \
            std::replace(s.begin(), s.end(), '_', '.'); \
            std::replace(sErr.begin(), sErr.end(), '_', '.'); \
            std::replace(sCount.begin(), sCount.end(), '_', '.'); \
            table->define ## Name (filter, schema[s], schema[sErr], schema[sCount]); \
        } else { \
            _fits->status = 0; \
        } \
        _fits->behavior |= lsst::afw::fits::Fits::AUTO_CHECK; \
    }


namespace except = lsst::pex::exceptions;

using lsst::daf::base::PropertyList;
using lsst::afw::table::BaseTable;
using lsst::afw::table::io::FitsReader;
using lsst::afw::table::io::FitsWriter;


namespace lsst { namespace ap { namespace cluster {

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
        virtual PTR(BaseTable) _clone() const {
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


// -- SourceClusterFitsWriter --------

namespace {

    // A custom FitsWriter for source cluster tables - this adds header keys that
    // define the slots. It also sets the AFW_TYPE key to SOURCE_CLUSTER, which
    // should ensure that SourceClusterFitsReader is used to read it.
    class SourceClusterFitsWriter : public FitsWriter {
    public:
        explicit SourceClusterFitsWriter(lsst::afw::fits::Fits * fits) :
            FitsWriter(fits) { }
    protected:
        virtual void _writeTable(CONST_PTR(BaseTable) const & table);
    };

    void SourceClusterFitsWriter::_writeTable(CONST_PTR(BaseTable) const & t) {
        CONST_PTR(SourceClusterTable) table =
            boost::dynamic_pointer_cast<SourceClusterTable const>(t);
        if (!table) {
            throw LSST_EXCEPT(except::LogicErrorException, "SourceClusterFitsWriter "
                              "can only write out SourceClusterTable instances!");
        }
        PTR(PropertyList) metadata = table->getMetadata();
        // exploit hole in constness of table and remove the FILTERS key if present
        if (metadata && metadata->exists("FILTERS")) {
            metadata->remove("FILTERS");
        }
        FitsWriter::_writeTable(table);
        metadata = boost::make_shared<PropertyList>();
        std::vector<std::string> const filters = table->getFilters();
        if (!filters.empty()) {
            metadata->set<std::string>("FILTERS", filters);
            _fits->writeMetadata(*metadata);
        }
        _fits->writeKey("AFW_TYPE", "SOURCE_CLUSTER",
                        "Tells lsst::afw to load this as a SourceClusterTable.");
        // save filter agnostic slots
        SAVE_SLOT(COORD_ERR, CoordErr)
        SAVE_SLOT(COORD2, WeightedCoord)
        SAVE_SLOT(COORD2_ERR, WeightedCoordErr)
        SAVE_SLOT(NUM_SOURCES, NumSources)
        SAVE_SLOT(TIME_MIN, TimeMin)
        SAVE_SLOT(TIME_MEAN, TimeMean)
        SAVE_SLOT(TIME_MAX, TimeMax)
        // save filters and filter-specific slots
        typedef std::vector<std::string>::const_iterator Iter;
        for (Iter i = filters.begin(), e = filters.end(); i != e; ++i) {
            std::string const f = *i;
            std::string const F = boost::to_upper_copy(f);
            SAVE_FILTER_SLOT(F, f, NUM_SOURCES, NumSources)
            SAVE_FILTER_SLOT(F, f, TIME_MIN, TimeMin)
            SAVE_FILTER_SLOT(F, f, TIME_MAX, TimeMax)
            SAVE_COMPOUND_SLOT(F, f, PSF_FLUX, PsfFlux)
            SAVE_COMPOUND_SLOT(F, f, MODEL_FLUX, ModelFlux)
            SAVE_COMPOUND_SLOT(F, f, AP_FLUX, ApFlux)
            SAVE_COMPOUND_SLOT(F, f, INST_FLUX, InstFlux)
            SAVE_COMPOUND_SLOT(F, f, SHAPE, Shape)
        }
    }

} // namespace <anonymous>


// -- SourceClusterFitsReader --------

namespace {

    // A custom FitsReader for source cluster tables - this reads header keys
    // defining slots. It is registered with the name SOURCE_CLUSTER, so will
    // be used whenever a table with AFW_TYPE set to that value is read.
    class SourceClusterFitsReader : public FitsReader {
    public:
        explicit SourceClusterFitsReader(lsst::afw::fits::Fits * fits) :
            FitsReader(fits) { }
    protected:
        virtual PTR(BaseTable) _readTable();
    };

    PTR(BaseTable) SourceClusterFitsReader::_readTable() {
        PTR(PropertyList) metadata = boost::make_shared<PropertyList>();
        _fits->readMetadata(*metadata, true);
        if (metadata->exists("AFW_TYPE")) {
            metadata->remove("AFW_TYPE");
        }
        lsst::afw::table::Schema schema(*metadata, true);
        PTR(SourceClusterTable) table = SourceClusterTable::make(
            schema, PTR(lsst::afw::table::IdFactory)());
        // read in filter agnostic slots
        LOAD_SLOT(COORD_ERR, CoordErr)
        LOAD_SLOT(COORD2, WeightedCoord)
        LOAD_SLOT(COORD2_ERR, WeightedCoordErr)
        LOAD_SLOT(NUM_SOURCES, NumSources)
        LOAD_SLOT(TIME_MIN, TimeMin)
        LOAD_SLOT(TIME_MEAN, TimeMean)
        LOAD_SLOT(TIME_MAX, TimeMax)
        // read in filter specific slots
        std::vector<std::string> filters;
        if (metadata->exists("FILTERS")) {
            filters = metadata->getArray<std::string>("FILTERS");
            metadata->remove("FILTERS");
        }
        typedef std::vector<std::string>::const_iterator Iter;
        for (Iter i = filters.begin(), e = filters.end(); i != e; ++i) {
            std::string const f = *i;
            std::string const F = boost::to_upper_copy(f);
            LOAD_FILTER_SLOT(F, f, NUM_SOURCES, NumSources)
            LOAD_FILTER_SLOT(F, f, TIME_MIN, TimeMin)
            LOAD_FILTER_SLOT(F, f, TIME_MAX, TimeMax)
            LOAD_COMPOUND_SLOT(F, f, PSF_FLUX, PsfFlux)
            LOAD_COMPOUND_SLOT(F, f, MODEL_FLUX, ModelFlux)
            LOAD_COMPOUND_SLOT(F, f, AP_FLUX, ApFlux)
            LOAD_COMPOUND_SLOT(F, f, INST_FLUX, InstFlux)
            LOAD_COMPOUND_SLOT(F, f, SHAPE, Shape)
        }
        _startRecords(*table);
        table->setMetadata(metadata);
        return table;
    }

    // registers the reader so FitsReader::make can use it.
    static FitsReader::FactoryT<SourceClusterFitsReader> sourceClusterFitsReaderFactory("SOURCE_CLUSTER");

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
    _keyWeightedCoord(),
    _keyWeightedCoordErr(),
    _keyNumSources(),
    _keyTimeMin(),
    _keyTimeMean(),
    _keyTimeMax(),
    _filterSlots()
{ }

SourceClusterTable::SourceClusterTable(SourceClusterTable const & other) :
    lsst::afw::table::SimpleTable(other),
    _keyCoordErr(other._keyCoordErr),
    _keyWeightedCoord(other._keyWeightedCoord),
    _keyWeightedCoordErr(other._keyWeightedCoordErr),
    _keyNumSources(other._keyNumSources),
    _keyTimeMin(other._keyTimeMin),
    _keyTimeMean(other._keyTimeMean),
    _keyTimeMax(other._keyTimeMax),
    _filterSlots(other._filterSlots)
{ }

SourceClusterTable::~SourceClusterTable() { }

std::vector<std::string> const SourceClusterTable::getFilters() const {
    typedef FilterSlotsMap::const_iterator Iter;
    std::vector<std::string> filters;
    for (Iter i = _filterSlots.begin(), e = _filterSlots.end(); i != e; ++i) {
        filters.push_back(i->first);
    }
    return filters;
}

PTR(lsst::afw::table::io::FitsWriter) SourceClusterTable::makeFitsWriter(
    lsst::afw::table::io::FitsWriter::Fits * fits) const
{
    return boost::make_shared<SourceClusterFitsWriter>(fits);
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
    keyTimeMin(),
    keyTimeMax(),
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


// -- UUtility function implementations --------

KeyTuple<lsst::afw::table::Shape> addShapeFields(
    lsst::afw::table::Schema & schema,
    std::string const & filter,
    std::string const & name,
    std::string const & doc)
{
    using lsst::afw::table::Shape;
    KeyTuple<Shape> kt;
    kt.mean = schema.addField<Shape::MeasTag>(
        filter + "." + name, doc, "rad^2");
    kt.err = schema.addField<Shape::ErrTag>(
        filter + "." + name + ".err",
        "covariance matrix for " + filter + "." + name,
        "rad^4");
    kt.count = schema.addField<int>(
        filter + "." + name + ".count",
        "Number of samples used to compute the " + filter + 
        "." + name + " mean");
    return kt;
}

KeyTuple<lsst::afw::table::Flux> addFluxFields(
    lsst::afw::table::Schema & schema,
    std::string const & filter,
    std::string const & name,
    std::string const & doc)
{
    using lsst::afw::table::Flux;
    KeyTuple<Flux> kt;
    kt.mean = schema.addField<Flux::MeasTag>(
        filter + "." + name, doc, "erg/s/cm^2/Hz");
    kt.err = schema.addField<Flux::ErrTag>(
        filter + "." + name + ".err",
        "uncertainty for " + filter + "." + name,
        "erg/s/cm^2/Hz");
    kt.count = schema.addField<int>(
        filter + "." + name + ".count",
        "Number of samples used to compute the " + filter +
        "." + name + " mean");
    return kt;
}

}}} // namespace lsst::ap::cluster

// -- Explicit instantiations

template class lsst::afw::table::CatalogT<lsst::ap::cluster::SourceClusterRecord>;
template class lsst::afw::table::CatalogT<lsst::ap::cluster::SourceClusterRecord const>;
template class lsst::afw::table::SimpleCatalogT<lsst::ap::cluster::SourceClusterRecord>;
template class lsst::afw::table::SimpleCatalogT<lsst::ap::cluster::SourceClusterRecord const>;

