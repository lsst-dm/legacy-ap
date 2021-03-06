changecom(`###')dnl
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

// THIS FILE IS AUTOMATICALLY GENERATED, AND WILL BE OVERWRITTEN IF EDITED MANUALLY.

/** @file
  * @brief Table and record classes for source cluster attributes.
  *
  * @ingroup ap
  * @author Serge Monkewitz
  */
define(`m4def', defn(`define'))dnl
m4def(`DECLARE_SLOT_ACCESSORS',
`/// @brief Get $3.
    $2 const get$1() const;
    /// @brief Set $3.
    void set$1($2 const & value);
')dnl
m4def(`DECLARE_FILTER_SLOT_ACCESSORS',
`/// @brief Get $3.
    $2 const get$1(std::string const & filter) const;
    /// @brief Set $3.
    void set$1(std::string const & filter, $2 const & value);
')dnl
m4def(`DECLARE_COMPOUND_SLOT_ACCESSORS',
`DECLARE_FILTER_SLOT_ACCESSORS(`$1$2', `lsst::afw::table::$2::MeasValue',
        `the inverse variance weighted mean of the $1$2 slot measurement in the given filter')
    DECLARE_FILTER_SLOT_ACCESSORS(`$1$2Err', `lsst::afw::table::$2::ErrValue',
        `the uncertainty of the $1$2 mean in the given filter')
    DECLARE_FILTER_SLOT_ACCESSORS(`$1$2Count', `int',
        `the number of measurements used to compute the $1$2 mean in the given filter')
')dnl
m4def(`DECLARE_FLUX_ACCESSORS', `DECLARE_COMPOUND_SLOT_ACCESSORS($1, `Flux')')dnl
m4def(`DECLARE_SHAPE_ACCESSORS', `DECLARE_COMPOUND_SLOT_ACCESSORS(`', `Shape')')dnl
m4def(`DEFINE_SLOT_ACCESSORS',
`inline $2 const SourceClusterRecord::get$1() const {
     return this->get(getTable()->get$1Key());
}
inline void SourceClusterRecord::set$1($2 const & value) {
     this->set(getTable()->get$1Key(), value);
}
')dnl
m4def(`DEFINE_FILTER_SLOT_ACCESSORS',
`inline $2 const SourceClusterRecord::get$1(std::string const & filter) const {
     return this->get(getTable()->get$1Key(filter));
}
inline void SourceClusterRecord::set$1(std::string const & filter, $2 const & value) {
     this->set(getTable()->get$1Key(filter), value);
}
')dnl
m4def(`DEFINE_COMPOUND_SLOT_ACCESSORS',
`DEFINE_FILTER_SLOT_ACCESSORS(`$1$2', `lsst::afw::table::$2::MeasValue')
DEFINE_FILTER_SLOT_ACCESSORS(`$1$2Err', `lsst::afw::table::$2::ErrValue')
DEFINE_FILTER_SLOT_ACCESSORS(`$1$2Count', `int')
')dnl
m4def(`DEFINE_FLUX_ACCESSORS', `DEFINE_COMPOUND_SLOT_ACCESSORS($1, `Flux')')dnl
m4def(`DEFINE_SHAPE_ACCESSORS', `DEFINE_COMPOUND_SLOT_ACCESSORS(`', `Shape')')dnl
m4def(`DECLARE_SLOT_DEFINERS',
`/// @brief Set the field used for the $1 slot using keys.
    void define$1(lsst::afw::table::Key<$2> const & key) {
        _key$1 = key;
    }
    /// @brief Set the field used for the $1 slot using a name.
    void define$1(std::string const & name) {
        lsst::afw::table::Schema schema = getSchema();
        _key$1 = schema[name];
    }
    /// @brief Get the name of the field used for the $1 slot.
    std::string const get$1Definition() const {
        return getSchema().find(_key$1).field.getName();
    }
    /// @brief Get the key used for the $1 slot.
    lsst::afw::table::Key<$2> const get$1Key() const {
        return _key$1;
    }
')dnl
m4def(`DECLARE_FILTER_SLOT_DEFINERS',
`/// @brief Set the field used for the $1 slot in the given filter using a Key.
    void define$1(std::string const & filter,
                  lsst::afw::table::Key<$2> const & key) {
        _filterSlots[filter].key$1 = key;
    }
    /// @brief Set the $1 slot in the given filter to the field named "<filter>.<name>".
    void define$1(std::string const & filter, std::string const & name) {
        lsst::afw::table::Schema schema = getSchema();
        _filterSlots[filter].key$1 = schema[filter][name];
    }
    /// @brief Return the name of the field used for the $1 slot in the given filter.
    std::string const get$1Definition(std::string const & filter) const {
        return getSchema().find(getFilterSlots(filter).key$1).field.getName();
    }
    /// @brief Return the key used for the $1 slot in the given filter.
    lsst::afw::table::Key<$2> const get$1Key(std::string const & filter) const {
        return getFilterSlots(filter).key$1;
    }
')dnl
m4def(`DECLARE_COMPOUND_SLOT_DEFINERS',
`/// @brief Set the fields used for the $1$2 slot in the given filter using Keys.
    void define$1$2(std::string const & filter,
                    lsst::afw::table::$2::MeasKey const & mean,
                    lsst::afw::table::$2::ErrKey const & err,
                    lsst::afw::table::Key<int> const & count) {
        _filterSlots[filter].key$1$2 = KeyTuple<lsst::afw::table::$2>(mean, err, count);
    }
    /// @brief Set the fields used for the $1$2 slot in the given filter to the fields named
    ///        "<filter>.<name>", "<filter>.<name>.err", and "<filter>.<name>.count".
    void define$1$2(std::string const & filter, std::string const & name) {
        lsst::afw::table::Schema schema = getSchema();
        _filterSlots[filter].key$1$2 = KeyTuple<lsst::afw::table::$2>(
            schema[filter][name],
            schema[filter][name]["err"],
            schema[filter][name]["count"]);
    }
    /// @brief Return the name of the field used for the $1$2 slot in the given filter.
    std::string const get$1$2Definition(std::string const & filter) const {
        return getSchema().find(getFilterSlots(filter).key$1$2.mean).field.getName();
    }
    /// @brief Return the key used for the $1$2 slot in the given filter.
    lsst::afw::table::$2::MeasKey const get$1$2Key(std::string const & filter) const {
        return getFilterSlots(filter).key$1$2.mean;
    }
    /// @brief Return the key used for $1$2 slot error or covariance.
    lsst::afw::table::$2::ErrKey const get$1$2ErrKey(std::string const & filter) const {
        return getFilterSlots(filter).key$1$2.err;
    }
    /// @brief Return the key used for the $1$2 slot count.
    lsst::afw::table::Key<int> const get$1$2CountKey(std::string const & filter) const {
        return getFilterSlots(filter).key$1$2.count;
    }
')dnl
m4def(`DECLARE_FLUX_DEFINERS', `DECLARE_COMPOUND_SLOT_DEFINERS($1, `Flux')')dnl
m4def(`DECLARE_SHAPE_DEFINERS', `DECLARE_COMPOUND_SLOT_DEFINERS(`', `Shape')')dnl
#ifndef LSST_AP_CLUSTER_SOURCECLUSTER_H
#define LSST_AP_CLUSTER_SOURCECLUSTER_H

#include "lsst/tr1/unordered_map.h"

#include "lsst/afw/geom.h"
#include "lsst/afw/coord.h"
#include "lsst/afw/table/Source.h"
#include "lsst/afw/table/IdFactory.h"
#include "lsst/afw/table/Catalog.h"
#include "lsst/afw/table/BaseColumnView.h"
#include "lsst/afw/table/io/FitsWriter.h"


namespace lsst { namespace ap { namespace cluster {

class SourceClusterTable;

template <typename RecordT> class SourceClusterColumnViewT;


#if !defined(SWIG)

/// @brief A three-element tuple of mean, uncertainty, and count keys.
template <typename MeasurementT>
struct KeyTuple {
    typename MeasurementT::MeasKey mean; ///< Key used for the mean measured value.
    typename MeasurementT::ErrKey err;   ///< Key used for the uncertainty.
    lsst::afw::table::Key<int> count;    ///< Key used for the sample count.

    /// Default-constructor; all keys will be invalid.
    KeyTuple() {}

    KeyTuple(
        typename MeasurementT::MeasKey const & mean_,
        typename MeasurementT::ErrKey const & err_,
        lsst::afw::table::Key<int> const & count_
    ) : mean(mean_), err(err_), count(count_) { }
};

/// Convenience function to setup fields for shapes.
KeyTuple<lsst::afw::table::Shape> addShapeFields(
    lsst::afw::table::Schema & schema,
    std::string const & filter,
    std::string const & name,
    std::string const & doc);

/// Convenience function to setup fields for fluxes.
KeyTuple<lsst::afw::table::Flux> addFluxFields(
    lsst::afw::table::Schema & schema,
    std::string const & filter,
    std::string const & name,
    std::string const & doc,
    std::string const & unit);

#endif // !SWIG


/** @brief Record class that contains measurement averages on clusters of
  *        single exposure sources.
  *
  * A system of aliases (called slots) in which a SourceClusterTable instance
  * stores keys for particular measurements is provided. SourceClusterRecord 
  * uses these keys to provide custom getters and setters.  These are not 
  * separate fields, but rather aliases that can point to custom fields.
  */
class SourceClusterRecord : public lsst::afw::table::SimpleRecord {
public:
    typedef SourceClusterTable Table;
    typedef SourceClusterColumnViewT<SourceClusterRecord> ColumnView;
    typedef lsst::afw::table::SortedCatalogT<SourceClusterRecord> Catalog;
    typedef lsst::afw::table::SortedCatalogT<SourceClusterRecord const> ConstCatalog;

    CONST_PTR(SourceClusterTable) getTable() const {
        return boost::static_pointer_cast<SourceClusterTable const>(
            lsst::afw::table::BaseRecord::getTable());
    }

    //@{
    /// @brief Convenience accessors for filter-agnostic keys.

    DECLARE_SLOT_ACCESSORS(`CoordErr', `Eigen::Matrix<float,2,2>',
        `the uncertainty of the sky-coordinates of the cluster')
    DECLARE_SLOT_ACCESSORS(`WeightedMeanCoord', `lsst::afw::coord::IcrsCoord',
        `the inverse variance weighted mean sky-coordinates of the cluster')
    DECLARE_SLOT_ACCESSORS(`WeightedMeanCoordErr', `Eigen::Matrix<float,2,2>',
        `the uncertainty of the WeightedMeanCoord slot')
    DECLARE_SLOT_ACCESSORS(`WeightedMeanCoordCount', `int',
        `the number of measurements used to compute the WeightedMeanCoord slot')
    DECLARE_SLOT_ACCESSORS(`NumSources', `int',
        `the number of sources in the cluster')
    DECLARE_SLOT_ACCESSORS(`TimeMin', `double',
        `the earliest observation time [MJD TAI] of sources in the cluster')
    DECLARE_SLOT_ACCESSORS(`TimeMean', `double',
        `the mean observation time [MJD TAI] of sources in the cluster')
    DECLARE_SLOT_ACCESSORS(`TimeMax', `double',
        `the latest observation time [MJD TAI] of sources in the cluster')
    //@}


    //@{
    /// @brief Convenience accessors for filter-specific keys.

    DECLARE_FILTER_SLOT_ACCESSORS(`NumSources', `int',
        `the number of sources in the given filter')
    DECLARE_FILTER_SLOT_ACCESSORS(`TimeMin', `double',
        `the earliest observation time [MJD TAI] of sources in the given filter')
    DECLARE_FILTER_SLOT_ACCESSORS(`TimeMax', `double',
        `the latest observation time [MJD TAI] of sources in the given filter')

    DECLARE_FLUX_ACCESSORS(`Psf')
    DECLARE_FLUX_ACCESSORS(`Model')
    DECLARE_FLUX_ACCESSORS(`Ap')
    DECLARE_FLUX_ACCESSORS(`Inst')
    DECLARE_SHAPE_ACCESSORS
    /// @brief Return the shape slot Ixx value for the given filter.
    double getIxx(std::string const & filter) const;

    /// @brief Return the shape slot Iyy value for the given filter.
    double getIyy(std::string const & filter) const;

    /// @brief Return the shape slot Ixy value for the given filter.
    double getIxy(std::string const & filter) const;

    //@}

protected:
    SourceClusterRecord(PTR(SourceClusterTable) const & table);
};


/** @brief Table class that contains measurement means on clusters of
  *        single exposure sources.
  *
  * Note that the minimal schema for a SourceClusterTable / 
  * SourceClusterRecord is identical to the minimal schema for
  * an lsst::afw::table::SimpleTable / lsst::afw::table::SimpleRecord.
  */
class SourceClusterTable : public lsst::afw::table::SimpleTable {
public:
    typedef SourceClusterRecord Record;
    typedef SourceClusterColumnViewT<SourceClusterRecord> ColumnView;
    typedef lsst::afw::table::SortedCatalogT<Record> Catalog;
    typedef lsst::afw::table::SortedCatalogT<Record const> ConstCatalog;

    /** @brief Construct a new table.
      *
      * @param[in] schema       Schema that defines the fields, offsets, and 
      *                         record size for the table.
      * @param[in] idFactory    Factory class to generate record IDs when they
      *                         are not explicitly given. If null, record IDs
      *                         will default to zero.
      *
      * Note that not passing an lsst::afw::table::IdFactory at all will call
      * the other override of make(), which will set the ID factory to
      * lsst::afw::table::IdFactory::makeSimple().
      */
    static PTR(SourceClusterTable) make(lsst::afw::table::Schema const & schema,
                                        PTR(lsst::afw::table::IdFactory) const & idFactory);

    /** @brief Construct a new table.
      *
      * @param[in] schema       Schema that defines the fields, offsets, and
      *                         record size for the table.
      *
      * The ID factory will be set to lsst::afw::table::IdFactory::makeSimple().
      */
    static PTR(SourceClusterTable) make(lsst::afw::table::Schema const & schema) {
        return make(schema, lsst::afw::table::IdFactory::makeSimple());
    }

    ~SourceClusterTable();

    /// @copydoc lsst::afw::table::BaseTable::clone
    PTR(SourceClusterTable) clone() const {
        return boost::static_pointer_cast<SourceClusterTable>(_clone());
    }

    /// @copydoc lsst::afw::table::BaseTable::makeRecord
    PTR(SourceClusterRecord) makeRecord() {
        return boost::static_pointer_cast<SourceClusterRecord>(_makeRecord());
    }

    /// @copydoc lsst::afw::table::BaseTable::copyRecord
    PTR(SourceClusterRecord) copyRecord(lsst::afw::table::BaseRecord const & other) {
        return boost::static_pointer_cast<SourceClusterRecord>(
            lsst::afw::table::BaseTable::copyRecord(other));
    }

    /// @copydoc lsst::afw::table::BaseTable::copyRecord
    PTR(SourceClusterRecord) copyRecord(lsst::afw::table::BaseRecord const & other,
                                        lsst::afw::table::SchemaMapper const & mapper) {
        return boost::static_pointer_cast<SourceClusterRecord>(
            lsst::afw::table::BaseTable::copyRecord(other, mapper));
    }

    //@{
    /// @brief Convenience definers for filter-agnostic keys.

    DECLARE_SLOT_DEFINERS(`CoordErr', `lsst::afw::table::Covariance<lsst::afw::table::Point<float> > ')
    DECLARE_SLOT_DEFINERS(`WeightedMeanCoord', `lsst::afw::coord::Coord')
    DECLARE_SLOT_DEFINERS(`WeightedMeanCoordErr', `lsst::afw::table::Covariance<lsst::afw::table::Point<float> > ')
    DECLARE_SLOT_DEFINERS(`WeightedMeanCoordCount', `int')
    DECLARE_SLOT_DEFINERS(`NumSources', `int')
    DECLARE_SLOT_DEFINERS(`TimeMin', `double')
    DECLARE_SLOT_DEFINERS(`TimeMean', `double')
    DECLARE_SLOT_DEFINERS(`TimeMax', `double')
    //@}

    /// @brief Get the lexicographically sorted list of filter names for which slots have been defined.
    std::vector<std::string> const getFilters() const;

    //@{
    /// @brief Convenience definers for filter-specific keys.

    DECLARE_FILTER_SLOT_DEFINERS(`NumSources', `int')
    DECLARE_FILTER_SLOT_DEFINERS(`TimeMin', `double')
    DECLARE_FILTER_SLOT_DEFINERS(`TimeMax', `double')
    DECLARE_FLUX_DEFINERS(`Psf')
    DECLARE_FLUX_DEFINERS(`Model')
    DECLARE_FLUX_DEFINERS(`Ap')
    DECLARE_FLUX_DEFINERS(`Inst')
    DECLARE_SHAPE_DEFINERS
    //@}

protected:
    SourceClusterTable(lsst::afw::table::Schema const & schema,
                       PTR(lsst::afw::table::IdFactory) const & idFactory);

    SourceClusterTable(SourceClusterTable const & other);

private:
    struct FilterSlots {
        lsst::afw::table::Key<double>     keyTimeMin;
        lsst::afw::table::Key<double>     keyTimeMax;
        lsst::afw::table::Key<int>        keyNumSources;
        KeyTuple<lsst::afw::table::Flux>  keyPsfFlux;
        KeyTuple<lsst::afw::table::Flux>  keyModelFlux;
        KeyTuple<lsst::afw::table::Flux>  keyApFlux;
        KeyTuple<lsst::afw::table::Flux>  keyInstFlux;
        KeyTuple<lsst::afw::table::Shape> keyShape;

        FilterSlots();
        ~FilterSlots();
    };

    typedef std::tr1::unordered_map<std::string, FilterSlots> FilterSlotsMap;

    // Return a writer object that knows how to save in FITS format.
    virtual PTR(lsst::afw::table::io::FitsWriter) makeFitsWriter(
        lsst::afw::table::io::FitsWriter::Fits * fits, int flags) const;

    FilterSlots const & getFilterSlots(std::string const & filter) const;

    lsst::afw::table::Key<lsst::afw::table::Covariance<lsst::afw::table::Point<float> > > _keyCoordErr;
    lsst::afw::table::Key<lsst::afw::coord::Coord> _keyWeightedMeanCoord;
    lsst::afw::table::Key<lsst::afw::table::Covariance<lsst::afw::table::Point<float> > > _keyWeightedMeanCoordErr;
    lsst::afw::table::Key<int> _keyWeightedMeanCoordCount;
    lsst::afw::table::Key<int> _keyNumSources;
    lsst::afw::table::Key<double> _keyTimeMin;
    lsst::afw::table::Key<double> _keyTimeMean;
    lsst::afw::table::Key<double> _keyTimeMax;

    FilterSlotsMap _filterSlots;

    friend class lsst::afw::table::io::FitsWriter;
};


template <typename RecordT>
class SourceClusterColumnViewT : public lsst::afw::table::ColumnViewT<RecordT> {
public:

    typedef RecordT Record;
    typedef typename RecordT::Table Table;

    ndarray::Array<int const,1> const getNumSources() const {
        return this->operator[](this->getTable()->getNumSourcesKey());
    }
    ndarray::Array<double const,1> const getTimeMin() const {
        return this->operator[](this->getTable()->getTimeMinKey());
    }
    ndarray::Array<double const,1> const getTimeMean() const {
        return this->operator[](this->getTable()->getTimeMeanKey());
    }
    ndarray::Array<double const,1> const getTimeMax() const {
        return this->operator[](this->getTable()->getTimeMaxKey());
    }

    ndarray::Array<int const,1> const getNumSources(std::string const & filter) const {
        return this->operator[](this->getTable()->getNumSourcesKey(filter));
    }
    ndarray::Array<double const,1> const getTimeMin(std::string const & filter) const {
        return this->operator[](this->getTable()->getTimeMinKey(filter));
    }
    ndarray::Array<double const,1> const getTimeMax(std::string const & filter) const {
        return this->operator[](this->getTable()->getTimeMaxKey(filter));
    }

    ndarray::Array<double const,1> const getPsfFlux(std::string const & filter) const {
        return this->operator[](this->getTable()->getPsfFluxKey(filter));
    }
    ndarray::Array<double const,1> const getApFlux(std::string const & filter) const {
        return this->operator[](this->getTable()->getApFluxKey(filter));
    }
    ndarray::Array<double const,1> const getModelFlux(std::string const & filter) const {
        return this->operator[](this->getTable()->getModelFluxKey(filter));
    }
    ndarray::Array<double const,1> const getInstFlux(std::string const & filter) const {
        return this->operator[](this->getTable()->getInstFluxKey(filter));
    }

    ndarray::Array<double const,1> const getIxx(std::string const & filter) const {
        return this->operator[](this->getTable()->getShapeKey(filter).getIxx());
    }
    ndarray::Array<double const,1> const getIyy(std::string const & filter) const {
        return this->operator[](this->getTable()->getShapeKey(filter).getIyy());
    }
    ndarray::Array<double const,1> const getIxy(std::string const & filter) const {
        return this->operator[](this->getTable()->getShapeKey(filter).getIxy());
    }

    /// @brief @copydoc lsst::afw::table::BaseColumnView::make
    template <typename InputIterator>
    static SourceClusterColumnViewT make(PTR(Table) const & table,
                                         InputIterator first,
                                         InputIterator last) {
        return SourceClusterColumnViewT(
            lsst::afw::table::BaseColumnView::make(table, first, last));
    }

protected:
    explicit SourceClusterColumnViewT(lsst::afw::table::BaseColumnView const & base) :
        lsst::afw::table::ColumnViewT<RecordT>(base) { }
};

typedef SourceClusterColumnViewT<SourceClusterRecord> SourceClusterColumnView;


/** Generates at most 2^32 - 1 contiguous record IDs. The upper 32 bits of the
  * record IDs generated by this class are fixed to a specific sky-tile ID.
  *
  * There is no support for notifying an instance that a given ID should not be
  * returned. An ID of 0 will never be returned.
  */
class SourceClusterIdFactory : public lsst::afw::table::IdFactory {
public:
    SourceClusterIdFactory(int skyTileId);
    virtual ~SourceClusterIdFactory();

    virtual lsst::afw::table::RecordId operator()();

    // Always throws!
    virtual void notify(lsst::afw::table::RecordId id);

    virtual PTR(lsst::afw::table::IdFactory) clone() const {
        return boost::make_shared<SourceClusterIdFactory>(*this);
    }

private:
    lsst::afw::table::RecordId _id;
    int _skyTileId;
};


#ifndef SWIG

typedef lsst::afw::table::SortedCatalogT<SourceClusterRecord> SourceClusterCatalog;
typedef lsst::afw::table::SortedCatalogT<SourceClusterRecord const> ConstSourceClusterCatalog;

DEFINE_SLOT_ACCESSORS(`CoordErr', `Eigen::Matrix<float,2,2>')
DEFINE_SLOT_ACCESSORS(`WeightedMeanCoord', `lsst::afw::coord::IcrsCoord')
DEFINE_SLOT_ACCESSORS(`WeightedMeanCoordErr', `Eigen::Matrix<float,2,2>')
DEFINE_SLOT_ACCESSORS(`WeightedMeanCoordCount', `int')
DEFINE_SLOT_ACCESSORS(`NumSources', `int')
DEFINE_SLOT_ACCESSORS(`TimeMin', `double')
DEFINE_SLOT_ACCESSORS(`TimeMean', `double')
DEFINE_SLOT_ACCESSORS(`TimeMax', `double')

DEFINE_FILTER_SLOT_ACCESSORS(`NumSources', `int')
DEFINE_FILTER_SLOT_ACCESSORS(`TimeMin', `double')
DEFINE_FILTER_SLOT_ACCESSORS(`TimeMax', `double')

DEFINE_FLUX_ACCESSORS(`Psf')
DEFINE_FLUX_ACCESSORS(`Model')
DEFINE_FLUX_ACCESSORS(`Ap')
DEFINE_FLUX_ACCESSORS(`Inst')
DEFINE_SHAPE_ACCESSORS
inline double SourceClusterRecord::getIxx(std::string const & filter) const {
    return get(getTable()->getShapeKey(filter).getIxx());
}
inline double SourceClusterRecord::getIyy(std::string const & filter) const {
    return get(getTable()->getShapeKey(filter).getIyy());
}
inline double SourceClusterRecord::getIxy(std::string const & filter) const {
    return get(getTable()->getShapeKey(filter).getIxy());
}

#endif // !SWIG

}}} // namespace lsst::ap::cluster

#endif // !LSST_AP_CLUSTER_SOURCECLUSTER_H
