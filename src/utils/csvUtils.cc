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
  * @brief Implementation of catalog to CSV file conversion.
  *
  * @ingroup ap
  * @author Serge Monkewitz
  */
#include "lsst/ap/utils/csvUtils.h"

#include "boost/ref.hpp"
#include "boost/type_traits/is_integral.hpp"
#include "boost/type_traits/is_same.hpp"

#include "lsst/afw/geom/Angle.h"
#include "lsst/afw/coord.h"

#include "lsst/ap/utils/Csv.h"

using lsst::afw::geom::Angle;

using lsst::afw::coord::Coord;
using lsst::afw::coord::IcrsCoord;

using lsst::afw::table::Array;
using lsst::afw::table::BaseCatalog;
using lsst::afw::table::BaseRecord;
using lsst::afw::table::Covariance;
using lsst::afw::table::Field;
using lsst::afw::table::Flag;
using lsst::afw::table::Key;
using lsst::afw::table::Moments;
using lsst::afw::table::Point;
using lsst::afw::table::Schema;
using lsst::afw::table::SchemaItem;


namespace lsst { namespace ap { namespace utils {

namespace {

    typedef std::vector<std::string>::const_iterator StringIter;

    uint8_t const SKIP_FIELD = 0x1;
    uint8_t const NULLABLE_INT = 0x2;
    uint8_t const RAD2 = 0x4;
    uint8_t const RAD4 = 0x8;

    // computed to 30 digit precision with Maple
    double const ARCSEC2_PER_RAD2 = 4.25451702961521995803707621197e10;
    double const ARCSEC4_PER_RAD4 = 1.81009151552859144174773079111e21;

    int nFlagInts(int count) {
        return (count + 62) / 63;
    }

    // Helper class for CsvConverter initialization
    class CsvConverterInit {
    public:
        CsvConverterInit(CsvConversionControl const & control, uint8_t * fieldFlags) :
            _control(control), _fieldFlags(fieldFlags), _curField(0) { }

        ~CsvConverterInit() { }

        template <typename T>
        void operator()(SchemaItem<T> const & item) const {
            std::string const name = item.field.getName();
            StringIter i = _control.ignoreFields.begin();
            StringIter e =  _control.ignoreFields.end();
            uint8_t c = 0;
            if (std::find(i, e, name) != e) {
                c |= SKIP_FIELD; 
            } else if (boost::is_integral<T>::value) {
                i = _control.nullableIntegers.begin();
                e = _control.nullableIntegers.end();
                if (std::find(i, e, name) != e) {
                    c |= NULLABLE_INT;
                }
            } else if (item.field.getUnits() == "rad^2" &&
                       (boost::is_same<T, Covariance<Point<float> > >::value ||
                        boost::is_same<T, Covariance<Point<double> > >::value ||
                        boost::is_same<T, Moments<float> >::value ||
                        boost::is_same<T, Moments<double> >::value)) {
                // Flag point covariance matrixes and moments with units of
                // rad^2 for unit conversion to arcsec^2.
                c |= RAD2;
            } else if (item.field.getUnits() == "rad^4" &&
                       (boost::is_same<T, Covariance<Moments<float> > >::value ||
                        boost::is_same<T, Covariance<Moments<double> > >::value)) {
                // Flag moment covariance matrixes with units of rad^4 for
                // unit conversion to arcsec^4.
                c |= RAD4;
            }
            _fieldFlags[_curField++] = c;
        }

    private:
        CsvConversionControl const & _control;
        uint8_t mutable * _fieldFlags;
        int mutable _curField;
    };

    // CSV conversion functor invoked for every item in a catalog's schema,
    // one catalog record at a time.
    class CsvConverter {
    public:
        CsvConverter(CsvWriter & writer,
                     Schema const & schema,
                     CsvConversionControl const & control) :
            _writer(&writer),
            _record(),
            _flagsAsBits(control.flagsAsBits),
            _nFlagInts(nFlagInts(schema.getFlagFieldCount())),
            _curField(0),
            _curFlag(0),
            _flagInts(new int64_t[_nFlagInts]),
            _fieldFlags(new uint8_t[schema.getFieldCount()]),
            _canonicalFlags()
        {
            std::memset(_fieldFlags.get(), 0,
                        schema.getFieldCount() * sizeof(uint8_t));
            CsvConverterInit init(control, _fieldFlags.get());
            schema.forEach(boost::ref(init));
            for (StringIter i = control.canonicalFlags.begin(),
                 e = control.canonicalFlags.end(); i != e; ++i) {
                _canonicalFlags.push_back(schema.find<Flag>(*i).key);
            }
        }

        ~CsvConverter() { }

        // Set the record to extract data from
        void startRecord(boost::shared_ptr<BaseRecord const> const &record) {
            _record = record;
            // zero out flag words
            if (!_flagsAsBits) {
                std::memset(_flagInts.get(), 0, _nFlagInts * sizeof(int64_t));
                _curFlag = 0;
            }
            _curField = 0;
        }

        // Finish conversion process for a record 
        void endRecord() {
            typedef std::vector<Key<Flag > >::const_iterator FlagIter;
            if (!_flagsAsBits) {
                // write our flag integers
                for (int i = 0; i < _nFlagInts; ++i) {
                    _writer->appendField(_flagInts[i]);
                }
                // write out canonical flag integers
                int64_t f = 0;
                int b = 0;
                for (FlagIter i = _canonicalFlags.begin(), e = _canonicalFlags.end();
                     i != e; ++i) {
                    if (_record->get(*i)) {
                        f |= static_cast<int64_t>(1) << b;
                    }
                    ++b;
                    if (b == 63) {
                        _writer->appendField(f);
                        f = 0;
                        b = 0;
                    }
                }
                if (b > 0) {
                    _writer->appendField(f);
                }
            }
            _writer->endRecord();
            _record.reset();
        }

        // -- Field conversion --------

        // basic types
        template <typename T>
        void operator()(SchemaItem<T> const & item) const {
            uint8_t const c = _fieldFlags[_curField++];
            if ((c & SKIP_FIELD) != 0) {
                return;
            }
            T val = _record->get(item.key);
            if (boost::is_integral<T>::value && (c & NULLABLE_INT) != 0 &&
                val == static_cast<T>(0)) {
                // T is integral, val is 0, and item is a nullable integer field
                _writer->appendNull();
            } else {
                _writer->appendField(val);
            }
        }

        // angles and coordinates
        void operator()(SchemaItem<Angle> const & item) const {
            if ((_fieldFlags[_curField++] & SKIP_FIELD) != 0) {
                return;
            }
            _writer->appendField(_record->get(item.key).asDegrees());
        }

        void operator()(SchemaItem<Coord> const & item) const {
            if ((_fieldFlags[_curField++] & SKIP_FIELD) != 0) {
                return;
            }
            Field<Coord>::Value c = _record->get(item.key);
            _writer->appendField(c.getLongitude().asDegrees());
            _writer->appendField(c.getLatitude().asDegrees());
        }

        // arrays
        template <typename T>
        void operator()(SchemaItem<Array<T> > const & item) const {
            if ((_fieldFlags[_curField++] & SKIP_FIELD) != 0) {
                return;
            }
            typename Field<Array<T> >::Value a = _record->get(item.key);
            for (int i = 0, e = item.field.getSize(); i != e; ++i) {
                _writer->appendField(a[i]);
            }
        }

        // points
        template <typename T>
        void operator()(SchemaItem<Point<T> > const & item) const {
            if ((_fieldFlags[_curField++] & SKIP_FIELD) != 0) {
                return;
            }
            typename Field<Point<T> >::Value p = _record->get(item.key);
            _writer->appendField(p.getX());
            _writer->appendField(p.getY());
        }

        // moments
        template <typename T>
        void operator()(SchemaItem<Moments<T> > const & item) const {
            uint8_t const c = _fieldFlags[_curField++];
            if ((c & SKIP_FIELD) != 0) {
                return;
            }
            typename Field<Moments<T> >::Value p = _record->get(item.key);
            double const scale = ((c & RAD2) == 0 ? 1.0 : ARCSEC2_PER_RAD2);
            _writer->appendField(scale * p.getIxx());
            _writer->appendField(scale * p.getIyy());
            _writer->appendField(scale * p.getIxy());
        }

        // covariances
        template <typename T>
        void operator()(SchemaItem<Covariance<T> > const & item) const {
            uint8_t const c = _fieldFlags[_curField++];
            if ((c & SKIP_FIELD) != 0) {
                return;
            }
            typename Field<Covariance<T> >::Value cov = _record->get(item.key);
            int const sz = item.field.getSize();
            double scale = 1.0;
            if ((c & RAD2) != 0) {
                scale = ARCSEC2_PER_RAD2;
            } else if ((c & RAD4) != 0) {
                scale = ARCSEC4_PER_RAD4;
            }
            for (int i = 0; i < sz; ++i) {
                for (int j = i; j < sz; ++j) {
                    _writer->appendField(scale * cov(i, j));
                }
            }
        }

        // flags
        void operator()(SchemaItem<Flag> const & item) const {
            if ((_fieldFlags[_curField++] & SKIP_FIELD) != 0) {
                return;
            }
            if (_flagsAsBits) {
                _writer->appendField(_record->get(item.key));
            } else {
                if (_record->get(item.key)) {
                    _flagInts[_curFlag / 63] |= static_cast<int64_t>(1) << (_curFlag % 63);
                }
                ++_curFlag;
            }
        }

    private:
        // disable copy construction/assignment
        CsvConverter(CsvConverter const &);
        CsvConverter& operator=(CsvConverter const &);

        CsvWriter mutable * _writer;
        boost::shared_ptr<BaseRecord const> _record;
        bool _flagsAsBits;
        int _nFlagInts;
        int mutable _curField;
        int mutable _curFlag;
        boost::scoped_array<int64_t> mutable _flagInts;
        boost::scoped_array<uint8_t> _fieldFlags;
        std::vector<Key<Flag> > _canonicalFlags;
    };

} // namespace <anonymous>


CsvConversionControl::CsvConversionControl() :
    flagsAsBits(true),
    ignoreFields(),
    nullableIntegers(),
    canonicalFlags()
{ }

CsvConversionControl::~CsvConversionControl() { }


void writeCsv(BaseCatalog const & catalog,
              CsvConversionControl const & cnvControl,
              CsvControl const & csvControl,
              std::string const & csvFile,
              bool truncate,
              bool append)
{
    CsvWriter writer(csvFile, csvControl, truncate, append);
    Schema schema = catalog.getSchema();
    CsvConverter cnv(writer, schema, cnvControl);
    boost::reference_wrapper<CsvConverter> cnvRef = boost::ref(cnv);
    for (BaseCatalog::const_iterator i = catalog.begin(), e = catalog.end(); i != e; ++i) {
        cnv.startRecord(i);
        schema.forEach(cnvRef);
        cnv.endRecord();
    }
    writer.flush();
}

}}} // namespace lsst:ap::utils
