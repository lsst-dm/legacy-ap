// -*- lsst-c++ -*-

/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
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
  * @brief  Inlines for CSV I/O classes.
  * @author Serge Monkewitz
  */
#ifndef LSST_AP_UTIL_CSV_CC
#define LSST_AP_UTIL_CSV_CC

#include "Csv.h"

#include <limits>

#include "boost/mpl/or.hpp"
#include "boost/static_assert.hpp"
#include "boost/type_traits/is_integral.hpp"
#include "boost/type_traits/is_floating_point.hpp"
#include "boost/type_traits/is_same.hpp"


namespace lsst { namespace ap { namespace util {

// -- CsvDialect inline members ----

/** Returns the explicit database NULL representation, or an empty string
  * if there is none (some dialects can still encode NULLs as the "\N"
  * escape sequence when this is the case).
  */
inline std::string const & CsvDialect::getNull() const {
    return _null;
}

/** Returns true if getNull() identifies the database NULL representation.
  * @sa isNullRecognizable() const
  */
inline bool CsvDialect::hasNull() const {
    return _hasNull;
}

/** Returns true if database NULLs are recognizable in this dialect.
  * This is the case when hasNull() is true, or when getEscapeChar() is not
  * '\\0' and standardEscapes() is true (in which case "\N" is recognized
  * as a NULL).
  */
inline bool CsvDialect::isNullRecognizable() const {
    return _hasNull || _standardEscapes;
}

/** Returns the dialects quoting style. Note that if this is equal to
  * QUOTE_NONE, no special handling of quote characters on reading is
  * performed.
  */
inline CsvDialect::Quoting CsvDialect::getQuoting() const {
    return _quoting;
}

/** Returns the character used to delimit fields. Embedded delimiters
  * are escaped with getEscapeChar(), unless the escape character is '\\0',
  * in which case an attempt to write a field with an embedded delimiter
  * will throw an exception. '\\[0nr]' are illegal delimiter characters.
  */
inline char CsvDialect::getDelimiter() const {
    return _delimiter;
}

/** Returns the character used to escape special characters. If the escape
  * character is '\\0', an attempt to write a field with an embedded special
  * character will throw an exception. The escape character is not allowed
  * to equal '\\[nr]' or the delimiter character - the only time it can equal
  * the quote character is if both are '\\0'.
  */
inline char CsvDialect::getEscapeChar() const {
    return _escapeChar;
}

/** Returns the character used to quote fields containing special characters.
  * Embedded quote characters are escaped with a leading quote (if
  * doubleQuote() returns true) or with getEscapeChar() otherwise. The quote
  * character is not allowed to equal '\\[nr]' or the delimiter character -
  * the only time it can equal the escape character is if both are '\\0'.
  */
inline char CsvDialect::getQuoteChar() const {
    return _quoteChar;
}

/** If true, whitespace immediately following the delimiter is ignored.
  */
inline bool CsvDialect::skipInitialSpace() const {
    return _skipInitialSpace;
}

/** If true, embedded quote characters are escaped with a leading quote
  * character. Otherwise the escape character is used. If escaping is
  * disabled and double-quoting is disabled, an attempt to write a field
  * with an embedded quote character will throw an exception.
  */
inline bool CsvDialect::doubleQuote() const {
    return _doubleQuote;
}

/** Should standard escapes be handled? If false, then the character sequence
  * "\C", where C is any character, is mapped to C (assuming '\\' is the
  * escape character). If true, behaviour is the the same except for the
  * following special cases:
  *
  * @li "\b" is mapped to BS - backspace (ASCII 8)
  * @li "\f" is mapped to FF - form feed (ASCII 12)
  * @li "\n" is mapped to NL - newline (ASCII 10)
  * @li "\r" is mapped to CR - carriage return (ASCII 13)
  * @li "\t" is mapped to TAB - horizontal tab (ASCII 9)
  * @li "\v" is mapped to VT - vertical tab (ASCII 11)
  * @li "\xD" and "\xDD", where D is a hexadecimal digit, is mapped to
  *     the character with that numeric code.
  * @li A field value of exactly "\N" (no quotes, whitespace, or other
  *     content) is treated as a database NULL.
  */
inline bool CsvDialect::standardEscapes() const {
    return _standardEscapes;
}

/** If true, then a trailing delimiter character is expected and written
  * at end of every record, immediately preceding the line terminator.
  */
inline bool CsvDialect::trailingDelimiter() const {
    return _trailingDelimiter;
}


// -- CsvReader inline/template members ----

/** Returns true if all records have been read and there is no current record.
  */
inline bool CsvReader::isDone() const {
    return _done;
}

/** Advances to the next record in the file. If all records have been read,
  * the function has no effect.
  *
  * @throw lsst::pex::exception::RuntimeErrorException
  *        If this exception is thrown, it is because the input file did not
  *        conform to the readers dialect. The current record will contain
  *        the fields succesfully read-in, but the last field may be
  *        incomplete or otherwise incorrectly decoded. The next call to
  *        nextRecord() will resume reading at the beginning of the next line
  *        in the file. If fields contain new-lines, this will not necessarily
  *        be at the start of a record!
  * @throw lsst::pex::exception::IoErrorException
  *        A system I/O call failed - one cannot recover in any general way.
  * @throw lsst::pex::exception::LogicErrorException
  *        There is a serious bug in the internal CSV parser. File a ticket!
  */
inline void CsvReader::nextRecord() {
    if (!_done) {
        _readRecord();
    }
}

/** Returns the number of lines read in. This is the 1-based index of the
  * last line in the current record, and is 0 only when the input stream
  * is empty.
  */
inline size_t CsvReader::getNumLines() const {
    return _numLines;
}

/** Returns the number of records read in. This is the 1-based index of the
  * current record, and is 0 only when the input stream is empty.
  */
inline size_t CsvReader::getNumRecords() const {
    return _numRecords;
}

/** Returns the field names defined for this CsvReader. The i-th name
  * identifies the i-th field in a record.
  */
inline std::vector<std::string> const & CsvReader::getFieldNames() const {
    return _names;
}

//@{
/** Returns the 0-based index (in a record) of the field with the specified
  * name, or -1 if the field name is not recognized.
  */
inline int CsvReader::getIndexOf(std::string const &name) const {
    FieldIndexes::const_iterator i = _indexes.find(name);
    return i == _indexes.end() ? -1 : static_cast<int>(i->second);
}
inline int CsvReader::getIndexOf(char const *name) const {
    return getIndexOf(std::string(name));
}
//@}

/** Returns the number of fields in the current record, or 0 if there is
  * no current record.
  */
inline int CsvReader::getNumFields() const {
    return static_cast<int>(_fields.size());
}

//@{
/** Returns true if the value of the given field is a database NULL.
  */
inline bool CsvReader::isNull(int i) const {
    if (_done) {
        _ioError("no current record (all records have been read)");
    }
    if (i < 0 || static_cast<size_t>(i) >= _fields.size()) {
        _ioError("attempt to access field with invalid field name or index");
    }
    return _fields[i] < 0;
}
inline bool CsvReader::isNull(std::string const &name) const {
    return isNull(getIndexOf(name));
}
inline bool CsvReader::isNull(char const *name) const {
    return isNull(getIndexOf(name));
}
//@}

//@{
/** Returns the value of a field as an instance of type T.
  */
template <typename T> inline T const CsvReader::get(int i) const {
    // force a compile time error for unsupported types
    static int const typeAllowed = boost::mpl::or_<
                                       boost::is_same<T, char const *>,
                                       boost::is_same<T, std::string>,
                                       boost::is_integral<T>,
                                       boost::is_floating_point<T>
                                   >::value;
    BOOST_STATIC_ASSERT(typeAllowed);
    if (isNull(i)) {
        return _null<T>();
    }
    // _get specializations live in the implementation file
    return _get<T>(_record.get() + _fields[i]);
}
template <typename T> inline T const CsvReader::get(std::string const &name) const {
    return get<T>(getIndexOf(name));
}
template <typename T> inline T const CsvReader::get(char const *name) const {
    return get<T>(getIndexOf(name));
}

/** Returns the value of a field as a std::string.
  */
inline std::string const CsvReader::get(int i) const {
    return get<std::string>(i);
}
inline std::string const CsvReader::get(std::string const &name) const {
    return get<std::string>(getIndexOf(name));
}
inline std::string const CsvReader::get(char const *name) const {
    return get<std::string>(getIndexOf(name));
}
//@}

// NULL values for various types
template <typename T> inline T CsvReader::_null() {
    BOOST_STATIC_ASSERT(sizeof(T) == 0);
}
template <> inline std::string CsvReader::_null<std::string>() {
    return std::string();
}
template <> inline char const * CsvReader::_null<char const *>() {
    return 0;
}
template <> inline bool CsvReader::_null<bool>() {
    return false;
}
template <> inline char CsvReader::_null<char>() {
    return '\0';
}

#define LSST_SPECIALIZE_NULL(T) \
    template <> inline T CsvReader::_null<T>() { \
        return std::numeric_limits<T>::min(); \
    }
LSST_SPECIALIZE_NULL(signed char)
LSST_SPECIALIZE_NULL(short)
LSST_SPECIALIZE_NULL(int)
LSST_SPECIALIZE_NULL(long)
LSST_SPECIALIZE_NULL(long long)
#undef LSST_SPECIALIZE_NULL

#define LSST_SPECIALIZE_NULL(T) \
    template <> inline T CsvReader::_null<T>() { \
        return std::numeric_limits<T>::max(); \
    }
LSST_SPECIALIZE_NULL(unsigned char)
LSST_SPECIALIZE_NULL(unsigned short)
LSST_SPECIALIZE_NULL(unsigned int)
LSST_SPECIALIZE_NULL(unsigned long)
LSST_SPECIALIZE_NULL(unsigned long long)
#undef LSST_SPECIALIZE_NULL

#define LSST_SPECIALIZE_NULL(T) \
    template <> inline T CsvReader::_null<T>() { \
        return std::numeric_limits<T>::quiet_NaN(); \
    }
LSST_SPECIALIZE_NULL(float)
LSST_SPECIALIZE_NULL(double)
LSST_SPECIALIZE_NULL(long double)
#undef LSST_SPECIALIZE_NULL

// Trivial _get implementations
template <>
inline char const * CsvReader::_get<char const *>(char const *field) const {
    return field;
}
template <>
inline std::string CsvReader::_get<std::string>(char const *field) const {
    return std::string(field);
}


// -- CsvWriter inlines ----

/** Forces a write of all user-space buffered data to the underlying stream.
  */
inline void CsvWriter::flush() {
    _out->flush();
}

/** Returns the number of lines written out.
  */
inline size_t CsvWriter::getNumLines() const {
    return _numLines;
}

/** Returns the number of records written out.
  */
inline size_t CsvWriter::getNumRecords() const {
    return _numRecords;
}

/** Returns the number of fields written to the current record.
  */
inline size_t CsvWriter::getNumFields() const {
    return _numFields;
}

inline void CsvWriter::appendField(std::string const &v) {
    _write(v.c_str());
}

inline void CsvWriter::appendField(char const *v) {
    if (v == 0) {
        appendNull();
    } else {
        _write(v);
    }
}

inline void CsvWriter::appendField(bool v) {
    _write(v ? "1" : "0");
}

inline void CsvWriter::appendField(char v) {
    char s[2];
    s[0] = v;
    s[1] = '\0';
    _write(s);
}

/** Raw write: the given string is written directly to the output stream
  * without being run through the quoting/escaping process. This can be useful
  * when repeatedly outputting the same records/fields.
  */
inline void CsvWriter::write(std::string const &s) {
    _out->write(s.c_str(), s.size());
}


// -- STL-style output for CsvWriter ----

/** Formatted output operator for CSV writers.
  */
template <typename T>
inline CsvWriter & operator<<(CsvWriter &w, T const &v) {
    w.appendField(v);
    return w;
}

/** Manipulator to end a CSV record.
  */
inline CsvWriter & endr(CsvWriter &w) {
    w.endRecord();
    return w;
}

/** Manipulator to flush a CSV writer.
  */
inline CsvWriter & flush(CsvWriter &w) {
    w.flush();
    return w;
}

/** Output operator that applies a manipulator to a CSV writer.
  */
template <typename T>
inline CsvWriter& operator<<(CsvWriter& w, CsvWriter& (*manip)(CsvWriter&)) {
    return manip(w);
}

}}} // namespace lsst::ap::util

#endif // LSST_AP_UTIL_CSV_CC
