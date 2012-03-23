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
  * @brief  Classes for CSV I/O.
  * @author Serge Monkewitz
  */
#ifndef LSST_AP_UTILS_CSV_H
#define LSST_AP_UTILS_CSV_H

#include <fstream>
#include <string>
#include <vector>

#include "lsst/tr1/unordered_map.h"

#include "boost/scoped_ptr.hpp"
#include "boost/scoped_array.hpp"

#include "lsst/pex/config.h"

#include "CsvControl.h"


namespace lsst { namespace ap { namespace utils {

/** A class for reading records from Character-Separated-Value files in
  * sequential order.
  *
  * @par Field Access
  *
  * Fields may be accessed either by name or by index - performance sensitive
  * code should favor the latter. Values can be retrieved as std::string
  * objects via @link get(int) const @endlink. Alternatively, to avoid
  * the cost of memory allocation, they can be retrieved via
  * <tt> get\<char const *\>(int) const</tt>, which returns a
  * null-terminated C string. But beware: the C string is located in the
  * readers internal decoded character buffer and is invalidated by the next
  * call to nextRecord()!
  *
  * @par
  * The <tt> get\<T\>(int) const </tt> method can be used to retrieve and
  * cast a field value simultaneously. The supported types are:
  * <dl>
  * <dt><b>bool</b></dt>
  * <dd>Values that are case insensitive matches to "y", "t", "yes", "true",
  *     or "1" are mapped to @c true. Case insensitive matches to "n", "f",
  *     "no", "false", "0" are mapped to @c false. Any other value results
  *     in an exception. Leading and trailing whitespace is permitted.
  * </dd>
  * <dt><b>char</b></dt>
  * <dd>Field values must contain exactly one character. Any other value
  *     results in an exception. Leading or trailing whitespace is illegal.
  * </dd>
  * <dt><b>integral types</b></dt>
  * <dd>If the field is a decimal representation of an integer that fits in
  *     the requested type, it is converted to binary value and returned.
  *     If the decimal value overflows the range of the type, cannot be
  *     converted, or contains extraneous characters other than leading or
  *     trailing whitespace, an exception is thrown.
  * </dd>
  * <dt><b>floating point types</b></dt>
  * <dd>If the field value is a floating point number, it is converted to 
  *     a binary value and returned. If the value overflows the range of the
  *     type, +/-INF is returned. If it underflows, 0.0 is returned. If the
  *     field cannot be converted or contains extraneous characters other
  *     than leading or trailing whitespace, an exception is thrown.
  * </dd>
  * </dl>
  *
  * @par NULL values
  *
  * When a field is recognized as a database NULL (e.g. \\N with the
  * CsvControl::standardEscapes == true), the field access methods return
  * the following values:
  *
  * <dl>
  * <dt><b> std::string </b></dt>  <dd> An empty string </dd>
  * <dt><b> char const * </b></dt> <dd> A null pointer </dd>
  * <dt><b> bool </b></dt>         <dd> @c false </dd>
  * <dt><b> char </b></dt>         <dd> '\\0' </dd>
  *
  * <dt><b> signed integral types   </b></dt><dd> The minimum representable value </dd>
  * <dt><b> unsigned integral types </b></dt><dd> The maximum representable value </dd>
  * <dt><b> floating point types    </b></dt><dd> A quiet NaN </dd>
  * </dl>
  *
  * To distinguish between a NULL field and one which happens to take one of
  * values above, use the @link isNull(int) const @endlink method.
  *
  * @par Exception safety
  *
  * All methods provide the strong exception safety guarantee, except for
  * nextRecord() - see the documentation of that method for details.
  *
  * @par Limitations
  *
  * @li CSV files containing embedded null characters are not read in properly.
  * @li The line terminator cannot be specified, and is always '\\n'.
  */
class CsvReader {
public:
    CsvReader(std::string const &path,
              CsvControl const &control,
              bool namesInFirstRecord=false);
    CsvReader(std::istream &in,
              CsvControl const &control,
              bool namesInFirstRecord=false);
    ~CsvReader();

    inline CsvControl const & getControl() const;

    // Return the number of lines/records read
    inline size_t getNumLines() const;
    inline size_t getNumRecords() const;

    // Get/set field names
    inline std::vector<std::string> const & getFieldNames() const;
    void setFieldNames(std::vector<std::string> const & names);
    void setFieldNames(std::string const & names,
                       std::string const & regex,
                       bool stripWhitespace=true);

    // Map field names to field indexes
    inline int getIndexOf(std::string const &name) const;
    inline int getIndexOf(char const *name) const;

    // Have all records been read?
    inline bool isDone() const;

    // Advance to the next record
    inline void nextRecord();

    ///@name Access fields in the current record
    //@{ 
    inline int getNumFields() const;

    inline bool isNull(int i) const;
    inline bool isNull(std::string const &name) const;
    inline bool isNull(char const *name) const;

    inline std::string const get(int i) const;
    inline std::string const get(std::string const &name) const;
    inline std::string const get(char const *name) const;

    // Get and convert a field value. Allowed types are std::string, bool,
    // built-in integral and floating point types, and char const *. The
    // pointers returned by get<char const *>() are invalidated by a call
    // to nextRecord().
    template <typename T> inline T const get(int i) const;
    template <typename T> inline T const get(std::string const &name) const;
    template <typename T> inline T const get(char const *name) const;
    //@}

private:
    typedef std::tr1::unordered_map<std::string, int> FieldIndexes;

    static std::string const WHITESPACE;
    static int const MAX_RECORD_LENGTH;
    static int const DEFAULT_CAPACITY;

    enum State {
        START_RECORD = 0,
        START_FIELD,
        IN_FIELD,
        IN_QUOTED_FIELD,
        INITIAL_ESCAPE,
        ESCAPE,
        EMBEDDED_QUOTE
    };

    // disable copy construction/assignment
    CsvReader(CsvReader const &);
    CsvReader & operator=(CsvReader const &);

    void _ioError(char const *msg) const;
    void _runtimeError(char const *msg) const;
    bool _readLine(int offset);
    void _readRecord();
    void _checkWhitespace(char const *s, char const *msg) const;
    template <typename T> static inline T _null();
    template <typename T> T _get(char const *value) const; 

    std::string _path;                  ///< File name.
    CsvControl _control;                ///< File format.
    std::vector<std::string> _names;    ///< Field names in order of occurence.
    FieldIndexes _indexes;              ///< Field name to index map.

    boost::scoped_ptr<std::ifstream> _stream; ///< File stream.
    std::istream *_in;                  ///< Input stream.
    size_t _numLines;                   ///< 1-based index of current line.
    size_t _numRecords;                 ///< 1-based index of current record.
    boost::scoped_array<char> _record;  ///< Data for a single record.
    int _capacity;                      ///< Capacity of _record.
    bool _done;                         ///< Finished reading file?
    std::vector<int> _fields;           ///< Index of first character for each
                                        ///  field in _record, -1 if NULL.
};

/** A class for writing fields and records to Character-Separated-Value
  * files in sequential order.
  *
  * @par Field Output
  *
  * Methods are provided to output strings, bool, signed/unsigned integers,
  * floating point types, and database NULLs. Integers are always output in
  * base 10, and booleans are output as "1" (true) or "0". The precision of
  * the decimal output for floating point numbers is such that the original
  * binary value can be recovered exactly if the systems float to decimal
  * conversion routines are correctly rounded and the input/output rounding
  * modes agree.
  *
  * @par Exception Safety
  *
  * This class provides the minimal exception safety guarantee - that is,
  * exceptions will not cause crashes or leaks, and one can technically
  * continue using the writer. However, the CSV files written under these
  * conditions are not guaranteed to be legal. For example, partially written
  * fields may be produced by a writer when raising an exception.
  *
  * @par Limitations
  *
  * @li Field values containing embedded null characters are not written out
  *     properly. Only the characters before the first null character will
  *     be written.
  * @li The line terminator cannot be specified, and is always '\\n'.
  */
class CsvWriter {
public:
    CsvWriter(std::string const &path,
              CsvControl const &control,
              bool truncate=false);
    CsvWriter(std::ostream &out,
              CsvControl const &control);
    ~CsvWriter();

    inline CsvControl const & getControl() const;

    // End the current record
    void endRecord();

    // Flush the underlying output stream without ending the current record
    inline void flush();

    // Return the number of lines/records/fields written
    inline size_t getNumLines() const;
    inline size_t getNumRecords() const;
    inline size_t getNumFields() const;

    /// @name Append field(s) to the end of a record
    //@{ 
    void appendFields(std::vector<std::string> const &fields);
    inline void appendField(std::string const &v);
    inline void appendField(char const *v);
    inline void appendField(bool v);
    inline void appendField(char v);
    void appendField(signed char v);
    void appendField(short v);
    void appendField(int v);
    void appendField(long v);
    void appendField(long long v);
    void appendField(unsigned char v);
    void appendField(unsigned short v);
    void appendField(unsigned int v);
    void appendField(unsigned long v);
    void appendField(unsigned long long v);
    void appendField(float v);
    void appendField(double v);
    void appendField(long double v);

    void appendNull();
    //@}

    // Raw writes. Beware: these bypass all output formatting!
    inline void write(std::string const &v);
    inline void write(char c);

private:
    // disable copy construction/assignment
    CsvWriter(CsvWriter const &);
    CsvWriter & operator=(CsvWriter const &);

    void _write(char const *s);
    void _writeQuoted(char const *s);
    void _writeUnquoted(char const *s);

    boost::scoped_ptr<std::ofstream> _stream; ///< Output file stream.
    std::ostream *_out;                       ///< Output stream.
    CsvControl _control;                      ///< File format.
    size_t _numRecords;                       ///< Number of records written.
    size_t _numLines;                         ///< Number of lines written.
    size_t _numFields;                        ///< Number of fields written.
};

// STL-style output for CsvWriter
inline CsvWriter & endr(CsvWriter &);
inline CsvWriter & flush(CsvWriter &);
inline CsvWriter & nullf(CsvWriter &);
template <typename T> inline CsvWriter & operator<<(CsvWriter &, T const &);
template <typename T> inline CsvWriter& operator<<(CsvWriter&, CsvWriter& (*)(CsvWriter&));

}}} // namespace lsst::ap::utils

#include "Csv.cc" // for inline/template members

#endif // LSST_AP_UTILS_CSV_H

