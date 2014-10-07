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
  * @brief  CSV I/O class implementations.
  * @author Serge Monkewitz
  */
#include "lsst/ap/utils/Csv.h"

#include <cerrno>
#include <cfloat>
#include <cstring>
#include <sstream>
#include <utility>

#include "boost/filesystem.hpp"
#include "boost/regex.hpp"

#include "lsst/utils/ieee.h"
#include "lsst/pex/exceptions.h"

using std::ifstream;
using std::ios;
using std::make_pair;
using std::memcpy;
using std::min;
using std::numeric_limits;
using std::ofstream;
using std::ostringstream;
using std::pair;
using std::string;
using std::strtof;
using std::strtod;
using std::strtold;
using std::strtol;
using std::strtoll;
using std::strtoul;
using std::strtoull;
using std::swap;
using std::vector;

namespace pexExcept = lsst::pex::exceptions;


namespace lsst { namespace ap { namespace utils {

// -- CsvReader implementation ----

std::string const CsvReader::WHITESPACE(" \t\f\v\n\r");
int const CsvReader::MAX_RECORD_LENGTH = 4*1024*1024;
int const CsvReader::DEFAULT_CAPACITY = 128*1024;

/** Creates a new CsvReader for a file and reads the first record.
  */
CsvReader::CsvReader(
    std::string const &path,   ///< Input file name.
    CsvControl const &control, ///< CSV format of the input file.
    bool namesInFirstRecord    ///< Set field names to the strings in the
                               ///  first record of the input file?
   ) :
    _path(path),
    _control(control),
    _names(),
    _indexes(),
    _stream(new ifstream(path.c_str(), ios::binary | ios::in)),
    _in(_stream.get()),
    _numLines(0),
    _numRecords(0),
    _record(new char[DEFAULT_CAPACITY]),
    _capacity(DEFAULT_CAPACITY),
    _done(false),
    _fields()
{
    if (!_stream->good()) {
        throw LSST_EXCEPT(pexExcept::IoError,
                          "failed to open file " + path + " for reading");
    }
    // exception mask for _stream is clear
    _readRecord();
    if (namesInFirstRecord && !isDone()) {
        vector<string> names;
        for (int i = 0; i < getNumFields(); ++i) {
            names.push_back(get(i));
        }
        setFieldNames(names);
        _readRecord();
    }
}

/** Creates a new CsvReader from an std::istream and reads the first record.
  *
  * This std::istream must be kept alive for the life-time of the CsvReader.
  * If it is externally modified while the reader is alive (e.g. its exception
  * mask is changed, or external reads are performed from it) then the
  * behaviour of the reader is undefined.
  */
CsvReader::CsvReader(
    std::istream &in,           ///< Input stream
    CsvControl const &control,  ///< CSV format of the input stream
    bool namesInFirstRecord     ///< Set field names to the strings in the
                                ///  first record read from the input stream?
   ) :
    _path(),
    _control(control),
    _names(),
    _indexes(),
    _stream(),
    _in(&in),
    _numLines(0),
    _numRecords(0),
    _record(new char[DEFAULT_CAPACITY]),
    _capacity(DEFAULT_CAPACITY),
    _done(false),
    _fields()
{
    if (!in.good()) {
        throw LSST_EXCEPT(pexExcept::IoError,
                          "std::istream not good() for reading");
    }
    // clear exception mask
    in.exceptions(ios::goodbit);
    _readRecord();
    if (namesInFirstRecord && !isDone()) {
        vector<string> names;
        for (int i = 0; i < getNumFields(); ++i) {
            names.push_back(get(i));
        }
        setFieldNames(names);
        _readRecord();
    }
}

CsvReader::~CsvReader() {
    _in = 0;
}

/** Associates field names with field indexes. The i-th field name in the
  * list corresponds to the i-th field in a record. Once field names have
  * been set, they can be used in lieu of field indexes to look up field
  * values.
  *
  * @par
  * Field names are case-sensitive and must all be distinct. The only
  * exception is that empty field names are skipped.
  */
void CsvReader::setFieldNames(std::vector<std::string> const &names) {
    if (names.size() > static_cast<size_t>(numeric_limits<int>::max())) {
        throw LSST_EXCEPT(pexExcept::InvalidParameterError,
                          "too many field names");
    }
    // create temporary name->index map.
    FieldIndexes indexes;
    for (vector<string>::size_type i = 0; i < names.size(); ++i) {
        if (names[i].empty()) {
            continue; // skip empty field names
        }
        pair<FieldIndexes::iterator, bool> p = indexes.insert(make_pair(
            names[i], static_cast<int>(i)));
        if (!p.second) {
            throw LSST_EXCEPT(pexExcept::InvalidParameterError,
                              "Duplicate field name: " + names[i]);
        }
    }
    // copy name list
    vector<string> copy(names);
    // commit state
    swap(copy, _names);
    swap(indexes, _indexes);
}

/** Splits a string containing a delimited list of field names using the given
  * regular expression and associates the resulting field names with field
  * indexes.
  *
  * @par
  * Field names are case-sensitive and must all be distinct. The only
  * exception is that empty field names are skipped.
  *
  * @sa setFieldNames(std::vector<std::string> const&)
  */
void CsvReader::setFieldNames(
    std::string const &names, ///< Delimited list of field names.
    std::string const &regex, ///< Regular expression matching delimiter.
    bool stripWhitespace      ///< Strip whitespace from field names?
) {
    vector<string> fieldNames;
    boost::regex re(regex);
    boost::sregex_token_iterator i(names.begin(), names.end(), re, -1);
    boost::sregex_token_iterator e;
    for (; i != e; ++i) {
        string name = i->str();
        if (stripWhitespace) {
            size_t w = name.find_last_not_of(WHITESPACE);
            if (w != string::npos) {
                name.erase(w + 1);
            } else {
                name.clear();
            }
            w = name.find_first_not_of(WHITESPACE);
            if (w != string::npos) {
                name.erase(0, w);
            } else {
                name.clear();
            }
        }
        fieldNames.push_back(name);
    }
    setFieldNames(fieldNames);
}

/** Throws an lsst::pex::exceptions::IoError with a file name,
  * line number and record number.
  */
void CsvReader::_ioError(char const *msg) const {
    ostringstream s;
    if (_path.empty()) {
        s << "CSV stream";
    } else {
        s << "CSV file " << _path;
    }
    s << " line " << _numLines << " record " << _numRecords << ": " << msg;
    throw LSST_EXCEPT(pexExcept::IoError, s.str());
}

/** Throws an lsst::pex::exceptions::RuntimeError with a file name,
  * line number and record number.
  */
void CsvReader::_runtimeError(char const *msg) const {
    ostringstream s;
    if (_path.empty()) {
        s << "CSV stream";
    } else {
        s << "CSV file " << _path;
    }
    s << " line " << _numLines << " record " << _numRecords << ": " << msg;
    throw LSST_EXCEPT(pexExcept::RuntimeError, s.str());
}

/** Reads a single line from the underlying stream.
  */
bool CsvReader::_readLine(int offset) {
    bool gotData = false;
    while (true) {
        _in->getline(_record.get() + offset, _capacity - 1 - offset);
        if (_in->eof()) {
            if (_in->gcount() > 0) {
                gotData = true;
                ++_numLines;
            }
            return gotData;
        } else if (_in->bad()) {
            _ioError("std::istream::getline() failed");
        } else if (_in->fail()) {
            // not enough space for line in buffer - expand it
            offset += _in->gcount();
            if (_capacity == MAX_RECORD_LENGTH) {
                _ioError("record too long");
            }
            int cap = min(_capacity*2, MAX_RECORD_LENGTH);
            boost::scoped_array<char> rec(new char[cap]);
            memcpy(rec.get(), _record.get(), static_cast<size_t>(_capacity));
            swap(rec, _record);
            _capacity = cap;
            continue;
        }
        offset += _in->gcount();
        ++_numLines;
        return true;
    }
}

/** Reads a single CSV record from the underlying stream.
  */
void CsvReader::_readRecord() {
    State state = START_RECORD;
    State popState = START_RECORD;
    int writeOffset = 0;
    int offset = 0;
    char c = 0, c2 = 0;

    _fields.clear();
    if (!_readLine(offset)) {
        _done = true;
        return;
    }
    ++_numRecords;

    // Run parsing state machine until we reach the end of a record,
    // then return.
    while (true) {
        c = _record[offset++];
        switch (state) {
            case START_RECORD:
                if (c == '\0') {
                    if (_control.trailingDelimiter) {
                        _runtimeError("expecting trailing delimiter "
                                      "at end of record");
                    } else if (_control.hasNull &&
                               _control.null.empty()) {
                       _fields.push_back(-1);
                    }
                    return; // finished record
                }
                state = START_FIELD;
                // Fall-through

            case START_FIELD:
                _fields.push_back(writeOffset);
                if (c == '\0') {
                    // finished empty field
                    _record[writeOffset++] = '\0';
                    if (_control.hasNull && _control.null.empty()) {
                        _fields.back() = -1; // NULL
                    }
                    if (_control.trailingDelimiter) {
                        // a trailing delimiter is expected and does not
                        // yield a trailing empty field
                        _fields.pop_back();
                    }
                    return; // finished record
                } else if (_control.quoting != "QUOTE_NONE" &&
                           c == _control.getQuoteChar()) {
                    // start quoted field
                    state = IN_QUOTED_FIELD;
                } else if (c == _control.getEscapeChar() &&
                           _control.getEscapeChar() != '\0') {
                    // escape at beginning of field
                    popState = IN_FIELD;
                    state = INITIAL_ESCAPE;
                } else if (c == _control.getDelimiter()) {
                    // finished empty field
                    _record[writeOffset++] = '\0';
                    if (_control.hasNull && _control.null.empty()) {
                        _fields.back() = -1; // NULL
                    }
                    state = START_FIELD;
                } else if (_control.skipInitialSpace &&
                           WHITESPACE.find(c) != string::npos) {
                    // eat initial whitespace
                    state = START_FIELD;
                } else {
                    // start unquoted field
                    _record[writeOffset++] = c;
                    state = IN_FIELD;
                }
                break;

            case IN_FIELD:
                if (c == '\0') {
                    // finished field
                    _record[writeOffset++] = '\0';
                    if (_control.trailingDelimiter) {
                        _runtimeError("expecting trailing delimiter "
                                      "at end of record");
                    }
                    if (_control.hasNull &&
                        _control.null == _record.get() + _fields.back()) {
                        _fields.back() = -1; // NULL
                    }
                    return; // finished record
                } else if (c == _control.getEscapeChar()) {
                    popState = IN_FIELD;
                    state = ESCAPE;
                } else if (c == _control.getDelimiter()) {
                    // finished field
                    _record[writeOffset++] = '\0';
                    if (_control.hasNull &&
                        _control.null == _record.get() + _fields.back()) {
                        _fields.back() = -1; // NULL
                    }
                    state = START_FIELD;
                } else {
                    _record[writeOffset++] = c;
                    state = IN_FIELD;
                }
                break;

            case IN_QUOTED_FIELD:
                if (c == '\0') {
                    // newline in quoted field
                    _record[writeOffset++] = '\n';
                    // keep reading field on next line
                    if (!_readLine(offset)) {
                        _record[writeOffset++] = '\0';
                        _runtimeError("expecting quote character "
                                      "at end of field");
                    }
                    state = IN_QUOTED_FIELD;
                } else if (c == '\0') {
                    _record[writeOffset++] = '\0';
                    _runtimeError("expecting quote character at end of field");
                } else if (c == _control.getEscapeChar()) {
                    popState = IN_QUOTED_FIELD;
                    state = ESCAPE;
                } else if (c == _control.getQuoteChar()) {
                    state = EMBEDDED_QUOTE;
                } else {
                    _record[writeOffset++] = c;
                    state = IN_QUOTED_FIELD;
                }
                break;

            case INITIAL_ESCAPE:
                if (c == 'N' && _control.standardEscapes) {
                    char c2 = _record[offset];
                    if (c2 == '\0') {
                        // finished NULL field
                        _fields.back() = -1;
                        if (_control.trailingDelimiter) {
                            _runtimeError("expecting trailing delimiter "
                                          "at end of record");
                        }
                        // finished record
                        return;
                    } else if (c2 == _control.getDelimiter()) {
                        // finished NULL field
                        _fields.back() = -1;
                        ++offset;
                        state = START_FIELD;
                        break;
                    }
                }
                // Fall-through

            case ESCAPE:
                if (c == '\0') {
                    _record[writeOffset++] = '\n';
                    // keep reading field on next line
                    if (!_readLine(offset)) {
                        _record[writeOffset++] = '\0';
                    }
                }
                if (_control.standardEscapes) {
                    // handle standard escape sequences
                    switch (c) {
                        case 'Z': c = 0x1a; break;
                        case 'b': c = 0x08; break;
                        case 'n': c = '\n'; break;
                        case 'r': c = '\r'; break;
                        case 't': c = '\t'; break;
                        case 'v': c = 0x0b; break;
                        case 'x':
                            // one or two digit hex-escape sequence
                            c2 = _record[offset++];
                            if (c2 >= '0' && c2 <= '9') {
                                c = c2 - '0';
                            } else if (c2 >= 'A' && c2 <= 'F') {
                                c = 10 + (c2 - 'A');
                            } else if (c2 >= 'a' && c2 <= 'f') {
                                c = 10 + (c2 - 'a');
                            } else {
                                _record[writeOffset++] = '\0';
                                _runtimeError("hex escape sequence must be "
                                              "followed by one or two hex "
                                              "digits");
                            }
                            c2 = _record[offset];
                            if (c2 >= '0' && c2 <= '9') {
                                c = c*16 + (c2 - '0');
                                ++offset;
                            } else if (c2 >= 'A' && c2 <= 'F') {
                                c = c*16 + 10 + (c2 - 'A');
                                ++offset;
                            } else if (c2 >= 'a' && c2 <= 'f') {
                                c = c*16 + 10 + (c2 - 'a');
                                ++offset;
                            }
                            break;
                        default:
                            break;
                    }
                }
                _record[writeOffset++] = c;
                state = popState;
                break;

            case EMBEDDED_QUOTE:
                if (c == '\0') {
                    // finished field
                    _record[writeOffset++] = '\0';
                    if (_control.trailingDelimiter) {
                        _runtimeError("expecting trailing delimiter "
                                      "at end of record");
                    }
                    return; // finished record
                } else if (c == _control.getQuoteChar() &&
                           _control.doubleQuote) {
                    // save "" as "
                    _record[writeOffset++] = c;
                    state = IN_QUOTED_FIELD;
                } else if (c == _control.getDelimiter()) {
                    // finished field
                    _record[writeOffset++] = '\0';
                    state = START_FIELD;
                } else {
                    _record[writeOffset++] = '\0';
                    _runtimeError("expecting delimiter after ending quote");
                }
                break;

            default:
                throw LSST_EXCEPT(pexExcept::LogicError,
                                  "CSV parser bug - state machine reached an "
                                  "illegal state");
        }
    }
}

/** Throws an exception if the given string contains anything but
  * whitespace.
  */
void CsvReader::_checkWhitespace(char const *s, char const *msg) const {
    for (; WHITESPACE.find(*s) != string::npos; ++s) { }
    if (*s != '\0') {
        _runtimeError(msg);
    }
}

template <> bool CsvReader::_get<bool>(char const *field) const {
    for (; WHITESPACE.find(*field) != string::npos; ++field) { }
    bool value = false;
    char c = *field++;
    switch (c) {
        case '1':
            value = true;
            break;
        case 't':
        case 'T':
            value = true;
            if ((field[0] == 'r' || field[0] == 'R') &&
                (field[1] == 'u' || field[1] == 'U') &&
                (field[2] == 'e' || field[2] == 'E')) {
                field += 3;
            }
            break;
        case 'y':
        case 'Y':
            value = true;
            if ((field[0] == 'e' || field[0] == 'E') &&
                (field[1] == 's' || field[1] == 'S')) {
                field += 2;
            }
            break;
        case '0':
            break;
        case 'f':
        case 'F':
            if ((field[0] == 'a' || field[0] == 'A') &&
                (field[1] == 'l' || field[1] == 'L') &&
                (field[2] == 's' || field[2] == 'S') &&
                (field[3] == 'e' || field[3] == 'E')) {
                field += 4;
            }
            break;
        case 'n':
        case 'N':
            if (field[0] == 'o' || field[0] == 'O') {
                ++field;
            }
            break;
        default:
            _runtimeError("failed to convert field value to bool");
    }
    _checkWhitespace(field, "failed to convert field value to bool");
    return value;
}

template <> char CsvReader::_get<char>(char const *field) const {
    if (field[0] == '\0') {
        _runtimeError("empty field");
    }
    if (field[1] != '\0') {
        _runtimeError("field value contains more than one character");
    }
    return field[0];
}

#define LSST_SPECIALIZE_GET(T, C, fun) \
    template <> T CsvReader::_get<T>(char const *field) const { \
        char *e; \
        errno = 0; \
        C v = fun(field, &e, 10); \
        if (e == field) { \
            _runtimeError("failed to convert field value to " #T); \
        } else if ((errno == ERANGE) || \
                   v > numeric_limits<T>::max() || \
                   v < numeric_limits<T>::min()) { \
            _runtimeError("field value overflow during conversion to " #T); \
        } \
        _checkWhitespace(e, "failed to convert field value to " #T); \
        return v; \
    }

    LSST_SPECIALIZE_GET(signed char, long, strtol)
    LSST_SPECIALIZE_GET(short, long, strtol)
    LSST_SPECIALIZE_GET(int, long, strtol)
    LSST_SPECIALIZE_GET(long, long, strtol)
    LSST_SPECIALIZE_GET(long long, long long, strtoll)
    LSST_SPECIALIZE_GET(unsigned char, unsigned long, strtoul)
    LSST_SPECIALIZE_GET(unsigned short, unsigned long, strtoul)
    LSST_SPECIALIZE_GET(unsigned int, unsigned long, strtoul)
    LSST_SPECIALIZE_GET(unsigned long, unsigned long, strtoul)
    LSST_SPECIALIZE_GET(unsigned long long, unsigned long long, strtoull)
#undef LSST_SPECIALIZE_GET

#define LSST_SPECIALIZE_GET(T, fun) \
    template <> T CsvReader::_get<T>(char const *field) const { \
        char *e; \
        errno = 0; \
        T v = fun(field, &e); \
        if (e == field) { \
            _runtimeError("failed to convert field value to " #T); \
        } \
        _checkWhitespace(e, "failed to convert field value to " #T); \
        return v; \
    }

    LSST_SPECIALIZE_GET(float, strtof)
    LSST_SPECIALIZE_GET(double, strtod)
    LSST_SPECIALIZE_GET(long double, strtold)
#undef LSST_SPECIALIZE_GET


// -- CsvWriter implementation ----

/** Creates a new CsvWriter that outputs to the given file.
  */
CsvWriter::CsvWriter(
    std::string const &path,   ///< Name of file to write to.
    CsvControl const &control, ///< CSV format of the output file.
    bool truncate,             ///< If true, an existing file is truncated.
                               ///  Otherwise, an attempt to create a writer
                               ///  for an existing file raises an exception
                               ///  unless append is true.
    bool append                ///  If true, append to an existing file.
                               ///  Has no effect if truncate is true.
) :
    _stream(),
    _out(0),
    _control(control),
    _numRecords(0),
    _numLines(0),
    _numFields(0)
{
    ios::openmode mode = ios::out | ios::binary;
    if (!truncate && boost::filesystem::exists(path)) {
        if (append) {
            mode |= ios::app;
        } else {
            throw LSST_EXCEPT(pexExcept::IoError,
                              "file " + path + " already exists");
        }
    }
    if (truncate) {
        mode |= ios::trunc;
    }
    _stream.reset(new ofstream(path.c_str(), mode));
    if (!_stream->good()) {
        throw LSST_EXCEPT(pexExcept::IoError,
                          "failed to open file " + path + " for writing");
    }
    _out = _stream.get();
    // throw on any kind of output error
    _stream->exceptions(ios::eofbit | ios::failbit | ios::badbit);
}

/** Creates a new CsvWriter that sends output to the given std::ostream.
  *
  * This std::ostream must be kept alive for the life-time of the CsvWriter.
  * If it is externally modified while the writer is alive (e.g. its exception
  * mask is changed, or external writes are performed on it), then the
  * behaviour and output of the writer is undefined.
  */
CsvWriter::CsvWriter(std::ostream &out, CsvControl const &control) :
    _stream(),
    _out(&out),
    _control(control),
    _numRecords(0),
    _numLines(0),
    _numFields(0)
{
    if (!out.good()) {
        throw LSST_EXCEPT(pexExcept::IoError,
                          "std::ostream not good() for writing");
    }
    out.exceptions(ios::eofbit | ios::failbit | ios::badbit);
}

CsvWriter::~CsvWriter() {
    _out = 0;
}

/** Terminates a record - subsequent apppends will be to a new record.
  */
void CsvWriter::endRecord() {
    if (_control.trailingDelimiter) {
        _out->put(_control.getDelimiter());
    }
    _out->put('\n');
    _numFields = 0;
    ++_numLines;
    ++_numRecords;
}

/** Appends a field for each string in the given list.
  */
void CsvWriter::appendFields(std::vector<std::string> const &fields) {
    for (vector<string>::const_iterator i = fields.begin(), e = fields.end();
         i != e; ++i) {
        appendField(*i);
    }
}

#define LSST_IMPLEMENT_APPEND_FIELD(T, C, fmt) \
    void CsvWriter::appendField(T v) { \
        char buf[64]; \
        int n = snprintf(buf, sizeof(buf), fmt, static_cast<C>(v)); \
        if (n <= 0) { \
            throw LSST_EXCEPT(pexExcept::RuntimeError, \
                              "failed to convert " #T " to a string"); \
        } else if (n >= static_cast<int>(sizeof(buf))) { \
            throw LSST_EXCEPT(pexExcept::LogicError, \
                              "internal buffer for string conversion too small"); \
        } \
        _write(buf); \
    }

    LSST_IMPLEMENT_APPEND_FIELD(short, short, "%hd")
    LSST_IMPLEMENT_APPEND_FIELD(int, int, "%d")
    LSST_IMPLEMENT_APPEND_FIELD(long, long, "%ld")
    LSST_IMPLEMENT_APPEND_FIELD(long long, long long, "%lld")
    LSST_IMPLEMENT_APPEND_FIELD(unsigned short, unsigned short, "%hu")
    LSST_IMPLEMENT_APPEND_FIELD(unsigned int, unsigned int, "%u")
    LSST_IMPLEMENT_APPEND_FIELD(unsigned long, unsigned long, "%lu")
    LSST_IMPLEMENT_APPEND_FIELD(unsigned long long, unsigned long long, "%llu")
    // Some platforms don't support %hhd/%hhu in snprintf
    LSST_IMPLEMENT_APPEND_FIELD(signed char, int, "%d")
    LSST_IMPLEMENT_APPEND_FIELD(unsigned char, unsigned int, "%u")
#undef LSST_IMPLEMENT_APPEND_FIELD


// For a floating point type F, gymnastics are required to determine the
// number of decimal digits to print such that F -> decimal -> F conversions
// are guaranteed lossless. Note that this guarantee can only be made if the
// platforms conversion functions (snprintf, strtof, strtod, strtold) are
// correctly rounded, and the same rounding mode is used on input and output.

void CsvWriter::appendField(float v) {
    char buf[64];
    char fmt[32];
    int n = 0;
    if (_control.nonfiniteAsNull && !lsst::utils::isfinite(v)) {
        appendNull();
        return;
    }
#if FLT_RADIX == 2
    if (FLT_MANT_DIG == 24) {
        // 32 bit IEEE fast path
        n = snprintf(buf, sizeof(buf), "%.9g", static_cast<double>(v));
    } else {
        // ceil(1 + FLT_MANT_DIG*log10(2))
        static unsigned long const ndig = 2 + (FLT_MANT_DIG*30103UL)/100000UL;
#elif defined(FLT_MAXDIG10)
        static unsigned long const ndig = FLT_MAXDIG10;
#elif defined(DECIMAL_DIG)
        static unsigned long const ndig = DECIMAL_DIG;
#else
#   error Unable to determine number of digits for lossless float->decimal->float conversion
#endif
        n = snprintf(fmt, sizeof(fmt), "%%.%lug", ndig);
        if (n <= 0 || n >= static_cast<int>(sizeof(fmt))) {
            throw LSST_EXCEPT(pexExcept::LogicError, \
                              "internal buffer for string conversion too small"); \
        }
        n = snprintf(buf, sizeof(buf), fmt, static_cast<double>(v));
#if FLT_RADIX == 2
    }
#endif
    if (n <= 0 || n >= static_cast<int>(sizeof(buf))) {
        throw LSST_EXCEPT(pexExcept::LogicError,
                          "snprintf() failed to convert float to string");
    }
    _write(buf);
}

void CsvWriter::appendField(double v) {
    char buf[64];
    char fmt[32];
    int n = 0;
    if (_control.nonfiniteAsNull && !lsst::utils::isfinite(v)) {
        appendNull();
        return;
    }
#if FLT_RADIX == 2
    if (DBL_MANT_DIG == 53) {
        // 64 bit IEEE fast path
        n = snprintf(buf, sizeof(buf), "%.17g", v);
    } else {
        // ceil(1 + DBL_MANT_DIG*log10(2))
        static unsigned long const ndig = 2 + (DBL_MANT_DIG*30103UL)/100000UL;
#elif defined(DBL_MAXDIG10)
        static unsigned long const ndig = DBL_MAXDIG10;
#elif defined(DECIMAL_DIG)
        static unsigned long const ndig = DECIMAL_DIG;
#else
#   error Unable to determine number of digits for lossless double->decimal->double conversion
#endif
        n = snprintf(fmt, sizeof(fmt), "%%.%lug", ndig);
        if (n <= 0 || n >= static_cast<int>(sizeof(fmt))) {
            throw LSST_EXCEPT(pexExcept::LogicError,
                              "snprintf() failed to produce format string");
        }
        n = snprintf(buf, sizeof(buf), fmt, v);
#if FLT_RADIX == 2
    }
#endif
    if (n <= 0 || n >= static_cast<int>(sizeof(buf))) {
        throw LSST_EXCEPT(pexExcept::LogicError,
                          "snprintf() failed to convert double to string");
    }
    _write(buf);
}

void CsvWriter::appendField(long double v) {
    char buf[64];
    char fmt[32];
    int n = 0;
    if (_control.nonfiniteAsNull && !lsst::utils::isfinite(v)) {
        appendNull();
        return;
    }
#if FLT_RADIX == 2
    if (LDBL_MANT_DIG == 64) {
        // 80bit IEEE long double fast path
        n = snprintf(buf, sizeof(buf), "%.21Lg", v);
    } else if (LDBL_MANT_DIG == 113) {
        // 128bit IEEE long double fast path
        n = snprintf(buf, sizeof(buf), "%.36Lg", v);
    } else {
        // ceil(1 + LDBL_MANT_DIG*log10(2))
        static unsigned long const ndig = 2 + (LDBL_MANT_DIG*30103UL)/100000UL;
#elif defined(DBL_MAXDIG10)
        static unsigned long const ndig = LDBL_MAXDIG10;
#elif defined(DECIMAL_DIG)
        static unsigned long const ndig = DECIMAL_DIG;
#else
#   error Unable to determine number of digits for lossless long double->decimal->long double conversion
#endif
        n = snprintf(fmt, sizeof(fmt), "%%.%luLg", ndig);
        if (n <= 0 || n >= static_cast<int>(sizeof(fmt))) {
            throw LSST_EXCEPT(pexExcept::LogicError,
                              "snprintf() failed to produce format string");
        }
        n = snprintf(buf, sizeof(buf), fmt, v);
#if FLT_RADIX == 2
    }
#endif
    if (n <= 0 || n >= static_cast<int>(sizeof(buf))) {
        throw LSST_EXCEPT(pexExcept::LogicError,
                          "snprintf() failed to convert long double to string");
    }
    _write(buf);
}

/** Appends a database NULL field to the current record. If there is no
  * NULL representation, an empty field is written.
  */
void CsvWriter::appendNull() {
    if (_control.hasNull) {
        // output leading delimiter except for the first field in a record
        if (_numFields > 0) {
            _out->put(_control.getDelimiter());
        }
        ++_numFields;
        // NULL is never quoted, and guaranteed not to require escaping
        _out->write(_control.null.c_str(), _control.null.size());
    } else if (_control.standardEscapes) {
        // output leading delimiter except for the first field in a record
        if (_numFields > 0) {
            _out->put(_control.getDelimiter());
        }
        ++_numFields;
        // write \N
        _out->put(_control.getEscapeChar());
        _out->put('N');
    } else {
        // write an empty field
        _write("");
    }
}

/** Appends a field to the current record.
  */
void CsvWriter::_write(char const *s) {
    // output leading delimiter except for the first field in a record
    if (_numFields > 0) {
        _out->put(_control.getDelimiter());
    }
    ++_numFields;
    if (_control.quoting == "QUOTE_NONE") {
        if (_control.hasNull && _control.null == s) {
            throw LSST_EXCEPT(pexExcept::RuntimeError,
                              "Field value coincides with NULL string "
                              "and quoting is disabled");
        }
        _writeUnquoted(s);
        return;
    } else if (_control.quoting == "QUOTE_ALL" ||
               (_control.hasNull && _control.null == s)) {
        _writeQuoted(s);
        return;
    }

    // minimal quoting mode - first determine whether to quote/escape
    size_t n = 0;
    bool wantEscape = false;
    bool wantQuote = false;
    for (char c = *s; c != '\0'; c = s[++n]) {
        if (c == '\n' || c == '\r') {
            if (_control.standardEscapes) {
                wantEscape = true;
            } else {
                wantQuote = true;
            }
        } else if (c == _control.getDelimiter()) {
            wantQuote = true;
        } else if (c == _control.getEscapeChar()) {
            wantQuote = true;
            wantEscape = true;
        } else if (c == _control.getQuoteChar()) {
            if (_control.doubleQuote) {
                wantQuote = true;
            }
            wantEscape = true;
        }
    }
    // if no escapes are necessary, blast chars to output stream, otherwise
    // delegate to _writeQuoted() and _writeUnquoted()
    if (wantQuote) {
        if (wantEscape) {
            _writeQuoted(s);
        } else {
            _out->put(_control.getQuoteChar());
            _out->write(s, n);
            _out->put(_control.getQuoteChar());
        }
    } else if (wantEscape) {
        _writeUnquoted(s);
    } else {
        _out->write(s, n);
    }
}

/** Writes quoted field data, handling any required escaping.
  */
void CsvWriter::_writeQuoted(char const *s) {
    _out->put(_control.getQuoteChar());
    size_t n = 0;
    while (true) {
        char const c = s[n];
        if (c == '\0') {
            _out->write(s, n);
            break;
        }
        if (c == _control.getQuoteChar()) {
            _out->write(s, n);
            s += n + 1;
            n = 0;
            if (_control.doubleQuote) {
                _out->put(c);
                _out->put(c);
            } else if (_control.getEscapeChar() == '\0') {
                _out->put(c);
                throw LSST_EXCEPT(pexExcept::InvalidParameterError,
                                  "Field value requires escaping, but "
                                  "no escape character is set");
            } else {
                _out->put(_control.getEscapeChar());
                _out->put(c);
            }
        } else if (c == _control.getEscapeChar()) {
            _out->write(s, n);
            s += n + 1;
            n = 0;
            _out->put(c);
            _out->put(c);
        } else if (_control.standardEscapes && (c == '\n' || c == '\r')) {
            _out->write(s, n);
            s += n + 1;
            n = 0;
            if (_control.getEscapeChar() == '\0') {
                _out->put(_control.getQuoteChar());
                throw LSST_EXCEPT(pexExcept::InvalidParameterError,
                                  "Field value requires escaping, but "
                                  "no escape character is set");
            }
            _out->put(_control.getEscapeChar());
            _out->put(c == '\n' ? 'n' : 'r');
        } else {
            ++n;
        }
    }
    _out->put(_control.getQuoteChar());
}

/** Writes unquoted field data, handling any required escaping.
  */
void CsvWriter::_writeUnquoted(char const *s) {
    size_t n = 0;
    while (true) {
        char const c = s[n];
        if (c == '\0') {
            _out->write(s, n);
            break;
        }
        if (c == '\n' ||
            c == '\r' ||
            c == _control.getDelimiter() ||
            c == _control.getEscapeChar()) {

            _out->write(s, n);
            s += n + 1;
            n = 0;
            if (_control.getEscapeChar() == '\0') {
                throw LSST_EXCEPT(pexExcept::InvalidParameterError,
                                  "Field value requires escaping, but "
                                  "no escape character is set");
            }
            _out->put(_control.getEscapeChar());
            if (c == '\n') {
                _out->put(_control.standardEscapes ? 'n' : c);
            } else if (c == '\r') {
                _out->put(_control.standardEscapes ? 'r' : c);
            } else {
                _out->put(c);
            }
            continue;
        } else if (c == _control.getQuoteChar() &&
                   _control.quoting != "QUOTE_NONE") {

            _out->write(s, n);
            s += n + 1;
            n = 0;
            if (_control.doubleQuote) {
                _out->put(c);
                _out->put(c);
            } else if (_control.getEscapeChar() == '\0') {
                throw LSST_EXCEPT(pexExcept::InvalidParameterError,
                                  "Field value requires escaping, but "
                                  "no escape character is set");
            } else {
               _out->put(_control.getEscapeChar());
               _out->put(c);
           }
        } else {
            ++n;
        }
    }
}

}}} // namespace lsst::ap::utils
