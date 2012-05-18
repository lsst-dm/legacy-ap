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
  * @brief  CSV format control.
  * @author Serge Monkewitz
  */
#ifndef LSST_AP_UTILS_CSVCONTROL_H
#define LSST_AP_UTILS_CSVCONTROL_H

#include <string>

#include "lsst/pex/config.h"


namespace lsst { namespace ap { namespace utils {

/** @brief Parameters that define a Character-Separated-Value dialect.
  *
  * These are intended to cover variations on the ubiquitous but often
  * inconsistently defined comma-separated-value and tab-separated-value
  * formats, and largely follow the parameters for dialect specification
  * in the python csv module. One important goal is to allow the CSV output
  * of various RDBMSes to be read and written. Accordingly, there are
  * parameters that allow database NULLs to be recognized.
  */
struct CsvControl {
    enum Quoting {
        QUOTE_NONE = 0, ///< Never quote fields.
        QUOTE_ALL,      ///< Always quote fields.
        QUOTE_MINIMAL   ///< Only quote fields when necessary.
    };

    CsvControl();
    ~CsvControl();

    LSST_CONTROL_FIELD(null, std::string,
        "String representation of NULL field values. Never quoted on output.\n"
        "If specified, the representation may not contain any delimiter,\n"
        "quote, escape or line terminator characters ('\\n'/'\\r').\n");

    LSST_CONTROL_FIELD(hasNull, bool,
        "Indicates whether the null string is valid. If set to false, the only\n"
        "way NULLs can be recognized/written is if standardEscapes is set to\n"
        "true (in which case '\\N' is mapped to NULL, assuming that '\\' is\n"
        "the escape character).\n");

    LSST_CONTROL_FIELD(quoting, std::string,
        "Field quoting style for CSV input/output. Legal values are:\n"
        "\n"
        "'QUOTE_MINIMAL': Only quote fields when necessary - for instance,\n"
        "                 when a field value contains a delimiter character.\n"
        "\n"
        "'QUOTE_NONE':    Never quote fields.\n"
        "\n"
        "'QUOTE_ALL':     Always quote fields.\n");

    LSST_CONTROL_FIELD(delimiter, std::string,
        "A one character string containing the field delimiter character.\n"
        "Values of '\\0', '\\n', or '\\r' are illegal.\n");

    LSST_CONTROL_FIELD(escapeChar, std::string,
        "A one character string containing the escape character. An empty\n"
        "string is mapped to an escape character of '\\0', which disables\n"
        "escaping. A value equal to the delimiter or quote character is\n"
        "illegal, with one exception: both the escape and quote characters\n"
        "can be '\\0'.  The '\\n' and '\\r' characters are also illegal.\n");

    LSST_CONTROL_FIELD(quoteChar, std::string,
        "A one character string containing the character used to quote fields\n"
        "when quoting is set to either 'QUOTE_ALL' or 'QUOTE_MINIMAL'. An\n"
        "empty string is mapped to a quote character of '\\0', which disables\n"
        "quoting (and is only legal when quoting is set to 'QUOTE_NONE'). A\n"
        "value equal to the delimiter or quote character is illegal, with one\n"
        "exception: both the escape and quote characters can be '\\0'. The\n"
        "'\\n' and '\\r' characters are also illegal.\n");

    LSST_CONTROL_FIELD(skipInitialSpace, bool,
        "If true, whitespace immediately following the delimiter is ignored.");

    LSST_CONTROL_FIELD(doubleQuote, bool,
        "If true, embedded quote characters are escaped with a leading quote\n"
        "character. Otherwise the escape character is used. If escaping and\n"
        "double-quoting are disabled, writing a field  with embedded quote\n"
        "character will raise an exception.\n");

    LSST_CONTROL_FIELD(standardEscapes, bool,
        "Flag indicating whether standard escape sequences should be handled.\n"
        "If false, then the character sequence '\\C', where C is any character,\n"
        "is mapped to C (assuming '\\' is the escape character). If true,\n"
        "the following special cases are handled differently:\n"
        "\n"
        "- '\\b' is mapped to BS - backspace (ASCII 8)\n"
        "- '\\f' is mapped to FF - form feed (ASCII 12)\n"
        "- '\\n' is mapped to NL - newline (ASCII 10)\n"
        "- '\\r' is mapped to CR - carriage return (ASCII 13)\n"
        "- '\\t' is mapped to TAB - horizontal tab (ASCII 9)\n"
        "- '\\v' is mapped to VT - vertical tab (ASCII 11)\n"
        "- '\\xD' and '\\xDD', where D is a hexadecimal digit, is mapped to\n"
        "  the character with that numeric code.\n"
        "- A field value of exactly '\\N' (no quotes, whitespace, or other\n"
        "  content) is treated as a NULL.\n");

    LSST_CONTROL_FIELD(trailingDelimiter, bool,
        "If true, then a trailing delimiter character is expected and written\n"
        "at end of every record, immediately preceding the line terminator.\n");

    LSST_CONTROL_FIELD(nonfiniteAsNull, bool,
        "If true, then non-finite (NaN, Inf, -Inf) floating point values are\n"
        "written out as NULL field values.\n");

    /** Returns true if database NULLs are recognizable in this dialect.
      * This is the case when hasNull is true, or when escapeChar is not
      * '\\0' and standardEscapes is true (in which case "\N" is recognized
      * as a NULL).
      */
    inline bool isNullRecognizable() const {
        return hasNull || standardEscapes;
    }

    /** Returns the dialects quoting style. Note that if this is equal to
      * QUOTE_NONE, no special handling of quote characters on reading is
      * performed.
      */
    inline Quoting getQuoting() const {
        if (quoting == "QUOTE_MINIMAL") {
            return QUOTE_MINIMAL;
        } else if (quoting == "QUOTE_ALL") {
            return QUOTE_ALL;
        }
        return QUOTE_NONE;
    }

    inline char getDelimiter() const {
        return delimiter.c_str()[0];
    }
    inline char getEscapeChar() const {
        return escapeChar.c_str()[0];
    }
    inline char getQuoteChar() const {
        return quoteChar.c_str()[0];
    }

    void validate() const;
};

}}} // namespace lsst::ap::utils

#endif // LSST_AP_UTILS_CSVCONTROL_H

