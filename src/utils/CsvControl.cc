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
  * @brief  CsvControl implementations.
  * @author Serge Monkewitz
  */
#include "lsst/ap/utils/CsvControl.h"

#include "lsst/pex/exceptions.h"


namespace pexExcept = lsst::pex::exceptions;
using std::string;

namespace lsst { namespace ap { namespace utils {

// -- CsvControl implementation ----

/** CsvControl default constructor.
  */
CsvControl::CsvControl() :
    null(),
    hasNull(false),
    quoting("QUOTE_MINIMAL"),
    delimiter(","),
    escapeChar("\\"),
    quoteChar("\""),
    skipInitialSpace(false),
    doubleQuote(false),
    standardEscapes(true),
    trailingDelimiter(false)
{
    validate();
}

CsvControl::~CsvControl() { }

void CsvControl::validate() const {
    if (delimiter.size() != 1) {
        throw LSST_EXCEPT(pexExcept::InvalidParameterException,
                          "delimiter must consist of a single character");
    }
    if (delimiter[0] == '\0' || delimiter[0] == '\n' || delimiter[0] == '\r') {
        throw LSST_EXCEPT(pexExcept::InvalidParameterException,
                          "delimiter equal to '\\[0nr]'");

    }
    if (escapeChar.size() > 1) {
        throw LSST_EXCEPT(pexExcept::InvalidParameterException,
                          "escapeChar string contains more than one character");

    }
    if (escapeChar == delimiter ||
        getEscapeChar() == '\n' || getEscapeChar() == '\r') {
        throw LSST_EXCEPT(pexExcept::InvalidParameterException,
                          "escapeChar equal to delimiter or '\\[nr]'.");
    }
    if (quoteChar.size() > 1) {
        throw LSST_EXCEPT(pexExcept::InvalidParameterException,
                          "quoteChar string contains more than one character");
    }
    if (quoteChar == delimiter ||
        getQuoteChar() == '\n' || getQuoteChar() == '\r') {
        throw LSST_EXCEPT(pexExcept::InvalidParameterException,
                          "quoteChar equal to delimiter or '\\[nr]'.");
    }
    if (escapeChar == quoteChar && getEscapeChar() != '\0') {
        throw LSST_EXCEPT(pexExcept::InvalidParameterException,
                          "escapeChar equal to quoteChar. Did you mean to use "
                          "the doubleQuote option instead?");
    }
    if (getEscapeChar() == '\0' && standardEscapes) {
        throw LSST_EXCEPT(pexExcept::InvalidParameterException,
                          "escapeChar set to '\\0', but standardEscapes = true");
    }
    if (getQuoteChar() == '\0' && quoting != "QUOTE_NONE") {
        throw LSST_EXCEPT(pexExcept::InvalidParameterException,
                          "quoteChar set to '\\0', but quoting is not "
                          "'QUOTE_NONE'");
    }
    if (null.find('\n') != string::npos ||
        null.find('\r') != string::npos ||
        null.find(getDelimiter()) != string::npos ||
        null.find(getEscapeChar()) != string::npos ||
        null.find(getQuoteChar()) != string::npos) {
        throw LSST_EXCEPT(pexExcept::InvalidParameterException,
                          "null string contains '\\n', '\\r', delimiter, "
                          "escapeChar, or quoteChar.");
    }
}

}}} // namespace lsst::ap::utils
