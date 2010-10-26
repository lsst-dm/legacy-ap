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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Csv

#include "boost/test/unit_test.hpp"

#include <cstdarg>
#include <sstream>

#include "lsst/pex/exceptions.h"
#include "lsst/ap/util/Csv.h"


using std::ostringstream;
using std::istringstream;
using std::vector;
using std::string;

using lsst::ap::util::CsvDialect;
using lsst::ap::util::CsvReader;
using lsst::ap::util::CsvWriter;
using lsst::pex::exceptions::Exception;


BOOST_AUTO_TEST_CASE(dialect) {
    // illegal delimiters should throw
    BOOST_CHECK_THROW(CsvDialect(CsvDialect::QUOTE_NONE, '\0'), Exception);
    BOOST_CHECK_THROW(CsvDialect(CsvDialect::QUOTE_NONE, '\n'), Exception);
    BOOST_CHECK_THROW(CsvDialect(CsvDialect::QUOTE_NONE, '\r'), Exception);
    // quote/escape char equal to delimiter should throw
    BOOST_CHECK_THROW(CsvDialect(CsvDialect::QUOTE_NONE, ',', ','), Exception);
    BOOST_CHECK_THROW(CsvDialect(CsvDialect::QUOTE_NONE, ',', '\0', ','), Exception);
    // identical quote/escape characters should throw
    BOOST_CHECK_THROW(CsvDialect(CsvDialect::QUOTE_NONE, ',', 'a', 'a'), Exception);
    // a quote char of '\0' requires QUOTE_NONE
    BOOST_CHECK_THROW(CsvDialect(CsvDialect::QUOTE_ALL, ',', '\\', '\0'), Exception);
    // an escape char of '\0' requires standardEscapes to be false
    BOOST_CHECK_THROW(CsvDialect(CsvDialect::QUOTE_ALL, ',', '\0', '"', false, false, true), Exception);
    // illegal escape characters should throw
    BOOST_CHECK_THROW(CsvDialect(CsvDialect::QUOTE_NONE, ',', '\n'), Exception);
    BOOST_CHECK_THROW(CsvDialect(CsvDialect::QUOTE_NONE, ',', '\r'), Exception);
    // illegal quote characters should throw
    BOOST_CHECK_THROW(CsvDialect(CsvDialect::QUOTE_NONE, ',', '\\', '\n'), Exception);
    BOOST_CHECK_THROW(CsvDialect(CsvDialect::QUOTE_NONE, ',', '\\', '\r'), Exception);
    // null string containing illegal characters should throw
    BOOST_CHECK_THROW(CsvDialect("\n", CsvDialect::QUOTE_NONE), Exception);
    BOOST_CHECK_THROW(CsvDialect("\r", CsvDialect::QUOTE_NONE), Exception);
    // null string containing delimter/escape/quote characters should throw
    BOOST_CHECK_THROW(CsvDialect(",", CsvDialect::QUOTE_NONE, ','), Exception);
    BOOST_CHECK_THROW(CsvDialect("a", CsvDialect::QUOTE_NONE, ',', 'a'), Exception);
    BOOST_CHECK_THROW(CsvDialect("b", CsvDialect::QUOTE_NONE, ',', 'a', 'b'), Exception);
}


static char const * const NULLF = 0;
static char const * const END = NULLF - 1;
static char const * const ENDR = NULLF - 2;

static void roundTrip(char const * const result,
                      CsvDialect const & dialect,
                      ...)
{
    typedef vector<vector<char const *> >::const_iterator RecIter;
    typedef vector<char const *>::const_iterator FieldIter;

    // read in the fields to output
    vector<vector<char const *> > records;
    std::va_list args;
    va_start(args, dialect);
    while (true) {
        vector<char const *> record;
        char const *s;
        while (true) {
            s = va_arg(args, char const *);
            if (s == ENDR || s == END) {
                break;
            } else {
                record.push_back(s);
            }
        }
        if (s != END || record.size() > 0) {
            records.push_back(record);
        }
        if (s == END) {
            break;
        }
    }

    // write out fields
    ostringstream oss;
    CsvWriter writer(oss, dialect);
    for (RecIter r = records.begin(), er = records.end(); r != er; ++r) {
        for (FieldIter f = r->begin(), ef = r->end(); f != ef; ++f) {
            writer.appendField(*f);
        }
        writer.endRecord();
    }

    // check that desired result is obtained
    BOOST_CHECK_EQUAL(oss.str(), result);

    // check that round tripping works
    istringstream iss(oss.str());
    CsvReader reader(iss, dialect);
    if (records.size() == 0) {
        BOOST_CHECK_EQUAL(reader.isDone(), true);
        return;
    }
    size_t nRec = 1;
    for (RecIter r = records.begin(), er = records.end(); r != er; ++r, ++nRec) {
        BOOST_CHECK_EQUAL(reader.isDone(), false);
        BOOST_CHECK_EQUAL(reader.getNumFields(), static_cast<int>(r->size()));
        BOOST_CHECK_EQUAL(reader.getNumRecords(), nRec);
        int nField = 0;
        for (FieldIter f = r->begin(), ef = r->end(); f != ef; ++f, ++nField) {
            if (*f == NULLF) {
                BOOST_CHECK_EQUAL(reader.isNull(nField), true);
            } else {
                BOOST_CHECK_EQUAL(reader.isNull(nField), false);
                BOOST_CHECK_EQUAL(reader.get(nField), *f);
            }
        }
        reader.nextRecord();
    }
    BOOST_CHECK_EQUAL(reader.isDone(), true);
}

BOOST_AUTO_TEST_CASE(quoting) {
    {
        CsvDialect dialect(CsvDialect::QUOTE_MINIMAL, ',', '\\', '\'');
        roundTrip("a,b,'c,d'\n", dialect, "a", "b", "c,d", END);
        roundTrip("a,b,'c,\\'d\\''\n", dialect, "a", "b", "c,'d'", END);
        roundTrip("\\'\n", dialect, "'", END);
        roundTrip("'a\nb',c\n", dialect, "a\nb", "c", END);
        roundTrip("'a\rb',c\n", dialect, "a\rb", "c", END);
    }
    {
        CsvDialect dialect(CsvDialect::QUOTE_MINIMAL, ',', '\0', '\'');
        roundTrip("a,b,'c,d'\n", dialect, "a", "b", "c,d", END);
        roundTrip("'a\nb,c'\n", dialect, "a\nb,c", END);
        roundTrip("'a\rb',c\n", dialect, "a\rb", "c", END);
        BOOST_CHECK_THROW(roundTrip("", dialect, "a", "b", "c,'d'", END), Exception);
    }
    {
        CsvDialect dialect(CsvDialect::QUOTE_NONE, ',', '\\', '\0');
        roundTrip("a,b,c\\,d\n", dialect, "a", "b", "c,d", END);
        roundTrip("'\n", dialect, "'", END);
        roundTrip("\\,\n", dialect, ",", END);
        roundTrip("\\\n\n", dialect, "\n", END);
    }
    {
        CsvDialect dialect(CsvDialect::QUOTE_NONE, ',', '\0', '\'');
        roundTrip("'\n", dialect, "'", END);
        BOOST_CHECK_THROW(roundTrip("", dialect, "\n", END), Exception);
        BOOST_CHECK_THROW(roundTrip("", dialect, "\r", END), Exception);
        BOOST_CHECK_THROW(roundTrip("", dialect, "a","b","c,d", END), Exception);
    }
    {
        CsvDialect dialect(CsvDialect::QUOTE_ALL, ',', '\\', '\'');
        roundTrip("'a','b','c,d'\n", dialect, "a", "b", "c,d", END);
        roundTrip("'a','b','c,\\'d\\''\n", dialect, "a", "b", "c,'d'", END);
        roundTrip("'a\n','b\r'\n", dialect, "a\n", "b\r", END);
        roundTrip("'\\\\'\n", dialect, "\\", END);
    }
    {
        CsvDialect dialect(CsvDialect::QUOTE_ALL, ',', '\0', '\'');
        BOOST_CHECK_THROW(roundTrip("", dialect, "'", END), Exception);
        roundTrip("','\n", dialect, ",", END);
        roundTrip("'\n'\n", dialect, "\n", END);
    }
}

BOOST_AUTO_TEST_CASE(escaping) {
    {
        CsvDialect dialect(CsvDialect::QUOTE_MINIMAL, ',', '\0', '\'');
        BOOST_CHECK_THROW(roundTrip("", dialect, "'", END), Exception);
    }
    {
        CsvDialect dialect(CsvDialect::QUOTE_MINIMAL, ',', '\\', '\'', false, true, true);
        roundTrip("\\n,\\r\n", dialect, "\n", "\r", END);
        roundTrip("'''',''''\n", dialect, "'", "'", END);
    }
    {
        CsvDialect dialect(CsvDialect::QUOTE_MINIMAL, ',', '\\', '\'', false, true);
        roundTrip("'\n','\r'\n", dialect, "\n", "\r", END);
        roundTrip("'''',''''\n", dialect, "'", "'", END);
    }
    {
        CsvDialect dialect(CsvDialect::QUOTE_NONE, ',', '\\', '\0', false, false, true);
        roundTrip("\\n,\\r\n", dialect, "\n", "\r", END);
        roundTrip("a,b,c\\,d\n", dialect, "a", "b", "c,d", END);
    }
}

BOOST_AUTO_TEST_CASE(null) {
    {
        CsvDialect dialect(CsvDialect::QUOTE_MINIMAL, ',', '\\', '\'', false, false, true);
        roundTrip("a,\\N,b\n", dialect, "a", NULLF, "b", END);
        roundTrip("\\N,\\N,\\N\n", dialect, NULLF, NULLF, NULLF, END);
    }
    {
        CsvDialect dialect("NULL", CsvDialect::QUOTE_ALL, ',', '\\', '\'');
        roundTrip("'a',NULL,'b'\n", dialect, "a", NULLF, "b", END);
        roundTrip("NULL\n", dialect, NULLF, END);
    }
    {
        CsvDialect dialect("NULL", CsvDialect::QUOTE_MINIMAL, ',', '\\', '\'');
        roundTrip("'NULL',NULL,'NULL',NULL\n", dialect, "NULL", NULLF, "NULL", NULLF, END);
        roundTrip("NULL\n", dialect, NULLF, END);
    }
    {
        CsvDialect dialect("", CsvDialect::QUOTE_MINIMAL, ',', '\\', '\'');
        roundTrip("a,,b,''\n", dialect, "a", NULLF, "b", "", END);
        roundTrip("\n\n", dialect, NULLF, ENDR, NULLF, END);
    }
}

BOOST_AUTO_TEST_CASE(trailingDelimiter) {
    CsvDialect dialect("", CsvDialect::QUOTE_NONE, '|', '\\', '\0', false, false, false, true);
    roundTrip("a|b|\n", dialect, "a", "b", END);
    roundTrip("a|b|\nc|d|\n", dialect, "a", "b", ENDR, "c", "d", END);
    roundTrip("||a||\n", dialect, NULLF, NULLF, "a", NULLF, END);
}

