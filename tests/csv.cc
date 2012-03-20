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

#include <cfloat>
#include <cmath>
#include <cstdarg>
#include <limits>
#include <sstream>

#include "lsst/utils/ieee.h"
#include "lsst/pex/exceptions.h"
#include "lsst/afw/math/Random.h"

#include "lsst/ap/utils/Csv.h"


using std::ostringstream;
using std::istringstream;
using std::numeric_limits;
using std::string;
using std::vector;

using lsst::pex::exceptions::Exception;
using lsst::afw::math::Random;
using lsst::ap::utils::CsvControl;
using lsst::ap::utils::CsvReader;
using lsst::ap::utils::CsvWriter;


BOOST_AUTO_TEST_CASE(control) {
    // illegal delimiters should throw
    CsvControl ctrl;
    ctrl.quoting = "QUOTE_NONE";
    ctrl.delimiter = "";
    BOOST_CHECK_THROW(ctrl.validate(), Exception);
    ctrl.delimiter = "\n";
    BOOST_CHECK_THROW(ctrl.validate(), Exception);
    ctrl.delimiter = "\r";
    BOOST_CHECK_THROW(ctrl.validate(), Exception);
    ctrl.delimiter = "ab";
    BOOST_CHECK_THROW(ctrl.validate(), Exception);

    // quote/escape char equal to delimiter should throw
    ctrl.delimiter = ",";
    ctrl.escapeChar = ",";
    BOOST_CHECK_THROW(ctrl.validate(), Exception);
    ctrl.escapeChar = "";
    ctrl.quoteChar = ",";
    BOOST_CHECK_THROW(ctrl.validate(), Exception);

    // identical quote/escape characters should throw
    ctrl.escapeChar = "a";
    ctrl.quoteChar = "a";
    BOOST_CHECK_THROW(ctrl.validate(), Exception);

    // quote/escape strings must be at most one character in length
    ctrl.escapeChar = "ab";
    BOOST_CHECK_THROW(ctrl.validate(), Exception);
    ctrl.escapeChar = "\\";
    ctrl.quoteChar = "quote";
    BOOST_CHECK_THROW(ctrl.validate(), Exception);

    // a quote char of '\0' requires QUOTE_NONE
    ctrl.quoting = "QUOTE_ALL";
    ctrl.quoteChar = "";
    BOOST_CHECK_THROW(ctrl.validate(), Exception);
    ctrl.quoting = "QUOTE_MINIMAL";
    BOOST_CHECK_THROW(ctrl.validate(), Exception);

    // an escape char of '\0' requires standardEscapes to be false
    ctrl.quoteChar = "\"";
    ctrl.escapeChar = "";
    ctrl.standardEscapes = true;
    BOOST_CHECK_THROW(ctrl.validate(), Exception);

    // illegal escape characters should throw
    ctrl.escapeChar = "\n";
    BOOST_CHECK_THROW(ctrl.validate(), Exception);
    ctrl.escapeChar = "\r";
    BOOST_CHECK_THROW(ctrl.validate(), Exception);

    // illegal quote characters should throw
    ctrl.escapeChar = "\\";
    ctrl.quoteChar = "\n";
    BOOST_CHECK_THROW(ctrl.validate(), Exception);
    ctrl.quoteChar = "\r";
    BOOST_CHECK_THROW(ctrl.validate(), Exception);

    // null string containing illegal characters should throw
    ctrl.quoteChar = "\"";
    ctrl.null = "\n";
    ctrl.hasNull = true;
    BOOST_CHECK_THROW(ctrl.validate(), Exception);
    ctrl.null = "\r";
    BOOST_CHECK_THROW(ctrl.validate(), Exception);

    // null string containing delimter/escape/quote characters should throw
    ctrl.null = ",";
    BOOST_CHECK_THROW(ctrl.validate(), Exception);
    ctrl.null = "\\";
    BOOST_CHECK_THROW(ctrl.validate(), Exception);
    ctrl.null = "\"";
    BOOST_CHECK_THROW(ctrl.validate(), Exception);
}


static char const * const NULLF = 0;
static char const * const END = NULLF - 1;
static char const * const ENDR = NULLF - 2;

static void roundTrip(char const * const result,
                      CsvControl const & control,
                      ...)
{
    typedef vector<vector<char const *> >::const_iterator RecIter;
    typedef vector<char const *>::const_iterator FieldIter;

    // read in the fields to output
    vector<vector<char const *> > records;
    std::va_list args;
    va_start(args, control);
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
    CsvWriter writer(oss, control);
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
    CsvReader reader(iss, control);
    if (records.size() == 0) {
        BOOST_CHECK(reader.isDone());
        return;
    }
    size_t nRec = 1;
    for (RecIter r = records.begin(), er = records.end(); r != er; ++r, ++nRec) {
        BOOST_CHECK(!reader.isDone());
        BOOST_CHECK_EQUAL(reader.getNumFields(), static_cast<int>(r->size()));
        BOOST_CHECK_EQUAL(reader.getNumRecords(), nRec);
        int nField = 0;
        for (FieldIter f = r->begin(), ef = r->end(); f != ef; ++f, ++nField) {
            if (*f == NULLF) {
                BOOST_CHECK(reader.isNull(nField));
            } else {
                BOOST_CHECK(!reader.isNull(nField));
                BOOST_CHECK_EQUAL(reader.get(nField), *f);
            }
        }
        reader.nextRecord();
    }
    BOOST_CHECK(reader.isDone());
}

BOOST_AUTO_TEST_CASE(quoting) {
    {
        CsvControl ctrl;
        ctrl.null = "";
        ctrl.hasNull = false;
        ctrl.quoting = "QUOTE_MINIMAL";
        ctrl.delimiter = ",";
        ctrl.escapeChar = "\\";
        ctrl.quoteChar = "'";
        ctrl.skipInitialSpace = false;
        ctrl.doubleQuote = false;
        ctrl.standardEscapes = false;
        ctrl.trailingDelimiter = false;
        ctrl.validate();
        roundTrip("a,b,'c,d'\n", ctrl, "a", "b", "c,d", END);
        roundTrip("a,b,'c,\\'d\\''\n", ctrl, "a", "b", "c,'d'", END);
        roundTrip("\\'\n", ctrl, "'", END);
        roundTrip("'a\nb',c\n", ctrl, "a\nb", "c", END);
        roundTrip("'a\rb',c\n", ctrl, "a\rb", "c", END);
    }
    {
        CsvControl ctrl;
        ctrl.null = "";
        ctrl.hasNull = false;
        ctrl.quoting = "QUOTE_MINIMAL";
        ctrl.delimiter = ","; 
        ctrl.escapeChar = "";
        ctrl.quoteChar = "'";
        ctrl.skipInitialSpace = false;
        ctrl.doubleQuote = false;
        ctrl.standardEscapes = false;
        ctrl.trailingDelimiter = false;
        ctrl.validate();
        roundTrip("a,b,'c,d'\n", ctrl, "a", "b", "c,d", END);
        roundTrip("'a\nb,c'\n", ctrl, "a\nb,c", END);
        roundTrip("'a\rb',c\n", ctrl, "a\rb", "c", END);
        BOOST_CHECK_THROW(roundTrip("", ctrl, "a", "b", "c,'d'", END), Exception);
    }
    {
        CsvControl ctrl;
        ctrl.null = "";
        ctrl.hasNull = false;
        ctrl.quoting = "QUOTE_NONE";
        ctrl.delimiter = ",";
        ctrl.escapeChar = "\\";
        ctrl.quoteChar = "";
        ctrl.skipInitialSpace = false;
        ctrl.doubleQuote = false;
        ctrl.standardEscapes = false;
        ctrl.trailingDelimiter = false;
        ctrl.validate();
        roundTrip("a,b,c\\,d\n", ctrl, "a", "b", "c,d", END);
        roundTrip("'\n", ctrl, "'", END);
        roundTrip("\\,\n", ctrl, ",", END);
        roundTrip("\\\n\n", ctrl, "\n", END);
    }
    {
        CsvControl ctrl;
        ctrl.null = "";
        ctrl.hasNull = false;
        ctrl.quoting = "QUOTE_NONE";
        ctrl.delimiter = ",";
        ctrl.escapeChar = "";
        ctrl.quoteChar = "'";
        ctrl.skipInitialSpace = false;
        ctrl.doubleQuote = false;
        ctrl.standardEscapes = false;
        ctrl.trailingDelimiter = false;
        ctrl.validate();
        roundTrip("'\n", ctrl, "'", END);
        BOOST_CHECK_THROW(roundTrip("", ctrl, "\n", END), Exception);
        BOOST_CHECK_THROW(roundTrip("", ctrl, "\r", END), Exception);
        BOOST_CHECK_THROW(roundTrip("", ctrl, "a","b","c,d", END), Exception);
    }
    {
        CsvControl ctrl;
        ctrl.null = "";
        ctrl.hasNull = false;
        ctrl.quoting = "QUOTE_ALL";
        ctrl.delimiter = ",";
        ctrl.escapeChar = "\\";
        ctrl.quoteChar = "'";
        ctrl.skipInitialSpace = false;
        ctrl.doubleQuote = false;
        ctrl.standardEscapes = false;
        ctrl.trailingDelimiter = false;
        ctrl.validate();
        roundTrip("'a','b','c,d'\n", ctrl, "a", "b", "c,d", END);
        roundTrip("'a','b','c,\\'d\\''\n", ctrl, "a", "b", "c,'d'", END);
        roundTrip("'a\n','b\r'\n", ctrl, "a\n", "b\r", END);
        roundTrip("'\\\\'\n", ctrl, "\\", END);
    }
    {
        CsvControl ctrl;
        ctrl.null = "";
        ctrl.hasNull = false;
        ctrl.quoting = "QUOTE_ALL";
        ctrl.delimiter = ",";
        ctrl.escapeChar = "";
        ctrl.quoteChar = "'";
        ctrl.skipInitialSpace = false;
        ctrl.doubleQuote = false;
        ctrl.standardEscapes = false;
        ctrl.trailingDelimiter = false;
        ctrl.validate();
        BOOST_CHECK_THROW(roundTrip("", ctrl, "'", END), Exception);
        roundTrip("','\n", ctrl, ",", END);
        roundTrip("'\n'\n", ctrl, "\n", END);
    }
}

BOOST_AUTO_TEST_CASE(escaping) {
    {
        CsvControl ctrl;
        ctrl.null = "";
        ctrl.hasNull = false;
        ctrl.quoting = "QUOTE_MINIMAL";
        ctrl.delimiter = ",";
        ctrl.escapeChar = "";
        ctrl.quoteChar = "'";
        ctrl.skipInitialSpace = false;
        ctrl.doubleQuote = false;
        ctrl.standardEscapes = false;
        ctrl.trailingDelimiter = false;
        ctrl.validate();
        BOOST_CHECK_THROW(roundTrip("", ctrl, "'", END), Exception);
    }
    {
        CsvControl ctrl;
        ctrl.null = "";
        ctrl.hasNull = false;
        ctrl.quoting = "QUOTE_MINIMAL";
        ctrl.delimiter = ",";
        ctrl.escapeChar = "\\";
        ctrl.quoteChar = "'";
        ctrl.skipInitialSpace = false;
        ctrl.doubleQuote = true;
        ctrl.standardEscapes = true;
        ctrl.trailingDelimiter = false;
        ctrl.validate();
        roundTrip("\\n,\\r\n", ctrl, "\n", "\r", END);
        roundTrip("'''',''''\n", ctrl, "'", "'", END);
    }
    {
        CsvControl ctrl;
        ctrl.null = "";
        ctrl.hasNull = false;
        ctrl.quoting = "QUOTE_MINIMAL";
        ctrl.delimiter = ",";
        ctrl.escapeChar = "\\";
        ctrl.quoteChar = "'";
        ctrl.skipInitialSpace = false;
        ctrl.doubleQuote = true;
        ctrl.standardEscapes = false;
        ctrl.trailingDelimiter = false;
        ctrl.validate();
        roundTrip("'\n','\r'\n", ctrl, "\n", "\r", END);
        roundTrip("'''',''''\n", ctrl, "'", "'", END);
    }
    {
        CsvControl ctrl;
        ctrl.null = "";
        ctrl.hasNull = false;
        ctrl.quoting = "QUOTE_NONE";
        ctrl.delimiter = ",";
        ctrl.escapeChar = "\\";
        ctrl.quoteChar = "";
        ctrl.skipInitialSpace = false;
        ctrl.doubleQuote = false;
        ctrl.standardEscapes = true;
        ctrl.trailingDelimiter = false;
        ctrl.validate();
        roundTrip("\\n,\\r\n", ctrl, "\n", "\r", END);
        roundTrip("a,b,c\\,d\n", ctrl, "a", "b", "c,d", END);
    }
}

BOOST_AUTO_TEST_CASE(null) {
    {
        CsvControl ctrl;
        ctrl.null = "";
        ctrl.hasNull = false;
        ctrl.quoting = "QUOTE_MINIMAL";
        ctrl.delimiter = ",";
        ctrl.escapeChar = "\\";
        ctrl.quoteChar = "'";
        ctrl.skipInitialSpace = false;
        ctrl.doubleQuote = false;
        ctrl.standardEscapes = true;
        ctrl.trailingDelimiter = false;
        ctrl.validate();
        roundTrip("a,\\N,b\n", ctrl, "a", NULLF, "b", END);
        roundTrip("\\N,\\N,\\N\n", ctrl, NULLF, NULLF, NULLF, END);
    }
    {
        CsvControl ctrl;
        ctrl.null = "NULL";
        ctrl.hasNull = true;
        ctrl.quoting = "QUOTE_ALL";
        ctrl.delimiter = ",";
        ctrl.escapeChar = "\\";
        ctrl.quoteChar = "'";
        ctrl.skipInitialSpace = false;
        ctrl.doubleQuote = false;
        ctrl.standardEscapes = false;
        ctrl.trailingDelimiter = false;
        ctrl.validate();
        roundTrip("'a',NULL,'b'\n", ctrl, "a", NULLF, "b", END);
        roundTrip("NULL\n", ctrl, NULLF, END);
    }
    {
        CsvControl ctrl;
        ctrl.null = "NULL";
        ctrl.hasNull = true;
        ctrl.quoting = "QUOTE_MINIMAL";
        ctrl.delimiter = ",";
        ctrl.escapeChar = "\\";
        ctrl.quoteChar = "'";
        ctrl.skipInitialSpace = false;
        ctrl.doubleQuote = false;
        ctrl.standardEscapes = false;
        ctrl.trailingDelimiter = false;
        ctrl.validate();
        roundTrip("'NULL',NULL,'NULL',NULL\n", ctrl, "NULL", NULLF, "NULL", NULLF, END);
        roundTrip("NULL\n", ctrl, NULLF, END);
    }
    {
        CsvControl ctrl;
        ctrl.null = "";
        ctrl.hasNull = true;
        ctrl.quoting = "QUOTE_MINIMAL";
        ctrl.delimiter = ",";
        ctrl.escapeChar = "\\";
        ctrl.quoteChar = "'";
        ctrl.skipInitialSpace = false;
        ctrl.doubleQuote = false;
        ctrl.standardEscapes = false;
        ctrl.trailingDelimiter = false;
        ctrl.validate();
        roundTrip("a,,b,''\n", ctrl, "a", NULLF, "b", "", END);
        roundTrip("\n\n", ctrl, NULLF, ENDR, NULLF, END);
    }
}

BOOST_AUTO_TEST_CASE(trailingDelimiter) {
    CsvControl ctrl;
    ctrl.null = "";
    ctrl.hasNull = true;
    ctrl.quoting = "QUOTE_NONE";
    ctrl.delimiter = "|";
    ctrl.escapeChar = "\\";
    ctrl.quoteChar = "";
    ctrl.skipInitialSpace = false;
    ctrl.doubleQuote = false;
    ctrl.standardEscapes = false;
    ctrl.trailingDelimiter = true;
    ctrl.validate();
    roundTrip("a|b|\n", ctrl, "a", "b", END);
    roundTrip("a|b|\nc|d|\n", ctrl, "a", "b", ENDR, "c", "d", END);
    roundTrip("||a||\n", ctrl, NULLF, NULLF, "a", NULLF, END);
}

template <typename T>
static void fpIoTest(T const tolerance) {
    static size_t const NREC = 1000;
    static size_t const NFIELD = 16;
    static unsigned long const EXP_RANGE = static_cast<unsigned long>(
        numeric_limits<T>::max_exponent - numeric_limits<T>::min_exponent);
    Random rng;
    ostringstream oss;

    CsvControl ctrl;
    ctrl.null = "";
    ctrl.hasNull = true;
    ctrl.quoting = "QUOTE_MINIMAL";
    ctrl.delimiter = "|";
    ctrl.escapeChar = "\\";
    ctrl.quoteChar = "'";
    ctrl.skipInitialSpace = false;
    ctrl.doubleQuote = false;
    ctrl.standardEscapes = false;
    ctrl.trailingDelimiter = false;
    ctrl.validate();

    CsvWriter writer(oss, ctrl);
    std::vector<T> input;
    for (size_t r = 0; r < NREC; ++r) {
        for (size_t i = 0; i < NFIELD; ++i) {
            int e = static_cast<int>(rng.uniformInt(EXP_RANGE)) +
                    numeric_limits<T>::min_exponent;
            T d = (2.0*static_cast<T>(rng.uniform()) - 1.0);
            for (; e > 0; --e) { d *= 2.0; }
            for (; e < 0; ++e) { d *= 0.5; }
            input.push_back(d);
            writer.appendField(d);
        }
        writer.endRecord();
    }
    // write out NaN, +-INF, and a numbers that underflow/overflow T
    writer.appendField(numeric_limits<T>::quiet_NaN());
    writer.appendField(numeric_limits<T>::infinity());
    writer.appendField(-numeric_limits<T>::infinity());
    writer.appendField("1.0e-50000");
    writer.appendField("1.0e50000");
    writer.appendField("-1.0e50000");
    writer.endRecord();
    // check insensitivity to whitespace, and that exceptions are thrown
    // for illegal values
    writer.appendField("1.0    ");
    writer.appendField("    1.0");
    writer.appendField("  0.0   ");
    writer.appendField("  a   ");
    writer.appendField("  1..0   ");
    writer.appendField("1.0   abc");
    writer.endRecord();

    istringstream iss(oss.str());
    CsvReader reader(iss, ctrl);
    bool warned = false;
    for (size_t r = 0; r < NREC; ++r) {
        BOOST_CHECK_EQUAL(r + 1, reader.getNumRecords());
        for (size_t f = 0; f < NFIELD; ++f) {
            T v1 = input[r*NFIELD + f];
            T v2 = reader.get<T>(static_cast<int>(f));
            if (v1 != v2) {
                BOOST_CHECK_CLOSE(v1, v2, tolerance);
                if (!warned) {
                    BOOST_WARN("The binary to decimal conversion routines of "
                               "your system do not appear to be correctly "
                               "rounded");
                }
                warned = true;
            }
        }
        BOOST_CHECK_EQUAL(reader.isDone(), false);
        reader.nextRecord();
    }

    // check specials
    BOOST_CHECK(!reader.isDone());
    BOOST_CHECK(lsst::utils::isnan(reader.get<T>(0)));
    BOOST_CHECK(lsst::utils::isinf(reader.get<T>(1)) && reader.get<T>(1) > 0);
    BOOST_CHECK(lsst::utils::isinf(reader.get<T>(2)) && reader.get<T>(2) < 0);
    BOOST_CHECK_EQUAL(reader.get<T>(3), 0.0);
    BOOST_CHECK(lsst::utils::isinf(reader.get<T>(4)) && reader.get<T>(4) > 0);
    BOOST_CHECK(lsst::utils::isinf(reader.get<T>(5)) && reader.get<T>(5) < 0);

    // check whitespace/illegal values
    reader.nextRecord();
    BOOST_CHECK(!reader.isDone());
    BOOST_CHECK_EQUAL(reader.get<T>(0), static_cast<T>(1));
    BOOST_CHECK_EQUAL(reader.get<T>(1), static_cast<T>(1));
    BOOST_CHECK_EQUAL(reader.get<T>(2), static_cast<T>(0));
    BOOST_CHECK_THROW(reader.get<T>(3), Exception);
    BOOST_CHECK_THROW(reader.get<T>(4), Exception);
    BOOST_CHECK_THROW(reader.get<T>(5), Exception);

    reader.nextRecord();
    BOOST_CHECK(reader.isDone());
}

BOOST_AUTO_TEST_CASE(floatIo) {
    fpIoTest(FLT_EPSILON*2.0f);
}
BOOST_AUTO_TEST_CASE(doubleIo) {
    fpIoTest(DBL_EPSILON*2.0);
}
BOOST_AUTO_TEST_CASE(longDoubleIo) {
    fpIoTest(LDBL_EPSILON*2.0l);
}

BOOST_AUTO_TEST_CASE(hardRounding) {
    static double const TWO_ULP = DBL_EPSILON*2.0;
    // Fred Tydeman's hard IEEE 754 double precision binary float to decimal
    // conversion cases, taken from a 1996 comp.arch.arithmetic post.
    static double const TYDEMAN_HARD[] = {
        9e0306, 4e-079, 7e-261, 6e-025,
        7e-161, 7e0289, 5e0079, 1e0080,
        7e-303, 5e0152, 5e0125, 2e0126,
        7e-141, 4e-192, 9e0043, 1e0303,

        95e-089, 85e0194, 69e0267, 97e-019,
        37e0046, 74e0046, 61e-099, 53e-208,
        93e-234, 79e-095, 87e-274, 83e0025,
        17e-036, 53e0033, 51e-074, 63e-022,

        839e0143, 749e-182, 999e-026, 345e0266,
        914e-102, 829e0102, 307e0090, 859e0182,
        283e0085, 589e0187, 302e0176, 604e0176,
        761e-244, 647e0230, 755e0174, 255e-075,

        3391e0055, 4147e-015, 3996e-026, 1998e-026,
        3338e-296, 1669e-296, 8699e-276, 5311e0243,
        7903e-096, 7611e-226, 3257e0058, 6514e0058,
        3571e0263, 7142e0263, 5311e0242, 1617e-063,

        51881e0037, 31441e-118, 30179e0079, 60358e0079,
        63876e-020, 31938e-020, 46073e-032, 32941e0051,
        82081e0041, 38701e-215, 62745e0047, 12549e0048,
        64009e-183, 89275e0261, 75859e0025, 57533e0287,

        584169e0229, 940189e-112, 416121e0197, 832242e0197,
        584738e0076, 933587e-140, 252601e0121, 358423e0274,
        892771e-213, 410405e0040, 928609e-261, 302276e-254,
        920657e-023, 609019e-025, 252601e0120, 654839e-060,

        8823691e0130, 2920845e0228, 9210917e0080, 5800419e-303,
        6119898e-243, 3059949e-243, 2572231e0223, 5444097e-021,
        5783893e-127, 3865421e-225, 4590831e0156, 9181662e0156,
        5906361e-027, 7315057e0235, 9088115e0106, 1817623e0107,

        44118455e0129, 35282041e0293, 31279898e-291, 15639949e-291,
        27966061e0145, 55932122e0145, 70176353e-053, 40277543e-032,
        50609263e0157, 66094077e0077, 84863171e0114, 89396333e0264,
        87575437e-309, 78693511e-044, 90285923e-206, 30155207e-030,

        245540327e0121, 263125459e0287, 566446538e-257, 283223269e-257,
        245540327e0122, 491080654e0122, 971212611e-126, 229058583e0052,
        325270231e0039, 989648089e-035, 653777767e0273, 923091487e0209,
        526250918e0288, 350301748e-309, 741111169e-203, 667284113e-240,

        1227701635e0120, 9981396317e-182, 5232604057e-298, 5572170023e-088,
        1964322616e0122, 3928645232e0122, 8715380633e-058, 4856063055e-127,
        8336960483e-153, 1007046393e-155, 5378822089e-176, 5981342308e-190,
        7214782613e-086, 5458466829e0142, 9078555839e-109, 6418488827e0079,

        65325840981e0069, 49573485983e0089, 46275205733e0074, 92550411466e0074,
        41129842097e-202, 93227267727e-049, 41297294357e0185, 41534892987e-067,
        42333842451e0201, 78564021519e-227, 53587107423e-061, 53827010643e-200,
        83356057653e0193, 45256834646e-118, 45392779195e-110, 23934638219e0291,

        995779191233e0113, 997422852243e-265,
        653532977297e-123, 938885684947e0147,
        619534293513e0124, 539879452414e-042,
        742522891517e0259, 254901016865e-022,
        685763015669e0280, 384865004907e-285,
        286556458711e0081, 573112917422e0081,
        769525178383e-150, 416780288265e0192,
        226963895975e-111, 665592809339e0063,

        3891901811465e0217, 4764593340755e0069,
        6336156586177e0269, 8233559360849e0095,
        3662265515198e-107, 1831132757599e-107,
        7812878489261e-179, 6363857920591e0145,
        8811915538555e0082, 9997878507563e-195,
        9224786422069e-291, 6284426329974e-294,
        9199302046091e-062, 6070482281213e-122,
        2780161250963e-301, 8233559360849e0094,

        72027097041701e0206, 97297545286625e0215,
        99021992302453e-025, 54104687080198e-022,
        33519685743233e0089, 67039371486466e0089,
        39064392446305e-180, 17796979903653e0261,
        28921916763211e0038, 87605699161665e0155,
        41921560615349e-067, 80527976643809e0061,
        72335858886654e-159, 52656615219377e0102,
        15400733123779e-072, 77003665618895e-073,

        475603213226859e-042, 972708181182949e0116,
        246411729980464e-071, 123205864990232e-071,
        609610927149051e-255, 475603213226859e-041,
        672574798934795e0065, 134514959786959e0066,
        294897574603217e-151, 723047919080275e0036,
        660191429952702e-088, 330095714976351e-088,
        578686871093232e-159, 144671717773308e-159,
        385018328094475e-074, 330095714976351e-089,

        2215901545757777e-212, 1702061899637397e-276,
        1864950924021923e0213, 3729901848043846e0213,
        7487252720986826e-165, 3743626360493413e-165,
        4988915232824583e0119, 3771476185376383e0277,
        6182410494241627e-119, 2572981889477453e0142,
        7793560217139653e0051, 9163942927285259e-202,
        6353227084707473e0155, 4431803091515554e-211,
        9324754620109615e0211, 8870461176410409e0263,

        90372559027740405e0143, 18074511805548081e0146,
        54897030182071313e0029, 76232626624829156e-032,
        59898021767894608e-165, 29949010883947304e-165,
        26153245263757307e0049, 27176258005319167e-261,
        18074511805548081e0147, 24691002732654881e-115,
        58483921078398283e0057, 64409240769861689e-159,
        94080055902682397e-242, 31766135423537365e0154,
        68985865317742005e0164, 13797173063548401e0165,

        902042358290366539e-281, 238296178309629163e0272,
        783308178698887621e0226, 439176241456570504e0029,
        899810892172646163e0283, 926145344610700019e-225,
        653831131593932675e0047, 130766226318786535e0048,
        557035730189854663e-294, 902042358290366539e-280,
        272104041512242479e0200, 544208083024484958e0200,
        680429695511221511e0192, 308975121073410857e0236,
        792644927852378159e0078, 783308178698887621e0223,

        8396094300569779681e-252, 3507665085003296281e-074,
        7322325862592278999e0074, 6014546754280072926e0209,
        7120190517612959703e0120, 3507665085003296281e-073,
        4345544743100783551e-218, 9778613303868468131e-090,
        7539204280836061195e-082, 7862637540082247119e-202,
        2176832332097939832e0200, 8643988913946659879e0115,
        5529436763613147623e0138, 6764958008109694533e-173,
        6802601037806061975e0197, 1360520207561212395e0198,

        62259110684423957791e0047, 88800290202542652011e-226,
        41010852717673354694e-221, 20505426358836677347e-221,
        66102447903809911604e0055, 35600952588064798515e0119,
        14371240869903838702e0205, 57500690832492901689e0043,
        23432630639573022093e-107, 62259110684423957791e0048,
        35620497849450218807e-306, 69658634627134074624e0200,
        99440755792436956989e-062, 55277197169490210673e0081,
        36992084760177624177e-318, 30888265282878466443e-111,
    };
    static size_t const NREC = sizeof(TYDEMAN_HARD)/(4*sizeof(double));

    ostringstream oss;
    CsvControl ctrl;
    ctrl.null = "";
    ctrl.hasNull = true;
    ctrl.quoting = "QUOTE_MINIMAL";
    ctrl.delimiter = "|";
    ctrl.escapeChar = "\\";
    ctrl.quoteChar = "'";
    ctrl.skipInitialSpace = false;
    ctrl.doubleQuote = false;
    ctrl.standardEscapes = false;
    ctrl.trailingDelimiter = false;
    ctrl.validate();

    CsvWriter writer(oss, ctrl);
    for (size_t r = 0; r < NREC; ++r) {
        for (size_t f = 0; f < 4u; ++ f) {
            writer.appendField(TYDEMAN_HARD[r*4 + f]);
        }
        writer.endRecord();
    }

    istringstream iss(oss.str());
    CsvReader reader(iss, ctrl);
    bool warned = false;
    for (size_t r = 0; r < NREC; ++r) {
        BOOST_CHECK_EQUAL(r + 1, reader.getNumRecords());
        for (size_t f = 0; f < 4u; ++ f) {
            double d1 = TYDEMAN_HARD[r*4 + f];
            double d2 = reader.get<double>(static_cast<int>(f));
            if (d1 != d2) {
                BOOST_CHECK_CLOSE(d1, d2, TWO_ULP);
                if (!warned) {
                    BOOST_WARN("The binary to decimal conversion routines of "
                               "your system do not appear to be correctly "
                               "rounded");
                }
                warned = true;
            }
        }
        BOOST_CHECK_EQUAL(reader.isDone(), false);
        reader.nextRecord();
    }
    BOOST_CHECK(reader.isDone());
}

template <typename T>
static void integerIoTest() {
    ostringstream oss;
    CsvControl ctrl;
    ctrl.null = "";
    ctrl.hasNull = true;
    ctrl.quoting = "QUOTE_MINIMAL";
    ctrl.delimiter = "\t";
    ctrl.escapeChar = "\\";
    ctrl.quoteChar = "\'";
    ctrl.skipInitialSpace = false;
    ctrl.doubleQuote = false;
    ctrl.standardEscapes = false;
    ctrl.trailingDelimiter = false;
    ctrl.validate();
    CsvWriter writer(oss, ctrl);
    for (int i = 0; i < 100; ++i) {
        writer.appendField(static_cast<T>(i));
        if (numeric_limits<T>::is_signed) {
            writer.appendField(-static_cast<T>(i));
        } else {
            writer.appendField(numeric_limits<T>::max() - static_cast<T>(i));
        }
    }
    writer.endRecord();
    // check limits, whitespace insensitivity, illegal values, overflow
    writer.appendField(numeric_limits<T>::min());
    writer.appendField(numeric_limits<T>::max());
    writer.appendField("   1");
    writer.appendField("1  ");
    writer.appendField(" 0 ");
    writer.appendField(" a ");
    writer.appendField("1 a");
    writer.appendField("1.0");
    writer.appendField("12345678901234567890123456789012345678901234567890123456789012345678901234567890");
    writer.appendField("-12345678901234567890123456789012345678901234567890123456789012345678901234567890");
    writer.endRecord();

    istringstream iss(oss.str());
    CsvReader reader(iss, ctrl);
    for (int i = 0; i < 100; ++i) {
        BOOST_CHECK_EQUAL(reader.get<T>(i*2), static_cast<T>(i));
        if (numeric_limits<T>::is_signed) {
            BOOST_CHECK_EQUAL(reader.get<T>(i*2 + 1), -static_cast<T>(i));
        } else {
            BOOST_CHECK_EQUAL(reader.get<T>(i*2 + 1),
                              numeric_limits<T>::max() - static_cast<T>(i));
        }
    }
    reader.nextRecord();
    BOOST_CHECK_EQUAL(reader.get<T>(0), numeric_limits<T>::min());
    BOOST_CHECK_EQUAL(reader.get<T>(1), numeric_limits<T>::max());
    BOOST_CHECK_EQUAL(reader.get<T>(2), static_cast<T>(1));
    BOOST_CHECK_EQUAL(reader.get<T>(3), static_cast<T>(1));
    BOOST_CHECK_EQUAL(reader.get<T>(4), static_cast<T>(0));
    for (int i = 5; i < 10; ++i) {
        BOOST_CHECK_THROW(reader.get<T>(i), Exception);
    }
    reader.nextRecord();
    BOOST_CHECK(reader.isDone());
}

BOOST_AUTO_TEST_CASE(signedCharIo)       { integerIoTest<signed char>();        }
BOOST_AUTO_TEST_CASE(unsignedCharIo)     { integerIoTest<unsigned char>();      }
BOOST_AUTO_TEST_CASE(shortIo)            { integerIoTest<short>();              }
BOOST_AUTO_TEST_CASE(unsignedShortIo)    { integerIoTest<unsigned short>();     }
BOOST_AUTO_TEST_CASE(intIo)              { integerIoTest<int>();                }
BOOST_AUTO_TEST_CASE(unsignedIntIo)      { integerIoTest<unsigned int>();       }
BOOST_AUTO_TEST_CASE(longIo)             { integerIoTest<long>();               }
BOOST_AUTO_TEST_CASE(unsignedLongIo)     { integerIoTest<unsigned long>();      }
BOOST_AUTO_TEST_CASE(longLongIo)         { integerIoTest<long long>();          }
BOOST_AUTO_TEST_CASE(unsignedLongLongIo) { integerIoTest<unsigned long long>(); }

BOOST_AUTO_TEST_CASE(boolIo) {
    ostringstream oss;
    CsvControl ctrl;
    ctrl.null = "";
    ctrl.hasNull = true;
    ctrl.quoting = "QUOTE_MINIMAL";
    ctrl.delimiter = ",";
    ctrl.escapeChar = "\\";
    ctrl.quoteChar = "'";
    ctrl.skipInitialSpace = false;
    ctrl.doubleQuote = false;
    ctrl.standardEscapes = false;
    ctrl.trailingDelimiter = false;
    ctrl.validate();
    CsvWriter writer(oss, ctrl);
    writer.appendField(true);
    writer.appendField(false);
    writer.appendField("2");
    writer.appendField("Foo");
    writer.appendField("falsely");
    writer.appendField("Ta");
    writer.appendField("Truest");
    writer.endRecord();
    writer.appendField("True  ");
    writer.appendField(" yes");
    writer.appendField("  1  ");
    writer.appendField("NO  ");
    writer.appendField(" f");
    writer.appendField(" 0  ");
    writer.endRecord();
    writer.appendField('t');
    writer.appendField('T');
    writer.appendField("true");
    writer.appendField("TRUE");
    writer.appendField('y');
    writer.appendField('Y');
    writer.appendField("yes");
    writer.appendField("YES");
    writer.endRecord();
    writer.appendField('f');
    writer.appendField('F');
    writer.appendField("false");
    writer.appendField("FALSE");
    writer.appendField('n');
    writer.appendField('N');
    writer.appendField("no");
    writer.appendField("NO");
    istringstream iss(oss.str());
    CsvReader reader(iss, ctrl);
    BOOST_CHECK_EQUAL(reader.get<char>(0), '1');
    BOOST_CHECK_EQUAL(reader.get<char>(1), '0');
    BOOST_CHECK(reader.get<bool>(0));
    BOOST_CHECK(!reader.get<bool>(1));
    for (int i = 2; i < 7; ++i) {
        BOOST_CHECK_THROW(reader.get<bool>(i), Exception);
    }
    reader.nextRecord();
    BOOST_CHECK(!reader.isDone());
    BOOST_CHECK(reader.get<bool>(0));
    BOOST_CHECK(reader.get<bool>(1));
    BOOST_CHECK(reader.get<bool>(2));
    BOOST_CHECK(!reader.get<bool>(3));
    BOOST_CHECK(!reader.get<bool>(4));
    BOOST_CHECK(!reader.get<bool>(5));
    reader.nextRecord();
    BOOST_CHECK(!reader.isDone());
    for (int i = 0; i < reader.getNumFields(); ++i) {
        BOOST_CHECK(reader.get<bool>(i));
    }
    reader.nextRecord();
    BOOST_CHECK(!reader.isDone());
    for (int i = 0; i < reader.getNumFields(); ++i) {
        BOOST_CHECK(!reader.get<bool>(i));
    }
    reader.nextRecord();
    BOOST_CHECK(reader.isDone());
}

BOOST_AUTO_TEST_CASE(charIo) {
    string characters("abcdefgHIJKLMNOP");
    ostringstream oss;
    CsvControl ctrl;
    ctrl.null = "";
    ctrl.hasNull = true;
    ctrl.quoting = "QUOTE_MINIMAL";
    ctrl.delimiter = "|";
    ctrl.escapeChar = "\\";
    ctrl.quoteChar = "'";
    ctrl.skipInitialSpace = false;
    ctrl.doubleQuote = false;
    ctrl.standardEscapes = false;
    ctrl.trailingDelimiter = false;
    ctrl.validate();

    CsvWriter writer(oss, ctrl);
    for (size_t i = 0; i < characters.size(); ++i) {
        writer.appendField(characters[i]);
    }
    writer.endRecord();
    writer.appendField("");
    writer.appendField("more than one character");
    writer.appendField('\\');
    writer.appendField('|');
    writer.appendField('\'');
    writer.endRecord();

    istringstream iss(oss.str());
    CsvReader reader(iss, ctrl);
    for (size_t i = 0; i < characters.size(); ++i) {
        BOOST_CHECK_EQUAL(reader.get<char>(static_cast<int>(i)), characters[i]);
    }
    reader.nextRecord();
    BOOST_CHECK(!reader.isDone());
    // check that retrieving empty or multi-character fields as a char throws
    BOOST_CHECK_THROW(reader.get<char>(0), Exception);
    BOOST_CHECK_THROW(reader.get<char>(1), Exception);
    BOOST_CHECK_EQUAL(reader.get<char>(2), '\\');
    BOOST_CHECK_EQUAL(reader.get<char>(3), '|');
    BOOST_CHECK_EQUAL(reader.get<char>(4), '\'');
    reader.nextRecord();
    BOOST_CHECK(reader.isDone());
}
