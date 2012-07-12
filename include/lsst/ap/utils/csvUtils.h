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
  * @brief  CSV utilities.
  * @author Serge Monkewitz
  */
#ifndef LSST_AP_UTILS_CSVUTILS_H
#define LSST_AP_UTILS_CSVUTILSL_H

#include <string>
#include <vector>

#include "lsst/pex/config.h"
#include "lsst/afw/table/Catalog.h"
#include "CsvControl.h"


namespace lsst { namespace ap { namespace utils {

/// @brief Catalog to CSV conversion parameters.
struct CsvConversionControl {
    CsvConversionControl();
    ~CsvConversionControl();

    LSST_CONTROL_FIELD(flagsAsBits, bool,
        "If true, then flag columns are mapped to 1 bit CSV columns. "
        "Otherwise, groups of 63 flag columns are packed into signed "
        "64 bit integer values. True by default.");

    LSST_CONTROL_FIELD(ignoreFields, std::vector<std::string>,
        "List of field names that should not appear in the CSV output.");

    LSST_CONTROL_FIELD(nullableIntegers, std::vector<std::string>,
        "List of names corresponding to nullable integer fields. "
        "A value of 0 in any of these fields will be mapped to a "
        "NULL during CSV conversion.");

    LSST_CONTROL_FIELD(canonicalFlags, std::vector<std::string>,
        "An order-sensitive list of flag field names corresponding to "
        "canonical flags. Ignored if flagsAsBits is set. Otherwise, these "
        "flags are packed into a separate set of 64 bit flag columns (63 "
        "flags per column), with the i-th entry of canonicalFlags "
        "corresponding to bit i % 63 of column i / 63.");
};


/// @brief Convert an afw Catalog to a CSV file.
void writeCsv(
    lsst::afw::table::BaseCatalog const & catalog, ///< @param[in] catalog    Catalog to convert.
    CsvConversionControl const & cnvControl,       ///< @param[in] cnvControl Conversion parameters.
    CsvControl const & csvControl,                 ///< @param[in] csvControl CSV dialect.
    std::string const & csvFile,                   ///< @param[in] csvFile    Name of file to write to
    bool truncate,                                 ///< @param[in] truncate   Truncate csvFile if it already exists?
    bool append                                    ///< @param[in] append     Append to csvFile if it exists and truncate is false?
);

}}} // namespace lsst::ap::utils

#endif // LSST_AP_UTILS_CSVUTILS_H
