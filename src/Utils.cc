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
 

/**
 * @file
 * @brief   Implementation of miscellaneous helper functions
 *
 * @ingroup ap
 */

#include <cerrno>
#include <sys/stat.h>
#include <unistd.h>

#include "lsst/pex/exceptions.h"

#include "lsst/ap/Utils.h"

namespace ex = lsst::pex::exceptions;


/**
 * Ensure that all directories along a path exist, creating them if necessary.
 *
 * @param[in] name  Pathname to file to be created
 */
void lsst::ap::verifyPathName(std::string const & name) {
    // Get the directory by stripping off anything after the last slash.
    std::string::size_type pos = name.find_last_of('/');
    if (pos == std::string::npos) return;
    std::string dirName = name.substr(0, pos);

    // Check to see if the directory exists.
    struct stat buf;
    int ret = ::stat(dirName.c_str(), &buf);

    if (ret == -1 && errno == ENOENT) {
        // It doesn't; check its parent and then create it.
        verifyPathName(dirName);

        ret = ::mkdir(dirName.c_str(), 0777);
        if (ret == -1) {
            throw LSST_EXCEPT(ex::IoErrorException, "Error creating directory " + dirName);
        }
    }
    else if (ret == -1) {
        // We couldn't read the (existing) directory for some reason.
        throw LSST_EXCEPT(ex::IoErrorException, "Error searching for directory " + dirName);
    }
    else if (!S_ISDIR(buf.st_mode)) {
        // It's not a directory.
        throw LSST_EXCEPT(ex::IoErrorException, "Non-directory in path " + dirName);
    }
}

