// -*- lsst-c++ -*-

/** 
 * @file
 * @brief Utility to dump object data from the database to binary files.
 *
 * @author $Author$
 * @version $Revision$
 * @date $Date$
 *
 * Contact: Kian-Tat Lim (ktl@slac.stanford.edu)
 *
 * @ingroup associate
 */

#ifndef __GNUC__
#  define __attribute__(x) /*NOTHING*/
#endif
static char const* SVNid __attribute__((unused)) = "$Id$";

#include <cstdlib>
#include <cerrno>
#include <fcntl.h>
#include <getopt.h>
#include <iostream>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>

#include <mysql/mysql.h>

#include "lsst/ap/Chunk.h"
#include "lsst/ap/Object.h"
#include "lsst/afw/image/Filter.h"
#include "lsst/pex/exceptions.h"
#include "lsst/daf/persistence/DbAuth.h"
#include "lsst/daf/persistence/DbStorageLocation.h"
#include "lsst/pex/policy/Policy.h"


// Command line option descriptors.
static struct option longopts[] = {
    {"minRA", required_argument, 0, 'r'},
    {"maxRA", required_argument, 0, 'R'},
    {"minDecl", required_argument, 0, 'd'},
    {"maxDecl", required_argument, 0, 'D'},
    {"url", required_argument, 0, 'u'},
    {"output", required_argument, 0, 'o'}
};

// Describe correct usage.
static void usage(char const* argv0) {
    std::cerr << "Usage: " << argv0 << " -o -d -D {-n | -r -R} [options]";
    std::cerr << std::endl << std::endl;
    std::cerr << "Required parameters:" << std::endl;
    std::cerr << "\t-o outputFile   Output filename (or base for appended chunk number)" << std::endl;
    std::cerr << "\t-d minDecl      Minimum declination [-90, 90)" << std::endl;
    std::cerr << "\t-D maxDecl      Maximum declination (-90, 90]" << std::endl;
    std::cerr << "Either -n or -r/-R must be specified:" << std::endl;
    std::cerr << "\t-n numChunks    Number of chunks in the stripe" << std::endl;
    std::cerr << "\t-r minRA        Minimum right ascension [0, 360)" << std::endl;
    std::cerr << "\t-R maxRA        Maximum right ascension (0, 360]" << std::endl;
    std::cerr << "Options:" << std::endl;
    std::cerr << "\t-u dbURL        Connection URL for database" << std::endl;
    std::cerr << "\t                (default: mysql://lsst10.ncsa.uiuc.edu:3306/test)" << std::endl;
    std::cerr << "\t-p policyFile   Filename of policy for database authentication" << std::endl;
    std::cerr << "\t                (default: none; falls back to file or environment)" << std::endl;
    std::exit(1);
}


// Queries to execute.

// Query for a single chunk.
static char const chunkQuery[] =
    "SELECT objectId, ra, decl,"
    " uVarProb, gVarProb, rVarProb, iVarProb, zVarProb, yVarProb"
    " FROM Object"
    " WHERE ra BETWEEN ? AND ? AND decl BETWEEN ? AND ?";
size_t const chunkQuerySize = sizeof(chunkQuery);

// Query for entire stripes.
// ORDER BY is required to avoid having many simultaneous open files.
static char const stripeQuery[] =
    "SELECT objectId, ra, decl,"
    " uVarProb, gVarProb, rVarProb, iVarProb, zVarProb, yVarProb"
    " FROM Object"
    " WHERE decl BETWEEN ? AND ? ORDER BY ra";
size_t const stripeQuerySize = sizeof(stripeQuery);


// Check a condition and throw an exception (with errno decoding) if failed.
static void ioAssert(bool cond, char const* msg) {
    if (cond) {
        return;
    }

    static char buf[64];
    std::snprintf(buf, sizeof(buf), "error %d", errno);
    throw lsst::pex::exceptions::Runtime(std::string(msg) + ": " + buf);
}


// Main program.
int main(int argc, char** argv) {

    // Database connection URL (with default).
    std::string url("mysql://lsst10.ncsa.uiuc.edu:3306/test");

    // Filename for output, or base if writing entire stripe.
    char const* fileBase = 0;

    // Optional name of policy file for database authentication.
    std::string policyFileName;

    // Number of chunks in this stripe.
    int numChunks = -1;

    // RA and declination minimum and maximum (initialized to bogus values).
    double raLimits[2] = {-1000.0, -1000.0};
    double declLimits[2] = {-1000.0, -1000.0};

    // Process command line arguments.
    int opt;
    while ((opt = getopt_long(argc, argv, "r:R:n:d:D:u:o:", longopts, 0)) != -1) {
        switch (opt) {
        case 'r':
            raLimits[0] = std::strtod(optarg, 0);
            break;
        case 'R':
            raLimits[1] = std::strtod(optarg, 0);
            break;
        case 'n':
            numChunks = std::strtol(optarg, 0, 10);
            break;
        case 'd':
            declLimits[0] = std::strtod(optarg, 0);
            break;
        case 'D':
            declLimits[1] = std::strtod(optarg, 0);
            break;
        case 'u':
            url = std::string(optarg);
            break;
        case 'o':
            fileBase = optarg;
            break;
        case 'p':
            policyFileName = std::string(optarg);
            break;
        default:
            usage(argv[0]);
            break;
        }
    }

    if (fileBase == 0) usage(argv[0]);
    int fileNameLen = strlen(fileBase) + 1 + 10 + 1;
    char* fileName = new char[fileNameLen];
    if (fileName == 0) {
        throw lsst::pex::exceptions::Runtime(
            "Unable to allocate memory for filename");
    }

    if (numChunks == -1) {
        // Check RA parameters.
        if (raLimits[0] < 0.0 || raLimits[1] > 360.0 ||
            raLimits[0] >= raLimits[1]) usage(argv[0]);
    }
    else if (numChunks < 1 ||
             raLimits[0] != -1000.0 || raLimits[1] != -1000.0) {
        usage(argv[0]);
    }

    if (declLimits[0] < -90.0 || declLimits[1] > 90.0 ||
        declLimits[0] > declLimits[1]) usage(argv[0]);

    argc -= optind;
    argv += optind;


    // Set up database authentication and the storage location.
    if (!policyFileName.empty()) {
        lsst::pex::policy::Policy::Ptr policyPtr(
            lsst::pex::policy::Policy::createPolicy(policyFileName));
        lsst::daf::persistence::DbAuth::setPolicy(policyPtr);
    }
    lsst::daf::persistence::DbStorageLocation dbLoc(url);

    // Set up the MySQL client library and connect to the database.
    // Use MySQL client directly for maximum performance.
    if (mysql_library_init(0, 0, 0)) {
        throw lsst::pex::exceptions::Runtime(
            "Unable to initialize mysql library");
    }
    MYSQL* db = mysql_init(0);
    if (db == 0) {
        throw lsst::pex::exceptions::Runtime(
            "Unable to create empty connection");
    }
    if (!mysql_real_connect(db, dbLoc.getHostname().c_str(),
                            dbLoc.getUsername().c_str(),
                            dbLoc.getPassword().c_str(),
                            dbLoc.getDbName().c_str(),
                            std::strtol(dbLoc.getPort().c_str(), 0, 10),
                            0, CLIENT_COMPRESS)) {
        throw lsst::pex::exceptions::Runtime(
            std::string("Unable to connect to database: ") + mysql_error(db));
    }

    // Prepare the query statement.
    MYSQL_STMT* stmt = mysql_stmt_init(db);
    if (stmt == 0) {
        throw lsst::pex::exceptions::Runtime(
            "Out of memory while preparing statement");
    }
    if (numChunks < 1) {
        if (mysql_stmt_prepare(stmt, chunkQuery, chunkQuerySize)) {
            throw lsst::pex::exceptions::Runtime(
                std::string("Unable to prepare statement: ") + mysql_error(db));
        }
    }
    else {
        if (mysql_stmt_prepare(stmt, stripeQuery, stripeQuerySize)) {
            throw lsst::pex::exceptions::Runtime(
                std::string("Unable to prepare statement: ") + mysql_error(db));
        }
    }

    // Set up the WHERE clause bindings.
    MYSQL_BIND bindArray[4];

    memset(bindArray, 0, sizeof(bindArray));

    if (numChunks < 1) {
        for (int i = 0; i < 4; ++i) {
            bindArray[i].buffer_type = MYSQL_TYPE_DOUBLE;
            bindArray[i].buffer_length = sizeof(double);
            bindArray[i].is_null = 0;
            if (i < 2) bindArray[i].buffer = &raLimits[i];
            else bindArray[i].buffer = &declLimits[i - 2];
        }
    }
    else {
        for (int i = 0; i < 2; ++i) {
            bindArray[i].buffer_type = MYSQL_TYPE_DOUBLE;
            bindArray[i].buffer_length = sizeof(double);
            bindArray[i].is_null = 0;
            bindArray[i].buffer = &declLimits[i];
        }
    }


    // Set up the result bindings.

    // Retrieved object.
    lsst::ap::Object obj;
    // Result binding array.
    MYSQL_BIND resultArray[3 + lsst::afw::image::Filter::NUM_FILTERS];
    // Null flags for object fields.
    my_bool isNull[3 + lsst::afw::image::Filter::NUM_FILTERS];
    // Error flags for object fields.
    my_bool error[3 + lsst::afw::image::Filter::NUM_FILTERS];

    // Initialize object to junk values to help detect bugs.
    obj._objectId = 0xcafefeeddeadbeefLL;
    obj._ra = -100.0;
    obj._decl = -1234567890.0;
    for (int i = 0; i < lsst::afw::image::Filter::NUM_FILTERS; ++i) {
        obj._varProb[i] = -1;
    }
    for (int i = 0; i < 3 + lsst::afw::image::Filter::NUM_FILTERS; ++i) {
        isNull[i] = false;
    }

    memset(resultArray, 0, sizeof(resultArray));

    // Set each result binding.
    resultArray[0].buffer_type = MYSQL_TYPE_LONGLONG;
    resultArray[0].buffer = &(obj._objectId);
    resultArray[0].buffer_length = sizeof(obj._objectId);
    resultArray[0].is_null = &isNull[0];
    resultArray[0].is_unsigned = false;
    resultArray[0].error = &error[0];

    resultArray[1].buffer_type = MYSQL_TYPE_DOUBLE;
    resultArray[1].buffer = &(obj._ra);
    resultArray[1].buffer_length = sizeof(obj._ra);
    resultArray[1].is_null = &isNull[1];
    resultArray[1].is_unsigned = false;
    resultArray[1].error = &error[1];

    resultArray[2].buffer_type = MYSQL_TYPE_DOUBLE;
    resultArray[2].buffer = &(obj._decl);
    resultArray[2].buffer_length = sizeof(obj._decl);
    resultArray[2].is_null = &isNull[2];
    resultArray[2].is_unsigned = false;
    resultArray[2].error = &error[2];

    for (int i = 0; i < lsst::afw::image::Filter::NUM_FILTERS; ++i) {
        resultArray[3 + i].buffer_type = MYSQL_TYPE_SHORT;
        resultArray[3 + i].buffer = &(obj._varProb[i]);
        resultArray[3 + i].buffer_length = sizeof(obj._varProb[i]);
        resultArray[3 + i].is_null = &isNull[3 + i];
        resultArray[3 + i].is_unsigned = false;
        resultArray[3 + i].error = &error[3 + i];
    }


    // Bind and execute the SQL statement.
    if (mysql_stmt_bind_param(stmt, bindArray)) {
        throw lsst::pex::exceptions::Runtime("Unable to bind parameters");
    }
    if (mysql_stmt_execute(stmt)) {
        throw lsst::pex::exceptions::Runtime("Unable to execute statement");
    }
    if (mysql_stmt_bind_result(stmt, resultArray)) {
        throw lsst::pex::exceptions::Runtime("Unable to bind results");
    }


    // Set up (most of) a chunk header.
    lsst::ap::BinChunkHeader header;
    header._recordSize = sizeof(obj);


    // Fetch rows from the database and write the resulting objects.

    // Error status return.
    int err = 0;
    // File descriptor for output file.
    int fd = -1;
    // Current chunk number.
    int chunkNum = -1;
    // Upper RA boundary of current chunk.
    double chunkBoundary = -1.0;

    while ((err = mysql_stmt_fetch(stmt)) == 0) {

        // Check for nulls.
        for (int i = 0; i < 9; ++i) {
            if (isNull[i]) {
                throw lsst::pex::exceptions::Runtime(
                    "Unexpected null value found");
            }
        }

        // Check for new chunk.
        if (obj._ra > chunkBoundary) {

            if (fd != -1) {
                // Rewrite the header, now that we know how many rows we have.
                off_t pos = lseek(fd, 0, SEEK_SET);
                if (pos == static_cast<off_t>(-1)) {
                    throw lsst::pex::exceptions::Runtime(
                        "Unable to rewind output file");
                }
                ssize_t bytes = write(fd, &header, sizeof(header));
                ioAssert(bytes == sizeof(header),
                         "Unable to write updated header");

                // Close the output file.
                err = close(fd);
                ioAssert(err == 0, "Unable to close output file");
            }

            if (numChunks < 1) {
                strcpy(fileName, fileBase);
                chunkBoundary = 1000.0;
            }
            else {
                chunkNum = static_cast<int>(floor(obj._ra * numChunks / 360.0));
                snprintf(fileName, fileNameLen, "%s.%d", fileBase, chunkNum);
                chunkBoundary = (360.0 * (chunkNum + 1)) / numChunks;
            }

            // Open the new output file.
            fd = open(fileName, O_WRONLY | O_CREAT | O_TRUNC, 0666);
            ioAssert(fd >= 0, "Unable to open output file");

            // Write the output file header.
            header._numRecords = 0;
            ssize_t bytes = write(fd, &header, sizeof(header));
            ioAssert(bytes == sizeof(header), "Unable to write header");

        }

        // Write the object and increment the count.
        ssize_t bytes = write(fd, &obj, sizeof(obj));
        ioAssert(bytes == sizeof(obj), "Unable to write object");
        ++header._numRecords;
    }
    if (err == 1) {
        throw lsst::pex::exceptions::Runtime(
            std::string("Error while fetching: ") + mysql_error(db));
    }
    else if (err == MYSQL_DATA_TRUNCATED) {
        throw lsst::pex::exceptions::Runtime(
            "Error while fetching: data truncated");
    }
    else if (err != MYSQL_NO_DATA) {
        throw lsst::pex::exceptions::Runtime("Error while fetching: unknown");
    }

    // Finish the final chunk.
    if (fd != -1) {
        // Rewrite the header, now that we know how many rows we have.
        off_t pos = lseek(fd, 0, SEEK_SET);
        if (pos == static_cast<off_t>(-1)) {
            throw lsst::pex::exceptions::Runtime(
                "Unable to rewind output file");
        }
        ssize_t bytes = write(fd, &header, sizeof(header));
        ioAssert(bytes == sizeof(header), "Unable to write updated header");

        // Close the output file.
        err = close(fd);
        ioAssert(err == 0, "Unable to close output file");
    }

    // Close down MySQL.
    if (mysql_stmt_close(stmt)) {
        throw lsst::pex::exceptions::Runtime("Error while closing statement");
    }
    mysql_close(db);
    mysql_library_end();

    delete[] fileName;

    return 0;
}
