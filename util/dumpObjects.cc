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
#include <iostream>
#include <cmath>
#include <unistd.h>

#include "mysql/mysql.h"

#include "boost/format.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/program_options.hpp"

#include "lsst/ap/Chunk.h"
#include "lsst/ap/Object.h"
#include "lsst/afw/image/Filter.h"
#include "lsst/pex/exceptions.h"
#include "lsst/daf/persistence/DbAuth.h"
#include "lsst/pex/policy/Policy.h"

using namespace boost::program_options;

namespace ex = lsst::pex::exceptions;

using lsst::pex::policy::Policy;
using lsst::daf::persistence::DbAuth;
using lsst::afw::image::Filter;
using lsst::ap::BinChunkHeader;
using lsst::ap::Object;


// Queries to execute.

// Query pattern for a single chunk.
static char const chunkQuery[] =
    "SELECT objectId, ra, decl"
    " FROM %1%"
    " WHERE decl BETWEEN ? AND ? AND ra BETWEEN ? AND ?";

// Query pattern for entire stripes.
// ORDER BY is required to avoid having many simultaneous open files.
static char const stripeQuery[] =
    "SELECT objectId, ra, decl"
    " FROM %1%"
    " WHERE decl BETWEEN ? AND ? ORDER BY ra";


// Check a condition and throw an exception (with errno decoding) if failed.
static void ioAssert(bool cond, char const* msg) {
    if (!cond) {
        throw LSST_EXCEPT(ex::RuntimeErrorException,
                          (boost::format("%1%: error %2%") % msg % errno).str());
    }
}

// Check whether two conflicting options were specified and complain if so
static void conflictingOptions(
    options_description const & desc,
    variables_map const & vm,
    char const * const opt1,
    char const * const opt2
) {
    if (vm.count(opt1) && !vm[opt1].defaulted() && vm.count(opt2) && !vm[opt2].defaulted()) {
        std::cerr << "Options --" << opt1 << " and --" << opt2 << " conflict" <<
                     std::endl << desc << std::endl;
        std::exit(1);
    }
}

static void requiredOption(
    options_description const & desc,
    variables_map const & vm,
    char const * const opt
) {
    if (vm.count(opt) != 1 && !vm[opt].defaulted()) {
        std::cerr << "--" << opt << " is a required parameter" <<
                     std::endl << desc << std::endl;
        std::exit(1);
    }
}

// Main program.
int main(int argc, char** argv) {

    // Number of chunks in this stripe.
    int numChunks = -1;
    // RA and declination minimum and maximum (initialized to bogus values)
    double minRa = -1000.0, maxRa = -1000.0;
    double minDecl = -1000.0, maxDecl = -1000.0;

    std::string fileBase;       // Filename for output, or base if writing entire stripe
    std::string policyFileName; // Optional name of policy file for database authentication
    std::string host;           // MySQL server hostname
    std::string port;           // MySQL server port number
    std::string database;       // MySQL database
    std::string table;          // Object table name

    options_description generalOptions("General options");
    generalOptions.add_options()
        ("help,h", "print usage help");

    options_description requiredOptions("Required (either -n OR -r/-R must be specified)");
    requiredOptions.add_options()
        ("min-decl,d", value(&minDecl), "Minimum declination [-90, 90)")
        ("max-decl,D", value(&maxDecl), "Maximum declination (-90, 90]")
        ("min-ra,r", value(&minRa), "Minimum right ascension [0, 360)")
        ("max-ra,R", value(&maxRa), "Maximum right ascension (0, 360]")
        ("num-chunks,n", value(&numChunks), "Number of chunks in the stripe")
        ("output,o", value(&fileBase), "Output filename (or base for appended chunk number)");

    options_description databaseOptions("Database options");
    databaseOptions.add_options()
        ("host,H", value(&host)->default_value("lsst10.ncsa.uiuc.edu"),
            "MySQL database server hostname")
        ("port,P", value(&port)->default_value("3306"),
            "MySQL database server port number")
        ("database,b", value(&database)->default_value("DC3a_catalogs"),
            "Database name")
        ("table,t", value(&table)->default_value("CFHTLSObject"),
            "Object table name")
        ("policy,p", value(&policyFileName),
            "Filename of policy for database authentication");

    options_description desc;
    desc.add(generalOptions).add(requiredOptions).add(databaseOptions);
    variables_map vm;
    store(parse_command_line(argc, argv, desc), vm);
    notify(vm);
    conflictingOptions(desc, vm, "num-chunks", "min-ra");
    conflictingOptions(desc, vm, "num-chunks", "max-ra");
    requiredOption(desc, vm, "min-decl");
    requiredOption(desc, vm, "max-decl");
    requiredOption(desc, vm, "output");

    if (vm.count("help")) {
        std::cout << desc << std::endl;
        return 0;
    }
    if (vm.count("num-chunks")) {
        // Check chunk count
        if (numChunks < 1) {
            std::cerr << "Number of chunks must be positive" <<
                         std::endl << desc << std::endl;
            return 1;
        }
    } else if (vm.count("min-ra") && vm.count("max-ra")) {
        // Check right ascensions limits
        if (minRa < 0.0 || maxRa > 360.0 || minRa >= maxRa) {
            std::cerr << "Illegal right ascension limits" <<
                         std::endl << desc << std::endl;
            return 1;
        }
    } else {
        std::cerr << "Right ascension limits or chunk count must be specified" <<
                     std::endl << desc << std::endl;
        return 1;
    }
    // Check declination limits
    if (minDecl < -90.0 || maxDecl > 90.0 || minDecl >= maxDecl) {
        std::cerr << "Illegal declination limits" <<
                     std::endl << desc << std::endl;
        return 1;
    }

    // Set up database authentication and the storage location.
    if (!policyFileName.empty()) {
        Policy::Ptr policy(new Policy(policyFileName));
        DbAuth::setPolicy(policy);
    }
    if (!DbAuth::available(host, port)) {
        std::cerr << "Database credentials for " << host << ":" << port << " not found." <<
                     " Check ~/.lsst/db-auth.paf or specify an alternate host." << std::endl;
    }

    // Set up the MySQL client library and connect to the database.
    // Use MySQL client directly for maximum performance.
    if (mysql_library_init(0, 0, 0)) {
        throw LSST_EXCEPT(ex::RuntimeErrorException, "Unable to initialize mysql library");
    }
    MYSQL* db = mysql_init(0);
    if (db == 0) {
        throw LSST_EXCEPT(ex::RuntimeErrorException, "Unable to create empty connection");
    }
    if (!mysql_real_connect(db, host.c_str(),
                            DbAuth::username(host, port).c_str(),
                            DbAuth::password(host, port).c_str(),
                            database.c_str(),
                            boost::lexical_cast<long>(port),
                            0, CLIENT_COMPRESS)) {
        throw LSST_EXCEPT(ex::RuntimeErrorException, 
                          std::string("Unable to connect to database: ") + mysql_error(db));
    }

    // Prepare the query statement.
    MYSQL_STMT* stmt = mysql_stmt_init(db);
    if (stmt == 0) {
        throw LSST_EXCEPT(ex::RuntimeErrorException, "Out of memory while preparing statement");
    }
    std::string query = (boost::format(numChunks < 1 ? chunkQuery : stripeQuery) % table).str();
    if (mysql_stmt_prepare(stmt, query.c_str(), query.size())) {
        throw LSST_EXCEPT(ex::RuntimeErrorException, 
                          std::string("Unable to prepare statement: ") + mysql_error(db));
    }

    // Set up the WHERE clause bindings.
    MYSQL_BIND bindArray[4];
    std::memset(bindArray, 0, sizeof(bindArray));
    for (int i = 0; i < 4; ++i) {
        bindArray[i].buffer_type = MYSQL_TYPE_DOUBLE;
        bindArray[i].buffer_length = sizeof(double);
        bindArray[i].is_null = 0;
    }
    bindArray[0].buffer = &minDecl;
    bindArray[1].buffer = &maxDecl;
    if (numChunks < 1) {
        bindArray[2].buffer = &minRa;
        bindArray[3].buffer = &maxRa;
    }

    // Set up the result bindings.
    Object obj;                                      // Retrieved object
    MYSQL_BIND resultArray[7 + Filter::NUM_FILTERS]; // Result binding array
    my_bool isNull[7 + Filter::NUM_FILTERS];         // Null flags for object fields
    my_bool error[7 + Filter::NUM_FILTERS];          // Error flags for object fields

    // Initialize object to junk values to help detect bugs.
    obj._objectId = 0xcafefeeddeadbeefLL;
    obj._ra = -100.0;
    obj._decl = -1234567890.0;

    // These fields aren't yet available in DC3a - set
    // proper motions and variability probabilities to zero
    obj._muRa = 0.0;
    obj._muDecl = 0.0;
    obj._parallax = 0.0;
    obj._radialVelocity = 0.0;
    for (int i = 0; i < Filter::NUM_FILTERS; ++i) {
        obj._varProb[i] = 0;
    }
    for (int i = 0; i < 7 + Filter::NUM_FILTERS; ++i) {
        isNull[i] = false;
    }
    std::memset(resultArray, 0, sizeof(resultArray));

    // Set each result binding.
    resultArray[0].buffer_type = MYSQL_TYPE_LONGLONG;
    resultArray[0].buffer = &(obj._objectId);
    resultArray[0].buffer_length = sizeof(obj._objectId);
    resultArray[0].is_null = &isNull[0];
    resultArray[0].is_unsigned = false;
    resultArray[0].error = &error[0];

    for (int i = 1; i < 3; ++i) {
        resultArray[i].buffer_type = MYSQL_TYPE_DOUBLE;
        resultArray[i].buffer_length = sizeof(double);
        resultArray[i].is_null = &isNull[i];
        resultArray[i].is_unsigned = false;
        resultArray[i].error = &error[i];
    }
    resultArray[1].buffer = &(obj._ra);
    resultArray[2].buffer = &(obj._decl);

    // Bind and execute the SQL statement.
    if (mysql_stmt_bind_param(stmt, bindArray)) {
        throw LSST_EXCEPT(ex::RuntimeErrorException, "Unable to bind parameters");
    }
    if (mysql_stmt_execute(stmt)) {
        throw LSST_EXCEPT(ex::RuntimeErrorException, "Unable to execute statement");
    }
    if (mysql_stmt_bind_result(stmt, resultArray)) {
        throw LSST_EXCEPT(ex::RuntimeErrorException, "Unable to bind results");
    }

    // Set up (most of) a chunk header.
    lsst::ap::BinChunkHeader header;
    header._recordSize = sizeof(obj);

    // Fetch rows from the database and write the resulting objects.

    std::string fileName;        // Filename for output
    int err = 0;                 // Error status return
    int fd = -1;                 // File descriptor for output file
    int chunkNum = -1;           // Current chunk number
    double chunkBoundary = -1.0; // Upper RA boundary of current chunk

    while ((err = mysql_stmt_fetch(stmt)) == 0) {
        // Check for nulls.
        for (int i = 0; i < 7 + Filter::NUM_FILTERS; ++i) {
            if (isNull[i]) {
                throw LSST_EXCEPT(ex::RuntimeErrorException, "Unexpected null value found");
            }
        }

        // Check for new chunk.
        if (obj._ra > chunkBoundary) {
            if (fd != -1) {
                // Rewrite the header, now that we know how many rows we have.
                off_t pos = lseek(fd, 0, SEEK_SET);
                if (pos == static_cast<off_t>(-1)) {
                    throw LSST_EXCEPT(ex::RuntimeErrorException, "Unable to rewind output file");
                }
                ssize_t bytes = write(fd, &header, sizeof(header));
                ioAssert(bytes == sizeof(header), "Unable to write updated header");

                // Close the output file.
                err = close(fd);
                ioAssert(err == 0, "Unable to close output file");
            }

            if (numChunks < 1) {
                fileName = fileBase;
                chunkBoundary = 1000.0;
            } else {
                chunkNum = static_cast<int>(floor(obj._ra * numChunks / 360.0));
                fileName = (boost::format("%1%_%2%.chunk") % fileBase % chunkNum).str();
                chunkBoundary = (360.0 * (chunkNum + 1)) / numChunks;
            }

            // Open the new output file.
            fd = open(fileName.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0666);
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
        throw LSST_EXCEPT(ex::RuntimeErrorException,
                          std::string("Error while fetching: ") + mysql_error(db));
    } else if (err == MYSQL_DATA_TRUNCATED) {
        throw LSST_EXCEPT(ex::RuntimeErrorException, "Error while fetching: data truncated");
    } else if (err != MYSQL_NO_DATA) {
        throw LSST_EXCEPT(ex::RuntimeErrorException, "Error while fetching: unknown");
    }

    // Finish the final chunk.
    if (fd != -1) {
        // Rewrite the header, now that we know how many rows we have.
        off_t pos = lseek(fd, 0, SEEK_SET);
        if (pos == static_cast<off_t>(-1)) {
            throw LSST_EXCEPT(ex::RuntimeErrorException, "Unable to rewind output file");
        }
        ssize_t bytes = write(fd, &header, sizeof(header));
        ioAssert(bytes == sizeof(header), "Unable to write updated header");

        // Close the output file.
        err = close(fd);
        ioAssert(err == 0, "Unable to close output file");
    }

    // Close down MySQL.
    if (mysql_stmt_close(stmt)) {
        throw LSST_EXCEPT(ex::RuntimeErrorException, "Error while closing statement");
    }
    mysql_close(db);
    mysql_library_end();

    return 0;
}
