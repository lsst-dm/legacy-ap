// -*- lsst-c++ -*-

/**
 * @file
 * @brief   Helper functions for storing file name/line number information in an ExceptionStack
 *
 * @ingroup associate
 */

#include <iostream>
#include <sstream>

#include <boost/any.hpp>

#include <lsst/daf/base/DataProperty.h>

#include <lsst/ap/Exceptions.h>


namespace lsst {
namespace ap {


/**
 * Adds a DataProperty named "origin" (containing a string describing the location
 * from which @a ex is to be thrown) to the given exception object.
 */
LSST_AP_API void setOrigin(
    lsst::pex::exceptions::ExceptionStack & ex,
    int  const         line,
    char const * const file,
    char const * const function
) throw() {
    if (file != 0) {
        try {
            std::ostringstream oss;
            oss << file << ':' << line;
            if (function != 0 && function[0] != 0) {
                oss << ": in function '" << function << '\'';
            }
            ex.getLast()->addProperty(
                lsst::daf::base::DataProperty("origin", boost::any(oss.str()))
            );
        } catch (...) {}
    }
}


/**
 * Adds a DataProperty named "error" containing a system error code
 * value (as obtained from e.g. errno) to the given exception object.
 */
LSST_AP_API void setSystemErrorCode(lsst::pex::exceptions::ExceptionStack & ex, int const errorCode) throw() {
    try {
        ex.getLast()->addProperty(
            lsst::daf::base::DataProperty("error", boost::any(errorCode))
        );
    } catch (...) {}
}


}} // end of namespace lsst::ap

