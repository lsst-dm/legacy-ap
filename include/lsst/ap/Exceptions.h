// -*- lsst-c++ -*-

/**
 * @file
 * @brief   Association pipeline exception classes and macros for throwing them with
 *          file name/line number information.
 *
 * @ingroup associate
 */

#ifndef LSST_AP_EXCEPTIONS_H
#define LSST_AP_EXCEPTIONS_H

#include <lsst/pex/exceptions.h>

#include "Common.h"
#include "CustomExceptions.h"


namespace lsst {
namespace ap {

using lsst::pex::exceptions::DomainError;
using lsst::pex::exceptions::InvalidParameter;
using lsst::pex::exceptions::LengthError;
using lsst::pex::exceptions::OutOfRange;
using lsst::pex::exceptions::Runtime;
using lsst::pex::exceptions::NotFound;
using lsst::pex::exceptions::Memory;


#define LSST_AP_JOIN_IMPL(a,b)  a##b
#define LSST_AP_JOIN(a,b)       LSST_AP_JOIN_IMPL(a,b)

/**
 * @def LSST_AP_THROW_TERSE
 *
 * Convenience macro for throwing an exception object of type @a Ex without an error message.
 * The line number, file name, and enclosing function name of the macro call are recorded in the
 * resulting exception.
 */
#define LSST_AP_THROW_TERSE(Ex) \
    do { \
        Ex LSST_AP_JOIN(exs_, __LINE__)(""); \
        setOrigin(LSST_AP_JOIN(exs_, __LINE__), __LINE__, __FILE__, LSST_AP_FUNC); \
        throw LSST_AP_JOIN(exs_, __LINE__); \
    } while (false)

/**
 * @def LSST_AP_THROW
 *
 * Convenience macro for throwing an exception object of type @a Ex with an error message @a msg.
 * The line number, file name, and enclosing function name of the macro call are recorded in the
 * resulting exception.
 */
#define LSST_AP_THROW(Ex, msg) \
    do { \
       Ex LSST_AP_JOIN(exs_, __LINE__)(msg); \
       setOrigin(LSST_AP_JOIN(exs_, __LINE__), __LINE__, __FILE__, LSST_AP_FUNC); \
       throw LSST_AP_JOIN(exs_, __LINE__); \
    } while (false)

/**
 * @def LSST_AP_THROW_ERR
 *
 * Convenience macro for throwing an exception object of type @a Ex with an error message @a msg.
 * The line number, file name, and enclosing function name of the macro call are recorded in the
 * resulting exception, as is an error code @a err.
 */
#define LSST_AP_THROW_ERR(Ex, msg, err) \
    do { \
       Ex LSST_AP_JOIN(exs_, __LINE__)(msg); \
       setOrigin(LSST_AP_JOIN(exs_, __LINE__), __LINE__, __FILE__, LSST_AP_FUNC); \
       setSystemErrorCode(LSST_AP_JOIN(exs_, __LINE__), err); \
       throw LSST_AP_JOIN(exs_, __LINE__); \
    } while (false)


LSST_AP_API void setOrigin(
    lsst::pex::exceptions::ExceptionStack & ex,
    int  const         line,
    char const * const file,
    char const * const function
) throw();


LSST_AP_API void setSystemErrorCode(lsst::pex::exceptions::ExceptionStack & ex, int const errorCode) throw();


}} // end of namespace lsst::ap

#endif // LSST_AP_EXCEPTIONS_H

