undefine(`format')dnl ' Stop m4 expanding format
// -*- lsst-c++ -*-

/**
 * @file
 * @brief   Association pipeline specific exception classes.
 *
 * This file is machine generated from CustomExceptions.h.m4. Please do not edit directly.
 *
 * @ingroup associate
 */

#ifndef LSST_AP_CUSTOM_EXCEPTIONS_H
#define LSST_AP_CUSTOM_EXCEPTIONS_H

#include <lsst/pex/exceptions.h>


namespace lsst {
namespace ap {
define(LSST_AP_NEW_EXCEPTION,
`/** @brief  Exception class thrown when $2. */
class LSST_AP_API $1 : public lsst::pex::exceptions::ExceptionStack  {
public:

    $1(std::string const & comment) throw() :
        lsst::pex::exceptions::ExceptionStack(std::string("$1"), comment)
    {
        lsst::pex::exceptions::ExceptionData eNode("Node$1");
        this->addExceptionData(eNode);
    }

    $1(boost::format const & comment) throw() :
        lsst::pex::exceptions::ExceptionStack(std::string("$1"), comment.str())
    {
        lsst::pex::exceptions::ExceptionData eNode("Node$1");
        this->addExceptionData(eNode);
    }

    $1(const $1 & orig) throw() : lsst::pex::exceptions::ExceptionStack(orig ) {}

    $1 & operator<<(lsst::pex::exceptions::ExceptionData rhs) throw() {
        this->getStack()->addProperty(rhs.getExceptionData());
        return *this;
    }

    $1 & operator<<(lsst::daf::base::DataProperty const & rhs) throw() {
        this->getLast()->addProperty(rhs);
        return *this;
    }

    $1 & operator<<(lsst::daf::base::DataProperty::PtrType const & rhs) throw() {
        this->getLast()->addProperty(rhs);
        return *this;
    }

    ~$1() throw() {}

    virtual char const * getPythonModule() const { return "lsst.ap.exceptions"; }
    virtual char const * getPythonClass()  const { return "Lsst$1"; }

}')

LSST_AP_NEW_EXCEPTION(Timeout, an operation times out);


LSST_AP_NEW_EXCEPTION(IoError, an input/output error occurs);


}} // end of namespace lsst::ap

#endif // LSST_AP_CUSTOM_EXCEPTIONS_H

