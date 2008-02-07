// -*- lsst-c++ -*-

/**
 * @file
 * @brief   Compilation unit for boost::test generated main function
 *
 * @ingroup associate
 */

// Use boost::test in header-only mode
#include <boost/version.hpp>
#if BOOST_VERSION < 103400
#   include <boost/test/included/unit_test_framework.hpp>
#   define BOOST_AUTO_TEST_MAIN "Association Pipeline Unit Tests"
#   include <boost/test/auto_unit_test.hpp>
#else
#   define BOOST_TEST_MODULE Association Pipeline Unit Tests
#   include <boost/test/included/unit_test.hpp>
#endif
#include <boost/test/unit_test_monitor.hpp>

#include <iostream>

#include <lsst/ap/Exceptions.h>


using namespace boost::unit_test;

namespace {

using lsst::mwi::exceptions::ExceptionStack;

void translator(ExceptionStack const & ex) {
    std::cerr << '\n' << ex.what() << '\n'
              << const_cast<ExceptionStack &>(ex).getStack()->toString("", true)
              << '\n' << std::endl;
    throw boost::execution_exception(
        boost::execution_exception::cpp_exception_error,
        "Caught lsst::mwi::exceptions::ExceptionStack"
    );
}

struct RegisterExceptionTranslators { RegisterExceptionTranslators(); };

RegisterExceptionTranslators::RegisterExceptionTranslators() {
    unit_test_monitor_t::instance().register_exception_translator<ExceptionStack>(&translator);
}

static RegisterExceptionTranslators translators;

}

