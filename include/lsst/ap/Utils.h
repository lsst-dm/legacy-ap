// -*- lsst-c++ -*-
//
//##====----------------                                ----------------====##/
//
//! \file   Utils.h
//! \brief  Miscellaneous helper methods.
//
//##====----------------                                ----------------====##/

#ifndef LSST_AP_UTILS_H
#define LSST_AP_UTILS_H

#include <boost/any.hpp>
#include <boost/format.hpp>
#include <boost/numeric/conversion/converter.hpp>
#include <boost/static_assert.hpp>
#include <boost/type_traits/is_integral.hpp>

#include <lsst/mwi/data/DataProperty.h>
#include <lsst/mwi/policy/Policy.h>

#include <lsst/ap/Common.h>


namespace lsst {
namespace ap {


LSST_AP_API lsst::mwi::data::DataProperty::PtrType extractRequired(
    lsst::mwi::data::DataProperty::PtrType const & properties,
    std::string                            const & key
);

LSST_AP_API std::string const getTableName(
    lsst::mwi::policy::Policy::Ptr         const & policy,
    lsst::mwi::data::DataProperty::PtrType const & properties
);

LSST_AP_API std::string const getTableTemplateName(
    lsst::mwi::policy::Policy::Ptr         const & policy,
    lsst::mwi::data::DataProperty::PtrType const & properties
);


/*!
    Extracts an integer of the specified type from the given boost::any. The extraction
    will succeed if and only if \a v contains an integer value (of built-in C++ type)
    that can be converted to an integer of the desired type without overflow.
 */
template <typename Target>
Target anyToInteger(boost::any const & v) {

    using namespace boost::numeric;

    BOOST_STATIC_ASSERT(boost::is_integral<Target>::value);

    std::type_info const & type = v.type();
    if (type == typeid(bool)) {
        typedef converter<Target, bool> Converter;
        return Converter::convert(boost::any_cast<bool>(v));
    } else if (type == typeid(char)) {
        typedef converter<Target, char> Converter;
        return Converter::convert(boost::any_cast<char>(v));
    } else if (type == typeid(signed char)) {
        typedef converter<Target, signed char> Converter;
        return Converter::convert(boost::any_cast<signed char>(v));
    } else if (type == typeid(unsigned char)) {
        typedef converter<Target, unsigned char> Converter;
        return Converter::convert(boost::any_cast<unsigned char>(v));
    } else if (type == typeid(short)) {
        typedef converter<Target, short> Converter;
        return Converter::convert(boost::any_cast<short>(v));
    } else if (type == typeid(unsigned short)) {
        typedef converter<Target, unsigned short> Converter;
        return Converter::convert(boost::any_cast<unsigned short>(v));
    } else if (type == typeid(int)) {
        typedef converter<Target, int> Converter;
        return Converter::convert(boost::any_cast<int>(v));
    } else if (type == typeid(unsigned int)) {
        typedef converter<Target, unsigned int> Converter;
        return Converter::convert(boost::any_cast<unsigned int>(v));
    } else if (type == typeid(long)) {
        typedef converter<Target, long> Converter;
        return Converter::convert(boost::any_cast<long>(v));
    } else if (type == typeid(unsigned long)) {
        typedef converter<Target, unsigned long> Converter;
        return Converter::convert(boost::any_cast<unsigned long>(v));
    } else if (type == typeid(long long)) {
        typedef converter<Target, long long> Converter;
        return Converter::convert(boost::any_cast<long long>(v));
    } else if (type == typeid(unsigned long long)) {
        typedef converter<Target, unsigned long long> Converter;
        return Converter::convert(boost::any_cast<unsigned long long>(v));
    }
    throw boost::bad_any_cast();
}


}} // end of namespace lsst::ap

#endif // LSST_AP_UTILS_H

