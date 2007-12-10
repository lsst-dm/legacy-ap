// -*- lsst-c++ -*-
//
//##====----------------                                ----------------====##/
//
//! \file   Utils.cc
//! \brief  Implementation of miscellaneous helper functions
//
//##====----------------                                ----------------====##/

#include <lsst/ap/Exceptions.h>
#include <lsst/ap/Utils.h>


namespace lsst {
namespace ap {


LSST_AP_API lsst::mwi::data::DataProperty::PtrType extractRequired(
    lsst::mwi::data::DataProperty::PtrType const & properties,
    std::string                            const & key
) {
    lsst::mwi::data::DataProperty::PtrType dp = properties->findUnique(key);
    if (!dp) {
        LSST_AP_THROW(NotFound, boost::format("\"%1%\" property not found") % key);
    }
    return dp;
}


static char const * const sDefaultVisitNamePat = "_visit%1%";


static int64_t extractVisitId(lsst::mwi::data::DataProperty::PtrType const & properties) {
    if (!properties) {
        LSST_AP_THROW(InvalidParameter, "null DataProperty");
    }
    lsst::mwi::data::DataProperty::PtrType dp = extractRequired(properties, "visitId");
    int64_t visitId = anyToInteger<int64_t>(dp->getValue());
    if (visitId < 0) {
        LSST_AP_THROW(Runtime, "\"visitId\" property value is negative");
    }
    return visitId;
}


static std::string const extractItemName(lsst::mwi::data::DataProperty::PtrType const & properties) {
    if (!properties) {
        LSST_AP_THROW(InvalidParameter, "null DataProperty");
    }
    lsst::mwi::data::DataProperty::PtrType dp = extractRequired(properties, "itemName");
    return boost::any_cast<std::string>(dp->getValue());
}


static std::string const extractPolicyString(
    lsst::mwi::policy::Policy::Ptr const & policy,
    std::string const & key,
    std::string const & def
) {
    if (policy) {
        return policy->getString(key, def);
    } else {
        return def;
    }
}


/*!
    Returns a visit specific table name, given a DataProperty that contains an
    integer-valued property named "visitId" and a string-valued property named "itemName".
 */
LSST_AP_API std::string const getTableName(
    lsst::mwi::policy::Policy::Ptr         const & policy,
    lsst::mwi::data::DataProperty::PtrType const & properties
) {
    std::string   itemName(extractItemName(properties));
    int64_t       visitId = extractVisitId(properties);
    boost::format fmt;
    fmt.exceptions(boost::io::all_error_bits);
    fmt.parse(extractPolicyString(
        policy,
        itemName + ".perVisitTableNamePattern",
        itemName + sDefaultVisitNamePat
    ));
    fmt % visitId;
    return fmt.str();
}


/*!
    Returns the name of the template table that should be used
    to create tables for a particular item.
 */
LSST_AP_API std::string const getTableTemplateName(
    lsst::mwi::policy::Policy::Ptr         const & policy,
    lsst::mwi::data::DataProperty::PtrType const & properties
) {
    std::string itemName(extractItemName(properties));
    return extractPolicyString(policy, itemName + ".templateTableName", itemName);
}


}} // end of namespace lsst::ap

