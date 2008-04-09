// -*- lsst-c++ -*-

/**
 * @file
 * @brief   Implementation of miscellaneous helper functions
 *
 * @ingroup associate
 */

#include <cerrno>
#include <sys/stat.h>
#include <unistd.h>

#include <lsst/ap/Exceptions.h>
#include <lsst/ap/Utils.h>


namespace lsst {
namespace ap {


LSST_AP_API lsst::daf::base::DataProperty::PtrType extractRequired(
    lsst::daf::base::DataProperty::PtrType const & properties,
    std::string                            const & key
) {
    lsst::daf::base::DataProperty::PtrType dp = properties->findUnique(key);
    if (!dp) {
        LSST_AP_THROW(NotFound, boost::format("\"%1%\" property not found") % key);
    }
    return dp;
}


static char const * const sDefaultVisitNamePat = "_visit%1%";


static int64_t extractVisitId(lsst::daf::base::DataProperty::PtrType const & properties) {
    if (!properties) {
        LSST_AP_THROW(InvalidParameter, "null DataProperty");
    }
    lsst::daf::base::DataProperty::PtrType dp = extractRequired(properties, "visitId");
    int64_t visitId = anyToInteger<int64_t>(dp->getValue());
    if (visitId < 0) {
        LSST_AP_THROW(Runtime, "\"visitId\" property value is negative");
    }
    return visitId;
}


static std::string const extractItemName(lsst::daf::base::DataProperty::PtrType const & properties) {
    if (!properties) {
        LSST_AP_THROW(InvalidParameter, "null DataProperty");
    }
    lsst::daf::base::DataProperty::PtrType dp = extractRequired(properties, "itemName");
    return boost::any_cast<std::string>(dp->getValue());
}


static std::string const extractPolicyString(
    lsst::pex::policy::Policy::Ptr const & policy,
    std::string const & key,
    std::string const & def
) {
    if (policy) {
        return policy->getString(key, def);
    } else {
        return def;
    }
}


/**
 * Returns a visit specific table name, given a DataProperty that contains an
 * integer-valued property named "visitId" and a string-valued property named "itemName".
 */
LSST_AP_API std::string const getTableName(
    lsst::pex::policy::Policy::Ptr         const & policy,
    lsst::daf::base::DataProperty::PtrType const & properties
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


/**
 * Returns the name of the template table that should be used
 * to create tables for a particular item.
 */
LSST_AP_API std::string const getTableTemplateName(
    lsst::pex::policy::Policy::Ptr         const & policy,
    lsst::daf::base::DataProperty::PtrType const & properties
) {
    std::string itemName(extractItemName(properties));
    return extractPolicyString(policy, itemName + ".templateTableName", itemName);
}


/**
 * Ensure that all directories along a path exist, creating them if necessary.
 *
 * @param[in] name  Pathname to file to be created
 */
LSST_AP_API void verifyPathName(std::string const & name) {
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
            LSST_AP_THROW(IoError, "Error creating directory: " + dirName);
        }
    }
    else if (ret == -1) {
        // We couldn't read the (existing) directory for some reason.
        LSST_AP_THROW(IoError, "Error searching for directory: " + dirName);
    }
    else if (!S_ISDIR(buf.st_mode)) {
        // It's not a directory.
        LSST_AP_THROW(IoError, "Non-directory in path: " + dirName);
    }
}


}} // end of namespace lsst::ap

