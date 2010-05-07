// -*- lsst-c++ -*-
/** @file
  * @brief Formatters for persistable classes in lsst::ap::cluster. 
  *
  * @ingroup ap
  * @author Serge Monkewitz
  */
#ifndef LSST_AP_CLUSTER_FORMATTERS_H
#define LSST_AP_CLUSTER_FORMATTERS_H

#include "lsst/daf/base/Persistable.h"
#include "lsst/daf/base/PropertySet.h"
#include "lsst/pex/policy/Policy.h"
#include "lsst/daf/persistence/Formatter.h"
#include "lsst/daf/persistence/Storage.h"

#include "../Common.h"


namespace lsst { namespace ap { namespace cluster {

/**
 * @brief  Formatter for PersistableSourceClusterVector.
 *
 * Supports persistence and retrieval with
 * @li lsst::daf::persistence::DbStorage
 * @li lsst::daf::persistence::DbTsvStorage
 * @li lsst::daf::persistence::BoostStorage
 */
class LSST_AP_API SourceClusterVectorFormatter :
    public lsst::daf::persistence::Formatter
{
public:
    virtual ~SourceClusterVectorFormatter();

    virtual void write(
        lsst::daf::base::Persistable const * persistable,
        lsst::daf::persistence::Storage::Ptr storage,
        lsst::daf::base::PropertySet::Ptr additionalData);

    virtual lsst::daf::base::Persistable * read(
        lsst::daf::persistence::Storage::Ptr storage,
        lsst::daf::base::PropertySet::Ptr additionalData);

    virtual void update(
        lsst::daf::base::Persistable * persistable,
        lsst::daf::persistence::Storage::Ptr storage,
        lsst::daf::base::PropertySet::Ptr additionalData);

    template <class Archive>
    static void delegateSerialize(
        Archive & ar,
        unsigned int const version,
        lsst::daf::base::Persistable * persistable);

private:
    lsst::pex::policy::Policy::Ptr _policy;

    static lsst::daf::persistence::FormatterRegistration registration;

    SourceClusterVectorFormatter(lsst::pex::policy::Policy::Ptr policy);

    static lsst::daf::persistence::Formatter::Ptr createInstance(
        lsst::pex::policy::Policy::Ptr policy);
};

}}} // namespace lsst::ap::cluster

#endif // LSST_AP_CLUSTER_FORMATTERS_H
