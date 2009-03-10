// -*- lsst-c++ -*-

/**
 * @file
 * @brief   Formatter subclasses for association pipeline result vectors.
 *
 * @ingroup associate
 */

#ifndef LSST_AP_IO_RESULT_FORMATTERS_H
#define LSST_AP_IO_RESULT_FORMATTERS_H

#include <string>
#include <vector>

#include "lsst/daf/base/Persistable.h"
#include "lsst/daf/base/PropertySet.h"
#include "lsst/pex/policy/Policy.h"
#include "lsst/daf/persistence/Formatter.h"
#include "lsst/daf/persistence/DbStorage.h"

#include "../Results.h"


namespace lsst { namespace ap { namespace io {

/**
 * @brief  Formatter for MatchPairVector instances.
 *
 * Supports persistence and retrieval with
 * - lsst::daf::persistence::DbStorage
 * - lsst::daf::persistence::DbTsvStorage
 * - lsst::daf::persistence::BoostStorage
 */
class LSST_AP_API MatchPairVectorFormatter : public lsst::daf::persistence::Formatter {
public:

    virtual ~MatchPairVectorFormatter();

    virtual void write(
        lsst::daf::base::Persistable const *,
        lsst::daf::persistence::Storage::Ptr,
        lsst::daf::base::PropertySet::Ptr
    );
    virtual lsst::daf::base::Persistable * read(
        lsst::daf::persistence::Storage::Ptr, 
        lsst::daf::base::PropertySet::Ptr
    );
    virtual void update(
        lsst::daf::base::Persistable *,
        lsst::daf::persistence::Storage::Ptr,
        lsst::daf::base::PropertySet::Ptr
    );

    template <class Archive>
    static void delegateSerialize(
        Archive &,
        unsigned int const,
        lsst::daf::base::Persistable *
    );

private:

    lsst::pex::policy::Policy::Ptr _policy;

    explicit MatchPairVectorFormatter(lsst::pex::policy::Policy::Ptr);

    static lsst::daf::persistence::Formatter::Ptr createInstance(lsst::pex::policy::Policy::Ptr);
    static lsst::daf::persistence::FormatterRegistration registration;
};


/**
 * @brief  Formatter for IdPairVector instances.
 *
 * Supports persistence and retrieval with
 * - lsst::daf::persistence::DbStorage
 * - lsst::daf::persistence::DbTsvStorage
 * - lsst::daf::persistence::BoostStorage
 */
class LSST_AP_API IdPairVectorFormatter : public lsst::daf::persistence::Formatter {
public:

    virtual ~IdPairVectorFormatter();

    virtual void write(
        lsst::daf::base::Persistable const *,
        lsst::daf::persistence::Storage::Ptr,
        lsst::daf::base::PropertySet::Ptr
    );
    virtual lsst::daf::base::Persistable * read(
        lsst::daf::persistence::Storage::Ptr, 
        lsst::daf::base::PropertySet::Ptr
    );
    virtual void update(
        lsst::daf::base::Persistable *,
        lsst::daf::persistence::Storage::Ptr,
        lsst::daf::base::PropertySet::Ptr
    );

    template <class Archive>
    static void delegateSerialize(
        Archive &,
        unsigned int const,
        lsst::daf::base::Persistable *
    );

private:

    lsst::pex::policy::Policy::Ptr _policy;

    explicit IdPairVectorFormatter(lsst::pex::policy::Policy::Ptr);

    static lsst::daf::persistence::Formatter::Ptr createInstance(lsst::pex::policy::Policy::Ptr);
    static lsst::daf::persistence::FormatterRegistration registration;
};


/**
 * @brief  Formatter for IdVector instances.
 *
 * Supports persistence and retrieval with
 * - lsst::daf::persistence::DbStorage
 * - lsst::daf::persistence::DbTsvStorage
 * - lsst::daf::persistence::BoostStorage
 */
class LSST_AP_API IdVectorFormatter : public lsst::daf::persistence::Formatter {
public:

    virtual ~IdVectorFormatter();

    virtual void write(
        lsst::daf::base::Persistable const *,
        lsst::daf::persistence::Storage::Ptr,
        lsst::daf::base::PropertySet::Ptr
    );
    virtual lsst::daf::base::Persistable * read(
        lsst::daf::persistence::Storage::Ptr, 
        lsst::daf::base::PropertySet::Ptr
    );
    virtual void update(
        lsst::daf::base::Persistable *,
        lsst::daf::persistence::Storage::Ptr,
        lsst::daf::base::PropertySet::Ptr
    );

    template <class Archive>
    static void delegateSerialize(
        Archive &,
        unsigned int const,
        lsst::daf::base::Persistable *
    );

private:

    lsst::pex::policy::Policy::Ptr _policy;

    explicit IdVectorFormatter(lsst::pex::policy::Policy::Ptr);

    static lsst::daf::persistence::Formatter::Ptr createInstance(lsst::pex::policy::Policy::Ptr);
    static lsst::daf::persistence::FormatterRegistration registration;
};


}}} // end of namespace lsst::ap::io

#endif // LSST_AP_IO_RESULT_FORMATTERS_H

