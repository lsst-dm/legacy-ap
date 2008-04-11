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

#include <lsst/daf/base/DataProperty.h>
#include <lsst/pex/policy/Policy.h>
#include <lsst/daf/persistence/Formatter.h>
#include <lsst/daf/persistence/DbStorage.h>

#include "../Results.h"


namespace lsst {
namespace ap {
namespace io {

using lsst::daf::persistence::Formatter;
using lsst::daf::persistence::FormatterRegistration;
using lsst::daf::base::Persistable;
using lsst::daf::persistence::Storage;
using lsst::daf::persistence::DbStorage;
using lsst::pex::policy::Policy;
using lsst::daf::base::DataProperty;


/**
 * @brief  Formatter for MatchPairVector instances.
 *
 * Supports persistence and retrieval with
 * - lsst::daf::persistence::DbStorage
 * - lsst::daf::persistence::DbTsvStorage
 * - lsst::daf::persistence::BoostStorage
 */
class LSST_AP_API MatchPairVectorFormatter : public Formatter {
public:

    virtual ~MatchPairVectorFormatter();

    virtual void write(Persistable const *, Storage::Ptr, DataProperty::PtrType);
    virtual Persistable * read(Storage::Ptr, DataProperty::PtrType);
    virtual void update(Persistable *, Storage::Ptr, DataProperty::PtrType);

    template <class Archive> static void delegateSerialize(Archive &, unsigned int const, Persistable *);

private:

    Policy::Ptr _policy;

    MatchPairVectorFormatter(Policy::Ptr const &);

    static Formatter::Ptr createInstance(Policy::Ptr);
    static FormatterRegistration registration;
};


/**
 * @brief  Formatter for IdPairVector instances.
 *
 * Supports persistence and retrieval with
 * - lsst::daf::persistence::DbStorage
 * - lsst::daf::persistence::DbTsvStorage
 * - lsst::daf::persistence::BoostStorage
 */
class LSST_AP_API IdPairVectorFormatter : public Formatter {
public:

    virtual ~IdPairVectorFormatter();

    virtual void write(Persistable const *, Storage::Ptr, DataProperty::PtrType);
    virtual Persistable * read(Storage::Ptr, DataProperty::PtrType);
    virtual void update(Persistable *, Storage::Ptr, DataProperty::PtrType);

    template <class Archive> static void delegateSerialize(Archive &, unsigned int const, Persistable *);

private:

    Policy::Ptr _policy;

    IdPairVectorFormatter(Policy::Ptr const &);

    static Formatter::Ptr createInstance(Policy::Ptr);
    static FormatterRegistration registration;
};


/**
 * @brief  Formatter for IdVector instances.
 *
 * Supports persistence and retrieval with
 * - lsst::daf::persistence::DbStorage
 * - lsst::daf::persistence::DbTsvStorage
 * - lsst::daf::persistence::BoostStorage
 */
class LSST_AP_API IdVectorFormatter : public Formatter {
public:

    virtual ~IdVectorFormatter();

    virtual void write(Persistable const *, Storage::Ptr, DataProperty::PtrType);
    virtual Persistable * read(Storage::Ptr, DataProperty::PtrType);
    virtual void update(Persistable *, Storage::Ptr, DataProperty::PtrType);

    template <class Archive> static void delegateSerialize(Archive &, unsigned int const, Persistable *);

private:

    Policy::Ptr _policy;

    IdVectorFormatter(Policy::Ptr const &);

    static Formatter::Ptr createInstance(Policy::Ptr);
    static FormatterRegistration registration;
};


}}} // end of namespace lsst::ap::io

#endif // LSST_AP_IO_RESULT_FORMATTERS_H

