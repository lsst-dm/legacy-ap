// -*- lsst-c++ -*-
//
//##====----------------                                ----------------====##/
//
//! \file   ResultFormatters.h
//! \brief  Formatter subclasses for association pipeline results.
//
//##====----------------                                ----------------====##/

#ifndef LSST_AP_IO_RESULT_FORMATTERS_H
#define LSST_AP_IO_RESULT_FORMATTERS_H

#include <string>
#include <vector>

#include <lsst/mwi/data/DataProperty.h>
#include <lsst/mwi/policy/Policy.h>
#include <lsst/mwi/persistence/Formatter.h>
#include <lsst/mwi/persistence/DbStorage.h>

#include "../Results.h"


namespace lsst {
namespace ap {
namespace io {

using lsst::mwi::persistence::Formatter;
using lsst::mwi::persistence::FormatterRegistration;
using lsst::mwi::persistence::Persistable;
using lsst::mwi::persistence::Storage;
using lsst::mwi::persistence::DbStorage;
using lsst::mwi::policy::Policy;
using lsst::mwi::data::DataProperty;


/*!
    \brief  Formatter for MatchPairVector instances.

    Supports persistence and retrieval with
    - lsst::mwi::persistence::DbStorage
    - lsst::mwi::persistence::DbTsvStorage
    - lsst::mwi::persistence::BoostStorage
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


/*!
    \brief  Formatter for IdPairVector instances.

    Supports persistence and retrieval with
    - lsst::mwi::persistence::DbStorage
    - lsst::mwi::persistence::DbTsvStorage
    - lsst::mwi::persistence::BoostStorage
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


/*!
    \brief  Formatter for IdVector instances.

    Supports persistence and retrieval with
    - lsst::mwi::persistence::DbStorage
    - lsst::mwi::persistence::DbTsvStorage
    - lsst::mwi::persistence::BoostStorage
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

