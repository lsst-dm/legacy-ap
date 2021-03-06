// -*- lsst-c++ -*-

/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
 * 
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the LSST License Statement and 
 * the GNU General Public License along with this program.  If not, 
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */
 

/**
 * @file
 * @brief   Formatter subclasses for association pipeline result vectors.
 *
 * @ingroup ap
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
class MatchPairVectorFormatter : public lsst::daf::persistence::Formatter {
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
class IdPairVectorFormatter : public lsst::daf::persistence::Formatter {
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
class IdVectorFormatter : public lsst::daf::persistence::Formatter {
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

