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
 * @brief   Utility class for automatically invoking a function when leaving a scope.
 *
 * @ingroup ap
 */

#ifndef LSST_AP_SCOPE_GUARD_H
#define LSST_AP_SCOPE_GUARD_H

#include "boost/noncopyable.hpp"
#include "boost/function.hpp"
#include "boost/scoped_ptr.hpp"

#include "Common.h"


namespace lsst { namespace ap {

/**
 * @brief   Utility class for automatically invoking a function when leaving a scope.
 *
 * See http://www.ddj.com/cpp/184403758 and http://www.zete.org/people/jlehrer/scopeguard.html
 * for details on the origins and uses of this class.
 */
class LSST_AP_API ScopeGuard : public boost::function<void ()>, private boost::noncopyable {

public :

    typedef boost::function<void ()> Base;

    template <typename F>
    explicit ScopeGuard(F const & f) : Base(f), _dismissed(false) {
        ++_numGuards; // stop compiler from optimizing away entire object instances
    }

    ~ScopeGuard() {
        --_numGuards;
        try {
            if (!_dismissed) {
                operator()();
            }
        } catch (...) {}
    }

    bool dismissed() const throw() {
        return _dismissed;
    }
    void dismiss() throw() {
        _dismissed = true;
    }

private :

    bool _dismissed;

    static int volatile _numGuards;
};

}} // end of namespace lsst::ap

#endif  // LSST_AP_SCOPE_GUARD_H

