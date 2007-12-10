// -*- lsst-c++ -*-
//
//##====----------------                                ----------------====##/
//
//! \file   ScopeGuard.h
//! \brief  Utility class for automatically invoking a function when leaving a scope.
//
//##====----------------                                ----------------====##/

#ifndef LSST_AP_SCOPE_GUARD_H
#define LSST_AP_SCOPE_GUARD_H

#include <boost/noncopyable.hpp>
#include <boost/function.hpp>

#include <lsst/ap/Common.h>


namespace lsst {
namespace ap {


/*!
    \brief  Utility class for automatically invoking a function when leaving a scope.

    See http://www.ddj.com/cpp/184403758 and http://www.zete.org/people/jlehrer/scopeguard.html
    for details on the origins and uses of this class.
 */
class LSST_AP_LOCAL ScopeGuard : public boost::function<void ()>, private boost::noncopyable {

public :

    typedef boost::function<void ()> BaseType;

    template <typename F>
    explicit ScopeGuard(F const & f) : BaseType(f) {}

    ~ScopeGuard() {
        try {
            if (!_dismissed) {
                operator()();
            }
        } catch (...) {}
    }

    bool dismissed() const { return _dismissed; }
    void dismiss()         { _dismissed = true; }

private :

    bool _dismissed;
};


}} // end of namespace lsst::ap

#endif  // LSST_AP_SCOPE_GUARD_H
