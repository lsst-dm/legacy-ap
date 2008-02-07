// -*- lsst-c++ -*-

/**
 * @file
 * @brief   Utility class for automatically invoking a function when leaving a scope.
 *
 * @ingroup associate
 */

#ifndef LSST_AP_SCOPE_GUARD_H
#define LSST_AP_SCOPE_GUARD_H

#include <boost/noncopyable.hpp>
#include <boost/function.hpp>
#include <boost/scoped_ptr.hpp>

#include "Common.h"


namespace lsst {
namespace ap {


/**
 * @brief   Utility class for automatically invoking a function when leaving a scope.
 *
 * See http://www.ddj.com/cpp/184403758 and http://www.zete.org/people/jlehrer/scopeguard.html
 * for details on the origins and uses of this class.
 */
class LSST_AP_API ScopeGuard : public boost::function<void ()>, private boost::noncopyable {

public :

    typedef boost::function<void ()> BaseType;

    template <typename F>
    explicit ScopeGuard(F const & f) : BaseType(f), _dismissed(false) {
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

    bool dismissed() const throw() { return _dismissed; }
    void dismiss()         throw() { _dismissed = true; }

private :

    bool _dismissed;

    static int volatile _numGuards;
};


}} // end of namespace lsst::ap

#endif  // LSST_AP_SCOPE_GUARD_H

