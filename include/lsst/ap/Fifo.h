// -*- lsst-c++ -*-

/**
 * @file
 * @brief   A fixed capacity FIFO buffer for integers.
 *
 * @ingroup associate
 */

#ifndef LSST_AP_FIFO_H
#define LSST_AP_FIFO_H

#include "boost/noncopyable.hpp"
#include "boost/static_assert.hpp"

#include "lsst/pex/exceptions.h"

#include "Common.h"


namespace lsst { namespace ap {

/** @brief  A First In, First Out (FIFO) queue of fixed capacity. */
template <int NumEntries>
class Fifo : private boost::noncopyable {

    // NumEntries must be a positive power of 2
    BOOST_STATIC_ASSERT(NumEntries > 0 && (NumEntries & (NumEntries - 1)) == 0);

public :

    /// Creates an empty Fifo.
    Fifo() { clear(); }

    /// Empties the Fifo.
    void clear() {
        _size  = 0;
        _back  = 0;
        _front = 0;
    }

    /// Returns @c true if the Fifo is empty
    bool empty() const {
        return _size == 0;
    }

    /// Returns @c true if the Fifo is full
    bool full() const {
        return _size == NumEntries;
    }

    /**
     * Inserts the given integer into the Fifo.
     *
     * @throw lsst::pex::exceptions::LengthError    Thrown if the Fifo is full.
     */
    void enqueue(boost::int64_t const elt) {
        int sz = _size;
        if (sz == NumEntries) {
            throw LSST_EXCEPT(lsst::pex::exceptions::LengthErrorException,
                              "unable to insert element into full FIFO");
        }
        int i = _back;
        _buffer[i] = elt;
        _back = (i + 1) & (NumEntries - 1);
        _size = sz + 1;
    }

    /**
     * Removes the least recently inserted integer from the Fifo.
     *
     * @throw lsst::pex::exceptions::LengthError   Thrown if the Fifo is empty.
     */
    boost::int64_t dequeue() {
        int sz = _size;
        if (sz == 0) {
            throw LSST_EXCEPT(lsst::pex::exceptions::LengthErrorException,
                              "unable to remove element from empty FIFO");
        }
        int i = _front;
        boost::int64_t elt = _buffer[i];
        _front = (i + 1) & (NumEntries - 1);
        _size = sz - 1;
        return elt;
    }

private :

    boost::int64_t _buffer[NumEntries];
    int _size;
    int _back;
    int _front;
};


}} // end of namespace lsst::ap

#endif // LSST_AP_FIFO_H
