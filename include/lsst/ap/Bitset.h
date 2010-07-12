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
 * @brief   A class for manipulating a fixed set of bits at the individual bit level.
 *
 * @ingroup ap
 */

#ifndef LSST_AP_BITSET_H
#define LSST_AP_BITSET_H

#include <cassert>
#include <cstring>

#include "boost/static_assert.hpp"

#include "Common.h"


namespace lsst { namespace ap { namespace detail {

template <typename WordT> struct BitTraits {
    static bool const IS_SPECIALIZED = false;
};

template <> struct LSST_AP_LOCAL BitTraits<uint8_t>  {
    static bool const IS_SPECIALIZED      = true;
    static int  const BITS_PER_WORD_LOG2  = 3;
    static int  const BITS_PER_WORD       = 8;
    static boost::uint8_t const WORD_MASK = 0xff;
};

template <> struct LSST_AP_LOCAL BitTraits<uint16_t> {
    static bool const IS_SPECIALIZED       = true;
    static int  const BITS_PER_WORD_LOG2   = 4;
    static int  const BITS_PER_WORD        = 16;
    static boost::uint16_t const WORD_MASK = 0xffff;
};

template <> struct LSST_AP_LOCAL BitTraits<uint32_t> {
    static bool const IS_SPECIALIZED       = true;
    static int  const BITS_PER_WORD_LOG2   = 5;
    static int  const BITS_PER_WORD        = 32;
    static boost::uint32_t const WORD_MASK = 0xffffffff;
};

template <> struct LSST_AP_LOCAL BitTraits<uint64_t> {
    static bool const IS_SPECIALIZED       = true;
    static int  const BITS_PER_WORD_LOG2   = 6;
    static int  const BITS_PER_WORD        = 64;
    static boost::uint64_t const WORD_MASK = UINT64_C(0xffffffffffffffff);
};

template <typename WordT>
inline int wordForBit(int const i) {
    return (i >> BitTraits<WordT>::BITS_PER_WORD_LOG2);
}

template <typename WordT>
inline WordT maskForBit(int const i) {
    return static_cast<WordT>(1) << (i & (BitTraits<WordT>::BITS_PER_WORD - 1));
}

template <typename WordT>
bool setBits(
    int * const indexes,
    WordT * const words,
    int const numBitsToSet,
    int const numBits
);

template <typename WordT>
void resetBits(
    WordT * const words,
    int const * const indexes,
    int const numBitsToReset,
    int const numBits
);

} // end of namespace detail


/** @brief  A fixed size set of bits. */
template <typename WordT, int NumBits>
class Bitset {

private :

    BOOST_STATIC_ASSERT(NumBits > 0);
    BOOST_STATIC_ASSERT(detail::BitTraits<WordT>::IS_SPECIALIZED);

public :

    static int const NUM_BITS  = NumBits;
    static int const NUM_WORDS = (NumBits + (detail::BitTraits<WordT>::BITS_PER_WORD - 1)) >>
                                 detail::BitTraits<WordT>::BITS_PER_WORD_LOG2;

    /** Clears all bits. */
    void reset() {
        std::memset(_bits, 0, sizeof(_bits));
    }

    /** Sets all bits to 1. */
    void set() {
        std::memset(_bits, -1, sizeof(_bits));
    }

    /** Sets the i-th bit in the set to zero. */
    void reset(int const i) {
        assert(i >= 0 && i < NumBits);
        _bits[detail::wordForBit<WordT>(i)] &= ~ detail::maskForBit<WordT>(i);
    }

    /** Sets the i-th bit in the set to one. */
    void set(int const i) {
        assert(i >= 0 && i < NumBits);
        _bits[detail::wordForBit<WordT>(i)] |= detail::maskForBit<WordT>(i);
    }

    /** Sets the i-th bit in the set to one if @a on is @c true and to zero otherwise. */
    void set(int const i, bool const on) {
        if (on) {
            set(i);
        } else {
            reset(i);
        }
    }

    /** Returns @c true if the i-th bit in the set is one and @c false otherwise. */
    bool test(int const i) const {
        assert(i >= 0 && i < NumBits);
        return _bits[detail::wordForBit<WordT>(i)] & detail::maskForBit<WordT>(i);
    }

    /**
     * If at least @a numBits zero bits are available in this Bitset, this function sets the first
     * @a numBits of them to one and returns @c true. Otherwise, @c false is returned.
     */
    bool set(int * const indexes, int const numBits) {
        return detail::setBits<WordT>(indexes, _bits, numBits, NumBits);
    }

    /** Sets @a numBits bits identified by the integers in @a indexes to zero. */
    void reset(int const * const indexes, int const numBits) {
        detail::resetBits<WordT>(_bits, indexes, numBits, NumBits);
    }

private :

    WordT _bits[NUM_WORDS];
};


}} // end of namespace lsst::ap

#endif // LSST_AP_BITSET_H
