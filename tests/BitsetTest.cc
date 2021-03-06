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
 * @brief   Tests for Bitset template and implementation details.
 *
 * @ingroup associate
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE BitsetTest
#include "boost/test/unit_test.hpp"

#include "lsst/ap/Common.h"
#include "lsst/ap/Bitset.h"


using namespace lsst::ap;


template <typename WordT, int NumBits>
void testBitset() {

    Bitset<WordT, NumBits> bs;
    bs.reset();
    for (int i = 0; i < NumBits; ++i) {
        BOOST_CHECK_EQUAL(bs.test(i), false);
        bs.set(i);
        BOOST_CHECK_EQUAL(bs.test(i), true);
        bs.reset(i);
        BOOST_CHECK_EQUAL(bs.test(i), false);
    }

    for (int i = 0; i < NumBits; ++i) { bs.set(i, (i & 3) == 0); }
    for (int i = 0; i < NumBits; ++i) { BOOST_CHECK_EQUAL(bs.test(i), (i & 3) == 0); }

    bs.set();
    for (int i = 0; i < NumBits; ++i) { BOOST_CHECK_EQUAL(bs.test(i), true); }
    bs.reset();
    for (int i = 0; i < NumBits; ++i) { BOOST_CHECK_EQUAL(bs.test(i), false); }
}


template <typename WordT, int NumBits>
void testAllocator() {
    Bitset<WordT, NumBits> bs;
    int indexes[NumBits];
    int nset = 0;

    bs.reset();
    for (int i = 0; i < NumBits; ++i) {
        bs.set(i, (i & 1) != 0);
        nset += (i & 1) != 0;
    }

    BOOST_CHECK_EQUAL(bs.set(indexes, NumBits - nset + 1), false);
    BOOST_CHECK_EQUAL(bs.set(indexes, NumBits - nset), true);
    for (int i = 0; i < NumBits; ++i) { BOOST_CHECK_EQUAL(bs.test(i), true); }

    bs.reset();
    nset = 0;
    for (int i = 0; i < NumBits; ++i) {
        bs.set(i, (i & 3) != 0);
        nset += (i & 3) != 0;
    }
    BOOST_CHECK_EQUAL(bs.set(indexes, NumBits - nset + 1), false);
    if (0 < ((NumBits - nset)>>1)) {
        BOOST_CHECK_EQUAL(bs.set(indexes, (NumBits - nset)>>1), true);
        bs.reset(indexes, (NumBits - nset)>>1);
        BOOST_CHECK_EQUAL(bs.set(indexes, (NumBits - nset)>>1), true);
        BOOST_CHECK_EQUAL(bs.set(indexes, (NumBits - nset) - ((NumBits - nset)>>1)), true);
        for (int i = 0; i < NumBits; ++i) { BOOST_CHECK_EQUAL(bs.test(i), true); }
    }
    bs.reset();
    nset = 0;
    for (int i = 0; i < NumBits; ++i) {
        bs.set(i, i%2 == 0 || i%3 == 0 || i%7 == 0 || i%13 == 0);
        nset += i%2 == 0 || i%3 == 0 || i%7 == 0 || i%13 == 0;
    }
    BOOST_CHECK_EQUAL(bs.set(indexes, NumBits - nset + 1), false);
    if (0 < ((NumBits - nset)>>2)) {
        BOOST_CHECK_EQUAL(bs.set(indexes, (NumBits - nset)>>2), true);
        bs.reset(indexes, (NumBits - nset)>>2);
        BOOST_CHECK_EQUAL(bs.set(indexes, (NumBits - nset)>>2), true);
        BOOST_CHECK_EQUAL(bs.set(indexes, (NumBits - nset)>>2), true);
        BOOST_CHECK_EQUAL(bs.set(indexes, (NumBits - nset)>>2), true);
        bs.reset(indexes, (NumBits - nset)>>2);
        BOOST_CHECK_EQUAL(bs.set(indexes, (NumBits - nset)>>2), true);
        BOOST_CHECK_EQUAL(bs.set(indexes, (NumBits - nset) - 3*((NumBits - nset)>>2)), true);
        BOOST_CHECK_EQUAL(bs.set(indexes, 1), false);
        for (int i = 0; i < NumBits; ++i) { BOOST_CHECK_EQUAL(bs.test(i), true); }
    }
}


BOOST_AUTO_TEST_CASE(bitsetTest) {

    BOOST_TEST_MESSAGE("    - testing Bitset template");

    testBitset<uint8_t,  8>();
    testBitset<uint8_t, 16>();
    testBitset<uint8_t, 40>();
    testBitset<uint8_t,  7>();
    testBitset<uint8_t, 15>();
    testBitset<uint8_t, 41>();
    testBitset<uint8_t, 20>();

    testBitset<uint16_t, 16>();
    testBitset<uint16_t, 32>();
    testBitset<uint16_t, 96>();
    testBitset<uint16_t, 15>();
    testBitset<uint16_t, 33>();
    testBitset<uint16_t, 95>();
    testBitset<uint16_t, 40>();

    testBitset<uint32_t, 32>();
    testBitset<uint32_t, 64>();
    testBitset<uint32_t, 96>();
    testBitset<uint32_t, 31>();
    testBitset<uint32_t, 65>();
    testBitset<uint32_t, 95>();
    testBitset<uint32_t, 37>();

    testBitset<uint64_t,  64>();
    testBitset<uint64_t, 128>();
    testBitset<uint64_t, 256>();
    testBitset<uint64_t,  63>();
    testBitset<uint64_t, 129>();
    testBitset<uint64_t, 255>();
    testBitset<uint64_t,   3>();
}


BOOST_AUTO_TEST_CASE(allocatorTest) {

    BOOST_TEST_MESSAGE("    - testing Bitset based memory allocation");

    testAllocator<uint8_t,  8>();
    testAllocator<uint8_t, 16>();
    testAllocator<uint8_t, 80>();
    testAllocator<uint8_t,  7>();
    testAllocator<uint8_t, 15>();
    testAllocator<uint8_t, 81>();
    testAllocator<uint8_t, 39>();

    testAllocator<uint16_t, 16>();
    testAllocator<uint16_t, 32>();
    testAllocator<uint16_t, 96>();
    testAllocator<uint16_t, 15>();
    testAllocator<uint16_t, 33>();
    testAllocator<uint16_t, 95>();
    testAllocator<uint16_t, 411>();

    testAllocator<uint32_t, 32>();
    testAllocator<uint32_t, 64>();
    testAllocator<uint32_t, 128>();
    testAllocator<uint32_t, 31>();
    testAllocator<uint32_t, 65>();
    testAllocator<uint32_t, 127>();
    testAllocator<uint32_t, 1024>();

    testAllocator<uint64_t,   64>();
    testAllocator<uint64_t,  128>();
    testAllocator<uint64_t,  256>();
    testAllocator<uint64_t,   63>();
    testAllocator<uint64_t,  129>();
    testAllocator<uint64_t,  255>();
    testAllocator<uint64_t,    3>();
    testAllocator<uint64_t, 5311>();
}

