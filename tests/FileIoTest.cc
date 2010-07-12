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
 * @brief   Tests for low-level file I/O classes.
 *
 * @ingroup associate
 */

#include <unistd.h>

#include <iostream>

#include "boost/bind.hpp"
#include "boost/scoped_array.hpp"
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE FileIoTest
#include "boost/test/unit_test.hpp"

#include "lsst/afw/math/Random.h"

#include "lsst/ap/Common.h"
#include "lsst/ap/ScopeGuard.h"
#include "lsst/ap/Time.h"
#include "lsst/ap/io/FileIo.h"


using std::size_t;
using boost::uint8_t;
using boost::uint32_t;
using lsst::afw::math::Random;

using namespace lsst::ap;


namespace {

Random & rng() {
    static Random * generator = 0;
    if (generator == 0) {
        TimeSpec ts;
        ts.systemTime();
        generator = new Random(Random::MT19937, static_cast<unsigned long>(ts.tv_sec + ts.tv_nsec));
        std::clog << "\n"
            << "     /\n"
            << "    | Note: Using random number seed " << generator->getSeed() << "\n"
            << "    |       and algorithm " << generator->getAlgorithmName() << "\n"
            << "     \\\n" << std::endl;
    }
    return *generator;
}


void doWrite(
    io::SequentialWriter & writer,
    uint8_t const * const  content,
    size_t  const          contentLength
) {
    writer.write(content, contentLength);
    writer.finish();
}


void doBlockedWrite(
    io::SequentialWriter & writer,
    uint8_t const * const  content,
    size_t  const          contentLength,
    size_t  const          blockSize
) {
    for (size_t i = 0; i*blockSize < contentLength; ++i) {
        size_t nb = contentLength - i*blockSize;
        writer.write(content + i*blockSize, nb < blockSize ? nb : blockSize);
    }
    writer.finish();
}


void doRead(
    io::SequentialReader & reader,
    uint8_t  * const       content,
    size_t const           contentLength
) {
    uint8_t * buf    = content;
    size_t    buflen = contentLength;
    size_t    nb;
    while ((nb = reader.read(buf, buflen)) > 0) {
        BOOST_REQUIRE(nb <= buflen);
        buf    += nb;
        buflen -= nb;
    }
    BOOST_CHECK_EQUAL(buflen, 0u);
    uint8_t slop[1];
    nb = reader.read(slop, 1);
    BOOST_CHECK_EQUAL(nb, 0u);
    BOOST_CHECK(reader.finished());
}


void doBlockedRead(
    io::SequentialReader & reader,
    uint8_t  * const       content,
    size_t const           contentLength,
    size_t const           blockSize
) {
    for (size_t i = 0; i*blockSize < contentLength; ++i) {
        uint8_t * buf    = content + i*blockSize;
        size_t    buflen = contentLength - i*blockSize;
        if (buflen > blockSize) {
            buflen = blockSize;
        }
        while (buflen > 0) {
            size_t nb = reader.read(buf, buflen);
            BOOST_REQUIRE(nb <= buflen);
            BOOST_REQUIRE(nb > 0);
            buf    += nb;
            buflen -= nb;
        }
    }
    uint8_t slop[1];
    size_t nb = reader.read(slop, 1);
    BOOST_CHECK_EQUAL(nb, 0u);
    BOOST_CHECK(reader.finished());
}


std::string const makeTempFile() {
    char name[64];
    std::strncpy(name, "/tmp/FileIoTest.XXXXXX", 63);
    name[63] = 0;
    int const fd = ::mkstemp(name);
    if (fd < 1) {
        BOOST_FAIL("Failed to create temporary file for testing purposes");
    }
    ::close(fd);
    return std::string(name);
}

} // end of anonymous namespace


BOOST_AUTO_TEST_CASE(sequentialIoTest1) {

    BOOST_TEST_MESSAGE("    - roundtrip file IO test, one-shot, no compression");
    std::string const name(makeTempFile());
    ScopeGuard fileGuard(boost::bind(::unlink, name.c_str()));

    size_t len = static_cast<size_t>(rng().flat(8192, 65536));
    boost::scoped_array<uint32_t> content(new uint32_t[len]);
    boost::scoped_array<uint32_t> contentCheck(new uint32_t[len]);
    for (uint32_t i = 0; i < len; ++i) { content[i] = i; }

    // write out, read back in, then compare contents
    io::SequentialFileWriter w(name, true);
    doWrite(w, reinterpret_cast<uint8_t *>(content.get()), len*sizeof(uint32_t));
    io::SequentialFileReader r(name);
    doRead(r, reinterpret_cast<uint8_t *>(contentCheck.get()), len*sizeof(uint32_t));
    int cmp = std::memcmp(contentCheck.get(), content.get(), len*sizeof(uint32_t));
    BOOST_CHECK_EQUAL(cmp, 0);
}


BOOST_AUTO_TEST_CASE(sequentialIoTest2) {

    BOOST_TEST_MESSAGE("    - roundtrip file IO test, blocked, no compression");
    std::string const name(makeTempFile());
    ScopeGuard  fileGuard(boost::bind(::unlink, name.c_str()));

    size_t len = static_cast<size_t>(rng().flat(3000, 65536));
    boost::scoped_array<uint32_t> content(new uint32_t[len]);
    boost::scoped_array<uint32_t> contentCheck(new uint32_t[len]);
    for (uint32_t i = 0; i < len; ++i) { content[i] = i; }

    // write out, read back in, then compare contents (no direct I/O)
    io::SequentialFileWriter w(name, true);
    doBlockedWrite(w, reinterpret_cast<uint8_t *>(content.get()), len*sizeof(uint32_t), 1024);
    io::SequentialFileReader r(name);
    doBlockedRead(r, reinterpret_cast<uint8_t *>(contentCheck.get()), len*sizeof(uint32_t), 2001);
    int cmp = std::memcmp(contentCheck.get(), content.get(), len*sizeof(uint32_t));
    BOOST_CHECK_EQUAL(cmp, 0);
}


BOOST_AUTO_TEST_CASE(sequentialIoTest3) {

    BOOST_TEST_MESSAGE("    - roundtrip file IO test, one-shot, with compression");
    std::string const name(makeTempFile());
    ScopeGuard  fileGuard(boost::bind(::unlink, name.c_str()));

    size_t len = static_cast<size_t>(rng().flat(512256, 2049024));
    boost::scoped_array<uint32_t> content(new uint32_t[len]);
    boost::scoped_array<uint32_t> contentCheck(new uint32_t[len]);
    for (uint32_t i = 0; i < len; ++i) { content[i] = i; }

    // write out, read back in, then compare contents
    io::CompressedFileWriter w(name, true);
    doBlockedWrite(w, reinterpret_cast<uint8_t *>(content.get()), len*sizeof(uint32_t), 1500);
    io::CompressedFileReader r(name);
    doBlockedRead(r, reinterpret_cast<uint8_t *>(contentCheck.get()), len*sizeof(uint32_t), 3100);
    int cmp = std::memcmp(contentCheck.get(), content.get(), len*sizeof(uint32_t));
    BOOST_CHECK_EQUAL(cmp, 0);
}


BOOST_AUTO_TEST_CASE(emptyFileTest) {

    BOOST_TEST_MESSAGE("    - empty file IO test");
    {
        std::string const name(makeTempFile());
        ScopeGuard  fileGuard(boost::bind(::unlink, name.c_str()));
        io::CompressedFileWriter w(name, true);
        doWrite(w, 0, 0);
        io::CompressedFileReader r(name);
        doRead(r, 0, 0);
    }
    {
        std::string const name(makeTempFile());
        ScopeGuard  fileGuard(boost::bind(::unlink, name.c_str()));
        io::CompressedFileWriter w(name, true);
        doWrite(w, 0, 0);
        io::CompressedFileReader r(name);
        doRead(r, 0, 0);
    }
}
