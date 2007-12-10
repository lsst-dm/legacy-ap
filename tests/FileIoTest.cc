// -*- lsst-c++ -*-
//
//##====----------------                                ----------------====##/
//
//! \file   FileIoTest.cc
//! \brief  Testing of file I/O classes.
//
//##====----------------                                ----------------====##/

#include <unistd.h>

#include <boost/bind.hpp>
#include <boost/scoped_array.hpp>
#include <boost/version.hpp>
#if BOOST_VERSION < 103400
#   include <boost/test/auto_unit_test.hpp>
#   define BOOST_TEST_MESSAGE BOOST_MESSAGE
#else
#   include <boost/test/unit_test.hpp>
#endif

#include <lsst/ap/Common.h>
#include <lsst/ap/Random.h>
#include <lsst/ap/ScopeGuard.h>
#include <lsst/ap/io/FileIo.h>


using namespace lsst::ap;


namespace {

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
    ScopeGuard  fileGuard(boost::bind(::unlink, name.c_str()));

    initRandom();
    size_t len = 16384;
    if (coinToss(0.5)) {
        len += static_cast<size_t>(uniformRandom()*16384.0) - 8192;
    }
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

    initRandom();
    size_t len = 15977;
    if (coinToss(0.5)) {
        len += static_cast<size_t>(uniformRandom()*16384.0) - 8192;
    }
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

    initRandom();
    size_t len = 1024512;
    if (coinToss(0.5)) {
        len += static_cast<size_t>(uniformRandom()*1024512.0) - 512256;
    }
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
