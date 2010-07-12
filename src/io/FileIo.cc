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
 * @brief   Implementation of low level file IO classes.
 *
 * @ingroup ap
 */

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/fcntl.h>
#include <aio.h>
#include <errno.h>

#include <zlib.h>
#if ZLIB_VERNUM < 0x123
#    warning Older version of zlib detected, upgrading to version 1.2.3 or later is recommended
#endif

#include <cstring>
#include <algorithm>

#include "boost/bind.hpp"
#include "boost/format.hpp"

#include "lsst/pex/exceptions.h"

#include "lsst/ap/Common.h"
#include "lsst/ap/ScopeGuard.h"
#include "lsst/ap/Time.h"
#include "lsst/ap/io/FileIo.h"

namespace ex = lsst::pex::exceptions;

namespace lsst { namespace ap { namespace io { namespace {

int openFile(
    std::string const & fileName,
    int         const   oflag,
    ::mode_t    const   mode = 0
) {
    int fd = ::open(fileName.c_str(), oflag, mode);
    if (fd == -1) {
        if (errno != ENOENT || (oflag & O_WRONLY) != 0) {
            throw LSST_EXCEPT(ex::IoErrorException,
                (boost::format("open() failed on file %1%, flags: %2%, errno: %3%") %
                    fileName % oflag % errno).str());
        }
        // file didn't exist (ok when opening for reading)
    }
    return fd;
}

}}}} // end of anonymous namespace


// -- SequentialIoBase ----------------

lsst::ap::io::SequentialIoBase::SequentialIoBase() : _state(IN_PROGRESS) {}

lsst::ap::io::SequentialIoBase::~SequentialIoBase() {}


// -- SequentialFileReader ----------------

lsst::ap::io::SequentialFileReader::SequentialFileReader(
    std::string const & fileName
) :
    _fd(-1)
{
    int const fd = openFile(fileName, false, O_RDONLY);
    if (fd == -1) {
        _state = FINISHED;
        return;
    }
    _fd = fd;
}


lsst::ap::io::SequentialFileReader::~SequentialFileReader() { cleanup(); }


void lsst::ap::io::SequentialFileReader::cleanup() {
    if (_fd != -1) {
        ::close(_fd);
        _fd = -1;
    }
}


std::size_t lsst::ap::io::SequentialFileReader::read(
    unsigned char * const buf,
    std::size_t const len
) {
    if (len == 0) {
        return 0;
    }
    if (buf == 0) {
        throw LSST_EXCEPT(ex::InvalidParameterException,
                          "null pointer to read destination");
    }
    if (_state == FAILED) {
        throw LSST_EXCEPT(ex::IoErrorException,
                          "read() called on a failed SequentialFileReader");
    } else if (_state == FINISHED) {
        return 0;
    }

    ::ssize_t n = ::read(_fd, buf, len);
    if (n < 0) {
        cleanup(FAILED);
        throw LSST_EXCEPT(ex::IoErrorException,
                          (boost::format("read() failed, errno: %1%") % errno).str());
    } else if (n == 0) {
        cleanup(FINISHED);
    }
    return static_cast<std::size_t>(n);
}


// -- SequentialFileWriter ----------------

lsst::ap::io::SequentialFileWriter::SequentialFileWriter(
    std::string const & fileName,
    bool        const   overwrite
) :
    _fd(-1)
{
    int const fd = openFile(
        fileName,
        O_WRONLY | O_CREAT | O_APPEND | (overwrite ? O_TRUNC : O_EXCL),
        S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH
    );
    assert(fd != -1);
    _fd = fd;
}


lsst::ap::io::SequentialFileWriter::~SequentialFileWriter() { cleanup(); }


void lsst::ap::io::SequentialFileWriter::cleanup() {
    if (_fd != -1) {
        ::close(_fd);
    }
}


void lsst::ap::io::SequentialFileWriter::write(
    unsigned char const * const buf,
    std::size_t const len
) {
    if (len == 0) {
        return;
    }
    if (buf == 0) {
        throw LSST_EXCEPT(ex::InvalidParameterException,
                          "null pointer to bytes to write");
    }
    if (_state != IN_PROGRESS) {
        throw LSST_EXCEPT(ex::IoErrorException,
                          "write() called on a finished or failed SequentialFileWriter");
    }

    unsigned char const * dst = buf;
    std::size_t nb = len;
    do {
        ::ssize_t n = ::write(_fd, dst, nb);
        if (n < 0) {
            cleanup(FAILED);
            throw LSST_EXCEPT(ex::IoErrorException,
                              (boost::format("write() failed, errno: %1%") % errno).str());
        }
        dst += n;
        nb  -= n;
    } while (nb > 0);
}


void lsst::ap::io::SequentialFileWriter::finish() {
    if (_state != IN_PROGRESS) {
        throw LSST_EXCEPT(ex::IoErrorException,
                          "finish() called on a finished or failed SequentialFileWriter");
    }

    // Flush userland and kernel buffers
    while (fsync(_fd) != 0) {
        if (errno != EINTR) {
            cleanup(FAILED);
            throw LSST_EXCEPT(ex::IoErrorException,
                              (boost::format("fsync() failed, errno: %1%") % errno).str());
        }
        errno = 0;
    }
    cleanup(FINISHED);
}


// -- CompressedFileReader ----------------

lsst::ap::io::CompressedFileReader::CompressedFileReader(
    std::string const & fileName,
    std::size_t const   blockSize
) :
    _memory(0),
    _buffers(0),
    _blockSize(blockSize),
    _fileSize(0),
    _remaining(0),
    _fd(-1)
{
    if (blockSize < 1024 || blockSize > 16777216 || (blockSize & (blockSize - 1)) != 0) {
        throw LSST_EXCEPT(ex::InvalidParameterException,
                          "I/O block size must be a power of 2 between 2^10 and 2^24 bytes");
    }

    std::memset(&_stream,  0, sizeof(::z_stream));
    std::memset(&_request, 0, sizeof(::aiocb));
    _request.aio_fildes = -1;

    // allocate IO buffers, store pointer to 8k aligned location inside buffer
    boost::scoped_array<unsigned char> mem(new unsigned char[2*blockSize + 8192]);
    std::size_t buf = reinterpret_cast<std::size_t>(mem.get());
    _buffers = reinterpret_cast<unsigned char *>((buf + 8191) & ~static_cast<std::size_t>(8191));

    // open file
    int const fd = openFile(fileName, O_RDONLY);
    if (fd == -1) {
        _state = FINISHED;
        return;
    }
    ScopeGuard fdGuard(boost::bind(::close, fd));

    // determine file size
    struct stat sb;
    if (::fstat(fd, &sb) != 0) {
        throw LSST_EXCEPT(ex::IoErrorException,
            (boost::format("fstat() failed on file %1%, errno: %2%") % fileName % errno).str());
    }
    _fileSize = static_cast<std::size_t>(sb.st_size);
    if (_fileSize == 0) {
        _state = FINISHED;
        return;
    }
    _remaining = _fileSize;

    // setup zlib (specify auto-detection of zlib/gzip header)
    if (inflateInit2(&_stream, MAX_WBITS + 32) != Z_OK) {
        throw LSST_EXCEPT(ex::RuntimeErrorException,
                          "[zlib] inflateInit2() failed to initialize decompression stream");
    }
    ScopeGuard zlibGuard(boost::bind(::inflateEnd, &_stream));

    // issue first async read
    std::size_t const nb = std::min(blockSize, _fileSize);
    _request.aio_fildes = fd;
    _request.aio_buf    = _buffers;
    _request.aio_nbytes = nb;
    if (::aio_read(&_request) != 0) {
        throw LSST_EXCEPT(ex::IoErrorException,
            (boost::format("aio_read() failed to enqueue IO request for file %1%, errno: %2%") %
                fileName % errno).str());
    }
    _fd = fd;
    using std::swap;
    swap(_memory, mem);
    zlibGuard.dismiss();
    fdGuard.dismiss();
}


lsst::ap::io::CompressedFileReader::~CompressedFileReader() {
    cleanup();
}


void lsst::ap::io::CompressedFileReader::cleanup() {
    if (_fd != -1) {
        if (_state == IN_PROGRESS) {
            ::aio_cancel(_fd, 0);
        }
        ::close(_fd);
        ::inflateEnd(&_stream);
        _fd = -1;
    }
}


std::size_t lsst::ap::io::CompressedFileReader::read(
    unsigned char * const buf,
    std::size_t const len
) {
    if (len == 0) {
        return 0;
    }
    if (buf == 0) {
        throw LSST_EXCEPT(ex::InvalidParameterException, "null pointer to read destination");
    }
    if (_state == FAILED) {
        throw LSST_EXCEPT(ex::IoErrorException, "read() called on a failed CompressedFileReader");
    } else if (_state == FINISHED) {
        return 0;
    }

    // if the current block has not been fully transferred to the user,
    // decompress as much as possible of what remains
    _stream.next_out  = buf;
    _stream.avail_out = len;
    if (_stream.avail_in > 0) {
        int const zret = ::inflate(&_stream, Z_SYNC_FLUSH);
        if (zret == Z_STREAM_END) {
            if (_remaining != 0) {
                cleanup(FAILED);
                throw LSST_EXCEPT(ex::RuntimeErrorException,
                                  "[zlib] inflate(): unexpected end of stream");
            }
            if (_stream.avail_in != 0) {
                cleanup(FAILED);
                throw LSST_EXCEPT(ex::RuntimeErrorException,
                                  "[zlib] inflate() failed to consume input");
            }
            cleanup(FINISHED);
            return len - _stream.avail_out;
        } else if (zret == Z_OK) {
            if (_stream.avail_out == 0) {
                return len;
            } else if (_remaining == 0) {
                cleanup(FAILED);
                throw LSST_EXCEPT(ex::RuntimeErrorException,
                                  "[zlib] inflate(): expected end of stream");
            } else if (_stream.avail_in != 0) {
                cleanup(FAILED);
                throw LSST_EXCEPT(ex::RuntimeErrorException,
                                  "[zlib] inflate() failed to consume input");
            }
            // user has space for even more data
        } else {
            cleanup(FAILED);
            throw LSST_EXCEPT(ex::RuntimeErrorException,
                (boost::format("[zlib] inflate() failed, return code: %1%") % zret).str());
        }
    }

    // Current block is consumed, so block on outstanding request
    while (::aio_error(&_request) == EINPROGRESS) {
        TimeSpec timeout;
        ::aiocb * const list = &_request;
        timeout.tv_sec = 10; // wait up to 10 seconds
        if (::aio_suspend(&list, 1, &timeout) != 0) {
            if (errno != EINTR) {
                if (errno == EAGAIN) {
                    cleanup(FAILED);
                    throw LSST_EXCEPT(ex::TimeoutException, "aio_read() timed out");
                } else {
                    cleanup(FAILED);
                    throw LSST_EXCEPT(ex::IoErrorException,
                        (boost::format("aio_suspend() failed, errno: %1%") % errno).str());
                }
            }
            errno = 0;
        }
    }

    // check status of the read that just completed
    int const nb = ::aio_return(&_request);
    if (nb < 0) {
        cleanup(FAILED);
        throw LSST_EXCEPT(ex::IoErrorException,
            (boost::format("aio_read() failed, errno: %1%") % ::aio_error(&_request)).str());
    } else if (nb == 0) {
        cleanup(FAILED);
        throw LSST_EXCEPT(ex::IoErrorException, "aio_read(): unexpected end of file reached");
    }

    unsigned char * ready = static_cast<unsigned char *>(const_cast<void *>(_request.aio_buf));
    _remaining -= static_cast<std::size_t>(nb);
    if (_remaining > 0) {
        // issue prefetch for the next block
        _request.aio_buf     = (ready == _buffers) ? ready + _blockSize : _buffers;
        _request.aio_offset  = _fileSize - _remaining;
        _request.aio_nbytes  = std::min(_blockSize, _remaining);
        if (::aio_read(&_request) != 0) {
            cleanup(FAILED);
            throw LSST_EXCEPT(ex::IoErrorException,
                (boost::format("aio_read() failed to enqueue IO request, errno: %1%") % errno).str());
        }
    }

    // Decompress into user-supplied buffer
    _stream.next_in  = ready;
    _stream.avail_in = nb;
    int const zret = ::inflate(&_stream, Z_SYNC_FLUSH);
    if (zret == Z_STREAM_END) {
        if (_remaining != 0) {
            cleanup(FAILED);
            throw LSST_EXCEPT(ex::RuntimeErrorException,
                              "[zlib] inflate(): unexpected end of stream");
        }
        if (_stream.avail_in != 0) {
            cleanup(FAILED);
            throw LSST_EXCEPT(ex::RuntimeErrorException,
                              "[zlib] inflate() failed to consume input");
        }
        cleanup(FINISHED);
    } else if (zret != Z_OK) {
        cleanup(FAILED);
        throw LSST_EXCEPT(ex::RuntimeErrorException,
            (boost::format("[zlib] inflate() failed, return code: %1%") % zret).str());
    } else if (_stream.avail_out != 0) {
        if (_remaining == 0) {
            cleanup(FAILED);
            throw LSST_EXCEPT(ex::RuntimeErrorException,
                              "[zlib] inflate(): expected end of stream");
        } else if (_stream.avail_in != 0) {
            cleanup(FAILED);
            throw LSST_EXCEPT(ex::RuntimeErrorException,
                              "[zlib] inflate() failed to consume input");
        }
    }
    return len - _stream.avail_out;
}


// -- CompressedFileWriter ----------------

lsst::ap::io::CompressedFileWriter::CompressedFileWriter(
    std::string const & fileName,
    bool        const   overwrite,
    std::size_t const   blockSize
) :
    _memory(0),
    _buffers(0),
    _blockSize(blockSize),
    _fd(-1),
    _started(false)
{
    if (blockSize < 1024 || blockSize > 16777216) {
        throw LSST_EXCEPT(ex::InvalidParameterException,
                          "I/O block size must be a power of 2 between 2^10 and 2^24 bytes");
    }

    std::memset(&_stream,  0, sizeof(::z_stream));
    std::memset(&_request, 0, sizeof(::aiocb));
    _request.aio_fildes = -1;

    // allocate IO buffers, store pointer to 8k aligned location inside buffer
    boost::scoped_array<unsigned char> mem(new unsigned char[2*blockSize + 8192]);
    std::size_t buf = reinterpret_cast<std::size_t>(mem.get());
    _buffers = reinterpret_cast<unsigned char *>((buf + 8191) & ~static_cast<std::size_t>(8191));

    // open file
    int const fd = openFile(
        fileName,
        O_WRONLY | O_CREAT | O_APPEND | (overwrite ? O_TRUNC : O_EXCL),
        S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH
    );
    assert(fd != -1);
    ScopeGuard fdGuard(boost::bind(::close, fd));

    // setup zlib (15 window bits, add 16 to indicate a gzip compatible header is desired)
    if (deflateInit2(&_stream, 1, Z_DEFLATED, MAX_WBITS + 16, MAX_MEM_LEVEL, Z_DEFAULT_STRATEGY) != Z_OK) {
        throw LSST_EXCEPT(ex::RuntimeErrorException,
                          "[zlib] deflateInit2() failed to initialize compression stream");
    }

    _request.aio_fildes = fd;
    _request.aio_buf    = _buffers;
    _request.aio_nbytes = blockSize;
    _fd                 = fd;
    using std::swap;
    swap(_memory, mem);
    fdGuard.dismiss();
}


lsst::ap::io::CompressedFileWriter::~CompressedFileWriter() { cleanup(); }


void lsst::ap::io::CompressedFileWriter::cleanup() {
    if (_fd != -1) {
        if (_started) {
            ::aio_cancel(_fd, 0);
        }
        ::close(_fd);
        ::deflateEnd(&_stream);            
        _fd = -1;
    }
}


void lsst::ap::io::CompressedFileWriter::write(
    unsigned char const * const buf,
    std::size_t const len
) {
    if (len == 0) {
        return;
    }
    if (buf == 0) {
        throw LSST_EXCEPT(ex::InvalidParameterException,
                          "null pointer to read destination");
    }
    if (_state != IN_PROGRESS) {
        throw LSST_EXCEPT(ex::IoErrorException,
                          "write() called on a finished or failed CompressedFileWriter");
    }

    _stream.next_in  = const_cast<Bytef *>(buf);
    _stream.avail_in = len;
    if (_stream.avail_out == 0) {
        unsigned char * buf = static_cast<unsigned char *>(const_cast<void *>(_request.aio_buf));
        _stream.next_out = (buf == _buffers) ? buf + _blockSize : _buffers;
        _stream.avail_out = _blockSize;
    }

    // deflate/write until the buffer passed in by the user has been consumed
    do {

        int const zret = ::deflate(&_stream, Z_NO_FLUSH);
        if (zret != Z_OK) {
            cleanup(FAILED);
            throw LSST_EXCEPT(ex::RuntimeErrorException,
                (boost::format("[zlib] deflate() failed, return code: %1%") % zret).str());
        }
        if (_stream.avail_out != 0) {
            if (_stream.avail_in != 0) {
                cleanup(FAILED);
                throw LSST_EXCEPT(ex::RuntimeErrorException,
                                  "[zlib] deflate() failed to consume input");
            }
            break;
        }

        if (_stream.avail_out == 0) {
            if (_started) {
                // wait for pending write to finish
                while (::aio_error(&_request) == EINPROGRESS) {
                    TimeSpec timeout;
                    ::aiocb * const list = &_request;
                    timeout.tv_sec = 10; // wait up to 10 seconds
                    if (::aio_suspend(&list, 1, &timeout) != 0) {
                        if (errno != EINTR) {
                            if (errno == EAGAIN) {
                                cleanup(FAILED);
                                throw LSST_EXCEPT(ex::TimeoutException, "aio_write() timed out");
                            } else {
                                cleanup(FAILED);
                                throw LSST_EXCEPT(ex::IoErrorException,
                                    (boost::format("aio_suspend() failed, errno: %1%") % errno).str());
                            }
                        }
                        errno = 0;
                    }
                }
                // check status of the write that just completed
                if (::aio_return(&_request) != static_cast<int>(_blockSize)) {
                    cleanup(FAILED);
                    throw LSST_EXCEPT(ex::IoErrorException,
                        (boost::format("aio_write() failed, errno: %1%") %
                            ::aio_error(&_request)).str());
                }
            }
            unsigned char * buf = static_cast<unsigned char *>(const_cast<void *>(_request.aio_buf));
            _request.aio_buf  = _stream.next_out - _blockSize;
            _stream.next_out  = buf;
            _stream.avail_out = _blockSize;
            // issue asynchronous write
            if (aio_write(&_request) != 0) {
                cleanup(FAILED);
                throw LSST_EXCEPT(ex::IoErrorException,
                    (boost::format("aio_write() failed to enqueue IO request") % errno).str());
            }
            _started = true;
        } else if (_stream.avail_in != 0) {
            throw LSST_EXCEPT(ex::RuntimeErrorException,
                              "[zlib] deflate() failed to consume input");
        }

    } while(_stream.avail_in != 0);
}


void lsst::ap::io::CompressedFileWriter::finish() {

    if (_state != IN_PROGRESS) {
        throw LSST_EXCEPT(ex::IoErrorException,
            "finish() called on a finished or failed CompressedFileWriter");
    }

    if (_stream.avail_out == 0) {
        unsigned char * buf = static_cast<unsigned char *>(const_cast<void *>(_request.aio_buf));
        _stream.next_out  = (buf == _buffers) ? buf + _blockSize : _buffers;
        _stream.avail_out = _blockSize;
    }

    while (true) {

        int const zret = ::deflate(&_stream, Z_FINISH);
        if (zret != Z_STREAM_END && zret != Z_OK) {
            cleanup(FAILED);
            throw LSST_EXCEPT(ex::RuntimeErrorException,
                (boost::format("[zlib] deflate() failed, return code: %1%") % zret).str());
        }
        if (_started) {
            while (::aio_error(&_request) == EINPROGRESS) {
                 TimeSpec timeout;
                 ::aiocb * const list = &_request;
                 timeout.tv_sec = 10; // wait up to 10 seconds
                 if (::aio_suspend(&list, 1, &timeout) != 0) {
                     if (errno != EINTR) {
                         if (errno == EAGAIN) {
                             cleanup(FAILED);
                             throw LSST_EXCEPT(ex::TimeoutException, "aio_write() timed out");
                         } else {
                             cleanup(FAILED);
                             throw LSST_EXCEPT(ex::IoErrorException,
                                 (boost::format("aio_suspend() failed, errno: %1%") % errno).str());
                         }
                     }
                     errno = 0;
                 }
            }
            // check status of the write that just completed
            if (::aio_return(&_request) != static_cast<int>(_blockSize)) {
                cleanup(FAILED);
                throw LSST_EXCEPT(ex::IoErrorException,
                    (boost::format("aio_write() failed, errno: %1%") %
                        ::aio_error(&_request)).str());
            }
        }
        if (zret == Z_STREAM_END) {
            break;
        }

        unsigned char * buf = static_cast<unsigned char *>(const_cast<void *>(_request.aio_buf));
        _stream.next_out  = buf;
        _stream.avail_out = _blockSize;
        // issue asynchronous write of the deflated data
        _request.aio_buf = (buf == _buffers) ? buf + _blockSize : _buffers;
        if (::aio_write(&_request) != 0) {
            cleanup(FAILED);
            throw LSST_EXCEPT(ex::IoErrorException,
                (boost::format("aio_write() failed to enqueue IO request") % errno).str());
        }
        _started = true;

    }

    // no outstanding IO remains, perform final (blocking) write
    std::size_t bytesRemaining = _blockSize - _stream.avail_out;
    unsigned char * buf = _stream.next_out - bytesRemaining;
    while (bytesRemaining > 0) {
        ::ssize_t nb = ::write(_fd, buf, bytesRemaining);
        if (nb < 0) {
            cleanup(FAILED);
            throw LSST_EXCEPT(ex::IoErrorException,
                              (boost::format("write() failed, errno: %1%") % errno).str());
        }
        bytesRemaining -= static_cast<std::size_t>(nb);
        buf += nb;
    }
    // flush both userland and kernel buffers
    while (::fsync(_fd) != 0) {
        if (errno != EINTR) {
            cleanup(FAILED);
            throw LSST_EXCEPT(ex::IoErrorException,
                              (boost::format("fsync() failed, errno: %1%") % errno).str());
        }
        errno = 0;
    }
    cleanup(FINISHED);
}

