// -*- lsst-c++ -*-

/**
 * @file
 * @brief   Implementation of low level file IO classes.
 *
 * @ingroup associate
 */

#include <sys/types.h>
#include <sys/stat.h>
#if LSST_AP_HAVE_DIRECTIO
#   include <sys/fcntl.h> // for ::directio()
#endif
#include <fcntl.h>
#include <aio.h>
#include <errno.h>

#include <zlib.h>
#if ZLIB_VERNUM < 0x123
#    warning Older version of zlib detected, upgrading to version 1.2.3 or later is recommended
#endif

#include <cstring>
#include <algorithm>

#include <boost/bind.hpp>

#include <lsst/ap/Common.h>
#include <lsst/ap/Exceptions.h>
#include <lsst/ap/ScopeGuard.h>
#include <lsst/ap/Time.h>
#include <lsst/ap/io/FileIo.h>


namespace lsst {
namespace ap {
namespace io {

namespace {

int openFile(
    std::string const & fileName,
    bool        const   direct,
    int         const   oflag,
    ::mode_t    const   mode = 0
) {

    int o = oflag;
//#if LSST_AP_HAVE_O_NOATIME // Linux: turn off last-access-time updates for the file (save on seeks)
//    o |= O_NOATIME;
//#endif

    // Open file
    int fd = ::open(fileName.c_str(), o, mode);

    if (fd == -1) {
        if (errno != ENOENT || (oflag & O_WRONLY) != 0) {
            LSST_AP_THROW_ERR(IoError, boost::format("open(): failed to open file %1%, flags: %2%, errno: %3%") % fileName % oflag % errno, errno);
        }
        // file didn't exist (ok when opening for reading)
        return -1;
    }
    ScopeGuard fdGuard(boost::bind(::close, fd));

#if LSST_AP_HAVE_DIRECTIO
    if (direct) {
        // Solaris: enable direct io
        if (::directio(fd, DIRECTIO_ON) != 0) {
            LSST_AP_THROW_ERR(
                IoError,
                boost::format("directio(): failed to enable direct IO on file %1%") % fileName,
                errno
            );
        }
    }
#elif LSST_AP_HAVE_F_NOCACHE
    if (direct) {
        // Darwin: disable file caching
        if (::fcntl(fd, F_NOCACHE, 1) == -1) {
            LSST_AP_THROW_ERR(
                IoError,
                boost::format("fcntl(): failed to disable caching for file %1%") % fileName,
                errno
            );
        }
    }
#endif
    fdGuard.dismiss();
    return fd;
}

} // end of anonymous namespace


// -- SequentialIoBase ----------------

SequentialIoBase::SequentialIoBase() : _state(IN_PROGRESS) {}

SequentialIoBase::~SequentialIoBase() {}


// -- SequentialFileReader ----------------

SequentialFileReader::SequentialFileReader(
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


SequentialFileReader::~SequentialFileReader() { cleanup(); }


void SequentialFileReader::cleanup() {
    if (_fd != -1) {
        ::close(_fd);
        _fd = -1;
    }
}


size_t SequentialFileReader::read(uint8_t * const buf, size_t const len) {
    if (len == 0) {
        return 0;
    }
    if (buf == 0) {
        LSST_AP_THROW(InvalidParameter, "null pointer to read destination");
    }
    if (_state == FAILED) {
        LSST_AP_THROW(IoError, "read() called on a failed SequentialFileReader");
    } else if (_state == FINISHED) {
        return 0;
    }

    ssize_t n = ::read(_fd, buf, len);
    if (n < 0) {
        cleanup(FAILED);
        LSST_AP_THROW_ERR(IoError, "read(): failed", errno);
    } else if (n == 0) {
        cleanup(FINISHED);
    }
    return static_cast<size_t>(n);
}


// -- SequentialFileWriter ----------------

SequentialFileWriter::SequentialFileWriter(
    std::string const & fileName,
    bool        const   overwrite
) :
    _fd(-1)
{
    int const fd = openFile(
        fileName,
        false,
        O_WRONLY | O_CREAT | O_APPEND | (overwrite ? O_TRUNC : O_EXCL),
        S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH
    );
    assert(fd != -1);
    _fd = fd;
}


SequentialFileWriter::~SequentialFileWriter() { cleanup(); }


void SequentialFileWriter::cleanup() {
    if (_fd != -1) {
        ::close(_fd);
    }
}


void SequentialFileWriter::write(uint8_t const * const buf, size_t const len) {
    if (len == 0) {
        return;
    }
    if (buf == 0) {
        LSST_AP_THROW(InvalidParameter, "null pointer to bytes to write");
    }
    if (_state != IN_PROGRESS) {
        LSST_AP_THROW(IoError, "write() called on a finished or failed SequentialFileWriter");
    }

    uint8_t const * dst = buf;
    size_t          nb  = len;
    do {
        ssize_t n = ::write(_fd, dst, nb);
        if (n < 0) {
            cleanup(FAILED);
            LSST_AP_THROW_ERR(IoError, "write(): failed", errno);
        }
        dst += n;
        nb  -= n;
    } while (nb > 0);
}


void SequentialFileWriter::finish() {
    if (_state != IN_PROGRESS) {
        LSST_AP_THROW(IoError, "finish() called on a finished or failed SequentialFileWriter");
    }

    // Flush userland and kernel buffers
    while (fsync(_fd) != 0) {
        if (errno != EINTR) {
            cleanup(FAILED);
            LSST_AP_THROW_ERR(IoError, "fsync(): failed", errno);
        }
        errno = 0;
    }
    cleanup(FINISHED);
}


// -- CompressedFileReader ----------------

CompressedFileReader::CompressedFileReader(
    std::string const & fileName,
    bool        const   direct,
    size_t      const   blockSize
) :
    _memory(0),
    _buffers(0),
    _blockSize(blockSize),
    _fileSize(0),
    _remaining(0),
    _fd(-1)
{
    if (blockSize < 1024 || blockSize > 16777216 || (blockSize & (blockSize - 1)) != 0) {
        LSST_AP_THROW(
            InvalidParameter,
            "I/O block size must be a power of 2 between 1024 and 16777216 bytes"
        );
    }

    std::memset(&_stream,  0, sizeof(::z_stream));
    std::memset(&_request, 0, sizeof(::aiocb));
    _request.aio_fildes = -1;

    // allocate IO buffers
    boost::scoped_array<uint8_t> mem(new uint8_t[2*blockSize + 8192]);
    size_t buf = reinterpret_cast<size_t>(mem.get());
    _buffers = reinterpret_cast<uint8_t *>((buf + 8191) & ~static_cast<size_t>(8191));

    // open file
    int const fd = openFile(fileName, direct, O_RDONLY);
    if (fd == -1) {
        _state = FINISHED;
        return;
    }
    ScopeGuard fdGuard(boost::bind(::close, fd));

    // determine file size
    struct stat sb;
    if (::fstat(fd, &sb) != 0) {
        LSST_AP_THROW_ERR(
            IoError,
            boost::format("fstat(): failed on file %1%") % fileName,
            errno
        );
    }
    _fileSize = static_cast<size_t>(sb.st_size);
    if (_fileSize == 0) {
        _state = FINISHED;
        return;
    }
    _remaining = _fileSize;

    // setup zlib (specify auto-detection of zlib/gzip header)
    if (inflateInit2(&_stream, MAX_WBITS + 32) != Z_OK) {
        LSST_AP_THROW(Runtime, "[zlib] inflateInit2(): failed to initialize decompression stream");
    }
    ScopeGuard zlibGuard(boost::bind(::inflateEnd, &_stream));

    // issue first async read
    size_t const nb = std::min(blockSize, _fileSize);
    _request.aio_fildes = fd;
    _request.aio_buf    = _buffers;
    _request.aio_nbytes = nb;
    if (aio_read(&_request) != 0) {
        LSST_AP_THROW_ERR(
            IoError,
            boost::format("aio_read(): failed to enqueue IO request for file %1%") % fileName,
            errno
        );
    }
    _fd = fd;
    using std::swap;
    swap(_memory, mem);
    zlibGuard.dismiss();
    fdGuard.dismiss();
}


CompressedFileReader::~CompressedFileReader() { cleanup(); }


void CompressedFileReader::cleanup() {
    if (_fd != -1) {
        if (_state == IN_PROGRESS) {
            ::aio_cancel(_fd, 0);
        }
        ::close(_fd);
        ::inflateEnd(&_stream);
        _fd = -1;
    }
}


size_t CompressedFileReader::read(uint8_t * const buf, size_t const len) {

    if (len == 0) {
        return 0;
    }
    if (buf == 0) {
        LSST_AP_THROW(InvalidParameter, "null pointer to read destination");
    }
    if (_state == FAILED) {
        LSST_AP_THROW(IoError, "read() called on a failed CompressedFileReader");
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
                LSST_AP_THROW(Runtime, "[zlib] inflate(): unexpected end of stream");
            }
            if (_stream.avail_in != 0) {
                cleanup(FAILED);
                LSST_AP_THROW(Runtime, "[zlib] inflate(): failed to consume input");
            }
            cleanup(FINISHED);
            return len - _stream.avail_out;
        } else if (zret == Z_OK) {
            if (_stream.avail_out == 0) {
                return len;
            } else if (_remaining == 0) {
                cleanup(FAILED);
                LSST_AP_THROW(Runtime, "[zlib] inflate(): expected end of stream");
            } else if (_stream.avail_in != 0) {
                cleanup(FAILED);
                LSST_AP_THROW(Runtime, "[zlib] inflate(): failed to consume input");
            }
            // user has space for even more data
        } else {
            cleanup(FAILED);
            LSST_AP_THROW(Runtime, boost::format("[zlib] inflate(): failed (error code %1%)") % zret);
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
                    LSST_AP_THROW(Timeout, "aio_read(): timed out");
                } else {
                    cleanup(FAILED);
                    LSST_AP_THROW_ERR(IoError, "aio_suspend(): failed", errno);
                }
            }
            errno = 0;
        }
    }

    // check status of the read that just completed
    int const nb = ::aio_return(&_request);
    if (nb < 0) {
        cleanup(FAILED);
        LSST_AP_THROW_ERR(IoError, "aio_read(): read failed", aio_error(&_request));
    } else if (nb == 0) {
        cleanup(FAILED);
        LSST_AP_THROW(IoError, "aio_read(): unexpected end of file reached");
    }

    uint8_t * ready = static_cast<uint8_t *>(const_cast<void *>(_request.aio_buf));
    _remaining -= static_cast<size_t>(nb);
    if (_remaining > 0) {
        // issue prefetch for the next block
        _request.aio_buf     = (ready == _buffers) ? ready + _blockSize : _buffers;
        _request.aio_offset  = _fileSize - _remaining;
        _request.aio_nbytes  = std::min(_blockSize, _remaining);
        if (aio_read(&_request) != 0) {
            cleanup(FAILED);
            LSST_AP_THROW_ERR(IoError, "aio_read(): failed to enqueue IO request", errno);
        }
    }

    // Decompress into user-supplied buffer
    _stream.next_in  = ready;
    _stream.avail_in = nb;
    int const zret = ::inflate(&_stream, Z_SYNC_FLUSH);
    if (zret == Z_STREAM_END) {
        if (_remaining != 0) {
            cleanup(FAILED);
            LSST_AP_THROW(Runtime, "[zlib] inflate(): unexpected end of stream");
        }
        if (_stream.avail_in != 0) {
            cleanup(FAILED);
            LSST_AP_THROW(Runtime, "[zlib] inflate(): failed to consume input");
        }
        cleanup(FINISHED);
    } else if (zret != Z_OK) {
        cleanup(FAILED);
        LSST_AP_THROW(Runtime, boost::format("[zlib] inflate(): failed (error code %1%)") % zret);
    } else if (_stream.avail_out != 0) {
        if (_remaining == 0) {
            cleanup(FAILED);
            LSST_AP_THROW(Runtime, "[zlib] inflate(): expected end of stream");
        } else if (_stream.avail_in != 0) {
            cleanup(FAILED);
            LSST_AP_THROW(Runtime, "[zlib] inflate(): failed to consume input");
        }
    }
    return len - _stream.avail_out;
}


// -- CompressedFileWriter ----------------

CompressedFileWriter::CompressedFileWriter(
    std::string const & fileName,
    bool        const   overwrite,
    bool        const   direct,
    size_t      const   blockSize
) :
    _memory(0),
    _buffers(0),
    _blockSize(blockSize),
    _fd(-1),
    _started(false)
{
    if (blockSize < 1024 || blockSize > 16777216) {
        LSST_AP_THROW(InvalidParameter, "I/O block size must be between 1024 and 16777216 bytes");
    }

    std::memset(&_stream,  0, sizeof(::z_stream));
    std::memset(&_request, 0, sizeof(::aiocb));
    _request.aio_fildes = -1;

    // allocate IO buffers
    boost::scoped_array<uint8_t> mem(new uint8_t[2*blockSize + 8192]);
    size_t buf = reinterpret_cast<size_t>(mem.get());
    _buffers = reinterpret_cast<uint8_t *>((buf + 8191) & ~static_cast<size_t>(8191));

    // open file
    int const fd = openFile(
        fileName,
        direct,
        O_WRONLY | O_CREAT | O_APPEND | (overwrite ? O_TRUNC : O_EXCL),
        S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH
    );
    assert(fd != -1);
    ScopeGuard fdGuard(boost::bind(::close, fd));

    // setup zlib (15 window bits, add 16 to indicate a gzip compatible header is desired)
    if (deflateInit2(&_stream, 1, Z_DEFLATED, MAX_WBITS + 16, MAX_MEM_LEVEL, Z_DEFAULT_STRATEGY) != Z_OK) {
        LSST_AP_THROW(Runtime, "[zlib] deflateInit2(): failed to initialize decompression stream");
    }

    _request.aio_fildes = fd;
    _request.aio_buf    = _buffers;
    _request.aio_nbytes = blockSize;
    _fd                 = fd;
    using std::swap;
    swap(_memory, mem);
    fdGuard.dismiss();
}


CompressedFileWriter::~CompressedFileWriter() { cleanup(); }


void CompressedFileWriter::cleanup() {
    if (_fd != -1) {
        if (_started) {
            ::aio_cancel(_fd, 0);
        }
        ::close(_fd);
        ::deflateEnd(&_stream);            
        _fd = -1;
    }
}


void CompressedFileWriter::write(uint8_t const * const buf, size_t const len) {

    if (len == 0) {
        return;
    }
    if (buf == 0) {
        LSST_AP_THROW(InvalidParameter, "null pointer to read destination");
    }
    if (_state != IN_PROGRESS) {
        LSST_AP_THROW(IoError, "write() called on a finished or failed CompressedFileWriter");
    }

    _stream.next_in   = const_cast<Bytef *>(buf);
    _stream.avail_in  = len;
    if (_stream.avail_out == 0) {
        uint8_t * buf     = static_cast<uint8_t *>(const_cast<void *>(_request.aio_buf));
        _stream.next_out  = (buf == _buffers) ? buf + _blockSize : _buffers;
        _stream.avail_out = _blockSize;
    }

    // deflate/write until the buffer passed in by the user has been consumed
    do {

        int const zret = ::deflate(&_stream, Z_NO_FLUSH);
        if (zret != Z_OK) {
            cleanup(FAILED);
            LSST_AP_THROW(Runtime, boost::format("[zlib] deflate(): failed (error code %1%)") % zret);
        }
        if (_stream.avail_out != 0) {
            if (_stream.avail_in != 0) {
                cleanup(FAILED);
                LSST_AP_THROW(Runtime, "[zlib] deflate(): failed to consume input");
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
                                LSST_AP_THROW(Timeout, "aio_write(): timed out");
                            } else {
                                cleanup(FAILED);
                                LSST_AP_THROW_ERR(IoError, "aio_suspend(): failed", errno);
                            }
                        }
                        errno = 0;
                    }
                }
                // check status of the write that just completed
                if (::aio_return(&_request) != static_cast<int>(_blockSize)) {
                    cleanup(FAILED);
                    LSST_AP_THROW_ERR(IoError, "aio_write(): write failed", aio_error(&_request));
                }
            }
            uint8_t * buf     = static_cast<uint8_t *>(const_cast<void *>(_request.aio_buf));
            _request.aio_buf  = _stream.next_out - _blockSize;
            _stream.next_out  = buf;
            _stream.avail_out = _blockSize;
            // issue asynchronous write
            if (aio_write(&_request) != 0) {
                cleanup(FAILED);
                LSST_AP_THROW_ERR(IoError, "aio_write() failed to enqueue IO request", errno);
            }
            _started = true;
        } else if (_stream.avail_in != 0) {
            LSST_AP_THROW(Runtime, "[zlib] deflate(): failed to consume input");
        }

    } while(_stream.avail_in != 0);
}


void CompressedFileWriter::finish() {

    if (_state != IN_PROGRESS) {
        LSST_AP_THROW(IoError, "finish() called on a finished or failed CompressedFileWriter");
    }

    if (_stream.avail_out == 0) {
        uint8_t * buf     = static_cast<uint8_t *>(const_cast<void *>(_request.aio_buf));
        _stream.next_out  = (buf == _buffers) ? buf + _blockSize : _buffers;
        _stream.avail_out = _blockSize;
    }

    while (true) {

        int const zret = ::deflate(&_stream, Z_FINISH);
        if (zret != Z_STREAM_END && zret != Z_OK) {
            cleanup(FAILED);
            LSST_AP_THROW(Runtime, boost::format("[zlib] deflate(): failed (error code %1%)") % zret);
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
                             LSST_AP_THROW(Timeout, "aio_write(): timed out");
                         } else {
                             cleanup(FAILED);
                             LSST_AP_THROW_ERR(IoError, "aio_suspend(): failed", errno);
                         }
                     }
                     errno = 0;
                 }
            }
            // check status of the write that just completed
            if (::aio_return(&_request) != static_cast<int>(_blockSize)) {
                cleanup(FAILED);
                LSST_AP_THROW_ERR(IoError, "aio_write(): write failed", aio_error(&_request));
            }
        }
        if (zret == Z_STREAM_END) {
            break;
        }

        uint8_t * buf = static_cast<uint8_t *>(const_cast<void *>(_request.aio_buf));
        _stream.next_out  = buf;
        _stream.avail_out = _blockSize;
        // issue asynchronous write of the deflated data
        _request.aio_buf = (buf == _buffers) ? buf + _blockSize : _buffers;
        if (aio_write(&_request) != 0) {
            cleanup(FAILED);
            LSST_AP_THROW_ERR(IoError, "aio_write() failed to enqueue IO request", errno);
        }
        _started = true;

    }

    // no outstanding IO remains, perform final (blocking) write
    size_t    bytesRemaining = _blockSize - _stream.avail_out;
    uint8_t * buf            = _stream.next_out - bytesRemaining;
    while (bytesRemaining > 0) {
        ssize_t nb = ::write(_fd, buf, bytesRemaining);
        if (nb < 0) {
            cleanup(FAILED);
            LSST_AP_THROW_ERR(IoError, "write(): failed", errno);
        }
        bytesRemaining -= static_cast<size_t>(nb);
        buf            += nb;
    }
    // flush both userland and kernel buffers
    while (fsync(_fd) != 0) {
        if (errno != EINTR) {
            cleanup(FAILED);
            LSST_AP_THROW_ERR(IoError, "fsync(): failed", errno);
        }
        errno = 0;
    }
    cleanup(FINISHED);
}


}}} // end of namespace lsst::ap::io
