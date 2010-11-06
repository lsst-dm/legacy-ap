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

/** @file
  * @brief  Single threaded arena (pool) memory allocator.
  */
#ifndef LSST_AP_UTILS_ARENA_H
#define LSST_AP_UTILS_ARENA_H

#include <cstdlib>
#include <vector>

#include "lsst/pex/exceptions.h"


namespace lsst { namespace ap { namespace utils {

/** A free-list based single-threaded arena allocator for objects of class T.
  *
  * @par
  * The arena grows in blocks of a run-time specified capacity, and never
  * shrinks during the lifetime of the arena. All blocks are allocated with
  * calloc(), and are freed when the arena is destroyed. Note that destroying
  * the arena does not result in destructor calls.
  *
  * @par
  * This class is quite low-level in that clients are required to explicitly
  * call destructors and deallocate from the arena. It is intended for use
  * in performance critical code only.
  */
template <typename T>
class Arena {
public:
    /** Creates a new Arena for objects of type T. The arena allocates
      * memory from the system in blocks, each of which is large enough
      * to contain @a blockCapacity objects.
      */
    Arena(size_t blockCapacity=262144/sizeof(T)) :
        _blocks(), _blockCapacity(blockCapacity), _free(0)
    {
        if (blockCapacity == 0) {
            throw LSST_EXCEPT(lsst::pex::exceptions::InvalidParameterException,
                              "Cannot create Arena with 0 capacity");
        }
        _grow();
    }

    ~Arena() {
        typedef std::vector<unsigned char *>::iterator Iter;
        for (Iter i = _blocks.begin(), e = _blocks.end(); i != e; ++i) {
            std::free(*i);
            *i = 0;
        }
    }

    /** Allocates and returns a pointer to a memory region large enough to
      * hold a single object of type T.
      */
    void *alloc() {
        if (_free == 0) {
            _grow();
        }
        unsigned char *next = _free;
        _free = *reinterpret_cast<unsigned char **>(next);
        return next;
    }

    /** Frees the memory at the given address.
      */
    void dealloc(void *ptr) {
        unsigned char *next = static_cast<unsigned char *>(ptr);
        *reinterpret_cast<unsigned char **>(next) = _free;
        _free = next;
    }

    /** Returns the total number of objects that can fit in the arena.
      */
    size_t capacity() const {
        return _blockCapacity*_blocks.size();
    }

    /** Returns the number of bytes in the arena.
      */
    size_t getNumBytes() const {
        return capacity()*SIZE;
    }

private:
    static const size_t REM = sizeof(T) % sizeof(unsigned char *);
    static const size_t PAD = (REM == 0) ? 0 : sizeof(unsigned char *) - REM;
    static const size_t SIZE = sizeof(T) + PAD;

    void _grow() {
        // pre-allocate space in block list
        _blocks.reserve(1);
        // allocate another block
        unsigned char *block = static_cast<unsigned char *>(
            std::calloc(_blockCapacity, SIZE));
        if (block == 0) {
            throw std::bad_alloc();
        }
        // what follows will never throw
        // initialize free-list
        unsigned char *item = block;
        for (size_t i = 0; i < _blockCapacity - 1; ++i, item += SIZE) {
            *reinterpret_cast<unsigned char **>(item) = item + SIZE;
        }
        // set _free to head of free-list and store new block pointer
        _free = block;
        _blocks.push_back(block);
    }

    std::vector<unsigned char *> _blocks; ///< List of memory blocks.
    size_t _blockCapacity; ///< Capacity of a memory block.
    unsigned char *_free;  ///< Head of linked free-list or 0 if the arena is full.
};

}}} // namespace lsst::ap::utils

template <typename T>
inline void * operator new(size_t, lsst::ap::utils::Arena<T> &arena) {
   return arena.alloc();
}

template <typename T>
inline void operator delete(void *ptr, lsst::ap::utils::Arena<T> &arena) {
   arena.dealloc(ptr);
}

#endif // LSST_AP_UTILS_ARENA_H

