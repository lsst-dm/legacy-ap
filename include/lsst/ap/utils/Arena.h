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
  * @author Serge Monkewitz
  */
#ifndef LSST_AP_UTILS_ARENA_H
#define LSST_AP_UTILS_ARENA_H

#include <cassert>
#include <cstdlib>
#include <algorithm>
#include <stdexcept>
#include <vector>

#include "boost/static_assert.hpp"

#include "lsst/pex/exceptions.h"


namespace lsst { namespace ap { namespace utils {

/** A free-list based single-threaded arena allocator for objects of class T.
  *
  * @par
  * The arena grows in blocks of a run-time specified capacity, and never
  * shrinks during the lifetime of the arena. All blocks are allocated with
  * malloc(), and are freed when the arena is destroyed. Note that destroying
  * the arena will also call the destructors of any live objects in the arena.
  *
  * @par
  * This class is quite low-level and is intended for use in performance
  * critical code only. 
  */
template <typename T>
class Arena {
public:
    Arena(size_t blockCapacity=262144/sizeof(T));
    ~Arena();

    inline void *alloc();
    inline void dealloc(void *ptr);
    inline void destroy(T *ptr);

    inline size_t capacity() const;
    inline size_t getNumBytes() const;
    inline size_t getBlockCapacity() const;

private:
    // T might be or contain a fixed size Eigen type
    static const size_t ALIGN = 16;
    static const size_t SIZE = (sizeof(T) + ALIGN - 1) & ~(ALIGN - 1);

    BOOST_STATIC_ASSERT((ALIGN & (ALIGN - 1)) == 0);

    void _grow();

    static inline unsigned char * _align(unsigned char *p);

    std::vector<unsigned char *> _blocks;   ///< List of memory blocks.
    std::vector<std::vector<bool> > _masks; ///< Per-block free bits - used
                                            ///  only in the arena destructor,
                                            ///  but is preallocated to avoid
                                            ///  std::bad_alloc therin.
    size_t const _blockCapacity;            ///< Capacity of a memory block.
    size_t _nFree;         ///< Number of free elements.
    unsigned char *_free;  ///< Head of linked free-list, 0 if arena is full.
};

}}} // namespace lsst::ap::utils

#include "Arena.cc"


template <typename T>
inline void * operator new(size_t, lsst::ap::utils::Arena<T> &arena) {
   return arena.alloc();
}

template <typename T>
inline void operator delete(void *ptr, lsst::ap::utils::Arena<T> &arena) {
   arena.dealloc(ptr);
}

#endif // LSST_AP_UTILS_ARENA_H

