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
#ifndef LSST_AP_UTILS_ARENA_CC
#define LSST_AP_UTILS_ARENA_CC

#include "Arena.h"

namespace lsst { namespace ap { namespace utils {

/** Creates a new Arena for objects of type T. The arena allocates
  * memory from the system in blocks, each of which is large enough
  * to contain @a blockCapacity objects.
  */
template <typename T>
Arena<T>::Arena(size_t blockCapacity) :
    _blocks(),
    _masks(),
    _blockCapacity(blockCapacity),
    _nFree(0), 
    _free(0)
{
    if (blockCapacity == 0) {
        throw LSST_EXCEPT(lsst::pex::exceptions::InvalidParameterException,
                          "Cannot create Arena with 0 capacity");
    }
    _grow();
}

template <typename T>
Arena<T>::~Arena() {
    typedef std::vector<unsigned char *>::iterator Iter;

    if (_nFree == 0) {
        assert(_free == 0 &&
               "No free elements, but free list head is non-null");
        // call destructor for every element of every block
        for (Iter i = _blocks.begin(), e = _blocks.end(); i != e; ++i) {
            unsigned char *block = *i;
            unsigned char *blockEnd = block + _blockCapacity*SIZE;
            for (; block < blockEnd; block += SIZE) {
                reinterpret_cast<T *>(block)->~T();
            }
        }
    } else if (_nFree < _blocks.size()*_blockCapacity) {
       assert(_free != 0 &&
              "Free count is non zero, but free list head is null");
       // need to call destructor for each allocated entry -
       // first compute free bits for every block.
       std::sort(_blocks.begin(), _blocks.end());
       unsigned char *ptr = _free;
       while (ptr != 0) {
           // figure out which block ptr is from and mark the element free
           Iter b = std::upper_bound(_blocks.begin(), _blocks.end(), ptr);
           assert(b != _blocks.begin() &&
                  "Free list contains pointer not belonging to arena");
           --b;
           assert((ptr - *b) % SIZE == 0 &&
                  "Free list contains unaligned pointer");
           _masks[b - _blocks.begin()][(ptr - *b)/SIZE] = true;
           --_nFree;
           ptr = *reinterpret_cast<unsigned char **>(ptr);
       }
       assert(_nFree == 0 && "Free list does not contain all free elements");
       // free any element not flagged as free in its occupancy mask
       for (size_t i = 0; i < _blocks.size(); ++i) {
           unsigned char *block = _blocks[i];
           std::vector<bool> const & free = _masks[i];
           for (size_t j = 0; j < _blocks.size(); block += SIZE, ++j) {
               if (!free[j]) {
                   reinterpret_cast<T *>(block)->~T();
               }
           }
       }
    } else {
       assert(_nFree == _blocks.size()*_blockCapacity &&
              "Number of free elements exceeds arena capacity");
    }
    // free the memory for each block
    for (Iter i = _blocks.begin(), e = _blocks.end(); i != e; ++i) {
        std::free(*i);
        *i = 0;
    }
}

/** Allocates and returns a pointer to a memory region large enough to
  * hold a single object of type T.
  */
template <typename T>
inline void *Arena<T>::alloc() {
    if (_free == 0) {
        _grow();
    }
    unsigned char *next = _free;
    _free = *reinterpret_cast<unsigned char **>(next);
    --_nFree;
    return next;
}

/** Frees the memory at the given address.
  */
template <typename T>
inline void Arena<T>::dealloc(void *ptr) {
    unsigned char *next = static_cast<unsigned char *>(ptr);
    *reinterpret_cast<unsigned char **>(next) = _free;
    _free = next;
    ++_nFree;
}

/** Destroys and frees the object at the given address.
  */
template <typename T>
inline void Arena<T>::destroy(T *ptr) {
    ptr->~T();
    dealloc(ptr);
}

/** Returns the total number of objects that can fit in the arena.
  */
template <typename T>
inline size_t Arena<T>::capacity() const {
    return _blockCapacity*_blocks.size();
}

/** Returns the number of bytes in the arena.
  */
template <typename T>
inline size_t Arena<T>::getNumBytes() const {
    return capacity()*SIZE;
}

/** Returns the number of objects that fit in a single arena block.
  */
template <typename T>
inline size_t Arena<T>::getBlockCapacity() const {
    return _blockCapacity;
}


template <typename T>
void Arena<T>::_grow() {
    // pre-allocate space in block and block mask list
    _blocks.reserve(_blocks.size() + 1);
    _masks.reserve(_masks.size() + 1);
    std::vector<bool> mask(_blockCapacity, false);
    // allocate another block
    unsigned char *block = static_cast<unsigned char *>(
        std::calloc(_blockCapacity, SIZE));
    if (block == 0) {
        throw std::bad_alloc();
    }
    // what follows will never throw
    _masks.push_back(std::vector<bool>());
    std::swap(_masks.back(), mask);
    // initialize free-list
    unsigned char *item = block;
    for (size_t i = 0; i < _blockCapacity - 1; ++i, item += SIZE) {
        *reinterpret_cast<unsigned char **>(item) = item + SIZE;
    }
    // set _free to head of free-list and store new block pointer
    _free = block;
    _nFree += _blockCapacity;
    _blocks.push_back(block);
}

}}} // namespace lsst::ap::utils

#endif // LSST_AP_UTILS_ARENA_CC

