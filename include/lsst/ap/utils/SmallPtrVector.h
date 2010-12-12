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
  * @brief  Pointer vector that avoids memory allocation for small vectors.
  * @author Serge Monkewitz
  */
#ifndef LSST_AP_UTILS_SMALLPTRVECTOR_H
#define LSST_AP_UTILS_SMALLPTRVECTOR_H


namespace lsst { namespace ap { namespace utils {

/** A class for storing vectors of pointers.
  * Up to N pointers can be held without incurring any memory allocations.
  * The interface is a subset of the std::vector interface - as a result a
  * std::vector can be substituted with minimal code changes.
  */
template <typename T, size_t N>
class LSST_AP_LOCAL SmallPtrVector {
public:
    typedef T ** iterator;
    typedef T const * const * const_iterator;
    typedef T *& reference;
    typedef T const * const & const_reference;
    typedef size_t size_type;
    typedef ptrdiff_t difference_type;
    typedef T * value_type;

    inline SmallPtrVector();
    inline ~SmallPtrVector(); 

    template <size_t N2>
    inline SmallPtrVector(SmallPtrVector<T, N2> const &v);

    template <size_t N2>
    SmallPtrVector & operator=(SmallPtrVector<T, N2> const &v);

    inline size_type size() const;
    inline size_type capacity() const;
    inline size_type max_size() const;
    inline bool empty() const;

    inline const_iterator begin() const;
    inline iterator begin();
    inline const_iterator end() const;
    inline iterator end();

    inline const_reference operator[](size_type i) const;
    inline reference operator[](size_type i);

    inline const_reference front() const;
    inline reference front();
    inline const_reference back() const;
    inline reference back();

    inline void clear();
    inline void pop_back();
    inline void push_back(value_type v);

    void reserve(size_type n);
    void swap(SmallPtrVector &v);

private:
    void _grow();

    T * _buf[N];
    T ** _beg;
    T ** _end;
    T ** _cap;
};

}}} // namespace lsst::ap::utils

#include "SmallPtrVector.cc"

#endif // LSST_AP_UTILS_SMALLPTRVECTOR_H

