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
  * @brief  SmallPtrVector implementation.
  * @author Serge Monkewitz
  */
#ifndef LSST_AP_UTILS_SMALLPTRVECTOR_CC
#define LSST_AP_UTILS_SMALLPTRVECTOR_CC

#include <cstring>
#include <algorithm>
#include <stdexcept>


namespace lsst { namespace ap { namespace utils {

template <typename T, size_t N>
inline SmallPtrVector<T, N>::SmallPtrVector() :
    _beg(_buf), _end(_buf), _cap(_buf + N) { }

template <typename T, size_t N>
    template <size_t N2>
inline SmallPtrVector<T, N>::SmallPtrVector(SmallPtrVector<T, N2> const &v) :
    _beg(_buf), _end(_buf), _cap(_buf + N)
{
    *this = v;
}

template <typename T, size_t N>
inline SmallPtrVector<T, N>::~SmallPtrVector() {
    if (_beg != _buf) {
        delete[] _beg;
    }
}

template <typename T, size_t N>
inline typename SmallPtrVector<T, N>::size_type SmallPtrVector<T, N>::size() const {
    return static_cast<size_type>(_end - _beg);
}

template <typename T, size_t N>
inline typename SmallPtrVector<T, N>::size_type SmallPtrVector<T, N>::capacity() const {
    return static_cast<size_type>(_cap - _beg);
}

template <typename T, size_t N>
inline typename SmallPtrVector<T, N>::size_type SmallPtrVector<T, N>::max_size() const {
    return static_cast<size_type>(std::numeric_limits<difference_type>::max()) /
           sizeof(value_type);
}

template <typename T, size_t N>
inline bool SmallPtrVector<T, N>::empty() const {
    return size() == 0;
}

template <typename T, size_t N>
inline typename SmallPtrVector<T, N>::const_iterator SmallPtrVector<T, N>::begin() const {
    return _beg;
}

template <typename T, size_t N>
inline typename SmallPtrVector<T, N>::iterator SmallPtrVector<T, N>::begin() {
    return _beg;
}

template <typename T, size_t N>
inline typename SmallPtrVector<T, N>::const_iterator SmallPtrVector<T, N>::end() const {
    return _end;
}

template <typename T, size_t N>
inline typename SmallPtrVector<T, N>::iterator SmallPtrVector<T, N>::end() {
    return _end;
}

template <typename T, size_t N>
inline typename SmallPtrVector<T, N>::const_reference SmallPtrVector<T, N>::operator[](
    typename SmallPtrVector<T, N>::size_type i) const
{
    return *(_beg + i);
}

template <typename T, size_t N>
inline typename SmallPtrVector<T, N>::reference SmallPtrVector<T, N>::operator[](
    typename SmallPtrVector<T, N>::size_type i)
{
    return *(_beg + i);
}

template <typename T, size_t N>
inline typename SmallPtrVector<T, N>::const_reference SmallPtrVector<T, N>::front() const {
    return *_beg;
}

template <typename T, size_t N>
inline typename SmallPtrVector<T, N>::reference SmallPtrVector<T, N>::front() {
    return *_beg;
}

template <typename T, size_t N>
inline typename SmallPtrVector<T, N>::const_reference SmallPtrVector<T, N>::back() const {
    return *(_end - 1);
}

template <typename T, size_t N>
inline typename SmallPtrVector<T, N>::reference SmallPtrVector<T, N>::back() {
    return *(_end - 1);
}

template <typename T, size_t N>
inline void SmallPtrVector<T, N>::clear() {
    _end = _beg;
}

template <typename T, size_t N>
inline void SmallPtrVector<T, N>::pop_back() {
    if (_end > _beg) {
        --end;
    }
}

template <typename T, size_t N>
inline void SmallPtrVector<T, N>::push_back(value_type v) {
    if (_end == _cap) {
        _grow();
    }
    *_end++ = v;
}


template <typename T, size_t N>
    template <size_t N2>
SmallPtrVector<T, N> & SmallPtrVector<T, N>::operator=(
    SmallPtrVector<T, N2> const &v)
{
    if (_beg == v._beg) {
        return *this;
    }
    clear();
    if (v.size() > 0) {
        reserve(v.size());
        std::memcpy(_beg, v._beg, v.size()*sizeof(value_type));
        _end = _beg + v.size();
    }
}

template <typename T, size_t N>
void SmallPtrVector<T, N>::reserve(typename SmallPtrVector<T, N>::size_type n) {
    if (n > max_size()) {
        throw std::bad_alloc();
    } else if (n > capacity()) {
        T **array = new T *[n];
        size_type sz = size();
        if (sz > 0) {
            std::memcpy(array, _beg, sz*sizeof(value_type));
        }
        std::swap(array, _beg);
        _end = _beg + sz;
        _cap = _beg + n;
        if (array != _buf) {
            delete[] array;
        }
    }
}

template <typename T, size_t N>
void SmallPtrVector<T, N>::swap(SmallPtrVector &v) {
    if (this == &v) {
        return;
    }
    if (_beg == _buf && v._beg == v._buf) {
        T *tmp[N];
        size_type sz = size();
        size_type vsz = v.size();
        std::memcpy(tmp, _beg, sz*sizeof(value_type));
        std::memcpy(_beg, v._beg, vsz*sizeof(value_type));
        std::memcpy(v._beg, tmp, sz*sizeof(value_type));
        _end = _beg + vsz;
        v._end = v._beg + sz;
    } else if (v._beg == v._buf) {
        size_type sz = v.size();
        std::memcpy(_buf, v._beg, sz*sizeof(value_type));
        v._beg = _beg;
        v._end = _end;
        v._cap = _cap;
        _beg = _buf;
        _end = _buf + sz;
        _cap = _buf + N;
    } else if (_beg == _buf) {
        size_type sz = size();
        std::memcpy(v._buf, _beg, sz*sizeof(value_type));
        _beg = v._beg;
        _end = v._end;
        _cap = v._cap;
        v._beg = v._buf;
        v._end = v._buf + sz;
        v._cap = v._buf + N;
    } else {
        std::swap(_beg, v._beg);
        std::swap(_end, v._end);
        std::swap(_cap, v._cap);
    }
}

template <typename T, size_t N>
void SmallPtrVector<T, N>::_grow() {
    size_type c = std::min(2*capacity(), max_size());
    if (c <= capacity()) {
        // overflow
        throw std::bad_alloc();
    }
    reserve(c);
}

}}} // namespace lsst::ap::utils

#endif // LSST_AP_UTILS_SMALLPTRVECTOR_CC

