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
 * @brief   Compile time constants related to data types to be stored in chunks.
 *
 * @ingroup ap
 */

#ifndef LSST_AP_DATA_TRAITS_H
#define LSST_AP_DATA_TRAITS_H

#include "Common.h"
#include "Object.h"


namespace lsst { namespace ap {

/**
 * @brief  Provides basic chunk parameters at compile time.
 *
 * Specializations of DataTraits must provide the following data type specific parameters:
 * <dl>
 * <dt><b> ENTRIES_PER_BLOCK_LOG2 </b></dt>
 * <dd> The base 2 logarithm of the number of entries in a single block,
 *      where blocks are the units of memory allocation. </dd>
 * <dt><b> MAX_BLOCKS_PER_CHUNK </b></dt>
 * <dd> The maximum number of blocks that can be allocated to a chunk.
 *      This should be determined from a very conservative (high) estimate
 *      of the worst case number of entries in a chunk. </dd>
 * <dt><b> MAX_CHUNKS_PER_FOV </b></dt>
 * <dd> The maximum number of chunks in a FOV, determined from the FOV size
 *      and partitioning granularity. </dd>
 * <dt><b> NUM_BLOCKS </b></dt>
 * <dd> Total number of allocateable blocks. Determined by calculating the
 *      necessary storage for 4 worst case FOVs. </dd>
 * </dl>
 */
template <typename D> struct DataTraits {};

template <> struct LSST_AP_LOCAL DataTraits<Object> {
    static int const ENTRIES_PER_BLOCK_LOG2 = 12;
    static int const MAX_BLOCKS_PER_CHUNK   = 128;
    static int const MAX_CHUNKS_PER_FOV     = 128;
    static int const NUM_BLOCKS             = 1024;
};


}} // end of namespace lsst::ap

#endif // LSST_AP_DATA_TRAITS_H
