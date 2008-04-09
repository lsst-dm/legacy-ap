// -*- lsst-c++ -*-

/**
 * @file
 * @brief   Compile time constants related to data types to be stored in chunks.
 *
 * @ingroup associate
 */

#ifndef LSST_AP_DATA_TRAITS_H
#define LSST_AP_DATA_TRAITS_H

#include <lsst/afw/detection/Source.h>

#include "Common.h"
#include "Object.h"


namespace lsst {
namespace ap {

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

template <> struct LSST_AP_LOCAL DataTraits<SimpleObject> {
    static uint32_t const ENTRIES_PER_BLOCK_LOG2 = 12;
    static uint32_t const MAX_BLOCKS_PER_CHUNK   = 128;
    static uint32_t const MAX_CHUNKS_PER_FOV     = 128;
    static uint32_t const NUM_BLOCKS             = 1024;
};


}} // end of namespace lsst::ap

#endif // LSST_AP_DATA_TRAITS_H
