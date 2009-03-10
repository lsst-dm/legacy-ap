// -*- lsst-c++ -*-

/**
 * @file
 * @brief   Implementation of all chunk to name mappings.
 *
 * @ingroup associate
 */

#include <stdexcept>
#include <sstream>

#include "lsst/ap/ChunkToNameMappings.h"


// -- ChunkToNameMapping ----------------

lsst::ap::ChunkToNameMapping::~ChunkToNameMapping() {}


// -- ChunkToFileNameMapping ----------------

lsst::ap::ChunkToFileNameMapping::ChunkToFileNameMapping(std::string const & pattern) :
    _format(pattern)
{
    // OK to have a pattern that does not use all arguments
    _format.exceptions(boost::io::all_error_bits ^ boost::io::too_many_args_bit);
}


lsst::ap::ChunkToFileNameMapping::~ChunkToFileNameMapping() {}


std::string const lsst::ap::ChunkToFileNameMapping::getName(
    std::string                  const & runId,
    ZoneStripeChunkDecomposition const & zsc,
    boost::int64_t               const   chunkId,
    int                          const   version
) {
    int const stripeId = ZoneStripeChunkDecomposition::chunkToStripe(chunkId);
    int const sequence = ZoneStripeChunkDecomposition::chunkToSequence(chunkId);
    _format.clear();
    _format % runId % stripeId % sequence % version;
    return _format.str();
}

