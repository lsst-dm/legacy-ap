// -*- lsst-c++ -*-
//
//##====----------------                                ----------------====##/
//
//! \file   ChunkToNameMappings.cc
//
//##====----------------                                ----------------====##/

#include <stdexcept>
#include <sstream>

#include <lsst/ap/ChunkToNameMappings.h>


namespace lsst {
namespace ap {


// -- ChunkToNameMapping ----------------

ChunkToNameMapping::~ChunkToNameMapping() {}


// -- ChunkToFileNameMapping ----------------

ChunkToFileNameMapping::ChunkToFileNameMapping(std::string const & pattern) :
    _format(pattern)
{
    // OK to have a pattern that does not use all arguments
    _format.exceptions(boost::io::all_error_bits ^ boost::io::too_many_args_bit);
}


ChunkToFileNameMapping::~ChunkToFileNameMapping() {}


std::string const ChunkToFileNameMapping::getName(
    ZoneStripeChunkDecomposition const & zsc,
    int64_t                      const   chunkId,
    int                          const   version
) {
    int32_t const stripeId = ZoneStripeChunkDecomposition::chunkToStripe(chunkId);
    int32_t const sequence = ZoneStripeChunkDecomposition::chunkToSequence(chunkId);
    _format.clear();
    _format % stripeId % sequence % version;
    return _format.str();
}


}}  // end of namespace lsst::ap

