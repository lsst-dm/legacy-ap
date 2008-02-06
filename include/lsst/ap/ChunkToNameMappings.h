// -*- lsst-c++ -*-
//
//##====----------------                                ----------------====##/
//
//! \file   ChunkToNameMappings.h
//! \brief  Classes for mapping chunk ids to names (of files, database tables, etc...).
//
//##====----------------                                ----------------====##/

#ifndef LSST_AP_CHUNK_TO_NAME_MAPPINGS_H
#define LSST_AP_CHUNK_TO_NAME_MAPPINGS_H

#include <string>

#include <boost/format.hpp>

#include "SpatialUtil.h"


namespace lsst {
namespace ap {    


/*! \brief  Interface for mapping chunk ids to names. */
class LSST_AP_LOCAL ChunkToNameMapping {
    
public :

    virtual ~ChunkToNameMapping();
    
    virtual std::string const getName(
        std::string                  const & runId,
        ZoneStripeChunkDecomposition const & zsc,
        int64_t                      const   chunkId,
        int                          const   version = 0
    ) = 0;
};


/*!
    \brief  Maps a chunk id to a file name according to a \c boost::format compatible pattern.

    The pattern can use any of 4 parameters, passed to the formatter as follows:
    <ol>
    <li>An identifier for the pipeline run</li>
    <li>The stripe id of the chunk</li>
    <li>The sequence number of the chunk (within its stripe)</li>
    <li>An integer file version number</li>
    </ol>
 */
class LSST_AP_LOCAL ChunkToFileNameMapping : public ChunkToNameMapping {

public :

    ChunkToFileNameMapping(std::string const & pattern);

    virtual ~ChunkToFileNameMapping();

    virtual std::string const getName(
        std::string                  const & runId,
        ZoneStripeChunkDecomposition const & zsc,
        int64_t                      const   chunkId,
        int                          const   version = 0
    );

private :

    boost::format _format;
};


}}  // end of namespace lsst::ap

#endif // LSST_AP_CHUNK_TO_NAME_MAPPINGS_H
