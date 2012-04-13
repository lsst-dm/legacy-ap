// -*- lsst-c++ -*-

/*
 * LSST Data Management System
 * Copyright 2012 LSST Corporation.
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
  * @brief  Clustering implementation.
  * @author Serge Monkewitz
  */
#include "lsst/ap/cluster/clustering.h"

#include "lsst/utils/ieee.h"
#include "lsst/pex/logging/Log.h"
#include "lsst/afw/image/ImageUtils.h"
#include "lsst/afw/table/SchemaMapper.h"
#include "lsst/ap/cluster/optics/Metrics.h"
#include "lsst/ap/cluster/optics/Optics.cc"

using lsst::pex::exceptions::InvalidParameterException;
using lsst::pex::logging::Log;
using lsst::afw::geom::AffineTransform;
using lsst::afw::geom::Angle;
using lsst::afw::geom::radians;
using lsst::afw::image::Filter;
using lsst::afw::table::Covariance;
using lsst::afw::table::Field;
using lsst::afw::table::Flag;
using lsst::afw::table::Flux;
using lsst::afw::table::IdFactory;
using lsst::afw::table::Key;
using lsst::afw::table::SourceCatalog;
using lsst::afw::table::SourceRecord;
using lsst::afw::table::SourceTable;
using lsst::afw::table::Schema;
using lsst::afw::table::SchemaMapper;
using lsst::afw::table::Shape;
using lsst::ap::match::ExposureInfo;
using lsst::ap::utils::PT1SkyTile;


namespace lsst { namespace ap { namespace cluster {

// -- Schema / table/ catalog setup --------

namespace {

    // Unary predicate that always returns true.
    struct True {
        template <typename T> bool operator()(T dummy) const {
            return true;
        }
    };

    void addFlux(Schema & schema,
                 Schema const & proto,
                 std::string const & filter,
                 Key<Flux::MeasTag> const & flux,
                 Key<Flux::ErrTag> const & fluxErr)
    {
        if (flux.isValid()) {
            std::string const name = proto.find(flux).field.getName();
            std::string doc(fluxErr.isValid() ? "inverse variance " : "un");
            doc += "weighted mean of " + filter + " filter source " + name;
            doc += " (" + proto.find(flux).field.getDoc() + ")";
            addFluxFields(schema, filter, name, doc);
        }
    }

    void addShape(Schema & schema,
                  Schema const & proto,
                  std::string const & filter,
                  Key<Shape::MeasTag> const & shape,
                  Key<Shape::ErrTag> const & shapeErr)
    {
        if (shape.isValid()) {
            std::string const name = proto.find(shape).field.getName();
            std::string doc(shapeErr.isValid() ? "inverse variance " : "un");
            doc += "weighted mean of " + filter + " filter source " + name;
            doc += " (" + proto.find(shape).field.getDoc() + ")";
            addShapeFields(schema, filter, name, doc);
        }
    }

}   // namespace <anonymous>


std::pair<boost::shared_ptr<SourceTable>, SchemaMapper> const makeOutputSourceTable(
    SourceTable const & prototype,
    SourceProcessingControl const & control)
{
    if (!prototype.getCoordKey().isValid() ||
        !prototype.getCentroidKey().isValid()) {
        throw LSST_EXCEPT(InvalidParameterException, "Prototypical "
            "SourceCatalog must have valid Coord and Centroid fields");
    }
    SchemaMapper mapper(prototype.getSchema());
    // add all fields in the input schema to the output schema.
    mapper.addMappingsWhere(True());
    // add coordinate error
    mapper.addOutputField(Field<Covariance<lsst::afw::table::Point<double> > >(
        "coord.err",
        "covariance matrix for source sky-coordinates",
        "rad^2"));
    if (!control.exposurePrefix.empty()) {
        // add exposure related fields
        mapper.addOutputField(Field<int64_t>(
            control.exposurePrefix + ".id",
            "ID of exposure source was measured on"));
        mapper.addOutputField(Field<int>(
            control.exposurePrefix + ".filter.id",
            "ID of filter for exposure"));
        mapper.addOutputField(Field<float>(
            control.exposurePrefix + ".time",
            "exposure time",
            "s"));
        mapper.addOutputField(Field<double>(
            control.exposurePrefix + ".time.mid",
            "middle of exposure time",
            "mjd"));
    }
    if (!control.clusterPrefix.empty()) {
        // add cluster related fields
        mapper.addOutputField(Field<int64_t>(
            control.clusterPrefix + ".id",
            "ID of cluster containing source; 0 if source is not assigned to any cluster."));
        mapper.addOutputField(Field<lsst::afw::table::Coord>(
            control.clusterPrefix + ".coord",
            "ICRS sky-coordinates of cluster containing source",
            "rad"));
    }
    boost::shared_ptr<SourceTable> table =
        SourceTable::make(mapper.getOutputSchema(), boost::shared_ptr<IdFactory>());
    // copy slot mappings from prototype
    table->definePsfFlux(prototype.getPsfFluxKey(),
                         prototype.getPsfFluxErrKey(),
                         prototype.getPsfFluxFlagKey());
    table->defineModelFlux(prototype.getModelFluxKey(),
                           prototype.getModelFluxErrKey(),
                           prototype.getModelFluxFlagKey());
    table->defineApFlux(prototype.getApFluxKey(),
                        prototype.getApFluxErrKey(),
                        prototype.getApFluxFlagKey());
    table->defineInstFlux(prototype.getInstFluxKey(),
                          prototype.getInstFluxErrKey(),
                          prototype.getInstFluxFlagKey());
    table->defineCentroid(prototype.getCentroidKey(),
                          prototype.getCentroidErrKey(),
                          prototype.getCentroidFlagKey());
    table->defineShape(prototype.getShapeKey(),
                       prototype.getShapeErrKey(),
                       prototype.getShapeFlagKey());
    return std::make_pair(table, mapper);
}


boost::shared_ptr<SourceClusterTable> const makeSourceClusterTable(
    SourceTable const & prototype,
    boost::shared_ptr<IdFactory> const & idFactory,
    SourceProcessingControl const & control)
{
    typedef std::vector<std::string>::const_iterator Iter;

    Schema schema = SourceClusterTable::makeMinimalSchema();
    // Add basic fields
    schema.addField<Covariance<lsst::afw::table::Point<double> > >(
        "coord.err",
        "sample covariance matrix for coord field (unweighted mean sky-coordinates)",
        "rad^2");
    schema.addField<int>(
        "obs.num",
        "number of sources in cluster");
    schema.addField<Flag>(
        "flag.bad",
        "set if cluster was created from a single bad source");
    schema.addField<Flag>(
        "flag.noise",
        "set if cluster was created from a single noise source");
    if (prototype.getCentroidKey().isValid() &&
        prototype.getCentroidErrKey().isValid()) {
        schema.addField<lsst::afw::table::Coord>(
            "coord2",
            "inverse variance weighted mean sky-coordinates (ICRS)",
            "rad");
        schema.addField<Covariance<lsst::afw::table::Point<double> > >(
            "coord2.err",
            "covariance matrix for coord2 field",
            "rad^2");
    }
    if (!control.exposurePrefix.empty()) {
        schema.addField<double>(
            "obs.time.min",
            "earliest observation time of sources in cluster (TAI)",
            "mjd");
        schema.addField<double>(
            "obs.time.mean",
            "mean observation time of sources in cluster (TAI)",
            "mjd");
        schema.addField<double>(
            "obs.time.max",
            "latest observation time of sources in cluster (TAI)",
            "mjd");
    }
    // Add filter specific fields
    std::vector<std::string> const filters = Filter::getNames();
    for (Iter f = filters.begin(), e = filters.end(); f != e; ++f) {
        schema.addField<int>(
            *f + ".obs.num",
            "number of " + *f + " filter sources in cluster");
        if (!control.exposurePrefix.empty()) {
            schema.addField<double>(
                *f + ".obs.time.min",
                "earliest observation time of " + *f + " filter sources in cluster (TAI)",
                "mjd");
            schema.addField<double>(
                *f + ".obs.time.max",
                "latest observation time of " + *f + " filter sources in cluster (TAI)",
                "mjd");
        }
        addFlux(schema, prototype.getSchema(), *f,
                prototype.getPsfFluxKey(), prototype.getPsfFluxErrKey());
        addFlux(schema, prototype.getSchema(), *f,
                prototype.getModelFluxKey(), prototype.getModelFluxErrKey());
        addFlux(schema, prototype.getSchema(), *f,
                prototype.getApFluxKey(), prototype.getApFluxErrKey());
        addFlux(schema, prototype.getSchema(), *f,
                prototype.getInstFluxKey(), prototype.getInstFluxErrKey());
        addShape(schema, prototype.getSchema(), *f,
                 prototype.getShapeKey(), prototype.getShapeErrKey());
    }

    // Create table
    boost::shared_ptr<SourceClusterTable> table =
        SourceClusterTable::make(schema, idFactory);

    // setup slot mappings
    table->defineCoordErr("coord.err");
    table->defineNumSources("obs.num");
    if (prototype.getCentroidKey().isValid() &&
        prototype.getCentroidErrKey().isValid()) {
        table->defineCoord2("coord2.err");
        table->defineCoord2Err("coord2.err");
    }
    if (!control.exposurePrefix.empty()) {
        table->defineTimeMin("obs.time.min");
        table->defineTimeMin("obs.time.mean");
        table->defineTimeMin("obs.time.max");
    }
    for (Iter f = filters.begin(), e = filters.end(); f != e; ++f) {
        table->defineNumSources(*f, "obs.num");
        if (!control.exposurePrefix.empty()) {
            table->defineTimeMin(*f, "obs.time.min");
            table->defineTimeMin(*f, "obs.time.max");
        }
        if (prototype.getPsfFluxKey().isValid()) {
            table->definePsfFlux(*f, prototype.getPsfFluxDefinition());
        }
        if (prototype.getModelFluxKey().isValid()) {
            table->defineModelFlux(*f, prototype.getModelFluxDefinition());
        }
        if (prototype.getApFluxKey().isValid()) {
            table->defineApFlux(*f, prototype.getApFluxDefinition());
        }
        if (prototype.getInstFluxKey().isValid()) {
            table->defineInstFlux(*f, prototype.getInstFluxDefinition());
        }
        if (prototype.getShapeKey().isValid()) {
            table->defineShape(*f, prototype.getShapeDefinition());
        }
    }
    return table;
}


// -- Source processing --------

namespace {

    // Return true if the slot mappings for the two catalogs are identical
    bool compareSlots(SourceTable const & a, SourceTable const & b) {
        return (
            a.getPsfFluxKey()       == b.getPsfFluxKey()       &&
            a.getPsfFluxErrKey()    == b.getPsfFluxErrKey()    &&
            a.getPsfFluxFlagKey()   == b.getPsfFluxFlagKey()   &&
            a.getModelFluxKey()     == b.getModelFluxKey()     &&
            a.getModelFluxErrKey()  == b.getModelFluxErrKey()  &&
            a.getModelFluxFlagKey() == b.getModelFluxFlagKey() &&
            a.getApFluxKey()        == b.getApFluxKey()        &&
            a.getApFluxErrKey()     == b.getApFluxErrKey()     &&
            a.getApFluxFlagKey()    == b.getApFluxFlagKey()    &&
            a.getInstFluxKey()      == b.getInstFluxKey()      &&
            a.getInstFluxErrKey()   == b.getInstFluxErrKey()   &&
            a.getInstFluxFlagKey()  == b.getInstFluxFlagKey()  &&
            a.getCentroidKey()      == b.getCentroidKey()      &&
            a.getCentroidErrKey()   == b.getCentroidErrKey()   &&
            a.getCentroidFlagKey()  == b.getCentroidFlagKey()  &&
            a.getShapeKey()         == b.getShapeKey()         &&
            a.getShapeErrKey()      == b.getShapeErrKey()      &&
            a.getShapeFlagKey()     == b.getShapeFlagKey()
        );
    }

}   // namespace <anonymous>


void processSources(
    SourceCatalog const & expSources,
    ExposureInfo const & expInfo,
    PT1SkyTile const & skyTile,
    SourceProcessingControl const & control,
    SchemaMapper const & mapper,
    SourceCatalog & sources,
    SourceCatalog & badSources,
    SourceCatalog & invalidSources)
{
    typedef std::vector<std::string>::const_iterator StringIter;
    typedef std::vector<Key<Flag> >::const_iterator FlagKeyIter;
    typedef SourceCatalog::const_iterator SourceIter;

    Log log(Log::getDefaultLog(), "lsst.ap.cluster");
    log.format(Log::INFO, "Processing sources from exposure %lld",
               static_cast<long long>(expInfo.getId()));

    // Validate that schemas and slots match as expected
    if (sources.getSchema() != badSources.getSchema() ||
        sources.getSchema() != invalidSources.getSchema()) {
        throw LSST_EXCEPT(InvalidParameterException, "output "
            "SourceCatalog schema mismatch");
    }
    if (!sources.getSchema().contains(expSources.getSchema())) {
        throw LSST_EXCEPT(InvalidParameterException, "output SourceCatalog "
            "schema does not contain input SourceCatalog schema");
    }
    if (mapper.getInputSchema() != expSources.getSchema()) {
        throw LSST_EXCEPT(InvalidParameterException, "SchemaMapper and "
            "input source catalog disagree on input schema.");
    }
    if (mapper.getOutputSchema() != sources.getSchema()) {
        throw LSST_EXCEPT(InvalidParameterException, "SchemaMapper and "
            "output source catalogs disagree on output schema.");
    }
    if (!compareSlots(*expSources.getTable(), *sources.getTable()) ||
        !compareSlots(*sources.getTable(), *badSources.getTable()) ||
        !compareSlots(*badSources.getTable(), *invalidSources.getTable())) {
        throw LSST_EXCEPT(InvalidParameterException, "slot mappings "
            "for input and output SourceCatalogs differ");
    }

    // Construct bad source flag key vector
    std::vector<Key<Flag> > badKeys;
    for (StringIter i = control.badFlagFields.begin(),
                    e = control.badFlagFields.end(); i != e; ++i) {
        badKeys.push_back(expSources.getSchema().find<Flag>(*i).key);
    }
    // extract exposure X/Y bounds
    double const xyMin = lsst::afw::image::indexToPosition(0) - 0.5;
    double const xMax  = lsst::afw::image::indexToPosition(expInfo.getWidth() - 1) + 0.5;
    double const yMax  = lsst::afw::image::indexToPosition(expInfo.getHeight() - 1) + 0.5;
    
    // Set up keys to additional fields
    Key<int64_t> expIdKey;
    Key<int> expFilterIdKey;
    Key<float> expTimeKey;
    Key<double> expTimeMidKey;
    Key<Covariance<lsst::afw::table::Point<double> > > centroidErrKey =
        sources.getTable()->getCentroidErrKey();
    Key<Covariance<lsst::afw::table::Point<double> > > coordErrKey;
    try {
        coordErrKey = sources.getSchema()["coord.err"];
    } catch (...) {
        // no easy way to ask whether a field exists by name?
    }
    if (!control.exposurePrefix.empty()) {
        Schema schema = mapper.getOutputSchema();
        expIdKey = schema[control.exposurePrefix + ".id"];
        expFilterIdKey = schema[control.exposurePrefix + ".filter.id"];
        expTimeKey = schema[control.exposurePrefix + ".time"];
        expTimeMidKey = schema[control.exposurePrefix + ".time.mid"];
    }
    size_t ngood = 0, nbad = 0, ninvalid = 0, noutside = 0;
    // Loop over input sources
    for (SourceIter s = expSources.begin(), es = expSources.end(); s != es; ++s) {
        // Check source validity
        bool invalid = false;
        double x = s->getX(), y = s->getY();
        if (lsst::utils::isnan(x) || lsst::utils::isnan(y)) {
            log.format(Log::WARN, "Centroid of source %lld contains NaNs",
                       static_cast<long long>(s->getId()));
            invalid = true;
        } else if (x < xyMin || x > xMax || y < xyMin || y > yMax) {
            log.format(Log::WARN, "Centroid of source %lld does not lie on exposure",
                       static_cast<long long>(s->getId()));
            invalid = true;
        } else if (lsst::utils::isnan(s->getRa().asRadians()) ||
                   lsst::utils::isnan(s->getDec().asRadians())) {
            log.format(Log::WARN, "Sky-coord of source %lld contains NaNs",
                       static_cast<long long>(s->getId()));
            invalid = true;
        }
        if (!invalid && !skyTile.contains(s->getCoord())) {
            // skip sources outside the sky-tilea
            noutside += 1;
            continue;
        }
        bool bad = false;
        // Check whether source is flagged as bad
        for (FlagKeyIter f = badKeys.begin(), ef = badKeys.end(); f != ef; ++f) {
            if (s->get(*f)) {
                bad = true;
                break;
            }
        }
        // Copy input source to appropriate output catalog
        boost::shared_ptr<SourceRecord> os;
        if (invalid) {
            ninvalid += 1;
            os = invalidSources.addNew();
        } else if (bad) {
            nbad += 1;
            os = badSources.addNew();
        } else {
            ngood += 1;
            os = sources.addNew();
        }
        os->assign(*s, mapper);

        // Add exposure parameters
        if (expIdKey.isValid()) {
            os->set(expIdKey, expInfo.getId());
        }
        if (expFilterIdKey.isValid()) {
            os->set(expFilterIdKey, expInfo.getFilter().getId());
        }
        if (expTimeKey.isValid()) {
            os->set(expTimeKey, expInfo.getExposureTime());
        }
        if (expTimeMidKey.isValid()) {
            os->set(expTimeMidKey, expInfo.getEpoch());
        }

        // Compute sky-coordinate errors
        if (!invalid && coordErrKey.isValid() && centroidErrKey.isValid()) {
            Eigen::Matrix2d m = expInfo.getWcs()->linearizePixelToSky(
                s->getCentroid(), radians).getLinear().getMatrix();
            Eigen::Matrix2d cov = s->getCentroidErr();
            os->set(coordErrKey, m * cov * m.transpose());
        }
    }
    log.format(Log::INFO, "processed %lld sources (invalid: %lld, "
               "outside sky-tile: %lld, bad: %lld, good: %lld)",
               static_cast<long long>(expSources.size()),
               static_cast<long long>(ninvalid),
               static_cast<long long>(noutside),
               static_cast<long long>(nbad),
               static_cast<long long>(ngood));
}


// -- Spatial clustering --------

namespace {

    typedef optics::Point<3, boost::shared_ptr<SourceRecord> > OpticsPoint;
    typedef optics::Optics<3, SourceRecord> Optics;

    /// @internal  Maximum number of sources that can be clustered at once.
    unsigned int const MAX_SOURCES =
        static_cast<unsigned int>(std::numeric_limits<int>::max());

} // namespace <anonymous>


std::vector<SourceCatalog> const cluster(
    SourceCatalog const & sources,
    ClusteringControl const & control)
{
    typedef SourceCatalog::const_iterator Iter;

    if (sources.size() > MAX_SOURCES) {
        throw LSST_EXCEPT(InvalidParameterException, "too many sources to cluster");
    }
    control.validate();
    boost::scoped_array<OpticsPoint> entries(new OpticsPoint[sources.size()]);
    std::vector<SourceCatalog> clusters;
    int i = 0;
    for (Iter s = sources.begin(), e = sources.end(); s != e; ++s, ++i) {
        entries[i].coords = s->getCoord().getVector().asEigen();
        entries[i].data = s;
    }
    if (i > 0) {
        // Convert epsilon and leafExtentThreshold to radians, and account
        // for the fact that our metric is the squared euclidian distance,
        // not angular separation.
        double eps = std::sin(0.5 * control.getEpsilon().asRadians());
        eps = 4.0 * eps * eps;
        double let = control.getLeafExtentThreshold().asRadians();
        if (let > 0.0) {
            let = std::sin(0.5 * let);
            let = 4.0 * let * let;
        }
        Optics optics(entries.get(), i, control.minNeighbors, eps, let, control.pointsPerLeaf);
        optics.run(sources.getTable(), clusters, optics::SquaredEuclidianDistanceOverSphere());
    }
    return clusters;
}

}}} // namespace lsst::ap::cluster

