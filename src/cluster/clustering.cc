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
#include "lsst/ap/cluster/detail/Metrics.h"
#include "lsst/ap/cluster/detail/Optics.cc"

using lsst::pex::exceptions::InvalidParameterError;
using lsst::pex::exceptions::NotFoundError;

using lsst::pex::logging::Log;

using lsst::afw::coord::IcrsCoord;

using lsst::afw::geom::AffineTransform;
using lsst::afw::geom::Angle;
using lsst::afw::geom::radians;

using lsst::afw::image::indexToPosition;
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
                 std::string const & name,
                 std::string const & unit)
    {
        std::string doc = "inverse variance weighted mean of " + filter +
                          "-filter source " + name + " (" +
                          proto.find<Flux::MeasTag>(name).field.getDoc() + ")";
        addFluxFields(schema, filter, name, doc, unit);
    }

    void addShape(Schema & schema,
                  Schema const & proto,
                  std::string const & filter,
                  std::string const & name)
    {
        std::string doc = "inverse variance weighted mean of " + filter +
                          "-filter source " + name + " (" +
                          proto.find<Shape::MeasTag>(name).field.getDoc() + ")";
        addShapeFields(schema, filter, name, doc);
    }

    bool hasFluxField(Schema const & schema, std::string const & name) {
        try {
            schema.find<Flux::MeasTag>(name);
        } catch (NotFoundError &) {
            return false;
        }
        return true;
    }

    bool hasShapeField(Schema const & schema, std::string const & name) {
        try {
            schema.find<Shape::MeasTag>(name);
        } catch (NotFoundError &) {
            return false;
        }
        return true;
    }

    struct FilterNameCmp {
        bool operator()(std::string const &, std::string const &) const; 
    };

    bool FilterNameCmp::operator()(std::string const & a, std::string const & b) const {
        return Filter(a, false).getId() < Filter(b, false).getId();
    }

}   // namespace <anonymous>


std::pair<boost::shared_ptr<SourceTable>, SchemaMapper> const makeOutputSourceTable(
    SourceTable const & prototype,
    SourceProcessingControl const & control)
{
    if (!prototype.getCoordKey().isValid() ||
        !prototype.getCentroidKey().isValid()) {
        throw LSST_EXCEPT(InvalidParameterError, "Prototypical "
            "SourceCatalog must have valid Coord and Centroid fields");
    }
    SchemaMapper mapper(prototype.getSchema());
    // add all fields in the input schema to the output schema.
    mapper.addMappingsWhere(True());
    // add coordinate error
    mapper.addOutputField(Field<Covariance<lsst::afw::table::Point<float> > >(
        "coord.err",
        "covariance matrix for source sky-coordinates",
        "rad^2"));
    if (!control.exposurePrefix.empty()) {
        // add exposure related fields
        mapper.addOutputField(Field<int64_t>(
            control.exposurePrefix + ".id",
            "ID of exposure source was measured on"));
        if (!control.multiBand) {
            mapper.addOutputField(Field<int>(
                control.exposurePrefix + ".filter.id",
                "ID of filter for exposure"));
        }
        if (!control.coadd) {
            mapper.addOutputField(Field<float>(
                control.exposurePrefix + ".time",
                "exposure time",
                "s"));
            mapper.addOutputField(Field<double>(
                control.exposurePrefix + ".time.mid",
                "middle of exposure time",
                "mjd"));
        }
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
    table->definePsfFlux(prototype.getPsfFluxDefinition());
    table->defineModelFlux(prototype.getModelFluxDefinition());
    table->defineApFlux(prototype.getApFluxDefinition());
    table->defineInstFlux(prototype.getInstFluxDefinition());
    table->defineCentroid(prototype.getCentroidDefinition());
    table->defineShape(prototype.getShapeDefinition());

    return std::make_pair(table, mapper);
}


boost::shared_ptr<SourceClusterTable> const makeSourceClusterTable(
    SourceTable const & prototype,
    boost::shared_ptr<IdFactory> const & idFactory,
    SourceProcessingControl const & control)
{
    typedef std::vector<std::string>::const_iterator Iter;

    Schema schema = SourceClusterTable::makeMinimalSchema();
    schema.setVersion(0);
    // Add basic fields
    schema.addField<Covariance<lsst::afw::table::Point<float> > >(
        "coord.err",
        "sample covariance matrix for coord field (unweighted mean sky-coordinates)",
        "rad^2");
    schema.addField<int>(
        "obs.count",
        "number of sources in cluster");
    schema.addField<Flag>(
        "flag.noise",
        "set if cluster was created from a single noise source");
    if (prototype.getCentroidKey().isValid() &&
        prototype.getCentroidErrKey().isValid()) {
        schema.addField<lsst::afw::table::Coord>(
            "coord.weightedmean",
            "inverse variance weighted mean sky-coordinates (ICRS)",
            "rad");
        schema.addField<Covariance<lsst::afw::table::Point<float> > >(
            "coord.weightedmean.err",
            "covariance matrix for coord.weightedmean field",
            "rad^2");
        schema.addField<int>(
            "coord.weightedmean.count",
            "number of samples used to compute coord.weightedmean");
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
    std::vector<std::string> filters = Filter::getNames();
    // Sort filter names by ID so that fields are added in a consistent order
    std::sort(filters.begin(), filters.end(), FilterNameCmp());
    // Add filter specific fields
    for (Iter filt = filters.begin(), eFilt = filters.end(); filt != eFilt; ++filt) {
        schema.addField<int>(
            *filt + ".obs.count",
            "number of " + *filt + "-filter sources in cluster");
        if (!control.exposurePrefix.empty()) {
            schema.addField<double>(
                *filt + ".obs.time.min",
                "earliest observation time of " + *filt + "-filter sources in cluster (TAI)",
                "mjd");
            schema.addField<double>(
                *filt + ".obs.time.max",
                "latest observation time of " + *filt + "-filter sources in cluster (TAI)",
                "mjd");
        }
        for (Iter flux = control.fluxFields.begin(), eFlux = control.fluxFields.end();
             flux != eFlux; ++flux) {
            if (hasFluxField(prototype.getSchema(), *flux)) {
                 addFlux(schema, prototype.getSchema(), *filt, *flux, control.fluxUnit);
            }
        }
        for (Iter shape = control.shapeFields.begin(), eShape = control.shapeFields.end();
             shape != eShape; ++shape) {
            if (hasShapeField(prototype.getSchema(), *shape)) {
                 addShape(schema, prototype.getSchema(), *filt, *shape);
            }
        }
    }

    // Create table
    boost::shared_ptr<SourceClusterTable> table =
        SourceClusterTable::make(schema, idFactory);

    // setup slot mappings
    table->defineCoordErr("coord.err");
    table->defineNumSources("obs.count");
    if (prototype.getCentroidKey().isValid() &&
        prototype.getCentroidErrKey().isValid()) {
        table->defineWeightedMeanCoord("coord.weightedmean");
        table->defineWeightedMeanCoordErr("coord.weightedmean.err");
        table->defineWeightedMeanCoordCount("coord.weightedmean.count");
    }
    if (!control.exposurePrefix.empty()) {
        table->defineTimeMin("obs.time.min");
        table->defineTimeMean("obs.time.mean");
        table->defineTimeMax("obs.time.max");
    }
    for (Iter filt = filters.begin(), eFilt = filters.end();
         filt != eFilt; ++filt) {
        table->defineNumSources(*filt, "obs.count");
        if (!control.exposurePrefix.empty()) {
            table->defineTimeMin(*filt, "obs.time.min");
            table->defineTimeMax(*filt, "obs.time.max");
        }
        Iter flux = control.fluxFields.begin(), eFlux = control.fluxFields.end();
        if (prototype.getPsfFluxKey().isValid() &&
            prototype.getPsfFluxErrKey().isValid()) {
            std::string def = prototype.getPsfFluxDefinition();
            if (std::find(flux, eFlux, def) != eFlux) {
                table->definePsfFlux(*filt, def);
            }
        }
        if (prototype.getModelFluxKey().isValid() &&
            prototype.getModelFluxErrKey().isValid()) {
            std::string def = prototype.getModelFluxDefinition();
            if (std::find(flux, eFlux, def) != eFlux) {
                table->defineModelFlux(*filt, def);
            }
        }
        if (prototype.getApFluxKey().isValid() &&
            prototype.getApFluxErrKey().isValid()) {
            std::string def = prototype.getApFluxDefinition();
            if (std::find(flux, eFlux, def) != eFlux) {
                table->defineApFlux(*filt, def);
            }
        }
        if (prototype.getInstFluxKey().isValid() &&
            prototype.getInstFluxErrKey().isValid()) {
            std::string def = prototype.getInstFluxDefinition();
            if (std::find(flux, eFlux, def) != eFlux) {
                table->defineInstFlux(*filt, def);
            }
        }
        if (prototype.getShapeKey().isValid() &&
            prototype.getShapeErrKey().isValid()) {
            std::string def = prototype.getShapeDefinition();
            Iter shape = control.shapeFields.begin(), eShape = control.shapeFields.end();
            if (std::find(shape, eShape, def) != eShape) {
                table->defineShape(*filt, def);
            }
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
    lsst::ap::utils::PT1SkyTile const * skyTile,
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
        throw LSST_EXCEPT(InvalidParameterError, "output "
            "SourceCatalog schema mismatch");
    }
    if (!sources.getSchema().contains(expSources.getSchema())) {
        throw LSST_EXCEPT(InvalidParameterError, "output SourceCatalog "
            "schema does not contain input SourceCatalog schema");
    }
    if (mapper.getInputSchema() != expSources.getSchema()) {
        throw LSST_EXCEPT(InvalidParameterError, "SchemaMapper and "
            "input source catalog disagree on input schema.");
    }
    if (mapper.getOutputSchema() != sources.getSchema()) {
        throw LSST_EXCEPT(InvalidParameterError, "SchemaMapper and "
            "output source catalogs disagree on output schema.");
    }
    if (!compareSlots(*expSources.getTable(), *sources.getTable()) ||
        !compareSlots(*sources.getTable(), *badSources.getTable()) ||
        !compareSlots(*badSources.getTable(), *invalidSources.getTable())) {
        throw LSST_EXCEPT(InvalidParameterError, "slot mappings "
            "for input and output SourceCatalogs differ");
    }

    // Construct bad source flag key vector
    std::vector<Key<Flag> > badKeys;
    for (StringIter i = control.badFlagFields.begin(),
                    e = control.badFlagFields.end(); i != e; ++i) {
        badKeys.push_back(expSources.getSchema().find<Flag>(*i).key);
    }
    // Extract exposure X/Y bounds relative to parent exposure (if there is one).
    // This is because centroids of coadd-sources are in the pixel coordinate system
    // of the tract exposure containing the patch they were measured on, NOT in
    // the pixel coordinate system of the patch.
    double const xMin = indexToPosition(expInfo.getX0()) - 0.5;
    double const yMin = indexToPosition(expInfo.getY0()) - 0.5;
    double const xMax = indexToPosition(expInfo.getX0() + expInfo.getWidth() - 1) + 0.5;
    double const yMax = indexToPosition(expInfo.getY0() + expInfo.getHeight() - 1) + 0.5;
    
    // Set up keys to additional fields
    Key<int64_t> expIdKey;
    Key<int> expFilterIdKey;
    Key<float> expTimeKey;
    Key<double> expTimeMidKey;
    afw::table::CovarianceMatrixKey<float,2> centroidErrKey =
        sources.getTable()->getCentroidErrKey();
    Key<Covariance<lsst::afw::table::Point<float> > > coordErrKey;
    Key<lsst::afw::coord::Coord> clusterCoordKey;
    try {
        coordErrKey = sources.getSchema()["coord.err"];
    } catch (NotFoundError &) {
        // no easy way to ask whether a field exists by name?
    }
    if (!control.exposurePrefix.empty()) {
        Schema schema = mapper.getOutputSchema();
        expIdKey = schema[control.exposurePrefix + ".id"];
        if (!control.multiBand) {
            expFilterIdKey = schema[control.exposurePrefix + ".filter.id"];
        }
        if (!control.coadd) {
            expTimeKey = schema[control.exposurePrefix + ".time"];
            expTimeMidKey = schema[control.exposurePrefix + ".time.mid"];
        }
    }
    if (!control.clusterPrefix.empty()) {
        Schema schema = mapper.getOutputSchema();
        clusterCoordKey = schema[control.clusterPrefix + ".coord"];
    }
    size_t ngood = 0, nbad = 0, ninvalid = 0, noutside = 0;
    // Loop over input sources
    for (SourceIter s = expSources.begin(), es = expSources.end(); s != es; ++s) {
        // Check source validity
        bool invalid = false;
        double x = s->getX();
        double y = s->getY();
        if (lsst::utils::isnan(x) || lsst::utils::isnan(y)) {
            log.format(Log::WARN, "Centroid of source %lld contains NaNs",
                       static_cast<long long>(s->getId()));
            invalid = true;
        } else if (x < xMin || x > xMax || y < yMin || y > yMax) {
            log.format(Log::WARN, "Centroid of source %lld does not lie on exposure",
                       static_cast<long long>(s->getId()));
            invalid = true;
        } else if (lsst::utils::isnan(s->getRa().asRadians()) ||
                   lsst::utils::isnan(s->getDec().asRadians())) {
            log.format(Log::WARN, "Sky-coord of source %lld contains NaNs",
                       static_cast<long long>(s->getId()));
            invalid = true;
        }
        if (!invalid && skyTile && !skyTile->contains(s->getCoord())) {
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
        // Source cluster coords default to source coordinates
        if (clusterCoordKey.isValid()) {
            os->set(clusterCoordKey, s->getCoord());
        }
        // Compute sky-coordinate errors
        if (!invalid && coordErrKey.isValid() && centroidErrKey.isValid()) {
            Eigen::Matrix2d cov = s->getCentroidErr().cast<double>();
            bool computeErr = true;
            if (lsst::utils::isnan(cov(0,0)) || lsst::utils::isnan(cov(1,1))) {
                computeErr = false;
            } else if (cov(0,1) != cov(1,0)) {
                // FIXME: for now, no measurement algorithms actually
                //        compute full covariance matrixes, and the
                //        off diagonal matrix elements are always NaN.
                //        Longer term, we will need some algorithm metadata
                //        to decide whether an off-diagonal NaN means a
                //        sample should be ignored because a computation failed,
                //        or whether it should be zeroed because the algorithm
                //        never computes it.
                if (lsst::utils::isnan(cov(0,1)) && lsst::utils::isnan(cov(1,0))) {
                    cov(0,1) = 0.0; cov(1,0) = 0.0;
                } else {
                    computeErr = false;
                }
            }
            if (computeErr) {
                Eigen::Matrix2d m = expInfo.getWcs()->linearizePixelToSky(
                    s->getCentroid(), radians).getLinear().getMatrix();
                os->set(coordErrKey, (m * cov * m.transpose()).cast<float>());
            }
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

    typedef detail::Point<3, boost::shared_ptr<SourceRecord> > OpticsPoint;
    typedef detail::Optics<3, SourceRecord> Optics;

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
        throw LSST_EXCEPT(InvalidParameterError, "too many sources to cluster");
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
        optics.run(sources.getTable(), clusters, detail::SquaredEuclidianDistanceOverSphere());
    }
    return clusters;
}


// -- Cluster attributes --------

void setClusterFields(
    SourceCatalog & sources,
    SourceClusterRecord const & record,
    SourceProcessingControl const & control)
{
    typedef SourceCatalog::iterator Iter;

    if (control.clusterPrefix.empty()) {
        return;
    }
    Key<int64_t> idKey;
    Key<lsst::afw::coord::Coord> coordKey;
    try {
        idKey = sources.getSchema().find<int64_t>(control.clusterPrefix + ".id").key;
        coordKey = sources.getSchema().find<lsst::afw::coord::Coord>(
            control.clusterPrefix + ".coord").key;
    } catch (NotFoundError &) {
        // could not find fields with the expected name
        return;
    }
    int64_t const id = record.getId();
    IcrsCoord const coord = record.getCoord();
    for (Iter i = sources.begin(), e = sources.end(); i != e; ++i) {
        i->set(idKey, id);
        i->set(coordKey, coord);
    }
}


}}} // namespace lsst::ap::cluster

