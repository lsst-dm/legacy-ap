#
# LSST Data Management System
# Copyright 2012 LSST Corporation.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#
import sys
import traceback

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.afw.image as afwImage
import lsst.afw.table as afwTable
import lsst.meas.algorithms as measAlgorithms
import lsst.skypix as skypix
import lsst.ap.match as apMatch
import lsst.ap.utils as apUtils
import lsst.ap.cluster as apCluster

from .sourceAssocArgumentParser import SourceAssocArgumentParser

__all__ = ["SourceAssocConfig", "SourceAssocTask"]


class SourceAssocConfig(pexConfig.Config):
    """Configuration parameters for SourceAssocTask.
    """
    inputLevel = pexConfig.Field(
        dtype=str, default="sensor",
        doc="""
            Level of input datasets identified by the
            inputSourceDataset and inputCalexpMetadataDataset
            configuration parameters.
            """)
    inputSourceDataset = pexConfig.Field(
        dtype=str, default="src",
        doc="Name of the butler dataset for input sources.")
    inputCalexpMetadataDataset = pexConfig.Field(
        dtype=str, default="calexp_md",
        doc="""
            Name of the butler dataset for calibrated exposure metadata.

            Note that this dataset must yield metadata for the calibrated
            exposures that the sources from inputSourceDataset were detected
            and measured on; otherwise, the SourceAssoc pipeline behavior 
            is undefined.
            """)

    sourceProcessing = pexConfig.ConfigField(
        dtype=apCluster.SourceProcessingConfig,
        doc="""
            Source processing parameters.

            To see their descriptions:

            >>> from lsst.ap.cluster import SourceProcessingConfig
            >>> help(SourceProcessingConfig)
            """)

    clustering = pexConfig.ConfigField(
        dtype=apCluster.ClusteringConfig,
        doc="""
            Source clustering parameters.

            To see their descriptions:

            >>> from lsst.ap.cluster import ClusteringConfig
            >>> help(ClusteringConfig)
            """)
    doCluster = pexConfig.Field(
        dtype=bool, default=True,
        doc="""
            If set to True, then "good" sources are clustered with the 
            OPTICS algorithm - this is an attempt to group sources from
            individual exposures which correspond to measurements of the
            same astronomical object.

            If set to False, then running SourceAssocTask reduces to simply
            processing sources - this involves adding exposure ID, filter,
            and middle-of-exposure time fields to each source, as well as
            computation of sky-coordinate errors from centroid errors. In
            other words, sources are prepared for database ingest into the
            LSST Source database table, but are not associated with one
            another.

            Note that a "good" source is one with valid sky-coordinates and
            which has not been identified as "bad" by one of the flag fields
            listed in the sourceProcessing.badFlagFields configuration
            parameter.
            """)
    doDiscardNoiseClusters = pexConfig.Field(
        dtype=bool, default=True,
        doc="Discard single source clusters?")
    doWriteClusters = pexConfig.Field(
        dtype=bool, default=True,
        doc="""
            Write source clusters to persistent storage via the butler?

            Source clusters are stored as lsst.ap.cluster.SourceClusterCatalog
            instances; one such catalog is stored per sky-tile. The
            corresponding butler dataset name is "object", so they can
            be retrieved using e.g.:

            >>> clusters = butler.get("object", skyTile=12345)

            Note that if no clusters were generated for a sky-tile, say
            because the doCluster configuration parameter was set to False,
            then nothing (not even an empty catalog!) is written to persistent
            storage. In other words:

            >>> butler.datasetExists("object", skyTile=12345)

            can return False, even after SourceAssocTask has been run on
            sky-tile 12345.
            """)
    algorithmFlags = pexConfig.DictField(
        keytype=str, itemtype=str,
        doc="""
            A dictionary mapping from algorithm names to strings containing
            comma separated lists of flag field names. If any flag is set for
            a source, then that source is ignored when computing the
            measurement mean of the corresponding algorithm.
            """)

    doWriteSources = pexConfig.Field(
        dtype=bool, default=True,
        doc="""
            Write processed "good" sources to persistent storage via the butler?

            A "good" source is one with valid coordinates and which has not
            been identified as "bad" by one of the flag fields listed in the
            sourceProcessing.badFlagFields configuration parameter. Sources are
            stored as lsst.afw.table.SourceCatalog instances; one such catalog
            is stored per sky-tile. The corresponding butler dataset name is
            "source", so they can be retrieved using e.g.:

            >>> sources = butler.get("source", skyTile=12345)

            Note that if no "good" sources were identified for a sky-tile, then
            nothing (not even an empty catalog!) is written to persistent
            storage. In other words:

            >>> butler.datasetExists("source", skyTile=12345)

            can return False. In this case, the clustering algorithm had no
            input, so no "object" dataset will have been written for the
            sky-tile either.
            """)
    doWriteBadSources = pexConfig.Field(
        dtype=bool, default=True,
        doc="""
            Write processed "bad" sources to persistent storage via the butler?

            A "bad" source is one with valid coordinates and for which at least
            one of the flag fields listed in the sourceProcessing.badFlagFields
            configuration parameter has been set. Bad sources are stored as
            lsst.afw.table.SourceCatalog instances; one such catalog is
            stored per sky-tile. The corresponding butler dataset name is
            "badSource", so they can be retrieved using e.g.:

            >>> badSources = butler.get("badSource", skyTile=12345)

            If no "bad" sources were identified in the sky-tile, no
            dataset is written (not even an empty catalog).
            """)
    doWriteInvalidSources = pexConfig.Field(
        dtype=bool, default=True,
        doc="""
            Write "invalid" sources to persistent storage via the butler?

            An "invalid" source is one with invalid coordinates/centroid,
            e.g. because the centroid algorithm failed. Invalid sources are
            stored as lsst.afw.table.SourceCatalog instances; one such catalog
            is stored per sky-tile. The corresponding butler dataset name is
            "invalidSource", so they can be retrieved using e,g,:

            >>> invalidSources =  butler.get("invalidSource", skyTile=12345)

            As for the other datasets, if no "invalid" sources were
            identified for the sky-tile, no dataset is written.
            """)

    sourceHistogramResolution = pexConfig.RangeField(
        dtype=int, default=2000, min=1,
        doc="X and Y resolution of 2D source position histograms.")
    doMakeSourceHistogram = pexConfig.Field(
        dtype=bool, default=True,
        doc="""
            Make 2D histogram of "good" source positions? If set, a square
            image covering the sky-tile being processed and with resolution
            equal to sourceHistogramResolution will be created. The value of
            each pixel in this image will be the number of "good" sources
            inside that pixel.
            """)
    doMakeBadSourceHistogram = pexConfig.Field(
        dtype=bool, default=True,
        doc="""
            Make 2D histogram of "bad" source positions? If set, a square
            image covering the sky-tile being processed and with resolution
            equal to sourceHistogramResolution will be created. The value of
            each pixel in this image will be the number of "bad" sources
            inside that pixel.
            """)

    doWriteSourceHistogram = pexConfig.Field(
        dtype=bool, default=True,
        doc="""
            Write "good" source histogram to persistent storage via the butler?

            If True, one histogram image is written per sky-tile containing
            at least one "good" source, via the butler. The corresponding
            dataset name is "sourceHist", and the type of these histograms is
            lsst.afw.image.DecoratedImageU. They can be retrieved using e.g.:

            >>> img = butler.get("sourceHist", skyTile=12345)
            """)
    doWriteBadSourceHistogram = pexConfig.Field(
        dtype=bool, default=True,
        doc="""
            Write "bad" source histogram to persistent storage via the butler?
            If True, one histogram image is written per sky-tile containing
            at least one "bad" source, via the butler. The corresponding
            dataset name is "badSourceHist", and the type of these histograms
            is lsst.afw.image.DecoratedImageU. They can be retrieved using e.g.:

            >>> img = butler.get("badSourceHist", skyTile=12345)
            """)

    measPrefix = pexConfig.Field(
        dtype=str, optional=True, default=None,
        doc="""
            Prefix for all source measurement fields. Must match the value
            of the lsst.meas.algorithms.SourceMeasurementConfig prefix
            configuration parameter, which is typically available as
            measurement.prefix in the CCD processing task configuration.
            """)
    measSlots = pexConfig.ConfigField(
        dtype=measAlgorithms.SourceSlotConfig,
        doc="""
            Mapping from algorithms to special aliases in
            lsst.afw.table.SourceTable. Must match the value of the
            lsst.meas.algorithms.SourceMeasurementConfig slots
            configuration parameter, which is typically available as
            measurement.slots in the CCD processing task configuration.
            For the details:

            >>> from lsst.meas.algorithms import SourceSlotConfig
            >>> help(SourceSlotConfig)
            """)

    def setDefaults(self):
        self.sourceProcessing.badFlagFields = ["flags.negative",
                                               "flags.pixel.edge",
                                               "shape.sdss.flags.unweightedbad",
                                              ]
        flags = ",".join(["flags.negative",
                          "flags.badcentroid",
                          "flags.pixel.edge",
                          "flags.pixel.interpolated.center",
                          "flags.pixel.saturated.center",
                         ])
        self.algorithmFlags = {
            "flux.gaussian": "flux.gaussian.flags," + flags,
            "flux.naive": "flux.naive.flags," + flags,
            "flux.psf": "flux.psf.flags," + flags,
            "flux.sinc": "flux.sinc.flags," + flags,
            "flux.kron": "flux.kron.flags,flux.kron.flags.aperture," + flags,
            "multishapelet.exp.flux":   "multishapelet.exp.flux.flags," + flags,
            "multishapelet.dev.flux":   "multishapelet.dev.flux.flags," + flags,
            "multishapelet.combo.flux": "multishapelet.combo.flux.flags," + flags,
            "shape.sdss": "shape.sdss.flags.unweightedbad," + flags,
        }


def _flagKeys(schema, config, alg):
    """Create an lsst.afw.table.FlagKeyVector identifying sources to
       ignore when computing measurement means for the given algorithm.
    """
    vec = afwTable.FlagKeyVector()
    if alg in config.algorithmFlags:
        flags = config.algorithmFlags[alg]
        for f in flags.split(","):
            f = f.strip()
            if len(f) == 0:
                continue
            si = schema.find(f)
            if si.field.getTypeString() != "Flag":
                raise TypeError(f + " field is not a Flag field")
            vec.append(si.key) 
    return vec


class SourceAssocTask(pipeBase.CmdLineTask):
    """Cluster the sources inside a sky-tile using the OPTICS algorithm - 
       this is an attempt to group sources from individual exposures that
       correspond to measurements of the same astronomical object. For each
       cluster, means of individual source measurements (currently including
       fluxes, shapes, and sky-coordinates) are computed; these cluster
       attributes are suitable for ingest into the LSST Object database table.

       For details on the clustering algorithm used, see:

           "OPTICS: Ordering Points To Identify the Clustering Structure".
           Mihael Ankerst, Markus M. Breunig, Hans-Peter Kriegel, Jorg Sander (1999).
           ACM SIGMOD international conference on Management of data.
           ACM Press. pp. 49-60.

           http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.129.6542
           http://en.wikipedia.org/wiki/OPTICS_algorithm

       The parameters of this algorithm can be adjusted via the 'clustering'
       attribute of the lsst.ap.tasks.SourceAssocConfig instance passed to
       to this task, e.g.

       >>> from lsst.ap.tasks import *
       >>> from lsst.ap.cluster import ClusteringConfig
       >>> # display documentation for top-level task configuration parameters
       ... help(SourceAssocConfig)
       >>> # display documentation for clustering parameters
       ... help(ClusteringConfig)    
       >>> # create a task configuration and fiddle with a value
       ... config = SourceAssocConfig()
       >>> config.clustering.minNeighbors = 2

       Some control over which sources participate in spatial clustering is
       provided via the 'sourceProcessing.badFlagFields' configuration 
       parameter - this is a list of flag field names. If at least one of
       these flags is set for an input source, then that source is deemed "bad"
       and does not participate in spatial clustering. Here's how one would
       set it (and discover other source processing parameters):

       >>> from lsst.ap.tasks import *
       >>> from lsst.ap.cluster import SourceProcessingConfig
       >>> # display documentation for source processing
       ... # configuration parameters
       ... help(SourceProcessingConfig)
       >>> config = SourceAssocConfig()
       >>> config.processing.badFlagFields = ["flags.evil", "flags.devilish"]

       Inputs and Assumptions:
       -----------------------

       This task retrieves "src", "calexp_md", and "processCcd_config" datasets
       from the butler (but note that the names of the input source and
       calibrated exposure metadata datasets can be changed via configuration
       parameters). It assumes that single frame measurement has been run
       (e.g. via lsst.pipe.tasks.processCcd.ProcessCcdTask, or some camera
       specific variant thereof). Note also that care is required when
       interleaving CCD processing and source association - this task does not
       support incremental updates to source clusters. In other words, one should
       ensure that all CCDs overlapping a sky-tile T have been processed before
       running this task on T. Otherwise, multiple runs on T will be required,
       where each run processes all the data for the sky-tile from scratch.

       Outputs:
       --------

       This task writes "object", "source", "badSource", "invalidSource",
       "sourceHist", and "badSourceHist" datasets via the butler. Whether any 
       of these is actually written for a sky-tile depends on configuration
       parameters and the input data.
    """
    ConfigClass = SourceAssocConfig
    _DefaultName = "sourceAssoc"

    def __init__(self, *args, **kwds):
        """Pass all parameters through to the the base class constructor."""
        pipeBase.Task.__init__(self, *args, **kwds)

    @pipeBase.timeMethod
    def cluster(self, skyTileId, butler):
        """Cluster sources falling inside the given sky-tile with the OPTICS
           algorithm.

           @param skyTileId: Integer sky-tile ID
           @param butler:    Butler responsible for retrieving calibrated
                             exposure metadata and associated sources

           @return A lsst.pipe.base.Struct with the following fields:

                   - sources:            Sources inside the sky-tile with valid
                                         positions and no "bad" flags set.
                   - badSources:         Sources inside the sky-tile with valid
                                         positions and at least one "bad" flag set.
                   - invalidSources:     Sources with invalid positions/centroids.
                   - clusters:           A list of lsst.afw.table.SourceCatalog objects,
                                         one per cluster generated. 
                   - exposures:          An lsst.ap.match.ExposureInfoMap,
                                         mapping calibrated exposure IDs to
                                         lsst.ap.match.ExposureInfo objects.
                   - sourceHistogram:    A 2D histogram of source positions.
                   - badSourceHistogram: A 2D histogram of bad source positions.

                   Note that any of these return values can legitimately be None
                   due to lack of inputs, problems reading them in, or the task
                   configuration.
        """
        # sky-tile setup
        qsp = skypix.createQuadSpherePixelization(butler.mapper.skypolicy)
        root, x, y = qsp.coords(skyTileId);
        skyTile = apUtils.PT1SkyTile(qsp.resolution, root, x, y, skyTileId)
        del root, x, y
        spControl = self.config.sourceProcessing.makeControl()
        sourceTable = None
        results = pipeBase.Struct(
            sources = None,
            badSources = None,
            invalidSources = None,
            clusters = None,
            exposures = apMatch.ExposureInfoMap(),
            sourceHistogram = None,
            badSourceHistogram = None
        )

        for dataRef in butler.subset(
            self.config.inputSourceDataset, self.config.inputLevel, skyTile=skyTileId):

            if not dataRef.datasetExists(self.config.inputSourceDataset):
                continue
            try:
                expMd = dataRef.get(self.config.inputCalexpMetadataDataset, immediate=True)
                expSources = dataRef.get(self.config.inputSourceDataset, immediate=True)
                expConfig = dataRef.get("processCcd_config", immediate=True)
                if (self.config.measPrefix != expConfig.measurement.prefix or
                    self.config.measSlots != expConfig.measurement.slots):
                    self.log.warn(str.format(
                        "skipping {} : processCcd measurement prefix and slot configuration "
                        "do not match sourceAssoc configuration", str(dataRef.dataId)))
                    continue
            except Exception:
                self.log.warn(str.format(
                    "skipping {} : failed to unpersist {}, {}, or processCcd_config dataset: {}",
                    str(dataRef.dataId), self.config.inputCalexpMetadataDataset,
                    self.config.inputSourceDataset, traceback.format_exc()))
                continue
            if sourceTable == None:
                # create output source table
                sourceTable, schemaMapper = apCluster.makeOutputSourceTable(
                    expSources.getTable(), spControl)
                # create output source catalogs
                results.sources = afwTable.SourceCatalog(sourceTable)
                results.badSources = afwTable.SourceCatalog(sourceTable)
                results.invalidSources = afwTable.SourceCatalog(sourceTable)
            # process sources: segregate into "good", "bad", and "invalid"
            # sets, discard sources outside sky-tile, and denormalize the
            # source schema.
            try:
                expInfo = apMatch.ExposureInfo(expMd)
            except:
                self.log.warn(str.format(
                    "skipping {} : failed to convert {} dataset to ExposureInfo",
                    str(dataRef.dataId), self.config.inputCalexpMetadataDataset))
                continue
            results.exposures.insert(expInfo)
            apCluster.processSources(
                expSources,
                expInfo,
                skyTile,
                spControl,
                schemaMapper,
                results.sources,
                results.badSources,
                results.invalidSources)
            # in hopes of freeing memory occupied by input sources
            del expSources, expMd

        if (sourceTable == None):
            return results # nothing to do
        # create clusters
        if self.config.doCluster and len(results.sources) > 0:
            results.clusters = apCluster.cluster(
                results.sources, self.config.clustering.makeControl())
        # create good/bad source histograms
        if self.config.doMakeSourceHistogram and len(results.sources) > 0:
            results.sourceHistogram, wcs = apUtils.createImageCoveringSkyTile(
                qsp, skyTileId, self.config.sourceHistogramResolution)
            apUtils.makeSourceHistogram(
                results.sourceHistogram.getImage(), results.sources, wcs, False)
        if self.config.doMakeBadSourceHistogram and len(results.badSources) > 0:
            results.badSourceHistogram, wcs = apUtils.createImageCoveringSkyTile(
                qsp, skyTileId, self.config.sourceHistogramResolution)
            apUtils.makeSourceHistogram(
                results.badSourceHistogram.getImage(), results.badSources, wcs, False)
        return results

    @pipeBase.timeMethod
    def attributes(self, skyTileId, clusters, exposures):
        """Compute source cluster attributes for a sky-tile.

           @param skyTileId: Integer sky-tile ID
           @param clusters:  List of lsst.afw.table.SourceCatalog objects,
                             each containing the sources for one cluster
           @param exposures: A lsst.ap.match.ExposureInfoMap object, mapping
                             calibrated exposure IDs to lsst.ap.match.ExposureInfo
                             objects.

           @return An lsst.ap.cluster.SourceClusterCatalog containing measurement
                   means for each cluster.
        """
        if len(clusters) == 0:
            return None
        self.log.info(str.format("Computing attributes for {} clusters", len(clusters)))
        spControl = self.config.sourceProcessing.makeControl()
        minNeighbors = self.config.clustering.makeControl().minNeighbors
        scTable = apCluster.makeSourceClusterTable(
            clusters[0].getTable(),
            apCluster.SourceClusterIdFactory(skyTileId),
            spControl)
        flagNoiseKey = scTable.getSchema().find("flag.noise").key
        scCat = apCluster.SourceClusterCatalog(scTable)
        algorithmFlags = dict()
        for alg in spControl.fluxFields:
            if alg in clusters[0].getSchema():
                algorithmFlags[alg] = _flagKeys(clusters[0].getSchema(), self.config, alg)
        for alg in spControl.shapeFields:
            if alg in clusters[0].getSchema():
                algorithmFlags[alg] = _flagKeys(clusters[0].getSchema(), self.config, alg)
        numNoise = 0
        for sources in clusters:
            if len(sources) == 1 and minNeighbors > 0:
                numNoise += 1
                if self.config.doDiscardNoiseClusters:
                    continue
                else:
                    sc = scCat.addNew()
                    sc.set(flagNoiseKey, True)
            else:
                sc = scCat.addNew()
            sev = apCluster.computeBasicAttributes(
                sc, sources, exposures, spControl.exposurePrefix)
            for alg in spControl.fluxFields:
                if alg in algorithmFlags:
                    apCluster.computeFluxMean(sc, sev, alg, algorithmFlags[alg],
                                              spControl.fluxScale)
            for alg in spControl.shapeFields:
                if alg in algorithmFlags:
                    apCluster.computeShapeMean(sc, sev, alg, algorithmFlags[alg]) 
            apCluster.setClusterFields(sources, sc, spControl)
        msg = "Computed attributes for {} clusters"
        if self.config.doDiscardNoiseClusters:
            msg += ", discarded {} noise clusters"
        else:
            msg += ", including {} noise clusters"
        self.log.info(str.format(msg, len(scCat), numNoise))
        return scCat

    @pipeBase.timeMethod
    def run(self, skyTileId, butler):
        """Run source association on a single sky-tile; return None.
        """
        self.log.info(str.format("Processing sky-tile {}", skyTileId))
        res = self.cluster(skyTileId, butler)
        if (self.config.doCluster and res.clusters != None and
            len(res.clusters) > 0):
            clusters = self.attributes(skyTileId, res.clusters, res.exposures)
            # persist clusters
            if self.config.doWriteClusters and len(clusters) > 0:
                butler.put(clusters, "object", skyTile=skyTileId)
        # persist sources
        if (self.config.doWriteSources and res.sources != None and
            len(res.sources) > 0):
            butler.put(res.sources, "source", skyTile=skyTileId)
        if (self.config.doWriteBadSources and res.badSources != None and
            len(res.badSources) > 0):
            butler.put(res.badSources, "badSource", skyTile=skyTileId)
        if (self.config.doWriteInvalidSources and res.invalidSources != None and
            len(res.invalidSources) > 0):
            butler.put(res.invalidSources, "invalidSource", skyTile=skyTileId)
        # persist source histograms
        if self.config.doWriteSourceHistogram and res.sourceHistogram != None:
            butler.put(res.sourceHistogram, "sourceHist", skyTile=skyTileId)
        if (self.config.doWriteBadSourceHistogram and
            res.badSourceHistogram != None):
            butler.put(res.badSourceHistogram, "badSourceHist", skyTile=skyTileId)

    @classmethod
    def _makeArgumentParser(cls):
        """Create an argument parser
        """
        return SourceAssocArgumentParser(name=cls._DefaultName,
                datasetType="source")

    @classmethod
    def parseAndRun(cls, args=None, config=None, log=None):
        """Parse argument list and run command.

            @param args:   list of command-line arguments;
                           if None use sys.argv
            @param config: config for task (instance of lsst.pex.config.Config);
                           if None use cls.ConfigClass()
            @param log:    log (instance of lsst.pex.logging.Log);
                           if None use the default log
        """
        argumentParser = cls._makeArgumentParser()
        if config is None:
            config = cls.ConfigClass()
        parsedCmd = argumentParser.parse_args(config=config, args=args, log=log)
        name = cls._DefaultName
        task = cls(name=name, config=parsedCmd.config, log=parsedCmd.log)
        if not hasattr(parsedCmd, "skyTileIds") or len(parsedCmd.skyTileIds) == 0:
            print >>sys.stderr, "Running on all sky-tiles"
            parsedCmd.skyTileIds = parsedCmd.butler.queryMetadata(
                "source", "skyTile")
        for skyTileId in parsedCmd.skyTileIds:
            parsedCmd.butler.put(
                parsedCmd.config, name + "_config", skyTile=skyTileId)
            try:
                task.run(skyTileId, parsedCmd.butler)
            except Exception, e:
                if parsedCmd.doraise:
                    raise
                task.log.fatal(str.format(
                    "Failed on skyTile {}: {}", skyTileId, e))
                traceback.print_exc(file=sys.stderr)
            parsedCmd.butler.put(
                task.getFullMetadata(), name + "_metadata", skyTile=skyTileId)

