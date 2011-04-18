#
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
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

import math
import operator
from textwrap import dedent

import lsst.daf.base as base
import lsst.pex.harness.stage as stage
import lsst.pex.policy as policy
import lsst.afw.detection as detection
import lsst.skypix as skypix
import lsst.ap.match as apMatch
import lsst.ap.utils as apUtils

import clusterLib

from lsst.pex.logging import Log


class SourceClusteringParallel(stage.ParallelProcessing):
    """
    Description:
        Discards input sources falling outside the current sky tile
        and runs the OPTICS clustering algorithm over the remaining
        sources.

    Policy Dictionary:
        policy/SourceClusteringStageDictionary.paf

    Clipboard Input:
        See inputKeys in policy/SourceClusteringStageDictionary.paf

    ClipboardOutput:
        See outputKeys in policy/SourceClusteringStageDictionary.paf
    """
    def setup(self):
        self.log = Log(self.log, "lsst.ap.cluster")
        policyFile = policy.DefaultPolicyFile(
            "ap", "SourceClusteringStageDictionary.paf", "policy")
        defaultPolicy = policy.Policy.createPolicy(
            policyFile, policyFile.getRepositoryPath())
        if self.policy is None:
            self.policy = policy.Policy()
        self.policy.mergeDefaults(defaultPolicy)

    def process(self, clipboard):
        # extract required clipboard data
        skyTileIdKey = self.policy.getString("inputKeys.skyTileId")
        jobIdentity = clipboard.get(
            self.policy.getString("inputKeys.jobIdentity"))
        qs = skypix.createQuadSpherePixelization(
            self.policy.getPolicy("quadSpherePolicy"))
        if isinstance(jobIdentity, base.PropertySet):
            skyTileId = jobIdentity.getInt(skyTileIdKey)
        else:
            skyTileId = jobIdentity[skyTileIdKey]
            if not isinstance(skyTileId, (int, long)):
                raise TypeError("Sky-tile id must be an integer")
        root, x, y = qs.coords(skyTileId);
        skyTile = apUtils.PT1SkyTile(qs.resolution, root, x, y, skyTileId)
        histogramRes = self.policy.getInt("debug.sourceHistogramResolution")
        badSourceMask = reduce(
            operator.__or__, self.policy.getIntArray('badSourceMask'), 0)
        computeSourceSkyCoords = self.policy.getBool('computeSourceSkyCoords')
        inputSources = clipboard.get(self.policy.getString("inputKeys.sources"))
        inputExposures = clipboard.get(
            self.policy.getString("inputKeys.exposures"))
        # Remove input sources and exposure metadata from clipboard
        clipboard.put(self.policy.getString("inputKeys.sources"), None)
        clipboard.put(self.policy.getString("inputKeys.exposures"), None)

        # turn input exposure metadata into a mapping from exposure ids to
        # ExposureInfo objects
        exposures = apMatch.ExposureInfoMap()
        if not isinstance(inputExposures, (list, tuple)):
            inputExposures = [inputExposures]
        for entry in inputExposures:
            exposures.insert(apMatch.ExposureInfo(entry))
        del inputExposures

        # turn input sources into a list of source sets
        sourceSets = []
        if not isinstance(inputSources, (list, tuple)):
            inputSources = [inputSources]
        for entry in inputSources:
             if isinstance(entry, detection.SourceSet):
                 sourceSets.append(entry)
             elif isinstance(entry, detection.PersistableSourceVector):
                 sourceSets.append(entry.getSources())
             else:
                 raise TypeError("Expecting lsst.afw.detection.SourceSet or " +
                                 "lsst.afw.detection.PersistableSourceVector")
        del inputSources

        # remove sources with invalid positions
        self.log.log(Log.INFO, "Computing source positions and/or removing " +
                     "invalid sources")
        invalidSources = detection.SourceSet()
        total = 0
        for ss in sourceSets:
            total += len(ss)
            if computeSourceSkyCoords:
                clusterLib.locateAndFilterSources(ss, invalidSources, exposures)
            else:
                clusterLib.segregateInvalidSources(ss, invalidSources)
        n = len(invalidSources)
        self.log.log(Log.INFO,
            "%d of %d sources are invalid; %d sources remain" %
            (n, total, total - n))
        if n > 0:
            outputInvalidSources = detection.PersistableSourceVector()
            outputInvalidSources.setSources(invalidSources)
            clipboard.put(self.policy.get("outputKeys.invalidSources"),
                          outputInvalidSources)

        # discard sources outside the current sky-tile
        self.log.log(Log.INFO, "Discarding sources lying outside the sky-tile")
        prunedSources = detection.SourceSet()
        n, total = 0, 0
        for ss in sourceSets:
            total += len(ss)
            skyTile.prune(ss)
            n += len(ss)
            prunedSources[len(prunedSources):] = ss
        del sourceSets
        self.log.log(Log.INFO, dedent("""\
            %d of %d valid sources are outside the current sky-tile;
            %d sources remain""") % (total - n, total, n))

        # do not feed "bad" sources to clustering algorithm
        self.log.log(Log.INFO, "Segregating bad sources")
        badSources = detection.SourceSet()
        clusterLib.segregateBadSources(prunedSources, badSources, badSourceMask)
        self.log.log(Log.INFO, dedent("""\
            %d of %d valid sources in the current sky-tile are bad; %d
            sources remain""") % (len(badSources), n, len(prunedSources)))
        if len(badSources) > 0:
            outputBadSources = detection.PersistableSourceVector()
            outputBadSources.setSources(badSources)
            clipboard.put(self.policy.get("outputKeys.badSources"),
                          outputBadSources)
            if self.policy.getBool("debug.createBadSourceHistogram"):
                self.log.log(Log.INFO, "Creating bad source histogram")
                hist, wcs = apUtils.createImageCoveringSkyTile(
                    qs, skyTileId, histogramRes)
                apUtils.makeSourceHistogram(
                    hist.getImage(), badSources, wcs, False)
                clipboard.put(
                    self.policy.getString("outputKeys.badSourceHistogram"),
                    hist)

        # cluster the remaining sources
        self.log.log(Log.INFO, "Clustering sources")
        sourceClusters = clusterLib.cluster(
            prunedSources, self.policy.getPolicy("sourceClusteringPolicy"))
        self.log.log(Log.INFO, "Finished clustering sources")
        # output exposures, clusters and good sources
        clipboard.put(self.policy.get("outputKeys.exposures"), exposures)
        clipboard.put(self.policy.get("outputKeys.sourceClusters"),
                      sourceClusters)
        if len(prunedSources) > 0:
            outputSources = detection.PersistableSourceVector()
            outputSources.setSources(prunedSources)
            clipboard.put(self.policy.get("outputKeys.sources"), outputSources)
            clipboard.put(self.policy.get("outputKeys.sourceClusteringPolicy"),
                          self.policy.getPolicy("sourceClusteringPolicy"))

            if self.policy.getBool("debug.createGoodSourceHistogram"):
                self.log.log(Log.INFO, "Creating good source histogram")
                hist, wcs = apUtils.createImageCoveringSkyTile(
                    qs, skyTileId, histogramRes)
                apUtils.makeSourceHistogram(
                    hist.getImage(), prunedSources, wcs, False)
                clipboard.put(
                    self.policy.getString("outputKeys.goodSourceHistogram"),
                    hist)

        # output invalid/bad sources
        if len(invalidSources) > 0:
            outputInvalidSources = detection.PersistableSourceVector()
            outputInvalidSources.setSources(invalidSources)
            clipboard.put(self.policy.get("outputKeys.invalidSources"),
                          outputInvalidSources)
        if len(badSources) > 0:
            outputBadSources = detection.PersistableSourceVector()
            outputBadSources.setSources(badSources)
            clipboard.put(self.policy.get("outputKeys.badSources"),
                          outputBadSources)


class SourceClusteringStage(stage.Stage):
    parallelClass = SourceClusteringParallel

