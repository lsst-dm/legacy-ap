import math
import operator
from textwrap import dedent

import lsst.pex.harness.stage as stage
import lsst.pex.policy as policy
import lsst.afw.detection as detection
import lsst.skypix as skypix
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
        event = clipboard.get(self.policy.getString("inputKeys.event"))
        qs = skypix.createQuadSpherePixelization(
            self.policy.getPolicy("quadSpherePolicy"))
        skyTileId = event.getInt("skyTileId");
        root, x, y = qs.coords(skyTileId);
        skyTile = clusterLib.PT1SkyTile(qs.resolution, root, x, y, skyTileId)
        inputSources = clipboard.get(self.policy.getString("inputKeys.sources"))
        badSourceMask = reduce(
            operator.__or__, self.policy.getIntArray('badSourceMask'), 0)

        # turn input into a list of source sets
        sourceSets = []
        if not isinstance(inputSources, (list, tuple)):
            inputSources = [ inputSources ]
        for entry in inputSources:
            if isinstance(entry, detection.SourceSet):
                sourceSets.append(entry)
            elif isinstance(entry, detection.PersistableSourceVector):
                sourceSets.append(entry.getSources())
            else:
                raise TypeError(
                    "Expecting sources in a [list/tuple of] %s or %s - got a %s" %
                    (t.__class__.__module__ + '.' + t.__class__.__name__ for t in
                     (detection.SourceSet, detection.PersistableSourceVector,
                      sources)))

        # remove sources with invalid positions
        self.log.log(Log.INFO, "Segregating sources with invalid positions")
        invalidSources = detection.SourceSet()
        total = 0
        for ss in sourceSets:
            total += len(ss)
            clusterLib.segregateInvalidSources(ss, invalidSources)
        n = len(invalidSources)
        self.log.log(Log.INFO,
            "%d of %d sources are invalid; %d sources remain" %
            (n, total, total - n))

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
        del inputSources
        self.log.log(Log.INFO, dedent("""\
            %d of %d valid sources are outside the current sky-tile;
            %d sources remain""") % (total - n, total, n))

        # do not feed "bad" sources to clustering algorithm
        self.log.log(Log.INFO, "Segregating bad sources")
        badSources = detection.SourceSet()
        clusterLib.segregateBadSources(prunedSources, badSources, badSourceMask)
        self.log.log(Log.INFO, dedent("""\
            %d of %d valid sources in the current sky-tile are bad; %d
            sources remain""") % (len(badSources), total, len(prunedSources)))

        # cluster the remaining sources
        self.log.log(Log.INFO, "Clustering sources")
        sourceClusters = clusterLib.cluster(
            prunedSources, self.policy.getPolicy("sourceClusteringPolicy"))
        self.log.log(Log.INFO, "Finished clustering sources")

        # output products
        clipboard.put(self.policy.get("outputKeys.sourceClusters"),
                      sourceClusters)
        outputSources = detection.PersistableSourceVector()
        outputSources.setSources(prunedSources)
        clipboard.put(self.policy.get("outputKeys.sources"), outputSources)
        clipboard.put(self.policy.get("outputKeys.sourceClusteringPolicy"),
                      self.policy.getPolicy("sourceClusteringPolicy"))
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

