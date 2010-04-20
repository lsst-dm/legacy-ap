#! /usr/bin/env python
import math
import lsst.pex.harness.stage as harnessStage

from lsst.pex.logging import Log

import lsst.pex.policy as policy
import lsst.afw.detection as detection
import lsst.skypix as skypix
import lsst.ap.cluster as cluster


class SourceClusteringParallel(harnessStage.ParallelProcessing):
    """
    Description:
        Discards input sources falling outside the current sky tile
        and runs the OPTICS clustering algorithm over the remaining
        sources.

    Policy Dictionary:
        policy/SourceClusteringStageDictionary.paf
        (references policy/SourceClusteringDictionary.paf)

    Clipboard Input:
        See inputKeys in policy/SourceClusteringStageDictionary.paf

    ClipboardOutput:
        See outputKeys in policy/SourceClusteringStageDictionary.paf
    """
    def setup(self):
        self.log = Log(self.log, "SourceClusteringParallel")
        policyFile = policy.DefaultPolicyFile(
            "ap", "SourceClusteringStageDictionary.paf", "policy")
        defaultPolicy = policy.Policy.createPolicy(
            policyFile, policyFile.getRepositoryPath())
        if self.policy is None:
            self.policy = policy.Policy()
        self.policy.mergeDefaults(defaultPolicy)

    def process(self, clipboard):
        # create a sky-tile from clipboard data
        event = clipboard.get(self.policy.getString("inputKeys.event"))
        qs = skypix.createQuadSpherePixelization(self.policy)
        skyTileId = event.get("skyTileId");
        root, x, y = qs.coords(skyTileId);
        skyTile = apCluster.PT1SkyTile(qs.resolution, root, x, y, skyTileId)
        inputSources = clipboard.get(self.policy.getString("inputKeys.sources"))

        # discard sources outside the current sky-tile
        self.log.log(Log.INFO, "discarding sources lying outside the sky-tile")
        if isinstance(inputSources, detection.SourceSet):
            sourceSets = [ inputSources ]
        elif isinstance(inputSources, detection.PersistableSourceVector):
            sourceSets = [ inputSources.getSources() ]
        elif isinstance(inputSources, detection.PersistableSourceVectorVector):
            sourceSets = [ psv.getSources() for psv in inputSources ]
        else:
            raise TypeError(
                "Expecting sources as a %s, %s, or %s - got a %s" %
                (t.__class__.__module__ + '.' + t.__class__.__name__ for t in
                 (detection.SourceSet, detection.PersistableSourceVector,
                  detection.PersistableSourceVectorVector, sources)))
        prunedSources = detection.SourceSet()
        sourcesInTile, totalSources = 0, 0
        for ss in sourceSets:
            totalSources += len(ss)
            skyTile.prune(ss)
            sourcesInTile += len(ss)
            prunedSources.append(ss.begin(), ss.end())
        del sourceSets
        del inputSources
        self.log.log(Log.INFO, "Discarded %d of %d sources; %d sources remain" %
                     (totalSources - sourcesInTile, totalSources, sourcesInTile))

        # cluster the remaining sources
        self.log.log(Log.INFO, "Clustering sources")
        sourceClusters = cluster.cluster(
            prunedSources, self.policy.getPolicy("sourceClusteringPolicy"))
        self.log.log(Log.INFO, "Finished clustering sources")

        # output products
        clipboard.put(self.policy.get("outputKeys.skyTile"), skyTile)
        clipboard.put(self.policy.get("outputKeys.sourceClusters"),
                      sourceClusters)
        outputSources = detection.PersistableSourceVector()
        outputSources.setSources(prunedSources)
        clipboard.put(self.policy.get("outputKeys.sources"), outputSources)


class SourceClusteringStage(harnessStage.Stage):
    parallelClass = SourceClusteringParallel

