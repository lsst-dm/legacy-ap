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

    Clipboard Input:
        Note that clipboard keys are configurable via policy. The
        following string-valued policy parameters must be present
        and set to the name of a clipboard item:

        "inputKeys.event":
            Key for a PropertySet containing properties for the
            pipeline trigger event. Must contain an integer-valued
            property named "skyTileId". 

        "inputKeys.sources":
            Key for the sources to cluster (read by an InputStage).
            Sources may be stored as a lsst.afw.detection.SourceSet,
            lsst.afw.detection.PersistableSourceVector, or a 
            lsst.afw.detection.PersistableSourceVectorVector.

    ClipboardOutput:
        Again, clipboard keys are configurable via policy. The following
        string-valued policy parameters must be present and set to the
        name of the desired output clipboard item:

        "outputKeys.skyTile":
        "outputKeys.clusters":
        "outputKeys.sources":
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
        sources = clipboard.get(self.policy.getString("inputKeys.sources"))

        # discard sources outside the current sky-tile
        self.log.log(Log.INFO, "discarding sources lying outside the sky-tile")
        if isinstance(sources, detection.SourceSet):
            sources = [ sources ]
        elif isinstance(sources, detection.PersistableSourceVector):
            sources = [ sources.getSources() ]
        elif isinstance(sources, detection.PersistableSourceVectorVector):
            sources = [ p.getSources() for p in sources ]
        else:
            raise TypeError(
                "Expecting sources as a %s, %s, or %s - got a %s" %
                (t.__class__.__module__ + '.' + t.__class__.__name__ for t in
                 (detection.SourceSet, detection.PersistableSourceVector,
                  detection.PersistableSourceVectorVector, sources)))
        prunedSources = detection.SourceSet()
        sourcesInTile, totalSources = 0, 0
        for s in sources:
            totalSources += len(s)
            skyTile.prune(s)
            sourcesInTile += len(s)
            prunedSources.append(s.begin(), s.end())
        del sources
        self.log.log(Log.INFO, "Discarded %d of %d sources; %d sources remain" %
                     (totalSources - sourcesInTile, totalSources, sourcesInTile))

        # Cluster the remaining sources
        self.log.log(Log.INFO, "Clustering sources")
        clusters = cluster.cluster(prunedSources)
        self.log.log(Log.INFO, "Finished clustering sources")

        # output products
        clipboard.put(self.policy.get("outputKeys.skyTile"), skyTile)
        clipboard.put(self.policy.get("outputKeys.clusters"), clusters)
        outputSources = detection.PersistableSourceVector()
        outputSources.setSources(prunedSources)
        clipboard.put(self.policy.get("outputKeys.sources"), outputSources)


class SourceClusteringStage(harnessStage.Stage):
    parallelClass = SourceClusteringParallel

