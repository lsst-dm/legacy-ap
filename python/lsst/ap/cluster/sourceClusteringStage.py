import math
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
        # create a sky-tile from clipboard data
        event = clipboard.get(self.policy.getString("inputKeys.event"))
        qs = skypix.createQuadSpherePixelization(
            self.policy.getPolicy("quadSpherePolicy"))
        skyTileId = event.get("skyTileId");
        root, x, y = qs.coords(skyTileId);
        skyTile = clusterLib.PT1SkyTile(qs.resolution, root, x, y, skyTileId)
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
            prunedSources[len(prunedSources):] = ss
        del sourceSets
        del inputSources
        self.log.log(Log.INFO, "Discarded %d of %d sources; %d sources remain" %
                     (totalSources - sourcesInTile, totalSources, sourcesInTile))

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


class SourceClusteringStage(stage.Stage):
    parallelClass = SourceClusteringParallel

