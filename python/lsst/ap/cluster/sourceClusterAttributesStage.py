from itertools import izip

import lsst.pex.harness.stage as stage
import lsst.pex.policy as policy
import lsst.afw.detection as detection
import clusterLib 

from lsst.pex.logging import Log


class SourceClusterAttributesParallel(stage.ParallelProcessing):
    """
    Description:
        Computes and outputs source cluster attributes, i.e. objects.

    Policy Dictionary:
        policy/SourceClusterAttributesStageDictionary.paf

    Clipboard Input:
        See inputKeys in policy/SourceClusterAttributesStageDictionary.paf

    ClipboardOutput:
        See outputKeys in policy/SourceClusterAttributesStageDictionary.paf
    """
    def setup(self):
        self.log = Log(self.log, "lsst.ap.cluster")
        policyFile = policy.DefaultPolicyFile(
            "ap", "SourceClusterAttributesStageDictionary.paf", "policy")
        defaultPolicy = policy.Policy.createPolicy(
            policyFile, policyFile.getRepositoryPath())
        if self.policy is None:
            self.policy = policy.Policy()
        self.policy.mergeDefaults(defaultPolicy)

    def process(self, clipboard):
        # retrieve clipboard data
        event = clipboard.get(self.policy.getString("inputKeys.event"))
        skyTileId = event.get("skyTileId");
        sourceClusters = clipboard.get(
            self.policy.getString("inputKeys.sourceClusters"))

        # create a SourceClusterAttributes object for each cluster
        self.log.log(Log.INFO, "Creating cluster attributes objects")
        sourceClusterAttributes = clusterLib.SourceClusterAttributesSet
        for sequenceNum, sources in enumerate(sourceClusters):
            sc = clusterLib.SourceClusterAttributes()
            clusterLib.setClusterId(sc, sources, sequenceNum + (skyTileId << 32))
            sourceClusterAttributes.append(sc)

        # compute position and proper motion estimate for each cluster
        self.log.log(Log.INFO, "Computing cluster positions and proper motions")
        for sc, sources in izip(sourceClusterAttributes, sourceClusters):
            clusterLib.computePositionAndVelocity(sc, sources)

        # compute per-filter properties
        for sc, sources in izip(sourceClusterAttributes, sourceClusters):
            # split sources into
            # compute average PSF magnitude
            # compute average shape measurements
            pass

        # output products
        clipboard.put(self.policy.get("outputKeys.sourceClusterAttributes"),
                      sourceClusterAttributes)


class SourceClusterAttributesStage(stage.Stage):
    parallelClass = SourceClusterAttributesParallel

