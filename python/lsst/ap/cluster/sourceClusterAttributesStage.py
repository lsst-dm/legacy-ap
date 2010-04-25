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

        # compute  SourceClusterAttributes for each cluster
        self.log.log(Log.INFO, "Computing source cluster attributes")
        sourceClusterAttributes = clusterLib.SourceClusterAttributesSet()
        for sequenceNum, sources in enumerate(sourceClusters):
            clusterId = sequenceNum + (skyTileId << 32)
            attributes = clusterLib.SourceClusterAttributes(sources, clusterId)
            clusterLib.updateSources(attributes, sources)
            sourceClusterAttributes.append(attributes)
        self.log.log(Log.INFO, "Finished computing source cluster attributes")

        # output products
        clipboard.put(self.policy.get("outputKeys.sourceClusterAttributes"),
                      sourceClusterAttributes)


class SourceClusterAttributesStage(stage.Stage):
    parallelClass = SourceClusterAttributesParallel

