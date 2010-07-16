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

from itertools import izip
import operator

import lsst.daf.base as base
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
        # extract input data and parameters
        skyTileIdKey = self.policy.getString("inputKeys.skyTileId")
        jobIdentity = clipboard.get(
            self.policy.getString("inputKeys.jobIdentity"))
        if isinstance(jobIdentity, base.PropertySet):
            skyTileId = jobIdentity.getInt(skyTileIdKey)
        else:
            skyTileId = jobIdentity[skyTileIdKey]
            if not isinstance(skyTileId, (int, long)):
                raise TypeError("Sky-tile id must be an integer")
        sourceClusters = clipboard.get(
            self.policy.getString("inputKeys.sourceClusters"))
        badSourcesKey = self.policy.getString("inputKeys.badSources")
        if clipboard.contains(badSourcesKey):
            badSources = clipboard.get(badSourcesKey).getSources()
        else:
            badSources = None
        fluxIgnoreMask = reduce(
            operator.__or__, self.policy.getArray("fluxIgnoreMask"), 0)
        ellipticityIgnoreMask = reduce(
            operator.__or__, self.policy.getArray("ellipticityIgnoreMask"), 0)
        createBadClusters = self.policy.getBool("createBadClusters")
        scpKey = self.policy.getString("inputKeys.sourceClusteringPolicy")
        if clipboard.contains(scpKey):
            scp = clipboard.get(scpKey)
        else:
            scp = policy.Policy()
        scp.mergeDefaults(self.policy.getPolicy("sourceClusteringPolicy"))
        minPoints = scp.getInt("minPoints")

        # compute SourceClusterAttributes for each cluster
        self.log.log(Log.INFO, "Computing source cluster attributes")
        scv = clusterLib.SourceClusterVector()
        sequenceNum = 0
        numNoise = 0
        for sources in sourceClusters:
            clusterId = sequenceNum + (skyTileId << 32)
            sequenceNum += 1
            sca = clusterLib.SourceClusterAttributes(
                sources, clusterId, fluxIgnoreMask, ellipticityIgnoreMask)
            if len(sources) == 1 and minPoints > 0:
                numNoise += 1
                sca.setFlags(sca.getFlags() |
                             clusterLib.SourceClusterAttributes.NOISE)
            clusterLib.updateSources(sca, sources)
            scv.append(sca)
        if len(sourceClusters) > 0:
            clipboard.put(self.policy.get("outputKeys.sourceClusterAttributes"),
                          clusterLib.PersistableSourceClusterVector(scv))
        self.log.log(Log.INFO,
            "Computed source cluster attributes for %d (%d noise) clusters" %
            (len(sourceClusters), numNoise))

        # create clusters from bad sources
        if badSources != None and len(badSources) > 0:
            if createBadClusters:
                self.log.log(Log.INFO,
                             "Creating source clusters for bad sources")
                badScv = clusterLib.SourceClusterVector()
                for source in badSources.getSources():
                    clusterId = sequenceNum + (skyTileId << 32)
                    sequenceNum += 1
                    badSca = clusterLib.SourceClusterAttributes(
                        source, clusterId, fluxIgnoreMask, ellipticityIgnoreMask)
                    badSca.setFlags(badSca.getFlags() |
                                    clusterLib.SourceClusterAttributes.BAD)
                    badScv.append(badSca)
                clipboard.put(
                    self.policy.get("outputKeys.badSourceClusterAttributes"),
                    clusterLib.PersistableSourceClusterVector(badScv))
                self.log.log(Log.INFO, "Created %d bad source clusters" %
                             len(badSources))
            else:
                clusterLib.updateBadSources(badSources)


class SourceClusterAttributesStage(stage.Stage):
    parallelClass = SourceClusterAttributesParallel

