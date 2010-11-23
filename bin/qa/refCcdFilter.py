#! /usr/bin/env python

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
import optparse
from textwrap import dedent

import lsst.daf.persistence as dafPersistence
import lsst.pex.logging as pexLogging
import lsst.pex.policy as pexPolicy
import lsst.ap.match as apMatch
from lsst.obs.lsstSim import LsstSimMapper


def getScienceCcdExposureId(visit, raft, sensor):
    r1, comma, r2 = raft
    s1, comma, s2 = sensor
    raftId = int(r1)*5 + int(r2)
    ccdNum = int(s1)*3 + int(s2)
    return (long(visit) << 9) + raftId*10 + ccdNum


def main():
    parser = optparse.OptionParser(dedent("""\
        %prog [options] <ref catalog> <output file> <root> [<registry>]

        Filters a reference catalog against the science CCD exposures in a run.
        The reference catalog must be in increasing declination order.

        Only reference catalog entries observable in at least one exposure are
        written out. For each entry written out, 6 fields are appended to the
        requested output fields: uCov, gCov, rCov, iCov, zCov, yCov. These give
        the number of exposures in each filter (u,g,r,i,z,y) containing the
        entry."""))
    parser.add_option(
        "-R", "--ref-policy", type="string", dest="refPolicy",
        help=dedent("""\
        Optional reference catalog policy file - defaults are taken from
        policy/ReferenceCatalogDictionary.paf in the ap package."""))
    parser.add_option(
        "-M", "--match-policy", type="string", dest="matchPolicy",
        help=dedent("""\
        Optional match policy file - defaults are taken from
        policy/ReferenceMatchDictionary.paf in the ap package."""))
    parser.add_option(
        "-p", "--parallaxThresh", type="float", dest="parallaxThresh",
        help=dedent("""\
        Parallax threshold (milliarcsec) above which a reduction for
        parallax from barycentric to geocentric place is applied to
        reference catalog positions. This option overrides the value
        stored in the match policy file."""))
    parser.add_option(
        "-s", "--no-ssb-to-geo", action="store_true", dest="noSsbToGeo",
        help=dedent("""\
        Disables reduction for parallax from barycentric to geocentric
        place."""))
    parser.add_option(
        "-o", "--output-fields", type="string", dest="outputFields",
        help=dedent("""\
        A comma separated list of reference catalog fields to output.
        Overrides the value of the "outputFields" policy parameter value in
        the reference catalog policy file. If this option is unspecified and
        the reference catalog policy file does not specify "outputFields",
        then all reference catalog fields are output. This can also be
        explictly requested with -f "*"."""))
         
    opts, args = parser.parse_args()
    if not len(args) in (3, 4):
        parser.error("Three or four arguments expected")
    if opts.refPolicy != None:
        refPolicy = pexPolicy.Policy(opts.refPolicy)
    else:
        refPolicy = pexPolicy.Policy()
    if opts.outputFields != None:
        refPolicy.remove("outputFields")
        for n in opts.outputFields.split(","):
            refPolicy.add("outputFields", n.strip())
    if not refPolicy.exists("outputFields"):
        refPolicy.set("outputFields", "*")
    if opts.matchPolicy != None:
        matchPolicy = pexPolicy.Policy(opts.matchPolicy)
    else:
        matchPolicy = pexPolicy.Policy()
    if opts.parallaxThresh != None:
        matchPolicy.set("parallaxThresh", opts.parallaxThresh)
    if opts.noSsbToGeo:
        matchPolicy.set("parallaxThresh", float("INF"))
    # Create butler for run
    root = args[2]
    registry = None
    if len(args) == 4:
        registry = args[3]
    mapper = LsstSimMapper(root=root, registry=registry)
    butler = dafPersistence.ButlerFactory(mapper=mapper).create()
    # Retrieve CCD metadata for all CCDs in run
    log = pexLogging.Log(pexLogging.getDefaultLog(), "lsst.ap.match")
    log.log(log.INFO, "Retrieving science CCD exposure metadata for all CCDs in " + root)
    exposures = apMatch.ExposureInfoVector()
    metadata = butler.queryMetadata("raw", "sensor", ("visit", "raft", "sensor"))
    for visit, raft, sensor in metadata:
        if butler.datasetExists("calexp_md", visit=visit, raft=raft, sensor=sensor):
            ps = butler.get("calexp_md", visit=visit, raft=raft, sensor=sensor)
            ps.set("scienceCcdExposureId", getScienceCcdExposureId(visit, raft, sensor))
            exposures.append(apMatch.ExposureInfo(ps))
    log.log(log.INFO, "Retrieved metadata for %d science CCDs" % len(exposures))
    apMatch.referenceFilter(args[0], args[1], exposures, refPolicy, matchPolicy)

if __name__ == "__main__":
    main()

