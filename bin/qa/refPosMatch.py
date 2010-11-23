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

import lsst.pex.policy as pexPolicy
import lsst.ap.match as apMatch


def main():
    parser = optparse.OptionParser(dedent("""\
        %prog [options] <ref catalog> <pos table> <output file>

        Matches a reference catalog against a table of positions.
        Both the reference catalog and position table must be
        in increasing declination order.
        """))
    parser.add_option(
        "-R", "--ref-policy", type="string", dest="refPolicy",
        help=dedent("""\
        Reference catalog policy file - defaults are taken from
        policy/ReferenceCatalogDictionary.paf in the ap package."""))
    parser.add_option(
        "-P", "--pos-policy", type="string", dest="posPolicy",
        help=dedent("""\
        Position table policy file - defaults are taken from
        policy/PositionTableDictionary.paf in the ap package."""))
    parser.add_option(
        "-M", "--match-policy", type="string", dest="matchPolicy",
        help=dedent("""\
        Match policy file - defaults are taken from
        policy/ReferenceMatchDictionary.paf in the ap package."""))
    parser.add_option(
        "-r", "--radius", type="float", dest="radius", help=dedent("""\
        Specifies the search radius (arcsec) for the reference catalog
        to position table match. Overrides the value stored in the match
        policy file."""))
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
        place. This reduction should be enabled for higher accuracy
        matching against sources, but not when matching against source
        clusters (the latter combine data from multiple epochs)."""))
    parser.add_option(
        "-f", "--fields", dest="fields", help=dedent("""\
        A comma separated list of the fields names in the position table.
        If omitted, the first line in the position table is expected to
        contain field names."""))

    opts, args = parser.parse_args()
    if len(args) != 3:
        parser.error("Three arguments expected")
    if opts.refPolicy != None:
        refPolicy = pexPolicy.Policy(opts.refPolicy)
    else:
        refPolicy = pexPolicy.Policy()
    if opts.posPolicy != None:
        posPolicy = pexPolicy.Policy(opts.posPolicy)
    else:
        posPolicy = pexPolicy.Policy()
    if opts.matchPolicy != None:
        matchPolicy = pexPolicy.Policy(opts.matchPolicy)
    else:
        matchPolicy = pexPolicy.Policy()
    if opts.radius != None:
        if opts.radius < 0.0:
            parser.error("Negative search radius")
        matchPolicy.set("radius", opts.radius)
    if opts.parallaxThresh != None:
        matchPolicy.set("parallaxThresh", opts.parallaxThresh)
    if opts.noSsbToGeo:
        matchPolicy.set("parallaxThresh", float("INF"))
    if opts.fields != None:
        posPolicy.remove("fieldNames")
        for n in opts.fields.split(","):
            posPolicy.add("fieldNames", n.strip())
    apMatch.referenceMatch(args[0], args[1], args[2],
                           refPolicy, posPolicy, matchPolicy)

if __name__ == "__main__":
    main()

