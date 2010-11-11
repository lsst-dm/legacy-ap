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
from __future__ import with_statement

from contextlib import nested
import csv
import math
import optparse
import os
import pdb
import random
import sys
from tempfile import NamedTemporaryFile
import unittest

import lsst.utils.tests as utilsTests
import lsst.geom as geom
import lsst.pex.policy as pexPolicy
import lsst.ap.match as apMatch


MAX_MATCHES = 8 # max no. of matches to generate per test position
MATCH_ACCURACY = 0.000001/3600.0 # 1 microarcsec in degrees


def pointInBox(minTheta, maxTheta, minPhi, maxPhi):
    """Returns spherical coordinates (theta, phi) for a point chosen
    approximately uniformly at random from within the specified
    longitude/latitude angle box.
    """
    # pick theta
    if minTheta <= maxTheta:
        theta = minTheta + random.random() * (maxTheta - minTheta)
        if theta > maxTheta:
            theta = maxTheta;
    else:
        # wrap-around
        m = minTheta - 360.0
        theta = m + random.random() * (maxTheta - m);
        if theta < 0:
            theta += 360.0
    # pick phi
    minZ = math.sin(math.radians(minPhi))
    maxZ = math.sin(math.radians(maxPhi))
    z = minZ + random.random() * (maxZ - minZ)
    phi = math.degrees(math.asin(z))
    if phi < minPhi:
        phi = minPhi
    elif phi > maxPhi:
        phi = maxPhi
    return theta, phi


def pointsOnCircle(c, r, n):
    """Generates an n-gon lying on the circle with center c and
    radius r. Vertices are equi-spaced.
    """
    points = []
    c = geom.cartesianUnitVector(c)
    north, east = geom.northEast(c)
    sr = math.sin(math.radians(r))
    cr = math.cos(math.radians(r))
    aoff = random.uniform(0.0, 2.0 * math.pi)
    for i in xrange(n):
        a = 2.0 * i * math.pi / n
        sa = math.sin(a + aoff)
        ca = math.cos(a + aoff)
        p = (ca * north[0] + sa * east[0],
             ca * north[1] + sa * east[1],
             ca * north[2] + sa * east[2])
        points.append(geom.sphericalCoords((cr * c[0] + sr * p[0],
                                            cr * c[1] + sr * p[1],
                                            cr * c[2] + sr * p[2])))
    return points


def buildPoints(refFile, posFile, radius):
    """Builds test data for point within circle matches.
    """
    assert radius > 0.0 and radius < 5.0

    # Divide unit sphere into latitude angle stripes
    phiMin = -90.0
    phiMax = 90.0
    i = 0
    deltaPhi = 4.0 * radius;
    phi = phiMin
    while phi < phiMax:
        centerPhi = geom.clampPhi(max(abs(phi), abs(phi + deltaPhi)))
        deltaTheta = geom.maxAlpha(4.0 * radius, centerPhi)
        theta = 0.0
        refOutput = []
        posOutput = []
        # Divide latitude angle stripes into boxes (by longitude angle)
        while theta < 360.0 - 2.0 * deltaTheta:
            # Create a random point inside a sub-region of each box
            # such that a circle of the given radius centered on that
            # point is guaranteed not to cross the box boundaries
            if theta == 0.0:
                # make sure longitude angle wrap-around is tested
                p = pointInBox(360.0 - 0.125 * deltaTheta,
                               0.125 * deltaTheta,
                               geom.clampPhi(phi + deltaPhi * 0.38),
                               geom.clampPhi(phi + deltaPhi * 0.62))
            else:
                p = pointInBox(theta + deltaTheta * 0.38,
                               theta + deltaTheta * 0.62,
                               geom.clampPhi(phi + deltaPhi * 0.38),
                               geom.clampPhi(phi + deltaPhi * 0.62))
            refId = i
            i += 1

            # Generate matches
            numMatches = random.randint(0, MAX_MATCHES)
            pIn = pointsOnCircle(p, radius - MATCH_ACCURACY, numMatches)
            pOut = pointsOnCircle(p, radius + MATCH_ACCURACY, numMatches)
            if len(pIn) > 1:
                # shift first point inwards towards p to make it the
                # closest match
                v1 = geom.cartesianUnitVector(p)
                v2 = geom.cartesianUnitVector(pIn[0])
                v3 = geom.normalize((v1[0] + v2[0], v1[1] + v2[1], v1[2] + v2[2]))
                pIn[0] = geom.sphericalCoords(v3)
            matches = ' '.join(map(str, xrange(i, i + len(pIn))))
            # Write out reference position
            refOutput.append(
                (p[1], "%d,%s,%s,%s\n" % (refId, repr(p[0]), repr(p[1]), matches)))
            # Write out matching positions
            for p in pIn:
                posOutput.append(
                    (p[1], "%d,%s,%s,%d\n" % (i, repr(p[0]), repr(p[1]), refId)))
                i += 1
            # Write out positions with no matches
            for p in pOut:
                posOutput.append(
                    (p[1], "%d,%s,%s,\n" % (i, repr(p[0]), repr(p[1]))))
                i += 1
            theta += deltaTheta
        posOutput.sort()
        refOutput.sort()
        for r in refOutput:
            refFile.write(r[1])
        for p in posOutput:
            posFile.write(p[1])
        phi += deltaPhi
    refFile.flush()
    posFile.flush()


class ReferenceMatchTestCase(unittest.TestCase):
    """Tests the reference matching implementation.
    """
    def setUp(self):
        self.radius = float(os.environ['TEST_RM_RADIUS'])
        self.refPolicy = pexPolicy.Policy()
        self.posPolicy = pexPolicy.Policy()
        for f in ("refObjectId", "ra", "decl", "refMatches"):
            self.refPolicy.add("fieldNames", f)
        self.refPolicy.set("outputFields", "refMatches")
        for f in ("posId", "ra", "decl", "posMatches"):
            self.posPolicy.add("fieldNames", f)
        self.posPolicy.set("outputFields", "posMatches")
        self.posPolicy.set("idColumn", "posId")
        self.matchPolicy = pexPolicy.Policy()
        self.matchPolicy.set("radius", self.radius*3600.0)

    def tearDown(self):
        # Without this, memory is leaked
        del self.refPolicy
        del self.posPolicy
        del self.matchPolicy

    def testMatch(self):
        random.seed(123456789)
        with NamedTemporaryFile() as matchFile:
            with nested(NamedTemporaryFile(),
                        NamedTemporaryFile()) as (refFile, posFile):
                buildPoints(refFile, posFile, self.radius)
                apMatch.referenceMatch(
                    refFile.name, posFile.name, matchFile.name,
                    self.refPolicy, self.posPolicy, self.matchPolicy)
                if False:
                    # print out the generated test data
                    for line in open(refFile.name, 'rb'):
                        sys.stdout.write(line)
                    print '\n---------\n'
                    for line in open(posFile.name, 'rb'):
                        sys.stdout.write(line)
                    print '\n---------\n'
                    for line in open(matchFile.name, 'rb'):
                        sys.stdout.write(line)
                    print '\n---------\n'

            # Validate that the output match file is correct
            with open(matchFile.name, 'rb') as checkFile:
                reader = csv.reader(checkFile,
                                    delimiter=',',
                                    escapechar=None,
                                    quoting=csv.QUOTE_NONE)
                refs = {}
                for row in reader:
                    refId = None
                    posId = None
                    angSep = None
                    nRefMatches = None
                    nPosMatches = None
                    isClosestMatchForRef = None
                    isClosestMatchForObj = None
                    refMatches = None
                    posMatches = None
                    if row[0] != '\\N':
                        refId = int(row[0])
                    if row[1] != '\\N':
                        posId = int(row[1])
                    if row[4] != '\\N':
                        angSep = float(row[4])
                    if row[5] != '\\N':
                        nRefMatches = int(row[5])
                    if row[6] != '\\N':
                        nPosMatches = int(row[6])
                    if row[7] != '\\N':
                        isClosestMatchForRef = row[7]
                    if row[8] != '\\N':
                        isClosestMatchForObj = row[8]
                    flags = row[9]
                    if refId in refs:
                        refMatches, alreadyFound = refs[refId]
                    elif row[10] != '\\N':
                        refMatches = []
                        if len(row[10]) > 0:
                            refMatches = map(int, row[10].split(' '))
                        alreadyFound = []
                    if row[11] != '\\N':
                        posMatches = []
                        if len(row[11]) > 0:
                            posMatches = map(int, row[11].split(' '))
                    if refId == None:
                        self.assertNotEqual(posId, None)
                        self.assertEqual(nRefMatches, None)
                        self.assertEqual(nPosMatches, 0)
                        self.assertEqual(angSep, None)
                        self.assertEqual(isClosestMatchForRef, None)
                        self.assertEqual(isClosestMatchForObj, None)
                        self.assertEqual(flags, '\\N')
                        self.assertEqual(refMatches, None)
                        self.assertEqual(posMatches, [])
                    elif posId == None:
                        self.assertNotEqual(refId, None)
                        self.assertEqual(nRefMatches, 0)
                        self.assertEqual(nPosMatches, None)
                        self.assertEqual(angSep, None)
                        self.assertEqual(isClosestMatchForRef, None)
                        self.assertEqual(isClosestMatchForObj, None)
                        self.assertEqual(flags, '\\N')
                        self.assertEqual(refMatches, [])
                        self.assertEqual(posMatches, None)
                    else:
                        # have a match
                        self.assertTrue(nRefMatches > 0 and nPosMatches == 1)
                        self.assertTrue(refMatches != None and posMatches != None)
                        self.assertEqual(len(posMatches), 1)
                        self.assertEqual(refId, posMatches[0])
                        self.assertTrue(angSep != None and angSep < self.radius*3600.0)
                        self.assertEqual(isClosestMatchForObj, "1")
                        self.assertEqual(flags, "0")
                        self.assertTrue(posId in refMatches)
                        self.assertFalse(posId in alreadyFound)
                        self.assertEqual(nRefMatches, len(refMatches))
                        if posId == refMatches[0]:
                            self.assertEqual(isClosestMatchForRef, "1")
                        alreadyFound.append(posId)
                        refs[refId] = (refMatches, alreadyFound)


def suite():
    """Returns a suite containing all the test cases in this module.
    """
    utilsTests.init()
    suites = map(unittest.makeSuite,
        [ReferenceMatchTestCase,
         utilsTests.MemoryTestCase
        ])
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests.
    """
    utilsTests.run(suite(), shouldExit)

if __name__ == '__main__':
    usage = """usage: %prog [options]

    Tests the reference matching implementation
    """
    parser = optparse.OptionParser(usage)
    parser.add_option(
        "-r", "--radius", type="string", dest="radius",
        help="Match radius (deg)")
    (options, args) = parser.parse_args()
    if options.radius:
        os.environ['TEST_RM_RADIUS'] = options.radius
    else:
        os.environ['TEST_RM_RADIUS'] = "1.0"
    run(True)

