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
import sys
from tempfile import NamedTemporaryFile
import unittest

import lsst.utils.tests as utilsTests
import lsst.daf.base as dafBase
import lsst.pex.policy as pexPolicy
import lsst.afw.coord as afwCoord
import lsst.ap.match as apMatch


def buildPoints(refFile, exposure):
    id = 0
    nMatches = 0
    points = []
    matches = {}
    for x in xrange(-100, exposure.getWidth() + 100, 50):
        for y in xrange(-100, exposure.getHeight() + 100, 50):
            x += 1.0e-5
            y += 1.0e-5
            if (x >= -0.5 and x <= exposure.getWidth() - 0.5 and
                y >= -0.5 and y <= exposure.getHeight() - 0.5):
                # (x y) falls on the exposure
                cov = 1
                nMatches += 1
            else:
                cov = 0
            p = exposure.getWcs().pixelToSky(x, y)
            entry= (id,
                    p.getLongitude(afwCoord.DEGREES),
                    p.getLatitude(afwCoord.DEGREES),
                    cov,
                    x,
                    y
                   )
            points.append(entry)
            if cov == 1:
                matches[id] = entry
            id += 1
    points.sort(key=lambda x: x[2])
    for p in points:
        refFile.write('%d,%r,%r,%d,%d,%d\n' % p)
    refFile.flush()
    return matches


class ReferenceFilterTestCase(unittest.TestCase):
    """Tests the reference catalog to exposure matching implementation.
    """
    def setUp(self):
        ps = dafBase.PropertySet()
        ps.setInt("scienceCcdExposureId", 0)
        ps.setString("RADESYS", "ICRS")
        ps.setString("FILTER", "u")
        ps.setInt("NAXIS1", 1000)
        ps.setInt("NAXIS2", 1000)
        ps.setDouble("EQUINOX", 2000.0)
        ps.setString("CTYPE1", "RA---TAN")
        ps.setString("CTYPE2", "DEC--TAN")
        ps.setString("CUNIT1", "deg")
        ps.setString("CUNIT2", "deg")
        ps.setDouble("CRPIX1", 500.5)
        ps.setDouble("CRPIX2", 500.5)
        ps.setDouble("CRVAL1", 0.0)
        ps.setDouble("CRVAL2", 0.0)
        ps.setDouble("CD1_1", 0.001)
        ps.setDouble("CD1_2", 0.0)
        ps.setDouble("CD2_1", 0.0)
        ps.setDouble("CD2_2", 0.001)
        ps.setString("TIME-MID", "2000-01-01T11:59:28.000000000Z")
        ps.setDouble("EXPTIME", 10.0)
        self.exposures = apMatch.ExposureInfoVector()
        self.exposures.append(apMatch.ExposureInfo(ps))
        self.refPolicy = pexPolicy.Policy()
        for f in ("refObjectId", "ra", "decl", "uCov", "x", "y"):
            self.refPolicy.add("fieldNames", f)
        for f in ("refObjectId", "ra", "decl", "uCov", "x", "y"):
            self.refPolicy.add("outputFields", f)

    def tearDown(self):
        # Without this, memory is leaked
        del self.refPolicy
        del self.exposures

    
    def testFilter(self):
        with NamedTemporaryFile() as filtFile:
            with NamedTemporaryFile() as refFile:
                matches = buildPoints(refFile, self.exposures[0])
                #import pdb
                #pdb.set_trace()
                apMatch.referenceFilter(refFile.name, filtFile.name,
                                        self.exposures, self.refPolicy)
                if False:
                    # print out the generated test data
                    for line in open(refFile.name, 'rb'):
                        sys.stdout.write(line)
                    print '\n---------\n'
                    for line in open(filtFile.name, 'rb'):
                        sys.stdout.write(line)
                    print '\n---------\n'
                    sys.stdout.flush()

            # Validate that the output match file is correct
            decl = float('-INF')
            with open(filtFile.name, 'rb') as checkFile:
                reader = csv.reader(checkFile,
                                    delimiter=',',
                                    escapechar=None,
                                    quoting=csv.QUOTE_NONE)
                for row in reader:
                    refId = int(row[0])
                    d = float(row[2])
                    # check that output is declination sorted
                    self.assertTrue(d >= decl)
                    decl = d
                    refCov = int(row[3])
                    filtCov = [int(row[i]) for i in xrange(6, 12)]
                    self.assertTrue(refId in matches)
                    if refId in matches:
                        del matches[refId]
                    self.assertTrue(all(filtCov[i] == 0 for i in xrange(1,6)))
                    self.assertEqual(refCov, filtCov[0])
                for m in matches:
                    print matches[m] 
                self.assertEqual(len(matches), 0)


def suite():
    """Returns a suite containing all the test cases in this module.
    """
    utilsTests.init()
    suites = map(unittest.makeSuite,
        [ReferenceFilterTestCase,
         utilsTests.MemoryTestCase
        ])
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests.
    """
    utilsTests.run(suite(), shouldExit)

if __name__ == '__main__':
    run(True)

