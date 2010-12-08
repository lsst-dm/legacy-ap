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
import getpass
from itertools import izip
import MySQLdb as sql
import optparse
import os, os.path
import pdb
from textwrap import dedent

import lsst.daf.base as dafBase
import lsst.daf.persistence as dafPersist
import lsst.pex.policy as pexPolicy
import lsst.afw.cameraGeom as cameraGeom
import lsst.afw.cameraGeom.utils as cameraGeomUtils
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.geom as geom
import lsst.skypix.quadsphere as qs
from lsst.obs.cfht import CfhtMapper
from lsst.obs.lsstSim import LsstSimMapper


def getAllSipWcs(cursor, qsp, kind):
    """Constructs a Wcs object from each entry in the
    Science_Ccd_Exposure table and returns a mapping from skytiles to the
    WCSes and associated identifiers for overlapping science CCDs.
    """
    wcsMap = {}
    wcsList = []
    qsp = qs.createQuadSpherePixelization()
    offset = 0
    blocksize = 1000
    # Getting this from actual calexp output files would be preferrable. 
    # However, runs often occur on clusters far far away and accessing
    # output images isn't always practical. Therefore, pull WCS related
    # cards from run output database.
    while True:
        cursor.execute(
            """SELECT scienceCcdExposureId, visit, raft, ccd, filterId,
                      raDeSys, equinox, ctype1, ctype2,
                      crpix1, crpix2, crval1, crval2,
                      cd1_1, cd1_2, cd2_1, cd2_2
               FROM Science_Ccd_Exposure LIMIT %d, %d
            """ % (offset, blocksize))
        rows = cursor.fetchall()
        if len(rows) == 0:
            break
        # For each CCD in the results, build up a metadata PropertySet
        # including SIP distortion terms and use it to obtain a WCS.
        for row in rows:
            ps = dafBase.PropertySet()
            ps.setString("RADESYS", row[5])
            ps.setDouble("EQUINOX", row[6])
            ps.setString("CTYPE1", row[7])
            ps.setString("CTYPE2", row[8])
            ps.setString("CUNIT1", "deg")
            ps.setString("CUNIT2", "deg")
            ps.setDouble("CRPIX1", row[9])
            ps.setDouble("CRPIX2", row[10])
            ps.setDouble("CRVAL1", row[11])
            ps.setDouble("CRVAL2", row[12])
            ps.setDouble("CD1_1", row[13])
            ps.setDouble("CD1_2", row[14])
            ps.setDouble("CD2_1", row[15])
            ps.setDouble("CD2_2", row[16])
            if row[7].endswith("-SIP") and row[8].endswith("-SIP"):
                cursor.execute("""
                    SELECT metadataKey, intValue
                    FROM Science_Ccd_Exposure_Metadata
                    WHERE scienceCcdExposureId = %d AND
                          intValue IS NOT NULL AND
                          metadataKey RLIKE '^[AB]P?_ORDER$'
                    """ % row[0])
                for k, v in cursor.fetchall():
                    ps.setInt(k, v)
                cursor.execute("""
                    SELECT metadataKey, doubleValue
                    FROM Science_Ccd_Exposure_Metadata
                    WHERE scienceCcdExposureId = %d AND
                          doubleValue IS NOT NULL AND
                          metadataKey RLIKE '^[AB]P?_[0-9]+_[0-9]+$'
                    """ % row[0])
                for k, v in cursor.fetchall():
                    ps.setDouble(k, v)
            wcs = afwImage.makeWcs(ps)
            data = (wcs, row[0], row[1], row[2], row[3], row[4])
            wcsList.append(data)
            if kind == 'imsim':
                poly = qs.imageToPolygon(wcs, 4072.0, 4000.0, 0.0)
            else:
                poly = qs.imageToPolygon(wcs, 2048.0, 4610.0, 0.0)
            pix = qsp.intersect(poly)
            for skyTile in pix:
                if skyTile in wcsMap:
                    wcsMap[skyTile].append(data)
                else:
                    wcsMap[skyTile] = [data]
        print "... %d" % (offset + len(rows))
        offset += blocksize
    return wcsMap, wcsList

def cornersFromBBox(bbox):
    return [ (bbox.getX0() - 0.5, bbox.getY0() - 0.5),
             (bbox.getX1() + 0.5, bbox.getY0() - 0.5),
             (bbox.getX1() + 0.5, bbox.getY1() + 0.5),
             (bbox.getX0() - 0.5, bbox.getY1() + 0.5) ]

def coordsToPolygon(coords):
    verts = [tuple(c.getVector()) for c in coords]
    convex, cc = geom.convex(verts)
    if not convex:
        raise RuntimeError('Image corners do not form a convex polygon: ' + cc)
    elif not cc:
        verts.reverse()
    return geom.SphericalConvexPolygon(verts)

lsstSimRafts = [       "0,1", "0,2", "0,3",
                "1,0", "1,1", "1,2", "1,3", "1,4",
                "2,0", "2,1", "2,2", "2,3", "2,4",
                "3,0", "3,1", "3,2", "3,3", "3,4",
                       "4,1", "4,2", "4,3"
               ]

def getAmps(ccdData, camera, ampDiskLayout, cursor, kind):
    """Gets a list of amps belonging to the specified CCD. For each amp, the following
    data is returned:
      - amp object (from cameraGeom)
      - amp WCS
      - datasec amp corners in the science CCD
      - datasec amp corners in the raw on-disk amp
    """
    visit, raftNum, sensorNum = ccdData[2:5]
    if kind == 'imsim':
        raft = lsstSimRafts[raftNum]
        s1 = sensorNum // 3
        s2 = sensorNum - 3 * s1
        ccdNum = int("%s%s%d%d" % (raft[0], raft[2], s1, s2))
        ccdName = "R:%s S:%d,%d" % (raft, s1, s2)
        ccdId = cameraGeom.Id(-1, ccdName)
    else:
        ccdName = "CFHT %d" % sensorNum
        ccdId = cameraGeom.Id(-1, ccdName)
    ccd = cameraGeomUtils.findCcd(camera, ccdId)
    cursor.execute("""
        SELECT r.amp, r.raDeSys, r.equinox, r.ctype1, r.ctype2,
               r.crpix1, r.crpix2, r.crval1, r.crval2,
               r.cd1_1, r.cd1_2, r.cd2_1, r.cd2_2
        FROM Raw_Amp_Exposure r, Raw_Amp_To_Science_Ccd_Exposure m
        WHERE r.rawAmpExposureId = m.rawAmpExposureId AND
              m.scienceCcdExposureId = %d AND m.snap = 0;
        """ % ccdData[1])
    results = []
    amps = {}
    for row in cursor.fetchall():
        amps[row[0]] = row
    for amp in ccd:
        ampSerial = amp.getId().getSerial()
        flipLR, flipTB = ampDiskLayout[ampSerial]
        chanX, chanY = amp.getId().getIndex()
        ampCorners = cornersFromBBox(amp.getDiskDataSec())
        ampCcdCorners = cornersFromBBox(amp.getDataSec(True))
        nQuarter = ccd.getOrientation().getNQuarter()
        # Deal with flips and rotations performed by CCD assembly.
        if flipLR:
            if flipTB: indexes = (2,3,0,1)
            else:      indexes = (1,0,3,2)
        elif flipTB:   indexes = (3,2,1,0)
        else:          indexes = (0,1,2,3)
        ampCorners = [ampCorners[(i - nQuarter) % 4] for i in indexes]
        # Construct raw WCS
        ps = dafBase.PropertySet()
        ps.setString("RADESYS", amps[ampSerial][1])
        ps.setDouble("EQUINOX", amps[ampSerial][2])
        ps.setString("CTYPE1", amps[ampSerial][3])
        ps.setString("CTYPE2", amps[ampSerial][4])
        ps.setString("CUNIT1", "deg")
        ps.setString("CUNIT2", "deg")
        ps.setDouble("CRPIX1", amps[ampSerial][5])
        ps.setDouble("CRPIX2", amps[ampSerial][6])
        ps.setDouble("CRVAL1", amps[ampSerial][7])
        ps.setDouble("CRVAL2", amps[ampSerial][8])
        ps.setDouble("CD1_1", amps[ampSerial][9])
        ps.setDouble("CD1_2", amps[ampSerial][10])
        ps.setDouble("CD2_1", amps[ampSerial][11])
        ps.setDouble("CD2_2", amps[ampSerial][12])
        results.append((amp, afwImage.makeWcs(ps), ampCcdCorners, ampCorners))
    return results

def verifySkyTiles(cursor, qsp, kind, wcsMap, wcsList, inButler):
    """Verifies that the sky-tile to CCD mapping in the registry used by inButler
    is correct. Also computes statistics on the angular separation between amp
    corners pre and post WCS determination.
    """
    if kind == 'imsim':
        geomPolicy = cameraGeomUtils.getGeomPolicy(pexPolicy.Policy.createPolicy(
            pexPolicy.DefaultPolicyFile("obs_lsstSim", "Full_STA_geom.paf", "description")))
    else:
        geomPolicy =  cameraGeomUtils.getGeomPolicy(pexPolicy.Policy.createPolicy(
            pexPolicy.DefaultPolicyFile("obs_cfht", "Full_Megacam_geom.paf", "megacam/description")))
    camera = cameraGeomUtils.makeCamera(geomPolicy)
    ampDiskLayout = {}
    for p in geomPolicy.get("CcdDiskLayout").getArray("Amp"):
        ampDiskLayout[p.get("serial")] = (p.get("flipLR"), p.get("flipTB"))
    actualTiles = set(wcsMap.keys())
    processedTiles = set(inButler.queryMetadata("raw", "skyTile"))
    emptyTiles = processedTiles - actualTiles
    missedTiles = actualTiles - processedTiles
    if len(emptyTiles) > 0:
        print dedent("""\
            %d tiles not overlapping any Science CCDs (calexps) were found. This can
            happen if the run didn't span all the input amps in the registry, or
            when a raw amp gets shifted off a tile after WCS determination. So this
            is normal, but note that empty sky-tiles should not contain any
            SourceAssoc output.""" % len(emptyTiles))
    for tile in emptyTiles:
        print "\t%d" % tile
    if len(missedTiles) > 0:
        print dedent("""\
            %d tiles overlapping Science CCDs (calexps) were not processed.
            This indicates the raw amp padding radius used by the registry
            creation script is too small!""" % len(missedTiles))
    for tile in missedTiles:
        print "\t%d" % tile
    print "----"
    missedAmps = set()
    for tile in actualTiles:
        tileGeom = qsp.getGeometry(tile)
        tileAmps = set()
        if kind == "imsim":
            for visit, raft, sensor, channel in inButler.queryMetadata(
                    "raw", "channel", ("visit", "raft", "sensor", "channel"),
                    skyTile=tile, snap=0):
                s1, comma, s2 = sensor
                c1, comma, c2 = channel
                raftNum = lsstSimRafts.index(raft)
                ccdNum = int(s1) * 3 + int(s2)
                ampNum = (int(c2) << 3) + int(c1)
                tileAmps.add((long(visit), raftNum, ccdNum, ampNum))
        else:
            for visit, ccd, amp in inButler.queryMetadata(
                    "raw", "amp", ("visit", "ccd", "amp"), skyTile=tile):
                tileAmps.add((long(visit), 0, int(ccd), int(amp)))
        for ccdData in wcsMap[tile]:
            amps = getAmps(ccdData, camera, ampDiskLayout, cursor, kind)
            for ampData in amps:
               sciCoords = [ccdData[0].pixelToSky(c[0], c[1]) for c in ampData[2]]
               # find amps that should have been included in tile
               ampSerial = ampData[0].getId().getSerial()
               if (tile in missedTiles or not
                   (ccdData[2], ccdData[3], ccdData[4], ampSerial) in tileAmps):
                   if tileGeom.intersects(coordsToPolygon(sciCoords)):
                       print ("Tile %d missing visit %s raft %s ccd %s amp %d" %
                             (tile, str(ccdData[1]), str(ccdData[2]), str(ccdData[3]), ampSerial))
                       missedAmps.add((ccdData[1], ampSerial))
    print "----"
    print "Computing statistics on angular separation between"
    print "science corners and raw corners/edges..."
    nSamples = 0
    maxDist, maxEdgeDist, maxMissEdgeDist = 0.0, 0.0, 0.0
    meanDist, meanEdgeDist, meanMissEdgeDist = 0.0, 0.0, 0.0
    minDist, minEdgeDist, minMissEdgeDist = 180.0, 180.0, 180.0
    for ccdData in wcsList:
        amps = getAmps(ccdData, camera, ampDiskLayout, cursor, kind)
        for ampData in amps:
           sciCoords = [ccdData[0].pixelToSky(c[0], c[1]) for c in ampData[2]]
           rawCoords = [ampData[1].pixelToSky(c[0], c[1]) for c in ampData[3]]
           rawPoly = coordsToPolygon(rawCoords)
           v = rawPoly.getVertices()
           e = rawPoly.getEdges()
           for raw, sci in izip(rawCoords, sciCoords):
              sciPos = tuple(sci.getVector())
              dist = geom.cartesianAngularSep(tuple(raw.getVector()),
                                              tuple(sciPos))
              minDist = min(dist, minDist)
              maxDist = max(dist, maxDist)
              meanDist += dist
              edgeDist = 180.0
              for i in xrange(len(v)):
                  edgeDist = min(edgeDist, geom.minEdgeSep(sciPos, e[i], v[i - 1], v[i]))
              minEdgeDist = min(edgeDist, minEdgeDist)
              maxEdgeDist = max(edgeDist, maxEdgeDist)
              meanEdgeDist += edgeDist
              nSamples += 1
              if (ccdData[1], ampData[0].getId().getSerial()) in missedAmps:
                  minMissEdgeDist = min(edgeDist, minMissEdgeDist)
                  maxMissEdgeDist = max(edgeDist, maxMissEdgeDist)
                  meanMissEdgeDist += edgeDist
    meanDist /= nSamples
    meanEdgeDist /= nSamples
    print "----"
    print "Number of amps examined: %d" % nSamples
    print "Min  raw to science corner distance: %.6f deg" % minDist
    print "Max  raw to science corner distance: %.6f deg" % maxDist
    print "Mean raw to science corner distance: %.6f deg" % meanDist
    print "Min  min science corner to raw amp edge distance: %.6f deg" % minEdgeDist
    print "Max  min science corner to raw amp edge distance: %.6f deg" % maxEdgeDist
    print "Mean min science corner to raw amp edge distance: %.6f deg" % meanEdgeDist
    print "Number of amps missed: %d" % len(missedAmps)
    if len(missedAmps) > 0:
        meanMissEdgeDist /= len(missedAmps)
        print "Min  min science corner to missed raw amp edge distance: %.6f deg" % minMissEdgeDist
        print "Max  min science corner to missed raw amp edge distance: %.6f deg" % maxMissEdgeDist
        print "Mean min science corner to missed raw amp edge distance: %.6f deg" % meanMissEdgeDist

def hostPort(sv):
    hp = sv.split(':')
    if len(hp) > 1:
        return (hp[0], int(hp[1]))
    else:
        return (hp[0], None)

def main():
    # Setup command line options
    usage = dedent("""\
    usage: %prog [options] <kind> <db> <inputRoot>

    Verifies sky-tile to raw amp mapping in registry

    <kind>:       Input dataset, one of 'imsim' or 'cfhtls'
    <db>:         Run database name
    <inputRoot>:  Input root
    """)
    parser = optparse.OptionParser(usage)
    parser.add_option(
        "-u", "--user", dest="user", default="serge",
        help="Database user name to use when connecting to MySQL.")
    parser.add_option(
        "-s", "--server", dest="server", default="lsst10.ncsa.uiuc.edu:3306",
        help="host:port of MySQL server to connect to; defaults to %default")
    parser.add_option(
        "-r", "--registry", dest="registry", default=None,
        help="Input registry for butler")
    opts, args = parser.parse_args()
    if len(args) != 3:
        parser.error("Invalid number of arguments")
    kind = args[0].lower()
    if kind not in ('imsim', 'cfhtls'):
        parser.error("Input dataset must be one of 'imsim' or 'cfhtls'")
    registry = opts.registry or os.path.join(args[2], "registry.sqlite3")
    if kind == 'imsim':
        mapper = LsstSimMapper(root=args[2], calibRoot=None, registry=registry)
    else:
        mapper = CfhtMapper(root=args[2], calibRoot=None, registry=registry)
    inButler = dafPersist.ButlerFactory(mapper=mapper).create()
    passwd = getpass.getpass()
    host, port = hostPort(opts.server)
    qsp = qs.createQuadSpherePixelization()
    conn = sql.connect(host=host, port=port, user=opts.user, passwd=passwd, db=args[1])
    try:
        cursor = conn.cursor()
        try:
            print "Reading WCSes for all Science CCDs in run"
            wcsMap, wcsList = getAllSipWcs(cursor, qsp, kind)
            print "Verifying sky-tile mapping in input registry :"
            print "=============================================="
            verifySkyTiles(cursor, qsp, kind, wcsMap, wcsList, inButler)
        finally:
            cursor.close()
    finally:
        conn.close()


if __name__ == "__main__":
    main()

