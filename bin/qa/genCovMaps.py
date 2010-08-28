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
import subprocess
from textwrap import dedent

import lsst.daf.base as dafBase
import lsst.daf.persistence as dafPersist
import lsst.pex.policy as pexPolicy
import lsst.geom as geom
import lsst.skypix.quadsphere as qs
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.ap.cluster as apCluster

def getAllSipWcs(cursor, qsp, kind):
    """Constructs a Wcs object from each entry in the Science_Ccd_Exposure
    table. SIP distortion parameters are read in from
    Science_Ccd_Exposure_Metadata. Returns a 3-tuple with the following
    contents:

    - a mapping from sky-tiles to a list of (WCS, filter) tuples for
      overlapping science CCDs
    - a list of all (WCS, filter) tuples read in
    - a bounding circle (lsst.geom.SphericalCircle) for all science CCDs.
    """
    wcsMap = {}
    wcsList = []
    offset = 0
    centers = []
    radius = 0.0
    n = 0
    blocksize = 1000
    # An alternative would be to create WCSes from FITS metadata or FITS
    # files themselves, but these aren't always available.
    while True:
        cursor.execute(
            """SELECT scienceCcdExposureId, filterId,
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
            ps.setString("RADESYS", row[2])
            ps.setDouble("EQUINOX", row[3])
            ps.setString("CTYPE1", row[4])
            ps.setString("CTYPE2", row[5])
            ps.setString("CUNIT1", "deg")
            ps.setString("CUNIT2", "deg")
            ps.setDouble("CRPIX1", row[6])
            ps.setDouble("CRPIX2", row[7])
            ps.setDouble("CRVAL1", row[8])
            ps.setDouble("CRVAL2", row[9])
            ps.setDouble("CD1_1", row[10])
            ps.setDouble("CD1_2", row[11])
            ps.setDouble("CD2_1", row[12])
            ps.setDouble("CD2_2", row[13])
            if row[4].endswith("-SIP") and row[5].endswith("-SIP"):
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
            wcsList.append([wcs, row[1]])
            if kind == 'imsim':
                poly = qs.imageToPolygon(wcs, 4072.0, 4000.0, 0.0)
            else:
                poly = qs.imageToPolygon(wcs, 2048.0, 4610.0, 0.0)
            bc = poly.getBoundingCircle()
            radius = max(radius, bc.getRadius())
            centers.append(geom.cartesianUnitVector(bc.getCenter()))
            pix = qsp.intersect(poly)
            for skyTile in pix:
                if skyTile in wcsMap:
                    wcsMap[skyTile].append([wcs, row[1]])
                else:
                    wcsMap[skyTile] = [[wcs, row[1]]]
        n = offset + len(rows)
        print "... %d" % n
        offset += blocksize
    x, y, z = 0.0, 0.0, 0.0
    for v in centers:
        x += v[0]; y += v[1]; z += v[2]
    cen = geom.normalize((x, y, z))
    r = max(geom.cartesianAngularSep(cen, v) for v in centers) + 2.0*radius
    return wcsMap, wcsList, geom.SphericalCircle(cen, r)


def persistCovMaps(covMaps, outputDir):
    """Persists coverage maps for up to 6 LSST filters (ugrizy) and
    produces an all-filter coverage map by summing the individual maps.
    Output files are automatically gzipped (they tend to compress well).
    """
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
    persistence = dafPersist.Persistence.getPersistence(pexPolicy.Policy())
    for covMap, filter in izip(covMaps, "ugrizy"):
        # Write out coverage map
        outPath = os.path.join(outputDir, '%sCovMap.fits' % filter)
        storageList = dafPersist.StorageList()
        storageList.append(persistence.getPersistStorage("FitsStorage",
            dafPersist.LogicalLocation(outPath)))
        persistence.persist(covMap, storageList, dafBase.PropertySet())
        subprocess.check_call(['gzip', '-f', outPath])
    allCovMap = covMaps.pop()
    allCovMapImg = allCovMap.getImage()
    for covMap in covMaps:
        allCovMapImg += covMap.getImage()
    outPath = os.path.join(outputDir, 'allCovMap.fits')
    storageList = dafPersist.StorageList()
    storageList.append(persistence.getPersistStorage("FitsStorage",
        dafPersist.LogicalLocation(outPath)))
    persistence.persist(allCovMap, storageList, dafBase.PropertySet())
    subprocess.check_call(['gzip', '-f', outPath])


def processSkyTile(wcsMap, outputDir, kind, skyTile, qsp, res):
    """Creates and persists coverage maps for a particular sky tile.
    """
    if skyTile not in wcsMap:
        return
    covMaps = []
    for i in xrange(6):
        covMap, tileWcs = apCluster.utils.createImageCoveringSkyTile(
            qsp, skyTile, res, afwImage.DecoratedImageF)
        covMaps.append(covMap)
    if kind == 'imsim':
        width, height = 4072, 4000
    else:
        width, height = 2048, 4610
    for wcs, filter in wcsMap[skyTile]:
        assert filter >= 0 and filter < 6
        apCluster.updateCoverageMap(covMaps[filter].getImage(), tileWcs,
                                    wcs, width, height, False)
    persistCovMaps(covMaps, os.path.join(outputDir, "st%d" % skyTile))

def processRun(wcsList, outputDir, kind, bc, res):
    """Creates and persists coverage maps for an entire run.
    """
    crpix = afwGeom.makePointD(0.5*(res + 1), 0.5*(res + 1))
    crval = afwGeom.makePointD(bc.getCenter()[0], bc.getCenter()[1])
    scale = bc.getRadius()*2.0/res
    runWcs = afwImage.createWcs(crval, crpix, scale, 0.0, 0.0, scale)
    covMaps = []
    for i in xrange(6):
        runImg = afwImage.DecoratedImageF(res, res)
        if hasattr(runImg, 'setWcs'):
            runImg.setWcs(runWcs)
        elif hasattr(runImg, 'setMetadata'):
            runImg.setMetadata(runWcs.getFitsMetadata())
        covMaps.append(runImg)
    if kind == 'imsim':
        width, height = 4072, 4000
    else:
        width, height = 2048, 4610
    for wcs, filter in wcsList:
        assert filter >= 0 and filter < 6
        apCluster.updateCoverageMap(covMaps[filter].getImage(), runWcs,
                                    wcs, width, height, True)
    persistCovMaps(covMaps, outputDir)

def hostPort(sv):
    hp = sv.split(':')
    if len(hp) > 1:
        return (hp[0], int(hp[1]))
    else:
        return (hp[0], None)

def main():
    usage = dedent("""\
    usage: %prog [options] <kind> <db> <output_dir>

    Generates coverage maps for an LSST pipeline run.

    <kind>:           Input dataset, one of 'imsim' or 'cfhtls'
    <db>:             Run database name
    <output_dir>:     Output directory.
    """)
    parser = optparse.OptionParser(usage)
    parser.add_option(
        "-d", "--database", dest="db", default="serge",
        help="Run database to extract image metadata from")
    parser.add_option(
        "-u", "--user", dest="user", default="serge",
        help="Database user name to use when connecting to MySQL.")
    parser.add_option(
        "-s", "--server", dest="server", default="lsst10.ncsa.uiuc.edu:3306",
        help="host:port of MySQL server to connect to; defaults to %default")
    parser.add_option(
        "-r", "--tile-res", type="int", dest="tileRes", default=2000,
        help=dedent("""\
            Tile coverage map resolution (default %default). Values less than 1
            will cause per-tile coverage map generation to be skipped"""))
    parser.add_option(
        "-R", "--run-res", type="int", dest="runRes", default=4000,
        help=dedent("""\
            Run coverage map resolution (default %default). Values less than 1
            will cause run coverage map generation to be skipped"""))
    opts, args = parser.parse_args()
    if len(args) != 3:
        parser.error("Invalid number of arguments")
    kind = args[0].lower()
    if kind not in ('imsim', 'cfhtls'):
        parser.error("Input kind must be one of 'imsim' or 'cfhtls'")
    if opts.tileRes <= 0 and opts.runRes <= 0:
        parser.error("Skipping both run and per sky-tile coverage map generation")
    # read in WCSes for all Science CCDs in run
    passwd = getpass.getpass()
    host, port = hostPort(opts.server)
    qsp = qs.createQuadSpherePixelization()
    conn = sql.connect(host=host, port=port, user=opts.user, passwd=passwd, db=args[1])
    try:
        cursor = conn.cursor()
        try:
            print "Reading WCSes for all Science CCDs in run"
            wcsMap, wcsList, bc = getAllSipWcs(cursor, qsp, kind)
        finally:
            cursor.close()
    finally:
        conn.close()
    # Build coverage maps
    print "Building run overview coverage maps"
    if opts.runRes > 0:
        processRun(wcsList, args[2], kind, bc, opts.runRes)
    if opts.tileRes > 0:
        for skyTile in wcsMap.keys():
            print "Processing sky-tile %d" % skyTile
            processSkyTile(wcsMap, args[2], kind, skyTile, qsp, opts.tileRes)


if __name__ == "__main__":
    main()

