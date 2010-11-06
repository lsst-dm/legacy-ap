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
import csv
import getpass
import MySQLdb as sql
import optparse
import pdb
import math
import numpy
# For now, use pywcs since lsst.afw.image.TanWcs doesn't expose a way
# to retrieve just SIP distortions
import pywcs
from textwrap import dedent

def maxDistort(transform, width, height, gridSize):
    a = []
    for x in xrange(0, width, gridSize):
        for y in xrange(0, height, gridSize):
            a.append([x, y])
        a.append([x, height])
    for y in xrange(0, height, gridSize):
        a.append([width, y])
    a.append([width, height])
    a = numpy.array(a, dtype=numpy.float64)
    a += 0.5
    deltas = transform.pix2foc(a, 1) - a
    maxd = 0.0
    for x, y in deltas:
        maxd = max(maxd, x*x + y*y)
    return math.sqrt(maxd)

def findBigSip(cursor, kind, gridSize, thresh, outFile):
    """Constructs a Wcs object from each entry in the Science_Ccd_Exposure
    table. SIP distortion parameters are read in from
    Science_Ccd_Exposure_Metadata. Returns a 3-tuple with the following
    contents:

    - a mapping from sky-tiles to a list of (WCS, filter) tuples for
      overlapping science CCDs
    - a list of all (WCS, filter) tuples read in
    - a bounding circle (lsst.geom.SphericalCircle) for all science CCDs.
    """
    offset = 0
    n = 0
    blocksize = 1000
    maxd = 0.0
    if kind == 'imsim':
        width, height = 4072, 4000
    else:
        width, height = 2048, 4610
    writer = csv.writer(open(outFile, "wb"))
    writer.writerow(["scienceCcdExposureId", "visit", "raft", "ccd", "filterId", "maxDistort"])
    # An alternative would be to create WCSes from FITS metadata or FITS
    # files themselves, but these aren't always available.
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
            if row[7].endswith("-SIP") and row[8].endswith("-SIP"):
                transform = pywcs.WCS()
                transform.wcs.ctype[0] = row[7]
                transform.wcs.ctype[1] = row[8]
                transform.wcs.crpix = [float(row[9]), float(row[10])]
                transform.wcs.crval = [float(row[11]), float(row[12])]
                cd = numpy.zeros(shape=(2,2), dtype=numpy.float64)
                cd[0,0] = float(row[13])
                cd[0,1] = float(row[14])
                cd[1,0] = float(row[15])
                cd[1,1] = float(row[16])
                transform.wcs.cd = cd
                cursor.execute("""
                    SELECT metadataKey, intValue
                    FROM Science_Ccd_Exposure_Metadata
                    WHERE scienceCcdExposureId = %d AND
                          intValue IS NOT NULL AND
                          metadataKey RLIKE '^[AB]P?_ORDER$'
                    """ % row[0])
                for k, v in cursor.fetchall():
                    if k == "A_ORDER":
                        a = numpy.zeros(shape=(v + 1, v + 1), dtype=numpy.float64)
                    elif k == "B_ORDER":
                        b = numpy.zeros(shape=(v + 1, v + 1), dtype=numpy.float64)
                    elif k == "AP_ORDER":
                        ap = numpy.zeros(shape=(v + 1, v + 1), dtype=numpy.float64)
                    elif k == "BP_ORDER":
                        bp = numpy.zeros(shape=(v + 1, v + 1), dtype=numpy.float64)
                cursor.execute("""
                    SELECT metadataKey, doubleValue
                    FROM Science_Ccd_Exposure_Metadata
                    WHERE scienceCcdExposureId = %d AND
                          doubleValue IS NOT NULL AND
                          metadataKey RLIKE '^[AB]P?_[0-9]+_[0-9]+$'
                    """ % row[0])
                for k, v in cursor.fetchall():
                    if k.startswith("A_"):    a[int(k[2]), int(k[4])] = v
                    elif k.startswith("B_"):  b[int(k[2]), int(k[4])] = v
                    elif k.startswith("AP_"): ap[int(k[3]), int(k[5])] = v
                    elif k.startswith("BP_"): bp[int(k[3]), int(k[5])] = v
                # make sure the SIP solution doesn't contain a translation
                assert a[0,0] == 0.0 and b[0,0] == 0.0
                transform.sip = pywcs.Sip(a, b, ap, bp, transform.wcs.crpix)
                d = maxDistort(transform, width, height, gridSize)
                maxd = max(maxd, d)
                if d >= thresh:
                    writer.writerow(map(str, row[0:5]) + [str(d)])
        offset += blocksize
    print "Maximum distortion: %.6f pixels" % maxd

def hostPort(sv):
    hp = sv.split(':')
    if len(hp) > 1:
        return (hp[0], int(hp[1]))
    else:
        return (hp[0], None)

def main():
    usage = dedent("""\
    usage: %prog [options] <kind> <db> <outfile>

    Prints a list of CCDs with big SIP distortions.

    <kind>:           Input dataset, one of 'imsim' or 'cfhtls'
    <db>:             Run database name
    <outputfile>:     Name of file to write output CSV to
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
        "-g", "--grid-size", type="int", dest="gridSize", default=64,
        help=dedent("""\
            Number of pixels in x,y separating adjacent test positions in grid
            used to esimate maximum SIP distortion"""))
    parser.add_option(
        "-t", "--threshold", type="float", dest="threshold", default=100.0,
        help=dedent("""\
            If the SIP distortion for a test position exceeds this threshold,
            it is reported by this program (default %default)"""))
    opts, args = parser.parse_args()
    if len(args) != 3:
        parser.error("Invalid number of arguments")
    kind = args[0].lower()
    if kind not in ('imsim', 'cfhtls'):
        parser.error("Input kind must be one of 'imsim' or 'cfhtls'")
    passwd = getpass.getpass()
    host, port = hostPort(opts.server)
    conn = sql.connect(host=host, port=port, user=opts.user, passwd=passwd, db=args[1])
    try:
        cursor = conn.cursor()
        try:
            findBigSip(cursor, kind, opts.gridSize, opts.threshold, args[2])
        finally:
            cursor.close()
    finally:
        conn.close()


if __name__ == "__main__":
    main()

