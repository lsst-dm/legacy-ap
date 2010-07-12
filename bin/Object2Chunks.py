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

"""
Convenience script that subdivides the sky into stripes and calls
dumpObjects to create chunk files for each one 
"""

import copy
import os, os.path
import optparse
import subprocess
import sys

import lsst.ap as ap


def main():
    defZpd = 180
    defZps = 63
    defHost = 'lsst10.ncsa.uiuc.edu'
    defPort = '3306'
    defDatabase = 'DC3a_catalogs'
    defTable = 'CFHTLSObject'
    parser = optparse.OptionParser()
    parser.add_option('-z', '--zones_per_degree', type=int, default=defZpd,
                      help='number of zones per degree (default: %d)' % defZpd)
    parser.add_option('-s', '--zones_per_stripe', type=int, default=defZps,
                      help='number of zones per stripe (default: %d)' % defZps)
    parser.add_option('-d', '--directory', help='Directory to create stripe directories in')
    parser.add_option('-H', '--host', default=defHost,
                      help='MySQL database server hostname (default: %s)' % defHost)
    parser.add_option('-P', '--port', type=int, default=defPort,
                      help='MySQL database server port number (default: %s)' % defPort)
    parser.add_option('-b', '--database', default=defDatabase,
                      help='database name (default: %s)' % defDatabase)
    parser.add_option('-t', '--table', default=defTable,
                      help='object table name (default: %s)' % defTable)
    parser.add_option('-p', '--policy',
                      help='filename of policy for database authentication')

    (options, args) = parser.parse_args()
    if not options.directory:
        print "Output directory must be specified"
        sys.exit(1)

    # create a decomposition of the unit sphere
    zsc = ap.ZoneStripeChunkDecomposition(options.zones_per_degree, options.zones_per_stripe, 1)
    if not os.path.isdir(options.directory):
        print "root directory %s doesn't exist" % options.directory
        sys.exit(1)

    fixedArgs = ['dumpObjects',
                 '-H', options.host,
                 '-P', str(options.port),
                 '-b', options.database,
                 '-t', options.table ]
    if options.policy:
        fixedArgs.extend(['-p', options.policy])

    # create stripe directories, run dumpObjects for each stripe
    for i in xrange(zsc.decToStripe(-90), zsc.decToStripe(90) + 1):
        stripeDir = os.path.join(options.directory, str(i))
        if not os.path.isdir(stripeDir):
            os.mkdir(stripeDir)
        args = copy.copy(fixedArgs)
        args.extend(['-o', os.path.join(stripeDir, 'ref'),
                     '-d', str(zsc.getStripeDecMin(i)),
                     '-D', str(zsc.getStripeDecMax(i)),
                     '-n', str(zsc.getNumChunksPerStripe(i))])
        exitcode = subprocess.Popen(args).wait()
        if exitcode != 0:
            print 'dumpObjects failed!'
            sys.exit(exitcode)


if __name__ == '__main__':
    main()

