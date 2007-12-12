#! /usr/bin/env python
"""
Convenience script that subdivides the sky into stripes and calls
dumpObjects to create chunk files for each one 
"""

import os
import os.path
import optparse
import subprocess
import sys

import lsst.ap.interface as ap


def main():
    defZpd   = 180
    defZps   = 63
    defDbUrl = 'mysql://lsst10.ncsa.uiuc.edu:3306/DC2'
    parser   = optparse.OptionParser()
    parser.add_option('-z', '--zones_per_degree', type=int,
                      default=defZpd, help='number of zones per degree (default: %d)' % defZpd)
    parser.add_option('-s', '--zones_per_stripe', type=int,
                      default=defZps, help='number of zones per stripe (default: %d)' % defZps)
    parser.add_option('-u', '--db_url',
                      default=defDbUrl,
                      help='Connection URL for database (default: %s)' % defDbUrl)
    parser.add_option('-d', '--directory', help='Directory to create stripe directories in')
    (options, args) = parser.parse_args()

    # create a decomposition of the unit sphere
    zsc = ap.ZoneStripeChunkDecomposition(options.zones_per_degree, options.zones_per_stripe, 1)
    dir = options.directory
    if not os.path.isdir(dir):
        print "root directory %s doesn't exist" % dir
        sys.exit(1)

    # create stripe directories, run dumpObjects for each stripe
    for i in xrange(zsc.decToStripe(-90), zsc.decToStripe(90) + 1):
        stripeDir = os.path.join(dir, 'stripe_%d' % i)
        if not os.path.isdir(stripeDir):
            os.mkdir(stripeDir)
        args = ['dumpObjects',
                '-u', options.db_url,
                '-o', os.path.join(stripeDir, 'objref_chunk'),
                '-d', str(zsc.getStripeDecMin(i)),
                '-D', str(zsc.getStripeDecMax(i)),
                '-n', str(zsc.getNumChunksPerStripe(i))]
        exitcode = subprocess.Popen(args).wait()
        if exitcode != 0:
            print 'dumpObjects failed!'
            sys.exit(exitcode)


if __name__ == '__main__':
    main()

