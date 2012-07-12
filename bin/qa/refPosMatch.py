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
import argparse

import lsst.ap.match as apMatch
from lsst.pipe.base import ConfigFileAction, ConfigValueAction

__all__ = ["referenceMatch"]

def referenceMatch(config, outputFile, referenceCatalog, positionCatalog):
    apMatch.referenceMatch(
        referenceCatalog,
        config.ref.makeControl(),
        config.refDialect.makeControl(),
        positionCatalog,
        config.pos.makeControl(),
        config.posDialect.makeControl(),
        outputFile,
        config.outDialect.makeControl(),
        config.getRadius(),
        config.getParallaxThresh(),
        config.doOutputRefExtras,
        True)

def main():
    parser = argparse.ArgumentParser(description=
        "Matches a reference catalog against a table of positions. "
        "Both the reference catalog and position table must be "
        "in increasing declination order. Matching is controlled by "
        "lsst.ap.match.ReferenceMatchConfig - supply overrides with -c or -C.")
    parser.add_argument("-c", "--config", nargs="*", action=ConfigValueAction,
        help="config override(s), e.g. -c foo=qux bar.baz=3",
        metavar="NAME=VALUE")
    parser.add_argument("-C", "--config-file", nargs="*",
        action=ConfigFileAction, help="config override file(s)")
    parser.add_argument("--camera", dest="camera", default="lsstSim",
        help="Name of desired camera (defaults to %(default)s)")
    parser.add_argument("outputFile", help="Name of output CSV file")
    parser.add_argument("referenceCatalog",
        help="Name of reference catalog CSV file")
    parser.add_argument("positionCatalog",
        help="Name of position (source/object) catalog CSV file")
    ns = argparse.Namespace()
    ns.config = apMatch.ReferenceMatchConfig()
    ns = parser.parse_args(namespace=ns)
    referenceMatch(ns.config,
                   ns.outputFile,
                   ns.referenceCatalog,
                   ns.positionCatalog)

if __name__ == "__main__":
    main()

