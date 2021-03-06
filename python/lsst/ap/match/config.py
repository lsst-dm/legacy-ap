#
# LSST Data Management System
# Copyright 2012 LSST Corporation.
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

import lsst.pex.config as pexConfig
import lsst.afw.geom as afwGeom

from lsst.ap.utils import CsvConfig
from .matchLib import CatalogControl

@pexConfig.wrap(CatalogControl)
class CatalogConfig(pexConfig.Config):
    pass

class ReferenceMatchConfig(pexConfig.Config):
    ref = pexConfig.ConfigField(
        dtype=CatalogConfig,
        doc="Columns and properties of the reference catalog")

    refDialect = pexConfig.ConfigField(
        dtype=CsvConfig,
        doc="CSV format of the reference catalog file (delimiter, quoting, etc...)")

    pos = pexConfig.ConfigField(
        dtype=CatalogConfig,
        doc="Columns and properties of the position catalog")

    posDialect = pexConfig.ConfigField(
        dtype=CsvConfig,
        doc="CSV format of the position catalog file (delimiter, quoting, etc...)")

    expDialect = pexConfig.ConfigField(
        dtype=CsvConfig,
        doc="CSV format of the exposure metadata file (delimiter, quoting, etc...)")

    outDialect = pexConfig.ConfigField(
        dtype=CsvConfig,
        doc="CSV format of the output file (delimiter, quoting, etc...)")

    radius = pexConfig.RangeField(
        dtype=float,
        doc="Match radius (arcsec)",
        default=2.0,
        min=0.0,
        inclusiveMin=False)

    parallaxThresh = pexConfig.RangeField(
        dtype=float,
        doc="Parallax threshold (milliarcsec). Positions of reference catalog "
            "entries with parallax below this value will not be subject to "
            "reduction from barycentric to geocentric place.",
        default=10.0,
        min=0.0)

    doOutputRefExtras = pexConfig.Field(
        dtype=bool,
        doc="If set, match records will contain the proper-motion/parallax "
            "corrected position of the reference object at the epoch of the "
            "matched position, as well as flags indicating which corrections "
            "were performed.",
        default=True)

    expIdKey = pexConfig.Field(
        dtype=str,
        default="Computed_ccdExposureId",
        doc="Name of metadata key corresponding to the exposure ID")

    def getRadius(self):
        """Return the radius parameter value as an lsst.afw.geom.Angle"""
        return self.radius * afwGeom.arcseconds

    def getParallaxThresh(self):
        """Return the self.parallaxThresh parameter value as an lsst.afw.geom.Angle""" 
        return self.parallaxThresh/1000.0 * afwGeom.arcseconds

