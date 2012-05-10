#! /usr/bin/env python

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

from lsst.afw.image import Filter

__all__ = ["makeMapper", "getFilterNames"]

_mapper = {
    "lsstsim": "lsst.obs.lsstSim.LsstSimMapper",
    "hscsim": "lsst.obs.hscSim.hscSimMapper.HscSimMapper",
    "suprimecam": "lsst.obs.suprimecam.suprimecamMapper.SuprimecamMapper",
    "test": "lsst.obs.test.TestMapper",
    "sdss": "lsst.obs.sdss.SdssMapper",
    "cfht": "lsst.obs.cfht.CfhtMapper",
}

def makeMapper(camera):
    """Return an instance of the lsst.daf.persistence.Mapper class
    to use for the camera of the given (case-insensitive) name.
    """
    camera = camera.lower()
    if camera not in _mapper:
        raise RuntimeError(str.format("{} is not a valid camera name", camera))
    name = _mapper[camera]
    try:
        pieces = name.split('.')
        cls = reduce(getattr, pieces[1:], __import__('.'.join(pieces[:-1])))
        mapper = cls()
    except:
        raise RuntimeError(str.format("Failed to construct a {} mapper", name))
    return mapper

def getFilterNames():
    """Return a list of filter names in increasing ID order. This assumes
    that an lsst.daf.persistence.Mapper which sets up filter definitions
    has been created, or that filters have been manually defined with e.g.
    lsst.afw.image.utils.defineFilter().
    """
    names = list(Filter.getNames())
    names.sort(key=lambda name: Filter(name, False).getId())
    return names

