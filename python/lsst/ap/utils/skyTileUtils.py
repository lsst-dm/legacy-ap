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

import math

import lsst.geom as geom
import lsst.skypix as skypix
import lsst.afw.geom as afwGeom
import lsst.afw.coord as afwCoord
import lsst.afw.image as afwImage

__all__ = ["findWcsCoveringSkyTile", "createImageCoveringSkyTile"]

def findWcsCoveringSkyTile(skyPixelization, skyTileId, imageRes):
    """Computes and returns a TAN WCS such that a 2D image with the
    given WCS and the following properties completely covers the
    sky-tile with the given pixel id:

    - NAXIS1/NAXIS2 >= imageRes
    - CRPIX1 = NAXIS1 / 2 + 0.5
    - CRPIX2 = NAXIS2 / 2 + 0.5
    """
    if not isinstance(imageRes, (int, long)):
        raise TypeError("Image resolution must be an integer")
    if imageRes < 1:
        raise RuntimeError("Image resolution must be at least 1")
    crpix = afwGeom.Point2D(0.5*(imageRes + 1), 0.5*(imageRes + 1))
    crval = geom.sphericalCoords(skyPixelization.getCenter(skyTileId))
    crval = afwCoord.makeCoord(afwCoord.ICRS, crval[0] * afwGeom.degrees, crval[1] * afwGeom.degrees)
    skyTile = skyPixelization.getGeometry(skyTileId)
    # Start with a huge TAN image centered at the sky-tile center,
    # then shrink it using binary search to determine suitable
    # CD matrix coefficients
    scale = 1000.0  # deg/pixel, ridiculously large
    delta = 0.5*scale
    frac = 0.01 # desired relative accuracy of CD matrix coeffs
    wcs = afwImage.makeWcs(crval, crpix, scale, 0.0, 0.0, scale)
    imagePoly = skypix.imageToPolygon(wcs, imageRes, imageRes)
    # Make sure the initial guess really is too large
    if not imagePoly.contains(skyTile):
        raise RuntimeError("Failed to construct image WCS covering sky-tile")
    # Search for a WCS with a tight fit to the sky-tile. Note that the
    # tightness of fit could be further improved by searching for a rotation
    # and not just a pixel scale.
    while delta >= frac * scale:
        tmp = scale - delta
        wcs = afwImage.makeWcs(crval, crpix, tmp, 0.0, 0.0, tmp)
        imagePoly = skypix.imageToPolygon(wcs, imageRes, imageRes)
        delta *= 0.5
        if imagePoly.contains(skyTile):
            scale = tmp
    return afwImage.makeWcs(crval, crpix, scale, 0.0, 0.0, scale)

def createImageCoveringSkyTile(skyPixelization, skyTileId, imageRes,
                               imageType=afwImage.DecoratedImageU):
    """Creates an image of the specified type and resolution along with
    a WCS that covers the given sky tile. If the requested image type
    is an lsst.afw.image.ExposureX, the exposure WCS is set automatically.
    Otherwise, the WCS is translated to a lsst.daf.base.PropertySet and set
    as lsst.afw.image.[Decorated]ImageX metadata.
    """
    wcs = findWcsCoveringSkyTile(skyPixelization, skyTileId, imageRes)
    img = imageType(afwGeom.Extent2I(imageRes, imageRes))
    if hasattr(img, 'setWcs'):
        img.setWcs(wcs)
    elif hasattr(img, 'setMetadata'):
        img.setMetadata(wcs.getFitsMetadata())
    return img, wcs

