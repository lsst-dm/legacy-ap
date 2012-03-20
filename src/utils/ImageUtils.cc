// -*- lsst-c++ -*-

/*
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
 *
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the LSST License Statement and
 * the GNU General Public License along with this program.  If not,
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */

/** @file
  * @brief Image utility method implementations.
  *
  * @ingroup ap
  * @author Serge Monkewitz
  */
#include "lsst/ap/utils/ImageUtils.h"

#include <algorithm>
#include <limits>
#include <utility>

#include "lsst/pex/exceptions.h"
#include "lsst/afw/geom/Point.h"


namespace except = lsst::pex::exceptions;
namespace detection = lsst::afw::detection;
namespace geom = lsst::afw::geom;
namespace image = lsst::afw::image;


namespace lsst { namespace ap { namespace utils {

namespace {

double orient2d(lsst::afw::geom::Point2D const &a,
                lsst::afw::geom::Point2D const &b,
                lsst::afw::geom::Point2D const &c)
{
    double detLeft = (a.getX() - c.getX()) * (b.getY() - c.getY());
    double detRight = (a.getY() - c.getY()) * (b.getX() - c.getX());
    double det = detLeft - detRight;
    double detSum;

    if (detLeft > 0.0) {
        if (detRight <= 0.0) {
            return det;
        } else {
            detSum = detLeft + detRight;
        }
    } else if (detLeft < 0.0) {
        if (detRight >= 0.0) {
            return det;
        } else {
            detSum = -detLeft - detRight;
        }
    } else {
        return det;
    }
    double epsilon = std::numeric_limits<double>::epsilon();
    double errBound = detSum * ((3.0 + 16.0 * epsilon) * epsilon);
    if (det >= errBound || -det >= errBound) {
        return det;
    }
    return 0.0;
}

/** A polygon edge, including pointers to vertices (ordered by y),
  * the corresponding line equation, and an embedded singly linked list.
  */
struct Edge {
    lsst::afw::geom::Point2D const *v1; /**< First (bottom) vertex */
    lsst::afw::geom::Point2D const *v2; /**< Second (top) vertex */
    double dx, dy;
    Edge *next;

    Edge(lsst::afw::geom::Point2D const & p1,
         lsst::afw::geom::Point2D const & p2) :
        v1(p1.getY() < p2.getY() ? &p1 : &p2),
        v2(p1.getY() < p2.getY() ? &p2 : &p1),
        next(0)
    {
        dy = v2->getY() - v1->getY();
        dx = v2->getX() - v1->getX();
    }

    // edges are compared by minimum y
    bool operator==(Edge const & e) const {
        return v1->getY() == e.v1->getY();
    }
    bool operator<(Edge const & e) const {
        return v1->getY() < e.v1->getY();
    }

    /** Returns true if this edge is to the right of e.
      */
    bool rightOf(Edge const &e) const {
        double det;
        if (v1 == e.v1) {
            det = orient2d(*e.v1, *e.v2, *v2);
        } else if (v2 == e.v2) {
            det = orient2d(*e.v1, *e.v2, *v1);
        } else {
            det = orient2d(*e.v1, *e.v2, *v1);
            if (det == 0.0) {
                det = orient2d(*e.v1, *e.v2, *v2);
            }
        }
        return det < 0.0;
    }

    /** Returns x-bounds of this edge within the given y-bounds.
      * Assumes that yMin < yMax, yMin >= v1->getY(), and yMax <= v2->getY().
      */
    std::pair<int, int> const getXBounds(double yMin, double yMax, int width) const {
        double xd1 = v1->getX() + ((yMin - v1->getY()) * dx) / dy;
        double xd2 = v1->getX() + ((yMax - v1->getY()) * dx) / dy;
        int x1 = width, x2 = width;
        if (dx > 0.0) {
            xd1 = xd1 < v1->getX() ? v1->getX() : (xd1 > v2->getX() ? v2->getX() : xd1);
            xd2 = xd2 < xd1 ? xd1 : (xd2 > v2->getX() ? v2->getX() : xd2);
            if (xd1 < image::PixelZeroPos - 0.5 + width) {
                x1 = (xd1 < image::PixelZeroPos - 0.5) ? 0 : image::positionToIndex(xd1);
            }
            if (xd2 < image::PixelZeroPos - 0.5 + width) {
                x2 = (xd2 < image::PixelZeroPos - 0.5) ? 0 : image::positionToIndex(xd2) + 1;
            }
        } else {
            xd2 = xd2 < v2->getX() ? v2->getX() : (xd2 > v1->getX() ? v1->getX() : xd2);
            xd1 = xd1 < xd2 ? xd2 : (xd1 > v1->getX() ? v1->getX() : xd1);
            if (xd2 < image::PixelZeroPos - 0.5 + width) {
                x1 = (xd2 < image::PixelZeroPos - 0.5) ? 0 : image::positionToIndex(xd2);
            }
            if (xd1 < image::PixelZeroPos - 0.5 + width) {
                x2 = (xd1 < image::PixelZeroPos - 0.5) ? 0 : image::positionToIndex(xd1) + 1;
            }
        }
        return std::make_pair(std::min(x1, x2), std::max(x1, x2));
    }

    /** Returns the area to the left of this edge within the specified pixel
        and y bounds.
      */
    double getCoverage(double const yMin, double const yMax, int x) const;
};


double Edge::getCoverage(double const yMin, double const yMax, int x) const {
   double const xp = image::indexToPosition(x);
   double const xMin = xp - 0.5;
   double const xMax = xp + 0.5;
   double cov = 0.0;
   if (dx == 0.0) { // vertical edge
       return (yMax - yMin) * (xMax - v1->getX());
   } else if (dx > 0.0) { // edge with positive slope
       double x1 = xMin;
       // compute intersection with left pixel boundary
       double y1 = v1->getY() + ((xMin - v1->getX()) * dy) / dx;
       if (y1 > yMin) {
           if (y1 >= yMax) {
               return yMax - yMin;
           }
           cov = y1 - yMin;
       } else {
           double x = v1->getX() + ((yMin - v1->getY()) * dx) / dy;
           x1 = (x < xMin) ? xMin : (x > xMax) ? xMax : x;
           y1 = yMin;
       }
       double x2 = xMax;
       // compute intersection with right pixel boundary
       double y2 = v1->getY() + ((xMax - v1->getX()) * dy) / dx;
       if (y2 < yMax) {
           if (y2 <= yMin) {
               return 0.0;
           }
       } else {
           double x = v1->getX() + ((yMax - v1->getY()) * dx) / dy;
           x2 = (x < xMin) ? xMin : (x > xMax) ? xMax : x;
           y2 = yMax;
       }
       if (y1 > y2) {
           y1 = y2;
       }
       if (x1 > x2) {
           x1 = x2;
       }
       cov += (xMax - x2 + 0.5 * (x2 - x1)) * (y2 - y1);
   } else { // edge with negative slope
       double x1 = xMax;
       // compute intersection with right pixel boundary
       double y1 = v1->getY() + ((xMax - v1->getX()) * dy) / dx;
       if (y1 > yMin) {
           if (y1 >= yMax) {
               return 0.0;
           }
       } else {
           double x = v1->getX() + ((yMin - v1->getY()) * dx) / dy;
           x1 = (x < xMin) ? xMin : (x > xMax) ? xMax : x;
           y1 = yMin;
       }
       double x2 = xMin;
       // compute intersection with left pixel boundary
       double y2 = v1->getY() + ((xMin - v1->getX()) * dy) / dx;
       if (y2 < yMax) {
           if (y2 <= yMin) {
               return yMax - yMin;
           }
           cov = yMax - y2;
       } else {
           double x = v1->getX() + ((yMax - v1->getY()) * dx) / dy;
           x2 = (x < xMin) ? xMin : (x > xMax) ? xMax : x;
           y2 = yMax;
       }
       if (y1 > y2) {
           y1 = y2;
       }
       if (x1 < x2) {
           x2 = x1;
       }
       cov += (xMax - x1 + 0.5 * (x2 - x1)) * (y2 - y1);
   }
   return cov;
}

/** A list of edges intersecting a horizontal sweep line, ordered in x.
  */
class SweepLine {
public:
    SweepLine(lsst::afw::image::Image<float>::Ptr img, double startY) :
       _cov(new double[img->getWidth()]),
       _img(img),
       _y(std::max(image::indexToPosition(0) - 0.5, startY)),
       _yIndex(image::positionToIndex(_y)),
       _head(0)
    {
        for (int i = 0; i < _img->getWidth(); ++i) {
            _cov[i] = 0.0;
        }
    }

    /** Adds the given edge to the sweep line.
      */
    void add(Edge &e) {
        if (_yIndex >= _img->getHeight()) {
            return;
        }
        // Insertion sort - in the typical case just 2 edges will
        // be in the sweep line for a given y. More complex data
        // structures are therefore unwarranted.
        Edge *after = 0;
        Edge *before = _head;
        while (before != 0 && e.rightOf(*before)) {
            after = before;
            before = before->next;
        }
        if (after != 0) {
            after->next = &e;
        } else {
            _head = &e;
        }
        e.next = before;
    }

    void advance(double y);
    void finish();

private:
    /** Removes edges with maximum y less than or equal to y
      * from the sweep line.
      */
    void removeBelow(double y) {
        Edge *prev = 0;
        Edge *e = _head;
        while (e != 0) {
            if (e->v2->getY() <= y) {
                // unlink e
                if (prev == 0) {
                    _head = e->next;
                } else {
                    prev->next = e->next;
                }
            } else {
                prev = e;
            }
            e = e->next;
        }
    }

    /** Adds coverage information for a single row of pixels to the given
        coverage map.
      */
    void updateCoverage(int y) {
        if (y < 0 || y >= _img->getHeight()) {
            return;
        }
        typedef image::Image<float>::x_iterator XIter;
        double *cov = _cov.get();
        for (XIter i = _img->row_begin(y), e = _img->row_end(y);
             i != e; ++i, ++cov) {
            double c = *cov;
            *cov = 0.0; // clear coverage information for next row
            if (c > 0.0) {
                *i += std::min(1.0f, static_cast<float>(c));
            }
        }
    }

    void rasterizeCoverage(double yMin, double yMax);

    boost::scoped_array<double> _cov;
    image::Image<float>::Ptr _img;
    double _y;
    int _yIndex;
    // linked edge list
    Edge *_head;
};

void SweepLine::rasterizeCoverage(double yMin, double yMax) {
    if (yMin >= yMax) {
        return;
    }
    // rasterize left/right edge pairs
    int const width = _img->getWidth();
    Edge *left = _head;
    while (left != 0) {
        Edge *right = left->next;
        if (right == 0) {
            break;
        }
        std::pair<int, int> xl = left->getXBounds(yMin, yMax, width);
        std::pair<int, int> xr = right->getXBounds(yMin, yMax, width);
        for (int x = xl.first; x < xl.second; ++x) {
            _cov[x] += left->getCoverage(yMin, yMax, x);
        }
        for (int x = xl.second; x < xr.second; ++x) {
            _cov[x] += yMax - yMin;
        }
        for (int x = xr.first; x < xr.second; ++x) {
            _cov[x] -= right->getCoverage(yMin, yMax, x);
        }
        left = right->next;
    }
}

void SweepLine::advance(double y) {
    if (_yIndex >= _img->getHeight() || _head == 0 || y <= _y) {
        return;
    }
    int const yIndex = image::positionToIndex(y);
    int const height = _img->getHeight();
    for (; _yIndex < yIndex && _yIndex < height; ++_yIndex) {
        double yMax = image::indexToPosition(_yIndex) + 0.5;
        rasterizeCoverage(_y, yMax);
        updateCoverage(_yIndex);
        _y = yMax;
    }
    if (yIndex < height) {
        rasterizeCoverage(_y, y);
        removeBelow(y);
        _y = y;
    } else {
        _head = 0;
    }
}

void SweepLine::finish() {
    double yMax = image::indexToPosition(_img->getHeight() - 1) + 0.5;
    while (_head != 0) {
        Edge *e = _head;
        double y = yMax;
        do {
            y = std::min(y, e->v2->getY());
            e = e->next;
        } while (e != 0);
        advance(y);
    }
    // flush coverage data to map and discard all edges
    updateCoverage(_yIndex);
    _head = 0;
}

} // namespace


/** Fills in a histogram of positions for the given sources. The value of each
  * histogram pixel is set to the number of sources falling within that pixel.
  * Sources falling outside the histogram image are either ignored or cause an
  * exception to be raised.
  *
  * @param[in,out] histogram  Histogram image to update.
  * @param[in] sources        Sources to generate a histogram for
  * @param[in] wcs            WCS of histogram image
  * @param[in] ignoreOffImage If true ignore off image sources, otherwise raise
  *                           an exception.
  */
void makeSourceHistogram(
    lsst::afw::image::Image<unsigned short>::Ptr histogram,
    lsst::afw::detection::SourceSet const & sources,
    lsst::afw::image::Wcs::Ptr wcs,
    bool ignoreOffImage)
{
    typedef detection::SourceSet::const_iterator SourceIter;
    for (SourceIter i = sources.begin(), e = sources.end(); i != e; ++i) {
        geom::Point2D xy = wcs->skyToPixel(*((*i)->getRaDec()));
       int x = histogram->positionToIndex(xy[0], image::X).first;
       int y = histogram->positionToIndex(xy[1], image::Y).first;
       if (x < 0 || x >= histogram->getWidth() ||
           y < 0 || y >= histogram->getHeight()) {
           if (!ignoreOffImage) {
               throw LSST_EXCEPT(except::RuntimeErrorException,
                                 "input SourceSet contains sources lying "
                                 "outside the histogram image");
           }
           continue;
       }
       histogram->operator()(x, y) += 1;
    }
}

/** Rasterizes the (not necessarily convex) polygon obtained by connecting the
  * the given vertices. For each pixel in the output image, the fraction of the
  * pixel area overlapping the input polygon is stored.
  *
  * @param[in,out] img  Image to rasterize to.
  * @param[in] verts    Polygon vertices
  */
void rasterizePolygon(
    std::vector<lsst::afw::geom::Point2D> const &verts,
    lsst::afw::image::Image<float>::Ptr img)
{
    typedef std::vector<geom::Point2D>::const_iterator VertexIter;
    typedef std::vector<Edge>::iterator EdgeIter;

    if (img->getWidth() <= 0 || img->getHeight() <= 0) {
        throw new LSST_EXCEPT(except::InvalidParameterException,
                              "image width/height must be at least 1");
    }
    if (verts.size() < 3) {
        throw new LSST_EXCEPT(except::InvalidParameterException,
                              "polygon must have at least 3 vertices");
    }
    std::vector<Edge> edges;
    edges.reserve(verts.size());
    double minY = image::indexToPosition(0) - 0.5;
    double maxY = image::indexToPosition(img->getHeight()) - 0.5;
    for (VertexIter v2 = verts.begin(), e = verts.end(), v1 = e - 1;
         v2 != e; ++v2) {
        VertexIter v = v1;
        v1 = v2;
        if (v->getY() == v2->getY()) {
            continue; // omit horizontal and degenerate edges
        } else if (v->getY() < v2->getY()) {
            if (v->getY() >= maxY || v2->getY() < minY) {
                continue; // edge not in y-range
            }
        } else {
            if (v2->getY() >= maxY || v->getY() < minY) {
                continue; // edge not in y-range
            }
        }
        edges.push_back(Edge(*v, *v2));
    }
    if (edges.size() > 0) {
        std::sort(edges.begin(), edges.end());
        SweepLine sweeper(img, edges.front().v1->getY());
        EdgeIter i = edges.begin(), e = edges.end();
        // add all edges intersecting minY to sweep line
        for (; i != e && i->v1->getY() <= minY; ++i) {
            sweeper.add(*i);
        }
        while (i != e) {
            double y = i->v1->getY();
            sweeper.advance(y);
            sweeper.add(*i);
            for (++i; i != e && i->v1->getY() == y; ++i) {
                sweeper.add(*i);
            }
        }
        sweeper.finish();
    }
}

/** Updates a coverage map to include contributions from the given image
  * (specified by a WCS and a pair of dimensions).
  *
  * @param[in,out] covMap   Coverage map to update.
  * @param[in] covMapWcs    WCS of coverage map.
  * @param[in] wcs          WCS of image to rasterize.
  * @param[in] width        Width of image to rasterize.
  * @param[in] height       Height of image to rasterize.
  * @param[in] step         If <= 0, the input image is rasterized as a
  *                         polygon formed by connecting the positions of
  *                         the 4 image corners with straight lines in
  *                         coverage-map pixel space. Otherwise, the boundary
  *                         of the input image is rasterize by connecting
  *                         step boundary pixels at a time (corners are
  *                         always included).
  */
void updateCoverageMap(
    lsst::afw::image::Image<float>::Ptr covMap,
    lsst::afw::image::Wcs::Ptr covMapWcs,
    lsst::afw::image::Wcs::Ptr wcs,
    int width,
    int height,
    int step)
{
    typedef std::vector<geom::Point2D>::iterator VertexIter;
    if (width <= 0 || height <= 0) {
        throw new LSST_EXCEPT(except::InvalidParameterException,
                              "Width/height of input image must be positive");
    }
    std::vector<geom::Point2D> v;
    if (step <= 0) {
        v.reserve(4);
        v.push_back(geom::Point2D(-0.5, -0.5));
        v.push_back(geom::Point2D(width - 0.5, -0.5));
        v.push_back(geom::Point2D(width - 0.5, height - 0.5));
        v.push_back(geom::Point2D(-0.5, height - 0.5));
    } else {
        v.reserve(2 * width + 2 * height + 2);
        for (int i = 0; i < width; i += step) {
            v.push_back(geom::Point2D(i - 0.5, -0.5));
        }
        for (int i = 0; i < height; i += step) {
            v.push_back(geom::Point2D(width - 0.5, i - 0.5));
        }
        for (int i = width; i > 0; i -= step) {
            v.push_back(geom::Point2D(i - 0.5, height - 0.5));
        }
        for (int i = height; i > 0; i -= step) {
            v.push_back(geom::Point2D(-0.5, i - 0.5));
        }
    }
    // Convert from image pixel space to coverage map pixel space.
    for (VertexIter i = v.begin(), e = v.end(); i != e; ++i) {
        *i = covMapWcs->skyToPixel(*wcs->pixelToSky(*i));
    }
    // rasterize polygon
    rasterizePolygon(v, covMap);
}

}}} // namespace lsst::ap::utils

