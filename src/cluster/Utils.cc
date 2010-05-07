// -*- lsst-c++ -*-
/** @file
  * @brief Utility method implementations.
  *
  * @ingroup ap
  * @author Serge Monkewitz
  */
#include "lsst/ap/cluster/Utils.h"

#include "lsst/pex/exceptions.h"
#include "lsst/afw/geom/Point.h"


namespace except = lsst::pex::exceptions;
namespace detection = lsst::afw::detection;
namespace geom = lsst::afw::geom;
namespace image = lsst::afw::image;


namespace lsst { namespace ap { namespace cluster {

/** Fills in a histogram of positions for the given sources. The value of each
  * histogram pixel is set to the number of sources falling within that pixel.
  * Sources falling outside the histogram image are either ignored or cause an
  * exception to be raised.
  */
LSST_AP_API void makeSourceHistogram(
    lsst::afw::detection::SourceSet const & sources,
    lsst::afw::image::Image<unsigned short>::Ptr histogram,
    lsst::afw::image::Wcs::Ptr wcs,
    bool ignoreOffImage)
{
    typedef detection::SourceSet::const_iterator SourceIter;
    for (SourceIter i = sources.begin(), e = sources.end(); i != e; ++i) {
       geom::Point2D xy = wcs->skyToPixel((*i)->getRa() * DEGREES_PER_RADIAN,
                                          (*i)->getDec() * DEGREES_PER_RADIAN);
       if (xy[0] < 0.0 || xy[0] + 0.5 >= histogram->getWidth() ||
           xy[1] < 0.0 || xy[1] + 0.5 >= histogram->getHeight()) {
           if (!ignoreOffImage) {
               throw LSST_EXCEPT(except::InvalidParameterException,
                                 "input SourceSet contains sources lying "
                                 "outside the histogram image");
           }
           continue;
       }
       int x = histogram->positionToIndex(xy[0], image::X).first;
       int y = histogram->positionToIndex(xy[1], image::Y).first;
       histogram->operator()(x, y) += 1;
    }
}

}}} // namespace lsst::ap::cluster

