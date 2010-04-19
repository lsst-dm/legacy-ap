// -*- lsst-c++ -*-
/** @file
  * @brief Class representing a PT1 sky tile.
  *
  * @ingroup ap
  * @author Serge Monkewitz
  */
#ifndef LSST_AP_CLUSTER_PT1SKYTILE_H
#define LSST_AP_CLUSTER_PT1SKYTILE_H

#include <vector>
#include <utility>

#include "lsst/afw/detection/Source.h"

#include "../Common.h"


namespace lsst { namespace ap { namespace cluster {

/** A PT1 sky-tile is a single pixel obtained from the quadrilateralized
  * spherical cube. Post-PT1, this code will be generalized and moved
  * to lsst::skypix. For now, this small class provides the minimal
  * functionality required by the PT1 SourceAssoc pipeline.
  */
class LSST_AP_API PT1SkyTile {
public:
    PT1SkyTile(int resolution, int root, int x, int y, int id);
    ~PT1SkyTile();

    int getId() const { return _id; }

    bool contains(double theta, double phi) const;

    std::pair<size_t, size_t> prune(
        lsst::afw::detection::SourceSet & sources) const;

private:
    int _resolution;
    int _root;
    int _x;
    int _y;
    int _id;
};

}}} // namespace lsst:ap::cluster

#endif // LSST_AP_CLUSTER_PT1SKYTILE_H
