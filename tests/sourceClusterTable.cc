/* 
 * LSST Data Management System
 * Copyright 2012 LSST Corporation.
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
#include <iostream>
#include <iterator>
#include <algorithm>
#include <map>

#include "boost/filesystem.hpp"
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE SourceClusterTable
#include "boost/test/unit_test.hpp"

#include "lsst/utils/ieee.h"
#include "lsst/ap/cluster/SourceCluster.h"

using std::vector;
using std::string;
using lsst::afw::geom::radians;
using lsst::afw::coord::IcrsCoord;


struct ExtractSchemaStrings {
    template <typename T>
    void operator()(lsst::afw::table::SchemaItem<T> const & item) const {
        names.push_back(item.field.getName());
        docs.push_back(item.field.getDoc());
        units.push_back(item.field.getUnits());
    }
    
    mutable vector<string> names;
    mutable vector<string> docs;
    mutable vector<string> units;
};


BOOST_AUTO_TEST_CASE(testSourceClusterTable) {
    using lsst::afw::table::Key;
    using lsst::afw::table::Point;
    using lsst::afw::table::Shape;
    using lsst::afw::table::Flux;
    using lsst::afw::table::Covariance;
    using lsst::afw::table::Schema;
    using lsst::afw::table::Flag;
    using lsst::afw::table::Coord;
    using namespace lsst::ap::cluster;


    string const filename = "tests/testTable.fits";

    Schema schema = SourceClusterTable::makeMinimalSchema();

    Key<int> numSources = schema.addField<int>("obs.num", "number of sources");
    Key<Flag> flag = schema.addField<Flag>("flag", "a mysterious boolean value");
    Key<double> timeMin = schema.addField<double>("obs.time.min", "earliest observation time", "mjd");
    Key<double> timeMean = schema.addField<double>("obs.time.mean", "mean observation time", "mjd");
    Key<double> timeMax = schema.addField<double>("obs.time.max", "latest observation time", "mjd");
    Key<Covariance<Point<float> > > coordErr = schema.addField<Covariance<Point<float> > >(
        "coord.err", "covariance matrix for coord field", "rad^2");
    Key<Coord> weightedCoord = schema.addField<Coord>("coord2", "another coordinate", "rad");
    Key<Covariance<Point<float> > > weightedCoordErr = schema.addField<Covariance<Point<float> > >(
        "coord2.err", "covariance matrix for coord2 field", "rad^2");
    Key<int> weightedCoordCount = schema.addField<int>("coord2.count", "sample count for coord2 field");
    Key<int> rNumSources = schema.addField<int>("r.obs.num", "number of sources");
    Key<double> rTimeMin = schema.addField<double>("r.obs.time.min", "earliest observation time", "mjd");
    Key<double> rTimeMax = schema.addField<double>("r.obs.time.max", "latest observation time", "mjd");

    KeyTuple<Flux> uModelFlux = addFluxFields(schema, "u", "flux.model", "model flux", "gremlins");
    KeyTuple<Flux> uApFlux = addFluxFields(schema, "u", "flux.ap", "aperture fluxs", "gremlins");
    KeyTuple<Flux> uPsfFlux = addFluxFields(schema, "u", "flux.psf", "PSF flux", "gremlins");
    KeyTuple<Shape> rShape = addShapeFields(schema, "r", "shape", "bogo-shape");
    KeyTuple<Flux> rFlux = addFluxFields(schema, "r", "flux.inst", "instrumental flux", "gremlins");

    SourceClusterCatalog outcat(SourceClusterTable::make(
        schema, boost::make_shared<SourceClusterIdFactory>(1)));

    outcat.getTable()->defineCoordErr(coordErr);
    outcat.getTable()->defineWeightedMeanCoord(weightedCoord);
    outcat.getTable()->defineWeightedMeanCoordErr(weightedCoordErr);
    outcat.getTable()->defineWeightedMeanCoordCount(weightedCoordCount);
    outcat.getTable()->defineNumSources(numSources);
    outcat.getTable()->defineTimeMin(timeMin);
    outcat.getTable()->defineTimeMean(timeMean);
    outcat.getTable()->defineTimeMax(timeMax);
    outcat.getTable()->defineModelFlux("u", uModelFlux.mean, uModelFlux.err, uModelFlux.count);
    outcat.getTable()->defineApFlux("u", uApFlux.mean, uApFlux.err, uApFlux.count);
    outcat.getTable()->definePsfFlux("u", uPsfFlux.mean, uPsfFlux.err, uPsfFlux.count);
    outcat.getTable()->defineNumSources("r", rNumSources);
    outcat.getTable()->defineTimeMin("r", rTimeMin);
    outcat.getTable()->defineTimeMax("r", rTimeMax);
    outcat.getTable()->defineInstFlux("r", rFlux.mean, rFlux.err, rFlux.count);
    outcat.getTable()->defineShape("r", rShape.mean, rShape.err, rShape.count);

    outcat.getTable()->setMetadata(boost::make_shared<lsst::daf::base::PropertyList>());
    outcat.getTable()->getMetadata()->add("BOGOSITY", 1.0, "bogosity factor");

    vector<string> outfilt = outcat.getTable()->getFilters();
    std::sort(outfilt.begin(), outfilt.end());
    BOOST_CHECK_EQUAL(outfilt.size(), 2u);
    BOOST_CHECK_EQUAL(outfilt[0], "r");
    BOOST_CHECK_EQUAL(outfilt[1], "u");

    {
        PTR(SourceClusterRecord) rec = outcat.getTable()->makeRecord();
        rec->set(flag, true);
        BOOST_CHECK_EQUAL(rec->get(flag), true);
        IcrsCoord c1(0.5*radians, 0.25*radians),
                  c2(1.0*radians, 0.5*radians);
        Eigen::Matrix2d cov1, cov2;
        cov1 << 1, 4,
                4, 2;
        cov2 << 0.5, 0.75,
                0.75, 1.25;
        rec->setCoord(c1);
        rec->setCoordErr(cov1);
        rec->setWeightedMeanCoord(c2);
        rec->setWeightedMeanCoordErr(cov2);
        rec->setWeightedMeanCoordCount(-999);
        rec->setNumSources(15);
        rec->setTimeMin(10.0);
        rec->setTimeMax(20.0);

        BOOST_CHECK_EQUAL(rec->getCoord().getLongitude(), c1.getLongitude());
        BOOST_CHECK_EQUAL(rec->getCoord().getLatitude(), c1.getLatitude());
        IcrsCoord c = rec->get(weightedCoord);
        BOOST_CHECK_EQUAL(c.getLongitude(), c2.getLongitude());
        BOOST_CHECK_EQUAL(c.getLatitude(), c2.getLatitude());
        Eigen::Matrix2d cov = rec->get(coordErr).cast<double>();
        BOOST_CHECK_EQUAL(cov, cov1);
        BOOST_CHECK_EQUAL(rec->getCoordErr(), cov1);
        cov = rec->get(weightedCoordErr).cast<double>();
        BOOST_CHECK_EQUAL(cov, cov2);
        BOOST_CHECK_EQUAL(rec->getWeightedMeanCoordErr().cast<double>(), cov2);
        BOOST_CHECK_EQUAL(rec->getWeightedMeanCoordCount(), -999);
        BOOST_CHECK_EQUAL(rec->getNumSources(), 15);
        BOOST_CHECK_EQUAL(rec->get(numSources), 15);
        BOOST_CHECK_EQUAL(rec->getTimeMin(), 10.0);
        BOOST_CHECK_EQUAL(rec->getTimeMax(), 20.0);
        BOOST_CHECK_EQUAL(rec->getTimeMin(), rec->get(timeMin));
        BOOST_CHECK_EQUAL(rec->getTimeMax(), rec->get(timeMax));

        rec->setNumSources("r", 8);
        BOOST_CHECK_EQUAL(rec->getNumSources("r"), rec->get(rNumSources));
        BOOST_CHECK_EQUAL(rec->getNumSources("r"), 8);
        rec->setTimeMin("r", 13.0);
        BOOST_CHECK_EQUAL(rec->getTimeMin("r"), 13.0);
        BOOST_CHECK_EQUAL(rec->getTimeMin("r"), rec->get(rTimeMin));
        rec->setTimeMax("r", 17.0);
        BOOST_CHECK_EQUAL(rec->getTimeMax("r"), 17.0);
        BOOST_CHECK_EQUAL(rec->getTimeMax("r"), rec->get(rTimeMax));

        rec->setPsfFlux("u", 100.0);
        BOOST_CHECK_EQUAL(rec->getPsfFlux("u"), 100.0);
        BOOST_CHECK_EQUAL(rec->getPsfFlux("u"), rec->get(uPsfFlux.mean));
        rec->setPsfFluxErr("u", 50.0);
        BOOST_CHECK_EQUAL(rec->getPsfFluxErr("u"), 50.0);
        BOOST_CHECK_EQUAL(rec->getPsfFluxErr("u"), rec->get(uPsfFlux.err));
        rec->setPsfFluxCount("u", 5);
        BOOST_CHECK_EQUAL(rec->getPsfFluxCount("u"), 5);
        BOOST_CHECK_EQUAL(rec->getPsfFluxCount("u"), rec->get(uPsfFlux.count));

        rec->setApFlux("u", 200.0);
        BOOST_CHECK_EQUAL(rec->getApFlux("u"), 200.0);
        BOOST_CHECK_EQUAL(rec->getApFlux("u"), rec->get(uApFlux.mean));
        rec->setApFluxErr("u", 60.0);
        BOOST_CHECK_EQUAL(rec->getApFluxErr("u"), 60.0);
        BOOST_CHECK_EQUAL(rec->getApFluxErr("u"), rec->get(uApFlux.err));
        rec->setApFluxCount("u", 6);
        BOOST_CHECK_EQUAL(rec->getApFluxCount("u"), 6);
        BOOST_CHECK_EQUAL(rec->getApFluxCount("u"), rec->get(uApFlux.count));

        rec->setModelFlux("u", 300.0);
        BOOST_CHECK_EQUAL(rec->getModelFlux("u"), 300.0);
        BOOST_CHECK_EQUAL(rec->getModelFlux("u"), rec->get(uModelFlux.mean));
        rec->setModelFluxErr("u", 70.0);
        BOOST_CHECK_EQUAL(rec->getModelFluxErr("u"), 70.0);
        BOOST_CHECK_EQUAL(rec->getModelFluxErr("u"), rec->get(uModelFlux.err));
        rec->setModelFluxCount("u", 7);
        BOOST_CHECK_EQUAL(rec->getModelFluxCount("u"), 7);
        BOOST_CHECK_EQUAL(rec->getModelFluxCount("u"), rec->get(uModelFlux.count));

        rec->setInstFlux("r", 400.0);
        BOOST_CHECK_EQUAL(rec->getInstFlux("r"), 400.0);
        BOOST_CHECK_EQUAL(rec->getInstFlux("r"), rec->get(rFlux.mean));
        rec->setInstFluxErr("r", 80.0);
        BOOST_CHECK_EQUAL(rec->getInstFluxErr("r"), 80.0);
        BOOST_CHECK_EQUAL(rec->getInstFluxErr("r"), rec->get(rFlux.err));
        rec->setInstFluxCount("r", 8);
        BOOST_CHECK_EQUAL(rec->getInstFluxCount("r"), 8);
        BOOST_CHECK_EQUAL(rec->getInstFluxCount("r"), rec->get(rFlux.count));

        outcat.push_back(rec);
    }
    {
        PTR(SourceClusterRecord) rec = outcat.getTable()->makeRecord();

        rec->setCoord(IcrsCoord(2.5*radians, -0.25*radians));
        rec->setNumSources(33);
        rec->setTimeMin(11.0);
        rec->setTimeMean(16.0);
        rec->setTimeMax(21.0);

        rec->setPsfFlux("u", 1000.0);
        rec->set(uApFlux.mean, 2000.0);
        rec->setModelFlux("u", 3000.0);
        rec->set(rFlux.err, 666.0); 

        outcat.push_back(rec);
    }

    outcat.writeFits(filename);

    SourceClusterCatalog incat = SourceClusterCatalog::readFits(filename);

    // check for schema and metadata equality
    {
        BOOST_CHECK_EQUAL(schema, incat.getSchema()); // only checks equality of keys
        BOOST_CHECK_EQUAL(incat.getTable()->getMetadata()->get<double>("BOGOSITY"),
                          outcat.getTable()->getMetadata()->get<double>("BOGOSITY"));
        BOOST_CHECK_EQUAL(incat.getTable()->getMetadata()->nameCount(),
                          outcat.getTable()->getMetadata()->nameCount());

        ExtractSchemaStrings func1;
        schema.forEach(boost::ref(func1));
        ExtractSchemaStrings func2;
        incat.getSchema().forEach(boost::ref(func2));
        BOOST_CHECK(func1.names == func2.names);
        BOOST_CHECK(func1.docs == func2.docs);
        BOOST_CHECK(func1.units == func2.units);
    }

    // check that filters and slots are preserved
    vector<string> infilt = outcat.getTable()->getFilters();
    std::sort(infilt.begin(), infilt.end());
    BOOST_CHECK_EQUAL(outfilt.size(), infilt.size());
    BOOST_CHECK(std::equal(outfilt.begin(), outfilt.end(), infilt.begin()));

    BOOST_CHECK_EQUAL(incat.getTable()->getCoordErrKey(),
                      outcat.getTable()->getCoordErrKey());
    BOOST_CHECK_EQUAL(incat.getTable()->getWeightedMeanCoordKey(),
                      outcat.getTable()->getWeightedMeanCoordKey());
    BOOST_CHECK_EQUAL(incat.getTable()->getWeightedMeanCoordErrKey(),
                      outcat.getTable()->getWeightedMeanCoordErrKey());
    BOOST_CHECK_EQUAL(incat.getTable()->getWeightedMeanCoordCountKey(),
                      outcat.getTable()->getWeightedMeanCoordCountKey());
    BOOST_CHECK_EQUAL(incat.getTable()->getNumSourcesKey(),
                      outcat.getTable()->getNumSourcesKey());
    BOOST_CHECK_EQUAL(incat.getTable()->getTimeMinKey(),
                      outcat.getTable()->getTimeMinKey());
    BOOST_CHECK_EQUAL(incat.getTable()->getTimeMeanKey(),
                      outcat.getTable()->getTimeMeanKey());
    BOOST_CHECK_EQUAL(incat.getTable()->getTimeMaxKey(),
                      outcat.getTable()->getTimeMaxKey());

    for (size_t i = 0; i < infilt.size(); ++i) {
        string f = infilt[i];
        BOOST_CHECK_EQUAL(incat.getTable()->getNumSourcesKey(f),
                          outcat.getTable()->getNumSourcesKey(f));
        BOOST_CHECK_EQUAL(incat.getTable()->getTimeMinKey(f),
                          outcat.getTable()->getTimeMinKey(f));
        BOOST_CHECK_EQUAL(incat.getTable()->getTimeMaxKey(f),
                          outcat.getTable()->getTimeMaxKey(f));

        BOOST_CHECK_EQUAL(incat.getTable()->getPsfFluxKey(f),
                          outcat.getTable()->getPsfFluxKey(f));
        BOOST_CHECK_EQUAL(incat.getTable()->getPsfFluxErrKey(f),
                          outcat.getTable()->getPsfFluxErrKey(f));
        BOOST_CHECK_EQUAL(incat.getTable()->getPsfFluxCountKey(f),
                          outcat.getTable()->getPsfFluxCountKey(f));

        BOOST_CHECK_EQUAL(incat.getTable()->getApFluxKey(f),
                          outcat.getTable()->getApFluxKey(f));
        BOOST_CHECK_EQUAL(incat.getTable()->getApFluxErrKey(f),
                          outcat.getTable()->getApFluxErrKey(f));
        BOOST_CHECK_EQUAL(incat.getTable()->getApFluxCountKey(f),
                          outcat.getTable()->getApFluxCountKey(f));

        BOOST_CHECK_EQUAL(incat.getTable()->getModelFluxKey(f),
                          outcat.getTable()->getModelFluxKey(f));
        BOOST_CHECK_EQUAL(incat.getTable()->getModelFluxErrKey(f),
                          outcat.getTable()->getModelFluxErrKey(f));
        BOOST_CHECK_EQUAL(incat.getTable()->getModelFluxCountKey(f),
                          outcat.getTable()->getModelFluxCountKey(f));

        BOOST_CHECK_EQUAL(incat.getTable()->getInstFluxKey(f),
                          outcat.getTable()->getInstFluxKey(f));
        BOOST_CHECK_EQUAL(incat.getTable()->getInstFluxErrKey(f),
                          outcat.getTable()->getInstFluxErrKey(f));
        BOOST_CHECK_EQUAL(incat.getTable()->getInstFluxCountKey(f),
                          outcat.getTable()->getInstFluxCountKey(f));

        BOOST_CHECK_EQUAL(incat.getTable()->getShapeKey(f),
                          outcat.getTable()->getShapeKey(f));
        BOOST_CHECK_EQUAL(incat.getTable()->getShapeErrKey(f),
                          outcat.getTable()->getShapeErrKey(f));
        BOOST_CHECK_EQUAL(incat.getTable()->getShapeCountKey(f),
                          outcat.getTable()->getShapeCountKey(f));
    }

    // compare data in records
    {
        SourceClusterRecord const & a = outcat[0];
        SourceClusterRecord const & b = incat[0];
        BOOST_CHECK_EQUAL(a.getCoord(), b.getCoord());
        BOOST_CHECK_EQUAL(a.getCoordErr(), b.getCoordErr());
        BOOST_CHECK_EQUAL(a.getWeightedMeanCoord(), b.getWeightedMeanCoord());
        BOOST_CHECK_EQUAL(a.getWeightedMeanCoordErr(), b.getWeightedMeanCoordErr());
        BOOST_CHECK_EQUAL(a.getWeightedMeanCoordCount(), b.getWeightedMeanCoordCount());
        BOOST_CHECK_EQUAL(a.getNumSources(), b.getNumSources());
        BOOST_CHECK_EQUAL(a.getTimeMin(), b.getTimeMin());
        BOOST_CHECK(lsst::utils::isnan(a.getTimeMean()) &&
                    lsst::utils::isnan(b.getTimeMean()));
        BOOST_CHECK_EQUAL(a.getTimeMax(), b.getTimeMax());
        BOOST_CHECK_EQUAL(a.getNumSources("r"), b.getNumSources("r"));
        BOOST_CHECK_EQUAL(a.getTimeMin("r"), b.getTimeMin("r"));
        BOOST_CHECK_EQUAL(a.getTimeMax("r"), b.getTimeMax("r"));
        BOOST_CHECK_EQUAL(a.getInstFlux("r"), b.getInstFlux("r"));
        BOOST_CHECK_EQUAL(a.getInstFluxErr("r"), b.getInstFluxErr("r"));
        BOOST_CHECK_EQUAL(a.getInstFluxCount("r"), b.getInstFluxCount("r"));
        BOOST_CHECK_EQUAL(a.getPsfFlux("u"), b.getPsfFlux("u"));
        BOOST_CHECK_EQUAL(a.getPsfFluxErr("u"), b.getPsfFluxErr("u"));
        BOOST_CHECK_EQUAL(a.getPsfFluxCount("u"), b.getPsfFluxCount("u"));
        BOOST_CHECK_EQUAL(a.getApFlux("u"), b.getApFlux("u"));
        BOOST_CHECK_EQUAL(a.getApFluxErr("u"), b.getApFluxErr("u"));
        BOOST_CHECK_EQUAL(a.getApFluxCount("u"), b.getApFluxCount("u"));
        BOOST_CHECK_EQUAL(a.getModelFlux("u"), b.getModelFlux("u"));
        BOOST_CHECK_EQUAL(a.getModelFluxErr("u"), b.getModelFluxErr("u"));
        BOOST_CHECK_EQUAL(a.getModelFluxCount("u"), b.getModelFluxCount("u"));
        BOOST_CHECK_EQUAL(a.get(flag), b.get(flag));
    }
    {
        SourceClusterRecord const & a = outcat[1];
        SourceClusterRecord const & b = incat[1];
        BOOST_CHECK_EQUAL(a.get(flag), b.get(flag));
        BOOST_CHECK_EQUAL(a.getCoord(), b.getCoord());
        BOOST_CHECK_EQUAL(a.getNumSources(), b.getNumSources());
        BOOST_CHECK_EQUAL(a.getTimeMin(), b.getTimeMin());
        BOOST_CHECK_EQUAL(a.getTimeMean(), b.getTimeMean());
        BOOST_CHECK_EQUAL(a.getTimeMax(), b.getTimeMax());
        BOOST_CHECK_EQUAL(a.getInstFluxErr("r"), b.getInstFluxErr("r"));
        BOOST_CHECK_EQUAL(a.getPsfFlux("u"), b.getPsfFlux("u"));
        BOOST_CHECK_EQUAL(a.getApFlux("u"), b.getApFlux("u"));
        BOOST_CHECK_EQUAL(a.getModelFlux("u"), b.getModelFlux("u"));
        Eigen::Matrix2d cova = a.getCoordErr(), covb = b.getCoordErr(),
                        cov2a = a.getWeightedMeanCoordErr(), cov2b = b.getWeightedMeanCoordErr();
        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                BOOST_CHECK(lsst::utils::isnan(cova(i,j)) &&
                            lsst::utils::isnan(covb(i,j)) &&
                            lsst::utils::isnan(cov2a(i,j)) &&
                            lsst::utils::isnan(cov2b(i,j)));
            }
        }
        IcrsCoord ca = a.getWeightedMeanCoord(), cb = b.getWeightedMeanCoord();
        BOOST_CHECK(lsst::utils::isnan(ca.getLongitude().asRadians()) &&
                    lsst::utils::isnan(ca.getLatitude().asRadians()) &&
                    lsst::utils::isnan(cb.getLongitude().asRadians()) &&
                    lsst::utils::isnan(cb.getLatitude().asRadians()));
        BOOST_CHECK(a.getWeightedMeanCoordCount() == 0 && b.getWeightedMeanCoordCount() == 0);
        BOOST_CHECK(a.getNumSources("r") == 0 && b.getNumSources("r") == 0);
        BOOST_CHECK(lsst::utils::isnan(a.getTimeMin("r")) && lsst::utils::isnan(b.getTimeMin("r")));
        BOOST_CHECK(lsst::utils::isnan(a.getTimeMax("r")) && lsst::utils::isnan(b.getTimeMax("r")));
        BOOST_CHECK(lsst::utils::isnan(a.getInstFlux("r")) && lsst::utils::isnan(b.getInstFlux("r")));
        BOOST_CHECK_EQUAL(a.getInstFluxErr("r"), b.getInstFluxErr("r"));
        BOOST_CHECK(a.getInstFluxCount("r") == 0 && b.getInstFluxCount("r") == 0);
        BOOST_CHECK(lsst::utils::isnan(a.getPsfFluxErr("u")) && lsst::utils::isnan(b.getPsfFluxErr("u")));
        BOOST_CHECK(a.getPsfFluxCount("u") == 0 && b.getPsfFluxCount("u") == 0);
        BOOST_CHECK(lsst::utils::isnan(a.getApFluxErr("u")) && lsst::utils::isnan(b.getApFluxErr("u")));
        BOOST_CHECK(a.getApFluxCount("u") == 0 && b.getApFluxCount("u") == 0);
        BOOST_CHECK(lsst::utils::isnan(a.getModelFluxErr("u")) && lsst::utils::isnan(b.getModelFluxErr("u")));
        BOOST_CHECK(a.getModelFluxCount("u") == 0 && b.getModelFluxCount("u") == 0);
    }

    boost::filesystem::remove(filename);
}

