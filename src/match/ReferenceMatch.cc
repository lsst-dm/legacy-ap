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
  * @brief  Reference matching implementation.
  * @author Serge Monkewitz
  */

#include "lsst/ap/match/ReferenceMatch.h"

#include <math.h>
#include <algorithm>
#include <limits>
#include <string>
#include <vector>

#include "boost/scoped_ptr.hpp"

#include "lsst/utils/ieee.h"

#include "lsst/pex/logging/Log.h"
#include "lsst/pex/policy/DefaultPolicyFile.h"
#include "lsst/pex/policy/Policy.h"
#include "lsst/afw/coord/Coord.h"
#include "lsst/afw/image/ImageUtils.h"

#include "lsst/ap/Common.h"
#include "lsst/ap/utils/Arena.h"
#include "lsst/ap/utils/Csv.h"
#include "lsst/ap/utils/SmallPtrVector.h"
#include "lsst/ap/utils/SpatialUtils.h"
#include "lsst/ap/match/BBox.h"
#include "lsst/ap/match/ExposureInfo.h"
#include "lsst/ap/match/ReferencePosition.h"
#include "lsst/ap/match/detail/SweepStructure.h"


using std::fabs;
using std::min;
using std::max;
using std::numeric_limits;
using std::vector;
using std::string;

using lsst::pex::logging::Log;
using lsst::pex::policy::Policy;
using lsst::pex::policy::DefaultPolicyFile;
using lsst::afw::image::PixelZeroPos;

using lsst::ap::match::detail::SphericalSweep;
using lsst::ap::utils::angularSeparation;
using lsst::ap::utils::Arena;
using lsst::ap::utils::cartesianToSpherical;
using lsst::ap::utils::clampPhi;
using lsst::ap::utils::CsvDialect;
using lsst::ap::utils::CsvReader;
using lsst::ap::utils::CsvWriter;
using lsst::ap::utils::degrees;
using lsst::ap::utils::maxAlpha;
using lsst::ap::utils::sphericalToCartesian;


namespace lsst { namespace ap { namespace match {

namespace {

// converts from milliarcsec/yr to rad/day
double const RADY_PER_MASD = 1.32734751843815467101961424328e-11;
double const RAD_PER_MAS = 4.84813681109535993589914102357e-9;

// -- Classes for match participants ----

class MatchablePos;
class MatchableRef;


/** @internal  A match between a reference position and a position.
  */
class RefPosMatch {
public:
     RefPosMatch(MatchableRef *ref,
                 MatchablePos *pos,
                 Eigen::Vector3d const &v,
                 double angularSeparation
                ) :
         _ref(ref),
         _pos(pos),
         _sc(cartesianToSpherical(v)),
         _angSep(angularSeparation),
         _refCount(2)
    { }

    ~RefPosMatch() {
#ifndef NDEBUG
        _ref = 0;
        _pos = 0;
#endif
    }

    MatchableRef const *getRef() const {
        return _ref;
    }
    MatchableRef *getRef() {
        return _ref;
    }
    MatchablePos const *getPos() const {
        return _pos;
    }
    MatchablePos *getPos() {
        return _pos;
    }
    double getRa() const {
        return _sc.x();
    }
    double getDecl() const {
        return _sc.y();
    }
    double getAngularSeparation() const {
        return _angSep;
    }

    int getRefCount() const {
        return _refCount;
    }
    int decrementRefCount() {
        return --_refCount;
    }

private:
    MatchableRef *_ref;
    MatchablePos *_pos;
    Eigen::Vector2d _sc;
    double _angSep;
    int _refCount;
};


/** @internal Base class for matchable entities.
  */
class Matchable : public BBox {
public:
    typedef lsst::ap::utils::SmallPtrVector<RefPosMatch, 3u> MatchVector;

    Matchable(std::string const &record) :
        _closestId(0),
        _closestSep(numeric_limits<double>::infinity()),
        _nMatches(0),
        _refCount(0),
        _matches(),
        _record(record)
    { }

    virtual ~Matchable() { }

    int64_t getClosestId() const {
        return _closestId;
    }
    int getNumMatches() const {
        return _nMatches;
    }
    std::string const & getRecord() const {
        return _record;
    }

    MatchVector & getMatches() {
        return _matches;
    }
    MatchVector const & getMatches() const {
        return _matches;
    }

    void appendMatch(int64_t id, RefPosMatch *m) {
        _matches.push_back(m);
        ++_refCount;
        ++_nMatches;
        if (m->getAngularSeparation() < _closestSep) {
            _closestId = id;
            _closestSep = m->getAngularSeparation();
        }
    }
    int decrementRefCount() {
        return --_refCount;
    }
    int getRefCount() const {
        return _refCount;
    }

private:
    int64_t _closestId;
    double _closestSep;
    int _nMatches;
    int _refCount;
    MatchVector _matches;
    std::string _record;
};


/** @internal  A matchable position and associated epoch.
  *
  * Coordinates are ICRS, except that the observer need not necessarily
  * be located at the SSB (i.e. may be topocentric). For galaxies and
  * distant stars, the distinction is moot. However, parallax will cause
  * barycentric and topocentric coordinates to disagree for near-by
  * stars. Since the positions being matched to reference objects have
  * no associated parallax/motion parameters, one cannot correct for this
  * effect here.
  */
class MatchablePos : public Matchable {
public:
    MatchablePos(int64_t id,
                 double ra,
                 double decl,
                 double epoch,
                 double radius,
                 std::string const &record);

    virtual ~MatchablePos() { }

    int64_t getId() const {
        return _id;
    }
    double getEpoch() const {
        return _epoch;
    }
    Eigen::Vector2d const & getSphericalCoordinates() const {
        return _sc;
    }
    Eigen::Vector3d const & getPosition() const {
        return _p;
    }
    double getRadius() const {
        return _radius;
    }

    // BBox API
    virtual double getMinCoord0() const {
        return _sc.x() - _alpha;
    }
    virtual double getMaxCoord0() const {
        return _sc.x() + _alpha;
    }
    virtual double getMinCoord1() const {
        return clampPhi(_sc.y() - _radius);
    }
    virtual double getMaxCoord1() const {
        return clampPhi(_sc.y() + _radius);
    }

private:
    int64_t _id;
    double _epoch;
    Eigen::Vector2d _sc;
    Eigen::Vector3d _p;
    double _radius;
    double _alpha;
};

MatchablePos::MatchablePos(
    int64_t id,    ///< Unique id.
    double ra,     ///< ICRS right ascension, rad.
    double decl,   ///< ICRS declination, rad
    double epoch,  ///< Epoch of position, MJD
    double radius, ///< Match radius, rad.
    std::string const &record ///< Ancillary match output data
) :
    Matchable(record),
    _id(id),
    _epoch(epoch),
    _sc(ra, decl),
    _p(sphericalToCartesian(ra, decl)),
    _radius(radius),
    _alpha(maxAlpha(radius, max(fabs(getMinCoord1()), fabs(getMaxCoord1()))))
{ }


/** @internal  A matchable reference position
  */
class MatchableRef : public Matchable {
public:
    MatchableRef(ReferencePosition const &rp,
                 std::string const &record);

    virtual ~MatchableRef() { }

    ReferencePosition & getReferencePosition() {
        return _rp;
    }
    ReferencePosition const & getReferencePosition() const {
        return _rp;
    }

    // BBox API
    virtual double getMinCoord0() const {
        return _rp.getMinCoord0();
    }
    virtual double getMaxCoord0() const {
        return _rp.getMaxCoord0();
    }
    virtual double getMinCoord1() const {
        return _rp.getMinCoord1();
    }
    virtual double getMaxCoord1() const {
        return _rp.getMaxCoord1();
    }

private:
    ReferencePosition _rp;
};

MatchableRef::MatchableRef(
    ReferencePosition const &rp,
    std::string const &record
) :
    Matchable(record),
    _rp(rp)
{ }

inline bool operator<(std::pair<double, MatchableRef *> const &r1,
                      std::pair<double, MatchableRef *> const &r2) {
    return r1.first > r2.first;
}


/** @internal  A reference position with coverage information.
  */
class RefWithCov : public BBox {
public:
    RefWithCov(ReferencePosition const &rp,
               std::string const &record);

    virtual ~RefWithCov() { }

    ReferencePosition & getReferencePosition() {
        return _rp;
    }
    ReferencePosition const & getReferencePosition() const {
        return _rp;
    }

    std::string const & getRecord() const {
        return _record;
    }

    void writeCoverage(lsst::ap::utils::CsvWriter &writer) const {
        for (size_t i = 0; i < sizeof(_filterCov)/sizeof(int); ++i) {
            writer.appendField(_filterCov[i]);
        }
    }

    int isCovered() const {
        return _covered;
    }

    void appendMatch(ExposureInfo *e) {
        _covered = true;
        ++_filterCov[e->getFilterId()];
    }

    // BBox API
    virtual double getMinCoord0() const {
        return _rp.getMinCoord0();
    }
    virtual double getMaxCoord0() const {
        return _rp.getMaxCoord0();
    }
    virtual double getMinCoord1() const {
        return _rp.getMinCoord1();
    }
    virtual double getMaxCoord1() const {
        return _rp.getMaxCoord1();
    }

private:
    ReferencePosition _rp;
    std::string _record;
    bool _covered;
    int _filterCov[6];
};

RefWithCov::RefWithCov(
    ReferencePosition const &rp,
    std::string const &record
) :
    _rp(rp),
    _record(record),
    _covered(false)
{
    for (size_t i = 0; i < sizeof(_filterCov)/sizeof(int); ++i) {
        _filterCov[i] = 0;
    }
}

inline bool operator<(std::pair<double, RefWithCov *> const &r1,
                      std::pair<double, RefWithCov *> const &r2) {
    return r1.first > r2.first;
}


// --  Classes for reading in reference objects and positions ----

/** Reads positions from a declination sorted CSV file.
  *
  * For now, matching works with a fixed radius, so this is equivalent
  * to producing positions sorted by minimum bounding box declination.
  */
class MatchablePosReader {
public:
    MatchablePosReader(std::string const &path,
                       lsst::pex::policy::Policy::Ptr policy,
                       CsvDialect const &outDialect,
                       double radius);

    ~MatchablePosReader();

    double getMinEpoch() const { return _minEpoch; }
    double getMaxEpoch() const { return _maxEpoch; }
    double getDefaultEpoch() const { return _defaultEpoch; }

    /** Returns a string consisting of pre-formatted all NULL
      * output for requested position table columns.
      */
    std::string const & getNullRecord() const {
        return _nullRecord;
    }

    /** Returns true if the positions read in have an associated output record.
      */
    bool haveOutputRecord() const {
        return !_columns.empty();
    }

    /** Returns @c true if all positions have been read and returned.
      */
    bool isDone() const {
        return _next == 0;
    }

    /** Returns the minimum declination of the next positions bounding box.
      * Assumes that <tt>isDone() == false</tt>.
      */
    double peek() const {
        return _next->getMinCoord1();
    }

    /** Returns the next MatchablePos. Assumes that <tt>isDone() == false</tt>.
      */
    MatchablePos *next() {
        MatchablePos *p = _next;
        _read();
        return p;
    }

    /** Frees a MatchablePos object created by this PosReader.
      */
    void destroy(MatchablePos *p) {
        _arena.destroy(p);
    }

private:
    void _read();
    void _scan(lsst::ap::utils::CsvReader &reader);

    Arena<MatchablePos> _arena;
    double _decl;
    double _defaultEpoch;
    double _minEpoch;
    double _maxEpoch;
    double _radius;
    MatchablePos *_next;
    boost::scoped_ptr<CsvReader> _reader;
    int _idCol;
    int _epochCol;
    int _raCol;
    int _declCol;
    double _raScale;
    double _declScale;
    std::vector<int> _columns;
    std::string _nullRecord;
    std::string _record;
    std::stringstream _buf;
    CsvWriter _writer;
};

MatchablePosReader::MatchablePosReader(
    std::string const &path,
    lsst::pex::policy::Policy::Ptr policy,
    CsvDialect const &outDialect,
    double radius
) :
    _arena(),
    _decl(-M_PI*0.5),
    _defaultEpoch(policy->getDouble("epoch")),
    _minEpoch(_defaultEpoch),
    _maxEpoch(_defaultEpoch),
    _radius(radius),
    _next(0),
    _reader(),
    _columns(),
    _nullRecord(),
    _record(),
    _buf(),
    _writer(_buf, outDialect)
{
    typedef std::vector<std::string>::const_iterator Iter;

    Log log(Log::getDefaultLog(), "lsst.ap.match");
    log.log(Log::INFO, "Opening position table " + path);

    // create CSV reader
    CsvDialect inDialect(policy->getPolicy("csvDialect"));
    _reader.reset(new CsvReader(path, inDialect, !policy->exists("fieldNames")));
    if (policy->exists("fieldNames")) {
        _reader->setFieldNames(policy->getStringArray("fieldNames"));
    }
    // get standard column info
    _idCol = _reader->getIndexOf(policy->getString("idColumn"));
    if (policy->exists("epochColumn")) {
        _epochCol = _reader->getIndexOf(policy->getString("epochColumn"));
    } else {
        _epochCol = -1;
    }
    _raCol = _reader->getIndexOf(policy->getString("raColumn"));
    _declCol = _reader->getIndexOf(policy->getString("declColumn"));
    _raScale = policy->getDouble("raScale");
    _declScale = policy->getDouble("declScale");
    if (_idCol < 0 || _raCol < 0 || _declCol < 0) {
        throw LSST_EXCEPT(pexExcept::RuntimeErrorException,
                          "Position table does not contain unique id, "
                          "right ascension, or declination column(s)");
    }
    // compute vector of output field indexes and NULL record
    if (policy->exists("outputFields")) { 
        vector<string> ofs = policy->getStringArray("outputFields");
        _writer.endRecord();
        _buf.str("");
        if (ofs.size() == 1 && ofs[0] == "*") {
            int nFields = static_cast<int>(_reader->getFieldNames().size());
            for (int i = 0; i < nFields; ++i) {
                _columns.push_back(i);
                _writer.appendNull();
            }
        } else {
            for (Iter i = ofs.begin(), e = ofs.end(); i != e; ++i) {
                int idx = _reader->getIndexOf(*i);
                if (idx < 0) {
                    throw LSST_EXCEPT(pexExcept::InvalidParameterException,
                                      "Position table has no column named " +
                                      *i);
                }
                _columns.push_back(idx);
                _writer.appendNull();
            }
        }
        _nullRecord = _buf.str();
    }
    if (_epochCol >= 0) {
        if (policy->exists("minEpoch") && policy->exists("maxEpoch")) {
            _minEpoch = policy->getDouble("minEpoch");
            _maxEpoch = policy->getDouble("maxEpoch");
        } else  {
            log.log(Log::INFO, "Scanning position table to determine min/max "
                    "epoch of input positions");
            CsvReader reader(path, _reader->getDialect(),
                             !policy->exists("fieldNames"));
            _scan(reader);
            log.format(Log::INFO, "Scanned %llu records",
                       static_cast<unsigned long long>(reader.getNumRecords()));
        }
    }
    log.format(Log::INFO, "  - time range is [%.3f, %.3f] MJD",
               _minEpoch, _maxEpoch);
    // read first record
    _read();
}

MatchablePosReader::~MatchablePosReader() { }

void MatchablePosReader::_read() {
    typedef std::vector<int>::const_iterator Iter;

    if (_reader->isDone()) {
        _next = 0;
        return;
    }
    // retrieve column values
    int64_t id = _reader->get<int64_t>(_idCol);
    double epoch = _defaultEpoch;
    if (_epochCol >= 0) {
        epoch = _reader->get<double>(_epochCol);
    }
    double ra = _reader->get<double>(_raCol)*_raScale;
    double decl = _reader->get<double>(_declCol)*_declScale;
    // check for NULLs and illegal values
    if (_reader->isNull(_idCol)) {
        throw LSST_EXCEPT(pexExcept::RuntimeErrorException,
                          "NULL unique id found in position table");
    }
    if (lsst::utils::isnan(ra) ||
        lsst::utils::isnan(decl) ||
        lsst::utils::isnan(epoch)) {
        throw LSST_EXCEPT(pexExcept::RuntimeErrorException,
                          "Position table contains NULL or NaN right "
                          "ascension, declination, or epoch");
    }
    if (decl < -M_PI*0.5 || decl > M_PI*0.5) {
        throw LSST_EXCEPT(pexExcept::RuntimeErrorException,
                          "Invalid declination found in position table");
    }
    if (epoch < J2000_MJD - 200.0*DAYS_PER_JY ||
        epoch > J2000_MJD + 200.0*DAYS_PER_JY) {
        throw LSST_EXCEPT(pexExcept::RuntimeErrorException,
                          "Position table epoch is not within 200 years of "
                          "J2000. Check your units - MJD required.");
    }
    // check that input file is declination sorted
    if (decl < _decl) {
        throw LSST_EXCEPT(pexExcept::RuntimeErrorException,
                          "Position table is not sorted by declination");
    }
    // Construct ancillary column output string
    if (!_columns.empty()) {
        _writer.endRecord();
        _buf.str("");
        _record.clear();
        for (Iter i = _columns.begin(), e = _columns.end(); i != e; ++i) {
            _writer.appendField(_reader->get<char const *>(*i));
        }
        _record = _buf.str();
    }
    _reader->nextRecord();
    // Create Pos object.
    _next = new (_arena) MatchablePos(id, ra, decl, epoch, _radius, _record);
    _decl = decl;
}

void MatchablePosReader::_scan(lsst::ap::utils::CsvReader &reader) {
    if (reader.isDone()) {
        return;
    }
    _minEpoch = reader.get<double>(_epochCol);
    _maxEpoch = _minEpoch;
    if (lsst::utils::isnan(_minEpoch)) {
        throw LSST_EXCEPT(pexExcept::RuntimeErrorException,
                          "Position table contains NULL or NaN epoch");
    }
    reader.nextRecord();
    while (!reader.isDone()) {
        double epoch = reader.get<double>(_epochCol);
        if (lsst::utils::isnan(epoch)) {
            throw LSST_EXCEPT(pexExcept::RuntimeErrorException,
                              "Position table contains NULL or NaN epoch");
        }
        if (epoch < _minEpoch) {
            _minEpoch = epoch;
        } else if (epoch > _maxEpoch) {
            _maxEpoch = epoch;
        }
        reader.nextRecord();
    }
}


/** @internal  Base class for declination sorted reference catalog readers.
  *
  * Since bounding boxes for reference positions have varying sizes, one
  * must read ahead in declination and reorder to produce reference positions
  * sorted by minimum bounding box declination. This class computes the 
  * read-ahead amount as part of its duties.
  */
class RefReaderBase {
public:
    RefReaderBase(std::string const &path,
                  lsst::pex::policy::Policy::Ptr inPolicy,
                  CsvDialect const &outDialect,
                  double minMatchEpoch,
                  double maxMatchEpoch,
                  double parallaxThresh);

    virtual ~RefReaderBase();

    double getMinEpoch() const { return _minEpoch; }
    double getMaxEpoch() const { return _maxEpoch; }
    double getMinMatchEpoch() const { return _minMatchEpoch; }
    double getMaxMatchEpoch() const { return _maxMatchEpoch; }
    double getMaxParallax() const { return _maxParallax; }
    double getMaxAngularVelocity() const { return _maxAngularVelocity; }
    double getParallaxThresh() const { return _parallaxThresh; }
    double getReadAhead() const { return _readAhead; }

    /** Returns a string consisting of pre-formatted all NULL
      * output for requested reference catalog columns.
      */
    std::string const & getNullRecord() const {
        return _nullRecord;
    }

    /** Returns true if the reference catalog entries read in have an
      * associated output record.
      */
    bool haveOutputRecord() const {
        return !_columns.empty();
    }

protected:
    ReferencePosition const * _readReferencePosition();

    std::string const & getRecord() const {
        return _record;
    }

    double _decl;
    double _readAhead;
    boost::scoped_ptr<CsvReader> _reader;

private:
    void _scan(lsst::ap::utils::CsvReader &reader,
               bool needEpochStats,
               bool needMotionStats);

    ReferencePosition _pos;
    double _defaultEpoch;
    double _minEpoch;
    double _maxEpoch;
    double _minMatchEpoch;
    double _maxMatchEpoch;
    double _maxParallax;
    double _maxAngularVelocity;
    double _parallaxThresh;
    int _idCol;
    int _epochCol;
    int _raCol;
    int _declCol;
    int _muRaCol;
    int _muDeclCol;
    int _vRadialCol;
    int _parallaxCol;
    double _raScale;
    double _declScale;
    double _muRaScale;
    double _muDeclScale;
    double _vRadialScale;
    double _parallaxScale;
    bool _muRaTrueAngle;
    std::vector<int> _columns;
    std::string _nullRecord;
    std::string _record;
    std::stringstream _buf;
    CsvWriter _writer;
};

RefReaderBase::RefReaderBase(
    std::string const &path,
    lsst::pex::policy::Policy::Ptr policy,
    lsst::ap::utils::CsvDialect const &outDialect,
    double minMatchEpoch,
    double maxMatchEpoch,
    double parallaxThresh
) :
    _decl(-M_PI*0.5),
    _readAhead(0.0),
    _reader(),
    _pos(0, 0.0, 0.0),
    _defaultEpoch(policy->getDouble("epoch")),
    _minEpoch(_defaultEpoch),
    _maxEpoch(_defaultEpoch),
    _minMatchEpoch(min(minMatchEpoch, maxMatchEpoch)),
    _maxMatchEpoch(max(minMatchEpoch, maxMatchEpoch)),
    _maxParallax(0.0),
    _maxAngularVelocity(0.0),
    _parallaxThresh(parallaxThresh),
    _columns(),
    _nullRecord(),
    _record(),
    _buf(),
    _writer(_buf, outDialect)
{
    typedef std::vector<std::string>::const_iterator Iter;

    Log log(Log::getDefaultLog(), "lsst.ap.match");
    log.log(Log::INFO, "Opening reference catalog " + path);

    // create CSV reader
    CsvDialect inDialect(policy->getPolicy("csvDialect"));
    _reader.reset(new CsvReader(path, inDialect, !policy->exists("fieldNames")));
    if (policy->exists("fieldNames")) {
        _reader->setFieldNames(policy->getStringArray("fieldNames"));
    }
    // get standard column info
    _idCol = _reader->getIndexOf(policy->getString("idColumn"));
    if (policy->exists("epochColumn")) {
        _epochCol = _reader->getIndexOf(policy->getString("epochColumn"));
    } else {
        _epochCol = -1;
    }
    _raCol = _reader->getIndexOf(policy->getString("raColumn"));
    _declCol = _reader->getIndexOf(policy->getString("declColumn"));
    _muRaCol = _reader->getIndexOf(policy->getString("muRaColumn"));
    _muDeclCol = _reader->getIndexOf(policy->getString("muDeclColumn"));
    _vRadialCol = _reader->getIndexOf(policy->getString("vRadialColumn"));
    _parallaxCol = _reader->getIndexOf(policy->getString("parallaxColumn"));
    _raScale = policy->getDouble("raScale");
    _declScale = policy->getDouble("declScale");
    _muRaScale = policy->getDouble("muRaScale");
    _muDeclScale = policy->getDouble("muDeclScale");
    _vRadialScale = policy->getDouble("vRadialScale");
    _parallaxScale = policy->getDouble("parallaxScale");
    _muRaTrueAngle = policy->getBool("muRaTrueAngle");
    if (_idCol < 0 || _raCol < 0 || _declCol < 0) {
        throw LSST_EXCEPT(pexExcept::RuntimeErrorException,
                          "Reference catalog doesn't contain unique id, "
                          "right ascension, or declination column(s)");
    }
    // compute vector of output field indexes and NULL record
    if (policy->exists("outputFields")) {
        vector<string> ofs = policy->getStringArray("outputFields");
        _writer.endRecord();
        _buf.str("");
        if (ofs.size() == 1 && ofs[0] == "*") {
            int nFields = static_cast<int>(_reader->getFieldNames().size());
            for (int i = 0; i < nFields; ++i) {
                _columns.push_back(i);
                _writer.appendNull();
            }
        } else {
            for (Iter i = ofs.begin(), e = ofs.end(); i != e; ++i) {
                int idx = _reader->getIndexOf(*i);
                if (idx < 0) {
                    throw LSST_EXCEPT(pexExcept::InvalidParameterException,
                                      "Reference catalog has no column "
                                      "named " + *i);
                }
                _columns.push_back(idx);
                _writer.appendNull();
            }
        }
        _nullRecord = _buf.str();
    }
    if (_epochCol >= 0 ||
        (_muRaCol >= 0 && _muDeclCol >= 0 && _parallaxCol >= 0)) {

        bool needEpochStats = false;
        bool needMotionStats = false;
        if (_epochCol >= 0) {
            if (policy->exists("minEpoch") && policy->exists("maxEpoch")) {
                _minEpoch = policy->getDouble("minEpoch");
                _maxEpoch = policy->getDouble("maxEpoch");
            } else {
                needEpochStats = true;
            }
        }
        if (_muRaCol >= 0 && _muDeclCol >= 0 && _parallaxCol >= 0) {
            if (policy->exists("maxParallax") &&
                policy->exists("maxAngularVelocity")) {
                _maxParallax = policy->getDouble("maxParallax")*RAD_PER_MAS;
                _maxAngularVelocity =
                    policy->getDouble("maxAngularVelocity")*RADY_PER_MASD;
            } else {
                needMotionStats = true;
            } 
        }

        if (needEpochStats || needMotionStats) {
            log.log(Log::INFO, "Scanning reference catalog to determine "
                    "time-range and/or maximum velocity/parallax");
            CsvReader reader(path, _reader->getDialect(),
                             !policy->exists("fieldNames"));
            _scan(reader, needEpochStats, needMotionStats);
            log.format(Log::INFO, "Scanned %llu records",
                       static_cast<unsigned long long>(reader.getNumRecords()));
 
        }
    }
    // determine read-ahead amount
    double dtMax = max(fabs(_maxEpoch - _minMatchEpoch),
                       fabs(_maxMatchEpoch - _minEpoch));
    _readAhead = 2.0*_maxParallax + _maxAngularVelocity*dtMax;
    log.format(Log::INFO, "  - time range is [%.3f, %.3f] MJD",
               _minEpoch, _maxEpoch);
    log.format(Log::INFO, "  - max parallax is %.3f milliarcsec",
               _maxParallax/RAD_PER_MAS);
    log.format(Log::INFO, "  - max angular velocity is %.3f milliarcsec/yr",
               _maxAngularVelocity/RADY_PER_MASD);
    log.format(Log::INFO, "  - read-ahead is %.3f milliarcsec",
               _readAhead/RAD_PER_MAS);
}

RefReaderBase::~RefReaderBase() { }

ReferencePosition const * RefReaderBase::_readReferencePosition() {
    typedef std::vector<int>::const_iterator Iter;

    if (_reader->isDone()) {
        return 0;
    }
    // retrieve column values
    int64_t id = _reader->get<int64_t>(_idCol);
    double epoch = _defaultEpoch;
    if (_epochCol >= 0) {
        epoch = _reader->get<double>(_epochCol);
    }
    double ra = _reader->get<double>(_raCol)*_raScale;
    double decl = _reader->get<double>(_declCol)*_declScale;
    // check for NULLs and illegal values
    if (_reader->isNull(_idCol)) {
        throw LSST_EXCEPT(pexExcept::RuntimeErrorException,
                          "NULL unique id found in reference catalog");
    }
    if (lsst::utils::isnan(ra) ||
        lsst::utils::isnan(decl) ||
        lsst::utils::isnan(epoch)) {
        throw LSST_EXCEPT(pexExcept::RuntimeErrorException,
                          "Reference catalog record contains NULL/NaN right "
                          "ascension, declination, or epoch");
    }
    if (decl < -M_PI*0.5 || decl > M_PI*0.5) {
        throw LSST_EXCEPT(pexExcept::RuntimeErrorException,
                          "Invalid declination found in reference catalog");
    }
    if (epoch < J2000_MJD - 200.0*DAYS_PER_JY ||
        epoch > J2000_MJD + 200.0*DAYS_PER_JY) {
        throw LSST_EXCEPT(pexExcept::RuntimeErrorException,
                          "Reference catalog epoch is not within 200 years "
                          "of J2000. Check your units - MJD required.");
    }
    // check that input file actually is declination sorted
    if (decl < _decl) {
        throw LSST_EXCEPT(pexExcept::RuntimeErrorException,
                          "Reference catalog is not sorted by declination");
    }
    // Construct ancillary column output string
    if (!_columns.empty()) {
        _writer.endRecord();
        _buf.str("");
        for (Iter i = _columns.begin(), e = _columns.end(); i != e; ++i) {
            _writer.appendField(_reader->get<char const *>(*i));
        }
        _record = _buf.str();
    }
    // create reference position
    _pos = ReferencePosition(id, ra, decl, epoch);
    if (_muRaCol >= 0 && _muDeclCol >= 0 && _parallaxCol >= 0) {
        // have motion parameters
        double muRa = _reader->get<double>(_muRaCol)*_muRaScale;
        double muDecl = _reader->get<double>(_muDeclCol)*_muDeclScale;
        double parallax = _reader->get<double>(_parallaxCol)*_parallaxScale;
        double vRadial = (_vRadialCol < 0) ? 0.0 :
                          _reader->get<double>(_vRadialCol)*_vRadialScale;
        if (parallax < 0.0) {
            throw LSST_EXCEPT(pexExcept::RuntimeErrorException,
                              "Reference catalog contains negative parallax");
        }
        if (!lsst::utils::isnan(muRa) &&
            !lsst::utils::isnan(muDecl) &&
            !lsst::utils::isnan(vRadial) &&
            !lsst::utils::isnan(parallax)) {
            // motion parameters are non-null
            _pos.setMotion(muRa, muDecl, parallax, vRadial,
                           _muRaTrueAngle, parallax > _parallaxThresh);
        }
    }
    _pos.setTimeRange(_minMatchEpoch, _maxMatchEpoch);
    _decl = decl;
    _reader->nextRecord();
    return &_pos;
}

void RefReaderBase::_scan(lsst::ap::utils::CsvReader &reader,
                          bool needEpochStats,
                          bool needMotionStats)
{
    if (reader.isDone()) {
        return;
    }
    if (needEpochStats) {
        _minEpoch = reader.get<double>(_epochCol);
        _maxEpoch = _minEpoch;
    }
    while (!reader.isDone()) {
        if (needEpochStats) {
            double epoch = reader.get<double>(_epochCol);
            if (lsst::utils::isnan(epoch)) {
                throw LSST_EXCEPT(pexExcept::RuntimeErrorException,
                                  "Record contains NULL or NaN epoch");
            }
            if (epoch < _minEpoch) {
                _minEpoch = epoch;
            } else if (epoch > _maxEpoch) {
                _maxEpoch = epoch;
            }
        }
        if (needMotionStats) {
            double parallax = reader.get<double>(_parallaxCol)*_parallaxScale;
            double muRa = reader.get<double>(_muRaCol)*_muRaScale;
            double muDecl = reader.get<double>(_muDeclCol)*_muDeclScale;
            if (!lsst::utils::isnan(parallax) &&
                !lsst::utils::isnan(muRa) &&
                !lsst::utils::isnan(muDecl)) {
                if (!_muRaTrueAngle) {
                    double decl = reader.get<double>(_declCol)*_declScale;
                    if (decl < -M_PI*0.5 || decl > M_PI*0.5 ||
                        lsst::utils::isnan(decl)) {
                        throw LSST_EXCEPT(pexExcept::RuntimeErrorException,
                                          "Invalid declination found in "
                                          "reference catalog");
                    }
                    muRa *= std::cos(decl);
                }
                if (parallax > _maxParallax) {
                    _maxParallax = parallax;
                }
                double v = sqrt(muRa*muRa + muDecl*muDecl);
                if (v > _maxAngularVelocity) {
                    _maxAngularVelocity = v;
                }
            }
        }
        reader.nextRecord();
    }
}


/** @internal  Reads matchable reference positions from a declination sorted
  *            CSV file.
  */
template <typename Ref>
class RefReader : public RefReaderBase {
public:
    RefReader(std::string const &path,
              lsst::pex::policy::Policy::Ptr inPolicy,
              CsvDialect const &outDialect,
              double minMatchEpoch,
              double maxMatchEpoch,
              double parallaxThresh);

    ~RefReader();

    /** Returns @c true if all positions have been read and returned.
      */
    bool isDone() const {
        return _heap.empty();
    }

    /** Returns the minimum declination of the next reference positions
      * bounding box.  Assumes that <tt>isDone() == false</tt>.
      */
    double peek() const {
        return _heap.front().first;
    }

    /** Returns the next MatchablePos. Assumes that <tt>isDone() == false</tt>.
      */
    Ref *next() {
        while (_decl - peek() <= _readAhead && !_reader->isDone()) {
            _read();
        }
        std::pop_heap(_heap.begin(), _heap.end());
        Ref *r = _heap.back().second;
        _heap.pop_back();
        return r;
    }

    /** Frees a MatchableRef object created by this RefReader.
      */
    void destroy(Ref *r) {
        _arena.destroy(r);
    }

private:
    void _read() {
        ReferencePosition const *p = _readReferencePosition();
        if (p != 0) {
            Ref *r = new (_arena) Ref(*p, getRecord());
            _heap.push_back(std::pair<double, Ref *>(r->getMinCoord1(), r));
            std::push_heap(_heap.begin(), _heap.end());
        }
    }

    Arena<Ref> _arena;
    std::vector<std::pair<double, Ref *> > _heap;
};

template <typename Ref>
RefReader<Ref>::RefReader(
    std::string const &path,
    lsst::pex::policy::Policy::Ptr policy,
    lsst::ap::utils::CsvDialect const &outDialect,
    double minMatchEpoch,
    double maxMatchEpoch,
    double parallaxThresh
) :
    RefReaderBase(path,
                  policy,
                  outDialect,
                  minMatchEpoch,
                  maxMatchEpoch,
                  parallaxThresh)
{
    _read();
}

template <typename Ref>
RefReader<Ref>::~RefReader() { }


// -- Matchers ----

/** @internal Reference position to position matcher.
  */
class RefPosMatcher {
public:
    RefPosMatcher();
    ~RefPosMatcher();

    // Sweep event call-backs
    void operator()(MatchablePos *p) {
        _finish(p);
    }
    void operator()(MatchableRef *r) {
        _finish(r);
    }
    void operator()(MatchableRef *r, MatchablePos *p) {
         _candidateMatch(r, p);
    }
    void operator()(MatchablePos *p, MatchableRef *r) {
         _candidateMatch(r, p);
    }

    void match(CsvWriter &writer,
               RefReader<MatchableRef> &refReader,
               MatchablePosReader &posReader);

private:
    void _writeMatch(RefPosMatch const *m);
    void _finish(MatchablePos *p);
    void _finish(MatchableRef *r);
    void _candidateMatch(MatchableRef *r, MatchablePos *p);

    // memory allocator for matches
    Arena<RefPosMatch> _arena;

    // data structures for matching
    SphericalSweep<MatchablePos> _posSweep;
    SphericalSweep<MatchableRef> _refSweep;

    CsvWriter *_writer;
    RefReader<MatchableRef> *_refReader;
    MatchablePosReader *_posReader;
};

RefPosMatcher::RefPosMatcher() :
    _arena(),
    _posSweep(),
    _refSweep(),
    _writer(0),
    _refReader(0),
    _posReader(0)
{ }

RefPosMatcher::~RefPosMatcher() {
    _writer = 0;
    _refReader = 0;
    _posReader = 0;
}

/** Matches a reference catalog to a position table.
  */
void RefPosMatcher::match(CsvWriter &writer,
                          RefReader<MatchableRef> &refReader,
                          MatchablePosReader &posReader)
{
    _writer = &writer;
    _refReader = &refReader;
    _posReader = &posReader;
    while (true) {
        if (_posReader->isDone()) {
            while (!_refReader->isDone()) {
                double refDecl = _refReader->peek();
                _posSweep.advance(refDecl, *this);
                MatchableRef *r = _refReader->next();
                _posSweep.search(r, *this);
                _finish(r);
            }
            break;
        } else if (_refReader->isDone()) {
            while (!_posReader->isDone()) {
                double posDecl = _posReader->peek();
                _refSweep.advance(posDecl, *this);
                MatchablePos *p = _posReader->next();
                _refSweep.search(p, *this);
                _finish(p);
            }
            break;
        }
        double posDecl = _posReader->peek();
        double refDecl = _refReader->peek();
        _refSweep.advance(posDecl, *this);
        _posSweep.advance(refDecl, *this);
        if (refDecl < posDecl) {
            MatchableRef *r = _refReader->next();
            _posSweep.search(r, *this);
            _refSweep.insert(r);
        } else {
            MatchablePos *p = _posReader->next();
            _refSweep.search(p, *this);
            _posSweep.insert(p);
        }
    }
    _refSweep.clear(*this);
    _posSweep.clear(*this); 
    _writer = 0;
    _refReader = 0;
    _posReader = 0;
}

/** Writes out a match pair.
  */
void RefPosMatcher::_writeMatch(RefPosMatch const *m) {
    MatchableRef const *r = m->getRef();
    MatchablePos const *p = m->getPos();
    _writer->appendField(r->getReferencePosition().getId());
    _writer->appendField(p->getId());
    _writer->appendField(degrees(m->getRa()));
    _writer->appendField(degrees(m->getDecl()));
    _writer->appendField(degrees(m->getAngularSeparation())*3600.0);
    _writer->appendField(r->getNumMatches());
    _writer->appendField(p->getNumMatches());
    _writer->appendField(r->getClosestId() == p->getId());
    _writer->appendField(p->getClosestId() == r->getReferencePosition().getId());
    _writer->appendField(r->getReferencePosition().getFlags());
    if (_refReader->haveOutputRecord()) {
        _writer->write(_writer->getDialect().getDelimiter());
        _writer->write(r->getRecord());
    }
    if (_posReader->haveOutputRecord()) {
        _writer->write(_writer->getDialect().getDelimiter());
        _writer->write(p->getRecord());
    }
    _writer->endRecord();
}

/** Called when @a p is removed from it's sweep structure, i.e.
  * when all matches (if any) for @a p have been found.
  */
void RefPosMatcher::_finish(MatchablePos *p) {
    typedef MatchablePos::MatchVector::iterator Iter;

    if (p->getNumMatches() == 0) {
        // p has no matches - write out p and delete it immediately.
        _writer->appendNull();
        _writer->appendField(p->getId());
        _writer->appendNull();
        _writer->appendNull();
        _writer->appendNull();
        _writer->appendNull();
        _writer->appendField(0);
        _writer->appendNull();
        _writer->appendNull();
        _writer->appendNull();
        if (_refReader->haveOutputRecord()) {
            _writer->write(_writer->getDialect().getDelimiter());
            _writer->write(_refReader->getNullRecord());
        }
        if (_posReader->haveOutputRecord()) {
            _writer->write(_writer->getDialect().getDelimiter());
            _writer->write(p->getRecord());
        }
        _writer->endRecord();
        _posReader->destroy(p);
        return;
    }
    // All matches for p have been found - write out all match
    // pairs (r, p) where all matches for r have also been found
    //
    // Note that writing (r, p) is delayed until this late stage
    // because the match table output includes flags that indicate
    // whether r is the closest match to p (and vice versa), as well
    // as match counts for r and p.
    for (Iter i = p->getMatches().begin(), e = p->getMatches().end();
         i != e; ++i) {
        RefPosMatch *m = *i;
        MatchableRef *r = m->getRef();
        if (m->decrementRefCount() == 0) {
            // p had the only reference to m
            _writeMatch(m);
            // delete m
            _arena.destroy(m);
            if (r->decrementRefCount() == 0) {
                // all matches involving r have been written
                _refReader->destroy(r);
            }
            p->decrementRefCount();
        }
    }
    if (p->getRefCount() == 0) {
        _posReader->destroy(p);
    }
}

/** Called when @a r is removed from it's sweep structure, i.e.
  * when all matches (if any) for @a r have been found.
  */
void RefPosMatcher::_finish(MatchableRef *r) {
    typedef MatchablePos::MatchVector::iterator Iter;

    if (r->getNumMatches() == 0) {
        // r has no matches - write out r and delete it immediately.
        _writer->appendField(r->getReferencePosition().getId());
        _writer->appendNull();
        _writer->appendNull();
        _writer->appendNull();
        _writer->appendNull();
        _writer->appendField(0);
        _writer->appendNull();
        _writer->appendNull();
        _writer->appendNull();
        _writer->appendNull();
        if (_refReader->haveOutputRecord()) {
            _writer->write(_writer->getDialect().getDelimiter());
            _writer->write(r->getRecord());
        }
        if (_posReader->haveOutputRecord()) {
            _writer->write(_writer->getDialect().getDelimiter());
            _writer->write(_posReader->getNullRecord());
        }
        _writer->endRecord();
        _refReader->destroy(r);
        return;
    }
    // All matches for r have been found - write out all match
    // pairs (r, p) where all matches for r have also been found
    //
    // Note that writing (r, p) is delayed until this late stage
    // because the match table output includes flags that indicate
    // whether r is the closest match to p (and vice versa), as well
    // as match counts for r and p.
    for (Iter i = r->getMatches().begin(), e = r->getMatches().end();
         i != e; ++i) {
        RefPosMatch *m = *i;
        MatchablePos *p = m->getPos();
        if (m->decrementRefCount() == 0) {
            // r had the only reference to m
            _writeMatch(m);
            // delete m
            _arena.destroy(m);
            if (p->decrementRefCount() == 0) {
                // all matches involving p have been written
                _posReader->destroy(p);
            }
            r->decrementRefCount();
        }
    }
    if (r->getRefCount() == 0) {
        _refReader->destroy(r);
    }
}

/** Called on candidate matches (r, p).
  */
void RefPosMatcher::_candidateMatch(MatchableRef *r, MatchablePos *p) {
    Eigen::Vector3d v = r->getReferencePosition().getPosition(p->getEpoch());
    double sep = angularSeparation(v, p->getPosition());
    if (sep <= p->getRadius()) {
        // got a match
        RefPosMatch *m = new (_arena) RefPosMatch(r, p, v, sep);
        r->appendMatch(p->getId(), m);
        p->appendMatch(r->getReferencePosition().getId(), m);
    }
}


/** @internal Matches reference positions to a list of exposures.
  */
class RefExpMatcher {
public:
    RefExpMatcher();
    ~RefExpMatcher();

    // Sweep event call-backs
    void operator()(ExposureInfo *e) { }

    void operator()(RefWithCov *r) {
        _finish(r);
    }
    void operator()(RefWithCov *r, ExposureInfo *e) {
         _candidateMatch(r, e);
    }
    void operator()(ExposureInfo *e, RefWithCov *r) {
         _candidateMatch(r, e);
    }

    void match(CsvWriter &writer,
               RefReader<RefWithCov> &refReader,
               std::vector<ExposureInfo::Ptr> &exposures);

private:
    struct ExposureInfoComparator {
        bool operator()(ExposureInfo::Ptr const &e1,
                        ExposureInfo::Ptr const &e2) const
        {
            return e1->getMinCoord1() < e2->getMinCoord1();
        }
    };

    void _finish(RefWithCov *r);
    void _candidateMatch(RefWithCov *r, ExposureInfo *e);

    // heap that reorders to produce declination sorted output
    std::vector<std::pair<double, RefWithCov *> > _heap;

    // data structures for matching
    SphericalSweep<ExposureInfo> _expSweep;
    SphericalSweep<RefWithCov> _refSweep;

    CsvWriter *_writer;
    RefReader<RefWithCov> *_refReader;
    double _maxDecl;
};

RefExpMatcher::RefExpMatcher() :
    _expSweep(),
    _refSweep(),
    _writer(0),
    _refReader(0),
    _maxDecl(-M_PI*0.5)
{ }

RefExpMatcher::~RefExpMatcher() {
    _writer = 0;
    _refReader = 0;
}

/** Matches a reference catalog to a list of exposures.
  */
void RefExpMatcher::match(CsvWriter &writer,
                          RefReader<RefWithCov> &refReader,
                          std::vector<ExposureInfo::Ptr> &exposures)
{
    typedef std::vector<ExposureInfo::Ptr>::const_iterator Iter;

    _writer = &writer;
    _refReader = &refReader;

    // sort exposures by minimum bounding box declination
    std::sort(exposures.begin(), exposures.end(), ExposureInfoComparator());
    Iter i = exposures.begin();
    Iter const end = exposures.end();
    // run the standard sweep line algorithm
    while (true) {
        if (i == end) {
            while (!_refReader->isDone()) {
                double refDecl = _refReader->peek();
                _expSweep.advance(refDecl, *this);
                RefWithCov *r = _refReader->next();
                _expSweep.search(r, *this);
                _finish(r);
            }
            break;
        } else if (_refReader->isDone()) {
            for (; i != end; ++i) {
                double expDecl = (*i)->getMinCoord1();
                _refSweep.advance(expDecl, *this);
                ExposureInfo *info = i->get();
                _refSweep.search(info, *this);
            }
            break;
        }
        double expDecl = (*i)->getMinCoord1();
        double refDecl = _refReader->peek();
        _refSweep.advance(expDecl, *this);
        _expSweep.advance(refDecl, *this);
        if (refDecl < expDecl) {
            RefWithCov *r = _refReader->next();
            _expSweep.search(r, *this);
            _refSweep.insert(r);
        } else {
            ExposureInfo *info = i->get();
            ++i;
            _refSweep.search(info, *this);
            _expSweep.insert(info);
        }
    }
    _refSweep.clear(*this);
    _expSweep.clear(*this); 
    while (!_heap.empty()) {
        // write out any reference catalog entries remaining in the
        // declination min-heap.
        std::pop_heap(_heap.begin(), _heap.end());
        RefWithCov *r = _heap.back().second;
        _heap.pop_back();
        if (_refReader->haveOutputRecord()) {
            _writer->write(r->getRecord());
            _writer->write(_writer->getDialect().getDelimiter());
        }
        r->writeCoverage(*_writer);
        _writer->endRecord();
        _refReader->destroy(r);
    }
    _writer = 0;
    _refReader = 0;
}

/** Called when r is removed from its sweep structure, i.e. when all
  * matches for r have been found.
  */
void RefExpMatcher::_finish(RefWithCov *r) {
    if (!r->isCovered()) {
        _refReader->destroy(r);
    } else {
        // insert r into declination heap, and update max declination
        // seen so far.
        double decl = r->getReferencePosition().getSphericalCoords().y();
        if (decl >= _maxDecl) {
            _maxDecl = decl;
        }
        _heap.push_back(std::pair<double, RefWithCov *>(decl, r));
        std::push_heap(_heap.begin(), _heap.end());
        while (!_heap.empty() &&
               _maxDecl - _heap.front().first > _refReader->getReadAhead()) {
            std::pop_heap(_heap.begin(), _heap.end());
            r = _heap.back().second;
            _heap.pop_back();
            if (_refReader->haveOutputRecord()) {
                _writer->write(r->getRecord());
                _writer->write(_writer->getDialect().getDelimiter());
            }
            r->writeCoverage(*_writer);
            _writer->endRecord();
            _refReader->destroy(r);
        }
    }
}

namespace {
    void _nullDeleter(void *) { }
}

/** Called on candidate matches (r, i).
  */
void RefExpMatcher::_candidateMatch(RefWithCov *r, ExposureInfo *e) {
    double epoch = e->getEpoch();
    Eigen::Vector3d v = r->getReferencePosition().getPosition(
        epoch, e->getEarthPosition());
    Eigen::Vector2d sc = cartesianToSpherical(v);
    lsst::afw::coord::IcrsCoord coord(degrees(sc.x()), degrees(sc.y()));
    // avoid at least one memory allocation
    boost::shared_ptr<lsst::afw::coord::Coord> sky(&coord, _nullDeleter);
    lsst::afw::geom::PointD p = e->getWcs()->skyToPixel(sky);
    if (p.getX() >= PixelZeroPos - 0.5 &&
        p.getX() <= e->getWidth() + PixelZeroPos - 0.5 &&
        p.getY() >= PixelZeroPos - 0.5 &&
        p.getY() <= e->getHeight() + PixelZeroPos - 0.5) {
        r->appendMatch(e);
    }
}

} // namespace lsst::ap::match::<anonymous>


/** Matches a declination sorted reference catalog (stored as a CSV file)
  * to a table of positions.
  */
LSST_AP_API void referenceMatch(
    std::string const &refPath,    ///< Path to declination sorted reference catalog CSV file.
    std::string const &posPath,    ///< Path to declination sorted position CSV file.
    std::string const &matchPath,  ///< Match output file path.
    lsst::pex::policy::Policy::Ptr refInPolicy, ///< Policy describing input reference catalog.
                                                ///  See policy/ReferenceCatalogDictionary 
                                                ///  for parameters and default values.
    lsst::pex::policy::Policy::Ptr posInPolicy, ///< Policy describing input position table.
                                                ///  See policy/PositionTableDictionary.paf
                                                ///  for parameters and default values.
    lsst::pex::policy::Policy::Ptr matchPolicy  ///< Policy describing match parameters.
                                                ///  See policy/ReferenceMatchDictionary.paf
                                                ///  for parameters and default values.
) {
    Log log(Log::getDefaultLog(), "lsst.ap.match");
    log.log(Log::INFO, "Matching reference catalog to position table...");

    // Merge input policies with defaults
    Policy::Ptr refPolicy;
    if (refInPolicy) {
        refPolicy.reset(new Policy(*refInPolicy, true));
    } else {
        refPolicy.reset(new Policy());
    }
    DefaultPolicyFile refDef("ap", "ReferenceCatalogDictionary.paf", "policy");
    refPolicy->mergeDefaults(Policy(refDef));
    Policy::Ptr posPolicy;
    if (posInPolicy) {
        posPolicy.reset(new Policy(*posInPolicy, true));
    } else {
        posPolicy.reset(new Policy());
    }
    DefaultPolicyFile posDef("ap", "PositionTableDictionary.paf", "policy");
    posPolicy->mergeDefaults(Policy(posDef));
    Policy::Ptr xPolicy;
    if (matchPolicy) {
        xPolicy.reset(new Policy(*matchPolicy, true));
    } else {
        xPolicy.reset(new Policy());
    }
    DefaultPolicyFile xDef("ap", "ReferenceMatchDictionary.paf", "policy");
    xPolicy->mergeDefaults(Policy(xDef));

    // Create readers, writer and matcher
    CsvDialect outDialect(xPolicy->getPolicy("csvDialect"));
    MatchablePosReader posReader(
        posPath, posPolicy, outDialect,
        xPolicy->getDouble("radius")*RADIANS_PER_ARCSEC);
    RefReader<MatchableRef> refReader(
        refPath, refPolicy, outDialect, posReader.getMinEpoch(),
        posReader.getMaxEpoch(), xPolicy->getDouble("parallaxThresh")*RAD_PER_MAS);
    CsvWriter writer(matchPath, outDialect);
    RefPosMatcher matcher;

    log.log(Log::INFO, "Starting reference catalog to position table match");
    matcher.match(writer, refReader, posReader);
    log.format(Log::INFO, "Wrote %llu records to output match table %s",
               static_cast<unsigned long long>(writer.getNumRecords()),
               matchPath.c_str());
}


/** Computes the number of times a reference catalog should have been observed in
  * each filter with an ideal observatory, given a set of exposures. The per-filter
  * observation counts are appended as columns [ugrizy]Cov. Reference catalog
  * entries not falling on any of the given exposures are dropped from the output.
  */
LSST_AP_API void referenceFilter(
    std::string const &refPath,  ///< Reference catalog path
    std::string const &filtPath, ///< Filtered output catalog path
    std::vector<ExposureInfo::Ptr> &exposures,  ///< Exposures to filter against - 
                                                ///  reordered by the call.
    lsst::pex::policy::Policy::Ptr refInPolicy, ///< Policy describing input reference catalog.
                                                ///  See policy/ReferenceCatalogDictionary 
                                                ///  for parameters and default values.
    lsst::pex::policy::Policy::Ptr matchPolicy  ///< Policy describing match parameters.
                                                ///  See policy/ReferenceMatchDictionary.paf
                                                ///  for parameters and default values.
) {
    typedef std::vector<ExposureInfo::Ptr>::const_iterator Iter;
    if (exposures.empty()) {
        throw LSST_EXCEPT(pexExcept::InvalidParameterException,
                          "no input exposure information");
    }

    Log log(Log::getDefaultLog(), "lsst.ap.match");
    log.log(Log::INFO, "Filtering out reference catalog entries not "
            "observable in any exposure...");

    // Merge input policies with defaults
    Policy::Ptr refPolicy;
    if (refInPolicy) {
        refPolicy.reset(new Policy(*refInPolicy, true));
    } else {
        refPolicy.reset(new Policy());
    }
    DefaultPolicyFile refDef("ap", "ReferenceCatalogDictionary.paf", "policy");
    refPolicy->mergeDefaults(Policy(refDef));
    Policy::Ptr xPolicy;
    if (matchPolicy) {
        xPolicy.reset(new Policy(*matchPolicy, true));
    } else {
        xPolicy.reset(new Policy());
    }
    DefaultPolicyFile xDef("ap", "ReferenceMatchDictionary.paf", "policy");
    xPolicy->mergeDefaults(Policy(xDef));

    // determine min/max epoch of exposures
    Iter i = exposures.begin();
    double minEpoch = (*i)->getEpoch();
    double maxEpoch = minEpoch;
    ++i;
    for (Iter e = exposures.end(); i != e; ++i) {
        double epoch = (*i)->getEpoch();
        if (epoch < minEpoch) {
            minEpoch = epoch;
        } else if (epoch > maxEpoch) {
            maxEpoch = epoch;
        }
    }
    log.format(Log::INFO, "Time range of exposures is [%.3f, %.3f] MJD",
               minEpoch, maxEpoch);

    // Create reader, writer and matcher
    CsvDialect outDialect(xPolicy->getPolicy("csvDialect"));
    RefReader<RefWithCov> refReader(
        refPath, refPolicy, outDialect, minEpoch, maxEpoch,
        xPolicy->getDouble("parallaxThresh")*RAD_PER_MAS);
    CsvWriter writer(filtPath, outDialect);
    RefExpMatcher matcher;

    log.log(Log::INFO, "Starting reference catalog to exposure match");
    matcher.match(writer, refReader, exposures);
    log.format(Log::INFO, "Wrote %llu records to output match table %s",
               static_cast<unsigned long long>(writer.getNumRecords()),
               filtPath.c_str());
}

}}} // namespace lsst::ap::match

