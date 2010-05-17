// -*- lsst-c++ -*-

/**
 * @file
 * @brief   Tests for association pipeline result IO.
 *
 * @ingroup associate
 */

#include <unistd.h>

#include <cmath>
#include <iostream>

#include "boost/bind.hpp"
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE ResultsIoTest
#include "boost/test/unit_test.hpp"

#include "lsst/daf/base/PropertySet.h"
#include "lsst/pex/policy/Policy.h"
#include "lsst/daf/persistence/DbAuth.h"
#include "lsst/daf/persistence/DbStorage.h"
#include "lsst/daf/persistence/Persistence.h"
#include "lsst/daf/persistence/LogicalLocation.h"

#include "lsst/afw/formatters/Utils.h"
#include "lsst/afw/math/Random.h"

#include "lsst/ap/Common.h"
#include "lsst/ap/Results.h"
#include "lsst/ap/ScopeGuard.h"
#include "lsst/ap/Time.h"


using lsst::daf::base::PropertySet;
using lsst::pex::policy::Policy;
using lsst::daf::persistence::DbAuth;
using lsst::daf::persistence::LogicalLocation;
using lsst::daf::persistence::Persistence;
using lsst::daf::base::Persistable;
using lsst::daf::persistence::Storage;
using lsst::daf::persistence::DbStorage;
using lsst::afw::math::Random;

using namespace lsst::ap;
namespace fmt = lsst::afw::formatters;


namespace {

Random & rng() {
    static Random * generator = 0;
    if (generator == 0) {
        TimeSpec ts;
        ts.systemTime();
        generator = new Random(Random::MT19937, static_cast<unsigned long>(ts.tv_sec + ts.tv_nsec));
        std::clog << "\n"
            << "     /\n"
            << "    | Note: Using random number seed " << generator->getSeed() << "\n"
            << "    |       and algorithm " << generator->getAlgorithmName() << "\n"
            << "     \\\n" << std::endl;
    }
    return *generator;
}


std::string const makeTempFile() {
    char name[64];
    std::strncpy(name, "/tmp/BoostArchive.XXXXXX", 63);
    name[63] = 0;
    int const fd = ::mkstemp(name);
    if (fd < 1) {
        BOOST_FAIL("Failed to create temporary file for testing purposes");
    }
    ::close(fd);
    return std::string(name);
}


void initTestData(IdVector & v) {
    int n = 32 + static_cast<int>(rng().flat(0, 32));
    for (int i = 0; i < n; ++i) {
        v.push_back(static_cast<boost::int64_t>(rng().flat(0, 2251799813685249.0)));
    }
}


void initTestData(IdPairVector & v) {
    int n = 32 + static_cast<int>(rng().flat(0, 32));
    IdPair ip;
    for (int i = 0; i < n; ++i) {
        ip.first  = static_cast<boost::int64_t>(rng().flat(0, 2251799813685249.0));
        ip.second = static_cast<boost::int64_t>(rng().flat(0, 2251799813685249.0));
        v.push_back(ip);
    }
}


void initTestData(MatchPairVector & v) {
    int n = 32 + static_cast<int>(rng().flat(0, 32));
    MatchPair mp;
    for (int i = 0; i < n; ++i) {
        mp.setFirst(static_cast<boost::int64_t>(rng().flat(0, 2251799813685249.0)));
        mp.setSecond(static_cast<boost::int64_t>(rng().flat(0, 2251799813685249.0)));
        mp.setDistance(rng().uniform());
        v.push_back(mp);
    }
}


// Make at least a token attempt at generating a unique visit id
// (in-db table name collisions could cause spurious testcase failures)
// This problem will go away once unit tests create their own sandbox databases
int createVisitId() {
    TimeSpec ts;
    ts.systemTime();
    return std::abs(static_cast<int>(ts.tv_sec) +
                    static_cast<int>(ts.tv_nsec) +
                    static_cast<int>(rng().flat(0, 1000)));
}


PropertySet::Ptr createDbTestProps(std::string const & itemName) {
    PropertySet::Ptr props(new PropertySet);
    props->add("visitId", createVisitId());
    props->add("itemName", itemName);
    return props;
}


template <typename VecT>
struct VecTraits {};

template <>
struct VecTraits<IdVector> {
    typedef IdVector Vector;
    typedef PersistableIdVector PersistableVector;
    static Vector & getVector(PersistableVector & pv) { return pv.getIds(); }
};

template <>
struct VecTraits<IdPairVector> {
    typedef IdPairVector Vector;
    typedef PersistableIdPairVector PersistableVector;
    static Vector & getVector(PersistableVector & pv) { return pv.getIdPairs(); }
};

template <>
struct VecTraits<MatchPairVector> {
    typedef MatchPairVector Vector;
    typedef PersistableMatchPairVector PersistableVector;
    static Vector & getVector(PersistableVector & pv) { return pv.getMatchPairs(); }
};


template <typename TraitsT>
void doTestBoost(std::string const & name) {
    typedef typename TraitsT::Vector Vector;
    typedef typename TraitsT::PersistableVector PersistableVector;

    std::string tempFile(makeTempFile());
    ScopeGuard fileGuard(boost::bind(::unlink, tempFile.c_str()));
    Policy::Ptr policy(new Policy);
    PropertySet::Ptr props(new PropertySet);
    LogicalLocation loc(tempFile);
    Persistence::Ptr pers(Persistence::getPersistence(policy));

    BOOST_TEST_MESSAGE("    - BoostStorage I/O test for " << name);

    PersistableVector pvec;
    Vector & vec = TraitsT::getVector(pvec);;
    initTestData(vec);
    // write out data
    {
        Storage::List storageList;
        storageList.push_back(pers->getPersistStorage("BoostStorage", loc));
        pers->persist(pvec, storageList, props);
    }
    // read in data
    {
        Storage::List storageList;
        storageList.push_back(pers->getRetrieveStorage("BoostStorage", loc));
        Persistable::Ptr p = pers->retrieve(name, storageList, props);
        BOOST_REQUIRE_MESSAGE(p.get() != 0, "Failed to retrieve Persistable");
        typename PersistableVector::Ptr pv = boost::dynamic_pointer_cast<PersistableVector, Persistable>(p);
        BOOST_REQUIRE_MESSAGE(pv, "Couldn't cast to " << typeid(PersistableVector).name());
        BOOST_CHECK_MESSAGE(
            TraitsT::getVector(*pv) == vec,
            "persist()/retrieve() resulted in " << typeid(PersistableVector).name() << " corruption"
        );
    }
}


struct LessThanHelper {
    bool operator()(MatchPair const & mp1, MatchPair const & mp2) {
        return mp1.getFirst() < mp2.getFirst() ||
               (mp1.getFirst() == mp2.getFirst() && mp1.getSecond() < mp2.getSecond());
    }

    bool operator()(IdPair const & ip1, IdPair const & ip2) {
        return ip1.first < ip2.first || (ip1.first == ip2.first && ip1.second < ip2.second);
    }

    bool operator()(int64_t const id1, int64_t const id2) { return id1 < id2; }
};


template <typename TraitsT>
void doTestDb(
    std::string const & storageType,
    std::string const & itemName,
    std::string const & name,
    std::string const & templateTable
) {
    typedef typename TraitsT::Vector Vector;
    typedef typename TraitsT::PersistableVector PersistableVector;

    if (!DbAuth::available("lsst10.ncsa.uiuc.edu", "3306")) {
        std::clog << "Skipping database tests - no authorization "
                     "credentials for lsst10.ncsa.uiuc.edu:3306" << std::endl;
        return;
    }

    // Create the required Policy and DataProperty
    Policy::Ptr policy(new Policy);
    std::string policyRoot("Formatter." + name);
    std::string itemRoot(policyRoot + "." + itemName);
    policy->set(itemRoot + ".templateTableName", templateTable);
    policy->set(itemRoot + ".tableNamePattern", "_tmp_test_" + itemName + "_v%(visitId)");

    Persistence::Ptr pers(Persistence::getPersistence(policy));
    LogicalLocation loc("mysql://lsst10.ncsa.uiuc.edu:3306/test");
    PropertySet::Ptr props(createDbTestProps(itemName));

    BOOST_TEST_MESSAGE("    - " << storageType << " I/O test for " << name);

    PersistableVector pvec;
    Vector & vec = TraitsT::getVector(pvec);
    initTestData(vec);
    // write out data
    {
        Storage::List storageList;
        storageList.push_back(pers->getPersistStorage(storageType, loc));
        pers->persist(pvec, storageList, props);
    }
    // and read it back in
    {
        Storage::List storageList;
        Storage::Ptr storage = pers->getRetrieveStorage(storageType, loc);
        storageList.push_back(storage);
        Persistable::Ptr p = pers->retrieve(name, storageList, props);
        BOOST_REQUIRE_MESSAGE(p.get() != 0, "Failed to retrieve Persistable");
        typename PersistableVector::Ptr pv = boost::dynamic_pointer_cast<PersistableVector, Persistable>(p);
        BOOST_REQUIRE_MESSAGE(pv, "Couldn't cast to " << typeid(PersistableVector).name());
        Vector & v = TraitsT::getVector(*pv);
        std::sort(v.begin(),  v.end(),  LessThanHelper());
        std::sort(vec.begin(), vec.end(), LessThanHelper());
        BOOST_CHECK_MESSAGE(
            v == vec,
            "persist()/retrieve() resulted in " << typeid(PersistableVector).name() << " corruption"
        );
        DbStorage * db = dynamic_cast<DbStorage *>(storage.get());
        BOOST_REQUIRE_MESSAGE(db != 0, "Didn't get DbStorage");
        db->startTransaction();
        db->dropTable(fmt::getTableName(policy->getPolicy(policyRoot), props));
        db->endTransaction();
    }
}


template void doTestBoost<VecTraits<IdVector> >(std::string const &);
template void doTestBoost<VecTraits<MatchPairVector> >(std::string const &);
template void doTestDb<VecTraits<IdVector> >(std::string const &, std::string const &,
                                             std::string const &, std::string const &);
template void doTestDb<VecTraits<IdPairVector> >(std::string const &, std::string const &,
                                                 std::string const &, std::string const &);
template void doTestDb<VecTraits<MatchPairVector> >(std::string const &, std::string const &,
                                                    std::string const &, std::string const &);


} // end of anonymous namespace


BOOST_AUTO_TEST_CASE(testBoost) {
    doTestBoost<VecTraits<IdVector> >("PersistableIdVector");
    doTestBoost<VecTraits<IdPairVector> >("PersistableIdPairVector");
    doTestBoost<VecTraits<MatchPairVector> >("PersistableMatchPairVector");
}


BOOST_AUTO_TEST_CASE(testDb) {
    doTestDb<VecTraits<IdVector> >("DbStorage", "Id", "PersistableIdVector", "_tmpl_Id");
    doTestDb<VecTraits<IdPairVector> >("DbStorage", "IdPair", "PersistableIdPairVector", "_tmpl_IdPair");
    doTestDb<VecTraits<MatchPairVector> >("DbStorage", "MatchPair", "PersistableMatchPairVector", "_tmpl_MatchPair");
}


BOOST_AUTO_TEST_CASE(testDbTsv) {
    doTestDb<VecTraits<IdVector> >("DbTsvStorage", "Id", "PersistableIdVector", "_tmpl_Id");
    doTestDb<VecTraits<IdPairVector> >("DbTsvStorage", "IdPair", "PersistableIdPairVector", "_tmpl_IdPair");
    doTestDb<VecTraits<MatchPairVector> >("DbTsvStorage", "MatchPair", "PersistableMatchPairVector", "_tmpl_MatchPair");
}

