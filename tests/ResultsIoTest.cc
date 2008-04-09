// -*- lsst-c++ -*-

/**
 * @file
 * @brief   Tests for association pipeline result IO.
 *
 * @ingroup associate
 */

#include <unistd.h>

#include <cmath>

#include <boost/version.hpp>
#include <boost/bind.hpp>
#if BOOST_VERSION < 103400
#   include <boost/test/auto_unit_test.hpp>
#   define BOOST_TEST_MESSAGE BOOST_MESSAGE
#else
#   include <boost/test/unit_test.hpp>
#endif

#include <lsst/daf/base/DataProperty.h>
#include <lsst/pex/policy/Policy.h>
#include <lsst/daf/persistence/Persistence.h>
#include <lsst/daf/persistence/LogicalLocation.h>

#include <lsst/ap/Common.h>
#include <lsst/ap/Exceptions.h>
#include <lsst/ap/Random.h>
#include <lsst/ap/Results.h>
#include <lsst/ap/ScopeGuard.h>
#include <lsst/ap/Time.h>
#include <lsst/ap/Utils.h>
#include <lsst/ap/io/ResultFormatters.h>


using lsst::daf::base::DataProperty;
using lsst::pex::policy::Policy;
using lsst::daf::persistence::LogicalLocation;
using lsst::daf::persistence::Persistence;
using lsst::daf::persistence::Persistable;
using lsst::daf::persistence::Storage;
using lsst::daf::persistence::DbStorage;

using namespace lsst::ap;


namespace {

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
    initRandom();
    int n = 32 + static_cast<int>(std::floor(uniformRandom()*32.0));
    for (int i = 0; i < n; ++i) {
        v.push_back(static_cast<int64_t>(std::floor(uniformRandom()*2251799813685249.0)));
    }
}


void initTestData(IdPairVector & v) {
    initRandom();
    int n = 32 + static_cast<int>(std::floor(uniformRandom()*32.0));
    IdPair ip;
    for (int i = 0; i < n; ++i) {
        ip.first  = static_cast<int64_t>(std::floor(uniformRandom()*2251799813685249.0));
        ip.second = static_cast<int64_t>(std::floor(uniformRandom()*2251799813685249.0));
        v.push_back(ip);
    }
}


void initTestData(MatchPairVector & v) {
    initRandom();
    int n = 32 + static_cast<int>(std::floor(uniformRandom()*32.0));
    MatchPair mp;
    for (int i = 0; i < n; ++i) {
        mp.setFirst(static_cast<int64_t>(std::floor(uniformRandom()*2251799813685249.0)));
        mp.setSecond(static_cast<int64_t>(std::floor(uniformRandom()*2251799813685249.0)));
        mp.setDistance(uniformRandom());
        v.push_back(mp);
    }
}


// Make at least a token attempt at generating a unique visit id
// (in-db table name collisions could cause spurious testcase failures)
int64_t createVisitId() {
    TimeSpec ts;
    ts.systemTime();
    return static_cast<int64_t>(ts.tv_sec )*INT64_C(1000000000) +
           static_cast<int64_t>(ts.tv_nsec/1000)*1000 +
           static_cast<int64_t>(std::floor(uniformRandom()*1000.0));
}


DataProperty::PtrType createDbTestProps(std::string const & itemName) {
    DataProperty::PtrType props = DataProperty::createPropertyNode("root");
    props->addProperty(DataProperty("visitId", createVisitId()));
    props->addProperty(DataProperty("itemName", itemName));
    return props;
}


template <typename VecT>
void doTestBoost(std::string const & name) {
    std::string           tempFile(makeTempFile());
    ScopeGuard            fileGuard(boost::bind(::unlink, tempFile.c_str()));
    Policy::Ptr           policy(new Policy);
    DataProperty::PtrType props(DataProperty::createPropertyNode("root"));
    LogicalLocation       loc(tempFile);
    Persistence::Ptr      pers(Persistence::getPersistence(policy));

    BOOST_TEST_MESSAGE("    - BoostStorage I/O test for " << name);

    VecT vec;
    initTestData(vec);
    // write out data
    {
        Storage::List storageList;
        storageList.push_back(pers->getPersistStorage("BoostStorage", loc));
        pers->persist(vec, storageList, props);
    }
    // read in data
    {
        Storage::List storageList;
        storageList.push_back(pers->getRetrieveStorage("BoostStorage", loc));
        Persistable::Ptr p = pers->retrieve(name, storageList, props);
        BOOST_REQUIRE_MESSAGE(p.get() != 0, "Failed to retrieve Persistable");
        typename VecT::Ptr v = boost::dynamic_pointer_cast<VecT, Persistable>(p);
        BOOST_REQUIRE_MESSAGE(v, "Couldn't cast to " << typeid(VecT).name());
        BOOST_CHECK_MESSAGE(
            *v == vec,
            "persist()/retrieve() resulted in " << typeid(VecT).name() << " corruption"
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


template <typename VecT>
void doTestDb(
    std::string const & storageType,
    std::string const & itemName,
    std::string const & name
) {
    // Create the required Policy and DataProperty
    Policy::Ptr           policy(new Policy);
    Persistence::Ptr      pers(Persistence::getPersistence(policy));
    LogicalLocation       loc("mysql://lsst10.ncsa.uiuc.edu:3306/test");
    DataProperty::PtrType props(createDbTestProps(itemName));

    BOOST_TEST_MESSAGE("    - " << storageType << " I/O test for " << name);

    VecT vec;
    initTestData(vec);
    // write out data
    {
        Storage::List storageList;
        storageList.push_back(pers->getPersistStorage(storageType, loc));
        pers->persist(vec, storageList, props);
    }
    // and read it back in
    {
        Storage::List storageList;
        Storage::Ptr storage = pers->getRetrieveStorage(storageType, loc);
        storageList.push_back(storage);
        Persistable::Ptr p = pers->retrieve(name, storageList, props);
        BOOST_REQUIRE_MESSAGE(p.get() != 0, "Failed to retrieve Persistable");
        typename VecT::Ptr v = boost::dynamic_pointer_cast<VecT, Persistable>(p);
        std::sort(v->begin(),  v->end(),  LessThanHelper());
        std::sort(vec.begin(), vec.end(), LessThanHelper());
        BOOST_REQUIRE_MESSAGE(v, "Couldn't cast to " << typeid(VecT).name());
        BOOST_CHECK_MESSAGE(
            *v == vec,
            "persist()/retrieve() resulted in " << typeid(VecT).name() << " corruption"
        );
        DbStorage * db = dynamic_cast<DbStorage *>(storage.get());
        BOOST_REQUIRE_MESSAGE(db != 0, "Didn't get DbStorage");
        std::string fmtPolicyName = "Formatter." + name;
        Policy::Ptr fmtPolicy;
        if (policy->exists(fmtPolicyName)) {
            fmtPolicy = policy->getPolicy(fmtPolicyName);
        }
        db->startTransaction();
        db->dropTable(getTableName(fmtPolicy, props));
        db->endTransaction();
    }
}


template void doTestBoost<IdVector>(std::string const &);
template void doTestBoost<MatchPairVector>(std::string const &);
template void doTestDb<IdVector>(std::string const &, std::string const &, std::string const &);
template void doTestDb<MatchPairVector>(std::string const &, std::string const &, std::string const &);


} // end of anonymous namespace


BOOST_AUTO_TEST_CASE(testBoost) {
    doTestBoost<IdVector>("IdVector");
    doTestBoost<IdPairVector>("IdPairVector");
    doTestBoost<MatchPairVector>("MatchPairVector");
}


BOOST_AUTO_TEST_CASE(testDb) {
    doTestDb<IdVector>("DbStorage", "Id", "IdVector");
    doTestDb<IdPairVector>("DbStorage", "IdPair", "IdPairVector");
    doTestDb<MatchPairVector>("DbStorage", "MatchPair", "MatchPairVector");
}


BOOST_AUTO_TEST_CASE(testDbTsv) {
    doTestDb<IdVector>("DbTsvStorage", "Id", "IdVector");
    doTestDb<IdPairVector>("DbTsvStorage", "IdPair", "IdPairVector");
    doTestDb<MatchPairVector>("DbTsvStorage", "MatchPair", "MatchPairVector");
}

