// -*- lsst-c++ -*-

/**
 * @file
 * @brief   Implementation of persistable association pipeline result vectors.
 *
 * @ingroup associate
 */

#include "lsst/ap/Results.h"


// -- PersistableMatchPairVector ----------------

lsst::ap::PersistableMatchPairVector::PersistableMatchPairVector() :
    lsst::daf::base::Citizen(typeid(*this)),
    _matchPairs()
{}


lsst::ap::PersistableMatchPairVector::PersistableMatchPairVector(lsst::ap::MatchPairVector const & matchPairs) :
    lsst::daf::base::Citizen(typeid(*this)),
    _matchPairs(matchPairs)
{}


lsst::ap::PersistableMatchPairVector::~PersistableMatchPairVector() {} 

bool lsst::ap::PersistableMatchPairVector::operator==(lsst::ap::MatchPairVector const & other) const {
    if (_matchPairs.size() != other.size()) {
        return false;
    }
    return std::equal(_matchPairs.begin(), _matchPairs.end(), other.begin());
}


// -- PersistableIdPairVector ----------------

lsst::ap::PersistableIdPairVector::PersistableIdPairVector() :
    lsst::daf::base::Citizen(typeid(*this)),
    _idPairs()
{}


lsst::ap::PersistableIdPairVector::PersistableIdPairVector(lsst::ap::IdPairVector const & idPairs) :
    lsst::daf::base::Citizen(typeid(*this)),
    _idPairs(idPairs)
{}


lsst::ap::PersistableIdPairVector::~PersistableIdPairVector() {}

bool lsst::ap::PersistableIdPairVector::operator==(lsst::ap::IdPairVector const & other) const {
    if (_idPairs.size() != other.size()) {
        return false;
    }
    return std::equal(_idPairs.begin(), _idPairs.end(), other.begin());
}


// -- IdVector ----------------

lsst::ap::PersistableIdVector::PersistableIdVector() :
    lsst::daf::base::Citizen(typeid(*this)),
    _ids()
{}

    
lsst::ap::PersistableIdVector::PersistableIdVector(lsst::ap::IdVector const & ids) :
    lsst::daf::base::Citizen(typeid(*this)),
    _ids(ids)
{}

    
lsst::ap::PersistableIdVector::~PersistableIdVector() {}

bool lsst::ap::PersistableIdVector::operator==(lsst::ap::IdVector const & other) const {
    if (_ids.size() != other.size()) {
        return false;
    }
    return std::equal(_ids.begin(), _ids.end(), other.begin());
}

