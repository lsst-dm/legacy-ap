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
 

/**
 * @file
 * @brief   Implementation of persistable association pipeline result vectors.
 *
 * @ingroup ap
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

