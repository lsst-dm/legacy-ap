// -*- lsst-c++ -*-

/**
 * @file
 * @brief   Implementation of persistable association pipeline result vectors.
 *
 * @ingroup associate
 */

#include <lsst/ap/Results.h>


namespace lsst {
namespace ap {


// -- MatchPairVector ----------------

MatchPairVector::MatchPairVector() :
    lsst::mwi::data::Citizen(typeid(*this)),
    _vec()
{}


MatchPairVector::MatchPairVector(size_type n) :
    lsst::mwi::data::Citizen(typeid(*this)),
    _vec(n)
{}


MatchPairVector::MatchPairVector(size_type n, value_type const & val) :
    lsst::mwi::data::Citizen(typeid(*this)),
    _vec(n, val)
{}


MatchPairVector::~MatchPairVector() {}


MatchPairVector::MatchPairVector(MatchPairVector const & v) :
    lsst::mwi::data::Citizen(typeid(*this)),
    _vec(v._vec)
{}


MatchPairVector::MatchPairVector(Vector const & v) :
    lsst::mwi::data::Citizen(typeid(*this)),
    _vec(v)
{}


MatchPairVector & MatchPairVector::operator=(MatchPairVector const & v) {
    if (this != &v) {
        _vec = v._vec;
    }
    return *this;
}

MatchPairVector & MatchPairVector::operator=(Vector const & v) {
    _vec = v;
    return *this;
}


// -- IdPairVector ----------------

IdPairVector::IdPairVector() :
    lsst::mwi::data::Citizen(typeid(*this)),
    _vec()
{}


IdPairVector::IdPairVector(size_type n) :
    lsst::mwi::data::Citizen(typeid(*this)),
    _vec(n)
{}


IdPairVector::IdPairVector(size_type n, value_type const & val) :
    lsst::mwi::data::Citizen(typeid(*this)),
    _vec(n, val)
{}


IdPairVector::~IdPairVector() {}


IdPairVector::IdPairVector(IdPairVector const & v) :
    lsst::mwi::data::Citizen(typeid(*this)),
    _vec(v._vec)
{}


IdPairVector::IdPairVector(Vector const & v) :
    lsst::mwi::data::Citizen(typeid(*this)),
    _vec(v)
{}


IdPairVector & IdPairVector::operator=(IdPairVector const & v) {
    if (this != &v) {
        _vec = v._vec;
    }
    return *this;
}

IdPairVector & IdPairVector::operator=(Vector const & v) {
    _vec = v;
    return *this;
}


// -- IdVector ----------------

IdVector::IdVector() :
    lsst::mwi::data::Citizen(typeid(*this)),
    _vec()
{}


IdVector::IdVector(size_type n) :
    lsst::mwi::data::Citizen(typeid(*this)),
    _vec(n)
{}


IdVector::IdVector(size_type n, value_type const & val) :
    lsst::mwi::data::Citizen(typeid(*this)),
    _vec(n, val)
{}


IdVector::~IdVector() {}


IdVector::IdVector(IdVector const & v) :
    lsst::mwi::data::Citizen(typeid(*this)),
    _vec(v._vec)
{}


IdVector::IdVector(Vector const & v) :
    lsst::mwi::data::Citizen(typeid(*this)),
    _vec(v)
{}


IdVector & IdVector::operator=(IdVector const & v) {
    if (this != &v) {
        _vec = v._vec;
    }
    return *this;
}

IdVector & IdVector::operator=(Vector const & v) {
    _vec = v;
    return *this;
}


}} // end of namespace lsst::ap

