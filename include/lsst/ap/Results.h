// -*- lsst-c++ -*-

/**
 * @file
 * @brief   Persistable vectors for association pipeline results.
 *
 * @ingroup ap
 */

#ifndef LSST_AP_RESULTS_H
#define LSST_AP_RESULTS_H

#include <utility>
#include <vector>

#include "boost/shared_ptr.hpp"

#include "lsst/daf/base/Citizen.h"
#include "lsst/daf/base/Persistable.h"

#include "Common.h"


/// @cond
namespace boost { namespace serialization {
    class access;
}}
/// @endcond


namespace lsst { namespace ap {

namespace io {
    // forward declarations for formatters
    class MatchPairVectorFormatter;
    class IdPairVectorFormatter;
    class IdVectorFormatter;
}


/** @brief  Holds a pair of ids and the distance between the corresponding positions on the unit sphere. */
class LSST_AP_API MatchPair {

public :

    MatchPair() : _first(-1), _second(-1), _distance(0.0) {}

    MatchPair(boost::int64_t const first, boost::int64_t const second, double const distance) :
        _first(first), _second(second), _distance(distance) {}

    boost::int64_t getFirst() const {
        return _first;
    }
    boost::int64_t getSecond() const {
        return _second;
    }
    double getDistance() const {
        return _distance;
    }

    void setFirst(boost::int64_t const first) {
        _first = first;
    }
    void setSecond(boost::int64_t const second) {
        _second = second;
    }
    void setDistance(double const distance) {
        _distance = distance;
    }

    bool operator==(MatchPair const & mp) const {
        return _first == mp.getFirst() && _second == mp.getSecond();
    }

    bool operator!=(MatchPair const & mp) const {
        return !operator==(mp);
    }

private :

    boost::int64_t _first;
    boost::int64_t _second;
    double _distance;

    template <typename Archive>
    void serialize(Archive & ar, unsigned int const version) {
        ar & _first;
        ar & _second;
        ar & _distance;
    }

    friend class boost::serialization::access;
    friend class io::MatchPairVectorFormatter;
};


/** @brief  Holds a pair of ids. */
typedef std::pair<boost::int64_t, boost::int64_t> IdPair;
/** @brief  A list of MatchPair instances. */
typedef std::vector<MatchPair> MatchPairVector;
/** @brief  A list of IdPair instances. */
typedef std::vector<IdPair> IdPairVector;
/** @brief  A list of integer ids. */
typedef std::vector<boost::int64_t> IdVector;


/** @brief  A persistable wrapper for a MatchPairVector. */
class LSST_AP_API PersistableMatchPairVector :
    public lsst::daf::base::Persistable,
    public lsst::daf::base::Citizen
{
public :
    typedef boost::shared_ptr<PersistableMatchPairVector> Ptr;

    PersistableMatchPairVector();
    explicit PersistableMatchPairVector(MatchPairVector const &);
    ~PersistableMatchPairVector();

    MatchPairVector & getMatchPairs() {
        return _matchPairs;
    }
    MatchPairVector const & getMatchPairs() const {
        return _matchPairs;
    }
    void setMatchPairs(MatchPairVector const & matchPairs) {
        _matchPairs = matchPairs;
    }

    bool operator==(MatchPairVector const & other) const;
    bool operator==(PersistableMatchPairVector const & other) const {
        return other == _matchPairs;
    }

private:
    LSST_PERSIST_FORMATTER(lsst::ap::io::MatchPairVectorFormatter);
    MatchPairVector _matchPairs;
};


/** @brief  A persistable wrapper for an IdPairVector. */
class LSST_AP_API PersistableIdPairVector :
    public lsst::daf::base::Persistable,
    public lsst::daf::base::Citizen
{
public :
    typedef boost::shared_ptr<PersistableIdPairVector> Ptr;

    PersistableIdPairVector();
    explicit PersistableIdPairVector(IdPairVector const &);
    ~PersistableIdPairVector();

    IdPairVector & getIdPairs() {
        return _idPairs;
    }
    IdPairVector const & getIdPairs() const {
        return _idPairs;
    }
    void setIdPairs(IdPairVector const & idPairs) {
        _idPairs = idPairs;
    }

    bool operator==(IdPairVector const & other) const;
    bool operator==(PersistableIdPairVector const & other) const {
        return other == _idPairs;
    }

private:
    LSST_PERSIST_FORMATTER(lsst::ap::io::IdPairVectorFormatter);
    IdPairVector _idPairs;
};


/** @brief  A persistable wrapper for an IdVector. */
class LSST_AP_API PersistableIdVector :
    public lsst::daf::base::Persistable,
    public lsst::daf::base::Citizen
{
public :
    typedef boost::shared_ptr<PersistableIdVector> Ptr;

    PersistableIdVector();
    explicit PersistableIdVector(IdVector const &);
    ~PersistableIdVector();

    IdVector & getIds() {
        return _ids;
    }
    IdVector const & getIds() const {
        return _ids;
    }
    void setIds(IdVector const & ids) {
        _ids = ids;
    }

    bool operator==(IdVector const & other) const;
    bool operator==(PersistableIdVector const & other) const {
        return other == _ids;
    }

private:
    LSST_PERSIST_FORMATTER(lsst::ap::io::IdVectorFormatter);
    IdVector _ids;
};


}}  // end of namespace lsst::ap

#endif // LSST_AP_RESULTS_H

