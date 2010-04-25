// -*- lsst-c++ -*-
/** @file
  * @brief Functions for clustering points (using the OPTICS algorithm)
  *        and for computing cluster attributes.
  *
  * @ingroup ap
  * @author Serge Monkewitz
  */
#ifndef LSST_AP_CLUSTER_SOURCECLUSTER_H
#define LSST_AP_CLUSTER_SOURCECLUSTER_H

#include <vector>

#include "boost/shared_ptr.hpp"

#include "lsst/tr1/unordered_map.h"
#include "lsst/daf/base/Persistable.h"
#include "lsst/pex/policy.h"
#include "lsst/afw/detection/Source.h"

#include "../Common.h"


namespace lsst { namespace ap { namespace cluster {

class SourceClusterAttributes;

LSST_AP_API std::vector<lsst::afw::detection::SourceSet> const cluster(
    lsst::afw::detection::SourceSet const & sources,
    lsst::pex::policy::Policy::Ptr policy);

LSST_AP_API void updateSources(
    SourceClusterAttributes const & cluster,
    lsst::afw::detection::SourceSet & sources);


/** Helper template for nullable scalars.
  */
template <typename ScalarT>
class Nullable {
public:
    Nullable() : _value(), _null(true) { }
    Nullable(ScalarT value) : _value(value), _null(false) { }

    bool isNull() const { return _null; }
    ScalarT get() const { return _value; }
    ScalarT operator()() const { return _value; }

    void setNull() {
        _null = true;
    }
    void set(ScalarT value) {
        _value = value;
        _null = false;
    }

private:
    ScalarT _value;
    bool _null;
};

// comparison operators
template <typename ScalarT>
inline bool operator==(Nullable<ScalarT> const & nullable, ScalarT value) {
    return !nullable.isNull() && nullable() == value;
}
template <typename ScalarT>
inline bool operator==(Nullable<ScalarT> const & n1, Nullable<ScalarT> const & n2) {
    return n1.isNull() ? n2.isNull() : n2 == n1();
}
template <typename ScalarT>
inline bool operator!=(Nullable<ScalarT> const & nullable, ScalarT value) {
    return !(nullable == value);
}
template <typename ScalarT>
inline bool operator!=(Nullable<ScalarT> const & n1, Nullable<ScalarT> const & n2) {
    return !(n1 == n2); 
}
template <typename ScalarT>
inline bool operator==(ScalarT value, Nullable<ScalarT> const & nullable) {
    return nullable == value;
}
template <typename ScalarT>
inline bool operator!=(ScalarT value, Nullable<ScalarT> const & nullable) {
    return nullable != value;
}


/** Per filter source cluster attributes.
  */
class LSST_AP_API PerFilterSourceClusterAttributes {
public:
    typedef boost::shared_ptr<PerFilterSourceClusterAttributes> Ptr;
    typedef boost::shared_ptr<PerFilterSourceClusterAttributes const> ConstPtr;

    PerFilterSourceClusterAttributes();
    PerFilterSourceClusterAttributes(lsst::afw::detection::SourceSet const & sources);
    ~PerFilterSourceClusterAttributes();

    int getFilterId() const { return _filterId; }
    int getNumObs() const { return _numObs; }
    double getEarliestObsTime() const { return _earliestObsTime; }
    double getLatestObsTime() const { return _latestObsTime; }
    double getFlux() const { return _flux; }
    double getFluxSigma() const { return _fluxSigma; }
    Nullable<double> const & getE1() const { return _e1; }
    Nullable<double> const & getE1Sigma() const { return _e1Sigma; }
    Nullable<double> const & getE2() const { return _e2; }
    Nullable<double> const & getE2Sigma() const { return _e2Sigma; }
    Nullable<double> const & getRadius() const { return _radius; }
    Nullable<double> const & getRadiusSigma() const { return _radiusSigma; }

    void setFilterId(int filterId) { _filterId = filterId; }
    void setNumObs(int numObs) { _numObs = numObs; }
    void setObsTimeRange(double earliest, double latest);
    void setFlux(double flux, double fluxSigma);
    void setEllipticities();
    void setEllipticities(Nullable<double> const & e1,
                          Nullable<double> const & e1Sigma,
                          Nullable<double> const & e2,
                          Nullable<double> const & e2Sigma,
                          Nullable<double> const & radius,
                          Nullable<double> const & radiusSigma);

private:
    int _filterId;
    int _numObs;
    double _earliestObsTime;
    double _latestObsTime;
    // flux
    double _flux;
    double _fluxSigma;
    // ellipticities
    Nullable<double> _e1;
    Nullable<double> _e1Sigma;
    Nullable<double> _e2;
    Nullable<double> _e2Sigma;
    Nullable<double> _radius;
    Nullable<double> _radiusSigma;

    void computeFlux(lsst::afw::detection::SourceSet const & sources);
    void computeEllipticities(lsst::afw::detection::SourceSet const & sources);
};


/** Attributes derived from a cluster of Sources; i.e. Sources that
  * have been determined to be observations of the same physical object
  * according to some clustering algorithm.
  */
class LSST_AP_API SourceClusterAttributes :
    public lsst::daf::base::Persistable
{
public:
    typedef boost::shared_ptr<SourceClusterAttributes> Ptr;
    typedef boost::shared_ptr<SourceClusterAttributes const> ConstPtr;
    typedef std::tr1::unordered_map<int, PerFilterSourceClusterAttributes>
        PerFilterAttributeMap;

    SourceClusterAttributes();
    SourceClusterAttributes(lsst::afw::detection::SourceSet const & sources,
                             int64_t id);
    ~SourceClusterAttributes();

    int64_t getClusterId() const { return _clusterId; }
    int getNumObs() const { return _numObs; }
    double getEarliestObsTime() const { return _earliestObsTime; }
    double getLatestObsTime() const { return _latestObsTime; }
    double getRa() const { return _ra; }
    double getDec() const { return _dec; }
    double getRaSigma() const { return _raSigma; }
    double getDecSigma() const { return _decSigma; }
    Nullable<double> const & getRaDecCov() const { return _raDecCov; }

    bool hasFilter(int filterId) const;
    std::vector<int> const getFilterIds() const;
    PerFilterAttributeMap const & getPerFilterAttributes() const {
        return _perFilterAttributes;
    }
    PerFilterSourceClusterAttributes const & getPerFilterAttributes(
        int filterId) const;

    void setClusterId(int clusterId) { _clusterId = clusterId; }
    void setNumObs(int numObs) { _numObs = numObs; }
    void setObsTimeRange(double earliest, double latest);
    void setPosition(double ra,
                     double dec,
                     double raSigma,
                     double decSigma,
                     Nullable<double> const & raDecCov);

    void clearPerFilterAttributes();
    bool setPerFilterAttributes(PerFilterSourceClusterAttributes const & attributes);
    bool removePerFilterAttributes(int filterId);

private:
    int64_t _clusterId;
    int _numObs;
    double _earliestObsTime;
    double _latestObsTime;
    // position
    double _ra;
    double _dec;
    double _raSigma;
    double _decSigma;
    Nullable<double> _raDecCov;
    // per filter propertiesa
    PerFilterAttributeMap _perFilterAttributes;
 
    void computePosition(lsst::afw::detection::SourceSet const & sources);
};

}}} // namespace lsst:ap::cluster

#endif // LSST_AP_CLUSTER_SOURCECLUSTER_H
