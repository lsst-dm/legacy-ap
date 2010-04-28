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

#include <math.h>
#include <limits>
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

LSST_AP_API void segregateInvalidSources(
    lsst::afw::detection::SourceSet & sources,
    lsst::afw::detection::SourceSet & badSources);

LSST_AP_API void segregateBadSources(
    lsst::afw::detection::SourceSet & sources,
    lsst::afw::detection::SourceSet & badSources,
    int badSourceMask);

LSST_AP_API void updateBadSources(
     lsst::afw::detection::SourceSet & badSources);


inline bool isNaN(double x) {
#if LSST_AP_HAVE_ISNAND
    return isnand(x);
#else
    return x != x;
#endif
}

inline bool isNaN(float x) {
#if LSST_AP_HAVE_ISNANF
    return isnanf(x);
#else
    return x != x;
#endif
}

/** Helper template for nullable floating point types. Nulls are equivalent
  */
template <typename FloatT>
class Nullable {
public:
    Nullable() {
        setNull();
    }
    Nullable(FloatT value) : _value(value) { }

    bool isNull() const {
        return isNaN(_value);
    }
    operator FloatT () const {
        return _value;
    }

    Nullable & operator=(Nullable const & n) {
        _value = n._value;
        return *this;
    }
    Nullable & operator=(FloatT const value) {
        _value = value;
        return *this;
    }
    void setNull() {
        _value = std::numeric_limits<FloatT>::quiet_NaN();
    }

private:
    FloatT _value;
};


/** Per filter source cluster attributes.
  */
class LSST_AP_API PerFilterSourceClusterAttributes {
public:
    typedef boost::shared_ptr<PerFilterSourceClusterAttributes> Ptr;
    typedef boost::shared_ptr<PerFilterSourceClusterAttributes const> ConstPtr;

    /** Number of flag bits used to store counts of the number of samples 
      * (sources) used to determine the flux and ellipticity parameter sample
      * means.
      */
    static int const NSAMPLE_BITS = 12;
    static int const NSAMPLE_MASK = (1 << NSAMPLE_BITS) - 1;
    /** The first of N_BITS flag bits used to store the number of sources
      * used to determine the PSF flux sample mean.
      */
    static int const FLUX_NSAMPLE_OFF = 0;
    /** The first of N_BITS flag bits used to store the number of sources
      * used to determine the ellipticity parameter sample means.
      */
    static int const ELLIPTICITY_NSAMPLE_OFF = FLUX_NSAMPLE_OFF + NSAMPLE_BITS;

    PerFilterSourceClusterAttributes();
    PerFilterSourceClusterAttributes(
        lsst::afw::detection::Source const & source,
        int fluxIgnoreMask,
        int ellipticityIgnoreMask);
    PerFilterSourceClusterAttributes(
        lsst::afw::detection::SourceSet const & sources,
        int fluxIgnoreMask,
        int ellipticityIgnoreMask);
    ~PerFilterSourceClusterAttributes();

    int getFilterId() const {
        return _filterId;
    }
    int getNumObs() const {
        return _numObs;
    }
    int getFlags() const {
        return _flags;
    }
    double getEarliestObsTime() const {
        return _earliestObsTime;
    }
    double getLatestObsTime() const {
        return _latestObsTime;
    }
    int getNumFluxSamples() const;
    Nullable<float> const & getFlux() const {
        return _flux;
    }
    Nullable<float> const & getFluxSigma() const {
        return _fluxSigma;
    }
    int getNumEllipticitySamples() const;
    Nullable<float> const & getE1() const {
        return _e1;
    }
    Nullable<float> const & getE1Sigma() const {
        return _e1Sigma;
    }
    Nullable<float> const & getE2() const {
        return _e2;
    }
    Nullable<float> const & getE2Sigma() const {
        return _e2Sigma;
    }
    Nullable<float> const & getRadius() const {
        return _radius;
    }
    Nullable<float> const & getRadiusSigma() const {
        return _radiusSigma;
    }

    void setFilterId(int filterId) {
        _filterId = filterId;
    }
    void setNumObs(int numObs) {
        _numObs = numObs;
    }
    void setFlags(int flags) {
        _flags = flags;
    }
    void setObsTimeRange(double earliest, double latest);
    void setNumFluxSamples(int samples);
    void setFlux(Nullable<float> const & flux,
                 Nullable<float> const & fluxSigma);
    void setNumEllipticitySamples(int samples);
    void setEllipticity();
    void setEllipticity(Nullable<float> const & e1,
                        Nullable<float> const & e2,
                        Nullable<float> const & radius,
                        Nullable<float> const & e1Sigma,
                        Nullable<float> const & e2Sigma,
                        Nullable<float> const & radiusSigma);

private:
    int _filterId;
    int _numObs;
    int _flags;
    double _earliestObsTime;
    double _latestObsTime;
    // flux
    Nullable<float> _flux;
    Nullable<float> _fluxSigma;
    // ellipticities
    Nullable<float> _e1;
    Nullable<float> _e2;
    Nullable<float> _radius;
    Nullable<float> _e1Sigma;
    Nullable<float> _e2Sigma;
    Nullable<float> _radiusSigma;

    void computeFlux(lsst::afw::detection::SourceSet const & sources,
                     int fluxIgnoreMask);
    void computeEllipticity(lsst::afw::detection::SourceSet const & sources,
                            int ellipticityIgnoreMask);
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

    /** Source cluster flag bits.
      */
    enum Flags {
        /** The cluster is a noise cluster, i.e. is a cluster generated for a
          * single noise source that is not a border source. Noise sources are
          * sources beneath the OPTICS core-source density threshold; border
          * sources are noise sources that have nevertheless been assigned to
          * a cluster.
          */
        NOISE = 0x1,
        /** The cluster was generated from a single bad source (where bad
          * sources are identified by a policy determined set of detection
          * flag bits).
          */
        BAD = 0x2,
    };

    SourceClusterAttributes();
    SourceClusterAttributes(lsst::afw::detection::Source const & source,
                            int64_t id,
                            int fluxIgnoreMask,
                            int ellipticityIgnoreMask);
    SourceClusterAttributes(lsst::afw::detection::SourceSet const & sources,
                            int64_t id,
                            int fluxIgnoreMask,
                            int ellipticityIgnoreMask);
    ~SourceClusterAttributes();

    int64_t getClusterId() const {
        return _clusterId;
    }
    int getNumObs() const {
        return _numObs;
    }
    int getFlags() const {
        return _flags;
    }
    double getEarliestObsTime() const {
        return _earliestObsTime;
    }
    double getLatestObsTime() const {
        return _latestObsTime;
    }
    double getRa() const {
        return _ra;
    }
    double getDec() const {
        return _dec;
    }
    Nullable<float> const & getRaSigma() const {
        return _raSigma;
    }
    Nullable<float> const & getDecSigma() const {
        return _decSigma;
    }
    Nullable<float> const & getRaDecCov() const {
        return _raDecCov;
    }

    bool hasFilter(int filterId) const;
    std::vector<int> const getFilterIds() const;
    PerFilterAttributeMap const & getPerFilterAttributes() const {
        return _perFilterAttributes;
    }
    PerFilterSourceClusterAttributes const & getPerFilterAttributes(
        int filterId) const;

    void setClusterId(int clusterId) {
        _clusterId = clusterId;
    }
    void setNumObs(int numObs) {
        _numObs = numObs;
    }
    void setFlags(int flags) {
        _flags = flags;
    }
    void setObsTimeRange(double earliest, double latest);
    void setPosition(double ra,
                     double dec,
                     Nullable<float> const & raSigma,
                     Nullable<float> const & decSigma,
                     Nullable<float> const & raDecCov);
    void clearPerFilterAttributes();
    bool setPerFilterAttributes(PerFilterSourceClusterAttributes const & attributes);
    bool removePerFilterAttributes(int filterId);

private:
    int64_t _clusterId;
    int _numObs;
    int _flags;
    double _earliestObsTime;
    double _latestObsTime;
    // position
    double _ra;
    double _dec;
    Nullable<float> _raSigma;
    Nullable<float> _decSigma;
    Nullable<float> _raDecCov;
    // per filter propertiesa
    PerFilterAttributeMap _perFilterAttributes;
 
    void computePosition(lsst::afw::detection::SourceSet const & sources);
};

}}} // namespace lsst:ap::cluster

#endif // LSST_AP_CLUSTER_SOURCECLUSTER_H
