// -*- lsst-c++ -*-
/** @file
  * @brief High-level API for sources clustering and cluster attributes.
  *
  * @ingroup ap
  * @author Serge Monkewitz
  */
#ifndef LSST_AP_CLUSTER_SOURCECLUSTER_H
#define LSST_AP_CLUSTER_SOURCECLUSTER_H

#include <limits>
#include <vector>

#include "boost/shared_ptr.hpp"

#include "lsst/tr1/unordered_map.h"
#include "lsst/daf/base/Persistable.h"
#include "lsst/pex/policy.h"
#include "lsst/afw/detection/Source.h"

#include "../Common.h"
#include "Formatters.h"


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


/** Helper template for possibly unset/invalid floating
    point types. Nulls are represented as NaNs.
  */
template <typename FloatT>
class NullOr {
public:
    NullOr() {
        setNull();
    }
    NullOr(FloatT value) : _value(value) { }

    bool isNull() const {
        return isNaN(_value);
    }
    operator FloatT() const {
        return _value;
    }

    NullOr & operator=(NullOr const & nullable) {
        _value = nullable._value;
        return *this;
    }
    NullOr & operator=(FloatT const value) {
        _value = value;
        return *this;
    }
    void setNull() {
        _value = std::numeric_limits<FloatT>::quiet_NaN();
    }

private:
    FloatT _value;

    template <typename Archive> void serialize(Archive &, unsigned int const); 
    friend class boost::serialization::access;
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
    /** The first of NSAMPLE_BITS flag bits used to store the number of sources
      * used to determine the PSF flux sample mean.
      */
    static int const FLUX_NSAMPLE_OFF = 0;
    /** The first of NSAMPLE_BITS flag bits used to store the number of sources
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

    // attribute access
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
    NullOr<float> const & getFlux() const {
        return _flux;
    }
    NullOr<float> const & getFluxSigma() const {
        return _fluxSigma;
    }
    int getNumEllipticitySamples() const;
    NullOr<float> const & getE1() const {
        return _e1;
    }
    NullOr<float> const & getE1Sigma() const {
        return _e1Sigma;
    }
    NullOr<float> const & getE2() const {
        return _e2;
    }
    NullOr<float> const & getE2Sigma() const {
        return _e2Sigma;
    }
    NullOr<float> const & getRadius() const {
        return _radius;
    }
    NullOr<float> const & getRadiusSigma() const {
        return _radiusSigma;
    }
    bool operator==(PerFilterSourceClusterAttributes const & attributes) const;

    // attribute modification
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
    void setFlux(NullOr<float> const & flux,
                 NullOr<float> const & fluxSigma);
    void setNumEllipticitySamples(int samples);
    void setEllipticity();
    void setEllipticity(NullOr<float> const & e1,
                        NullOr<float> const & e2,
                        NullOr<float> const & radius,
                        NullOr<float> const & e1Sigma,
                        NullOr<float> const & e2Sigma,
                        NullOr<float> const & radiusSigma);

private:
    int _filterId;
    int _numObs;
    int _flags;
    double _earliestObsTime;
    double _latestObsTime;
    // flux
    NullOr<float> _flux;
    NullOr<float> _fluxSigma;
    // ellipticity parameters
    NullOr<float> _e1;
    NullOr<float> _e2;
    NullOr<float> _radius;
    NullOr<float> _e1Sigma;
    NullOr<float> _e2Sigma;
    NullOr<float> _radiusSigma;

    void computeFlux(lsst::afw::detection::SourceSet const & sources,
                     int fluxIgnoreMask);
    void computeEllipticity(lsst::afw::detection::SourceSet const & sources,
                            int ellipticityIgnoreMask);

    template <typename Archive> void serialize(Archive &, unsigned int const); 
    friend class boost::serialization::access;
};

inline bool operator!=(PerFilterSourceClusterAttributes const & lhs,
                       PerFilterSourceClusterAttributes const & rhs) {
    return !(lhs == rhs);
}


/** Attributes derived from a cluster of Sources; i.e. Sources that
  * have been determined to be observations of the same physical object
  * according to some clustering algorithm.
  */
class LSST_AP_API SourceClusterAttributes
{
public:
    typedef boost::shared_ptr<SourceClusterAttributes> Ptr;
    typedef boost::shared_ptr<SourceClusterAttributes const> ConstPtr;
    typedef std::tr1::unordered_map<int, PerFilterSourceClusterAttributes>
        PerFilterAttributesMap;

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

    // attribute access
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
    NullOr<float> const & getRaSigma() const {
        return _raSigma;
    }
    NullOr<float> const & getDecSigma() const {
        return _decSigma;
    }
    NullOr<float> const & getRaDecCov() const {
        return _raDecCov;
    }

    bool hasFilter(int filterId) const;
    std::vector<int> const getFilterIds() const;

    PerFilterAttributesMap const & getPerFilterAttributes() const {
        return _perFilterAttributes;
    }

    PerFilterSourceClusterAttributes const & getPerFilterAttributes(
        int filterId) const;

    bool operator==(SourceClusterAttributes const & attributes) const;

    // attribute modification
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
                     NullOr<float> const & raSigma,
                     NullOr<float> const & decSigma,
                     NullOr<float> const & raDecCov);
    void clearPerFilterAttributes();
    bool setPerFilterAttributes(
        PerFilterSourceClusterAttributes const & attributes);
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
    NullOr<float> _raSigma;
    NullOr<float> _decSigma;
    NullOr<float> _raDecCov;
    // per filter propertiesa
    PerFilterAttributesMap _perFilterAttributes;
 
    void computePosition(lsst::afw::detection::SourceSet const & sources);

    template <typename Archive> void serialize(Archive &, unsigned int const);
    friend class boost::serialization::access;
};

inline bool operator!=(SourceClusterAttributes const & lhs,
                       SourceClusterAttributes const & rhs) {
    return !(lhs == rhs);
}


typedef std::vector<SourceClusterAttributes::Ptr> SourceClusterVector;

/** A lsst::daf::base::Persistable wrapper for vectors of shared pointers to
  * SourceClusterAttributes objects.
  */
class LSST_AP_API PersistableSourceClusterVector :
    public lsst::daf::base::Persistable
{
public:
    typedef boost::shared_ptr<PersistableSourceClusterVector> Ptr;

    PersistableSourceClusterVector() : _clusters() { }
    PersistableSourceClusterVector(SourceClusterVector const & clusters) :
        _clusters(clusters)
    { }

    virtual ~PersistableSourceClusterVector() { }

    size_t size() const {
        return _clusters.size();
    }
    SourceClusterVector const & getClusters() const {
        return _clusters;
    }

    void setClusters(SourceClusterVector const & clusters) {
        _clusters = clusters;
    }

private:
    LSST_PERSIST_FORMATTER(lsst::ap::cluster::SourceClusterVectorFormatter);
    SourceClusterVector _clusters;
};

}}} // namespace lsst:ap::cluster

#endif // LSST_AP_CLUSTER_SOURCECLUSTER_H
