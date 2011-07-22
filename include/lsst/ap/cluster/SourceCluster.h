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
#include "lsst/utils/ieee.h"
#include "lsst/daf/base/Persistable.h"
#include "lsst/pex/policy.h"
#include "lsst/afw/geom/AffineTransform.h"
#include "lsst/afw/detection/Source.h"

#include "../Common.h"
#include "../match/ExposureInfo.h"
#include "Formatters.h"


namespace lsst { namespace ap { namespace cluster {

class SourceClusterAttributes;

LSST_AP_API std::vector<lsst::afw::detection::SourceSet> const cluster(
    lsst::afw::detection::SourceSet const & sources,
    lsst::pex::policy::Policy::Ptr policy);

LSST_AP_API void updateSources(
    SourceClusterAttributes const & cluster,
    lsst::afw::detection::SourceSet & sources);

LSST_AP_API void locateAndFilterSources(
    lsst::afw::detection::SourceSet & sources,
    lsst::afw::detection::SourceSet & invalidSources,
    lsst::ap::match::ExposureInfoMap const & exposures);

LSST_AP_API void segregateInvalidSources(
    lsst::afw::detection::SourceSet & sources,
    lsst::afw::detection::SourceSet & invalidSources);

LSST_AP_API void segregateBadSources(
    lsst::afw::detection::SourceSet & sources,
    lsst::afw::detection::SourceSet & badSources,
    int badSourceMask);

LSST_AP_API void updateUnclusteredSources(
     lsst::afw::detection::SourceSet & badSources);


/** Helper template for possibly unset/invalid floating
  * point types. Nulls are represented as NaNs.
  */
template <typename FloatT>
class NullOr {
public:
    NullOr() {
        setNull();
    }
    NullOr(FloatT value) : _value(value) { }

    bool isNull() const {
        return lsst::utils::isnan(_value);
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

    bool operator==(NullOr const & n) const {
        return n == _value;
    }
    bool operator!=(NullOr const & n) const {
        return !(n == _value);
    }
    bool operator==(FloatT const & value) const {
        return isNull() ? lsst::utils::isnan(value) : value == _value;
    }
    bool operator!=(FloatT const & value) const {
        return isNull() ? !lsst::utils::isnan(value) : value != _value;
    }

private:
    FloatT _value;

    template <typename Archive> void serialize(Archive &, unsigned int const); 
    friend class boost::serialization::access;
};

template <typename FloatT>
inline bool operator==(FloatT const & value, NullOr<FloatT> const & n) {
    return n == value;
}
template <typename FloatT>
inline bool operator!=(FloatT const & value, NullOr<FloatT> const & n) {
    return n != value;
}


/** Helper class that bundles a source, it's originating exposure, and an
  * affine transform from the pixel space of that exposure to the pixel space
  * of a tangent plane projection centered at the fiducial cluster position
  * and with the standard north, east basis.
  */
class LSST_AP_API SourceAndExposure {
public:
    SourceAndExposure(
        boost::shared_ptr<lsst::afw::detection::Source const> source,
        lsst::ap::match::ExposureInfo::ConstPtr exposure,
        lsst::afw::geom::AffineTransform const & transform) :
        _source(source), _exposure(exposure), _transform(transform) { }

    ~SourceAndExposure() { }

    boost::shared_ptr<lsst::afw::detection::Source const> getSource() const {
        return _source;
    }
    lsst::ap::match::ExposureInfo::ConstPtr getExposureInfo() const {
        return _exposure;
    }
    lsst::afw::geom::AffineTransform const & getTransform() const {
        return _transform;
    }

private:
    boost::shared_ptr<lsst::afw::detection::Source const> _source;
    lsst::ap::match::ExposureInfo::ConstPtr _exposure;
    lsst::afw::geom::AffineTransform _transform;
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
    static int const NSAMPLE_BITS = 8;
    /** Integer with NSAMPLE_BITS least significant bits set to 1.
      */
    static int const NSAMPLE_MASK = (1 << NSAMPLE_BITS) - 1;
    /** The first of NSAMPLE_BITS flag bits used to store the number of sources
      * used to determine the point source (PSF) flux sample mean.
      */
    static int const FLUX_PS_NSAMPLE_OFF = 0;
    /** The first of NSAMPLE_BITS flag bits used to store the number of sources
      * used to determine the experimental small galaxy model flux sample mean.
      */
    static int const FLUX_SG_NSAMPLE_OFF = FLUX_PS_NSAMPLE_OFF + NSAMPLE_BITS;
    /** The first of NSAMPLE_BITS flag bits used to store the number of sources
      * used to determine the elliptical gaussian model flux sample mean.
      */
    static int const FLUX_GAUSS_NSAMPLE_OFF = FLUX_SG_NSAMPLE_OFF +
                                              NSAMPLE_BITS;
    /** The first of NSAMPLE_BITS flag bits used to store the number of sources
      * used to determine the ellipticity parameter sample means.
      */
    static int const ELLIPTICITY_NSAMPLE_OFF = FLUX_GAUSS_NSAMPLE_OFF +
                                               NSAMPLE_BITS;

    PerFilterSourceClusterAttributes();
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

    int getNumPsFluxSamples() const;
    int getNumSgFluxSamples() const;
    int getNumGaussianFluxSamples() const;

    NullOr<float> const & getPsFlux() const {
        return _psFlux;
    }
    NullOr<float> const & getPsFluxSigma() const {
        return _psFluxSigma;
    }
    NullOr<float> const & getSgFlux() const {
        return _sgFlux;
    }
    NullOr<float> const & getSgFluxSigma() const {
        return _sgFluxSigma;
    }
    NullOr<float> const & getGaussianFlux() const {
        return _gaussianFlux;
    }
    NullOr<float> const & getGaussianFluxSigma() const {
        return _gaussianFluxSigma;
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

    void setNumPsFluxSamples(int samples);
    void setNumSgFluxSamples(int samples);
    void setNumGaussianFluxSamples(int samples);
    void setPsFlux(NullOr<float> const & flux,
                   NullOr<float> const & fluxSigma);
    void setSgFlux(NullOr<float> const & flux,
                   NullOr<float> const & fluxSigma);
    void setGaussianFlux(NullOr<float> const & flux,
                         NullOr<float> const & fluxSigma);

    void setNumEllipticitySamples(int samples);
    void setEllipticity();
    void setEllipticity(NullOr<float> const & e1,
                        NullOr<float> const & e2,
                        NullOr<float> const & radius,
                        NullOr<float> const & e1Sigma,
                        NullOr<float> const & e2Sigma,
                        NullOr<float> const & radiusSigma);

    void computeFlux(SourceAndExposure const & source,
                     double fluxScale,
                     int psFluxIgnoreMask,
                     int sgFluxIgnoreMask,
                     int gaussianFluxIgnoreMask);

    void computeFlux(std::vector<SourceAndExposure> const & sources,
                     double fluxScale,
                     int psFluxIgnoreMask,
                     int sgFluxIgnoreMask,
                     int gaussianFluxIgnoreMask);

    void computeEllipticity(SourceAndExposure const & source,
                            int ellipticityIgnoreMask);

    void computeEllipticity(std::vector<SourceAndExposure> const & sources,
                            int ellipticityIgnoreMask);

private:
    int _filterId;
    int _numObs;
    int _flags;
    double _earliestObsTime;
    double _latestObsTime;
    // fluxes
    NullOr<float> _psFlux;
    NullOr<float> _psFluxSigma;
    NullOr<float> _sgFlux;
    NullOr<float> _sgFluxSigma;
    NullOr<float> _gaussianFlux;
    NullOr<float> _gaussianFluxSigma;
    // ellipticity parameters
    NullOr<float> _e1;
    NullOr<float> _e2;
    NullOr<float> _radius;
    NullOr<float> _e1Sigma;
    NullOr<float> _e2Sigma;
    NullOr<float> _radiusSigma;

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
    SourceClusterAttributes(int64_t id);
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
    double getMeanObsTime() const {
        return _meanObsTime;
    }

    double getRaPs() const {
        return _raPs;
    }
    double getDecPs() const {
        return _decPs;
    }
    NullOr<float> const & getRaPsSigma() const {
        return _raPsSigma;
    }
    NullOr<float> const & getDecPsSigma() const {
        return _decPsSigma;
    }
    NullOr<float> const & getRaDecPsCov() const {
        return _raDecPsCov;
    }

    NullOr<double> const & getRaSg() const {
        return _raSg;
    }
    NullOr<double> const & getDecSg() const {
        return _decSg;
    }
    NullOr<float> const & getRaSgSigma() const {
        return _raSgSigma;
    }
    NullOr<float> const & getDecSgSigma() const {
        return _decSgSigma;
    }
    NullOr<float> const & getRaDecSgCov() const {
        return _raDecSgCov;
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
    void setObsTime(double earliest, double latest, double mean);
    void setPsPosition(double ra,
                       double dec,
                       NullOr<float> const & raSigma,
                       NullOr<float> const & decSigma,
                       NullOr<float> const & raDecCov);
    void setSgPosition(NullOr<double> const & ra,
                       NullOr<double> const & dec,
                       NullOr<float> const & raSigma,
                       NullOr<float> const & decSigma,
                       NullOr<float> const & raDecCov);
    void clearPerFilterAttributes();
    bool setPerFilterAttributes(
        PerFilterSourceClusterAttributes const & attributes);
    bool removePerFilterAttributes(int filterId);

    void computeAttributes(
        boost::shared_ptr<lsst::afw::detection::Source const> source,
        lsst::ap::match::ExposureInfoMap const & exposures,
        double fluxScale,
        int psFluxIgnoreMask,
        int sgFluxIgnoreMask,
        int gaussianFluxIgnoreMask,
        int ellipticityIgnoreMask);

    void computeAttributes(
        lsst::afw::detection::SourceSet const & sources,
        lsst::ap::match::ExposureInfoMap const & exposures,
        double fluxScale,
        int psFluxIgnoreMask,
        int sgFluxIgnoreMask,
        int gaussianFluxIgnoreMask,
        int ellipticityIgnoreMask);

private:
    int64_t _clusterId;
    int _numObs;
    int _flags;
    double _earliestObsTime;
    double _latestObsTime;
    double _meanObsTime;

    // position, point source model (unweighted mean)
    double _raPs;
    double _decPs;
    NullOr<float> _raPsSigma;
    NullOr<float> _decPsSigma;
    NullOr<float> _raDecPsCov;
    // position, small galaxy model (inverse variance weighted mean)
    NullOr<double> _raSg;
    NullOr<double> _decSg;
    NullOr<float> _raSgSigma;
    NullOr<float> _decSgSigma;
    NullOr<float> _raDecSgCov;

    // per filter propertiesa
    PerFilterAttributesMap _perFilterAttributes;
 
    void _computePsPosition(lsst::afw::detection::SourceSet const & sources);

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
