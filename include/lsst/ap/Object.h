// -*- lsst-c++ -*-

/**
 * @file
 * @brief   C++ representations of an LSST Object.
 *
 * @ingroup associate
 */

#ifndef LSST_AP_OBJECT_H
#define LSST_AP_OBJECT_H

#include <lsst/daf/base/DateTime.h>

#include <lsst/afw/image/Filter.h>

#include "Common.h"
#include "Bitset.h"


namespace lsst {
namespace ap {


/**
 * @brief   A simplified Object representation containing only id, position and
 *          per-filter variability probabilities.
 *
 * This is sufficient for performing spatial crosss matches. None of these fields may ever be NULL.
 */
class LSST_AP_API SimpleObject {

public :

    int64_t _objectId;
    double  _ra;
    double  _decl;
    int16_t  _varProb[lsst::afw::image::Filter::NUM_FILTERS];

    int64_t getId()  const { return _objectId; }
    double  getRa()  const { return _ra;       }
    double  getDec() const { return _decl;     }

    int16_t  getVarProb(lsst::afw::image::Filter const f) const { return _varProb[f]; }
};

LSST_AP_API bool operator==(SimpleObject const & o1, SimpleObject const & o2);

inline bool operator!=(SimpleObject const & o1, SimpleObject const & o2) {
    return !(o1 == o2);
}


/** @brief  A well named class. */
class LSST_AP_LOCAL PerFilterObjectData {

public :

    /** @brief  Nullable field ordinals, counting from 0. */
    enum NullableField {
        MAG = 0,
        MAG_ERR,
        ERR_A,
        ERR_B,
        ERR_THETA,
        TIMESCALE,
        AP_MAG,
        AP_MAG_ERR,
        ISO_AREA_IMAGE,
        MU_MAX,
        FLUX_RADIUS,
        PETRO_MAG,
        PETRO_MAG_ERR,
        AMPLITUDE,
        PERIOD,
        IXX,
        IXX_ERR,
        IYY,
        IYY_ERR,
        IXY,
        IXY_ERR,
        SCALEGRAM,
        NUM_OBS,
        FLAGS,
        VAR_PROB,
        NUM_NULLABLE_FIELDS
    };

    double  _mag;           // DOUBLE   NULL
    double  _magErr;        // DOUBLE   NULL
    double  _errA;          // DOUBLE   NULL
    double  _errB;          // DOUBLE   NULL
    double  _errTheta;      // DOUBLE   NULL
    double  _timescale;     // DOUBLE   NULL
    double  _apMag;         // DOUBLE   NULL
    double  _apMagErr;      // DOUBLE   NULL
    double  _isoAreaImage;  // DOUBLE   NULL
    double  _muMax;         // DOUBLE   NULL
    double  _fluxRadius;    // DOUBLE   NULL
    float   _petroMag;      // FLOAT(0) NULL
    float   _petroMagErr;   // FLOAT(0) NULL
    float   _amplitude;     // FLOAT(0) NULL
    float   _period;        // FLOAT(0) NULL
    float   _Ixx;           // FLOAT(0) NULL
    float   _IxxErr;        // FLOAT(0) NULL
    float   _Iyy;           // FLOAT(0) NULL
    float   _IyyErr;        // FLOAT(0) NULL
    float   _Ixy;           // FLOAT(0) NULL
    float   _IxyErr;        // FLOAT(0) NULL
    float   _scalegram[25]; // FLOAT(0) NULL
    int32_t _numObs;        // INTEGER  NULL
    int32_t _flags;         // INTEGER  NULL
    int16_t _varProb;       // SMALLINT NULL

    Bitset<uint8_t, NUM_NULLABLE_FIELDS> _nulls;

    bool    isNull    (NullableField const f) const            { return _nulls.test(f); }
    void    setNull   (NullableField const f)                  { _nulls.set(f);         }
    void    setNotNull(NullableField const f)                  { _nulls.reset(f);       }
    void    setNull   (NullableField const f, bool const null) { _nulls.set(f, null);   }
    void    setNull   ()                                       { _nulls.set();          }
    void    setNotNull()                                       { _nulls.reset();        }
};


/**
 * @brief  Contains data for Object records.
 *
 * This is a memory critical structure - for some AP designs, tens of millions of instances
 * can be held in memory.
 *
 * This structure is derived from the <b>DC2</b> version of the LSST MySQL schema, where the Object and
 * ObjectPhotoZ schema have been merged (storage for two sets of redshift parameters is provided).
 *
 * The C++ version consists of a per-filter structure and a common header, and columns are sorted
 * by type size. This minimizes the number of padding bytes the compiler must insert to meet type
 * alignment requirements.
 */
class LSST_AP_LOCAL Object {

public :

    /** @brief  Nullable field ordinals for header (non filter specific) fields. */
    enum NullableField {
        MU_RA = 0,
        MU_DECL,
        MU_RA_ERR,
        MU_DEC_ERR,
        UG_COLOR,
        GR_COLOR,
        RI_COLOR,
        IZ_COLOR,
        ZY_COLOR,
        CX,
        CX_ERR,
        CY,
        CY_ERR,
        CZ,
        CZ_ERR,
        EARLIEST_OBS_TIME,
        LATEST_OBS_TIME,
        PARALLAX,
        PARALLAX_ERR,
        REDSHIFT,
        REDSHIFT_ERR,
        PHOTO_Z1,
        PHOTO_Z1_ERR,
        PHOTO_Z2,
        PHOTO_Z2_ERR,
        PHOTO_Z1_OUTLIER,
        PHOTO_Z2_OUTLIER,
        PROC_HISTORY_ID,
        FLAG_4_STAGE1,
        FLAG_4_STAGE2,
        FLAG_4_STAGE3,
        IS_PROVISIONAL,
        PROBABILITY,
        NUM_NULLABLE_FIELDS
    };

    int64_t _objectId;       // BIGINT  NOT NULL
    double  _ra;             // DOUBLE  NOT NULL
    double  _decl;           // DOUBLE  NOT NULL
    double  _muRA;           // DOUBLE      NULL
    double  _muDecl;         // DOUBLE      NULL
    double  _muRAErr;        // DOUBLE      NULL
    double  _muDecErr;       // DOUBLE      NULL
    double  _ugColor;        // DOUBLE      NULL
    double  _grColor;        // DOUBLE      NULL
    double  _riColor;        // DOUBLE      NULL
    double  _izColor;        // DOUBLE      NULL
    double  _zyColor;        // DOUBLE      NULL
    double  _cx;             // DOUBLE      NULL
    double  _cxErr;          // DOUBLE      NULL
    double  _cy;             // DOUBLE      NULL
    double  _cyErr;          // DOUBLE      NULL
    double  _cz;             // DOUBLE      NULL
    double  _czErr;          // DOUBLE      NULL
    lsst::daf::base::DateTime _earliestObsTime; // DATETIME NULL
    lsst::daf::base::DateTime _latestObsTime;   // DATETIME NULL
    float   _parallax;       // FLOAT(0)    NULL
    float   _parallaxErr;    // FLOAT(0)    NULL
    float   _redshift;       // FLOAT(0)    NULL
    float   _redshiftErr;    // FLOAT(0)    NULL
    float   _photoZ1;        // FLOAT(0)    NULL
    float   _photoZ1Err;     // FLOAT(0)    NULL
    float   _photoZ2;        // FLOAT(0)    NULL
    float   _photoZ2Err;     // FLOAT(0)    NULL
    float   _photoZ1Outlier; // FLOAT(0)    NULL
    float   _photoZ2Outlier; // FLOAT(0)    NULL
    int32_t _procHistoryId;  // INTEGER     NULL
    int32_t _flag4stage1;    // INTEGER     NULL
    int32_t _flag4stage2;    // INTEGER     NULL
    int32_t _flag4stage3;    // INTEGER     NULL
    bool    _isProvisional;  // BOOL        NULL DEFAULT FALSE
    int8_t  _probability;    // TINYINT     NULL

    Bitset<uint8_t, NUM_NULLABLE_FIELDS> _nulls;

    PerFilterObjectData _filters[lsst::afw::image::Filter::NUM_FILTERS];

    // methods

    int64_t getId() const  { return _objectId; }
    double  getRa() const  { return _ra;       }
    double  getDec() const { return _decl;     }

    int8_t  getVarProb(lsst::afw::image::Filter const f) const { return _filters[f]._varProb; }

    PerFilterObjectData const & getFilter(lsst::afw::image::Filter const f) const {
        assert(f >= 0 && f < lsst::afw::image::Filter::NUM_FILTERS);
        return _filters[f];
    }

    PerFilterObjectData & getFilter(lsst::afw::image::Filter const f) {
        assert(f >= 0 && f < lsst::afw::image::Filter::NUM_FILTERS);
        return _filters[f];
    }

    bool    isNull    (NullableField const f) const            { return _nulls.test(f); }
    void    setNull   (NullableField const f)                  { _nulls.set(f);         }
    void    setNotNull(NullableField const f)                  { _nulls.reset(f);       }
    void    setNull   (NullableField const f, bool const null) { _nulls.set(f, null);   }
    void    setNull   ()                                       { _nulls.set();          }
    void    setNotNull()                                       { _nulls.reset();        }
};


}}  // end of namespace lsst::ap

#endif // LSST_AP_OBJECT_H
