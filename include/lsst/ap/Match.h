// -*- lsst-c++ -*-

/**
 * @file
 * @brief   Spatial crossmatch algorithms.
 *
 * The spatial crossmatch routines currently implemented take match list (or pair) processors
 * as parameters. In this model, a processor could for example be some part of the LSST alert
 * generation pipeline. More likely, it will be something that sets various flags on entities
 * involved in the match and/or records matches in some form (e.g. as an id-id pair) for later
 * use.
 *
 * Note that the routines currently don't do anything special (like match-processing) for cases
 * where no matches to a given entity are found. The routines could be extended with a
 * no-match-processor parameter, or the existing processors could simply be called with empty
 * match lists. Finally, the user could write a match processor that flags all entities involved
 * in a match, and then post-process over unflagged entities.
 *
 * Note that a single match pair/list processor instance is used by all threads involved
 * in a parallel match. Ensuring that this doesn't cause problems is the responsibility of the
 * match pair/list processor author.
 *
 * @ingroup associate
 */

#ifndef LSST_AP_MATCH_H
#define LSST_AP_MATCH_H

#if LSST_AP_HAVE_OPEN_MP
#   include <omp.h>
#endif

#include <stdexcept>
#include <vector>

#include "Common.h"
#include "EllipseTypes.h"
#include "SpatialUtil.h"
#include "ZoneTypes.h"


namespace lsst { namespace ap {


/** @brief  A default "let everything through" filter implementation. */
template <typename T>
struct PassthroughFilter {
    bool operator() (T const &) { return true; }
};


/** @brief  Contains a pointer to a match and an associated distance. */
template <typename T>
struct MatchWithDistance {

    MatchWithDistance() {}

    MatchWithDistance(T * const match, double const d2) :
        _match(match),
        _distance(2.0*std::asin(0.5*std::sqrt(d2)))
    {}

    T & operator*()  const { return *_match; }
    T * operator->() const { return  _match; }

    T *    _match;
    double _distance;
};


/** @brief Contains a pointer to a match. */
template <typename T>
struct MatchWithoutDistance {

    MatchWithoutDistance() {}

    MatchWithoutDistance(T * const match, double const) : _match(match) {}

    T & operator*()  const { return *_match; }
    T * operator->() const { return  _match; }

    T * _match;
};


/** @brief  A default "do nothing" match list processing implementation. */
template <typename F, typename M>
struct EmptyMatchListProcessor {

    typedef M Match;
    typedef typename std::vector<M>::const_iterator MatchIterator;

    inline void operator() (F & entity, MatchIterator begin, MatchIterator end) {}
};


/** @brief  A default "do nothing" match pair processing implementation. */
template <typename F, typename S>
struct EmptyMatchPairProcessor {
    inline void operator() (F & first, S & second) {}
};


/**
 * Spatial cross-match routine -- finds match pairs in a first and a second set of
 * entities (both subject to filtering), where both sets consist of points. An entity from
 * the second set is deemed a match to an entity from the first set if the two are within
 * the given angle of eachother. All matches for a given entity in the first set are found
 * at once, and sent off to a match list processor for further inspection.
 *
 * This routine is optimized for the case where few matches are expected for any given entity.
 *
 * @pre @code radius >= 0 @endcode
 *
 * @param[in] first                 A first set of entities.
 * @param[in] second                A second set of entities.
 * @param[in] radius                The match-radius (in degrees).
 * @param[in] firstFilter           A filter on the first set of entities.
 * @param[in] secondFilter          A filter on the second set of entities.
 * @param[in] matchListProcessor    A processor for match lists.
 * @return                          The number of match pairs found.
 */
template <
    typename FirstEntryT,
    typename SecondEntryT,
    typename FirstFilterT,        // = PassthroughFilter<FirstEntryT>,
    typename SecondFilterT,       // = PassthroughFilter<SecondEntryT>,
    typename MatchListProcessorT  // = EmptyMatchListProcessor<FirstEntryT, MatchWithoutDistance<SecondEntryT> >
>
std::size_t distanceMatch(
    ZoneIndex<FirstEntryT>  & first,
    ZoneIndex<SecondEntryT> & second,
    double const              radius,
    FirstFilterT            & firstFilter,
    SecondFilterT           & secondFilter,
    MatchListProcessorT     & matchListProcessor
) {
    typedef typename ZoneIndex<FirstEntryT>::Zone  FirstZone;
    typedef typename ZoneIndex<SecondEntryT>::Zone SecondZone;
    typedef typename MatchListProcessorT::Match    Match;

    if (radius < 0) {
        throw LSST_EXCEPT(lsst::pex::exceptions::RangeErrorException,
                          "match radius must be greater than or equal to zero degrees");
    }

    double const shr     = std::sin(radians(radius*0.5));
    double const d2Limit = 4.0*shr*shr;
    int const minZone = first.getMinZone();
    int const maxZone = first.getMaxZone();
    boost::int32_t const deltaDec = deltaDecToScaledInteger(radius);

    std::size_t numMatchPairs = 0;

    second.computeMatchParams(radius);

    // loop over first set of zones in parallel
#if LSST_AP_HAVE_OPEN_MP
#   pragma omp parallel default(shared)
#endif
    {
        // allocate per-thread data structures
        SecondZone * zones[2048];
        int limits[2048];
        std::vector<Match> matches;
        matches.reserve(32);

        // loop over the first set of zones in parallel, assigning batches of
        // adjacent zones to threads (for cache coherency)
#if LSST_AP_HAVE_OPEN_MP
#   pragma omp for \
               reduction(+:numMatchPairs) \
               schedule(static,8)
#endif
        for (int fzi = minZone; fzi <= maxZone; ++fzi) {

            FirstZone   * const __restrict fz = first.getZone(fzi);
            FirstEntryT * const __restrict fze = fz->_entries;
            int const nfze = fz->_size;
            if (nfze <= 0) {
                continue; // no entries in first zone
            }

            // populate secondary zone array with potentially matching zones
            int nsz = 0;
            {
                double d = first.getDecomposition().getZoneDecMin(fz->_zone) - radius;
                int const minz = second.getDecomposition().decToZone(d <= -90.0 ? -90.0 : d);
                d = first.getDecomposition().getZoneDecMax(fz->_zone) + radius;
                int const maxz = second.getDecomposition().decToZone(d >= 90.0 ? 90.0 : d);
                // A search circle should never cover more than 2048 zones
                assert(maxz - minz + 1 <= 2048 && "match radius too large");

                SecondZone *       __restrict sz    = second.firstZone(minz, maxz);
                SecondZone * const __restrict szend = second.endZone(minz, maxz);

                for ( ; sz < szend; ++sz) {
                    int const nsze = sz->_size;
                    if (nsze > 0) {
                        zones[nsz] = sz;
                        if ((nsze >> 4) > nfze) {
                            // second set much larger than first, use binary search in inner loop
                            limits[nsz] = -1;
                        } else {
                            // use linear walk in inner loop, find starting point now
                            limits[nsz] = sz->findGte(fze[0]._ra - sz->_deltaRa);
                        }
                        ++nsz;
                    }
                }
            }

            if (nsz == 0) {
                // no entries in any potentially matching zones
                continue;
            }

            // loop over entries in first zone
            for (int fe = 0; fe < nfze; ++fe) {

                if (!firstFilter(fze[fe])) {
                    continue; // entry was filtered out
                }

                matches.clear();

                boost::uint32_t const ra  = fze[fe]._ra;
                boost::int32_t  const dec = fze[fe]._dec;
                double const fx = fze[fe]._x;
                double const fy = fze[fe]._y;
                double const fz = fze[fe]._z;

                // loop over all potentially matching zones.
                for (int szi = 0; szi < nsz; ++szi) {

                    SecondZone   * const __restrict sz  = zones[szi];
                    SecondEntryT * const __restrict sze = sz->_entries;

                    int const seWrap = sz->_size;
                    boost::uint32_t const deltaRa = sz->_deltaRa;
                    boost::uint32_t const deltaRaWrap = -deltaRa;

                    int start = limits[szi];
                    int se;
                    boost::uint32_t dra;

                    if (start >= 0) {

                        // use linear walk from last starting point
                        // to get to first point in ra range
                        se  = start;
                        dra = ra - sze[se]._ra;
                        if (dra > deltaRa && dra < deltaRaWrap) {

                            bool cont = false;
                            do {
                                ++se;
                                if (se == seWrap) {
                                    se = 0; // ra wrap around
                                }
                                if (se == start) {
                                    cont = true;
                                    break; // avoid infinite loops
                                }
                                dra = ra - sze[se]._ra;
                            } while (dra > deltaRa && dra < deltaRaWrap);

                            if (cont) {
                                continue; // no starting point found
                            }
                            // found starting point -- remember it
                            limits[szi] = se;
                            start = se;
                        }

                    } else {

                        // use binary search to find starting point
                        start = sz->findGte(ra - deltaRa);
                        se    = start;
                        dra   = ra - sze[se]._ra;
                        if (dra > deltaRa && dra < deltaRaWrap) {
                            continue; // no starting point found
                        }

                    }

                    // At this point, entry se is within ra range of fe --
                    // loop over all zone entries within ra range
                    do {

                        // test whether entry se is within dec range
                        boost::int32_t ddec = dec - sze[se]._dec;
                        // C standard says: result of a right shifting a negative integer
                        // is implementation defined. If a signed right shift is available,
                        // use it to perform a branchless abs().
#if LSST_AP_HAVE_SIGNED_RSHIFT
                        boost::int32_t sgn  = ddec >> 31;
                        ddec = (ddec ^ sgn) - sgn; // abs(dec - sze[se]._dec)
#else
                        if (ddec < 0) {
                            ddec = -ddec;
                        }
#endif
                        if (ddec <= deltaDec) {
                            // yes -- perform detailed distance test
                            double xd = (fx - sze[se]._x);
                            double yd = (fy - sze[se]._y);
                            double zd = (fz - sze[se]._z);
                            double d2 = xd*xd + yd*yd + zd*zd;
                            if (d2 < d2Limit) {
                                // Note: this isn't necessarily the best place for the second filter test...
                                if (secondFilter(sze[se])) {
                                    // found a match, record it
                                    matches.push_back(Match(&sze[se], d2));
                                }
                            }
                        }
                        ++se;
                        if (se == seWrap) {
                            se = 0; // ra wrap around
                        }
                        if (se == start) {
                            break; // avoid infinite loops
                        }
                        dra = ra - sze[se]._ra;

                    } while (dra <= deltaRa || dra >= deltaRaWrap);

                } // end of loop over potentially matching zones

                // All matches (if any) for fe are found
                std::size_t nm = matches.size();
                if (nm > 0) {
                    // pass them on to the match processor
                    numMatchPairs = numMatchPairs + nm;
                    matchListProcessor(fze[fe], matches.begin(), matches.end());
                }

            } // end of loop over entries in first zone

        } // end of omp for
    } // end of omp parallel

    return numMatchPairs;
}


/**
 * Spatial cross-match routine -- finds match pairs in a first and a second set of entities
 * (both subject to filtering), where the first set consists of ellipses and the second of
 * points. An entity in the second set is deemed a match for an entity in the first set if it
 * is within the ellipse defined by the first entity. Note -- this particular routine does
 * not package up all matches for a given ellipse before sending them off to a match list
 * processor. Instead, it reports match pairs to a match pair processor in no particular order.
 * The routine is also currently single threaded. For a parallel routine that uses a
 * match list processor, see ellipseGroupedMatch(). The advantage of this routine is that it
 * should have much better memory access patterns (smaller cache footprint, more sequential
 * access) in the presence of ellipses with wildly varying size.
 *
 * @param[in] first                 A first set of entities.
 * @param[in] second                A second set of entities.
 * @param[in] firstFilter           A filter on the first set of entities.
 * @param[in] secondFilter          A filter on the second set of entities.
 * @param[in] matchPairProcessor    A processor for match pairs.
 * @return                          The number of match pairs found.
 */
template <
    typename FirstEntryT,
    typename SecondEntryT,
    typename FirstFilterT,        // = PassthroughFilter<FirstEntryT>,
    typename SecondFilterT,       // = PassthroughFilter<SecondEntryT>,
    typename MatchPairProcessorT  // = EmptyMatchPairProcessor<FirstEntryT, SecondEntryT>
>
std::size_t ellipseMatch(
    EllipseList<FirstEntryT> & first,
    ZoneIndex<SecondEntryT>  & second,
    FirstFilterT             & firstFilter,
    SecondFilterT            & secondFilter,
    MatchPairProcessorT      & matchPairProcessor
) {

    typedef Ellipse<FirstEntryT> Ellipse;
    typedef typename ZoneIndex<SecondEntryT>::Zone SecondZone;

    int const minZone = second.getMinZone();
    int const maxZone = second.getMaxZone();
    std::size_t numMatchPairs = 0;

    first.prepareForMatch(second.getDecomposition());

    Ellipse * activeHead = 0;
    Ellipse * activeTail = 0;
    Ellipse * searchHead = &(*first.begin());
    Ellipse * const end = &(*first.end());

    // initialize the linked list of active ellipses (those that intersect minZone)
    while (searchHead < end && searchHead->_minZone <= minZone) {
        if (searchHead->_maxZone >= minZone) {
            if (firstFilter(*searchHead)) {
                if (activeTail == 0) {
                    activeHead = searchHead;
                } else {
                    activeTail->_next = searchHead;
                }
                activeTail = searchHead;
            }
        }
        ++searchHead;
    }

    // loop over the set of zones in the second entity set
    int szi = minZone;
    while (true) {

        SecondZone   * const __restrict sz  = second.getZone(szi);
        SecondEntryT * const __restrict sze = sz->_entries;
        int const seWrap = sz->_size;

        ++szi;

        // Traverse the active ellipse list
        Ellipse * active = activeHead;
        Ellipse * prev   = 0;

        while (active != 0) {

#if defined(__GNUC__)
            __builtin_prefetch(active->_next, 0, 3);
#endif
            if (seWrap > 0) {
                // find entries within ra range of the active ellipse
                boost::uint32_t const ra          = active->_ra;
                boost::uint32_t const deltaRa     = active->_deltaRa;
                boost::uint32_t const deltaRaWrap = -deltaRa;
                int const start = sz->findGte(ra - deltaRa);

                int se  = start;
                boost::uint32_t dra = ra - sze[se]._ra;
                while (dra <= deltaRa || dra >= deltaRaWrap) {

                    // check whether entry is within dec range
                    boost::int32_t const dec = sze[se]._dec;
                    if (dec >= active->_minDec && dec <= active->_maxDec) {
                        // perform detailed in ellipse test
                        if (active->contains(sze[se]._x, sze[se]._y, sze[se]._z)) {
                            if (secondFilter(sze[se])) {
                                // process match pair
                                matchPairProcessor(*active, sze[se]);
                                ++numMatchPairs;
                            }
                        }
                    }
                    ++se;
                    if (se == seWrap) {
                        se = 0; // ra wrap around
                    }
                    if (se == start) {
                        break; // avoid infinite loops
                    }
                    dra = ra - sze[se]._ra;
                }
            }

            if (active->_maxZone < szi) {
                // ellipse no longer active in following zone -- remove it from the active list
                active = active->_next;
                if (prev == 0) {
                    activeHead = active;
                } else {
                    prev->_next = active;
                }
            } else {
                // keep ellipse in the active list
                prev   = active;
                active = active->_next;
            }

        } // end of loop over active ellipses

        if (szi > maxZone) {
            break;
        }

        // Add ellipses with minimum zone equal to szi (the zone that will be searched
        // in the next loop iteration) to the active list
        while (searchHead < end && searchHead->_minZone <= szi) {
            if (firstFilter(*searchHead)) {
                if (prev == 0) {
                    activeHead = searchHead;
                } else {
                    prev->_next = searchHead;
                }
                prev = searchHead;
            }
            ++searchHead;
        }

    } // end of loop over second set of zones

    return numMatchPairs;
}


/**
 * Spatial cross-match routine -- finds match pairs in a first and a second set of entities
 * (both subject to filtering), where the first set consists of ellipses and the second of
 * points. An entity in the second set is deemed a match for an entity in the first set if it
 * is within the ellipse defined by the first entity. All matches for a given ellipse are found
 * at once and sent off to a match list processor for further inspection.
 *
 * @param[in] first                 A first set of entities.
 * @param[in] second                A second set of entities.
 * @param[in] firstFilter           A filter on the first set of entities.
 * @param[in] secondFilter          A filter on the second set of entities.
 * @param[in] matchListProcessor    A processor for match lists.
 * @return                          The number of match pairs found.
 */
template <
    typename FirstEntryT,
    typename SecondEntryT,
    typename FirstFilterT,       // = PassthroughFilter<FirstEntryT>,
    typename SecondFilterT,      // = PassthroughFilter<SecondEntryT>,
    typename MatchListProcessorT // = EmptyMatchListProcessor<FirstEntryT, SecondEntryT>
>
std::size_t ellipseGroupedMatch(
    EllipseList<FirstEntryT> & first,
    ZoneIndex<SecondEntryT>  & second,
    FirstFilterT             & firstFilter,
    SecondFilterT            & secondFilter,
    MatchListProcessorT      & matchListProcessor
) {
    typedef typename EllipseList<FirstEntryT>::Ellipse Ellipse;
    typedef typename ZoneIndex<SecondEntryT>::Zone     SecondZone;
    typedef typename MatchListProcessorT::Match        Match;

    std::size_t const numEllipses   = first.size();
    std::size_t       numMatchPairs = 0;

    first.prepareForMatch(second.getDecomposition());

#if LSST_AP_HAVE_OPEN_MP
#   pragma omp parallel default(shared)
#endif
    {
        // allocate per-thread match list
        std::vector<Match> matches;
        matches.reserve(2048);

        // loop over the list of ellipses in parallel
#if LSST_AP_HAVE_OPEN_MP
#   pragma omp for \
               reduction(+:numMatchPairs) \
               schedule(static,128)
#endif
        for (std::size_t i = 0; i < numEllipses; ++i) {

            if (!firstFilter(first[i])) {
                continue; // ellipse was filtered out
            }

            Ellipse    * const __restrict ell   = &first[i];
            SecondZone *       __restrict sz    = second.firstZone(ell->_minZone, ell->_maxZone);
            SecondZone * const __restrict szend = second.endZone(ell->_minZone, ell->_maxZone);

            boost::uint32_t const ra          = ell->_ra;
            boost::uint32_t const deltaRa     = ell->_deltaRa;
            boost::uint32_t const deltaRaWrap = -deltaRa;

            matches.clear();

            // loop over zones covered by the ellipse
            for ( ; sz < szend; ++sz) {

                int const seWrap = sz->_size;
                if (seWrap == 0) {
                    continue;
                }

                // find starting point within ra range of the ellipse
                SecondEntryT * const __restrict sze = sz->_entries;
                int const start = sz->findGte(ra - deltaRa);
                int se = start;
                boost::uint32_t dra = ra - sze[se]._ra;

                while (dra <= deltaRa || dra >= deltaRaWrap) {

                    // check whether entry is within dec range
                    boost::int32_t const dec = sze[se]._dec;
                    if (dec >= ell->_minDec && dec <= ell->_maxDec) {
                        // perform detailed in ellipse test
                        if (ell->contains(sze[se]._x, sze[se]._y, sze[se]._z)) {
                            if (secondFilter(sze[se])) {
                                // record match
                                matches.push_back(Match(&sze[se]));
                            }
                        }
                    }
                    ++se;
                    if (se == seWrap) {
                        se = 0; // ra wrap around
                    }
                    if (se == start) {
                        break; // avoid infinite loops
                    }
                    dra = ra - sze[se]._ra;
                }
            }

            // All matches (if any) for ell are found
            std::size_t nm = matches.size();
            if (nm > 0) {
                // pass them on to the match processor
                numMatchPairs = numMatchPairs + nm;
                matchListProcessor(*ell, matches.begin(), matches.end());
            }

        } // end of omp for
    } // end of omp parallel

    return numMatchPairs;
}


}}  // end of namespace lsst::ap

#endif // LSST_AP_MATCH_H
