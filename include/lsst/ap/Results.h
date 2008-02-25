// -*- lsst-c++ -*-

/**
 * @file
 * @brief   Persistable vectors for association pipeline results.
 *
 * @ingroup associate
 */

#ifndef LSST_AP_RESULTS_H
#define LSST_AP_RESULTS_H

#include <utility>
#include <vector>

#include <boost/shared_ptr.hpp>

#include <lsst/mwi/data/Citizen.h>
#include <lsst/mwi/persistence/Persistable.h>

#include "Common.h"


/// @cond
namespace boost {
namespace serialization {
    class access;
}}
/// @endcond


namespace lsst {
namespace ap {

// forward declarations for formatters
namespace io {
    class MatchPairVectorFormatter;
    class IdPairVectorFormatter;
    class IdVectorFormatter;
}


/** @brief  Holds a pair of ids and the distance between the corresponding positions on the unit sphere. */
class LSST_AP_API MatchPair {

public :

    MatchPair() : _first(-1), _second(-1), _distance(0.0) {}

    MatchPair(int64_t const first, int64_t const second, double const distance) :
        _first(first), _second(second), _distance(distance) {}

    int64_t getFirst()    const { return _first;    }
    int64_t getSecond()   const { return _second;   }
    double  getDistance() const { return _distance; }

    void setFirst   (int64_t const first)    { _first    = first;    }
    void setSecond  (int64_t const second)   { _second   = second;   }
    void setDistance(double  const distance) { _distance = distance; }

private :

    int64_t _first;
    int64_t _second;
    double  _distance;

    template <typename Archive> void serialize(Archive & ar, unsigned int const version) {
        ar & _first;
        ar & _second;
        ar & _distance;
    }

    friend class boost::serialization::access;
    friend class io::MatchPairVectorFormatter;
};

inline bool operator==(MatchPair const & mp1, MatchPair const & mp2) {
    return mp1.getFirst() == mp2.getFirst() && mp1.getSecond() == mp2.getSecond();
}

inline bool operator!=(MatchPair const & d1, MatchPair const & d2) {
    return !(d1 == d2);
}


// Classes that require special handling in the SWIG interface file follow
#ifndef SWIG

/** @brief  Holds a pair of ids. */
typedef std::pair<int64_t, int64_t> IdPair;


/** @brief  A persistable container of MatchPair instances, implemented using std::vector. */
class LSST_AP_API MatchPairVector :
    public lsst::mwi::persistence::Persistable,
    public lsst::mwi::data::Citizen
{
public :

    typedef boost::shared_ptr<MatchPairVector> Ptr;
    typedef std::vector<MatchPair>             Vector;

    typedef Vector::allocator_type         allocator_type;
    typedef Vector::iterator               iterator;
    typedef Vector::const_iterator         const_iterator;
    typedef Vector::reverse_iterator       reverse_iterator;
    typedef Vector::const_reverse_iterator const_reverse_iterator;
    typedef Vector::size_type              size_type;
    typedef Vector::difference_type        difference_type;
    typedef Vector::reference              reference;
    typedef Vector::const_reference        const_reference;
    typedef Vector::value_type             value_type;

    MatchPairVector();
    explicit MatchPairVector(size_type sz);
    MatchPairVector(size_type sz, value_type const & val);

    template <typename InputIteratorT>
    MatchPairVector(InputIteratorT beg, InputIteratorT end) :
        lsst::mwi::data::Citizen(typeid(*this)),
        _vec(beg, end)
    {}

    virtual ~MatchPairVector();

    MatchPairVector(MatchPairVector const & vec);
    explicit MatchPairVector(Vector const & vec);
    MatchPairVector & operator=(MatchPairVector const & vec);
    MatchPairVector & operator=(Vector const & vec);

    void swap(MatchPairVector & v) { using std::swap; swap(_vec, v._vec); }
    void swap(Vector & v)          { using std::swap; swap(_vec, v);      }

    size_type size()     const { return _vec.size();     }
    size_type max_size() const { return _vec.max_size(); }
    bool      empty()    const { return _vec.empty();    }
    size_type capacity() const { return _vec.capacity(); }

    void reserve(size_type const n) { _vec.reserve(n); }

    template <typename InputIteratorT>
    void assign(InputIteratorT beg, InputIteratorT end)    { _vec.assign(beg, end); }
    void assign(size_type const n, value_type const & val) { _vec.assign(n, val);   }

    reference       at        (size_type const i)       { return _vec.at(i); }
    const_reference at        (size_type const i) const { return _vec.at(i); }
    reference       operator[](size_type const i)       { return _vec[i];    }
    const_reference operator[](size_type const i) const { return _vec[i];    }

    reference       front()       { return _vec.front(); }
    const_reference front() const { return _vec.front(); }
    reference       back ()       { return _vec.back();  }
    const_reference back () const { return _vec.back();  }

    iterator               begin ()       { return _vec.begin();  }
    const_iterator         begin () const { return _vec.begin();  }
    reverse_iterator       rbegin()       { return _vec.rbegin(); }
    const_reverse_iterator rbegin() const { return _vec.rbegin(); }
    iterator               end   ()       { return _vec.end();    }
    const_iterator         end   () const { return _vec.end();    }
    reverse_iterator       rend  ()       { return _vec.rend();   }
    const_reverse_iterator rend  () const { return _vec.rend();   }

    void push_back (value_type const & value) { _vec.push_back(value);  }

    void pop_back () { _vec.pop_back();  }
    void clear()     { _vec.clear();     }

    template <typename InputIteratorT>
    void     insert(iterator pos, InputIteratorT beg, InputIteratorT end) { _vec.insert(pos, beg, end);      }
    iterator insert(iterator pos, value_type const & val)                 { return _vec.insert(pos, val);    }
    void     insert(iterator pos, size_type n, value_type const & val)    { return _vec.insert(pos, n, val); }

    iterator erase(iterator pos)               { return _vec.erase(pos);      }
    iterator erase(iterator beg, iterator end) { return _vec.erase(beg, end); }

    void resize(size_type n)                   { _vec.resize(n);        }
    void resize(size_type n, value_type value) { _vec.resize(n, value); }

    bool operator==(MatchPairVector const & v) { return _vec == v._vec; }
    bool operator!=(MatchPairVector const & v) { return _vec != v._vec; }

private :

    LSST_PERSIST_FORMATTER(io::MatchPairVectorFormatter);

    Vector _vec;
};


/** @brief  A persistable container of IdPair identifiers, implemented using std::vector. */
class LSST_AP_API IdPairVector :
    public lsst::mwi::persistence::Persistable,
    public lsst::mwi::data::Citizen
{
public :

    typedef boost::shared_ptr<IdPairVector> Ptr;
    typedef std::vector<IdPair>             Vector;

    typedef Vector::allocator_type         allocator_type;
    typedef Vector::iterator               iterator;
    typedef Vector::const_iterator         const_iterator;
    typedef Vector::reverse_iterator       reverse_iterator;
    typedef Vector::const_reverse_iterator const_reverse_iterator;
    typedef Vector::size_type              size_type;
    typedef Vector::difference_type        difference_type;
    typedef Vector::reference              reference;
    typedef Vector::const_reference        const_reference;
    typedef Vector::value_type             value_type;

    IdPairVector();
    explicit IdPairVector(size_type sz);
    IdPairVector(size_type sz, value_type const & val);

    template <typename InputIteratorT>
    IdPairVector(InputIteratorT beg, InputIteratorT end) :
        lsst::mwi::data::Citizen(typeid(*this)),
        _vec(beg, end)
    {}

    virtual ~IdPairVector();

    IdPairVector(IdPairVector const & vec);
    explicit IdPairVector(Vector const & vec);
    IdPairVector & operator=(IdPairVector const & vec);
    IdPairVector & operator=(Vector const & vec);

    void swap(IdPairVector & v) { using std::swap; swap(_vec, v._vec); }
    void swap(Vector & v)       { using std::swap; swap(_vec, v);      }

    size_type size()     const { return _vec.size();     }
    size_type max_size() const { return _vec.max_size(); }
    bool      empty()    const { return _vec.empty();    }
    size_type capacity() const { return _vec.capacity(); }

    void reserve(size_type const n) { _vec.reserve(n); }

    template <typename InputIteratorT>
    void assign(InputIteratorT beg, InputIteratorT end)    { _vec.assign(beg, end); }
    void assign(size_type const n, value_type const & val) { _vec.assign(n, val);   }

    reference       at        (size_type const i)       { return _vec.at(i); }
    const_reference at        (size_type const i) const { return _vec.at(i); }
    reference       operator[](size_type const i)       { return _vec[i];    }
    const_reference operator[](size_type const i) const { return _vec[i];    }

    reference       front()       { return _vec.front(); }
    const_reference front() const { return _vec.front(); }
    reference       back ()       { return _vec.back();  }
    const_reference back () const { return _vec.back();  }

    iterator               begin ()       { return _vec.begin();  }
    const_iterator         begin () const { return _vec.begin();  }
    reverse_iterator       rbegin()       { return _vec.rbegin(); }
    const_reverse_iterator rbegin() const { return _vec.rbegin(); }
    iterator               end   ()       { return _vec.end();    }
    const_iterator         end   () const { return _vec.end();    }
    reverse_iterator       rend  ()       { return _vec.rend();   }
    const_reverse_iterator rend  () const { return _vec.rend();   }

    void push_back (value_type const & value) { _vec.push_back(value);  }

    void pop_back () { _vec.pop_back();  }
    void clear()     { _vec.clear();     }

    template <typename InputIteratorT>
    void     insert(iterator pos, InputIteratorT beg, InputIteratorT end) { _vec.insert(pos, beg, end);      }
    iterator insert(iterator pos, value_type const & val)                 { return _vec.insert(pos, val);    }
    void     insert(iterator pos, size_type n, value_type const & val)    { return _vec.insert(pos, n, val); }

    iterator erase(iterator pos)               { return _vec.erase(pos);      }
    iterator erase(iterator beg, iterator end) { return _vec.erase(beg, end); }

    void resize(size_type n)                   { _vec.resize(n);        }
    void resize(size_type n, value_type value) { _vec.resize(n, value); }

    bool operator==(IdPairVector const & v) { return _vec == v._vec; }
    bool operator!=(IdPairVector const & v) { return _vec != v._vec; }

private :

    LSST_PERSIST_FORMATTER(io::IdPairVectorFormatter);

    Vector _vec;
};


/** @brief  A persistable container of integer (int64_t) identifiers, implemented using std::vector. */
class LSST_AP_API IdVector :
    public lsst::mwi::persistence::Persistable,
    public lsst::mwi::data::Citizen
{
public :

    typedef boost::shared_ptr<IdVector>    Ptr;
    typedef std::vector<int64_t>           Vector;

    typedef Vector::allocator_type         allocator_type;
    typedef Vector::iterator               iterator;
    typedef Vector::const_iterator         const_iterator;
    typedef Vector::reverse_iterator       reverse_iterator;
    typedef Vector::const_reverse_iterator const_reverse_iterator;
    typedef Vector::size_type              size_type;
    typedef Vector::difference_type        difference_type;
    typedef Vector::reference              reference;
    typedef Vector::const_reference        const_reference;
    typedef Vector::value_type             value_type;

    IdVector();
    explicit IdVector(size_type sz);
    IdVector(size_type sz, value_type const & val);

    template <typename InputIteratorT>
    IdVector(InputIteratorT beg, InputIteratorT end) :
        lsst::mwi::data::Citizen(typeid(*this)),
        _vec(beg, end)
    {}

    virtual ~IdVector();

    IdVector(IdVector const & vec);
    explicit IdVector(Vector const & vec);
    IdVector & operator=(IdVector const & vec);
    IdVector & operator=(Vector const & vec);

    void swap(IdVector & v) { using std::swap; swap(_vec, v._vec); }
    void swap(Vector & v)   { using std::swap; swap(_vec, v);      }

    size_type size()     const { return _vec.size();     }
    size_type max_size() const { return _vec.max_size(); }
    bool      empty()    const { return _vec.empty();    }
    size_type capacity() const { return _vec.capacity(); }

    void reserve(size_type const n) { _vec.reserve(n); }

    template <typename InputIteratorT>
    void assign(InputIteratorT beg, InputIteratorT end)    { _vec.assign(beg, end); }
    void assign(size_type const n, value_type const & val) { _vec.assign(n, val);   }

    reference       at        (size_type const i)       { return _vec.at(i); }
    const_reference at        (size_type const i) const { return _vec.at(i); }
    reference       operator[](size_type const i)       { return _vec[i];    }
    const_reference operator[](size_type const i) const { return _vec[i];    }

    reference       front()       { return _vec.front(); }
    const_reference front() const { return _vec.front(); }
    reference       back ()       { return _vec.back();  }
    const_reference back () const { return _vec.back();  }

    iterator               begin ()       { return _vec.begin();  }
    const_iterator         begin () const { return _vec.begin();  }
    reverse_iterator       rbegin()       { return _vec.rbegin(); }
    const_reverse_iterator rbegin() const { return _vec.rbegin(); }
    iterator               end   ()       { return _vec.end();    }
    const_iterator         end   () const { return _vec.end();    }
    reverse_iterator       rend  ()       { return _vec.rend();   }
    const_reverse_iterator rend  () const { return _vec.rend();   }

    void push_back (value_type const & value) { _vec.push_back(value);  }

    void pop_back () { _vec.pop_back();  }
    void clear()     { _vec.clear();     }

    template <typename InputIteratorT>
    void     insert(iterator pos, InputIteratorT beg, InputIteratorT end) { _vec.insert(pos, beg, end);      }
    iterator insert(iterator pos, value_type const & val)                 { return _vec.insert(pos, val);    }
    void     insert(iterator pos, size_type n, value_type const & val)    { return _vec.insert(pos, n, val); }

    iterator erase(iterator pos)               { return _vec.erase(pos);      }
    iterator erase(iterator beg, iterator end) { return _vec.erase(beg, end); }

    void resize(size_type n)                   { _vec.resize(n);        }
    void resize(size_type n, value_type value) { _vec.resize(n, value); }

    bool operator==(IdVector const & v) { return _vec == v._vec; }
    bool operator!=(IdVector const & v) { return _vec != v._vec; }

private :

    LSST_PERSIST_FORMATTER(io::IdVectorFormatter);

    Vector _vec;
};

#endif // SWIG


}}  // end of namespace lsst::ap

#endif // LSST_AP_RESULTS_H

