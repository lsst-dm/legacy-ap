// -*- lsst-c++ -*-
%define ap_DOCSTRING
"
Access to association pipeline persistable result objects and implementation methods.
"
%enddef

%feature("autodoc", "1");
%module(package="lsst.ap", docstring=ap_DOCSTRING) interface

%{
#include "lsst/ap/Exceptions.h"
#include "lsst/ap/Results.h"
#include "lsst/ap/io/ResultFormatters.h"
#include "lsst/ap/Stages.h"
#include "lsst/ap/Utils.h"
%}

%init %{
%}

%pythoncode %{
import lsst.ap.exceptions
import lsst.fw.exceptions
%}

%include "lsst/mwi/p_lsstSwig.i"
%include "lsst/mwi/persistenceMacros.i"

%import "lsst/mwi/data/Citizen.h"
%import "lsst/mwi/persistence/Persistable.h"
%import "lsst/mwi/data/DataProperty.h"
%import "lsst/mwi/policy/Policy.h"
%import "lsst/mwi/persistence/LogicalLocation.h"
%import "lsst/mwi/persistence/Persistence.h"
%import "lsst/mwi/persistence/Storage.h"

%include <stdint.i>
%include <typemaps.i>
%include <std_pair.i>
%include <std_vector.i>

%import  "lsst/ap/Common.h"

%rename(IdVec)        lsst::ap::IdVector;
%rename(IdPairVec)    lsst::ap::IdPairVector;
%rename(MatchPairVec) lsst::ap::MatchPairVector;

%include "lsst/ap/Results.h"


%define %lsst_idpair(Type)
    %template(IdPair) std::pair<Type, Type>;

    %extend std::pair<Type, Type> {
        std::string toString() {
            std::ostringstream os;
            os << "IdPair (" << $self->first << ", " << $self->second << ")";
            return os.str();
        }
    };

    %pythoncode %{
    IdPair.__str__ = IdPair.toString
    %}
%enddef

#if defined(SWIGWORDSIZE64)
    %lsst_idpair(long);
#else
    %lsst_idpair(long long);
#endif

%extend lsst::ap::MatchPair {
    std::string toString() {
        std::ostringstream os;
        os << "MatchPair (" << $self->getFirst() << ", " << $self->getSecond() << ", " <<
              $self->getDistance() << ")";
        return os.str();
    }
};

%pythoncode %{
MatchPair.__str__ = MatchPair.toString
%}


%define %lsst_vector_traits(UnqualifiedType, Type)
    %fragment(#UnqualifiedType "Traits","header",fragment="StdSequenceTraits")
    %{
    namespace swig {
        template <>
        struct traits_asptr<Type>  {
            static int asptr(PyObject *obj, Type **vec) {
                return traits_asptr_stdseq<Type>::asptr(obj, vec);
            }
        };

        template <>
        struct traits_from<Type> {
            static PyObject *from(const Type& vec) {
                return traits_from_stdseq<Type>::from(vec);
            }
        };
    }
    %}
%enddef

%lsst_vector_traits(MatchPairVector, lsst::ap::MatchPairVector);
%lsst_vector_traits(IdPairVector,    lsst::ap::IdPairVector);
%lsst_vector_traits(IdVector,        lsst::ap::IdVector);


// An analogue to %std_vector_methods() with support for comparison operators
// and without support for std::vector::get_allocator()
%define %lsst_vector_methods(Type...)
    Type();
    Type(const Type&);
    Type(size_type size);
    Type(size_type size, const value_type& value);

    bool empty() const;
    size_type size() const;
    void clear();

    void swap(Type& v);

    #ifdef SWIG_EXPORT_ITERATOR_METHODS
    class iterator;
    class reverse_iterator;
    class const_iterator;
    class const_reverse_iterator;

    iterator begin();
    const_iterator begin() const;
    iterator end();
    const_iterator end() const;
    reverse_iterator rbegin();
    const_reverse_iterator rbegin() const;
    reverse_iterator rend();
    const_reverse_iterator rend() const;

    iterator erase(iterator pos);
    iterator erase(iterator first, iterator last);

    iterator insert(iterator pos, const value_type& x);
    void insert(iterator pos, size_type n, const value_type& x);
    #endif

    void pop_back();
    void push_back(const value_type& x);

    const value_type& front() const;
    const value_type& back() const;

    void assign(size_type n, const value_type& x);
    void resize(size_type new_size);
    void resize(size_type new_size, const value_type& x);

    void reserve(size_type n);
    size_type capacity() const;

    bool operator==(Type const & v);
    bool operator!=(Type const & v);
%enddef


// Applies SWIG std::vector machinery to result vectors. These classes contain a std::vector and look
// like a std::vector, but do not derive from one (since std::vector has a non-virtual destructor).
%define %lsst_persistable_vector(UnqualifiedType, Type, ValueType...)
    namespace lsst {
    namespace ap {

    class UnqualifiedType : public lsst::mwi::persistence::Persistable {
    public:
        typedef size_t size_type;
        typedef ptrdiff_t difference_type;
        typedef ValueType value_type;
        typedef value_type* pointer;
        typedef const value_type* const_pointer;
        typedef ValueType& reference;
        typedef const ValueType& const_reference;
        typedef std::allocator<ValueType > allocator_type;

        %traits_swigtype(ValueType);

        %fragment(SWIG_Traits_frag(Type), "header",
                  fragment=SWIG_Traits_frag(ValueType),
                  fragment=#UnqualifiedType "Traits") {
            namespace swig {
                template <>  struct traits<Type> {
                    typedef pointer_category category;
                    static const char* type_name() {
                        return #Type ;
                    }
                };
            }
        }

        %typemap_traits_ptr(SWIG_TYPECHECK_VECTOR, Type);

        #ifdef %swig_vector_methods
        // Add swig/language extra methods
        %swig_vector_methods(Type);
        #endif

        %lsst_vector_methods(Type);
    };

    }}
%enddef


%lsst_persistable_vector(MatchPairVector, lsst::ap::MatchPairVector, lsst::ap::MatchPair);
#if defined(SWIGWORDSIZE64)
    %lsst_persistable_vector(IdPairVector, lsst::ap::IdPairVector, std::pair<long, long>);
    %lsst_persistable_vector(IdVector, lsst::ap::IdVector, long);
#else
    %lsst_persistable_vector(IdPairVector, lsst::ap::IdPairVector, std::pair<long long, long long>);
    %lsst_persistable_vector(IdVector, lsst::ap::IdVector, long long);
#endif


// Make sure SWIG generates type information for certain types that are wrapped in other modules
%types(boost::shared_ptr<lsst::mwi::persistence::Persistable> *);

namespace lsst {
namespace fw {
    class DiaSourceVector;
    class MovingObjectPredictionVector;
}}

%types(boost::shared_ptr<lsst::fw::DiaSourceVector> *);
%types(boost::shared_ptr<lsst::fw::MovingObjectPredictionVector> *);


// Export instantiations of boost::shared_ptr for persistable data vectors
%lsst_persistable_shared_ptr(IdVecSharedPtr,        lsst::ap::IdVector);
%lsst_persistable_shared_ptr(IdPairVecSharedPtr,    lsst::ap::IdPairVector);
%lsst_persistable_shared_ptr(MatchPairVecSharedPtr, lsst::ap::MatchPairVector);

%pythoncode %{
def IdVecPtr(*args):
    """Creates an IdVec from the given arguments and returns an IdVecSharedPtr that owns it"""
    v = IdVec(*args)
    out = IdVecSharedPtr(v)
    return out

def IdPairVecPtr(*args):
    """Creates an IdPairVec from the given arguments and returns an IdPairVecSharedPtr that owns it"""
    v = IdPairVec(*args)
    out = IdPairVecSharedPtr(v)
    return out

def MatchPairVecPtr(*args):
    """Creates an MatchPairVec from the given arguments and returns an MatchPairVecSharedPtr that owns it"""
    v = MatchPairVec(*args)
    out = MatchPairVecSharedPtr(v)
    return out
%}


namespace lsst {
namespace ap {
    std::string const getTableName(
        boost::shared_ptr<lsst::mwi::policy::Policy> const &,
        boost::shared_ptr<lsst::mwi::data::DataProperty> const &
    );
    std::string const getTableTemplateName(
        boost::shared_ptr<lsst::mwi::policy::Policy> const &,
        boost::shared_ptr<lsst::mwi::data::DataProperty> const &
    );
}}


%include "lsst/ap/Stages.h"

