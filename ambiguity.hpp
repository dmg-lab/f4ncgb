#ifndef AMBIGUITY_H
#define AMBIGUITY_H

#include <algorithm>
#include <concepts>
#include <boost/functional/hash.hpp>

namespace kommunopp {

template<typename I>
size_t compute_hash(const std::array<I,7>& values) {
    size_t seed = 0;
    for (I value : values)
         boost::hash_combine(seed, value);
    return seed;
}

// container to store when two monomials overlaps
template<typename I>
struct ambiguity {
    std::array<I, 7> values;
    size_t hash_;
    
    ambiguity(I d, I i, I j, I ai, I ci, I aj, I cj) {
        values[0] = d;
        values[1] = i; values[2] = j;
        values[3] = ai; values[4] = ci;
        values[5] = aj; values[6] = cj;
        hash_ = compute_hash(values);
    }
    
    inline I degree() const {return values[0];}
    inline I i() const {return values[1];}
    inline I j() const {return values[2];}
    inline I ai() const {return values[3];}
    inline I ci() const {return values[4];}
    inline I aj() const {return values[5];}
    inline I cj() const {return values[6];}
    
    // Equality operator for unordered_set
    bool operator==(const ambiguity& other) const {
        return (hash_ == other.hash_) and (values == other.values);
    }
};

// Custom hash function
template<typename I>
struct ambiguity_hash {
    size_t operator()(const ambiguity<I>& a) const { return a.hash_; }
};

}

#endif // AMBIGUITY_H


