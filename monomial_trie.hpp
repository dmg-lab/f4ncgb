#ifndef MONOMIAL_TRIE_H
#define MONOMIAL_TRIE_H

#include <algorithm>
#include <concepts>
#include <cstdint>
#include <limits>
#include <span>
#include <vector>
#include <unordered_set>

#include "./aho_corasick/src/aho_corasick/aho_corasick.hpp"

#include "kommunopp.hpp"
#include "ambiguity.hpp"

namespace kommunopp {

template<internal::metadata_concept M = internal::metadata<uint8_t>,
         internal::value_concept V = uint8_t,
         typename I = uint32_t>
class monomial_trie {
    using monomial = std::span<const V>;
    using monomial_store = internal::monomial_store<M, V, I>;
    using ambiguity = ambiguity<I>;
    using amb_hash = ambiguity_hash<I>;
        
    public:
    aho_corasick::trie<V> T;
        
    monomial_trie(){}
    ~monomial_trie() = default;

    inline void insert_monomial(const monomial& m) {
        std::vector<V> vec(m.begin(), m.end());
        T.insert(vec);
    }
    // -----------------------------------------------------------------
    inline void insert_monomial_reversed(const monomial& m) {
        std::vector<V> vec(m.rbegin(), m.rend());
        T.insert(vec,true);
    }
    // -----------------------------------------------------------------
    std::unordered_set<I> compute_divisors(monomial_store& s, I a) {
    monomial m = s[a];
    auto results = T.parse_text(m.begin(), m.end());        
    
    std::unordered_set<I> divisors;
    divisors.reserve(results.size());
    
    for (auto r : results ) {
        auto vec = r.get_keyword();
        divisors.insert( s.getid(vec) );
    }    
    return divisors;
    }
    // -----------------------------------------------------------------
    std::unordered_set<ambiguity, amb_hash> compute_overlaps(monomial_store& s, I i) {
    /*
    * Find all overlaps ABC with a = AB
    */
    monomial ab = s[i];
    std::unordered_set<ambiguity, amb_hash> overlaps;
    
    for (I id = 1; id < ab.size(); id++) {     
        auto results = T.prefixes(ab.begin() + id, ab.end());                 
        
        // compute the ambiguities
        for (auto &  r : results) {                
            std::vector<V> v = r.get_keyword();
            monomial bc (v.begin(), v.size());
            
            I d = id + v.size();
            I j = s.getid(bc);
            I aj = s.getid(ab.first(id));
            I ci = s.getid(bc.last(r.get_end()));
            
            ambiguity a (d,i,j,0,ci,aj,0);
            overlaps.insert(a);
        }      
    }               
    return overlaps;
    }
    // -----------------------------------------------------------------
    std::unordered_set<ambiguity, amb_hash> compute_overlaps_reversed(monomial_store& s, I i) {
    /*
    * Find all overlaps ABC with a = BC
    */
    monomial bc = s[i];
    std::unordered_set<ambiguity, amb_hash> overlaps;
    
    std::cout << "Computing REVERSED overlaps\n";
    
    for (I id = 1; id < bc.size(); id++) {     
        auto results = T.prefixes(bc.rbegin() + id, bc.rend());                 
        
        // compute the ambiguities
        for (auto &  r : results) {                
            std::vector<V> v = r.get_keyword();
            monomial ab (v.begin(), v.size());
            
            std::cout << "AB = ";
            for (const auto & ch : ab)
                std::cout << static_cast<int>(ch) << ", ";
            std::cout << "\n";
            
            
            I d = id + v.size();
            I j = s.getid(ab);
            I ai = s.getid(ab.first(r.get_end()));
            I cj = s.getid(bc.last(id));
            
            ambiguity a (d,j,i,0,cj,ai,0);
            overlaps.insert(a);
        }      
    }               
    return overlaps;
}
    // -----------------------------------------------------------------    
    std::unordered_set<ambiguity, amb_hash> compute_inclusions(monomial_store& s, I i) {
    monomial m = s[i]; 

    // Find all inclusions ABC with m = ABC
    std::unordered_set<ambiguity, amb_hash> inclusions;
    auto results = T.parse_text(m.begin(), m.end());
    
    for (const auto & r : results) {
        std::vector<V> v = r.get_keyword();
        monomial b (v.begin(), v.size());
        // pos say where A ends in ABC
        I pos = r.get_start();                                   
        I d = m.size();
        I j = s.getid(b);            
        I aj = s.getid(m.first(pos));
        I cj = s.getid(m.last(d - pos - b.size()));

        ambiguity a (d,j,i,aj,cj,0,0);
        inclusions.insert(a);
    } 
            
    // Find all inclusions ABC with m = B   
    results = T.inclusions(m.begin(), m.end()); 
            
    for (const auto & r : results) {
        std::vector<V> v = r.get_keyword();            
        monomial abc (v.begin(), v.size());
        // pos say how long C is in ABC
        I pos = r.get_end();                    
        I d = abc.size();
        I j = s.getid(abc);        
        I ai = s.getid(abc.first(d - pos - m.size()));
        I ci = s.getid(abc.last(pos));
        
        ambiguity a (d,i,j,ai,ci,0,0);
        inclusions.insert(a);
    }    
                
    return inclusions;
}
};
}

#endif // MONOMIAL_TRIE_H

// -----------------------------------------------------------------
//     const std::unordered_set<ambiguity,amb_hash> compute_ambiguities(
//                                                 const monomial_trie& prefix, 
//                                                 const monomial_trie& suffix, 
//                                                 I i) {
//         std::unordered_set<ambiguity, amb_hash> amb;
//         
//         // compute all overlaps ABC where a = AB
//         amb = prefix.compute_overlaps(*this, i);
//         
//         // compute all overlaps ABC where a = BC
//         amb.merge( suffix.compute_overlaps_reversed(*this, i) );
//         
//         // compute all inclusions
//         amb.merge( prefix.compute_inclusions(*this, i) );
//                 
//         return amb;                                             
//     }