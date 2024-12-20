#pragma once

#include <algorithm>
#include <concepts>
#include <cstdint>
#include <limits>
#include <span>
#include <vector>

#include <boost/align/align_down.hpp>
#include <boost/align/align_up.hpp>
#include <boost/container/small_vector.hpp>
#include <boost/multiprecision/gmp.hpp>
#include <boost/unordered/unordered_flat_map.hpp>

#include "ambiguity.hpp"
#include "kommunopp.hpp"
#include "monomial_trie.hpp"
#include "signal_statistics.hpp"




namespace kommunopp {

template<internal::metadata_concept MM = internal::metadata<uint8_t>,
         internal::metadata_concept PM = internal::metadata<uint8_t>,
         internal::value_concept V = uint8_t,
         typename I = uint32_t,
         typename C = boost::multiprecision::gmp_rational>
struct impl {
  using var = V;
  using idx = I;
  using coefficient = C;
  using monomial_store = internal::monomial_store<MM, V, I>;
  using polynomial_store = internal::polynomial_store<PM, MM, V, I, C>;
  // using monomial_trie = 
};

template<size_t N>
struct metadata_monomial {
  uint8_t length = 0;
  uint16_t t[N];
};

template<size_t N>
struct metadata_polynomial {
  uint8_t length = 0;
  uint16_t t[N];
};

template<size_t N,
        internal::value_concept V = uint8_t,
        typename I = uint32_t,
        typename C = boost::multiprecision::gmp_rational>
struct data {
    using impl = kommunopp::impl<metadata_monomial<N>, metadata_polynomial<1>, V, I, C>;
    using var = V;
    using idx = I;
    using coefficient = C;
    using monomial_store = impl::monomial_store;
    using polynomial_store = impl::polynomial_store;
    // using monomial_trie = impl::monomial_trie;
    
    public:
        monomial_store mons;
        polynomial_store poly;
        // monomial_trie lm;
        // monomial_trie lm_reversed;
    
    data(): mons(), poly(mons) {}
    
};
    




}
