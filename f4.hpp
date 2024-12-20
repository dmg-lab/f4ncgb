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
#include <boost/unordered_set.hpp>

#include "ambiguity.hpp"
#include "kommunopp.hpp"
#include "monomial_trie.hpp"
#include "signal_statistics.hpp"

namespace kommunopp {

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

//------------------------------------------------------------------------------

template<size_t N,
         internal::value_concept V = uint8_t,
         typename I = uint32_t,
         typename C = boost::multiprecision::gmp_rational>
struct f4 {
  using MM = metadata_monomial<N>;
  using PM = metadata_polynomial<1>;
  using coefficient = C;
  using monomial_store = internal::monomial_store<MM, V, I>;
  using polynomial_store = internal::polynomial_store<PM, MM, V, I, C>;
  using monomial_trie = monomial_trie<MM, V, I>;

  using monomial = std::span<const V>;
  using mon_id = I;
  using poly_id = I;

  using ambiguity = ambiguity<I>;
  using amb_hash = ambiguity_hash<I>;
  using crit_pair = std::pair<I, I>;

  public:
  monomial_store mons;
  polynomial_store poly;
  monomial_trie lm;
  monomial_trie lm_reversed;
  std::vector<poly_id> basis;
  std::map<size_t, boost::unordered_set<ambiguity>> amb;
  std::set<crit_pair> crit_pairs;
  std::map<mon_id, poly_id> lm_to_poly;

  f4()
    : mons()
    , poly(mons)
    , lm()
    , lm_reversed() {}

  //------------------------------------------------------------------------------

  std::vector<poly_id> compute_basis(size_t maxiter) {
    msg("Computing Gröbner basis");

    poly_id p1 = poly.add_polynomial({ { 1l, { 2, 3 } }, { 2l, { 2, 3 } } });
    poly_id p2 = poly.add_polynomial({ { 3l, { 3, 4 } }, { 4l, { 5, 6 } } });

    poly.print(std::cout);

    std::vector<poly_id> input = { p1, p2 };

    // add something at 0th position
    // so that index 0 remains free
    basis.push_back(0);

    // add input to critical pairs
    for(const auto& p : input)
      crit_pairs.insert(std::make_pair(p, p));

    // main loop
    size_t iter = 0;
    while(!amb.empty() and crit_pairs.empty() and iter < maxiter) {

      stage_crit_pairs();

      auto new_elements = reduction();

      update_basis(new_elements);
    }

    return basis;
  }

  //------------------------------------------------------------------------------
  void stage_crit_pairs() {
    if(amb.empty())
      return;

    auto minimal_amb = amb.begin();
    size_t d = minimal_amb->first;
    for(const auto& a : minimal_amb->second) {
      // crit_pair p = a.to_crit_pair();
      // crit_pairs.insert(p);
    }
    amb.erase(d);
  }
  //------------------------------------------------------------------------------
  void compute_ambiguities(mon_id i) {
    boost::unordered_set<ambiguity&, amb_hash> new_amb;

    // compute all overlaps ABC where a = AB
    new_amb = lm.compute_overlaps(mons, i);

    // compute all overlaps ABC where a = BC
    new_amb.merge(lm_reversed.compute_overlaps_reversed(mons, i));

    // compute all inclusions
    new_amb.merge(lm.compute_inclusions(mons, i));

    for(const auto& a : new_amb)
      amb[a.degree()].insert(a);
  }

  //------------------------------------------------------------------------------
  std::vector<poly_id>& reduction() {}

  //------------------------------------------------------------------------------
  void update_basis(std::vector<poly_id>& new_elements) {

    size_t n = basis.size();
    for(const poly_id& p_id : new_elements) {
      // update lm data
      mon_id m_id = poly.get_lm_id(p_id);
      lm_to_poly[m_id] = p_id;

      // update tries
      monomial m = mons[m_id];
      lm.insert_monomial(m);
      lm_reversed.insert_monomial_reversed(m);

      // update basis
      basis.push_back(p_id);
    }
  }
};
}
