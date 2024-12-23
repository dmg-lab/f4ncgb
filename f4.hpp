#pragma once

#include <algorithm>
#include <concepts>
#include <cstdint>
#include <limits>
#include <span>
#include <utility>
#include <vector>

#include <boost/align/align_down.hpp>
#include <boost/align/align_up.hpp>
#include <boost/container/small_vector.hpp>
#include <boost/multiprecision/gmp.hpp>
#include <boost/unordered/unordered_flat_map.hpp>
#include <unordered_set>

#include "ambiguity.hpp"
#include "gmp.h"
#include "kommunopp.hpp"
#include "monomial_trie.hpp"
#include "signal_statistics.hpp"
#include "sparse_rref/sparse_mat.h"
#include "linear_algebra.hpp"

using namespace boost::multiprecision;

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

  using ambiguity = ambiguity<mon_id>;
  using amb_hash = ambiguity_hash<mon_id>;
  using crit_pair = std::pair<poly_id, poly_id>;

  public:
  monomial_store mons;
  polynomial_store poly;
  monomial_trie lm;
  monomial_trie lm_reversed;
  std::vector<poly_id> basis;
  std::map<size_t, std::unordered_set<ambiguity, amb_hash>> amb;
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

    poly_id p1 = poly.add_polynomial({ { 5l, { 2, 3 } }, { 14l, { 2, 3 } } });
    poly_id p2 = poly.add_polynomial({ { 3l, { 3, 4 } }, { 4l, { 5, 6 } } });

    C src = (poly.get_coefficients(p1)[0]);
    C src2 = (poly.get_coefficients(p1)[1]);
    gmp_int i(src);
    gmp_int j(src2);

    sparse_mat_t<gmp_int> mat;
    sparse_mat_init(mat, 5, 5);

    auto row = sparse_mat_row(mat, 0);
    _sparse_vec_set_entry(row, 3, &i);
    _sparse_vec_set_entry(row, 1, &j);
    auto [idxs, entries] = multimodular_rref(mat);
    sparse_mat_clear(mat);

    std::vector<poly_id> input = { p1, p2 };

    // add something at 0th position
    // so that index 0 remains free
    basis.push_back(0);

    // add input to critical pairs
    for(const auto& p : input) {
      crit_pair c(p, p);
      crit_pairs.insert(c);
    }

    // main loop
    size_t iter = 0;
    while(!amb.empty() and !crit_pairs.empty() and iter < maxiter) {

      stage_crit_pairs();

      msg("Reducing %d critical pairs", crit_pairs.size());
      auto new_elements = reduction();

      update_basis_and_amb(new_elements);

      iter++;
    }

    return basis;
  }

  //------------------------------------------------------------------------------
  inline crit_pair to_crit_pair(const ambiguity& a) {
    poly_id i = lm_to_poly[a.i()];
    poly_id j = lm_to_poly[a.j()];

    mon_id ai = a.ai();
    mon_id ci = a.ci();
    mon_id aj = a.aj();
    mon_id cj = a.cj();

    poly_id f = poly.multiply_front_and_back(ai, i, ci);
    poly_id g = poly.multiply_front_and_back(aj, j, cj);
    crit_pair c(f, g);
    return c;
  }
  //------------------------------------------------------------------------------

  inline void stage_crit_pairs() {
    if(amb.empty())
      return;

    auto minimal_amb = amb.begin();
    size_t d = minimal_amb->first;
    for(const auto& a : minimal_amb->second)
      crit_pairs.insert(to_crit_pair(a));
    amb.erase(d);
  }
  //------------------------------------------------------------------------------
  void compute_ambiguities(mon_id i) {
    std::unordered_set<ambiguity, amb_hash> new_amb;

    // compute all overlaps ABC where a = AB
    new_amb = lm.compute_overlaps(mons, i);

    // compute all overlaps ABC where a = BC
    new_amb.merge(lm_reversed.compute_overlaps_reversed(mons, i));

    // compute all inclusions
    new_amb.merge(lm.compute_inclusions(mons, i));

    for(const auto& a : new_amb)
      amb[a.degree()].insert(a);
  }

  void symbolic_preprocessing() {
    std::set<mon_id> todo;
    std::set<mon_id> done;
    std::set<poly_id> reducers;
    for(const auto& [f, g] : crit_pairs) {
      auto mons = poly.get_monomial_ids(f);
      done.insert(mons[0]);
      todo.insert(mons.begin() + 1, mons.end());
      mons = poly.get_monomial_ids(g);
      done.insert(mons[0]);
      todo.insert(mons.begin() + 1, mons.end());
    }
   while (!todo.empty()) {
     //   mon_id m = todo.pop();


} 


  }

  //------------------------------------------------------------------------------
  std::vector<poly_id>& reduction() {
    // symbolic preprocessing
    symbolic_preprocessing();
    // set up matrix

    // reduction

    // identify new elements
  }

  //------------------------------------------------------------------------------
  void update_basis_and_amb(std::vector<poly_id>& new_elements) {

    for(const poly_id& p_id : new_elements) {
      // update lm data
      mon_id m_id = poly.get_lm_id(p_id);
      lm_to_poly[m_id] = p_id;

      // compute ambiguities
      compute_ambiguities(m_id);

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
