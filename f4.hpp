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

#include "parser.hpp"
#include "ambiguity.hpp"
#include "gmp.h"
#include "kommunopp.hpp"
#include "linear_algebra.hpp"
#include "monomial_trie.hpp"
#include "signal_statistics.hpp"
#include "sparse_rref/sparse_mat.h"

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
  inline void read_input(parser_context & context) {
    parse_res res = parse_rest_into_polynomial_store(context, poly);
    }

  //------------------------------------------------------------------------------

  std::vector<poly_id> compute_basis(size_t maxiter) {
    msg("Computing Gröbner basis");

    // add something at 0th position
    // so that index 0 remains free
    basis.push_back(0);

    std::vector<poly_id> input(poly.begin(), poly.end());     
    update_basis_and_amb(input);

    mons.print(std::cout);
    poly.print(std::cout);

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

  //------------------------------------------------------------------------------
  std::unordered_set<poly_id> symbolic_preprocessing() {
    std::unordered_set<mon_id> todo;
    std::unordered_set<mon_id> done;
    std::unordered_set<poly_id> rows;

    msg("============ Sym Pre ===================");

    for(const auto& [f, g] : crit_pairs) {
      auto mons = poly.get_monomial_ids(f);
      done.insert(mons[0]);
      todo.insert(mons.begin() + 1, mons.end());
      mons = poly.get_monomial_ids(g);
      done.insert(mons[0]);
      todo.insert(mons.begin() + 1, mons.end());
      rows.insert(f);
      rows.insert(g);
    }
    while(!todo.empty()) {
      die(-1,"");
      std::cout << "Todo length = " << todo.size() << "\n";
      for (auto c : todo)
        std::cout << (int)c << ", ";
      std::cout << "\n";
      
      mon_id m = *todo.begin();
      todo.erase(todo.begin());
      done.insert(m);
      poly_id reducer = find_reducer(m);
      if(!reducer)
        continue;
      assert(m == poly.get_lm_id(reducer));

      std::cout << "Found reducer " << reducer << "\n";

      mons.print(std::cout);
      std::cout << "---------------------------\n";
      poly.print(std::cout);
      
      rows.insert(reducer);
      for(auto rm : poly.get_monomial_ids(reducer))
        if(!done.count(rm))
          todo.insert(rm);
    }

        msg("============ Sym Pre finished  ===================");

    return rows;
  }

  //------------------------------------------------------------------------------
  poly_id find_reducer(mon_id m, bool strategy = false) {

        msg("============ Find reducer  ===================");
    
    auto reducers = lm.compute_divisors(mons, m);
    if(reducers.empty())
      return 0;

    monomial mm = mons[m];

    poly_id g;
    // strategy 1 : the last one
    if(strategy)
      g = *std::max_element(reducers.begin(), reducers.end());
    // strategy 2 : the one with smallest lm
    else if(true)
      g = *std::max_element(
        reducers.begin(), reducers.end(), [this](mon_id a, mon_id b) {
          return this->mons.cmp(b, a);
        });
    // strategy 3 : the one with largest lm
    else
      g = *std::max_element(
        reducers.begin(), reducers.end(), [this](mon_id a, mon_id b) {
          return this->mons.cmp(a, b);
        });

    monomial lm_g = poly.get_lm(g);

    auto it = std::search(mm.begin(), mm.end(), lm_g.begin(), lm_g.end());
    assert(it != lm_g.end());
    mon_id a = mons.getid(std::span(mm.begin(), it));
    mon_id b = mons.getid(std::span(it + lm_g.size(), mm.end()));

    poly_id res = poly.multiply_front_and_back(a, g, b);

    msg("============ Find reducer finished  ===================");

    
    return res;
  }
  //------------------------------------------------------------------------------
  struct compare_lm {

    monomial_store * mons;
    polynomial_store * poly;

    compare_lm(monomial_store *s,  polynomial_store * p) : mons(s), poly(p) {}
    
    bool operator()(const poly_id a, const poly_id b) {
      return mons->cmp(poly->get_lm_id(a), poly->get_lm_id(b));
    }
  };
  //------------------------------------------------------------------------------
  std::vector<poly_id>& reduction() {
    // symbolic preprocessing
    std::unordered_set<poly_id> rows = symbolic_preprocessing();

    // make columns
    std::set<mon_id, compare_lm> columns (compare_lm(&mons,&poly));
    for(const auto r : rows) { auto m = poly.get_monomial_ids(r);
      columns.insert(m.begin(), m.end());
    }

    for(auto c : columns)
      mons.print_monomial(c, std::cout);
    // set up matrix
    

    // reduction

    // identify new elements
  }

  //------------------------------------------------------------------------------
  void update_basis_and_amb(std::vector<poly_id>& new_elements) {

    for(const poly_id p_id : new_elements) {
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
