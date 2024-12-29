#pragma once

#include <algorithm>
#include <chrono>
#include <concepts>
#include <cstdint>
#include <limits>
#include <span>
#include <string>
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
#include "linear_algebra.hpp"
#include "monomial_trie.hpp"
#include "parser.hpp"
#include "signal_statistics.hpp"
#include "sparse_rref/sparse_mat.h"
#include "sparse_rref/sparse_vec.h"

using namespace boost::multiprecision;

extern double amb_time, crit_pair_time, sym_pre_time, reduction_time,
  reset_time;
extern double overlap_time, inclusion_time, overlap_map_time;

std::chrono::time_point<std::chrono::high_resolution_clock> start;
std::chrono::time_point<std::chrono::high_resolution_clock> end;
std::chrono::duration<double> elapsed;

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

  template<size_t N, size_t Nvars,
         internal::value_concept V = uint8_t,
         typename I = uint32_t,
         typename C = boost::multiprecision::gmp_rational>
struct f4 {
  using MM = metadata_monomial<N>;
  using PM = metadata_polynomial<1>;
  using coefficient = C;
  using monomial_store = internal::monomial_store<MM, V, I>;
  using polynomial_store = internal::polynomial_store<PM, MM, V, I, C>;

  using monomial = std::span<const V>;
  using mon_id = I;
  using poly_id = I;

  using ambiguity = ambiguity<mon_id>;
  using amb_hash = ambiguity_hash<mon_id>;
  using crit_pair = std::pair<poly_id, poly_id>;

  using monomial_trie = monomial_trie<Nvars, V, I>;

  // Custom hash function for std::span
  struct monomial_hash {
    size_t operator()(const monomial& m) const {
      size_t seed = 0;
      for(auto it = m.begin(); it != m.end(); it++)
        boost::hash_combine(seed, *it);
      return seed;
    }
  };

  // Custom equality comparator for std::span
  struct monomial_equal {
    bool operator()(const monomial& m1, const monomial& m2) const {
      return m1.size() == m2.size()
             && std::equal(m1.begin(), m1.end(), m2.begin());
    }
  };

  public:
  monomial_store mons;
  polynomial_store poly;
  std::vector<poly_id> basis;
  std::map<size_t, std::unordered_set<ambiguity, amb_hash>> amb;
  std::set<crit_pair> crit_pairs;
  std::unordered_map<mon_id, poly_id> lm_to_poly;
  monomial_trie prefix_trie;
  monomial_trie suffix_trie;

  size_t maxdeg;

  f4()
    : mons()
    , poly(mons)
    , prefix_trie()
    , suffix_trie() {}

  //------------------------------------------------------------------------------
  inline void read_input(parser_context& context) {
    parse_res res = parse_rest_into_polynomial_store(context, poly);
  }

  //------------------------------------------------------------------------------

  std::vector<poly_id> compute_basis(size_t maxiter, size_t maxdeg_ = INT_MAX) {
    msg("Computing Gröbner basis");

    maxdeg = maxdeg_;

    // add something at 0th position
    // so that index 0 remains free
    basis.push_back(0);
    poly_id p1
      = poly.add_polynomial({ { 1l, { 3, 2, 3 } }, { -1l, { 2, 1, 2 } } });
    poly_id p2
      = poly.add_polynomial({ { 1l, { 3, 1, 2 } }, { -1l, { 1, 2, 1 } } });
    poly_id p3
      = poly.add_polynomial({ { 1l, { 3, 1, 3 } }, { -1l, { 2, 3, 1 } } });
    poly_id p4 = poly.add_polynomial({ { 1l, { 3, 3, 3 } },
                                       { 1l, { 2, 2, 2 } },
                                       { 1l, { 1, 2, 3 } },
                                       { 1l, { 1, 1, 1 } } });
    std::vector<poly_id> input = { p1, p2, p3, p4 };

    /* poly_id p1 = poly.add_polynomial({ { 1l, { 1, 1, 1 } }, { -1l, {} } });
     */
    /* poly_id p2 = poly.add_polynomial({ { 1l, { 2, 2, 2 } }, { -1l, {} } });
     */
    /* poly_id p3 = poly.add_polynomial( */
    /*   { { 1l, { 2, 1, 2, 1, 1, 2, 1, 2, 1, 1 } }, { -1l, {} } }); */
    /* std::vector<poly_id> input = { p1, p2, p3 }; */

    /* poly_id p1 = poly.add_polynomial({ { 1l, { 3, 2 } }, { 1l, { 2, 1 } } });
     */
    /* poly_id p2 = poly.add_polynomial({ { 1l, { 3, 3 } }, */
    /*                                    { 1l, { 3, 2 } }, */
    /*                                    { -1l, { 2, 3 } }, */
    /*                                    { -1l, { 2, 2 } } }); */
    /* std::vector<poly_id> input = { p1, p2 }; */

    // add input to critical pairs
    for(const auto& p : input) {
      crit_pair c(p, p);
      crit_pairs.insert(c);
    }

    // main loop
    size_t iter = 0;
    while((!amb.empty() or !crit_pairs.empty()) and iter < maxiter) {

      start = std::chrono::high_resolution_clock().now();
      stage_crit_pairs();
      end = std::chrono::high_resolution_clock().now();
      elapsed = end - start;
      crit_pair_time += elapsed.count();

      msg("Reducing %d critical pairs.", crit_pairs.size());
      std::vector<poly_id> new_elements = reduction();
      msg("Adding %d new elements to basis.", new_elements.size());

      update_basis_and_amb(new_elements);

      msg("Iteration %d has finished. Basis has now %d elements",
          iter,
          basis.size() - 1);
      iter++;
    }

    /* for (auto i: basis){ */
    /*   auto coeffs = poly.get_coefficients(i); */
    /*   size_t k = 0; */
    /*   for (auto it = poly[i].begin(); it != poly[i].end(); it++) { */
    /*     std::cout << mpq_rational(coeffs[k++]) << "*"; */
    /*     mons.print_monomial(*it, std::cout); */
    /*     std::cout << " "; */
    /*   } */
    /*   std::cout << "\n"; */
    /*   } */

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
    for(const auto& a : minimal_amb->second) {
      crit_pairs.insert(to_crit_pair(a));
    }
    amb.erase(d);
  }
  //------------------------------------------------------------------------------
  void compute_ambiguities(mon_id i) {
    std::unordered_set<ambiguity, amb_hash> new_amb;

    monomial m = mons[i];
    monomial ab, b, bc, abc;
    std::vector<std::pair<mon_id, size_t>> overlaps;
    std::vector<std::pair<mon_id, size_t>> inclusions;

    auto s = std::chrono::high_resolution_clock().now();
    prefix_trie.overlaps_and_inclusions(m, overlaps, inclusions);

    // overlaps with m = AB
    ab = m;
    // last k elements of m = AB form overlap B
    for(auto [j, k] : overlaps) {
      bc = mons[j];
      I d = ab.size() + bc.size() - k;
      if(d > maxdeg)
        continue;
      I aj = mons.getid(ab.first(ab.size() - k));
      I ci = mons.getid(bc.last(bc.size() - k));
      ambiguity a(d, i, j, 0, ci, aj, 0);
      new_amb.insert(a);
    }
    overlaps.clear();

    // overlaps with m = BC
    suffix_trie.overlaps_rev(m, overlaps);
    bc = m;
    // k determines where B starts in m = BC
    for(auto [j, k] : overlaps) {
      ab = mons[j];
      I d = ab.size() + bc.size() - k;
      if(d > maxdeg)
        continue;
      I ai = mons.getid(ab.first(ab.size() - k));
      I cj = mons.getid(bc.last(bc.size() - k));
      ambiguity a(d, j, i, 0, cj, ai, 0);
      new_amb.insert(a);
    }
    auto e = std::chrono::high_resolution_clock().now();
    std::chrono::duration<double> elapsed = e - s;
    overlap_time += elapsed.count();

    s = std::chrono::high_resolution_clock().now();
    // inclusions with m = ABC
    I d = m.size();
    // k determines where B starts in m = ABC
    for(auto [j, k] : inclusions) {
      if(i == j)
        continue;
      b = mons[j];
      I aj = mons.getid(m.first(k));
      I cj = mons.getid(m.last(d - k - b.size()));
      ambiguity a(d, j, i, aj, cj, 0, 0);
      new_amb.insert(a);
    }
    inclusions.clear();

    // inclusions with m = B
    b = m;
    prefix_trie.inclusions(m, inclusions);
    // last k elements in ABC form C
    for(auto [j, k] : inclusions) {
      if(i == j)
        continue;
      abc = mons[j];
      I d = abc.size();
      if(d > maxdeg)
        continue;
      I ai = mons.getid(abc.first(abc.size() - b.size() - k));
      I ci = mons.getid(abc.last(k));
      ambiguity a(d, i, j, ai, ci, 0, 0);
      new_amb.insert(a);
    }
    e = std::chrono::high_resolution_clock().now();
    elapsed = e - s;
    inclusion_time += elapsed.count();

    for(const auto& a : new_amb)
      amb[a.degree()].insert(a);
  }

  //------------------------------------------------------------------------------
  std::unordered_set<poly_id> symbolic_preprocessing() {
    std::unordered_set<mon_id> todo;
    std::unordered_set<mon_id> done;
    std::unordered_set<poly_id> rows;

    for(const auto& [f, g] : crit_pairs) {
      auto mon_it = poly[f];
      done.insert(*mon_it.begin());
      todo.insert(++mon_it.begin(), mon_it.end());
      mon_it = poly[g];
      done.insert(*mon_it.begin());
      todo.insert(++mon_it.begin(), mon_it.end());
      rows.insert(f);
      rows.insert(g);
    }

    crit_pairs.clear();

    while(!todo.empty()) {
      mon_id m = *todo.begin();
      todo.erase(todo.begin());
      done.insert(m);
      poly_id reducer = find_reducer(m);
      if(!reducer)
        continue;
      assert(m == poly.get_lm_id(reducer));
      rows.insert(reducer);
      for(auto it = poly[reducer].begin(); it < poly[reducer].end(); it++) {
        if(!done.count(*it))
          todo.insert(*it);
      }
    }
    return rows;
  }

  //------------------------------------------------------------------------------
  poly_id find_reducer(mon_id m, bool strategy = false) {

    auto reducers = prefix_trie.divisors(mons[m]);
    if(reducers.empty())
      return 0;

    std::pair<mon_id, size_t> match;
    // strategy 1 : the last one
    if(strategy)
      match = *std::max_element(reducers.begin(), reducers.end());
    // strategy 2 : the one with smallest lm
    else if(true)
      match = *std::max_element(
        reducers.begin(), reducers.end(), [this](auto a, auto b) {
          return this->mons.cmp(b.first, a.first);
        });
    // strategy 3 : the one with largest lm
    else
      match = *std::max_element(
        reducers.begin(), reducers.end(), [this](auto a, auto b) {
          return this->mons.cmp(a.first, b.first);
        });

    monomial mm = mons[m];
    monomial lm = mons[match.first];
    mon_id a = mons.getid(mm.first(match.second));
    mon_id b = mons.getid(mm.last(mm.size() - match.second - lm.size()));

    poly_id res = poly.multiply_front_and_back(a, lm_to_poly[match.first], b);

    return res;
  }
  //------------------------------------------------------------------------------
  std::vector<poly_id> identify_new_elements(
    std::vector<std::pair<size_t, size_t>>& idxs,
    std::vector<C>& coeffs,
    std::vector<mon_id>& columns) {

    std::vector<poly_id> res;
    std::vector<std::pair<C, mon_id>> p;

    std::unordered_map<size_t, std::vector<std::pair<size_t, C>>> poly_map;
    size_t k = 0;
    for(auto [i, j] : idxs) {
      poly_map[i].emplace_back(j, coeffs[k++]);
    }

    for(auto it = poly_map.begin(); it != poly_map.end(); it++) {
      auto poly_data = it->second;
      size_t j0 = poly_data[0].first;
      // test if leading monomial already exists
      if(!prefix_trie.divisors(mons[columns[j0]]).empty()) {
        continue;
      }
      // new leading monomial, actually make polynomial
      for(auto [j, c] : poly_data)
        p.emplace_back(c, columns[j]);

      res.push_back(poly.add_polynomial(p));
      p.clear();
    }
    return res;
  }
  //------------------------------------------------------------------------------
  std::vector<poly_id> reduction() {
    // symbolic preprocessing

    std::cout << "Symbolic preprocessing\n";
    start = std::chrono::high_resolution_clock::now();
    std::unordered_set<poly_id> rows = symbolic_preprocessing();
    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    sym_pre_time += elapsed.count();
    std::cout << "Symbolic preprocessing done\n";

    // make columns
    // columns are sorted in DESCENDING order
    std::unordered_set<mon_id> col_set;
    for(const auto r : rows) {
      auto p = poly[r];
      col_set.insert(p.begin(), p.end());
    }
    std::vector<mon_id> columns(col_set.begin(), col_set.end());
    auto cmp
      = [this](const mon_id a, const mon_id b) { return this->mons.cmp(b, a); };
    std::sort(columns.begin(), columns.end(), cmp);

    // set up matrix
    std::cout << "Setting  up matrix\n";
    sfmpz_mat_t mat;
    set_up_matrix(mat, rows, columns);
    std::cout << "Setting up matrix done\n";

    // reduction
    start = std::chrono::high_resolution_clock::now();
    auto [idxs, entries] = multimodular_rref(mat);
    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    reduction_time += elapsed.count();

    std::cout << "RRef done\n";

    // identify new elements
    std::vector<poly_id> new_elements
      = identify_new_elements(idxs, entries, columns);

    std::cout << "Computed new elements\n";

    return new_elements;
  }
  //------------------------------------------------------------------------------
  inline void get_common_denom(fmpz_t denom, fmpz_t tmp, std::span<C>& coeffs) {
    fmpz_set_ui(denom, 1);
    for(auto& c : coeffs) {
      fmpz_set_mpz(tmp, mpq_denref(c.data()));
      fmpz_lcm(denom, denom, tmp);
    }
  }
  //------------------------------------------------------------------------------
  void set_up_matrix(sfmpz_mat_t mat,
                     std::unordered_set<poly_id>& rows,
                     std::vector<mon_id>& columns) {

    std::unordered_map<mon_id, size_t> col_to_id;
    size_t i = 0;
    for(auto c : columns)
      col_to_id[c] = i++;

    msg("Setting up matrix of size (%d, %d)", rows.size(), columns.size());

    // initialize matrix
    sparse_mat_init(mat, rows.size(), columns.size());

    // set all entries
    fmpz_t denom;
    fmpz_t tmp;
    fmpz_t c;
    fmpz_init(c);
    fmpz_init(denom);
    fmpz_init(tmp);
    i = 0;

    for(auto r : rows) {
      auto row = sparse_mat_row(mat, i++);
      std::span<C> coeffs = poly.get_coefficients(r);

      // compute common denominator so that we can normalize row
      get_common_denom(denom, tmp, coeffs);

      size_t k = 0;
      for(auto j = poly[r].begin(); j != poly[r].end(); j++) {
        auto cc = coeffs[k++].data();
        fmpz_set_mpz(tmp, mpq_denref(cc));
        assert(fmpz_divisible(denom, tmp));
        fmpz_divexact(tmp, denom, tmp);
        fmpz_set_mpz(c, mpq_numref(cc));
        fmpz_mul(c, c, tmp);
        _sparse_vec_set_entry(row, col_to_id[*j], c);
      }
    }

    fmpz_clear(denom);
    fmpz_clear(tmp);
    fmpz_clear(c);
  }

  //------------------------------------------------------------------------------
  void update_basis_and_amb(std::vector<poly_id>& new_elements) {

    for(const poly_id p_id : new_elements) {
      // update lm data
      mon_id m_id = poly.get_lm_id(p_id);
      lm_to_poly[m_id] = p_id;

      // update tries
      monomial m = mons[m_id];
      prefix_trie.insert(m, m_id);
      suffix_trie.insert_rev(m, m_id);

      // compute ambiguities
      start = std::chrono::high_resolution_clock().now();
      compute_ambiguities(m_id);
      end = std::chrono::high_resolution_clock().now();
      elapsed = end - start;
      amb_time += elapsed.count();

      // update basis
      basis.push_back(p_id);
    }
  }
};
}
