#pragma once

#include <algorithm>
#include <chrono>
#include <concepts>
#include <coroutine>
#include <cstdint>
#include <limits>
#include <map>
#include <memory>
#include <ostream>
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
extern double overlap_time, inclusion_time, new_elements_time;

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

template<size_t N,
         size_t Nvars,
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

  size_t characteristic = 0;
  size_t iter = 0;
  size_t maxiter = UINT_MAX;
  size_t maxdeg = UINT_MAX;
  size_t threads = 1;

  f4(size_t prime_, size_t maxiter_, size_t maxdeg_, size_t threads_)
    : mons()
    , poly(mons)
    , prefix_trie()
    , suffix_trie()
    , characteristic(prime_)
    , maxiter(maxiter_)
    , maxdeg(maxdeg_)
    , threads(threads_) {}

  //------------------------------------------------------------------------------
  inline parse_res read_input(parser_context& context) {
    parse_res res = parse_rest_into_polynomial_store(context, poly);
    // Print in case of doubt.
    // context.to_msolve(std::cout, poly);
    return res;
  }
  //------------------------------------------------------------------------------
  void interreduce_and_add_to_basis(std::vector<poly_id> input) {

    msg("Linearly interreducing input of size %d.", input.size());
    for(const auto& p : input) {
      crit_pair c(p, p);
      crit_pairs.insert(c);
    }
    start = std::chrono::high_resolution_clock().now();
    stage_crit_pairs();
    end = std::chrono::high_resolution_clock().now();
    elapsed = end - start;
    crit_pair_time += elapsed.count();

    std::vector<poly_id> new_elements = reduction(true);
    msg("Adding %d input elements to basis.", new_elements.size());

    update_basis_and_amb(new_elements);
  }
  //------------------------------------------------------------------------------

  std::vector<poly_id> compute_basis() {
    // add something at 0th position
    // so that index 0 remains free
    basis.push_back(0);
    std::vector<poly_id> input(poly.begin() + 1, poly.end());

    // add input to critical pairs
    interreduce_and_add_to_basis(input);

    // main loop
    iter = 0;
    while((!amb.empty() or !crit_pairs.empty()) and iter <= maxiter) {

      start = std::chrono::high_resolution_clock().now();
      stage_crit_pairs();
      end = std::chrono::high_resolution_clock().now();
      elapsed = end - start;
      crit_pair_time += elapsed.count();

      msg("Reducing %d critical pairs.", crit_pairs.size());
      std::vector<poly_id> new_elements = reduction();
      msg("Adding %d new elements to basis.", new_elements.size());

      update_basis_and_amb(new_elements);

      msg("==== Iteration %d has finished. Basis has now %d elements ====",
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
  /*
   * assumes that self.i() == other.i()
   * returns
   *   -1 if self is not divisible by other
   *   0 if self == other
   *   1 if self is properly divisble by other
   **/
  int inline divisible_by(const ambiguity& self, const ambiguity& other) const {
    assert(self.i() == other.i());

    std::span<const V> s_ai = mons[self.ai()];
    std::span<const V> s_ci = mons[self.ci()];
    std::span<const V> o_ai = mons[other.ai()];
    std::span<const V> o_ci = mons[other.ci()];

    if(endswith(s_ai, o_ai) and startswith(s_ci, o_ci))
      if(s_ai.size() == o_ai.size() and s_ci.size() == o_ci.size())
        return 0;
      else
        return 1;
    else
      return -1;
  }
  //------------------------------------------------------------------------------
  void gebauer_moeller(std::unordered_set<ambiguity, amb_hash>& new_amb) {
    // first index is always the newer polynomial

    auto cmp = [this](ambiguity& a, ambiguity& b) {
      if(a.degree() != b.degree())
        return a.degree() < b.degree();
      if(a.j() != b.j())
        return a.j() < b.j();
      return this->mons.cmp(a.aj(), b.aj());
    };

    std::vector<ambiguity> tmp(new_amb.begin(), new_amb.end());
    std::sort(tmp.begin(), tmp.end(), cmp);

    bool* to_remove = new bool[tmp.size()];
    for(size_t i = 0; i < tmp.size(); i++)
      to_remove[i] = false;

    for(size_t i = 0; i < tmp.size(); i++) {
      auto a = tmp[i];
      if(to_remove[i]) {
        new_amb.erase(a);
        continue;
      }
      // use a to remove other ambiguities
      for(size_t j = i + 1; j < tmp.size(); j++) {
        if(to_remove[j])
          continue;
        int d = divisible_by(tmp[j], a);
        // if divisible, remove
        // ordering ensures that other conditions are satisfied
        if(d >= 0)
          to_remove[j] = true;
      }
    }
    delete[] to_remove;
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
      ambiguity a(d, i, j, ai, 0, 0, cj);
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
      ambiguity a(d, i, j, 0, 0, aj, cj);
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

    gebauer_moeller(new_amb);

    for(const auto& a : new_amb)
      amb[a.degree()].insert(a);
  }

  //------------------------------------------------------------------------------
  std::unordered_set<poly_id> symbolic_preprocessing() {
    std::unordered_set<mon_id> todo;
    std::unordered_set<mon_id> done;
    std::unordered_set<poly_id> rows;

    for(const auto& [f, g] : crit_pairs) {
      // add monomials to corresponding sets
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

    size_t k = 0;
    while(!todo.empty()) {
      auto s = std::chrono::high_resolution_clock().now();
      mon_id m = *todo.begin();
      todo.erase(todo.begin());
      done.insert(m);
      std::vector<poly_id> reducer = find_reducer(m);
      if(reducer.size() != 3)
        continue;

      auto a = reducer[0];
      auto g = reducer[1];
      auto b = reducer[2];

      poly_id agb = poly.multiply_front_and_back(a, lm_to_poly[g], b);
      assert(m == poly.get_lm_id(agb));
      rows.insert(agb);
      for(auto mm : poly[agb]) {
        if(!done.count(mm))
          todo.insert(mm);
      }
      auto e = std::chrono::high_resolution_clock().now();
      std::chrono::duration<double> el = e - s;

      auto n = todo.begin();
      while(n != todo.end()) {
        if(mons.is_divisible(*n, g)) {
          k++;
          auto s = std::chrono::high_resolution_clock().now();
          done.insert(*n);
          auto nn = mons[*n];
          auto i = std::distance(nn.begin(),
                                 std::ranges::search(nn, mons[g]).begin());
          auto a = mons.getid(nn.first(i));
          auto b = mons.getid(nn.last(nn.size() - i - mons.get_length(g)));
          poly_id agb = poly.multiply_front_and_back(a, lm_to_poly[g], b);
          assert(*n == poly.get_lm_id(agb));
          rows.insert(agb);
          n = todo.erase(n);
          for(auto mm : poly[agb]) {
            if(!done.count(mm))
              todo.insert(mm);
          }
          auto e = std::chrono::high_resolution_clock().now();
          std::chrono::duration<double> el = e - s;
        } else
          n++;
      }
    }

    std::cout << "Saved = " << k << "\n";
    return rows;
  }

  //------------------------------------------------------------------------------
  std::vector<poly_id> find_reducer(mon_id m, bool strategy = false) {

    auto reducers = prefix_trie.divisors(mons[m]);
    if(reducers.empty())
      return std::vector<poly_id>();

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

    std::vector<poly_id> res = { a, match.first, b };

    return res;
  }
  //------------------------------------------------------------------------------

  std::unordered_set<poly_id> symbolic_preprocessing_orig() {
    std::unordered_set<mon_id> todo;
    std::unordered_set<mon_id> done;
    std::unordered_set<poly_id> rows;

    for(const auto& [f, g] : crit_pairs) {
      // add monomials to corresponding sets
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
      poly_id reducer = find_reducer_orig(m);
      if(reducer == 0)
        continue;
      assert(m == poly.get_lm_id(reducer));
      rows.insert(reducer);
      for(auto mm : poly[reducer]) {
        if(!done.count(mm))
          todo.insert(mm);
      }
    }
    return rows;
  }

  //------------------------------------------------------------------------------
  poly_id find_reducer_orig(mon_id m, bool strategy = false) {

    auto& reducers = prefix_trie.divisors(mons[m]);
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
  std::vector<poly_id> compute_new_polynomials(
    std::vector<std::pair<size_t, size_t>>& idxs,
    std::vector<C>& coeffs,
    std::vector<mon_id>& columns) {

    std::vector<poly_id> res;
    std::vector<std::pair<C, mon_id>> p;

    if(idxs.size() == 0)
      return res;

    size_t k = 0;
    size_t cur_i = idxs[0].first;
    for(auto [i, j] : idxs) {
      // a new polynomial starts
      if(i != cur_i) {
        res.push_back(poly.add_polynomial(p));
        p.clear();
        cur_i = i;
      }
      p.emplace_back(coeffs[k++], columns[j]);
    }
    // don't forget to add last element
    res.push_back(poly.add_polynomial(p));

    return res;
  }
  //------------------------------------------------------------------------------
  std::vector<poly_id> reduction(bool interreduce = false) {
    // symbolic preprocessing
    start = std::chrono::high_resolution_clock::now();
    auto rows = symbolic_preprocessing_orig();
    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    sym_pre_time += elapsed.count();

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
    sfmpz_mat_t mat;
    set_up_matrix(mat, rows, columns);

    if(mat->nrow < 20) {
      std::stable_sort(mat->rows, mat->rows + mat->nrow, [](auto a, auto b) {
        auto idx_a = a.indices[0];
        auto idx_b = b.indices[0];
        if(idx_a != idx_b)
          return idx_a > idx_b;
        auto nnz_a = a.nnz;
        auto nnz_b = b.nnz;
        return nnz_a < nnz_b;
      });

      sparse_mat_write(mat, std::cout);
    }

    // reduction
    start = std::chrono::high_resolution_clock::now();
    auto [idxs, entries] = multimodular_gauss_elim(mat, interreduce);
    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    reduction_time += elapsed.count();

    // compute new elements
    start = std::chrono::high_resolution_clock().now();
    std::vector<poly_id> new_elements
      = compute_new_polynomials(idxs, entries, columns);
    end = std::chrono::high_resolution_clock().now();
    elapsed = end - start;
    new_elements_time += elapsed.count();

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
    fmpz_init(denom);
    fmpz_init(tmp);

    i = 0;
    for(auto r : rows) {
      auto row = sparse_mat_row(mat, i++);
      std::span<C> coeffs = poly.get_coefficients(r);

      // compute common denominator so that we can normalize row
      get_common_denom(denom, tmp, coeffs);

      auto p = poly[r];
      auto nnz = p.size();
      sparse_vec_realloc(row, nnz);
      row->nnz = nnz;

      size_t k = 0;
      for(auto it = p.begin(); it != p.end(); it++) {
        auto cc = coeffs[k].data();
        fmpz_set_mpz(tmp, mpq_denref(cc));
        assert(fmpz_divisible(denom, tmp));
        fmpz_divexact(tmp, denom, tmp);
        fmpz_set_mpz(row->entries + k, mpq_numref(cc));
        fmpz_mul(row->entries + k, row->entries + k, tmp);
        row->indices[k] = col_to_id[*it];
        k++;
      }
    }

    fmpz_clear(denom);
    fmpz_clear(tmp);
  }

  //------------------------------------------------------------------------------
  void update_basis_and_amb(std::vector<poly_id>& new_elements) {

    for(const poly_id p_id : new_elements) {
      // only do all of this if we don't termiante next iteration
      if(iter < maxiter) {
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
      }

      // update basis
      basis.push_back(p_id);
    }
  }
};
}
