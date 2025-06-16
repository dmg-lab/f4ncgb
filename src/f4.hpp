#pragma once

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <ostream>
#include <set>
#include <span>
#include <utility>
#include <vector>

#include <boost/align/align_down.hpp>
#include <boost/align/align_up.hpp>
#include <boost/container/small_vector.hpp>
#include <boost/multiprecision/gmp.hpp>
#include <boost/unordered/unordered_map.hpp>
#include <boost/unordered/unordered_set.hpp>

#include "ambiguity.hpp"
#include "f4ncgb.hpp"
#include "gmp.h"
#include "linear_algebra.hpp"
#include "parser.hpp"
#include "profiling.hpp"
#include "signal_statistics.hpp"
#include "sparse_rref/sparse_mat.h"
#include "sparse_rref/sparse_vec.h"
#include "sparse_rref/thread_pool.hpp"

#include "monomial_trie.hpp"

extern template struct f4ncgb::monomial_trie<uint8_t, uint32_t>;
extern int verbose;
extern int proof;

using namespace boost::multiprecision;

namespace f4ncgb {

struct metadata_monomial {
  uint32_t length = 0;
};

struct metadata_polynomial {
  uint32_t length = 0;
};

//------------------------------------------------------------------------------

template<size_t Nblocks,
         internal::value_concept V = uint8_t,
         typename I = uint32_t,
         typename C = boost::multiprecision::gmp_rational>
struct f4 {
  using MM = metadata_monomial;
  using PM = metadata_polynomial;
  using coefficient = C;
  using monomial_store = internal::monomial_store<MM, V, I, Nblocks>;
  using polynomial_store = internal::polynomial_store<PM, MM, V, I, C, Nblocks>;

  using monomial = std::span<const V>;
  using mon_id = I;
  using poly_id = I;

  using ambiguity_ = ambiguity<mon_id>;
  using amb_hash = ambiguity_hash<mon_id>;
  using crit_pair = std::pair<poly_id, poly_id>;

  using monomial_trie_ = monomial_trie<V, I>;

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
  parser_context& context;

  monomial_store mons;
  polynomial_store poly;
  std::vector<poly_id> basis;
  std::map<size_t, boost::unordered_set<ambiguity_, amb_hash>> amb;
  std::set<crit_pair> crit_pairs;
  boost::unordered_map<mon_id, poly_id> lm_to_poly;
  monomial_trie_ prefix_trie;
  monomial_trie_ suffix_trie;
  std::unique_ptr<BS::thread_pool<BS::none>> pool;

  size_t characteristic = 0;
  size_t iter = 0;
  size_t maxiter = UINT_MAX;
  size_t maxdeg = UINT_MAX;
  static constexpr bool block_order = Nblocks > 0;

  std::ofstream proof_file;

  f4(parser_context& context_,
     size_t nvars,
     size_t characteristic_,
     size_t maxiter_,
     size_t maxdeg_,
     size_t num_threads,
     const std::string& proof_file_)
    : context(context_)
    , mons()
    , poly(mons)
    , prefix_trie(nvars)
    , suffix_trie(nvars)
    , characteristic(characteristic_)
    , maxiter(maxiter_)
    , maxdeg(maxdeg_) {

    mons.set_blocks(context.block_sizes());

    if(num_threads > 1)
      pool = std::make_unique<BS::thread_pool<BS::none>>(num_threads);

    if(proof_file_ != "") {
      proof_file.open(proof_file_, std::ios_base::trunc);
      if(!proof_file)
        die(19, "Failed to open proof file.");
    }
  }
  //------------------------------------------------------------------------------
  struct triplet {
    mon_id a;
    int i;
    mon_id b;

    triplet(mon_id a_, int i_, mon_id b_)
      : a(a_)
      , i(i_)
      , b(b_) {}
  };

  struct cofactor {
    C c;
    triplet t;

    cofactor(C c_, mon_id a_, int i_, mon_id b_)
      : c(c_)
      , t(a_, i_, b_) {}

    inline mon_id a() { return t.a; }
    inline int i() { return t.i; }
    inline mon_id b() { return t.b; }

    inline void log(std::ostream& o,
                    parser_context& context,
                    polynomial_store& poly,
                    bool first = false) {
      if(!first)
        o << (mpq_sgn(c.data()) > 0 ? " + " : " ");
      o << c;
      if(t.a != 0)
        context.to_msolve_mon(o, poly, t.a, false);
      if(t.i >= 0)
        o << "*[" << t.i << "]";
      else
        o << "*[i" << -t.i - 1 << "]";
      if(t.b != 0)
        context.to_msolve_mon(o, poly, t.b, false);
    }
  };

  std::vector<triplet> extended_rows;
  std::vector<std::vector<cofactor>> cofactors;
  //------------------------------------------------------------------------------
  inline parse_res read_input(parser_context& context) {
    parse_res res = parse_rest_into_polynomial_store(context, poly);
    // Print in case of doubt.
    // context.to_msolve(std::cout, poly);
    return res;
  }
  //------------------------------------------------------------------------------
  inline void write_basis(std::ostream& o) {
    bool first = true;
    for(size_t n = 1; n < basis.size(); n++) {
      context.to_msolve_poly(o, poly, basis[n], first);
      first = false;
    }
    o << "\n";
  }
  //------------------------------------------------------------------------------
  void interreduce_and_add_to_basis(std::vector<poly_id> input) {

    if(verbose > 1)
      msg("Linearly interreducing input of size %d.", input.size());

    size_t i = 0;
    for(const auto& p : input) {
      extended_rows.emplace_back(0, i, 0);
      extended_rows.emplace_back(0, i++, 0);
      crit_pair c(p, p);
      crit_pairs.insert(c);
    }

    {
      F4NCGB_TIME(crit_pair);
      stage_crit_pairs();
    }

    // to leave 0th position open; just like in basis
    cofactors.emplace_back();

    std::vector<poly_id> new_elements = reduction(true);
    if(verbose > 1)
      msg("Adding %d input elements to basis.", new_elements.size());

    update_basis_and_amb(new_elements);
  }
  //------------------------------------------------------------------------------

  std::vector<poly_id> input;
  void compute_basis() {
    // add something at 0th position
    // so that index 0 remains free
    basis.push_back(0);
    input.clear();
    input.resize(
      static_cast<size_t>(std::distance(poly.begin() + 1, poly.end())));
    std::copy(poly.begin() + 1, poly.end(), input.begin());

    // add input to critical pairs
    interreduce_and_add_to_basis(input);

    // main loop
    iter = 0;
    while((!amb.empty() or !crit_pairs.empty()) and iter < maxiter) {
      {
        F4NCGB_TIME(crit_pair);
        stage_crit_pairs();
      }

      if(verbose > 1)
        msg("Reducing %d critical pairs.", crit_pairs.size());
      std::vector<poly_id> new_elements = reduction();

      if(verbose > 1)
        msg("Adding %d new elements to basis.", new_elements.size());

      update_basis_and_amb(new_elements);

      iter++;
      if(verbose > 0)
        msg("==== Iteration %d has finished. Basis has now %d elements ====",
            iter,
            basis.size() - 1);
    }
  }
  //------------------------------------------------------------------------------
  inline crit_pair to_crit_pair(const ambiguity_& a) {
    poly_id i = lm_to_poly[a.i()];
    poly_id j = lm_to_poly[a.j()];

    mon_id ai = a.ai();
    mon_id ci = a.ci();
    mon_id aj = a.aj();
    mon_id cj = a.cj();

    poly_id f = poly.multiply_front_and_back(ai, i, ci);
    poly_id g = poly.multiply_front_and_back(aj, j, cj);

    if(proof > 0) {
      assert(poly.get_idx(i) > 0);
      assert(poly.get_idx(j) > 0);
      extended_rows.emplace_back(ai, poly.get_idx(i), ci);
      extended_rows.emplace_back(aj, poly.get_idx(j), cj);
    }
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
  int inline divisible_by(const ambiguity_& self,
                          const ambiguity_& other) const {
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
  std::vector<ambiguity_> tmp;
  std::vector<char> to_remove;
  void gebauer_moeller(boost::unordered_set<ambiguity_, amb_hash>& new_amb) {
    // first index is always the newer polynomial

    auto cmp = [this](ambiguity_& a, ambiguity_& b) {
      if(a.degree() != b.degree())
        return a.degree() < b.degree();
      if(a.j() != b.j())
        return a.j() < b.j();
      return this->mons.template cmp<block_order>(a.aj(), b.aj());
    };

    tmp.clear();
    tmp.reserve(
      static_cast<size_t>(std::distance(new_amb.begin(), new_amb.end())));
    std::copy(new_amb.begin(), new_amb.end(), std::back_inserter(tmp));
    std::sort(tmp.begin(), tmp.end(), cmp);

    to_remove.clear();
    to_remove.resize(tmp.size());
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
  }
  //------------------------------------------------------------------------------
  std::vector<std::pair<mon_id, size_t>> overlaps;
  std::vector<std::pair<mon_id, size_t>> inclusions;
  boost::unordered_set<ambiguity_, amb_hash> new_amb;
  void compute_ambiguities(mon_id i) {
    monomial m = mons[i];
    monomial ab, b, bc, abc;

    new_amb.clear();
    overlaps.clear();
    inclusions.clear();

    {
      F4NCGB_TIME(overlap);
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
        ambiguity_ a(d, i, j, 0, ci, aj, 0);
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
        ambiguity_ a(d, i, j, ai, 0, 0, cj);
        new_amb.insert(a);
      }
    }

    {
      F4NCGB_TIME(inclusion);
      // inclusions with m = ABC
      I d = m.size();
      // k determines where B starts in m = ABC
      for(auto [j, k] : inclusions) {
        if(i == j)
          continue;
        b = mons[j];
        I aj = mons.getid(m.first(k));
        I cj = mons.getid(m.last(d - k - b.size()));
        ambiguity_ a(d, i, j, 0, 0, aj, cj);
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
        ambiguity_ a(d, i, j, ai, ci, 0, 0);
        new_amb.insert(a);
      }
    }

    gebauer_moeller(new_amb);

    for(const auto& a : new_amb)
      amb[a.degree()].insert(a);
  }

  //------------------------------------------------------------------------------
  std::vector<poly_id> symbolic_preprocessing() {
    F4NCGB_TIME(sym_pre);
    boost::unordered_set<mon_id> todo;
    boost::unordered_set<mon_id> done;
    std::vector<poly_id> rows;

    for(const auto& [f, g] : crit_pairs) {
      // add monomials to corresponding sets
      auto mon_it = poly[f];
      done.insert(*mon_it.begin());
      todo.insert(++mon_it.begin(), mon_it.end());

      mon_it = poly[g];
      done.insert(*mon_it.begin());
      todo.insert(++mon_it.begin(), mon_it.end());

      rows.push_back(f);
      rows.push_back(g);
    }

    while(!todo.empty()) {
      mon_id m = *todo.begin();
      todo.erase(todo.begin());
      done.insert(m);

      poly_id reducer = find_reducer(m);
      if(reducer == 0)
        continue;
      assert(m == poly.get_lm_id(reducer));
      rows.push_back(reducer);
      for(auto mm : poly[reducer]) {
        if(!done.count(mm))
          todo.insert(mm);
      }
    }
    return rows;
  }

  //------------------------------------------------------------------------------
  poly_id find_reducer(mon_id m) {
    auto& reducers = prefix_trie.divisors(mons[m]);
    if(reducers.empty())
      return 0;

    // strategy  : the one with smallest lm
    std::pair<mon_id, size_t> match = *std::max_element(
      reducers.begin(), reducers.end(), [this](auto a, auto b) {
        return this->mons.template cmp<block_order>(b.first, a.first);
      });

    monomial mm = mons[m];
    monomial lm = mons[match.first];
    mon_id a = mons.getid(mm.first(match.second));
    mon_id b = mons.getid(mm.last(mm.size() - match.second - lm.size()));
    poly_id g = lm_to_poly[match.first];

    if(proof > 0) {
      assert(poly.get_idx(g) > 0);
      extended_rows.emplace_back(a, poly.get_idx(g), b);
    }

    poly_id res = poly.multiply_front_and_back(a, g, b);

    return res;
  }

  //------------------------------------------------------------------------------
  void log_cofactors() {

    // compute expanded proofs
    if(proof > 1 and basis.size() > 1) {
      std::vector<cofactor> expanded;
      for(size_t n = basis.size(); n < cofactors.size(); n++) {
        expanded.clear();
        for(auto& cofactor : cofactors[n]) {
          C& c = cofactor.c;
          mon_id a = cofactor.a();
          mon_id b = cofactor.b();
          for(auto& cofactor_i : cofactors[(size_t)cofactor.i()]) {
            C cc;
            mpq_mul(cc.data(), c.data(), cofactor_i.c.data());
            mon_id aa = mons.get_product_id(a, cofactor_i.a());
            mon_id bb = mons.get_product_id(cofactor_i.b(), b);
            expanded.emplace_back(cc, aa, cofactor_i.i(), bb);
          }
        }
        cofactors[n] = std::move(expanded);
      }
      // non-expanded proof
    } else if(proof == 1) {
      // mark input
      if(basis.size() == 1)
        for(size_t n = basis.size(); n < cofactors.size(); n++)
          for(auto& cofactor : cofactors[n])
            cofactor.t.i = -cofactor.t.i - 1;
      // not input -> shift all indices down by one
      else
        for(size_t n = basis.size(); n < cofactors.size(); n++)
          for(auto& cofactor : cofactors[n])
            cofactor.t.i--;
    }

    // write to file
    bool first = true;
    for(size_t n = basis.size(); n < cofactors.size(); n++) {
      first = true;
      for(auto& cofactor : cofactors[n]) {
        cofactor.log(proof_file, context, poly, first);
        first = false;
      }
      proof_file << std::endl;
    }
  }
  //------------------------------------------------------------------------------
  std::vector<poly_id> res;
  std::vector<std::pair<C, mon_id>> p;
  std::vector<cofactor> current_cofactors;
  const std::vector<poly_id>& compute_new_polynomials(
    std::vector<std::pair<size_t, size_t>>& idxs,
    std::vector<C>& coeffs,
    std::vector<mon_id>& columns) {

    res.clear();
    p.clear();
    current_cofactors.clear();

    poly.reset();

    if(idxs.size() == 0)
      return res;

    size_t n = columns.size();
    size_t k = 0;
    size_t cur_i = idxs[0].first;
    for(auto [i, j] : idxs) {
      // a new polynomial starts
      if(i != cur_i) {
        res.push_back(poly.add_polynomial(p));
        p.clear();
        cur_i = i;
        if(proof > 0) {
          cofactors.emplace_back(std::move(current_cofactors));
          current_cofactors.clear();
        }
      }
      if(j < n) {
        assert(mpq_rational(coeffs[k]) != 0);
        p.emplace_back(coeffs[k], columns[j]);
      } else {
        current_cofactors.emplace_back(coeffs[k],
                                       extended_rows[j - n].a,
                                       extended_rows[j - n].i,
                                       extended_rows[j - n].b);
      }
      k++;
    }
    // don't forget to add last element
    res.push_back(poly.add_polynomial(p));
    if(proof > 0)
      cofactors.emplace_back(std::move(current_cofactors));

    extended_rows.clear();

    return res;
  }
  //------------------------------------------------------------------------------
  boost::unordered_set<mon_id> col_set;
  std::vector<mon_id> columns;
  const std::vector<poly_id>& reduction(bool interreduce = false) {
    // symbolic preprocessing
    auto rows = symbolic_preprocessing();
    crit_pairs.clear();

    col_set.clear();

    // make columns
    // columns are sorted in DESCENDING order
    for(const auto r : rows) {
      auto p = poly[r];
      col_set.insert(p.begin(), p.end());
    }
    columns.clear();
    columns.resize(
      static_cast<size_t>(std::distance(col_set.begin(), col_set.end())));
    std::move(col_set.begin(), col_set.end(), columns.begin());
    auto cmp = [this](const mon_id a, const mon_id b) {
      return this->mons.template cmp<block_order>(b, a);
    };
    std::sort(columns.begin(), columns.end(), cmp);

    // set up matrix
    sfmpz_mat_t mat;
    set_up_matrix(mat, rows, columns);

    // reduction
    F4NCGB_PROFILE(auto timer = gstats.time(gstats.reduction));
    auto [idxs, entries]
      = linear_algebra(mat, characteristic, pool, interreduce);

    // compute new elements
    F4NCGB_PROFILE(auto timer2 = gstats.time(gstats.new_elements));
    auto& new_elements = compute_new_polynomials(idxs, entries, columns);

    sparse_mat_clear(mat);

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
  boost::unordered_map<mon_id, size_t> col_to_id;
  void set_up_matrix(sfmpz_mat_t mat,
                     std::vector<poly_id>& rows,
                     std::vector<mon_id>& columns) {

    col_to_id.clear();
    size_t i = 0;
    for(auto c : columns)
      col_to_id[c] = i++;

    size_t m = rows.size();
    size_t n = columns.size();

    // initialize matrix
    if(proof > 0)
      sparse_mat_init(mat, m, m + n);
    else
      sparse_mat_init(mat, m, n);

    if(verbose > 2)
      msg("Setting up matrix of size (%d, %d)", mat->nrow, mat->ncol);

    // set all entries
    fmpz_t denom;
    fmpz_t tmp;
    fmpz_init(denom);
    fmpz_init(tmp);

    i = 0;
    for(auto r : rows) {
      auto row = sparse_mat_row(mat, i);
      std::span<C> coeffs = poly.get_coefficients(r);
      // compute common denominator so that we can normalize row
      get_common_denom(denom, tmp, coeffs);

      auto p = poly[r];
      auto nnz = p.size();
      if(proof > 0)
        nnz += 1;// for transformation matrix
      sparse_vec_realloc(row, nnz);
      row->nnz = nnz;

      // insert poly
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

      // insert transformation matrix - if required
      if(proof > 0) {
        row->indices[k] = n + i;
        fmpz_set(row->entries + k, denom);
      }
      i++;
    }

    fmpz_clear(denom);
    fmpz_clear(tmp);
  }
  //------------------------------------------------------------------------------
  void update_basis_and_amb(std::vector<poly_id>& new_elements) {

    // log cofactors
    if(proof > 0) {
      F4NCGB_TIME(other);
      log_cofactors();
    }

    uint16_t n = basis.size();
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
        F4NCGB_TIME(amb);
        compute_ambiguities(m_id);
      }

      // update basis
      poly.set_idx(p_id, n++);
      basis.push_back(p_id);
    }
  }
};
}
