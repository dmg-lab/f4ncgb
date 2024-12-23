#pragma once

#include <algorithm>
#include <concepts>
#include <cstdint>
#include <limits>
#include <utility>
#include <vector>

#include <boost/multiprecision/gmp.hpp>
#include <unordered_set>

#include "gmp.h"
#include "kommunopp.hpp"
#include "sparse_rref/scalar.h"
#include "sparse_rref/sparse_vec.h"
#include "sparse_rref/sparse_mat.h"
#include "signal_statistics.hpp"

using namespace boost::multiprecision;

namespace kommunopp {


std::vector<ulong> primes = { 5, 7, 23, 17, 31, 37 };

typedef sparse_vec_t<gmp_int> gmp_int_vec_t;
typedef sparse_mat_t<gmp_int> gmp_int_mat_t;

inline ulong
mod_coeff(const gmp_int* c, uint p) {
  return mpz_fdiv_ui(c->data(), p);
}

//------------------------------------------------------------------------------

inline void
vec_mod(snmod_vec_t vec, const gmp_int_vec_t src, uint p) {
  sparse_vec_realloc(vec, src->nnz);
  vec->alloc = src->nnz;
  vec->nnz = 0;
  for(size_t i = 0; i < src->nnz; i++) {
    ulong val = mod_coeff(src->entries + i, p);
    _sparse_vec_set_entry(vec, src->indices[i], &val);
  }
}

//------------------------------------------------------------------------------

inline void
mat_mod(snmod_mat_t mat, const gmp_int_mat_t src, uint p) {
  for(size_t i = 0; i < src->nrow; i++)
    vec_mod(sparse_mat_row(mat, i), sparse_mat_row(src, i), p);
}
//------------------------------------------------------------------------------

using pivots = std::vector<std::pair<slong, slong>>;
inline int
cmp_pivots(pivots x, pivots y) {
  if(x.size() < y.size())
    return -1;
  else if(x.size() > y.size())
    return 1;
  else if(x == y)
    return 0;
  else if(x < y)
    return 1;

  return -1;
}
//------------------------------------------------------------------------------

inline gmp_int
height(gmp_int_mat_t& mat) {
  gmp_int h;
  for(size_t i = 0; i < mat->nrow; i++) {
    auto row = sparse_mat_row(mat, i);
    for(size_t j = 0; j < row->nnz; j++)
      if(mpz_cmpabs((row->entries + j)->data(), h.data()))
        mpz_abs(h.data(), (row->entries + j)->data());
  }
  return h;
}
//------------------------------------------------------------------------------

inline std::vector<mpz_int>
crt_basis(const std::vector<ulong>& moduli, const mpz_int& M) {
  std::vector<mpz_int> res;
  mpz_int m, mm, g, s, t;

  if(moduli.size() == 0)
    return res;
  for(const auto& m_ulong : moduli) {
    m = m_ulong;
    mm = M / m;
    mpz_gcdext(g.backend().data(),
               s.backend().data(),
               t.backend().data(),
               m.backend().data(),
               mm.backend().data());

    assert(g == 1);
    res.push_back((t * mm) % M);
  }
  return res;
}
//------------------------------------------------------------------------------

void
reconstruct_mat(std::vector<gmp_int>& entries,
                std::vector<std::pair<size_t, size_t>>& idxs,
                std::vector<snmod_mat_t*> rrefs,
                std::vector<ulong> primes,
                mpz_int prod) {
  if(rrefs.size() == 0)
    return;

  size_t n_rows = (*rrefs[0])->nrow;
  entries.reserve(sparse_mat_nnz(*rrefs[0]));
  idxs.reserve(sparse_mat_nnz(*rrefs[0]));

  // get all (i,j) where at least one rref is nonzero
  std::map<size_t, std::set<size_t>> nnz_pos;
  for(size_t i = 0; i < n_rows; i++) {
    std::set<size_t> nnz_pos_row;
    for(auto& rref : rrefs) {
      auto row = sparse_mat_row(*rref, i);
      for(size_t k = 0; k < row->nnz; k++)
        nnz_pos_row.insert(row->indices[k]);
    }
    for(auto& j : nnz_pos_row)
      nnz_pos[i].insert(j);
  }

  // actual CRT
  auto crt_base = crt_basis(primes, prod);
  for(size_t i = 0; i < n_rows; i++) {
    for(auto j : nnz_pos[i]) {
      mpz_int tmp = 0;
      for(size_t k = 0; k < rrefs.size(); k++)
        tmp = tmp + crt_base[k] * (*sparse_mat_entry(*rrefs[k], i, j));
      idxs.emplace_back(i, j);
      entries.push_back(tmp.backend());
    }
  }
}
//------------------------------------------------------------------------------

inline bool
verify_result() {
  return true;
}
//------------------------------------------------------------------------------

struct ratrec_data {
  mpz_t mod;
  mpz_t N;
  mpz_t r0;
  mpz_t r1;
  mpz_t t0;
  mpz_t t1;
  mpz_t q;
  mpz_t tmp;
  mpz_t n;
  mpz_t d;

  // Constructor to initialize mpz_t members
  ratrec_data() {
    mpz_init(mod);
    mpz_set_ui(mod, 0);
    mpz_init(N);
    mpz_set_ui(N, 0);
    mpz_init(r0);
    mpz_set_ui(r0, 0);
    mpz_init(r1);
    mpz_set_ui(r1, 0);
    mpz_init(t0);
    mpz_set_ui(t0, 0);
    mpz_init(t1);
    mpz_set_ui(t1, 0);
    mpz_init(q);
    mpz_set_ui(q, 0);
    mpz_init(tmp);
    mpz_set_ui(tmp, 0);
    mpz_init(n);
    mpz_set_ui(n, 0);
    mpz_init(d);
    mpz_set_ui(d, 0);
  }

  // Destructor to clear mpz_t members
  ~ratrec_data() {
    mpz_clear(mod);
    mpz_clear(N);
    mpz_clear(r0);
    mpz_clear(r1);
    mpz_clear(t0);
    mpz_clear(t1);
    mpz_clear(q);
    mpz_clear(tmp);
    mpz_clear(n);
    mpz_clear(d);
  }
};
//------------------------------------------------------------------------------

inline bool
ratrecon(gmp_rational& res, mpz_t u, ratrec_data* data) {

  bool success = false;

  std::cout << "========= Reconstructing entry  ==============" << "\n";
  mpz_int p1(u);
  mpz_int p2(data->mod);
  std::cout << p1 << " % " << p2 << "\n";

  while(mpz_cmp_ui(u, 0) < 0) {
    mpz_add(u, u, data->mod);
  }

  mpz_set(data->r0, data->mod);
  mpz_set_ui(data->t0, 0);

  mpz_set(data->r1, u);
  mpz_set_ui(data->t1, 1);

  while(mpz_cmp(data->r1, data->N) > 0) {
    mpz_fdiv_q(data->q, data->r0, data->r1);

    mpz_mul(data->tmp, data->q, data->r1);
    mpz_sub(data->tmp, data->r0, data->tmp);
    mpz_swap(data->r0, data->r1);
    mpz_swap(data->r1, data->tmp);

    mpz_mul(data->tmp, data->q, data->t1);
    mpz_sub(data->tmp, data->t0, data->tmp);
    mpz_swap(data->t0, data->t1);
    mpz_swap(data->t1, data->tmp);
  }
  mpz_set(data->n, data->r1);
  mpz_set(data->d, data->t1);

  if(mpz_sgn(data->d) < 0) {
    mpz_neg(data->n, data->n);
    mpz_neg(data->d, data->d);
  }
  mpz_gcd(data->q, data->n, data->d);
  if(mpz_cmp(data->d, data->N) <= 0 && mpz_cmp_ui(data->q, 1) == 0) {
    success = true;
    mpz_set(mpq_numref(res.data()), data->n);
    mpz_set(mpq_denref(res.data()), data->d);
    std::cout << "Reconstructed " << mpq_rational(res) << "\n";

  } else
    msg("Rational reconstruction does not exist");

  return success;
}
//------------------------------------------------------------------------------

inline bool
rational_reconstruction(std::vector<gmp_rational>& entries,
                        std::vector<gmp_int>& crt_entries,
                        mpz_int& prod) {
  ratrec_data data;
  mpz_set(data.mod, prod.backend().data());
  // N = floor(sqrt(m/2))
  mpz_fdiv_q_2exp(data.N, data.mod, 1);
  mpz_sqrt(data.N, data.N);

  std::cout << "========= Reconstructing matrix ==============" << "\n";

  entries.reserve(crt_entries.size());
  for(auto& e : crt_entries) {
    gmp_rational r;
    bool success = ratrecon(r, e.data(), &data);
    if(success)
      entries.push_back(r);
    else
      return false;
  }
  return true;
}
//------------------------------------------------------------------------------

std::pair<std::vector<std::pair<size_t, size_t>>, std::vector<gmp_rational>>
multimodular_rref(gmp_int_mat_t& mat, bool proof = true) {
  field_t F;
  rref_option_t opt;
  opt->verbose = true;
  opt->is_back_sub = false;
  opt->print_step = 100;
  opt->pivot_dir = false;
  opt->search_depth = INT_MAX;

  // TODO : adapt
  BS::thread_pool pool(2);

  std::vector<snmod_mat_t*> rrefs;
  pivots best_piv;
  std::vector<pivots> pivs;
  std::vector<ulong> used_primes;
  std::vector<snmod_mat_t*> good_rrefs;
  std::vector<pivots> good_pivs;
  std::vector<ulong> good_primes;

  std::vector<gmp_rational> rat_entries;
  std::vector<std::pair<size_t, size_t>> idxs;

  size_t i = 0;

  mpz_int h = mpz_int(height(mat));
  mpz_int prod = 1;
  mpz_int M = mat->ncol * 100000 * (h + 100) * h + 1;

  // TODO : remove !!!!
  M = 25;

  size_t MAX_PRIMES = primes.size();
  ulong p;

  while(true) {
    while(prod < M) {
      if(i >= MAX_PRIMES)
        die(-1, "Multimodular rref is not converging");
      p = primes[i++];

      msg("Computing mod %d", p);

      field_init(F, FIELD_Fp, std::vector<ulong>{ p });

      snmod_mat_t nmod_mat;
      sparse_mat_init(nmod_mat, mat->nrow, mat->ncol);
      mat_mod(nmod_mat, mat, p);
      pivots piv = sparse_mat_rref(nmod_mat, F, pool, opt);

      sparse_mat_write(nmod_mat, std::cout);

      if(cmp_pivots(best_piv, piv) <= 0) {
        best_piv = piv;
        pivs.push_back(piv);
        rrefs.push_back(&nmod_mat);
        used_primes.push_back(p);
        prod = prod * p;
      } else {
        msg("Excluding prime %d (bad pivots)", p);
        sparse_mat_clear(nmod_mat);
      }
    }
    prod = 1;
    good_primes.clear();
    good_rrefs.clear();
    good_pivs.clear();
    for(size_t r = 0; r < rrefs.size(); r++) {
      if(cmp_pivots(best_piv, pivs[r]) <= 0) {
        good_primes.push_back(used_primes[r]);
        good_rrefs.push_back(rrefs[r]);
        good_pivs.push_back(pivs[r]);
        prod = prod * used_primes[r];
      }
    }

    std::vector<gmp_int> crt_entries;
    reconstruct_mat(crt_entries, idxs, good_rrefs, good_primes, prod);
    bool success = rational_reconstruction(rat_entries, crt_entries, prod);

    if(!success) {
      msg("Reconstruction unsuccessfull. Increasing bound.");
      M = prod * p * p * p;
    }
    if(!proof or verify_result())
      break;
  }
  for(const auto& rref : rrefs)
    sparse_mat_clear(*rref);

  return std::make_pair(idxs, rat_entries);
}
}
