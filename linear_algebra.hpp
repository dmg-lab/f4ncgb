#pragma once

#include <algorithm>
#include <chrono>
#include <concepts>
#include <cstdint>
#include <limits>
#include <utility>
#include <vector>

#include <boost/multiprecision/gmp.hpp>
#include <unordered_set>

#include "gmp.h"
#include "kommunopp.hpp"
#include "signal_statistics.hpp"
#include "sparse_rref/scalar.h"
#include "sparse_rref/sparse_mat.h"
#include "sparse_rref/sparse_vec.h"

extern double crt_time, ratrec_time, rref_time;

using namespace boost::multiprecision;

namespace kommunopp {

/* std::vector<ulong> primes = { 5, 7, 23, 17, 31, 37 }; */

std::vector<ulong> primes
  = { 2147483629, 2147483587, 2147483579, 2147483563, 2147483549, 2147483543,
      2147483497, 2147483489, 2147483477, 2147483423, 2147483399, 2147483353,
      2147483323, 2147483269, 2147483249, 2147483237, 2147483179, 2147483171,
      2147483137, 2147483123, 2147483077, 2147483069, 2147483059, 2147483053,
      2147483033, 2147483029, 2147482951, 2147482949, 2147482943, 2147482937,
      2147482921 };

typedef sparse_vec_t<fmpz> sfmpz_vec_t;
typedef sparse_mat_t<fmpz> sfmpz_mat_t;

//------------------------------------------------------------------------------

inline void
vec_mod(snmod_vec_t vec, const sfmpz_vec_t src, nmod_t mod) {
  sparse_vec_realloc(vec, src->nnz);
  vec->alloc = src->nnz;
  vec->nnz = 0;
  for(size_t i = 0; i < src->nnz; i++) {
    ulong val = fmpz_get_nmod(src->entries + i, mod);
    _sparse_vec_set_entry(vec, src->indices[i], &val);
  }
}

//------------------------------------------------------------------------------

inline void
mat_mod(snmod_mat_t mat, const sfmpz_mat_t src, ulong p) {
  nmod_t mod;
  nmod_init(&mod, p);
  for(size_t i = 0; i < src->nrow; i++)
    vec_mod(sparse_mat_row(mat, i), sparse_mat_row(src, i), mod);
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
height(sfmpz_mat_t mat) {
  gmp_int h;
  mpz_t tmp;
  mpz_init(tmp);
  for(size_t i = 0; i < mat->nrow; i++) {
    auto row = sparse_mat_row(mat, i);
    for(size_t j = 0; j < row->nnz; j++) {
      fmpz_get_mpz(tmp, row->entries + j);
      if(mpz_cmpabs(tmp, h.data()))
        mpz_abs(h.data(), tmp);
    }
  }
  mpz_clear(tmp);

  return h;
}
//------------------------------------------------------------------------------

void
reconstruct_mat(fmpz*& entries,
                std::vector<std::pair<size_t, size_t>>& idxs,
                std::vector<snmod_mat_t*> rrefs,
                std::vector<ulong> primes,
                mpz_int prod) {
  if(rrefs.size() == 0)
    return;

  size_t n_rows = (*rrefs[0])->nrow;

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
  size_t nnz = 0;
  for(const auto& [key, values] : nnz_pos)
    nnz += values.size();

  size_t len = primes.size();
  idxs.reserve(sparse_mat_nnz(*rrefs[0]));
  fmpz* moduli = new fmpz[len];
  for(size_t i = 0; i < len; ++i)
    fmpz_init_set_ui(moduli + i, primes[i]);
  entries = new fmpz[nnz];
  for(size_t i = 0; i < nnz; i++)
    fmpz_init(entries + i);

  fmpz_multi_CRT_t crt_base;
  fmpz_multi_CRT_init(crt_base);
  int res = fmpz_multi_CRT_precompute(crt_base, moduli, len);
  if(!res)
    die(12, "Problem with CRT");

  fmpz* inputs = new fmpz[rrefs.size()];
  for(size_t i = 0; i < rrefs.size(); i++)
    fmpz_init(inputs + i);

  size_t idx = 0;
  for(size_t i = 0; i < n_rows; i++) {
    for(auto j : nnz_pos[i]) {
      idxs.emplace_back(i, j);
      for(size_t k = 0; k < rrefs.size(); k++)
        fmpz_set_ui(inputs + k, *sparse_mat_entry(*rrefs[k], i, j));
      fmpz_multi_CRT_precomp(entries + idx++, crt_base, inputs, 0);
    }
  }

  fmpz_multi_CRT_clear(crt_base);
  for(size_t i = 0; i < primes.size(); i++)
    fmpz_clear(moduli + i);
  for(size_t i = 0; i < rrefs.size(); i++)
    fmpz_clear(inputs + i);
  free(moduli);
  free(inputs);
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

  mpz_int p1(u);
  mpz_int p2(data->mod);
  // std::cout << p1 << " % " << p2 << "\n";

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

  } else
    msg("Rational reconstruction does not exist");

  return success;
}
//------------------------------------------------------------------------------

inline bool
rational_reconstruction(std::vector<gmp_rational>& entries,
                        fmpz * crt_entries,
                        mpz_int& prod,
                        size_t N) {
  ratrec_data data;
  mpz_set(data.mod, prod.backend().data());
  // N = floor(sqrt(m/2))
  mpz_fdiv_q_2exp(data.N, data.mod, 1);
  mpz_sqrt(data.N, data.N);

  bool res = true;

  std::cout << "========= Reconstructing matrix ==============" << "\n";
  entries.reserve(N);
  mpz_t u;
  mpz_init(u);
  for(size_t i = 0; i < N; i++) {
    gmp_rational r;
    fmpz_get_mpz(u, crt_entries + i);
    bool success = ratrecon(r, u, &data);
    if(success)
      entries.push_back(r);
    else {
      res = false;
      break;
    }
  }
  mpz_clear(u);
  return res;
}
//------------------------------------------------------------------------------
template<typename T>
std::vector<std::pair<slong, slong>>
my_sparse_mat_rref(sparse_mat_t<T> mat,
                   field_t F,
                   BS::thread_pool& pool,
                   rref_option_t opt) {
  // first canonicalize, sort and compress the matrix
  sparse_mat_compress(mat);

  T scalar[1];
  scalar_init(scalar);

  slong r;
  slong start_row = 0;
  slong min_row;
  slong min;

  std::vector<std::pair<slong, slong>> pivots;

  for(slong c = 0; c < mat->ncol; c++) {

    min_row = -1;
    min = mat->ncol + 1;
    for(r = start_row; r < mat->nrow; r++) {
      auto therow = sparse_mat_row(mat, r);
      if(therow->nnz > 0 and therow->nnz < min)
        if(therow->indices[0] == c) {
          min_row = r;
          min = therow->nnz;
        }
    }
    if(min_row == -1)
      continue;

    // will use row r to reduce column c
    r = min_row;

    // rescale row
    scalar_inv(scalar, sparse_mat_entry(mat, r, c, true), F);
    sparse_vec_rescale(sparse_mat_row(mat, r), scalar, F);

    // swap row to top
    std::swap(mat->rows[start_row], mat->rows[r]);
    pivots.emplace_back(start_row, c);

    // eliminate
    auto therow = sparse_mat_row(mat, start_row);
    for(size_t i = 0; i < mat->nrow; i++) {
      if(i == start_row)
        continue;
      auto row_i = sparse_mat_row(mat, i);
      auto b = sparse_vec_entry(row_i, c);
      if(b != NULL) {
        sparse_vec_sub_mul(sparse_mat_row(mat, i), therow, b, F);
      }
    }
    start_row++;
  }

  return pivots;
}
//------------------------------------------------------------------------------

std::pair<std::vector<std::pair<size_t, size_t>>, std::vector<gmp_rational>>
multimodular_rref(sfmpz_mat_t& mat, bool proof = true) {
  field_t F;
  rref_option_t opt;
  opt->verbose = true;
  opt->is_back_sub = true;
  opt->print_step = 100;
  opt->pivot_dir = true;
  opt->search_depth = INT_MAX;

  // TODO : adapt
  BS::thread_pool pool(1);

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

  size_t MAX_PRIMES = primes.size();
  ulong p;

  while(true) {
    while(prod < M) {
      if(i >= MAX_PRIMES)
        die(-1, "Multimodular rref is not converging");
      p = primes[i++];

      msg("Computing mod %ul", p);

      field_init(F, FIELD_Fp, std::vector<ulong>{ p });

      snmod_mat_t* nmod_mat = new snmod_mat_t[1];
      sparse_mat_init(*nmod_mat, mat->nrow, mat->ncol);
      mat_mod(*nmod_mat, mat, p);

      auto start = std::chrono::high_resolution_clock().now();
      pivots piv = my_sparse_mat_rref(*nmod_mat, F, pool, opt);
      auto end = std::chrono::high_resolution_clock().now();
      std::chrono::duration<double> elapsed = end - start;
      rref_time += elapsed.count();

      std::cout << "Computed rref\n";

      if(cmp_pivots(best_piv, piv) <= 0) {
        best_piv = piv;
        pivs.push_back(piv);
        rrefs.push_back(nmod_mat);
        used_primes.push_back(p);
        prod = prod * p;
      } else {
        msg("Excluding prime %d (bad pivots)", p);
        sparse_mat_clear(*nmod_mat);
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

    // Initialize crt_entries in reconstruction
    // and clear here
    fmpz* crt_entries;
    auto start = std::chrono::high_resolution_clock().now();
    reconstruct_mat(crt_entries, idxs, good_rrefs, good_primes, prod);
    auto end = std::chrono::high_resolution_clock().now();
    std::chrono::duration<double> elapsed = end - start;
    crt_time += elapsed.count();

    start = std::chrono::high_resolution_clock().now();
    bool success = rational_reconstruction(rat_entries, crt_entries, prod, idxs.size());
    for(size_t i = 0; i < idxs.size(); i++)
      fmpz_clear(crt_entries + i);
    free(crt_entries);
    end = std::chrono::high_resolution_clock().now();
    elapsed = end - start;
    ratrec_time += elapsed.count();

    if(!success) {
      msg("Reconstruction unsuccessfull. Increasing bound.");
      M = prod * p * p;
      continue;
    }

    if(!proof or verify_result())
      break;
  }
  for(const auto& rref : rrefs)
    sparse_mat_clear(*rref);

  return std::make_pair(idxs, rat_entries);
}

}
