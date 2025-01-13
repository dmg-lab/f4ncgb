#pragma once

#include <algorithm>
#include <cstdint>
#include <flint/nmod.h>
#include <map>
#include <set>
#include <sys/errno.h>
#include <utility>
#include <vector>

#include <boost/multiprecision/gmp.hpp>
#include <boost/unordered_set.hpp>

#include "gmp.h"
#include "signal_statistics.hpp"

#include "sparse_rref/sparse_mat.h"
#include "sparse_rref/sparse_vec.h"

#include "primes.hpp"
#include "profiling.hpp"

using namespace boost::multiprecision;

namespace kommunopp {

extern std::vector<uint32_t> PRIMES;
static size_t const MAX_PRIMES = PRIMES.size();

typedef sparse_vec_t<fmpz> sfmpz_vec_t;
typedef sparse_mat_t<fmpz> sfmpz_mat_t;
typedef sparse_vec_t<uint32_t> uint32_vec_t;
typedef sparse_mat_t<uint32_t> uint32_mat_t;

//------------------------------------------------------------------------------

inline bool
vec_mod(uint32_vec_t vec, const sfmpz_vec_t src, nmod_t mod) {
  assert(vec->nnz <= src->nnz);
  auto nnz = src->nnz;
  sparse_vec_realloc(vec, nnz);
  vec->nnz = nnz;
  std::copy(src->indices, src->indices + nnz, vec->indices);
  for(size_t i = 0; i < nnz; i++) {
    uint32_t val = fmpz_get_nmod(src->entries + i, mod);
    if(i == 0 and val == 0)
      return false;
    vec->entries[i] = val;
  }
  return true;
}

//------------------------------------------------------------------------------

inline bool
mat_mod(uint32_mat_t mat, const sfmpz_mat_t src, nmod_t mod, bool* trace) {
  bool res = true;
  for(size_t i = 0; i < src->nrow; i++) {
    if(!trace[i]) {
      res = vec_mod(sparse_mat_row(mat, i), sparse_mat_row(src, i), mod);
      if(!res)
        break;
    }
  }
  return res;
}
//------------------------------------------------------------------------------

using pivots = std::vector<size_t>;
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
crt_reconstruction(fmpz*& entries,
                   std::vector<std::pair<size_t, size_t>>& idxs,
                   std::vector<sparse_mat_struct<uint32_t>*>& rrefs,
                   std::vector<ulong>& primes,
                   std::vector<ulong>& relevant_rows) {
  if(rrefs.size() == 0)
    return;

  // get all (i,j) where at least one rref is nonzero
  std::map<size_t, std::set<size_t>> nnz_pos;
  for(size_t i : relevant_rows) {
    auto& nnz_pos_row = nnz_pos[i];
    for(auto& rref : rrefs) {
      auto row = sparse_mat_row(rref, i);
      nnz_pos_row.insert(row->indices, row->indices + row->nnz);
    }
  }
  size_t nnz = 0;
  for(const auto& [key, values] : nnz_pos)
    nnz += values.size();
  idxs.reserve(nnz);

  size_t len = primes.size();
  fmpz* moduli = new fmpz[len];
  entries = new fmpz[nnz];
  for(size_t i = 0; i < len; i++)
    fmpz_init_set_ui(moduli + i, primes[i]);
  for(size_t i = 0; i < nnz; i++)
    fmpz_init(entries + i);

  fmpz_multi_CRT_t crt_base;
  fmpz_multi_CRT_init(crt_base);
  int res = fmpz_multi_CRT_precompute(crt_base, moduli, static_cast<long>(len));
  if(!res)
    die(12, "Problem with CRT");

  fmpz* inputs = new fmpz[rrefs.size()];
  for(size_t i = 0; i < rrefs.size(); i++)
    fmpz_init(inputs + i);

  size_t idx = 0;
  for(size_t i : relevant_rows) {
    for(auto j : nnz_pos[i]) {
      idxs.emplace_back(i, j);
      for(size_t k = 0; k < rrefs.size(); k++) {
        auto c = sparse_mat_entry(rrefs[k], i, j);
        if(c == nullptr)
          fmpz_set_ui(inputs + k, 0);
        else
          fmpz_set_ui(inputs + k, *c);
      }
      fmpz_multi_CRT_precomp(entries + idx++, crt_base, inputs, 0);
    }
  }

  fmpz_multi_CRT_clear(crt_base);
  for(size_t i = 0; i < primes.size(); i++)
    fmpz_clear(moduli + i);
  for(size_t i = 0; i < rrefs.size(); i++)
    fmpz_clear(inputs + i);
  delete[] moduli;
  delete[] inputs;
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
                        fmpz* crt_entries,
                        mpz_int& prod,
                        size_t N) {
  ratrec_data data;
  mpz_set(data.mod, prod.backend().data());
  // N = floor(sqrt(m/2))
  mpz_fdiv_q_2exp(data.N, data.mod, 1);
  mpz_sqrt(data.N, data.N);

  entries.reserve(N);

  bool res = true;
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
void inline copy_to_buffer(T* buffer, uint32_vec_t vec) {
  for(size_t i = 0; i < vec->nnz; i++)
    buffer[vec->indices[i]] = vec->entries[i];
}
//------------------------------------------------------------------------------
void inline copy_from_buffer_and_clear(int64_t* buffer,
                                       std::vector<size_t>& buffer_ids,
                                       uint32_vec_t vec) {
  size_t nnz = buffer_ids.size();

  sparse_vec_clear(vec);
  sparse_vec_realloc(vec, nnz);
  vec->nnz = nnz;

  size_t j = 0;
  for(size_t i : buffer_ids) {
    vec->indices[j] = i;
    vec->entries[j] = (uint32_t)buffer[i];
    j++;
    buffer[i] = 0;
  }
}
//------------------------------------------------------------------------------
void inline normalize_row(uint32_vec_t vec, nmod_t mod) {
  uint32_t c = vec->entries[0];
  uint32_t inv = nmod_inv(c, mod);
  for(size_t i = 0; i < vec->nnz; i++) {
    uint32_t v = nmod_mul(inv, vec->entries[i], mod);
    vec->entries[i] = v;
  }
}
//------------------------------------------------------------------------------
// Compute x - ay mod p
// but leave out the 0th entry of y
// because that will be zero anyway
void inline xmay(int64_t* x, int64_t a, uint32_vec_t y, int64_t p2) {
  size_t j;
  int64_t t;
  for(size_t i = 1; i < y->nnz; i++) {
    j = y->indices[i];
    t = x[j];
    t -= a * y->entries[i];
    t += (t >> 63) & p2;
    x[j] = t;
  }
}
//------------------------------------------------------------------------------
pivots
gauss_elim(uint32_mat_t mat, nmod_t mod, bool* trace) {
  // first canonicalize and compress the matrix
  sparse_mat_compress(mat);

  slong* pivots = new slong[mat->ncol];
  std::fill(pivots, pivots + mat->ncol, -1);

  int64_t* buffer = new int64_t[mat->ncol];
  std::fill(buffer, buffer + mat->ncol, 0);
  uint64_t p = mod.n;
  int64_t p2 = static_cast<int64_t>(p * p);

  std::vector<size_t> buffer_ids;
  buffer_ids.reserve(32);

  for(size_t r = 0; r < mat->nrow; r++) {

    if(trace[r])
      continue;

    auto row = sparse_mat_row(mat, r);
    auto c = row->indices[0];

    // we found a new pivot => rescale and insert in pivots
    if(pivots[c] < 0) {
      normalize_row(row, mod);
      pivots[c] = static_cast<slong>(r);
      continue;
    }

    // we already have a pivot => reduce this row by all pivots
    copy_to_buffer(buffer, row);
    buffer_ids.clear();

    // reduce current row with all pivots
    int64_t cc;
    slong rr;
    for(size_t i = c; i < mat->ncol; i++) {
      cc = buffer[i];
      if(cc == 0)
        continue;
      cc %= p;
      buffer[i] = cc;
      if(cc == 0)
        continue;
      rr = pivots[i];
      if(rr < 0) {
        buffer_ids.push_back(i);
        continue;
      }
      buffer[i] = 0;
      xmay(buffer, cc, sparse_mat_row(mat, rr), p2);
    }

    // we have a zero row
    if(buffer_ids.empty()) {
      sparse_vec_clear(row);
      trace[r] = true;
      continue;
    }

    copy_from_buffer_and_clear(buffer, buffer_ids, row);
    normalize_row(row, mod);
    pivots[row->indices[0]] = static_cast<slong>(r);
  }

  std::vector<size_t> piv;
  for(size_t i = 0; i < mat->ncol; i++)
    if(pivots[i] >= 0)
      piv.push_back(i);

  delete[] buffer;
  delete[] pivots;
  return piv;
}
//------------------------------------------------------------------------------
std::vector<size_t> inline compute_relevant_rows(
  sfmpz_mat_t mat,
  std::vector<sparse_mat_struct<uint32_t>*>& rrefs) {

  std::vector<size_t> relevant_rows;
  boost::unordered_set<ulong> old_pivot_columns;

  for(size_t i = 0; i < mat->nrow; i++)
    old_pivot_columns.insert(sparse_mat_row(mat, i)->indices[0]);

  for(size_t i = 0; i < mat->nrow; i++) {
    for(size_t k = 0; k < rrefs.size(); k++) {
      auto row_ik = sparse_mat_row(rrefs[k], i);
      if(row_ik->nnz == 0)
        continue;
      // test if we have a new leading monomial
      if(old_pivot_columns.find(row_ik->indices[0])
         == old_pivot_columns.end()) {
        relevant_rows.push_back(i);
        break;
      }
    }
  }
  return relevant_rows;
}
//------------------------------------------------------------------------------
std::vector<size_t> inline compute_relevant_rows(sfmpz_mat_t mat,
                                                 uint32_mat_t rref) {

  std::vector<size_t> relevant_rows;
  boost::unordered_set<ulong> old_pivot_columns;

  for(size_t i = 0; i < mat->nrow; i++)
    old_pivot_columns.insert(sparse_mat_row(mat, i)->indices[0]);

  for(size_t i = 0; i < rref->nrow; i++) {
    auto row = sparse_mat_row(rref, i);
    if(row->nnz == 0)
      continue;
    // test if we have a new leading monomial
    if(old_pivot_columns.find(row->indices[0]) == old_pivot_columns.end())
      relevant_rows.push_back(i);
  }
  return relevant_rows;
}
//------------------------------------------------------------------------------
std::pair<std::vector<std::pair<size_t, size_t>>, std::vector<gmp_rational>>
multimodular_gauss_elim(sfmpz_mat_t mat,
                        bool interreduce = false,
                        bool proof = true) {

  std::vector<std::unique_ptr<sparse_mat_struct<uint32_t>>> rrefs;
  pivots best_piv;
  std::vector<pivots> pivs;
  std::vector<ulong> used_primes;
  std::vector<sparse_mat_struct<uint32_t>*> good_rrefs;
  std::vector<pivots> good_pivs;
  std::vector<ulong> good_primes;

  std::vector<gmp_rational> rat_entries;
  std::vector<std::pair<size_t, size_t>> idxs;

  bool* trace = new bool[mat->nrow];
  std::fill(trace, trace + mat->nrow, false);

  size_t i = 0;

  mpz_int h = mpz_int(height(mat));
  mpz_int prod = 1;
  mpz_int M = mat->ncol * 10000000 * (h + 100) * h + 1;

  uint32_t p;
  nmod_t mod;
  bool res;

  while(true) {
    while(prod < M) {
      if(i >= MAX_PRIMES)
        die(-1, "Multimodular Gaussian elimination is not converging");
      p = PRIMES[i++];

      msg("Computing mod %lu", p);

      std::unique_ptr<sparse_mat_struct<uint32_t>> nmod_mat
        = std::make_unique<sparse_mat_struct<uint32_t>>();
      sparse_mat_init(nmod_mat.get(), mat->nrow, mat->ncol);
      nmod_init(&mod, p);
      res = mat_mod(nmod_mat.get(), mat, mod, trace);

      // a pivot was set to zero -- we don't want that
      if(!res) {
        msg("Excluding prime %lu (bad pivots)", p);
        sparse_mat_clear(nmod_mat.get());
        continue;
      }

      KOMMUNOPP_TIME(rref);
      pivots piv(gauss_elim(nmod_mat.get(), mod, trace));
      KOMMUNOPP_PROFILE(timer.~adding_timer());

      if(cmp_pivots(best_piv, piv) <= 0) {
        best_piv = piv;
        pivs.push_back(piv);
        rrefs.push_back(std::move(nmod_mat));
        used_primes.push_back(p);
        prod = prod * p;
      } else {
        msg("Excluding prime %lu (bad pivots)", p);
        sparse_mat_clear(nmod_mat.get());
      }
    }
    prod = 1;
    good_primes.clear();
    good_rrefs.clear();
    good_pivs.clear();
    for(size_t r = 0; r < rrefs.size(); r++) {
      if(cmp_pivots(best_piv, pivs[r]) <= 0) {
        good_primes.push_back(used_primes[r]);
        good_rrefs.push_back(rrefs[r].get());
        good_pivs.push_back(pivs[r]);
        prod = prod * used_primes[r];
      }
    }

    // before reconstruction, compute rows with new leading terms
    std::vector<size_t> relevant_rows;
    if(interreduce)
      for(size_t i = 0; i < mat->nrow; i++)
        relevant_rows.push_back(i);
    else
      relevant_rows = compute_relevant_rows(mat, good_rrefs);

    // Initialize crt_entries in reconstruction
    // and clear here
    fmpz* crt_entries;
    idxs.clear();
    rat_entries.clear();

    {
      KOMMUNOPP_TIME(crt);
      crt_reconstruction(
        crt_entries, idxs, good_rrefs, good_primes, relevant_rows);
    }

    KOMMUNOPP_PROFILE(auto timer = gstats.time(gstats.ratrec));
    bool success
      = rational_reconstruction(rat_entries, crt_entries, prod, idxs.size());
    KOMMUNOPP_PROFILE(timer.~adding_timer());

    for(size_t i = 0; i < idxs.size(); i++)
      fmpz_clear(crt_entries + i);
    delete[] crt_entries;

    if(!success) {
      msg("Reconstruction unsuccessfull. Increasing bound.");
      M = prod * p * p;
      // reset trace
      std::fill(trace, trace + mat->nrow, false);
      continue;
    }

    if(!proof or verify_result())
      break;
  }
  for(const auto& rref : rrefs)
    sparse_mat_clear(rref.get());

  delete[] trace;

  return std::make_pair(idxs, rat_entries);
}
//------------------------------------------------------------------------------

std::pair<std::vector<std::pair<size_t, size_t>>,
          std::vector<gmp_rational>> inline nmod_gauss_elim(sfmpz_mat_t mat,
                                                            size_t p,
                                                            bool interreduce
                                                            = false) {

  // not needed but use it to avoid code duplication
  bool* trace = new bool[mat->nrow];
  std::fill(trace, trace + mat->nrow, false);

  nmod_t mod;
  nmod_init(&mod, p);

  uint32_mat_t nmod_mat;
  sparse_mat_init(nmod_mat, mat->nrow, mat->ncol);
  mat_mod(nmod_mat, mat, mod, trace);

  KOMMUNOPP_TIME(rref);
  gauss_elim(nmod_mat, mod, trace);
  KOMMUNOPP_PROFILE(timer.~adding_timer());

  // compute rows with new leading terms
  std::vector<size_t> relevant_rows;
  if(interreduce)
    for(size_t i = 0; i < mat->nrow; i++)
      relevant_rows.push_back(i);
  else
    relevant_rows = compute_relevant_rows(mat, nmod_mat);

  // compute nonzero indices & entries
  std::vector<std::pair<size_t, size_t>> idxs;
  std::vector<gmp_rational> entries;
  size_t nnz = sparse_mat_nnz(nmod_mat);
  idxs.reserve(nnz);
  entries.reserve(nnz);
  for(size_t i : relevant_rows) {
    auto row = sparse_mat_row(nmod_mat, i);
    for(size_t j = 0; j < row->nnz; j++) {
      idxs.emplace_back(i, row->indices[j]);
      gmp_rational r;
      mpq_set_ui(r.data(), row->entries[j], 1UL);
      entries.push_back(r);
    }
  }

  sparse_mat_clear(nmod_mat);
  delete[] trace;

  return std::make_pair(idxs, entries);
}
//------------------------------------------------------------------------------

std::pair<std::vector<std::pair<size_t, size_t>>,
          std::vector<gmp_rational>> inline linear_algebra(sfmpz_mat_t mat,
                                                           size_t
                                                             characteristic,
                                                           bool interreduce
                                                           = false,
                                                           bool proof = true) {

  // sort rows by first index and nnz
  std::sort(mat->rows, mat->rows + mat->nrow, [](auto& r1, auto& r2) {
    auto id1 = r1.indices[0];
    auto id2 = r2.indices[0];
    if(id1 != id2)
      return id1 > id2;
    return r1.nnz < r2.nnz;
  });

  if(characteristic == 0)
    return multimodular_gauss_elim(mat, interreduce, proof);
  else
    return nmod_gauss_elim(mat, characteristic, interreduce);
}
}
