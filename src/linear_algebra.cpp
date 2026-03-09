#include <algorithm>
#include <atomic>
#include <boost/unordered/unordered_map.hpp>
#include <cstddef>
#include <cstdint>
#include <map>
#include <ranges>
#include <set>
#include <span>
#include <utility>
#include <vector>

#include "coeff.hpp"
#include "fast_div.hpp"
#include "primes.hpp"
#include "profiling.hpp"
#include "signal_statistics.hpp"
#include "sparse_rref/sparse_mat.h"
#include "store.hpp"

#include "linear_algebra.hpp"

using namespace boost::multiprecision;

extern int verbose;

namespace f4ncgb {

static bool use_trace;
uint32_mat_ro_t red_mat;
boost::unordered_map<const coeff*, size_t> ref_rows;

void inline set_up_reducer_mat(const std::vector<std::vector<size_t>>& idxs,
                               const std::vector<std::span<coeff>>& entries) {
  const size_t m = idxs.size();

  ref_rows.clear();

  // first pass: compute total nnz and unique_nnz
  // also remember first row for each unique span
  size_t nnz = 0;
  size_t unique_nnz = 0;
  for(size_t i = 0; i < m; i++) {
    nnz += idxs[i].size();
    const coeff* ptr = entries[i].data();
    if(ref_rows.find(ptr) == ref_rows.end()) {
      ref_rows[ptr] = i;
      unique_nnz += idxs[i].size();
    }
  }

  sparse_mat_ro_init(red_mat, m, nr_cols, nnz, unique_nnz);
  double ratio = nnz > 0 ? (100. * unique_nnz) / nnz : 0;
  if(verbose > 2)
    msg("Reducer matrix has size (%d, %d) (%.1f%% entries actually stored)",
        red_mat->nrow,
        red_mat->ncol,
        ratio);

  size_t idx_k = 0;
  size_t entry_k = 0;
  uint32_t* indices = red_mat->indices;
  uint32_t* index_start = red_mat->index_start;
  uint32_t* row_start = red_mat->row_start;
  uint32_t* row_nnz = red_mat->row_nnz;

  for(size_t i = 0; i < m; i++) {
    index_start[i] = idx_k;
    row_nnz[i] = idxs[i].size();

    // TODO: transformation matrix
    // store indices
    for(auto j : idxs[i])
      indices[idx_k++] = j;

    // store coeff offset
    const coeff* ptr = entries[i].data();
    size_t ref_row = ref_rows[ptr];
    if(i == ref_row) {
      row_start[i] = entry_k;
      entry_k += idxs[i].size();// TODO: transformation matrix
    } else {
      row_start[i] = row_start[ref_row];
    }
  }
}
//------------------------------------------------------------------------------
bool inline fill_reducer_mat(const std::vector<std::span<coeff>>& entries,
                             nmod_t mod) {

  uint32_t* mat_entries = red_mat->entries;

  // only fill entries for first occurrence of each span
  for(size_t i = 0; i < red_mat->nrow; i++) {
    const auto& coeffs = entries[i];
    const coeff* ptr = coeffs.data();

    // just a duplicate -> nothing to do
    if(ref_rows[ptr] != i)
      continue;

    // check leading coeff for zero
    uint32_t cc = fmpz_get_nmod(coeffs[0].value, mod);
    if(!cc)
      return false;

    // only write entries for first occurrence
    size_t k = red_mat->row_start[i];
    mat_entries[k++] = 1;
    uint32_t inv = nmod_inv(cc, mod);
    for(size_t j = 1; j < coeffs.size(); j++) {
      cc = fmpz_get_nmod(coeffs[j].value, mod);
      mat_entries[k++] = nmod_mul(inv, cc, mod);
    }
  }
  return true;
}

//------------------------------------------------------------------------------
std::vector<int32_t> red_pivots;
void inline set_reducer_pivots() {
  // clear
  red_pivots.resize(red_mat->ncol);
  std::fill(red_pivots.begin(), red_pivots.end(), -1);

  for(uint32_t i = 0; i < red_mat->nrow; i++) {
    uint32_t j = red_mat->indices[red_mat->index_start[i]];
    red_pivots[j] = static_cast<int32_t>(i);
  }
}

//------------------------------------------------------------------------------

bool
set_up_matrix(uint32_mat_t mat,
              nmod_t mod,
              const std::vector<std::vector<size_t>>& idxs,
              const std::vector<std::span<coeff>>& entries) {
  const size_t m = idxs.size();

  // ---- initialize matrix ----
  if(proof_level > 0)
    sparse_mat_init(mat, m, nr_cols + m);
  else
    sparse_mat_init(mat, m, nr_cols);

  if(verbose > 2)
    msg("S-Pol matrix has size (%d, %d)", mat->nrow, mat->ncol);

  // ---- fill rows ----
  for(size_t i = 0; i < m; i++) {

    if(red_mat->trace[i]) {
      // dummy init to be safe
      sparse_mat_row_init(mat, i, 0);
      continue;
    }

    const auto& coeffs = entries[i];
    const auto& cols = idxs[i];

    const size_t max_nnz = coeffs.size() + (proof_level > 0 ? 1 : 0);

    auto row = sparse_mat_row_init(mat, i, max_nnz);

    size_t nnz = 0;

    // --- insert reduced entries ---
    for(size_t k = 0; k < coeffs.size(); k++) {
      uint32_t c = fmpz_get_nmod(coeffs[k].value, mod);
      if(c != 0) {
        row->indices[nnz] = cols[k];
        row->entries[nnz] = c;
        nnz++;
      }
      // ---- detect zero leading entry ----
      if(nnz == 0)
        return false;
    }

    // ---- transformation block ----
    if(proof_level > 0) {
      row->indices[nnz] = nr_cols + i;
      // insert denom from input polynomial
      if(interreduce) {
        uint32_t c = fmpz_get_nmod(input_denoms[i].value, mod);
        row->entries[nnz] = c;
      } else
        row->entries[nnz] = 1;
      nnz++;
    }

    row->nnz = nnz;
  }
  return true;
}
//------------------------------------------------------------------------------

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

static inline coeff
height(const std::vector<std::span<coeff>>& entries_spol,
       const std::vector<std::span<coeff>>& entries_red) {

  coeff h(0L);
  for(const auto& row : entries_spol) {
    for(const auto& c : row) {
      if(fmpz_cmpabs(c.value, h.value) > 0)
        fmpz_abs(h.value, c.value);
    }
  }
  for(const auto& row : entries_red) {
    for(const auto& c : row) {
      if(fmpz_cmpabs(c.value, h.value) > 0)
        fmpz_abs(h.value, c.value);
    }
  }

  return h;
}
//------------------------------------------------------------------------------

static void
crt_reconstruction(fmpz*& entries,
                   std::vector<std::pair<size_t, size_t>>& idxs,
                   std::vector<sparse_mat_struct<uint32_t>*>& rrefs,
                   std::vector<uint32_t>& primes,
                   size_t n_piv) {
  if(rrefs.size() == 0)
    return;

  // get all (i,j) where at least one rref is nonzero
  std::map<size_t, std::set<size_t>> nnz_pos;
  for(size_t i = 0; i < n_piv; i++) {
    auto& nnz_pos_row = nnz_pos[i];
    for(auto& rref : rrefs) {
      auto row = sparse_mat_row(rref, i);
      // include row only if polynomial part is nonzero
      if(proof_level == 0 or row->nnz > 0)
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
  for(size_t i = 0; i < n_piv; i++) {
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
  fmpz_cleanup(moduli, primes.size());
  fmpz_cleanup(inputs, rrefs.size());
}
//------------------------------------------------------------------------------

inline bool
verify_result(fmpz* nums,
              fmpz* denoms,
              size_t len,
              coeff height,
              size_t n,
              coeff P) {
  fmpz_t d;
  fmpz_init_set_ui(d, 1);
  for(size_t i = 0; i < len; i++)
    fmpz_lcm(d, d, denoms + i);

  // height_rref = max | (d / denom_i) * num_i |
  fmpz_t tmp;
  fmpz_t height_rref;
  fmpz_init(tmp);
  fmpz_init_set_ui(height_rref, 0);

  for(size_t i = 0; i < len; i++) {
    fmpz_divexact(tmp, d, denoms + i);
    fmpz_mul(tmp, tmp, nums + i);
    fmpz_abs(tmp, tmp);
    if(fmpz_cmp(tmp, height_rref) > 0)
      fmpz_set(height_rref, tmp);
  }

  // Now compute lhs = height_rref * height * n
  fmpz_t lhs_fmpz;
  fmpz_init(lhs_fmpz);

  fmpz_mul(lhs_fmpz, height_rref, height.value);
  fmpz_mul_ui(lhs_fmpz, lhs_fmpz, n);

  bool result = (fmpz_cmp(lhs_fmpz, P.value) < 0);

  // cleanup
  fmpz_clear(d);
  fmpz_clear(tmp);
  fmpz_clear(height_rref);
  fmpz_clear(lhs_fmpz);

  return result;
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
ratrecon(fmpz_t num, fmpz_t den, mpz_t u, ratrec_data* data) {
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
    fmpz_set_mpz(num, data->n);
    fmpz_set_mpz(den, data->d);

  } else if(verbose > 2)
    msg("Rational reconstruction does not exist");

  return success;
}
//------------------------------------------------------------------------------

static inline bool
rational_reconstruction(fmpz*& nums,
                        fmpz*& denoms,
                        fmpz* crt_entries,
                        coeff& prod,
                        size_t N) {
  ratrec_data data;
  fmpz_get_mpz(data.mod, prod.value);
  // N = floor(sqrt(m/2))
  mpz_fdiv_q_2exp(data.N, data.mod, 1);
  mpz_sqrt(data.N, data.N);

  nums = new fmpz[N];
  denoms = new fmpz[N];
  for(size_t i = 0; i < N; i++) {
    fmpz_init(nums + i);
    fmpz_init(denoms + i);
  }

  bool res = true;
  mpz_t u;
  mpz_init(u);
  for(size_t i = 0; i < N; i++) {
    fmpz_get_mpz(u, crt_entries + i);
    bool success = ratrecon(nums + i, denoms + i, u, &data);
    if(!success) {
      res = false;
      break;
    }
  }
  mpz_clear(u);
  return res;
}
//------------------------------------------------------------------------------

void
clear_denoms(std::vector<std::pair<size_t, size_t>>& idxs,
             fmpz* nums,
             fmpz* denoms) {

  if(idxs.empty())
    return;

  fmpz_t L;
  fmpz_t tmp;
  fmpz_init(L);
  fmpz_init(tmp);

  size_t cur_row = idxs[0].first;
  size_t row_start = 0;

  while(row_start < idxs.size()) {

    cur_row = idxs[row_start].first;
    fmpz_set_ui(L, 1);

    // find row range and compute denom-lcm
    size_t row_end = row_start;
    while(row_end < idxs.size() && idxs[row_end].first == cur_row) {
      fmpz_lcm(L, L, denoms + row_end);
      row_end++;
    }

    // scale nums if necessary
    if(!fmpz_is_one(L)) {
      for(size_t i = row_start; i < row_end; i++) {
        fmpz_divexact(tmp, L, denoms + i);
        fmpz_mul(nums + i, nums + i, tmp);
      }
    }

    row_start = row_end;
  }

  fmpz_clear(L);
  fmpz_clear(tmp);
}
//------------------------------------------------------------------------------

template<typename T>
static void inline copy_to_buffer(T& buffer, uint32_vec_t vec) {
  for(size_t i = 0; i < vec->nnz; i++)
    buffer[vec->indices[i]] = vec->entries[i];
}
//------------------------------------------------------------------------------

static void inline copy_from_buffer_and_clear(int64_t* buffer,
                                              std::vector<uint32_t>& buffer_ids,
                                              uint32_vec_t vec) {
  size_t nnz = buffer_ids.size();

  sparse_vec_realloc(vec, nnz);
  vec->nnz = nnz;

  std::move(buffer_ids.begin(), buffer_ids.end(), vec->indices);
  size_t j = 0;
  for(size_t i : buffer_ids) {
    assert(buffer[i] != 0);
    vec->entries[j++] = static_cast<uint32_t>(buffer[i]);
    buffer[i] = 0;
  }
}
//------------------------------------------------------------------------------

static void inline normalize_row(uint32_vec_t vec, nmod_t mod) {
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
static inline void
xmay(int64_t* __restrict x,
     int64_t a,
     const uint32_t* __restrict idx,
     const uint32_t* __restrict val,
     uint32_t nnz,
     int64_t p2) {
  uint32_t i = 1;
  for(; i + 3 < nnz; i += 4) {
    uint32_t j0 = idx[i];
    uint32_t j1 = idx[i + 1];
    uint32_t j2 = idx[i + 2];
    uint32_t j3 = idx[i + 3];

    int64_t v0 = x[j0];
    int64_t v1 = x[j1];
    int64_t v2 = x[j2];
    int64_t v3 = x[j3];

    int64_t t0 = v0 - a * val[i];
    int64_t t1 = v1 - a * val[i + 1];
    int64_t t2 = v2 - a * val[i + 2];
    int64_t t3 = v3 - a * val[i + 3];

    t0 += (t0 >> 63) & p2;
    t1 += (t1 >> 63) & p2;
    t2 += (t2 >> 63) & p2;
    t3 += (t3 >> 63) & p2;

    x[j0] = t0;
    x[j1] = t1;
    x[j2] = t2;
    x[j3] = t3;
  }

  // Handle remaining elements
  for(; i < nnz; i++) {
    uint32_t j = idx[i];
    int64_t t = x[j] - a * val[i];
    t += (t >> 63) & p2;
    x[j] = t;
  }
}
//------------------------------------------------------------------------------

template<bool mersenne = false>
static pivots
reverse_solve(uint32_mat_t mat, nmod_t mod) {
  uint64_t p = mod.n;
  int64_t p2 = static_cast<int64_t>(p * p);

  std::vector<int32_t> piv_array(mat->ncol, -1);
  std::vector<int64_t> buffer(mat->ncol, 0);
  std::vector<uint32_t> buffer_ids;
  int64_t* buff = buffer.data();
  buffer_ids.reserve(32);

  // sort rows by first index and nnz
  // assume: mat is in ref
  std::sort(mat->rows, mat->rows + mat->nrow, [](auto& r1, auto& r2) {
    // move zero rows down
    if(r1.nnz == 0 or r2.nnz == 0)
      return r1.nnz > r2.nnz;
    // sort by pivots
    return r1.indices[0] > r2.indices[0];
  });

  // collect pivots -- assume: mat is in ref
  pivots piv;
  piv.reserve(64);
  for(size_t r = 0; r < mat->nrow; r++) {
    auto row = sparse_mat_row(mat, r);
    if(row->nnz == 0)
      break;
    piv_array[row->indices[0]] = static_cast<int32_t>(r);
    piv.push_back(row->indices[0]);
  }

  // reduce the new pivot rows fully
  for(size_t r = 0; r < piv.size(); r++) {
    auto row = sparse_mat_row(mat, r);
    // skip rows with only one entry
    if(row->nnz < 2)
      continue;

    copy_to_buffer(buff, row);
    buffer_ids.clear();
    buffer_ids.push_back(row->indices[0]);
    int64_t cc;
    int32_t rr;
    uint32_t i = row->indices[1];
    while(i < mat->ncol) {
#ifdef F4NCGB_FAST_MERSENNE_PRIME_MODULO
      if constexpr(mersenne) {
        cc = mersenne_mod(buff[i]);
      } else {
#endif
        cc = (int64_t)(((uint64_t)buff[i]) % p);
#ifdef F4NCGB_FAST_MERSENNE_PRIME_MODULO
      }
#endif
      buff[i] = cc;
      if(cc != 0) {
        rr = piv_array[i];
        if(rr < 0) {
          buffer_ids.push_back(i);
        } else {
          assert(sparse_mat_row(mat, rr)->indices[0] == i);
          assert(sparse_mat_row(mat, rr)->entries[0] == 1);
          buff[i] = 0;
          auto row_rr = sparse_mat_row(mat, rr);
          xmay(buff, cc, row_rr->indices, row_rr->entries, row_rr->nnz, p2);
        }
      }
      for(i++; i < mat->ncol; i++)
        if(buff[i] != 0)
          break;
    }
    copy_from_buffer_and_clear(buff, buffer_ids, row);
  }

  return piv;
}
//------------------------------------------------------------------------------

template<bool mersenne = false>
pivots
gauss_elim(uint32_mat_t mat,
           nmod_t mod,
           std::unique_ptr<BS::thread_pool<BS::none>>& pool) {

  uint64_t p = mod.n;
  int64_t p2 = static_cast<int64_t>(p * p);

  thread_local std::vector<int64_t> buffer_local;
  thread_local std::vector<uint32_t> buffer_ids_local;

  std::vector<std::atomic_int32_t> atomic_pivots(mat->ncol);
  std::fill(atomic_pivots.begin(), atomic_pivots.end(), -1);

  for(size_t r = 0; r < mat->nrow; r++) {
    if(red_mat->trace[r])
      continue;

    // reduce row with all already known pivots
    auto task = [r, &mat, &atomic_pivots, p, p2, &mod]() {
      F4NCGB_TIME(elim_task_cpu);
      auto row = sparse_mat_row(mat, r);
      int32_t rr;
      int64_t cc;
      int32_t expected = -1;
      uint32_t* indices;
      uint32_t* entries;
      uint32_t nnz;

      buffer_local.resize(mat->ncol, 0);
      int64_t* buff = buffer_local.data();

      do {
        copy_to_buffer(buff, row);
        buffer_ids_local.clear();
        uint32_t i = row->indices[0];
        uint32_t max_col = row->indices[row->nnz - 1];
        while(i <= max_col) {
          // find next nonzero entry
          for(; i <= max_col; i++)
            if(buff[i] != 0)
              break;

          if(i > max_col)
            break;

          assert(buff[i] > 0);
          // v must be smaller than 2^2b, i.e. 2^62
          assert(buff[i] < 4611686018427387904);
#ifdef F4NCGB_FAST_MERSENNE_PRIME_MODULO
          if constexpr(mersenne) {
            (void)p;// p is not used in this case.
            cc = mersenne_mod(buff[i]);
          } else
#endif
            cc = (int64_t)(((uint64_t)buff[i]) % p);
          buff[i] = cc;
          if(cc == 0)
            continue;

          // first check red pivots
          rr = red_pivots[i];
          if(rr >= 0) {
            buff[i] = 0;
            indices = red_mat_indices(red_mat, rr);
            entries = red_mat_entries(red_mat, rr);
            nnz = red_mat_nnz(red_mat, rr);
            max_col = std::max(max_col, indices[nnz - 1]);
            xmay(buff, cc, indices, entries, nnz, p2);
            i++;
            continue;
          }

          // then check atomic pivots
          rr = atomic_pivots[i];
          if(rr >= 0) {
            buff[i] = 0;
            auto row = sparse_mat_row(mat, rr);
            indices = row->indices;
            entries = row->entries;
            nnz = row->nnz;
            max_col = std::max(max_col, indices[nnz - 1]);
            xmay(buff, cc, indices, entries, nnz, p2);
            i++;
            continue;
          }

          // no reducer found
          buffer_ids_local.push_back(i);
          i++;
        }

        // we have a zero row
        if(buffer_ids_local.empty()
           or (proof_level > 0
               and buffer_ids_local[0] >= mat->ncol - mat->nrow)) {
          sparse_vec_clear(row);
          for(auto id : buffer_ids_local)
            buff[id] = 0;
          red_mat->trace[r] = use_trace;// only set them if we use tracer
          return;
        }
        copy_from_buffer_and_clear(buff, buffer_ids_local, row);
        normalize_row(row, mod);
        expected = -1;
      } while(!atomic_pivots[row->indices[0]].compare_exchange_weak(
        expected, static_cast<int32_t>(r)));
    };
    if(pool)
      pool->detach_task(task);
    else
      task();
  }
  if(pool)
    pool->wait();

  if(mod.n == PRIMES[0]) {
    return reverse_solve<true>(mat, mod);
  } else {
    return reverse_solve<false>(mat, mod);
  }
}
//------------------------------------------------------------------------------

std::pair<std::vector<std::pair<size_t, size_t>>, fmpz*>
multimodular_gauss_elim(std::vector<std::vector<size_t>>& idxs_spol,
                        std::vector<std::span<coeff>>& entries_spol,
                        std::vector<std::span<coeff>>& entries_red,
                        std::unique_ptr<BS::thread_pool<BS::none>>& pool) {

  std::vector<std::unique_ptr<sparse_mat_struct<uint32_t>>> rrefs;
  pivots best_piv;
  std::vector<pivots> pivs;
  std::vector<uint32_t> used_primes;
  std::vector<sparse_mat_struct<uint32_t>*> good_rrefs;
  std::vector<pivots> good_pivs;
  std::vector<uint32_t> good_primes;

  std::vector<std::pair<size_t, size_t>> idxs;
  fmpz* nums = nullptr;
  fmpz* denoms = nullptr;

  size_t nrow = entries_spol.size() + entries_red.size();
  size_t ncol = 0;

  size_t i = 0;

  coeff h = height(entries_spol, entries_red);
  coeff prod(1L);
  coeff M = 20000000 * nrow * (h + 100) * h + 1;

  uint32_t p = PRIMES[0];
  nmod_t mod;
  bool res;

  while(true) {
    while(prod < M) {
      if(i >= MAX_PRIMES)
        die(-1, "Multimodular Gaussian elimination is not converging");
      p = PRIMES[i++];
      nmod_init(&mod, p);

      if(verbose > 2)
        msg("Computing mod %lu", p);

      // set up the matrices
      res = fill_reducer_mat(entries_red, mod);
      std::unique_ptr<sparse_mat_struct<uint32_t>> nmod_mat
        = std::make_unique<sparse_mat_struct<uint32_t>>();
      res &= set_up_matrix(nmod_mat.get(), mod, idxs_spol, entries_spol);

      // a pivot was set to zero -- we don't want that
      if(!res) {
        if(verbose > 2)
          msg("Excluding prime %lu (bad pivots)", p);
        sparse_mat_clear(nmod_mat.get());
        continue;
      }

      ncol = nmod_mat.get()->ncol;

      F4NCGB_TIME(rref);
      auto gauss_elim_wrapper = [&]() -> pivots {
        if(mod.n == PRIMES[0]) {
          return gauss_elim<true>(nmod_mat.get(), mod, pool);
        } else {
          return gauss_elim<false>(nmod_mat.get(), mod, pool);
        }
      };
      pivots piv(gauss_elim_wrapper());

      if(cmp_pivots(best_piv, piv) <= 0) {
        best_piv = piv;
        pivs.push_back(piv);
        rrefs.push_back(std::move(nmod_mat));
        used_primes.push_back(p);
        prod = prod * p;
      } else {
        if(verbose > 2)
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
    // Initialize crt_entries in reconstruction
    // and clear here
    fmpz* crt_entries = nullptr;
    idxs.clear();

    {
      F4NCGB_TIME(crt);
      crt_reconstruction(
        crt_entries, idxs, good_rrefs, good_primes, best_piv.size());
    }

    bool success;
    {
      F4NCGB_PROFILE(auto timer = gstats.time(gstats.ratrec));
      success
        = rational_reconstruction(nums, denoms, crt_entries, prod, idxs.size());
    }

    fmpz_cleanup(crt_entries, idxs.size());

    if(!success) {
      if(verbose > 2)
        msg("Reconstruction unsuccessful. Increasing bound.");
      M = prod * p * p;
      // reset trace and cleanup
      fmpz_cleanup(nums, idxs.size());
      fmpz_cleanup(denoms, idxs.size());
      continue;
    }

    if(verify_result(nums, denoms, idxs.size(), h, ncol, prod))
      break;
  }

  for(const auto& rref : rrefs)
    sparse_mat_clear(rref.get());
  clear_denoms(idxs, nums, denoms);
  fmpz_cleanup(denoms, idxs.size());

  return std::make_pair(idxs, nums);
}
//------------------------------------------------------------------------------

std::pair<std::vector<std::pair<size_t, size_t>>, fmpz*>
nmod_gauss_elim(std::vector<std::vector<size_t>>& idxs_spol,
                std::vector<std::span<coeff>>& entries_spol,
                std::vector<std::span<coeff>>& entries_red,
                size_t p,
                std::unique_ptr<BS::thread_pool<BS::none>>& pool) {
  nmod_t mod;
  nmod_init(&mod, p);

  uint32_mat_t nmod_mat;
  set_up_matrix(nmod_mat, mod, idxs_spol, entries_spol);
  fill_reducer_mat(entries_red, mod);

  F4NCGB_TIME(rref);
  auto gauss_elim_wrapper = [&]() -> pivots {
    if(mod.n == PRIMES[0]) {
      return gauss_elim<true>(nmod_mat, mod, pool);
    } else {
      return gauss_elim<false>(nmod_mat, mod, pool);
    }
  };
  pivots piv(gauss_elim_wrapper());

  // compute nonzero indices & entries
  size_t nnz = sparse_mat_nnz(nmod_mat);
  std::vector<std::pair<size_t, size_t>> idxs;
  idxs.reserve(nnz);
  fmpz* entries = new fmpz[nnz];

  size_t k = 0;
  for(size_t i = 0; i < nmod_mat->nrow; i++) {
    auto row = sparse_mat_row(nmod_mat, i);
    for(size_t j = 0; j < row->nnz; j++) {
      idxs.emplace_back(i, row->indices[j]);
      fmpz_set_ui(entries + k++, row->entries[j]);
    }
  }

  sparse_mat_clear(nmod_mat);

  return std::make_pair(idxs, entries);
}
//------------------------------------------------------------------------------

std::pair<std::vector<std::pair<size_t, size_t>>, fmpz*>
linear_algebra(std::vector<std::vector<size_t>>& idxs_spol,
               std::vector<std::vector<size_t>>& idxs_red,
               std::vector<std::span<coeff>>& entries_spol,
               std::vector<std::span<coeff>>& entries_red,
               size_t characteristic,
               std::unique_ptr<BS::thread_pool<BS::none>>& pool,
               bool tracer) {

  use_trace = tracer;

  set_up_reducer_mat(idxs_red, entries_red);
  set_reducer_pivots();
  sparse_mat_init_trace(red_mat, idxs_spol.size());

  std::pair<std::vector<std::pair<size_t, size_t>>, fmpz*> res;
  if(characteristic == 0)
    res = multimodular_gauss_elim(idxs_spol, entries_spol, entries_red, pool);

  else
    res = nmod_gauss_elim(
      idxs_spol, entries_spol, entries_red, characteristic, pool);

  sparse_mat_ro_clear(red_mat);
  return res;
}
}
