#include <algorithm>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <map>
#include <set>
#include <sys/errno.h>
#include <utility>
#include <vector>

#include <boost/multiprecision/gmp.hpp>
#include <boost/unordered_set.hpp>

#include "fast_div.hpp"
#include "primes.hpp"
#include "profiling.hpp"
#include "signal_statistics.hpp"
#include "sparse_rref/thread_pool.hpp"

#include "linear_algebra.hpp"

using namespace boost::multiprecision;

namespace kommunopp {

static inline bool
vec_mod(uint32_vec_t vec,
        const sfmpz_vec_t src,
        nmod_t mod,
        std::vector<size_t>& indices,
        std::vector<uint32_t>& entries) {
  assert(vec->nnz <= src->nnz);

  indices.clear();
  entries.clear();

  for(size_t i = 0; i < src->nnz; i++) {
    uint32_t val = fmpz_get_nmod(src->entries + i, mod);
    if(val != 0) {
      indices.push_back(src->indices[i]);
      entries.push_back(val);
    }
  }

  auto nnz = indices.size();
  if(nnz == 0 or indices[0] != src->indices[0])
    return false;

  sparse_vec_realloc(vec, nnz);
  vec->nnz = nnz;
  std::move(indices.begin(), indices.end(), vec->indices);
  std::move(entries.begin(), entries.end(), vec->entries);

  return true;
}

static inline bool
mat_mod(uint32_mat_t mat, const sfmpz_mat_t src, nmod_t mod, bool* trace) {
  bool res = true;
  std::vector<size_t> indices;
  std::vector<uint32_t> entries;
  for(size_t i = 0; i < src->nrow; i++) {
    if(!trace[i]) {
      res = vec_mod(
        sparse_mat_row(mat, i), sparse_mat_row(src, i), mod, indices, entries);
      if(!res)
        break;
    }
  }
  return res;
}

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

static inline gmp_int
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

static void
crt_reconstruction(fmpz*& entries,
                   std::vector<std::pair<size_t, size_t>>& idxs,
                   std::vector<sparse_mat_struct<uint32_t>*>& rrefs,
                   std::vector<ulong>& primes,
                   size_t n_piv) {
  if(rrefs.size() == 0)
    return;

  // get all (i,j) where at least one rref is nonzero
  std::map<size_t, std::set<size_t>> nnz_pos;
  for(size_t i = 0; i < n_piv; i++) {
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
  for(size_t i = 0; i < primes.size(); i++)
    fmpz_clear(moduli + i);
  for(size_t i = 0; i < rrefs.size(); i++)
    fmpz_clear(inputs + i);
  delete[] moduli;
  delete[] inputs;
}

inline bool
verify_result() {
  return true;
}
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

static inline bool
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

template<typename T>
static void inline copy_to_buffer(T* buffer, uint32_vec_t vec) {
  for(size_t i = 0; i < vec->nnz; i++)
    buffer[vec->indices[i]] = vec->entries[i];
}

static void inline copy_from_buffer_and_clear(int64_t* buffer,
                                              std::vector<size_t>& buffer_ids,
                                              uint32_vec_t vec) {
  size_t nnz = buffer_ids.size();

  sparse_vec_clear(vec);
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

static void inline normalize_row(uint32_vec_t vec, nmod_t mod) {
  uint32_t c = vec->entries[0];
  uint32_t inv = nmod_inv(c, mod);
  for(size_t i = 0; i < vec->nnz; i++) {
    uint32_t v = nmod_mul(inv, vec->entries[i], mod);
    vec->entries[i] = v;
  }
}

// Compute x - ay mod p
// but leave out the 0th entry of y
// because that will be zero anyway
static void inline xmay(int64_t* x, int64_t a, uint32_vec_t y, int64_t p2) {
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

static inline bool
is_ref(uint32_mat_t mat) {

  std::unordered_map<size_t, size_t> pivs;

  for(size_t i = 0; i < mat->nrow; i++) {
    auto row = sparse_mat_row(mat, i);
    if(row->nnz == 0)
      continue;
    pivs[row->indices[0]] += 1;
  }

  for(auto [k, v] : pivs)
    if(v > 1)
      return false;

  return true;
}

static inline bool
is_rref(uint32_mat_t mat) {

  if(!is_ref(mat))
    return false;

  std::set<size_t> pivs;

  for(size_t i = 0; i < mat->nrow; i++) {
    auto row = sparse_mat_row(mat, i);
    if(row->nnz == 0)
      continue;
    pivs.insert(row->indices[0]);
  }

  for(size_t i = 0; i < mat->nrow; i++) {
    auto row = sparse_mat_row(mat, i);
    if(row->nnz == 0 or !row->is_new_piv)
      continue;
    for(size_t j = 1; j < row->nnz; j++)
      if(pivs.contains(row->indices[j]))
        return false;
  }

  return true;
}

pivots
reverse_solve(uint32_mat_t mat, nmod_t mod) {
  uint64_t p = mod.n;
  int64_t p2 = static_cast<int64_t>(p * p);

  int64_t* piv_array = new int64_t[mat->ncol];
  std::fill(piv_array, piv_array + mat->ncol, -1);

  int64_t* buffer = new int64_t[mat->ncol];
  std::fill(buffer, buffer + mat->ncol, 0);

  std::vector<size_t> buffer_ids;
  buffer_ids.reserve(32);

  // sort new pivot rows up -- assume: maat is in ref
  // sort rows by first index and nnz
  std::sort(mat->rows, mat->rows + mat->nrow, [](auto& r1, auto& r2) {
    // move new pivot rows up
    if(r1.is_new_piv != r2.is_new_piv)
      return r1.is_new_piv;
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
    piv_array[row->indices[0]] = static_cast<int64_t>(r);
    if(row->is_new_piv)
      piv.push_back(row->indices[0]);
  }

  // reduce the new pivot rows fully
  for(size_t r = 0; r < piv.size(); r++) {
    auto row = sparse_mat_row(mat, r);
    // skip rows with only one entry
    if(row->nnz < 2)
      continue;

    copy_to_buffer(buffer, row);
    buffer_ids.clear();
    buffer_ids.push_back(row->indices[0]);
    int64_t cc;
    slong rr;
    size_t i = row->indices[1];
    while(i < mat->ncol) {
      buffer[i] %= p;
      cc = buffer[i];
      if(cc != 0) {
        rr = piv_array[i];
        if(rr < 0) {
          buffer_ids.push_back(i);
        } else {
          assert(sparse_mat_row(mat, rr)->indices[0] == i);
          buffer[i] = 0;
          xmay(buffer, cc, sparse_mat_row(mat, rr), p2);
        }
      }
      i++;
      while(i < mat->ncol and buffer[i] == 0)
        i++;
    }
    copy_from_buffer_and_clear(buffer, buffer_ids, row);
  }

  assert(is_rref(mat));

  delete[] buffer;
  delete[] piv_array;
  return piv;
}

pivots
gauss_elim(uint32_mat_t mat, nmod_t mod, size_t num_threads, bool* trace) {
  uint64_t p = mod.n;
  int64_t p2 = static_cast<int64_t>(p * p);

  thread_local int64_t* buffer_local;
  thread_local std::vector<size_t> buffer_ids_local;

  std::vector<std::atomic_int64_t> atomic_pivots(mat->ncol);
  std::fill(atomic_pivots.begin(), atomic_pivots.end(), -1);

  BS::thread_pool pool(num_threads, [n = mat->ncol] {
    buffer_local = new int64_t[n];
    std::fill(buffer_local, buffer_local + n, 0);
    buffer_ids_local.reserve(32);
  });
  pool.set_cleanup_func([]() { delete[] buffer_local; });

  for(size_t r = 0; r < mat->nrow; r++) {
    if(trace[r])
      continue;

    auto row = sparse_mat_row(mat, r);
    size_t c = row->indices[0];

    // we found a new pivot => rescale and insert in pivots
    int64_t rr = atomic_pivots[c];
    if(rr < 0) {
      normalize_row(row, mod);
      assert(atomic_pivots[c] == -1);
      atomic_pivots[c] = static_cast<int64_t>(r);
      continue;
    }

    // reduce row with all already known pivots
    if(p == 2147483647) {
      // This might be cleaned up a bit, but the mersenne prime case in general
      // can be optimized way better.
      pool.detach_task([r, &mat, &atomic_pivots, &trace, &p2, &mod]() {
        KOMMUNOPP_TIME(elim_task_cpu);
        auto row = sparse_mat_row(mat, r);

        int64_t rr;
        int64_t cc;
        int64_t expected = -1;
        do {
          copy_to_buffer(buffer_local, row);
          buffer_ids_local.clear();
          size_t i = row->indices[0];
          while(i < mat->ncol) {
            assert(buffer_local[i] > 0);
            // v must be smaller than 2^2b, i.e. 2^62
            assert(buffer_local[i] < 4611686018427387904);
            buffer_local[i] = v_mod_2_31_1((uint64_t)buffer_local[i]);
            cc = buffer_local[i];
            if(cc != 0) {
              rr = atomic_pivots[i];
              if(rr < 0) {
                buffer_ids_local.push_back(i);
              } else {
                buffer_local[i] = 0;
                xmay(buffer_local, cc, sparse_mat_row(mat, rr), p2);
              }
            }
            i++;
            while(i < mat->ncol and buffer_local[i] == 0)
              i++;
          }
          // we have a zero row
          if(buffer_ids_local.empty()) {
            sparse_vec_clear(row);
            trace[r] = true;
            return;
          }

          copy_from_buffer_and_clear(buffer_local, buffer_ids_local, row);
          normalize_row(row, mod);
          expected = -1;
        } while(!atomic_pivots[row->indices[0]].compare_exchange_weak(
          expected, static_cast<int64_t>(r)));
        row->is_new_piv = true;
      });
    } else {
      pool.detach_task([r, &mat, &atomic_pivots, &trace, &p, &p2, &mod]() {
        KOMMUNOPP_TIME(elim_task_cpu);
        auto row = sparse_mat_row(mat, r);

        int64_t rr;
        int64_t cc;
        int64_t expected = -1;
        do {
          copy_to_buffer(buffer_local, row);
          buffer_ids_local.clear();
          size_t i = row->indices[0];
          while(i < mat->ncol) {
            assert(buffer_local[i] != 0);
            buffer_local[i] %= p;
            cc = buffer_local[i];
            if(cc != 0) {
              rr = atomic_pivots[i];
              if(rr < 0) {
                buffer_ids_local.push_back(i);
              } else {
                buffer_local[i] = 0;
                xmay(buffer_local, cc, sparse_mat_row(mat, rr), p2);
              }
            }
            i++;
            while(i < mat->ncol and buffer_local[i] == 0)
              i++;
          }
          // we have a zero row
          if(buffer_ids_local.empty()) {
            sparse_vec_clear(row);
            trace[r] = true;
            return;
          }

          copy_from_buffer_and_clear(buffer_local, buffer_ids_local, row);
          normalize_row(row, mod);
          expected = -1;
        } while(!atomic_pivots[row->indices[0]].compare_exchange_weak(
          expected, static_cast<int64_t>(r)));
        row->is_new_piv = true;
      });
    }
  }
  pool.wait();

  return reverse_solve(mat, mod);
}

std::pair<std::vector<std::pair<size_t, size_t>>, std::vector<gmp_rational>>
multimodular_gauss_elim(sfmpz_mat_t mat,
                        size_t num_threads,
                        bool interreduce,
                        bool proof) {
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
      pivots piv(gauss_elim(nmod_mat.get(), mod, num_threads, trace));
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
    // Initialize crt_entries in reconstruction
    // and clear here
    fmpz* crt_entries;
    idxs.clear();
    rat_entries.clear();

    // the first n_piv rows will be reconstructed
    size_t n_piv;
    if(interreduce)
      n_piv = mat->nrow;
    else
      n_piv = best_piv.size();

    {
      KOMMUNOPP_TIME(crt);
      crt_reconstruction(crt_entries, idxs, good_rrefs, good_primes, n_piv);
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

std::pair<std::vector<std::pair<size_t, size_t>>, std::vector<gmp_rational>>
nmod_gauss_elim(sfmpz_mat_t mat,
                size_t p,
                size_t num_threads,
                bool interreduce) {
  // not needed but use it to avoid code duplication
  bool* trace = new bool[mat->nrow];
  std::fill(trace, trace + mat->nrow, false);

  nmod_t mod;
  nmod_init(&mod, p);

  uint32_mat_t nmod_mat;
  sparse_mat_init(nmod_mat, mat->nrow, mat->ncol);
  mat_mod(nmod_mat, mat, mod, trace);

  KOMMUNOPP_TIME(rref);
  pivots piv(gauss_elim(nmod_mat, mod, num_threads, trace));
  KOMMUNOPP_PROFILE(timer.~adding_timer());

  // compute rows with new leading terms
  size_t n_piv;
  if(interreduce)
    n_piv = mat->nrow;
  else
    n_piv = piv.size();

  // compute nonzero indices & entries
  std::vector<std::pair<size_t, size_t>> idxs;
  std::vector<gmp_rational> entries;
  size_t nnz = sparse_mat_nnz(nmod_mat);
  idxs.reserve(nnz);
  entries.reserve(nnz);
  for(size_t i = 0; i < n_piv; i++) {
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

std::pair<std::vector<std::pair<size_t, size_t>>, std::vector<gmp_rational>>
linear_algebra(sfmpz_mat_t mat,
               size_t characteristic,
               size_t num_threads,
               bool interreduce,
               bool proof) {
  // sort rows by first index and nnz
  std::sort(mat->rows, mat->rows + mat->nrow, [](auto& r1, auto& r2) {
    auto id1 = r1.indices[0];
    auto id2 = r2.indices[0];
    if(id1 != id2)
      return id1 > id2;
    return r1.nnz < r2.nnz;
  });

  if(characteristic == 0)
    return multimodular_gauss_elim(mat, num_threads, interreduce, proof);
  else
    return nmod_gauss_elim(mat, num_threads, characteristic, interreduce);
}
}
