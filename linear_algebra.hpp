#pragma once

#include <cstddef>
#include <vector>

#include "gmp.h"
#include <boost/multiprecision/gmp.hpp>
#include <flint/nmod.h>

#include "sparse_rref/sparse_mat.h"
#include "sparse_rref/sparse_vec.h"
#include "sparse_rref/thread_pool.hpp"

namespace kommunopp {

typedef sparse_vec_t<fmpz> sfmpz_vec_t;
typedef sparse_mat_t<fmpz> sfmpz_mat_t;
typedef sparse_vec_t<uint32_t> uint32_vec_t;
typedef sparse_mat_t<uint32_t> uint32_mat_t;
using pivots = std::vector<size_t>;

pivots
reverse_solve(uint32_mat_t mat, nmod_t mod);

pivots
gauss_elim(uint32_mat_t mat, nmod_t mod, size_t num_threads, bool* trace);

std::pair<std::vector<std::pair<size_t, size_t>>,
          std::vector<boost::multiprecision::gmp_rational>>
multimodular_gauss_elim(sfmpz_mat_t mat,
                        std::unique_ptr<BS::thread_pool<BS::none>>& pool,
                        bool interreduce = false,
                        bool verify = true);

std::pair<std::vector<std::pair<size_t, size_t>>,
          std::vector<boost::multiprecision::gmp_rational>>
nmod_gauss_elim(sfmpz_mat_t mat,
                size_t p,
                std::unique_ptr<BS::thread_pool<BS::none>>& pool,
                bool interreduce);

std::pair<std::vector<std::pair<size_t, size_t>>,
          std::vector<boost::multiprecision::gmp_rational>>
linear_algebra(sfmpz_mat_t mat,
               size_t characteristic,
               std::unique_ptr<BS::thread_pool<BS::none>>& pool,
               bool interreduce = false,
               bool verify = true);
}
