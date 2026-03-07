#pragma once

#include <cstddef>
#include <vector>

#include <flint/fmpz.h>
#include <flint/nmod.h>

#include "coeff.hpp"
#include "sparse_rref/sparse_mat.h"
#include "sparse_rref/sparse_vec.h"
#include "sparse_rref/thread_pool.hpp"

namespace f4ncgb {

typedef sparse_vec_t<uint32_t> uint32_vec_t;
typedef sparse_mat_t<uint32_t> uint32_mat_t;
using pivots = std::vector<size_t>;


std::pair<std::vector<std::pair<size_t, size_t>>, fmpz*>
linear_algebra(std::vector<std::vector<size_t>>& idxs_in,
               std::vector<std::span<coeff>>& entries_in,
               size_t characteristic,
               std::unique_ptr<BS::thread_pool<BS::none>>& pool,
               bool tracer = true);
}
