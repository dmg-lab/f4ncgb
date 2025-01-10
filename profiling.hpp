#pragma once

#include <chrono>
#include <cstddef>

#ifdef KOMUNOPP_ENABLE_PROFILING
#define KOMMUNOPP_PROFILE(EXPR) EXPR
#else
#define KOMMUNOPP_PROFILE(EXPR)
#endif

namespace kommunopp {
#ifdef KOMMUNOPP_ENABLE_PROFILING
struct statistics {
  size_t hashmap_hits = 0, hashmap_calls = 0;
  size_t monomial_hits = 0, monomial_calls = 0;

  // Times
  double init = 0, other = 0, new_elements = 0;
  double amb = 0, crit_pair = 0, sym_pre = 0, reduction = 0;
  double overlap = 0, inclusion = 0, crt = 0, ratrec = 0, rref = 0;

  struct timer {
    timer() = default;
    ~timer() = default;

    std::chrono::high_resolution_clock start
      (std::chrono::high_resolution_clock().now());
  }
};
extern statistics gstats;
#endif
}
