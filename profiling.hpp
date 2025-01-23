#pragma once

#include <chrono>
#include <cstddef>

#ifdef KOMMUNOPP_ENABLE_PROFILING
#define KOMMUNOPP_PROFILE(EXPR) EXPR
#define KOMMUNOPP_TIME(TGT) \
  auto timer = ::kommunopp::gstats.time(::kommunopp::gstats.TGT)
#else
#define KOMMUNOPP_PROFILE(EXPR)
#define KOMMUNOPP_TIME(EXPR)
#endif

namespace kommunopp {
#ifdef KOMMUNOPP_ENABLE_PROFILING
struct statistics {
#ifdef KOMMUNOPP_USE_MONOMIAL_PRODUCTS_MAP
  size_t hashmap_hits = 0, hashmap_calls = 0;
#endif

  size_t store_find_calls = 0, store_find_hits = 0;

  // Times
  double init = 0, other = 0, new_elements = 0;
  double amb = 0, crit_pair = 0, sym_pre = 0, reduction = 0;
  double overlap = 0, inclusion = 0, crt = 0, ratrec = 0, rref = 0;
  double binarysearch = 0;
  double elim_task_cpu = 0;

  // Automated timing to add time (in seconds) to some target.
  struct adding_timer {
    adding_timer(double& tgt)
      : target(&tgt) {
      start = std::chrono::high_resolution_clock().now();
    }
    ~adding_timer() {
      if(!target)
        return;

      // The addition is only carried out once! The destructor may be called
      // multiple times.
      std::chrono::time_point<std::chrono::high_resolution_clock> end
        = std::chrono::high_resolution_clock().now();
      std::chrono::duration<double> elapsed = end - start;
      *target += elapsed.count();
      target = nullptr;
    };

    adding_timer(const adding_timer& o) = delete;
    adding_timer& operator=(const adding_timer& o) = delete;

    private:
    std::chrono::time_point<std::chrono::high_resolution_clock> start;
    double* target;
  };

  adding_timer time(double& tgt) { return adding_timer(tgt); }

  void print(size_t threads = 1);
};
extern statistics gstats;
#endif
}
