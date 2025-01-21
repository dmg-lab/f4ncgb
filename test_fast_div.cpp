#include <boost/test/unit_test.hpp>

#include "fast_div.hpp"
#include "profiler.h"
#include "profiling.hpp"

using namespace kommunopp;

#define TEST(v) BOOST_REQUIRE(v_mod_2_31_1(v) == (v % 2'147'483'647))

BOOST_AUTO_TEST_CASE(test_v_mod_2_31_1) {
  TEST(10);
  TEST(12313);
  TEST(121);
  TEST(0);
  TEST(1111111);
  TEST(123211);
  TEST(12124124);
  TEST(12124124);
  TEST(1212124124);
  TEST(2'147'483'646);

#ifdef KOMMUNOPP_ENABLE_PROFILING
  volatile size_t N = 4'000'000;
  double mod_time = 0;
  double fast_mod_time = 0;

  volatile uint32_t val;
  volatile uint32_t res;

  BOOST_TEST_MESSAGE("Profile fast mod using mersenne prime.");
  {
    auto timer = statistics::adding_timer(fast_mod_time);
    for(size_t i = 0; i < N; ++i) {
      val = static_cast<uint32_t>(N);
      res = v_mod_2_31_1(val);
    }
  }

  BOOST_TEST_MESSAGE("Profile regular mod using %");
  {
    auto timer = statistics::adding_timer(mod_time);
    for(size_t i = 0; i < N; ++i) {
      val = static_cast<uint32_t>(N);
      res = val % 2'147'483'647;
    }
  }
  BOOST_REQUIRE(res != 0);

  double speedup = mod_time / fast_mod_time;

  BOOST_TEST_MESSAGE("N=" << N << ", mod time: " << mod_time
                          << ", fast mod time: " << fast_mod_time
                          << ", speedup: " << speedup);
#endif
}

#undef TEST

#define TEST(a, b, x)        \
  BOOST_REQUIRE(             \
    mult_mod_2_31_1(a, b, x) \
    == (((uint64_t)a * (uint64_t)x + (uint64_t)b) % 2'147'483'647))

BOOST_AUTO_TEST_CASE(test_mult_mod_2_31_1) {
  TEST(12, 4, 55);
  TEST(252, 5235, 55);
  TEST(1, 0, 2425);
  TEST(0, 0, 2425);
  TEST(425, 635, 2425);

#ifdef KOMMUNOPP_ENABLE_PROFILING
  volatile size_t N = 4'000'000;
  double mod_time = 0;
  double fast_mod_time = 0;

  volatile uint32_t val;
  volatile uint32_t res;

  BOOST_TEST_MESSAGE("Profile fast mult mod using mersenne prime.");
  {
    auto timer = statistics::adding_timer(fast_mod_time);
    for(size_t i = 0; i < N; ++i) {
      val = static_cast<uint32_t>(N);
      res = mult_mod_2_31_1(val, val + 1, val + 2);
    }
  }

  BOOST_TEST_MESSAGE("Profile regular mult mod using %");
  {
    auto timer = statistics::adding_timer(mod_time);
    for(size_t i = 0; i < N; ++i) {
      val = static_cast<uint32_t>(N);
      res = ((uint64_t)val * ((uint64_t)val + 2) + ((uint64_t)val + 1))
            % 2'147'483'647;
    }
  }
  BOOST_REQUIRE(res != 0);

  double speedup = mod_time / fast_mod_time;

  BOOST_TEST_MESSAGE("N=" << N << ", mod time: " << mod_time
                          << ", fast mod time: " << fast_mod_time
                          << ", speedup: " << speedup);
#endif
}

#undef TEST
