#include <boost/program_options/value_semantic.hpp>
#include <iostream>

#include <boost/program_options.hpp>

#include <config.hpp>

namespace po = boost::program_options;

#include "fast_div.hpp"
#include "profiling.hpp"

using namespace kommunopp;

static void
bench_mod_mersenne() {
  volatile uint32_t res = 0;
  for(uint32_t i = 0; i < 2'147'483'646; ++i) {
    res = v_mod_2_31_1(i);
  }
}

static void
bench_mod_arbitrary(uint32_t p) {
  volatile uint32_t res = 0;
  for(uint32_t i = 0; i < p; ++i) {
    res = v_mod_p(i, p);
  }
}

template<uint32_t p>
static void
bench_mod_tmpl() {
  volatile uint32_t res = 0;
  for(uint32_t i = 0; i < p; ++i) {
    res = v_mod_tmpl<p>(i);
  }
}

static void
bench_mod_standard(uint32_t p) {
  volatile uint32_t res = 0;
  for(uint32_t i = 0; i < p; ++i) {
    res = i % p;
  }
}

static void
bench_mod(uint32_t p) {
  if(p == 2'147'483'647) {
    bench_mod_mersenne();
  } else if(p == 2147483629) {
    bench_mod_tmpl<2147483629>();
  } else {
    bench_mod_arbitrary(p);
  }
}

inline static void
check(uint32_t res, uint32_t i, uint32_t p) {
  if(res != i % p) {
    std::cout << "ERROR " << i << " % " << p << " = " << res
              << " and not the expected " << (i % p) << std::endl;
  }
}

static void
test_mod_mersenne() {
  volatile uint32_t res = 0;
  for(uint32_t i = 0; i < 2'147'483'646; ++i) {
    res = v_mod_2_31_1(i);
    check(res, i, 2'147'483'646);
  }
}

static void
test_mod_arbitrary(uint32_t p) {
  volatile uint32_t res = 0;
  for(uint32_t i = 0; i < p; ++i) {
    res = v_mod_p(i, p);
    check(res, i, 2'147'483'646);
  }
}

static void
test_mod(uint32_t p) {
  if(p == 2'147'483'647) {
    test_mod_mersenne();
  } else {
    test_mod_arbitrary(p);
  }
}

template<typename Functor>
static double
test_function(Functor f, uint32_t p) {
  double wall_time = 0;
  {
    statistics::adding_timer t(wall_time);
    f(p);
  }
  return wall_time;
}

int
main(int argc, char* argv[]) {
  uint32_t prime;

  bool test = false;
  bool bench = false;

  po::options_description desc(
    "mod_prime_tester is a small tool that tests and benchmarks alternative\n"
    "methods of calculating modulo using the nice\n"
    "properties of Mersenne primes. After it runs it prints some statistics\n"
    "that may be analyzed.\n\n"
    "Part of freegb version " FREEGB_VERSION);

  // clang-format off
  desc.add_options()
    ("help,h", "produce help message")
    ("version", "produce version message")
    ("benchmark,b", po::bool_switch(&bench), "run benchmarks")
    ("test,t", po::bool_switch(&test), "run tests")
    ("prime,p", po::value(&prime)->default_value(2'147'483'647), "the prime number to test")
  ;

  po::positional_options_description positional_desc;
  positional_desc.add("prime", 1);
  // clang-format on

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv)
              .options(desc)
              .positional(positional_desc)
              .run(),
            vm);
  po::notify(vm);

  if(vm.count("help")) {
    desc.print(std::cout);
    return EXIT_SUCCESS;
  }

  if(vm.count("version")) {
    std::cerr << FREEGB_VERSION << std::endl;
    return EXIT_SUCCESS;
  }

  if(vm.count("version")) {
    std::cerr << FREEGB_VERSION << std::endl;
    return EXIT_SUCCESS;
  }

  std::cout << std::fixed;

  if(bench) {
    double std = test_function(bench_mod_standard, prime);
    std::cout << "benchmark-standard\t" << prime << "\t" << std << std::endl;
    double opt = test_function(bench_mod, prime);
    std::cout << "benchmark-optimized\t" << prime << "\t" << opt << "\t"
              << (std / opt) << std::endl;
  }
  if(test) {
    std::cout << "test\t\t\t" << prime << "\t" << test_function(test_mod, prime)
              << std::endl;
  }
}
