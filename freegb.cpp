#include <boost/program_options/options_description.hpp>
#include <boost/program_options/positional_options.hpp>
#include <boost/program_options/variables_map.hpp>
#include <cstdlib>
#include <filesystem>
#include <string.h>

#include "freegb.hpp"
#include "parser.hpp"
#include "profiling.hpp"
#include "signal_statistics.hpp"
#include <flint/ulong_extras.h>

#include <boost/program_options.hpp>
#include <config.hpp>

namespace po = boost::program_options;

// / Name of the input file
static std::string input_name = "";

// / \brief
// / Name of output files
static std::string output_name = "";
static std::string proof_file = "";

// / Selected prime, maxiter, maxdeg, and number of threads
static size_t maxiter = 0;
static size_t maxdeg = 0;
static size_t threads = 0;

// / Verified multimodular computations and proof logging
static bool verified_algebra = true;
bool proof = false;
bool tracer = true;
/*------------------------------------------------------------------------*/
// ERROR CODES:

static int err_no_file = 10;// no input file given
static int err_prime_too_big = 11;
static int err_no_prime = 12;

static bool print_problem = false;

using namespace kommunopp;

/**
    Main Function of freegb.

    Prints statistics to stdout after finishing.
*/
int
main(int argc, char** argv) {
  // clang-format off
  po::options_description desc("freegb, version " FREEGB_VERSION "\n"
                               "Copyright(C) 2025 Clemens Hofstadler, Maximilian Heisinger\n"
                               "JKU Linz, Austria\n"
                               "USAGE");
  desc.add_options()
    ("version", "produce version message")
    ("help,h", "produce help message")
    ("input,i", po::value<std::string>(&input_name)->default_value(""), "Set the input file (either msolve, poly, or sympoly; switched according to content).")
    ("maxiter,m", po::value<size_t>(&maxiter)->default_value(10), "Maximal number of iterations of the F4-algorithm to be performed.")
    ("maxdeg,d", po::value<size_t>(&maxdeg)->default_value(UINT_MAX), "Maximal degree of ambiguities that are considered.")
    ("output,o", po::value<std::string>(&output_name)->default_value(""), "Set the output file.")
    ("print-problem", po::value<bool>(&print_problem)->default_value(false), "Re-print the problem after parsing.")
    ("proof,p",  po::value<std::string>(&proof_file)->default_value(""), "Proof logging.")
    ("threads,T", po::value<size_t>(&threads)->default_value(1), "Number of threads to be used.")
    ("tracer,t", po::value<bool>(&tracer)->default_value(true), "Whether computations with the first prime shall be traced. Speeds up the computation, but yields the correct result only with high probability.")
    ("verbosity,V", po::value<int>(&verbose)->default_value(1), "Set the verbosity level.")
    ("verify-algebra,v", po::value<bool>(&verified_algebra)->default_value(true), "Whether the multimodular linear algebra shall be verified (otherwise, the result is only correct with high probability).")
  ;

  po::positional_options_description positional_desc;
  positional_desc.add("input", 1);
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
    std::cout << FREEGB_VERSION << std::endl;
    return EXIT_SUCCESS;
  }

  if(!std::filesystem::exists(input_name))
    die(err_no_file, "no input file given(try '-h')");

  int res = 0;

  std::filesystem::path in = input_name;

  parser_context context;
  parse_res r;
  {
    KOMMUNOPP_TIME(parse);
    context.open(in);
    r = context.parse_header();
  }
  if(r.has_value()) {
    std::cerr << *r << std::endl;
    die(17, "Error in parsing input file.");
  }

  size_t nvars = context.num_vars();
  size_t nblocks = 0;
  if(context.num_blocks() > 1)
    nblocks = context.num_blocks();
  if(nblocks > KOMMUNOPP_MAX_BLOCKS)
    die(4, "More blocks than current compilation allows\n");

  size_t characteristic = context.characteristic();
  if(characteristic > 2147483647l)// 2^31 -1
    die(err_prime_too_big,
        "Provided characteristic %lu is too large. Only p < 2^31 supported",
        characteristic);
  if(characteristic != 0 and !n_is_prime(characteristic))
    die(err_no_prime,
        "Provided nonzero characteristic %lu is not prime.",
        characteristic);

  proof = (proof_file != "");

  res = kommunopp::kommunopp_main(context,
                                  nblocks,
                                  nvars,
                                  characteristic,
                                  maxiter,
                                  maxdeg,
                                  threads,
                                  verified_algebra,
                                  output_name,
                                  proof_file,
                                  print_problem);

  KOMMUNOPP_PROFILE(gstats.print(threads));
  return res;
}
