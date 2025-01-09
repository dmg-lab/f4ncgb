#include <boost/program_options/options_description.hpp>
#include <boost/program_options/positional_options.hpp>
#include <concepts>
#include <cstdlib>
#include <filesystem>
#include <memory>
#include <set>
#include <sstream>
#include <string.h>
#include <vector>

#include "f4.hpp"
#include "kommunopp.hpp"
#include "parser.hpp"
#include "signal_statistics.hpp"

#include <boost/program_options.hpp>
#include <config.hpp>

namespace po = boost::program_options;


static const size_t MAX_VARS = 10;
static const size_t MAX_BLOCKS = 3;

// / Name of the input file
static std::string input_name = "";

// / \brief
// / Name of output file
static std::string output_name = "";

// / Selected prime, maxiter, maxdeg, and number of threads
static int64_t prime = -1;
static size_t maxiter = 0;
static size_t maxdeg = 0;
static size_t threads = 0;
/*------------------------------------------------------------------------*/
// ERROR CODES:

static int err_no_file = 10;// no input file given
static int err_char_sel
  = 11;// characteristic has already been selected/not selected
static int err_maxiter_sel
  = 12;                        // maxiter has already been selected/not selected
static int err_maxdeg_sel = 13;// maxdeg has already been selected/not selected
static int err_thread_sel
  = 14;// number of threads has already been selected/not selected
static int err_wrong_arg = 15;// wrong number of arguments given

/*------------------------------------------------------------------------*/
/**
    Calls the deallocaters of the involved data types
    @see reset_all_signal_handlers()
    @see delete_gates()
    @see deallocate_terms()
    @see deallocate_mstack()
    @see clear_mpz()
*/
static void
reset_all() {
  reset_all_signal_handlers();
  //   delete_gates();
  //   deallocate_terms();
  //   deallocate_mstack();
  //   clear_mpz();
}
/*------------------------------------------------------------------------*/

using namespace kommunopp;

/**
    Main Function of freegb.

    Prints statistics to stdout after finishing.
*/
int
main(int argc, char** argv) {
  // clang-format off
  po::options_description desc("freegb, version " FREEGB_VERSION "\n"
                               "Copyright(C) 2024 Clemens Hofstadler, Maximilian Heisinger\n"
                               "JKU Linz, Austria\n"
                               "USAGE");
  desc.add_options()
    ("help,h", "produce help message")
    ("version", "produce version message")
    ("input,i", po::value<std::string>(&input_name)->default_value(""), "set the input file (either msolve, poly, or sympoly; switched according to content)")
    ("output,o", po::value<std::string>(&output_name)->default_value(""), "set the output file")
    ("verbosity,v", po::value<int>(&verbose)->default_value(1), "set the verbosity level")
    ("prime,p", po::value<int64_t>(&prime)->default_value(0), "characteristic of the coefficient field. Either 0 (computation over QQ) or a prime < 2^31")
    ("maxiter,m", po::value<size_t>(&maxiter)->default_value(10), "Maximal number of iterations of the F4-algorithm to be performed.")
    ("maxdeg,d", po::value<size_t>(&maxdeg)->default_value(UINT_MAX), "Maximal degree of ambiguities that are considered.")
    ("threads,t", po::value<size_t>(&threads)->default_value(1), "Number of threads to be used.")
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
    std::cerr << FREEGB_VERSION << std::endl;
    return EXIT_SUCCESS;
  }

  if(!std::filesystem::exists(input_name))
    die(err_no_file, "no input file given(try '-h')");

  int res = 0;
  init_all_signal_handers();

  std::filesystem::path in = input_name;

  parser_context context;
  context.open(in);
  parse_res r = context.parse_header();
  if(r.has_value()) {
    std::cerr << *r << std::endl;
    die(17, "Error in parsing input file.");
  }

  size_t nvars = context.num_vars();
  size_t nblocks = 0;
  if(context.num_blocks() > 1)
    nblocks = context.num_blocks();

  if(nvars > MAX_VARS)
    die(4, "More variables than current compilation allows\n");
  if(nblocks > MAX_BLOCKS)
    die(4, "More blocks than current compilation allows\n");

  boost::mp11::mp_with_index<MAX_VARS>(nvars, [&context,nblocks] (auto Nvars) {
    boost::mp11::mp_with_index<MAX_BLOCKS>(nblocks, [&context, Nvars](auto Nblocks) {
    f4<Nvars, Nblocks> algo((size_t)prime, maxiter, maxdeg, threads);
    if(auto err = algo.read_input(context)) {
      std::cerr << *err << std::endl;
      die(17, "Error in parsing body of input file.");
    }

    if(verbose > 1 or true) {
      msg("==== Input Parameters ====");
      if(output_name == "")
        msg("No output file specified. Writing output to console.");
      msg("Computing in characteristic %lu.", prime);
      msg("Executing at most %lu iterations.", maxiter);
      msg("Considering ambiguities up to degree %lu.", maxdeg);
      msg("Using %lu threads.", threads);
      msg("==== Starting Gröbner Basis Computation ====");
    }

    algo.compute_basis();
  });});

  print_statistics();
  reset_all();

  return res;
}
