#include <concepts>
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

#define VERSION "0.1"
/*------------------------------------------------------------------------*/
// / Manual, will be printed with command line '-h'
static const char* USAGE
  = "\n"
    "### USAGE ###\n"
    "usage : freegb <input file> [<output file>] [<option> ...] \n"
    "\n"
    "\n"
    "<input file> = path to an input file in either msolve format or in "
    "(sym)poly format\n"
    "<output file> = path to an output file that will contain the (partial) "
    "Gröbner basis\n"
    "<option> = the following options are available \n"
    "       -h | --help      print this command line summary \n"
    "       -v<1,2,3,4>      different levels of verbosity \n"
    "                        Default: -v1. \n"
    "       -p PRIME         characteristic of the coefficient field.  \n"
    "                        Either 0 (computation over QQ) or a prime < "
    "2^31.\n"
    "                        Default: -p 0. \n"
    "       -m MAXITER       Maximal number of iterations of the F4-algorithm. "
    "\n"
    "                        to be performed. \n"
    "                        Default: -m 10. \n"
    "       -d MAXDEG        Maximal degree of ambiguities that are "
    "considered. \n"
    "                        Default: -d UINT_MAX. \n"
    "       -t  THR          Number of threads to be used. \n"
    "                        Default: -t 1. \n";

// / Name of the input file
static const char* input_name = 0;

// / \brief
// / Name of output file
static const char* output_name = 0;

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

  msg("freegb Version " VERSION);
  msg("");
  msg("Copyright(C) 2024 Clemens Hofstadler, Maximilian Heisinger");
  msg("JKU Linz, Austria");

  for(int i = 1; i < argc; i++) {
    if(!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) {
      fputs(USAGE, stdout);
      fflush(stdout);
      exit(0);
    } else if(!strcmp(argv[i], "-v1")) {
      verbose = 1;
    } else if(!strcmp(argv[i], "-v2")) {
      verbose = 2;
    } else if(!strcmp(argv[i], "-v3")) {
      verbose = 3;
    } else if(!strcmp(argv[i], "-v4")) {
      verbose = 4;
    } else if(!strcmp(argv[i], "-p")) {
      if(prime >= 0) {
        die(err_char_sel,
            "Characteristic has alreday been selected (try '-h')");
      } else if(i + 1 >= argc or !strncmp(argv[i + 1], "-", 1)) {
        die(err_char_sel, "No characteristic provided (try '-h')");
      } else {
        prime = std::atoi(argv[i + 1]);
        i++;
        continue;
      }
    } else if(!strcmp(argv[i], "-m")) {
      if(maxiter > 0) {
        die(err_maxiter_sel, "Maxiter has alreday been selected (try '-h')");
      } else if(i + 1 >= argc or !strncmp(argv[i + 1], "-", 1)) {
        die(err_maxiter_sel, "No maxiter value provided (try '-h')");
      } else {
        maxiter = std::strtoul(argv[i + 1], 0, 10);
        i++;
        continue;
      }
    } else if(!strcmp(argv[i], "-d")) {
      if(maxdeg > 0) {
        die(err_maxdeg_sel, "Maxdeg has alreday been selected (try '-h')");
      } else if(i + 1 >= argc or !strncmp(argv[i + 1], "-", 1)) {
        die(err_maxdeg_sel, "No maxdeg value provided (try '-h')");
      } else {
        maxdeg = std::strtoul(argv[i + 1], 0, 10);
        i++;
        continue;
      }
    } else if(!strcmp(argv[i], "-t")) {
      if(threads > 0) {
        die(err_thread_sel, "Threads has alreday been selected (try '-h')");
      } else if(i + 1 >= argc or !strncmp(argv[i + 1], "-", 1)) {
        die(err_thread_sel, "No number of threads provided (try '-h')");
      } else {
        threads = std::strtoul(argv[i + 1], 0, 10);
        i++;
        continue;
      }
    } else if(output_name) {
      die(err_wrong_arg,
          "too many arguments '%s', '%s', and '%s'(try '-h')",
          input_name,
          output_name,
          argv[i]);
    } else if(input_name) {
      output_name = argv[i];
    } else {
      input_name = argv[i];
    }
  }

  if(!input_name || !std::filesystem::exists(input_name))
    die(err_no_file, "no input file given(try '-h')");

  if(prime < 0)
    prime = 0;
  if(maxiter == 0)
    maxiter = 10;
  if(maxdeg == 0)
    maxdeg = UINT_MAX;
  if(threads == 0)
    threads = 1;

  int res = 0;
  init_all_signal_handers();

  std::filesystem::path in = input_name;

  parser_context context;
  context.open(in);
  parse_res r = context.parse_header();
  if(r.has_value())
    die(17, "Error in parsing input file.");

  size_t n = 0;
  if(context.num_blocks() > 1)
    n = context.num_blocks();

  boost::mp11::mp_with_index<6>(n, [&context](auto N) {
    f4<N, 5> algo((size_t)prime, maxiter, maxdeg, threads);
    if(auto err = algo.read_input(context)) {
      die(17, "Error in parsing body of input file.");
    }

    if(verbose > 1 or true) {
      msg("==== Input Parameters ====");
      if(!output_name)
        msg("No output file specified. Writing output to console.");
      msg("Computing in characteristic %lu.", prime);
      msg("Executing at most %lu iterations.", maxiter);
      msg("Considering ambiguities up to degree %lu.", maxdeg);
      msg("Using %lu threads.", threads);
      msg("==== Starting Gröbner Basis Computation ====");
    }

    algo.compute_basis();
  });

  print_statistics();
  reset_all();

  return res;
}
