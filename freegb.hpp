#ifndef FREEGB_H
#define FREEGB_H

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
    "usage : freeag <mode> <input.aig> <output files> [<option> ...] \n"
    "\n"
    "Depending on the <mode> the <output files> and <options> have to be set:"
    "\n"
    "\n"
    "<mode> = -substitute:\n"
    "    <output files> =  2 output files passed in the following order \n"
    "      <out.cnf>:        CNF miter for correctness of adder substitution \n"
    "      <out.aig>:        rewritten aiger is stored in this file \n"
    "\n"
    "    <option> = the following options are available \n"
    "      -h | --help       print this command line summary \n"
    "      -v<1,2,3,4>       different levels of verbosity(see below) \n"
    "      -signed           option for non-negative integer multipliers \n"
    "\n"
    "\n"
    "<mode> = -verify:\n"
    "    <output files> =  no output files are required \n"
    "     \n "
    "    <option> = the following options are available \n"
    "       -h | --help           print this command line summary \n"
    "       -v<1,2,3,4>           different levels of verbosity(default -v1) \n"
    "       -signed               option for non-negative integer multipliers "
    "\n"
    "       -no-counter-examples  do not generate and write counter examples\n"
    "     \n"
    "     \n"
    "<mode> = -certify:\n"
    "    <output files> =  3 output files passed in the following order\n"
    "      <out.polys>:      initial polynomial set \n"
    "      <out.proof>:      proof rules \n"
    "      <out.spec> :      spec which should be checked \n"
    "     \n "
    "    <option> = the following options are available \n"
    "       -h | --help      print this command line summary \n"
    "       -v<1,2,3,4>      different levels of verbosity(default -v1) \n"
    "       -signed          option for non-negative integer multipliers \n"
    "       -no-counter-examples  do not generate and write counter examples\n"
    "\n"
    "       -p1          expanded proof \n"
    "       -p2          middle condensed proof(some linear combinations "
    "occur, default)\n"
    "       -p3          condensed proof(one single linear combination)\n";

// / Name of the input file
static const char* input_name = 0;

// / \brief
// / Name of output file
static const char* output_name = 0;

// / Selected mode, '-verify' = 1, '-certify' = 2
static int mode;
/*------------------------------------------------------------------------*/
// ERROR CODES:

static int err_no_file = 10;  // no input file given
static int err_mode_sel = 11; // mode has already been selected/not selected
static int err_wrong_arg = 12;// wrong number of arguments given

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
  msg("Copyright(C) 2024 Clemens Hofstadler, JKU Linz, Austria");

  for(int i = 1; i < argc; i++) {
    if(!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) {
      fputs(USAGE, stdout);
      fflush(stdout);
      exit(0);
    } else if(!strcmp(argv[i], "-v0")) {
      verbose = 0;
    } else if(!strcmp(argv[i], "-v1")) {
      verbose = 1;
    } else if(!strcmp(argv[i], "-v2")) {
      verbose = 2;
    } else if(!strcmp(argv[i], "-v3")) {
      verbose = 3;
    } else if(!strcmp(argv[i], "-v4")) {
      verbose = 4;
    } else if(!strcmp(argv[i], "-verify")) {
      if(!mode) {
        msg("selected mode: verification");
        mode = 1;
      } else {
        die(err_mode_sel, "mode has alreday been selected(try '-h')");
      }
    } else if(!strcmp(argv[i], "-certify")) {
      if(!mode) {
        msg("selected mode: verification + certificates");
        mode = 2;
      } else {
        die(err_mode_sel, "mode has alreday been selected(try '-h')");
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

  if(!input_name)
    die(err_no_file, "no input file given(try '-h')");
  // if(!output_name)
  //   die(err_no_file, "no output file given(try '-h')");

  int res = 0;
  init_all_signal_handers();

  std::filesystem::path in = std::filesystem::current_path()
                             / std::filesystem::path("test_inputs/braid3.ms");
  
  parser_context context;
  context.open(in);
  parse_res r = parse_header(context);
  if(r.has_value())
    die(17, "Error in parsing input file.");

  size_t maxiter = 10;
  size_t n = 0;
  if(context.num_blocks > 1)
    n = context.num_blocks;

  boost::mp11::mp_with_index<6>(n, [&maxiter, &context](auto N) {
    f4<N> algo;
    // algo.read_input(context);
    algo.compute_basis(maxiter, 11);
  });

  print_statistics();

  return res;
}

#endif// FREEGB_H
