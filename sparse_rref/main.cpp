// // use mimalloc to replace the default malloc
// #include <cstdlib>
// #include "mimalloc-override.h"
// #include "mimalloc-new-delete.h"

#include <filesystem>
#include <fstream>
#include <iostream>

#include "argparse.hpp"
#include "sparse_mat.h"
#include "sparse_vec.h"

#define printtime(str)                                                  \
  std::cout << (str) << " spent " << std::fixed << std::setprecision(6) \
            << sparse_base::usedtime(start, end) << " seconds." << std::endl

#define printmatinfo(mat)                             \
  std::cout << "nnz: " << sparse_mat_nnz(mat) << " "; \
  std::cout << "nrow: " << (mat)->nrow << " ";        \
  std::cout << "ncol: " << (mat)->ncol << std::endl

int
main(int argc, char** argv) {
  std::string version = "v0.2.1";
  argparse::ArgumentParser program("sparserref", version);
  program.set_usage_max_line_width(80);
  program.add_description("(exact) Sparse Reduced Row Echelon Form " + version);
  program.add_argument("input_file").help("input file in matrix market format");
  program.add_argument("-o", "--output")
    .help("output file in matrix market format")
    .default_value("input_file.rref")
    .nargs(1);
  program.add_usage_newline();
  program.add_argument("-k", "--kernel")
    .default_value(false)
    .help("output the kernel")
    .implicit_value(true)
    .nargs(0);
  program.add_argument("--output-pivots")
    .help("output pivots")
    .default_value(false)
    .implicit_value(true)
    .nargs(0);
  program.add_usage_newline();
  program.add_argument("-F", "--field")
    .default_value("QQ")
    .help("QQ: rational field\nZp or Fp: Z/p for a prime p")
    .nargs(1);
  program.add_argument("-p", "--prime")
    .default_value("34534567")
    .help("a prime number, only vaild when field is Zp ")
    .nargs(1);
  program.add_argument("-t", "--threads")
    .help("the number of threads ")
    .default_value(1)
    .nargs(1)
    .scan<'i', int>();
  program.add_argument("-pd", "--pivot_direction")
    .help("the direction to select pivots")
    .default_value("row")
    .nargs(1);
  program.add_argument("-sd", "--search_depth")
    .help("the depth of search, default is the max of int ")
    .default_value(0)
    .nargs(1)
    .scan<'i', int>();
  program.add_usage_newline();
  program.add_argument("-V", "--verbose")
    .default_value(false)
    .help("prints information of calculation")
    .implicit_value(true)
    .nargs(0);
  program.add_argument("-ps", "--print_step")
    .default_value(100)
    .help("print step when --verbose is enabled")
    .nargs(1)
    .scan<'i', int>();
  program.add_usage_newline();
  program.add_argument("--no-backward-substitution")
    .help("no backward substitution")
    .default_value(false)
    .implicit_value(true)
    .nargs(0);

  try {
    program.parse_args(argc, argv);
  } catch(const std::exception& err) {
    std::cerr << err.what() << std::endl;
    std::cerr << program;
    std::exit(1);
  }

  ulong prime;
  nmod_t p;
  if(program.get<std::string>("--field") == "QQ") {
    prime = 0;
  } else {
    prime = 65521;
    nmod_init(&p, prime);
    std::cout << "Using prime: " << prime << std::endl;
  } 

  int nthread = 4;
  BS::thread_pool pool(nthread);

  field_t F;
  if(prime == 0)
    field_init(F, FIELD_QQ, 1, NULL);
  else
    field_init(F, FIELD_Fp, std::vector<ulong>{ prime });

  sparse_mat_t<fmpq> mat_Q;
  sparse_mat_t<ulong> mat_Zp;

  auto start = sparse_base::clocknow();
  auto input_file = program.get<std::string>("input_file");
  std::filesystem::path filePath = input_file;
  if(!std::filesystem::exists(filePath)) {
    std::cerr << "File does not exist: " << filePath << std::endl;
    return 1;
  }

  std::ifstream file(filePath);
  sfmpq_mat_read(mat_Q, file);
  if(prime != 0) {
    sparse_mat_init(mat_Zp, mat_Q->nrow, mat_Q->ncol);
    snmod_mat_from_sfmpq(mat_Zp, mat_Q, p);
    sparse_mat_clear(mat_Q);
  }
  file.close();

  auto end = sparse_base::clocknow();
  std::cout << "-------------------" << std::endl;
  printtime("read");

  if(prime == 0) {
    printmatinfo(mat_Q);
  } else {
    printmatinfo(mat_Zp);
  }

  std::cout << "-------------------" << std::endl;
  std::cout << "RREFing: " << std::endl;

  rref_option_t opt;
  opt->verbose = (program["--verbose"] == true);
  opt->is_back_sub = (program["--no-backward-substitution"] == false);
  opt->print_step = program.get<int>("--print_step");
  opt->search_depth = (ulong)program.get<int>("--search_depth");
  opt->pivot_dir = (program.get<std::string>("--pivot_direction") == "row");
  if(opt->search_depth == 0)
    opt->search_depth = INT_MAX;

  start = sparse_base::clocknow();
  std::vector<std::pair<slong, slong>> pivots;
  if(prime == 0) {
    pivots = sparse_mat_rref(mat_Q, F, pool, opt);
  } else {
    pivots = sparse_mat_rref(mat_Zp, F, pool, opt);
  }

  end = sparse_base::clocknow();
  std::cout << "-------------------" << std::endl;
  printtime("RREF");

  std::cout << "rank: " << pivots.size() << " ";
  if(prime == 0) {
    printmatinfo(mat_Q);
  } else {
    printmatinfo(mat_Zp);
  }

  start = sparse_base::clocknow();
  std::ofstream file2;
  std::string outname, outname_add("");
  if(program.get<std::string>("--output") == "input_file.rref")
    outname = input_file;
  else
    outname = program.get<std::string>("--output");

  if(program["--output-pivots"] == true) {
    outname_add = ".piv";
    file2.open(outname + outname_add);
    for(auto& ii : pivots)
      file2 << ii.first + 1 << ", " << ii.second + 1 << '\n';
    file2.close();
  }

  if(outname == input_file)
    outname_add = ".rref";
  else
    outname_add = "";

  // // print matrix
  // for(size_t i = 0; i < mat_Zp->nrow; i++) {
  //   auto row = sparse_mat_row(mat_Zp, i);
  //   for(size_t j = 0; j < 10; j++) {
  //     auto c = sparse_vec_entry(row, j);
  //     std::cout << c << ", ";
  //   }
  //   std::cout << "\n";
  // }

  // file2.open(outname + outname_add);
  // // file2.open(outname + outname_add, std::ios::binary);
  // if (prime == 0) {
  // 	sparse_mat_write(mat_Q, file2);
  // }
  // else {
  // 	sparse_mat_write(mat_Zp, file2);
  // 	// auto buffer = snmod_mat_to_binary(mat_Zp);
  // 	// file2.write(buffer.second, buffer.first * sizeof(char));
  // 	// s_free(buffer.second);
  // }
  // file2.close();

  end = sparse_base::clocknow();
  printtime("write files");

  field_clear(F);

  // clean is very expansive, leave to OS :(
  sparse_mat_clear(mat_Q);
  sparse_mat_clear(mat_Zp);
  return 0;
}
