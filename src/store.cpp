#include "store.hpp"
#include "f4.hpp"
#include "f4ncgb.hpp"
#include "parser.hpp"
#include "store.hpp"

#include <boost/multiprecision/detail/default_ops.hpp>
#include <boost/multiprecision/gmp.hpp>
#include <fstream>

extern bool tracer;

namespace f4ncgb {
int
coefficient_sign(const boost::multiprecision::backends::gmp_rational& r) {
  using namespace boost::multiprecision;
  mpq_rational rr(r);
  return boost::multiprecision::sign(rr);
}

boost::multiprecision::mpq_rational
coefficient_abs(const boost::multiprecision::backends::gmp_rational& r) {
  using namespace boost::multiprecision;
  mpq_rational rr(r);
  return boost::multiprecision::abs(rr);
}

bool
coefficient_is_posneg_neutral(
  const boost::multiprecision::backends::gmp_rational& r) {
  using namespace boost::multiprecision;
  mpq_rational rr(r);
  return (numerator(rr) == -1l || numerator(rr) == 1l) && denominator(rr) == 1l;
}

int
f4ncgb_main(parser_context& context,
            size_t nblocks,
            size_t nvars,
            size_t characteristic,
            size_t maxiter,
            size_t maxdeg,
            size_t threads,
            const std::string& output_name,
            const std::string& proof_file,
            std::function<void()>* stats_print_function,
            bool leak_memory,
            bool print_read_problem,
            void* userdata,
            f4ncgb_add_cb add_cb,
            f4ncgb_end_poly_cb end_poly_cb) {
  return boost::mp11::mp_with_index<F4NCGB_MAX_BLOCKS>(
    nblocks,
    [&context,
     nvars,
     characteristic,
     maxiter,
     maxdeg,
     threads,
     output_name,
     proof_file,
     stats_print_function,
     leak_memory,
     print_read_problem,
     userdata,
     add_cb,
     end_poly_cb](auto Nblocks) -> int {
      std::unique_ptr<f4<Nblocks>> algo_ptr = std::make_unique<f4<Nblocks>>(
        context, nvars, characteristic, maxiter, maxdeg, threads, proof_file);

      auto& algo = *algo_ptr;
      {
        F4NCGB_TIME(parse);
        if(auto err = algo.read_input(context)) {
          std::cerr << *err << std::endl;
          die(17, "Error in parsing body of input file.");
        }
      }

      if(print_read_problem) {
        context.to_msolve(std::cout, algo.poly);
      }

#ifdef F4NCGB_ENABLE_STORE_DUMP
      algo.mons.set_binary_dump_path(
        std::filesystem::path(context.get_filename()).filename().string()
        + "_mons.bin");
      algo.poly.set_binary_dump_path(
        std::filesystem::path(context.get_filename()).filename().string()
        + "_poly.bin");
#endif

      if(verbose > 0) {
        msg("==== Input Parameters ====");
        std::string out_name = "Output file:       ";
        if(output_name == "")
          out_name += "None. Writing output to console.";
        else
          out_name += output_name;

        std::ostringstream mon_order;
        mon_order << "Monomial order:    ";
        bool first = true;
        if(Nblocks == 0) {
          for(size_t v = 1; v <= nvars; v++) {
            if(!first)
              mon_order << " < ";
            context.var_to_ostream(mon_order, v);
            first = false;
          }
        } else {
          for(size_t i = 0; i < context.num_blocks(); i++) {
            first = true;
            for(auto v : context.block(i)) {
              if(!first)
                mon_order << " < ";
              context.var_to_ostream(mon_order, v);
              first = false;
            }
            if(i < context.num_blocks() - 1)
              mon_order << " << ";
          }
        }

        std::string proof_str = proof ? "on " : "off";
        if(proof > 1)
          proof_str += "(expanded) ";
        if(proof > 0)
          proof_str += "(writing to " + proof_file + ")";

        msg("Characteristic:    %lu", characteristic);
        msg("Max. Iterations:   %lu", maxiter);
        msg("Max. amb. degree:  %lu", maxdeg);
        // msg("Number of blocks:  %lu", nblocks);
        // msg("Max nr. of blocks: %lu", F4NCGB_MAX_BLOCKS);
        msg(mon_order.str().c_str());
        msg("Nr. threads:       %lu", threads);
        msg(out_name.c_str());
        msg("Proof logging:     %s", proof_str.c_str());
        msg("Tracer:            %s", tracer ? "on" : "off");
        msg("PID of F4NCGB:     %lu", getpid());
        msg("==== Starting Gröbner Basis Computation ====");
      }

      if(stats_print_function) {
        *stats_print_function = [&algo]() {
          msg("Current basis:");
          algo.write_basis(std::cerr);
        };
      }

      algo.compute_basis();

      msg("==== Basis computation finished ====");

      if(add_cb && end_poly_cb) {
        algo.write_basis(userdata, add_cb, end_poly_cb);
      } else {
        std::ostream* out_file = &std::cout;
        std::ofstream filestream;
        if(output_name != "") {
          filestream.open(output_name, std::ios_base::trunc);
          if(!filestream)
            die(18, "Failed to open output file.");
          out_file = &filestream;
        }
        algo.write_basis(*out_file);
      }

      if(stats_print_function) {
        *stats_print_function = []() {};
      }

    // This is the main function. We do not need to clean up
    // usually. This optimization is only done when not doing
    // Address Sanitizing and when building without assertions.
    // Store dumping also disables this optimization.

#if defined(__has_feature) && NDEBUG && !defined(F4NCGB_ENABLE_STORE_DUMP)
#if !__has_feature(address_sanitizer)
      if(leak_memory)
        algo_ptr.release();
#else
      (void)leak_memory;
#endif
#else
      (void)leak_memory;
#endif
      return 0;
    });
}
}

std::ostream&
operator<<(std::ostream& o,
           const boost::multiprecision::backends::gmp_rational& r) {
  using namespace boost::multiprecision;
  mpq_rational rr(r);
  return o << rr;
}
