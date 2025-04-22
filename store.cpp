#include "store.hpp"
#include "f4.hpp"
#include "store.hpp"
#include "parser.hpp"

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
               bool verified_algebra,
               const std::string& output_name,
               const std::string& proof_file,
               bool print_read_problem) {
  return boost::mp11::mp_with_index<F4NCGB_MAX_BLOCKS>(
    nblocks,
    [&context,
     nvars,
     characteristic,
     maxiter,
     maxdeg,
     threads,
     verified_algebra,
     output_name,
     proof_file,
     print_read_problem](auto Nblocks) -> int {
      std::unique_ptr<f4<Nblocks>> algo_ptr
        = std::make_unique<f4<Nblocks>>(context,
                                        nvars,
                                        characteristic,
                                        maxiter,
                                        maxdeg,
                                        threads,
                                        verified_algebra,
                                        proof_file);
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
        msg("Verified algebra:  %s", verified_algebra ? "on" : "off");
        msg("==== Starting Gröbner Basis Computation ====");
      }

      algo.compute_basis();

      msg("==== Basis computation finished ====");

      std::ostream* out_file = &std::cout;
      std::ofstream filestream;
      if(output_name != "") {
        filestream.open(output_name);
        if(!filestream)
          die(18, "Failed to open output file.");
        out_file = &filestream;
      }
      algo.write_basis(*out_file);

    // This is the main function. We do not need to clean up
    // usually. This optimization is only done when not doing
    // Address Sanitizing and when building without assertions.

#if defined(__has_feature) && NDEBUG
#if !__has_feature(address_sanitizer)
      algo_ptr.release();
#endif
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
