#include "kommunopp.hpp"
#include "f4.hpp"
#include "freegb.hpp"
#include "parser.hpp"

#include <boost/multiprecision/detail/default_ops.hpp>
#include <boost/multiprecision/gmp.hpp>

extern bool tracer;

namespace kommunopp {
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
kommunopp_main(parser_context& context,
               size_t nblocks,
               size_t nvars,
               size_t characteristic,
               size_t maxiter,
               size_t maxdeg,
               size_t threads,
               bool verified_algebra,
               const std::string& output_name,
               bool print_read_problem) {
  boost::mp11::mp_with_index<MAX_BLOCKS>(
    nblocks,
    [&context,
     nvars,
     characteristic,
     maxiter,
     maxdeg,
     threads,
     verified_algebra,
     output_name,
     print_read_problem](auto Nblocks) {
      f4<Nblocks> algo(
        nvars, characteristic, maxiter, maxdeg, threads, verified_algebra);
      {
        KOMMUNOPP_TIME(parse);
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
        std::string out = "Output file:       ";
        if(output_name == "")
          out += "None. Writing output to console.";
        else
          out += output_name;

        msg(out.c_str());
        msg("Characteristic:    %lu", characteristic);
        msg("Max. Iterations:   %lu", maxiter);
        msg("Max. amb. degree:  %lu", maxdeg);
        msg("Nr. threads:       %lu", threads);
        msg("Proof logging:     %d",  proof);
        msg("Tracer:            %d",  tracer);
        msg("Verified algebra:  %d",  verified_algebra);
        msg("==== Starting Gröbner Basis Computation ====");
      }

      algo.compute_basis();
      // msg("Success");
    });
  return 0;
}
}

std::ostream&
operator<<(std::ostream& o,
           const boost::multiprecision::backends::gmp_rational& r) {
  using namespace boost::multiprecision;
  mpq_rational rr(r);
  return o << rr;
}
