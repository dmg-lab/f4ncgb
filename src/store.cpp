#include "store.hpp"
#include "f4.hpp"
#include "f4ncgb.hpp"
#include "parser.hpp"
#include "store.hpp"

#include <boost/multiprecision/detail/default_ops.hpp>
#include <boost/multiprecision/gmp.hpp>
#include <fstream>

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

struct f4_base {
  virtual ~f4_base() = default;
  virtual parse_res read_input(parser_context& context) = 0;
  virtual void compute_basis() = 0;
  virtual void reduced_form() = 0;
  virtual void write_basis(std::ostream& os) = 0;
  virtual void write_basis(void* userdata,
                           f4ncgb_add_cb add_cb,
                           f4ncgb_end_poly_cb end_poly_cb)
    = 0;
  virtual void print_read_problem(parser_context& context) = 0;
  virtual void print_parameters(parser_context& context,
                                const std::string& output_name)
    = 0;

#ifdef F4NCGB_ENABLE_STORE_DUMP
  virtual void set_dump_paths(std::string path) = 0;
#endif
};

template<std::size_t Nblocks, typename V>
struct f4_wrapper : f4_base {
  f4<Nblocks, V> algo;

  f4_wrapper(parser_context& context)
    : algo(context) {}

  parse_res read_input(parser_context& context) override {
    return algo.read_input(context);
  }
  void compute_basis() override { algo.compute_basis(); }
  void reduced_form() override { algo.reduced_form(); }
  void write_basis(std::ostream& os) override { algo.write_basis(os); }
  void write_basis(void* userdata,
                   f4ncgb_add_cb add_cb,
                   f4ncgb_end_poly_cb end_poly_cb) override {
    algo.write_basis(userdata, add_cb, end_poly_cb);
  }
  void print_read_problem(parser_context& context) override {
    context.to_msolve(std::cout, algo.poly);
  }

  void print_parameters(parser_context& context,
                        const std::string& output_name) override {
    msg("==== Input Parameters ====");
    std::string out_name = "Output file:       ";
    if(output_name == "")
      out_name += "None. Writing output to console.";
    else
      out_name += output_name;

    std::ostringstream mon_order;
    mon_order << "Monomial order:    ";
    bool first = true;
    if(context.num_blocks() == 0) {
      for(size_t v = 1; v <= context.num_vars(); v++) {
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

    std::string proof_str = context.proof_level() > 0 ? "on " : "off";
    if(context.proof_level() > 1)
      proof_str += "(expanded) ";
    if(context.proof_level() > 0)
      proof_str += "(writing to " + context.proof_file() + ")";

    msg("Characteristic:    %lu", context.characteristic());
    msg("Max. Iterations:   %lu", context.maxiter());
    msg("Max. amb. degree:  %lu", context.maxdeg());
    msg(mon_order.str().c_str());
    msg("Nr. threads:       %lu", context.threads());
    msg(out_name.c_str());
    msg("Proof logging:     %s", proof_str.c_str());
    msg("Tracer:            %s", context.tracer() ? "on" : "off");
    msg("PID of F4NCGB:     %lu", getpid());
    msg("==== Starting Gröbner Basis Computation ====");
  }

#ifdef F4NCGB_ENABLE_STORE_DUMP
  void set_dump_paths(std::string path) override {
    algo.mons.set_binary_dump_path(std::filesystem::path(path + "_mons.bin"));
    algo.poly.set_binary_dump_path(std::filesystem::path(path + "_poly.bin"));
  }
#endif
};

int
f4ncgb_main(parser_context& context,
            const std::string& output_name,
            std::function<void()>* stats_print_function,
            bool leak_memory,
            bool print_read_problem,
            void* userdata,
            f4ncgb_add_cb add_cb,
            f4ncgb_end_poly_cb end_poly_cb) {

  return boost::mp11::mp_with_index<F4NCGB_MAX_BLOCKS>(
    context.num_blocks(),
    [&context,
     output_name,
     stats_print_function,
     leak_memory,
     print_read_problem,
     userdata,
     add_cb,
     end_poly_cb](auto Nblocks) -> int {
      std::unique_ptr<f4_base> algo_ptr;
      size_t nvars = context.num_vars();

      if(nvars < std::numeric_limits<uint8_t>::max()) {
        algo_ptr = std::make_unique<f4_wrapper<Nblocks, uint8_t>>(context);
      } else if(nvars < std::numeric_limits<uint16_t>::max()) {
        algo_ptr = std::make_unique<f4_wrapper<Nblocks, uint16_t>>(context);
      } else {
        die(31, "Too many variables. At most 2^16-1 supported.");
      }

      auto& algo = *algo_ptr;

      {
        F4NCGB_TIME(parse);
        if(auto err = algo.read_input(context)) {
          std::cerr << *err << std::endl;
          die(17, "Error in parsing body of input file.");
        }
      }

      if(print_read_problem) {
        algo.print_read_problem(context);
      }

#ifdef F4NCGB_ENABLE_STORE_DUMP
      algo.set_dump_paths(context.get_filename()).filename().string());
#endif

      if(stats_print_function) {
        *stats_print_function = [&algo]() {
          msg("Current basis:");
          algo.write_basis(std::cerr);
        };
      }

      if(verbose > 0)
        algo.print_parameters(context, output_name);

      if(context.reduce())
        algo.reduced_form();
      else
        algo.compute_basis();

      if(verbose > 0)
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
