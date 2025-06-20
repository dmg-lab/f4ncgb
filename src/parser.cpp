#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "parser.hpp"
#include "store.hpp"
#include "zip_utils.hpp"

#define xstr(s) str(s)
#define str(s) #s

#define RETURN_ERROR(MSG)                                       \
  do {                                                          \
    std::source_location loc = std::source_location::current(); \
    return error_at_current(                                    \
      std::format("{}\n{}:{}:{}: location of last error",       \
                  MSG,                                          \
                  loc.file_name(),                              \
                  loc.line(),                                   \
                  loc.column()));                               \
  } while(false)

#define EXPECT(C)                                                 \
  do {                                                            \
    if(current_c() != C) {                                        \
      std::source_location loc = std::source_location::current(); \
      return error_at_current(                                    \
        std::format("{} '{}'\n{}:{}:{}: location of last error",  \
                    "Expected " str(C) " but received",           \
                    static_cast<char>(current_c()),               \
                    loc.file_name(),                              \
                    loc.line(),                                   \
                    loc.column()));                               \
    }                                                             \
  } while(false)

static inline bool
isdigit_(char c) {
  return c >= '0' && c <= '9';
}

#define EXPECT_DIGIT()                                            \
  do {                                                            \
    if(!isdigit_(current_c())) {                                  \
      std::source_location loc = std::source_location::current(); \
      return error_at_current(                                    \
        std::format("{} {}\n{}:{}:{}: location of last error",    \
                    "Expected a digit but received",              \
                    static_cast<char>(current_c()),               \
                    loc.file_name(),                              \
                    loc.line(),                                   \
                    loc.column()));                               \
    }                                                             \
  } while(false)

#define EXPECT_POSNEGDIGIT()                                      \
  do {                                                            \
    if(!isdigit_(current_c()) && c != '-' && c != '+') {          \
      std::source_location loc = std::source_location::current(); \
      return error_at_current(                                    \
        std::format("{} {}\n{}:{}:{}: location of last error",    \
                    "Expected a digit but received",              \
                    static_cast<char>(current_c()),               \
                    loc.file_name(),                              \
                    loc.line(),                                   \
                    loc.column()));                               \
    }                                                             \
  } while(false)

#define READ(C) \
  c = getc();   \
  EXPECT(C);

#define ADD(V) add_cb(add_cb_userdata, V)

#define BOUNDARY(N, D, R) boundary_cb(boundary_cb_userdata, N, D, R)

namespace f4ncgb {
parse_res
parser_context::open(std::filesystem::path p) {
  if(!std::filesystem::exists(p)) {
    return std::format("File \"{}\" does not exist!", p.string());
  }

  FILE* f_ptr = f4ncgb::open_compressed_or_direct(p);

  if(!f_ptr) {
    return std::format(
      "Could not open file \"{}\", error: {}", p.string(), strerror(errno));
  }
  f.reset(f_ptr);
  filename = p.string();
  return std::nullopt;
}

parse_res
parser_context::impl_msolve(parse_add_cb add_cb,
                            void* add_cb_userdata,
                            parse_monomial_boundary_cb boundary_cb,
                            void* boundary_cb_userdata) {
  int c = 0;
  // Read polynomials.

  // Each polynomial is a line ending with ',' except the last one. They are
  // made up of Monomials. Monomials always start with some coefficient
  // (either rational or integer), optionally followed by variables.

  while(true) {
    // Parse a polynomial.

    // Parse a monomial.
    do {

      // Parse the coefficient.
      bool rational = false;
      long numerator = 1;
      size_t denominator = 1;
      c = swallow_whitespace_no_newline();

      if(c == '-') {
        numerator = -1;
        c = getc();
      } else if(c == '+') {
        numerator = 1;
        c = getc();
      }

      c = swallow_whitespace_no_newline();

      if(isdigit_(c)) {
        EXPECT_POSNEGDIGIT();
        numerator = read_positive_or_negative_int() * numerator;
        c = swallow_whitespace_no_newline();
        denominator = 1;
        if(c == '/') {
          c = getc();
          denominator = read_positive_int();
          rational = true;
        }
        c = swallow_whitespace_no_newline();
      } else {
        // Trick the rest of the parser into thinking it read a coefficient and
        // a *, but it actually read nothing.
        ungetc(c);
        c = '*';
      }

      // Parse the variables
      while(c == '*') {
        c = getc();
        c = swallow_whitespace_no_newline();

        const std::string& ident
          = read_ident(parser_context::msolve_ident_filter);
        auto id = str_to_id(ident);

        c = swallow_whitespace_no_newline();

        if(c == '^') {
          getc();
          c = swallow_whitespace_no_newline();
          EXPECT_DIGIT();
          size_t sup = read_positive_int();
          if(sup > 1000) {
            return "Exponent too high";
          }
          c = swallow_whitespace_no_newline();
          for(size_t i = 0; i < sup; ++i) {
            ADD(id);
          }
        } else {
          ADD(id);
        }
      }

      // Submit the current monomial.
      BOUNDARY(numerator, static_cast<long>(denominator), rational);
    } while(c == '+' || c == '-');

    // Maybe another polynomial, or break.
    if(c == ',') {
      c = getc();
      if(c == '\r') {
        c = getc();
      }
      if(c == '\n') {
        c = getc();
      }

      // Boundary between polynomials.
      ADD(0);
    } else if(c == '\n' || c == '\r') {
      // Very last polynomial.
      if(c == '\r') {
        c = getc();
      }
      ADD(0);
      break;
    } else {
      // Also very last polynomial with some weird ending.
      ADD(0);
      break;
    }
  }

  return std::nullopt;
}

parse_res
parser_context::impl_msolve_header(char first_char, char second_char) {
  // Handle the first two characters first, then process the rest of the input.

  // Msolve always has 1 block.
  num_blocks_ = 1;

  if(first_char == ',') {
    RETURN_ERROR("First character of msolve must be some ident, not ','");
  }
  if(parser_context::msolve_ident_filter(first_char)
     && (second_char == ',' || second_char == ' ' || second_char == '\n')) {
    // First char was some variable that has to be counted.
    std::string ident;
    ident += first_char;
    str_to_id(ident);
    ++num_vars_;
  } else if(parser_context::msolve_ident_filter(first_char)
            && parser_context::msolve_ident_filter(second_char)) {
    c = getc();
    if(c == ',') {
      std::string ident;
      ident += first_char;
      ident += second_char;
      str_to_id(ident);
      // Some variable name was started and finished with the third char.
      ++num_vars_;
    } else if(parser_context::msolve_ident_filter(c)) {
      // Some variable name was started but not finished yet, read until finish.
      std::string ident(read_ident(parser_context::msolve_ident_filter));
      std::string prefix;
      prefix += first_char;
      prefix += second_char;
      ident.insert(0, prefix);
      str_to_id(ident);
      ++num_vars_;
    } else {
      RETURN_ERROR(
        std::format("Unexpected third character for msolve format: '{}'", c));
    }
  } else if(first_char == ' ' && second_char == ' ') {
    // Swallow up to the first new character.
    c = swallow_whitespace_and_newline();
  } else {
    RETURN_ERROR(std::format(
      "Unexpected start for msolve format: \"{}{}\"", first_char, second_char));
  }

  // Now, we have the first ident and are either at a ',' or at the end of the
  // variables with '\n'
  assert(num_vars_ == 1);
  if(c != '\n' && c != ',' && c != ' ') {
    RETURN_ERROR(std::format(
      "Unexpected msolve parse state after first variable, have c: '{}'", c));
  }

  blocks_.emplace_back();
  blocks_[0].emplace_back(1);

  auto* current_block = &blocks_[0];

  while(true) {
    c = swallow_whitespace_no_newline();

    if(c == '\r') {
      c = getc();
    }

    // No more variables.
    if(c == '\n') {
      break;
    } else if(c == ',') {
      getc();
      swallow_whitespace_no_newline();
      if(c == '\n') {
        getc();
        blocks_.emplace_back();
        current_block = &blocks_[blocks_.size() - 1];
        // Swallow as much as needed until there is something again.
        swallow_whitespace_and_newline();
      }
      auto id = str_to_id(read_ident(parser_context::msolve_ident_filter));
      current_block->emplace_back(id);
      ++num_vars_;
    } else {
      RETURN_ERROR(std::format("Unexpected character '{}'", c));
    }
  }

  // Step over newline.
  if(c == '\r') {
    c = getc();
  }
  EXPECT('\n');
  c = getc();

  // Read the characteristic.
  EXPECT_DIGIT();
  characteristic_ = read_positive_int();

  // Step over newline.
  if(c == '\r') {
    c = getc();
  }
  EXPECT('\n');
  c = getc();

  impl = [this](parse_add_cb add_cb,
                void* add_cb_userdata,
                parse_monomial_boundary_cb boundary_cb,
                void* boundary_cb_userdata) {
    return impl_msolve(
      add_cb, add_cb_userdata, boundary_cb, boundary_cb_userdata);
  };

  return std::nullopt;
}

template<parser_symbolic_context::id (*getV)(parser_context&)>
parse_res
parser_context::impl_poly_gen(parse_add_cb add_cb,
                              void* add_cb_userdata,
                              parse_monomial_boundary_cb boundary_cb,
                              void* boundary_cb_userdata) {
  // Variable Blocks, then polynomials.

  // The number of variables was already given in the header. Now, the blocks
  // are read.
  size_t block = 0;
  while(block < num_blocks_) {
    if constexpr(std::is_same_v<typeof(getV), typeof(read_numeric_var)>) {
      EXPECT_DIGIT();
    }
    parser_symbolic_context::id v = getV(*this);
    if(v == 0) {
      ++block;
    } else {
      blocks_[block].push_back(v);
    }
    swallow_whitespace_and_newline();
  }

  // Now, read polynomials. Each polynomial is separated by a 0, with each
  // monomial within a polynomial being separated by a +. Negative coefficients
  // are always separated using a + and then just start with a -.

  while(c != EOF) {
    parser_symbolic_context::id v;

    // Polynomial
    do {
      // Monomial
      do {
        // Monomial Coefficient
        swallow_whitespace_and_newline();
        EXPECT_POSNEGDIGIT();
        long numerator = read_positive_or_negative_int();
        swallow_whitespace_and_newline();
        EXPECT_DIGIT();
        long denominator = static_cast<long>(read_positive_int());

        if(denominator == 0) {
          RETURN_ERROR("Denominator must not be 0!");
        }

        swallow_whitespace_and_newline();

        // Monomial Vars
        while(c != '+' && c != '0' && c != EOF) {
          v = getV(*this);
          ADD(v);
          swallow_whitespace_and_newline();
        }

        BOUNDARY(numerator, denominator, numerator == denominator);

        // Prepare next monomial. If '0', nothing has to be done.
        if(c == '+') {
          getc();
          swallow_whitespace_and_newline();
        } else {
          break;
        }
      } while(c != EOF);

      // Prepare next polynomial.
      if(c == '0') {
        getc();
        swallow_whitespace_and_newline();
        ADD(0);
      } else {
        break;
      }
    } while(c != EOF);
  }

  return std::nullopt;
}

template parse_res
parser_context::impl_poly_gen<&parser_context::read_numeric_var>(
  parse_add_cb add_cb,
  void* add_cb_userdata,
  parse_monomial_boundary_cb boundary_cb,
  void* boundary_cb_userdata);

template parse_res
parser_context::impl_poly_gen<&parser_context::read_symbolic_var>(
  parse_add_cb add_cb,
  void* add_cb_userdata,
  parse_monomial_boundary_cb boundary_cb,
  void* boundary_cb_userdata);

parse_res
parser_context::impl_poly(parse_add_cb add_cb,
                          void* add_cb_userdata,
                          parse_monomial_boundary_cb boundary_cb,
                          void* boundary_cb_userdata) {
  return impl_poly_gen<&parser_context::read_numeric_var>(
    add_cb, add_cb_userdata, boundary_cb, boundary_cb_userdata);
}

parse_res
parser_context::impl_sympoly(parse_add_cb add_cb,
                             void* add_cb_userdata,
                             parse_monomial_boundary_cb boundary_cb,
                             void* boundary_cb_userdata) {
  return impl_poly_gen<&parser_context::read_symbolic_var>(
    add_cb, add_cb_userdata, boundary_cb, boundary_cb_userdata);
}

parse_res
parser_context::parse_header() {
  // Parsing procedure:
  //
  // 1. Decide on format to use based on first few read characters.
  // 2. Call the respective implementation.
  // 3. Let the implementation call the callbacks.
  //
  // Return eventual errors as strings.

  getc();
  if(c != 'p') {
    init_symbols();
    char first = c;
    return impl_msolve_header(first, getc());
  }

  // File reads "p"
  c = getc();
  if(c != ' ') {
    init_symbols();
    return impl_msolve_header('p', c);
  }

  // File reads "p "

  using impl_fun_type
    = std::function<parse_res(parse_add_cb add_cb,
                              void* add_cb_userdata,
                              parse_monomial_boundary_cb boundary_cb,
                              void* boundary_cb_userdata)>;
  impl_fun_type impl_fun = nullptr;

  // Now the file must be one of our own new format types!
  getc();
  if(c == 's') {
    // File should read "p sympoly"
    READ('y');
    READ('m');
    READ('p');
    READ('o');
    READ('l');
    READ('y');
    READ(' ');

    impl_fun = [this](parse_add_cb add_cb,
                      void* add_cb_userdata,
                      parse_monomial_boundary_cb boundary_cb,
                      void* boundary_cb_userdata) {
      return impl_sympoly(
        add_cb, add_cb_userdata, boundary_cb, boundary_cb_userdata);
    };
    init_symbols();
  } else if(c == 'p') {
    // File should read "p poly"
    READ('o');
    READ('l');
    READ('y');
    READ(' ');

    impl_fun = [this](parse_add_cb add_cb,
                      void* add_cb_userdata,
                      parse_monomial_boundary_cb boundary_cb,
                      void* boundary_cb_userdata) {
      return impl_poly(
        add_cb, add_cb_userdata, boundary_cb, boundary_cb_userdata);
    };
  } else {
    return std::format("character {} is not a valid start for problem "
                       "specifier, must be either s or p",
                       c);
  }

  getc();
  EXPECT_DIGIT();

  num_vars_ = read_positive_int();
  EXPECT(' ');
  c = getc();
  EXPECT_DIGIT();
  num_blocks_ = read_positive_int();
  if(num_blocks_ > 100)
    return "too many blocks";
  blocks_.resize(num_blocks_);
  EXPECT(' ');
  c = getc();
  EXPECT_DIGIT();
  characteristic_ = read_positive_int();
  EXPECT('\n');
  getc();// Swallow the newline.

  impl = impl_fun;
  return std::nullopt;
}

#ifdef __linux__
void
dummy_parse_string(std::string input) {
  (void)input;
  int add_cb_count = 0;
  int boundary_cb_count = 0;

  auto add_cb = [&add_cb_count](uint32_t i) {
    (void)i;
    if(add_cb_count++ > 10000000) {
      assert(false);
    }
    return std::nullopt;
  };
  auto boundary_cb
    = [&boundary_cb_count](long numerator, long denominator, bool is_rational) {
        (void)numerator;
        (void)denominator;
        (void)is_rational;
        if(boundary_cb_count++ > 10000000) {
          assert(false);
        }
        return std::nullopt;
      };

  FILE* f = fmemopen((void*)input.c_str(), input.size(), "r");

  if(!f) {
    std::cout << "Could not make a memory-backed file. Error: "
              << strerror(errno) << std::endl;
    exit(1);
  }

  parser_context ctx(f);

  parse(ctx, add_cb, boundary_cb);
}
#else
void
dummy_parse_string(std::string input) {}
#endif

}
