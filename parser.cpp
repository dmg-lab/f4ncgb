#include <cerrno>
#include <cstdio>
#include <cstring>

#include "parser.hpp"

#define xstr(s) str(s)
#define str(s) #s

#define RETURN_ERROR(MSG)                                       \
  do {                                                          \
    std::source_location loc = std::source_location::current(); \
    return ctx.error_at_current(                                \
      std::format("{}\n{}:{}:{}: location of last error",       \
                  MSG,                                          \
                  loc.file_name(),                              \
                  loc.line(),                                   \
                  loc.column()));                               \
  } while(false)

#define EXPECT(C)                                                 \
  do {                                                            \
    if(ctx.current_c() != C) {                                    \
      std::source_location loc = std::source_location::current(); \
      return ctx.error_at_current(                                \
        std::format("{} '{}'\n{}:{}:{}: location of last error",  \
                    "Expected " str(C) " but received",           \
                    static_cast<char>(ctx.current_c()),           \
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
    if(!isdigit_(ctx.current_c())) {                              \
      std::source_location loc = std::source_location::current(); \
      return ctx.error_at_current(                                \
        std::format("{} {}\n{}:{}:{}: location of last error",    \
                    "Expected a digit but received",              \
                    static_cast<char>(ctx.current_c()),           \
                    loc.file_name(),                              \
                    loc.line(),                                   \
                    loc.column()));                               \
    }                                                             \
  } while(false)

#define EXPECT_POSNEGDIGIT()                                      \
  do {                                                            \
    if(!isdigit_(ctx.current_c()) && c != '-' && c != '+') {      \
      std::source_location loc = std::source_location::current(); \
      return ctx.error_at_current(                                \
        std::format("{} {}\n{}:{}:{}: location of last error",    \
                    "Expected a digit but received",              \
                    static_cast<char>(ctx.current_c()),           \
                    loc.file_name(),                              \
                    loc.line(),                                   \
                    loc.column()));                               \
    }                                                             \
  } while(false)

#define READ(C)   \
  c = ctx.getc(); \
  EXPECT(C);

#define ADD(V) add_cb(add_cb_userdata, V)

#define BOUNDARY(N, D, R) boundary_cb(boundary_cb_userdata, N, D, R)

namespace kommunopp {
parse_res
parser_context::open(std::filesystem::path p) {
  if(!std::filesystem::exists(p)) {
    return std::format("File \"{}\" does not exist!", p.string());
  }
  FILE* f_ptr = std::fopen(p.c_str(), "r");
  if(!f_ptr) {
    return std::format(
      "Could not open file \"{}\", error: {}", p.string(), strerror(errno));
  }
  f.reset(f_ptr);
  filename = p.string();
  return std::nullopt;
}

static parse_res
parse_impl_msolve(parser_context& ctx,
                  parse_add_cb add_cb,
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
      c = ctx.swallow_whitespace_no_newline();

      if(c == '+' || c == '-' || isdigit_(c)) {
        EXPECT_POSNEGDIGIT();
        numerator = ctx.read_positive_or_negative_int();
        c = ctx.swallow_whitespace_no_newline();
        denominator = 1;
        if(c == '/') {
          c = ctx.getc();
          denominator = ctx.read_positive_int();
          rational = true;
        }
        c = ctx.swallow_whitespace_no_newline();
      } else {
        // Trick the rest of the parser into thinking it read a coefficient and
        // a *, but it actually read nothing.
        ctx.ungetc(c);
        c = '*';
      }

      // Parse the variables
      while(c == '*') {
        c = ctx.getc();
        c = ctx.swallow_whitespace_no_newline();

        const std::string& ident
          = ctx.read_ident(parser_context::msolve_ident_filter);
        auto id = ctx.str_to_id(ident);

        ADD(id);

        c = ctx.swallow_whitespace_no_newline();
      }

      // Submit the current monomial.
      BOUNDARY(numerator, static_cast<long>(denominator), rational);
    } while(c == '+' || c == '-');

    // Maybe another polynomial, or break.
    if(c == ',') {
      READ('\n');
      c = ctx.getc();

      // Boundary between polynomials.
      ADD(0);
    } else if(c == '\n') {
      // Very last polynomial.
      ADD(0);
      break;
    } else {
      RETURN_ERROR(
        std::format("Unexpected line termination: '{}'", static_cast<char>(c)));
    }
  }

  return std::nullopt;
}

static parse_res
parse_impl_msolve_header(parser_context& ctx,
                         char first_char,
                         char second_char) {
  // Handle the first two characters first, then process the rest of the input.

  int c = ctx.current_c();

  // Msolve always has 1 block.
  ctx.num_blocks = 1;

  if(first_char == ',')
    return "First character of msolve must be some ident, not ','";
  if(parser_context::msolve_ident_filter(first_char) && second_char == ',') {
    // First char was some variable that has to be counted.
    ++ctx.num_vars;
  } else if(parser_context::msolve_ident_filter(first_char)
            && parser_context::msolve_ident_filter(second_char)) {
    c = ctx.getc();
    if(c == ',') {
      // Some variable name was started and finished with the third char.
      ++ctx.num_vars;
    } else if(parser_context::msolve_ident_filter(c)) {
      // Some variable name was started but not finished yet, read until finish.
      ctx.read_ident(parser_context::msolve_ident_filter);
    } else {
      RETURN_ERROR(
        std::format("Unexpected third character for msolve format: '{}'", c));
    }
  } else {
    RETURN_ERROR(std::format(
      "Unexpected start for msolve format: \"{}{}\"", first_char, second_char));
  }

  // Now, we have the first ident and are either at a ',' or at the end of the
  // variables with '\n'
  assert(ctx.num_vars == 1);
  if(c != '\n' && c != ',' && c != ' ') {
    RETURN_ERROR(std::format(
      "Unexpected msolve parse state after first variable, have c: '{}'", c));
  }

  while(true) {
    c = ctx.swallow_whitespace_no_newline();

    // No more variables.
    if(c == '\n') {
      break;
    } else if(c == ',') {
      ctx.swallow_whitespace_no_newline();
      ctx.read_ident(parser_context::msolve_ident_filter);
      ++ctx.num_vars;
    } else {
      RETURN_ERROR(std::format("Unexpected character '{}'", c));
    }
  }

  // Step over newline.
  EXPECT('\n');
  c = ctx.getc();

  // Read the characteristic.
  EXPECT_DIGIT();
  ctx.characteristic = ctx.read_positive_int();

  // Step over newline.
  EXPECT('\n');
  c = ctx.getc();

  ctx.impl = [](parser_context& ctx,
                parse_add_cb add_cb,
                void* add_cb_userdata,
                parse_monomial_boundary_cb boundary_cb,
                void* boundary_cb_userdata) {
    return parse_impl_msolve(
      ctx, add_cb, add_cb_userdata, boundary_cb, boundary_cb_userdata);
  };

  return std::nullopt;
}

static parse_res
parse_impl_poly(parser_context& ctx,
                parse_add_cb add_cb,
                void* add_cb_userdata,
                parse_monomial_boundary_cb boundary_cb,
                void* boundary_cb_userdata) {
  // TODO
  return std::nullopt;
}

static parse_res
parse_impl_sympoly(parser_context& ctx,
                   parse_add_cb add_cb,
                   void* add_cb_userdata,
                   parse_monomial_boundary_cb boundary_cb,
                   void* boundary_cb_userdata) {
  // TODO
  return std::nullopt;
}

parse_res
parse_header(parser_context& ctx) {
  // Parsing procedure:
  //
  // 1. Decide on format to use based on first few read characters.
  // 2. Call the respective implementation.
  // 3. Let the implementation call the callbacks.
  //
  // Return eventual errors as strings.

  int c = ctx.getc();
  if(c != 'p') {
    ctx.init_symbols();
    return parse_impl_msolve_header(ctx, 'p', ctx.getc());
  }

  // File reads "p"
  c = ctx.getc();
  if(c != ' ') {
    ctx.init_symbols();
    return parse_impl_msolve_header(ctx, 'p', c);
  }

  // File reads "p "

  using impl_fun_type = std::add_pointer<decltype(parse_impl_sympoly)>::type;
  impl_fun_type impl_fun = nullptr;

  // Now the file must be one of our own new format types!
  c = ctx.getc();
  if(c == 's') {
    // File should read "p sympoly"
    READ('y');
    READ('m');
    READ('p');
    READ('o');
    READ('l');
    READ('y');
    READ(' ');

    impl_fun = parse_impl_sympoly;
    ctx.init_symbols();
  } else if(c == 'p') {
    // File should read "p poly"
    READ('o');
    READ('l');
    READ('y');
    READ(' ');

    impl_fun = parse_impl_poly;
  } else {
    return std::format("character {} is not a valid start for problem "
                       "specifier, must be either s or p",
                       c);
  }

  c = ctx.getc();
  EXPECT_DIGIT();

  ctx.num_vars = ctx.read_positive_int();
  c = ctx.current_c();
  EXPECT(' ');
  c = ctx.getc();
  EXPECT_DIGIT();
  ctx.num_blocks = ctx.read_positive_int();
  c = ctx.current_c();
  EXPECT(' ');
  c = ctx.getc();
  EXPECT_DIGIT();
  ctx.characteristic = ctx.read_positive_int();
  c = ctx.current_c();
  EXPECT('\n');

  ctx.impl = impl_fun;
  return std::nullopt;
}
}
