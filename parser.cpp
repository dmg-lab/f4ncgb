#include <cerrno>
#include <cstdio>
#include <cstring>

#include "parser.hpp"

#define xstr(s) str(s)
#define str(s) #s

#define EXPECT(C)                                                            \
  do {                                                                       \
    if(c != C) {                                                             \
      return ctx.error_at_current(std::format(                               \
        "{} {}", "Expected " str(C) " but received", static_cast<char>(c))); \
    }                                                                        \
  } while(false)

static inline bool
isdigit_(char c) {
  return c >= '0' && c <= '9';
}

#define EXPECT_DIGIT()                                             \
  do {                                                             \
    if(!isdigit_(c)) {                                             \
      return ctx.error_at_current(                                 \
        std::format("{} {}", "Expected a digit but received", c)); \
    }                                                              \
  } while(false)

#define READ(C)   \
  c = ctx.getc(); \
  EXPECT(C);

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
                  void* boundary_cb_userdata,
                  char first_char,
                  char second_char) {
  // TODO
  return std::nullopt;
}

static parse_res
parse_impl_poly(parser_context& ctx,
                parse_add_cb add_cb,
                void* add_cb_userdata,
                parse_monomial_boundary_cb boundary_cb,
                void* boundary_cb_userdata,
                size_t var_count,
                size_t var_blocks,
                size_t characteristic) {
  // TODO
  return std::nullopt;
}

static parse_res
parse_impl_sympoly(parser_context& ctx,
                   parse_add_cb add_cb,
                   void* add_cb_userdata,
                   parse_monomial_boundary_cb boundary_cb,
                   void* boundary_cb_userdata,
                   size_t var_count,
                   size_t var_blocks,
                   size_t characteristic) {
  // TODO
  return std::nullopt;
}

parse_res
parse_impl(parser_context& ctx,
           parse_add_cb add_cb,
           void* add_cb_userdata,
           parse_monomial_boundary_cb boundary_cb,
           void* boundary_cb_userdata) {
  // Parsing procedure:
  //
  // 1. Decide on format to use based on first few read characters.
  // 2. Call the respective implementation.
  // 3. Let the implementation call the callbacks.
  //
  // Return eventual errors as strings.

  int c = ctx.getc();
  if(c != 'p') {
    int second = ctx.getc();
    return parse_impl_msolve(ctx,
                             add_cb,
                             add_cb_userdata,
                             boundary_cb,
                             boundary_cb_userdata,
                             c,
                             second);
  }

  // File reads "p"
  c = ctx.getc();
  if(c != ' ') {
    return parse_impl_msolve(
      ctx, add_cb, add_cb_userdata, boundary_cb, boundary_cb_userdata, 'p', c);
  }

  // File reads "p "

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
  } else if(c == 'p') {
    // File should read "p poly"
    READ('o');
    READ('l');
    READ('y');
    READ(' ');
  }

  return std::nullopt;
}
}
