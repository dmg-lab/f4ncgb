#pragma once

#include <filesystem>
#include <format>
#include <memory>
#include <optional>
#include <string>
#include <type_traits>

#include <boost/bimap.hpp>
#include <boost/bimap/unordered_set_of.hpp>

namespace kommunopp {
using parse_res = std::optional<std::string>;
using parse_add_cb = std::add_pointer<parse_res(void*, uint32_t)>::type;

// Also gives the coefficent of the monomial as either a rational or an integer.
// The bool is true if it is a rational number. The first long is the numerator,
// the second long is the denominator. If integer, the denominator is
// always 1.
using parse_monomial_boundary_cb
  = std::add_pointer<parse_res(void*, long, long, bool)>::type;

class parser_symbolic_context {
  public:
  using id = uint32_t;

  private:
  using map = boost::bimap<boost::bimaps::unordered_set_of<std::string>,
                           boost::bimaps::unordered_set_of<id>>;
  map m;
  id counter = 1;

  public:
  parser_symbolic_context() = default;
  ~parser_symbolic_context() = default;

  inline std::string_view id_to_str(id i) { return m.right.at(i); }
  inline id str_to_id(const std::string& s) {
    auto it = m.left.find(s);
    if(it != m.left.end()) {
      return it->second;
    }
    m.insert(map::value_type(s, counter));
    return counter++;
  }
};

class parser_context {
  std::unique_ptr<parser_symbolic_context> symbols;
  struct FILE_deleter {
    void operator()(FILE* f) { fclose(f); }
  };
  std::unique_ptr<FILE, FILE_deleter> f;

  size_t last_line = 0;
  size_t last_col = 0;
  size_t line = 1;
  size_t col = 1;
  std::string filename = "";
  std::string ident = "";
  int c = 0;

  public:
  parser_context(FILE* in)
    : f(in) {}
  parser_context() = default;
  ~parser_context() = default;

  size_t characteristic = 0;
  size_t num_blocks = 0;
  size_t num_vars = 0;
  std::function<parse_res(parser_context& ctx,
                          parse_add_cb add_cb,
                          void* add_cb_userdata,
                          parse_monomial_boundary_cb boundary_cb,
                          void* boundary_cb_userdata)>
    impl;

  parse_res open(std::filesystem::path p);

  inline int current_c() const { return c; }

  inline int getc() {
    assert(f);
    last_line = line;
    last_col = col;

    int c = fgetc(f.get());

    // Debug aid:
    // std::printf("Read '%c'\n", c);

    this->c = c;
    if(c == '\n') {
      ++line;
      col = 1;
    } else {
      ++col;
    }
    return c;
  }

  inline void ungetc(int c) { ::ungetc(c, f.get()); }

  inline size_t read_positive_int() {
    assert(c >= '0' && c <= '9');
    size_t num = 0;
    do {
      num *= 10;
      num += (static_cast<size_t>(c) - '0');
      getc();
    } while(c >= '0' && c <= '9');
    return num;
  }

  inline long read_positive_or_negative_int() {
    assert((c >= '0' && c <= '9') || c == '-' || c == '+');
    long num = 0;
    bool neg = false;
    if(c == '-') {
      getc();
      neg = true;
    } else if(c == '+') {
      getc();
      neg = false;
    }
    do {
      num *= 10;
      num += (static_cast<size_t>(c) - '0');
      getc();
    } while(c >= '0' && c <= '9');
    return num * (neg ? -1 : 1);
  }

  inline int swallow_whitespace_no_newline() {
    while(c == ' ' || c == '\t') {
      getc();
    }
    return c;
  }

  static inline bool msolve_ident_filter(char c) {
    return c != '-' && c != '/' && c != '*' && c != '+' && c != ',' && c != ' '
           && c != '\n' && c != '\t';
  }
  static inline bool sympoly_ident_filter(char c) {
    return c != ' ' && c != '\n' && c != '\t';
  }

  template<typename Filter>
  inline const std::string& read_ident(Filter c_allowed) {
    ident = "";
    do {
      ident.push_back(c);
      getc();
    } while(c_allowed(c));
    return ident;
  }

  std::string error_at_current(std::string_view msg) {
    if(filename != "") {
      return std::format(
        "{}:{}:{}: error: {}", filename, last_line, last_col, msg);
    } else {
      return std::format("{}:{}: error: {}", last_line, last_col, msg);
    }
  }

  inline std::string_view id_to_str(parser_symbolic_context::id i) {
    return symbols->id_to_str(i);
  }
  inline parser_symbolic_context::id str_to_id(const std::string& s) {
    return symbols->str_to_id(s);
  }

  void init_symbols() { symbols = std::make_unique<parser_symbolic_context>(); }
};

parse_res
parse_header(parser_context& ctx);

template<typename AddCB, typename BoundaryCB>
parse_res
parse_rest(parser_context& ctx, AddCB add_cb, BoundaryCB boundary_cb) {
  parse_add_cb add_cb_impl = [](void* userdata, uint32_t i) -> parse_res {
    AddCB* cb = reinterpret_cast<AddCB*>(userdata);
    return (*cb)(i);
  };
  parse_monomial_boundary_cb boundary_cb_impl
    = [](void* userdata,
         long numerator,
         long denominator,
         bool is_rational) -> parse_res {
    BoundaryCB* cb = reinterpret_cast<BoundaryCB*>(userdata);
    return (*cb)(numerator, denominator, is_rational);
  };

  return ctx.impl(ctx, add_cb_impl, &add_cb, boundary_cb_impl, &boundary_cb);
}

template<typename AddCB, typename BoundaryCB>
parse_res
parse(parser_context& ctx, AddCB add_cb, BoundaryCB boundary_cb) {
  parse_res err = parse_header(ctx);
  if(err)
    return err;
  return parse_rest(ctx, add_cb, boundary_cb);
}

template<class PS>
parse_res
parse_rest_into_polynomial_store(parser_context& ctx, PS& s) {
  using monomial_store = PS::monomial_store_;
  using monomial_ref = PS::value_type;
  using coefficient = PS::coefficient;
  using var = PS::monomial;
  std::vector<var> tmp_monomial;
  std::vector<std::pair<coefficient, monomial_ref>> tmp_polynomial;

  auto add_cb = [&s, &tmp_monomial, &tmp_polynomial](uint32_t id) {
    if(id) {
      tmp_monomial.push_back(id);
    } else {
      s.add_polynomial(tmp_polynomial);
      tmp_polynomial.clear();
    }
    return std::nullopt;
  };

  auto boundary_cb = [&s, &tmp_monomial, &tmp_polynomial](
                       long numerator, long denominator, bool is_rational) {
    (void)is_rational;

    monomial_store& m = s.get_monomial_store();
    auto monomial_id = m.getid(tmp_monomial);
    tmp_polynomial.push_back(
      std::make_pair(coefficient(numerator, denominator), monomial_id));

    tmp_monomial.clear();
    return std::nullopt;
  };

  return parse_rest(ctx, add_cb, boundary_cb);
}

}
