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
using parse_add_cb = std::add_pointer<parse_res(void*, int)>::type;
using parse_monomial_boundary_cb = std::add_pointer<parse_res(void*)>::type;

class parser_symbolic_context {
  public:
  using id = uint32_t;

  private:
  using map = boost::bimap<boost::bimaps::unordered_set_of<std::string>,
                           boost::bimaps::unordered_set_of<id>>;
  map m;
  id counter = 0;

  public:
  parser_symbolic_context() = default;
  ~parser_symbolic_context() = default;

  inline std::string_view id_to_str(id i) { return m.right.at(i); }
  inline id str_to_id(const std::string& s) { return m.left.at(s); }
  inline id insert_str_to_id(const std::string& s) {
    m.insert(map::value_type(s, counter++));
    return m.left.at(s);
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
  size_t line = 0;
  size_t col = 0;
  std::string filename = "";

  public:
  parser_context(FILE* in)
    : f(in) {}
  parser_context() = default;
  ~parser_context() = default;

  parse_res open(std::filesystem::path p);

  inline int getc() {
    last_line = line;
    last_col = col;

    int c = fgetc(f.get());
    if(c == '\n') {
      ++line;
      col = 0;
    } else {
      ++col;
    }
    return c;
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
};

parse_res
parse_impl(parser_context& ctx,
           parse_add_cb add_cb,
           void* add_cb_userdata,
           parse_monomial_boundary_cb boundary_cb,
           void* boundary_cb_userdata);

template<typename AddCB, typename BoundaryCB>
parse_res
parse(parser_context& ctx, AddCB add_cb, BoundaryCB boundary_cb) {
  parse_add_cb add_cb_impl = [](void* userdata, int i) -> parse_res {
    AddCB* cb = reinterpret_cast<AddCB*>(userdata);
    return (*cb)(i);
  };
  parse_monomial_boundary_cb boundary_cb_impl
    = [](void* userdata) -> parse_res {
    BoundaryCB* cb = reinterpret_cast<BoundaryCB*>(userdata);
    return (*cb)();
  };

  return parse_impl(ctx, add_cb_impl, &add_cb, boundary_cb_impl, &boundary_cb);
}
}
