#include <array>
#include <cstdio>
#include <filesystem>
#include <string>

#include <boost/test/unit_test.hpp>

#include "parser.hpp"

using namespace kommunopp;

std::optional<std::filesystem::path>
get_parse_tests_location() {
  std::array<std::string, 2> candidates = { "parse_tests/", "../parse_tests/" };

  auto curr = std::filesystem::current_path();

  for(auto& c : candidates) {
    auto path = curr / std::filesystem::path(c);
    if(std::filesystem::exists(path)) {
      return path;
    }
  }
  return std::nullopt;
}

// These tests only test the parser code itself, not the mechanics underneath.
// They are fed using some input files in the parse_tests subdirectory.

static std::optional<std::string>
null_parser(std::filesystem::path p) {
  parser_context ctx;
  if(auto err = ctx.open(p); err) {
    return err;
  }
  return parse(
    ctx, [](int) { return std::nullopt; }, []() { return std::nullopt; });
}

BOOST_AUTO_TEST_CASE(parse_files) {
  auto dir_optional = get_parse_tests_location();
  BOOST_REQUIRE(dir_optional.has_value());
  std::filesystem::path dir = *dir_optional;

  for(auto& e : std::filesystem::directory_iterator(dir)) {
    BOOST_REQUIRE(e.is_regular_file());
    if(e.path().filename().string().starts_with("fail_")) {
      // Expect some parse error.
      auto ret = null_parser(e);
      BOOST_CHECK_MESSAGE(ret.has_value(),
                          "parsing the input "
                            << e.path().filename()
                            << " did not lead to an error as expected!");
    } else {
      // Expect the parsing to be successful.
      auto ret = null_parser(e);
      BOOST_CHECK_MESSAGE(!ret.has_value(),
                          "parsing the input " << e.path().filename()
                                               << " did lead to an error:\n"
                                               << *ret);
    }
  }
}
