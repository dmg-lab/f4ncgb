#include <boost/test/unit_test.hpp>

#include "./aho_corasick/src/aho_corasick/aho_corasick.hpp"
#include "kommunopp.hpp"
#include "monomial_trie.hpp"
#include "parser.hpp"

using namespace kommunopp;

std::optional<std::filesystem::path>
get_test_input_files_location() {
  std::array<std::string, 2> candidates = { "test_inputs/", "../test_inputs/" };

  auto curr = std::filesystem::current_path();

  for(auto& c : candidates) {
    auto path = curr / std::filesystem::path(c);
    if(std::filesystem::exists(path)) {
      return path;
    }
  }
  return std::nullopt;
}

template<internal::metadata_concept MM = internal::metadata<uint8_t>,
         internal::metadata_concept PM = internal::metadata<uint8_t>,
         internal::value_concept V = uint8_t,
         typename I = uint32_t,
         typename C = boost::multiprecision::gmp_rational>
struct impl {
  using var = V;
  using idx = I;
  using coefficient = C;
  using monomial_store = internal::monomial_store<MM, V, I>;
  using polynomial_store = internal::polynomial_store<PM, MM, V, I, C>;
};

BOOST_AUTO_TEST_CASE(std_span_behavior) {
  std::vector<uint8_t> v{ 1, 2, 3, 1, 2, 3 };
  std::span<uint8_t> s1{ v.begin(), 3 };
  std::span<uint8_t> s2{ v.begin() + 3, 3 };

  BOOST_TEST(std::equal(s1.begin(), s1.end(), s2.begin(), s2.end()));
}

BOOST_AUTO_TEST_CASE(simple_store) {
  using I = impl<>;

  I::monomial_store s;
  auto m1_idx_1 = s.getid({ 1, 2, 3 });
  auto m1_idx_2 = s.getid({ 1, 2, 3 });

  BOOST_TEST(m1_idx_1 == m1_idx_2);
  BOOST_TEST(m1_idx_1 > 0);

  auto m2_idx_1 = s.getid({ 1, 2, 3, 3 });

  BOOST_TEST(m2_idx_1 != m1_idx_1);

  auto m1_idx_3 = s.getid({ 1, 2, 3 });

  BOOST_TEST(m1_idx_1 == m1_idx_3);
}

BOOST_AUTO_TEST_CASE(simple_product) {
  using I = impl<>;

  I::monomial_store s;

  I::idx m1_idx = s.getid({ 1 });
  I::idx m23_idx = s.getid({ 2, 3 });
  I::idx m123_idx = s.getid({ 1, 2, 3 });
  I::idx m123123_idx = s.getid({ 1, 2, 3, 1, 2, 3 });

  I::idx prod_idx = s.get_product_id(m1_idx, m23_idx);

  BOOST_TEST(m123_idx == prod_idx);

  I::idx prod2_idx = s.get_product_id(m123_idx, m123_idx);
  BOOST_TEST(m123123_idx == prod2_idx);
}

BOOST_AUTO_TEST_CASE(simple_cmp) {
  using I = impl<>;

  I::monomial_store s;

  I::idx m123_idx = s.getid({ 1, 2, 3 });
  I::idx m123_idx2 = s.getid({ 1, 2, 3 });
  I::idx m132_idx = s.getid({ 1, 3, 2 });
  I::idx m123123_idx = s.getid({ 1, 2, 3, 1, 2, 3 });

  BOOST_TEST(!s.cmp(m123_idx, m123_idx2));
  BOOST_TEST(s.cmp(m123_idx, m132_idx));
  BOOST_TEST(!s.cmp(m123123_idx, m123_idx));
}

BOOST_AUTO_TEST_CASE(simple_polynomials) {
  using I = impl<>;
  I::monomial_store s;
  I::polynomial_store p(s);

  I::idx p1 = p.add_polynomial({ { 1l, { 2, 3 } }, { 2l, { 2, 3 } } });
  I::idx p2 = p.add_polynomial({ { 3l, { 3, 4 } }, { 4l, { 5, 6 } } });

  BOOST_TEST(p1 == 1);
  BOOST_TEST(p2 != p1);

  BOOST_TEST(p.get_coefficients(p1)[0].compare(I::coefficient(1l)) == 0);
  BOOST_TEST(p.get_coefficients(p1)[1].compare(I::coefficient(2l)) == 0);
  BOOST_TEST(p.get_coefficients(p2)[0].compare(I::coefficient(3l)) == 0);
  BOOST_TEST(p.get_coefficients(p2)[1].compare(I::coefficient(4l)) == 0);

  I::idx p3 = p.multiply_front(s.getid({ 1, 2 }), p1);
  I::idx p4 = p.multiply_back(p2, s.getid({ 1, 2 }));

  // Printing is also possible:

  //   std::cout << "p1: " << p1 << ", p2: " << p2 << ", p3: " << p3
  //             << ", p4: " << p4 << std::endl;
  //
  //   s.print(std::cout);
  //   p.print(std::cout);

  BOOST_TEST(p.get_monomial_id(p3, 0) == s.getid({ 1, 2, 2, 3 }));
  BOOST_TEST(p.get_monomial_id(p4, 0) == s.getid({ 3, 4, 1, 2 }));
}

template<size_t N>
struct metadata_with_tuple {
  uint8_t length = 0;
  uint16_t t[N];
};

template<size_t N>
void
run(/* some input parameters */) {
  using I
    = impl<metadata_with_tuple<N>, metadata_with_tuple<1>, uint8_t, uint16_t>;

  // Do things with impl.
  typename I::monomial_store s;

  s.getid({ 1, 2, 3 });
}

BOOST_AUTO_TEST_CASE(select_impl_from_number) {
  size_t n = 1;

  // Generate a switch over n with at most 60 variables. The implementation is
  // duplicated internally for each N.
  //
  // The implementation could also be a more complex struct or something else.
  // The run function is only a demonstration.
  boost::mp11::mp_with_index<60>(n, [](auto N) { run<N>(); });
}

BOOST_AUTO_TEST_CASE(simple_with_extended_metadata) {
  struct metadata_with_extra {
    uint8_t length = 0;
    uint16_t m = 1;
  };
  static_assert(sizeof(metadata_with_extra) == 4);
  static_assert(alignof(metadata_with_extra) == 2);

  using I = impl<metadata_with_extra, metadata_with_extra, uint8_t, uint16_t>;
  I::monomial_store s;
  I::idx m1_idx = s.getid({ 1 });
  metadata_with_extra& m1_metadata = s.get_metadata(m1_idx);

  BOOST_TEST(m1_metadata.m == 1);

  I::idx m23_idx = s.getid({ 2, 3 });
  metadata_with_extra& m23_metadata = s.get_metadata(m23_idx);

  BOOST_TEST(m23_metadata.m == 1);

  metadata_with_extra& m1_metadata_ = s.get_metadata(m1_idx);

  BOOST_TEST(m1_metadata.m == 1);
  BOOST_TEST(m1_metadata_.m == 1);
  BOOST_TEST(&m1_metadata == &m1_metadata_);

  I::idx prod_id = s.get_product_id(m1_idx, m23_idx);
  metadata_with_extra& prod_metadata = s.get_metadata(prod_id);
  BOOST_TEST(prod_metadata.m == 1);
  I::idx prod_id_expected = s.getid({ 1, 2, 3 });

  BOOST_TEST(prod_id == prod_id_expected);

  I::idx long_id = s.getid({ 1, 1, 3, 3 });
  metadata_with_extra& long_metadata = s.get_metadata(long_id);
  BOOST_TEST(long_metadata.m == 1);

  std::vector<uint8_t> v(s[long_id].begin(), s[long_id].end());
  BOOST_TEST(v[0] == 1);
  BOOST_TEST(v[1] == 1);
  BOOST_TEST(v[2] == 3);
  BOOST_TEST(v[3] == 3);
}


BOOST_AUTO_TEST_CASE(parse_small_ms_into_polynomial) {
  auto dir_optional = get_test_input_files_location();
  BOOST_REQUIRE(dir_optional.has_value());
  std::filesystem::path dir = *dir_optional;

  std::filesystem::path input = dir / "small.ms";

  using I = impl<>;
  I::monomial_store ms;
  I::polynomial_store ps(ms);

  parser_context ctx;
  ctx.open(input);
  parse_res r = parse_header(ctx);
  BOOST_REQUIRE(!r.has_value());

  BOOST_CHECK(ctx.characteristic == 0);
  BOOST_CHECK(ctx.num_vars == 2);
  BOOST_CHECK(ctx.num_blocks == 1);

  r = parse_rest_into_polynomial_store(ctx, ps);
  BOOST_REQUIRE(!r.has_value());

  std::vector<I::idx> polys(ps.begin(), ps.end());
  BOOST_REQUIRE(polys[0] == 0);
  BOOST_REQUIRE(polys[1] == 1);
  BOOST_REQUIRE(polys[2] == 17);
  std::vector<I::idx> poly1(ps[polys[1]].begin(), ps[polys[1]].end());

  BOOST_TEST(poly1[0] == 1);
}
