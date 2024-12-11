#include <boost/test/unit_test.hpp>

#include "kommunopp.hpp"

using namespace kommunopp;

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

BOOST_AUTO_TEST_CASE(simple_polynomial) {
  using I = impl<>;

  I::monomial_store s;

  I::polynomial p1(s);
  p1.append({ 1, 2 });
  p1.append(std::initializer_list<I::var>{ 1 });

  I::polynomial p2(s);
  p2.append({ 1, 2 });
  p2.append(std::initializer_list<I::var>{ 1 });

  I::idx m12 = s.getid({ 1, 2 });
  I::idx m112 = s.getid({ 1, 1, 2 });
  I::idx m1212 = s.getid({ 1, 2, 1, 2 });
  I::idx m121 = s.getid({ 1, 2, 1 });

  I::idx p1m0 = p1.getid(0);

  BOOST_TEST(m12 == p1m0);

  p1.multiply_back(m12);
  p2.multiply_front(m12);

  p1m0 = p1.getid(0);
  I::idx p1m1 = p1.getid(1);

  BOOST_TEST(m1212 == p1m0);
  BOOST_TEST(m112 == p1m1);

  I::idx p2m0 = p2.getid(0);
  I::idx p2m1 = p2.getid(1);

  BOOST_TEST(m1212 == p2m0);
  BOOST_TEST(m121 == p2m1);
}

BOOST_AUTO_TEST_CASE(simple_with_extended_metadata) {
  struct metadata_with_extra {
    uint8_t length = 0;
    uint8_t extra1 = 1;
    uint8_t extra2 = 2;
  };
  static_assert(sizeof(metadata_with_extra) == 3);
  static_assert(alignof(metadata_with_extra) == 1);

  using I = impl<metadata_with_extra, uint8_t, uint16_t>;
  I::monomial_store s;
  I::idx m1_idx = s.getid({ 1 });
  metadata_with_extra& m1_metadata = s.get_metadata(m1_idx);

  BOOST_TEST(m1_metadata.extra1 == 1);
  BOOST_TEST(m1_metadata.extra2 == 2);

  I::idx m23_idx = s.getid({ 2, 3 });
  metadata_with_extra& m23_metadata = s.get_metadata(m23_idx);

  BOOST_TEST(m23_metadata.extra1 == 1);
  BOOST_TEST(m23_metadata.extra2 == 2);

  metadata_with_extra& m1_metadata_ = s.get_metadata(m1_idx);

  BOOST_TEST(m1_metadata.extra1 == 1);
  BOOST_TEST(m1_metadata.extra2 == 2);
  BOOST_TEST(m1_metadata_.extra1 == 1);
  BOOST_TEST(m1_metadata_.extra2 == 2);
  BOOST_TEST(&m1_metadata == &m1_metadata_);

  I::idx prod_id = s.get_product_id(m1_idx, m23_idx);
  metadata_with_extra& prod_metadata = s.get_metadata(prod_id);
  BOOST_TEST(prod_metadata.extra1 == 1);
  BOOST_TEST(prod_metadata.extra2 == 2);
  I::idx prod_id_expected = s.getid({ 1, 2, 3 });

  BOOST_TEST(prod_id == prod_id_expected);
}
