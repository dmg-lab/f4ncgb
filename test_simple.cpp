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
  auto m1_idx_1 = s.getidx({ 1, 2, 3 });
  auto m1_idx_2 = s.getidx({ 1, 2, 3 });

  BOOST_TEST(m1_idx_1 == m1_idx_2);
  BOOST_TEST(m1_idx_1 > 0);

  auto m2_idx_1 = s.getidx({ 1, 2, 3, 3 });

  BOOST_TEST(m2_idx_1 != m1_idx_1);

  auto m1_idx_3 = s.getidx({ 1, 2, 3 });

  BOOST_TEST(m1_idx_1 == m1_idx_3);
}

BOOST_AUTO_TEST_CASE(simple_product) {
  using I = impl<>;

  I::monomial_store s;

  I::idx m1_idx = s.getidx({ 1 });
  I::idx m23_idx = s.getidx({ 2, 3 });
  I::idx m123_idx = s.getidx({ 1, 2, 3 });
  I::idx m123123_idx = s.getidx({ 1, 2, 3, 1, 2, 3 });

  I::idx prod_idx = s.getproductidx(m1_idx, m23_idx);

  BOOST_TEST(m123_idx == prod_idx);

  I::idx prod2_idx = s.getproductidx(m123_idx, m123_idx);
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

  I::idx m12 = s.getidx({ 1, 2 });
  I::idx m112 = s.getidx({ 1, 1, 2 });
  I::idx m1212 = s.getidx({ 1, 2, 1, 2 });
  I::idx m121 = s.getidx({ 1, 2, 1 });

  I::idx p1m0 = p1.getidx(0);

  BOOST_TEST(m12 == p1m0);

  p1.multiply_back(m12);
  p2.multiply_front(m12);

  p1m0 = p1.getidx(0);
  I::idx p1m1 = p1.getidx(1);

  BOOST_TEST(m1212 == p1m0);
  BOOST_TEST(m112 == p1m1);

  I::idx p2m0 = p2.getidx(0);
  I::idx p2m1 = p2.getidx(1);

  BOOST_TEST(m1212 == p2m0);
  BOOST_TEST(m121 == p2m1);
}
