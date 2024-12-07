#include <boost/test/unit_test.hpp>

#include "kommunopp.hpp"

using namespace kommunopp;

BOOST_AUTO_TEST_CASE(simple_store) {
  using I=impl<8>;

  I::monomial_store s;
  I::idx i = s.find(m12_1)
}
