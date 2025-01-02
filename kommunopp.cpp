#include "kommunopp.hpp"

#include <boost/multiprecision/detail/default_ops.hpp>
#include <boost/multiprecision/gmp.hpp>

namespace kommunopp {
int
coefficient_sign(const boost::multiprecision::backends::gmp_rational& r) {
  using namespace boost::multiprecision;
  mpq_rational rr(r);
  return boost::multiprecision::sign(rr);
}

boost::multiprecision::mpq_rational
coefficient_abs(const boost::multiprecision::backends::gmp_rational& r) {
  using namespace boost::multiprecision;
  mpq_rational rr(r);
  return boost::multiprecision::abs(rr);
}

bool
coefficient_is_posneg_neutral(
  const boost::multiprecision::backends::gmp_rational& r) {
  using namespace boost::multiprecision;
  mpq_rational rr(r);
  return (numerator(rr) == -1l || numerator(rr) == 1l) && denominator(rr) == 1l;
}
}

std::ostream&
operator<<(std::ostream& o,
           const boost::multiprecision::backends::gmp_rational& r) {
  using namespace boost::multiprecision;
  mpq_rational rr(r);
  return o << rr;
}
