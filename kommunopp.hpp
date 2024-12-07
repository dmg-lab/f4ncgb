#pragma once

#include <algorithm>
#include <cstdint>
#include <limits>
#include <ostream>
#include <vector>

#include <boost/unordered/unordered_flat_map.hpp>

namespace kommunopp {
namespace internal {
template<typename V = uint8_t>
class monomial {
  public:
  using var = V;
  using self = monomial<V>;

  private:
};

template<class I, class M>
class monomial_store {
  public:
  using monomial = M;

  const monomial& operator[](I idx) const {
    assert(idx > 0);
    return (monomials_.get())[idx];
  }

  private:
  I capacity_ = std::numeric_limits<I>::max();
  std::unique_ptr<typename monomial::var[]> monomials_;
  boost::unordered_flat_map<std::reference_wrapper<monomial>, I> map_;

  I find(const monomial& m) {
  }
};

template<class I, class M>
class polynomial {
  public:
  using monomial = M;

  void add(I idx);

  private:
  // List of indices into the global monomial vector.
  std::vector<I> monomials;
};
}

template<typename I = uint32_t, typename V = uint8_t>
struct impl {
  using idx = I;
  using monomial = internal::monomial<V>;
  using monomial_store = internal::monomial_store<I, monomial>;
  using polynomial = internal::polynomial<I, monomial>;
};
}
