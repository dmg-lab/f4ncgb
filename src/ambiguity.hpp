#ifndef AMBIGUITY_H
#define AMBIGUITY_H

#include <algorithm>
#include <array>
#include <boost/functional/hash.hpp>
#include <boost/unordered/unordered_map.hpp>
#include <span>

namespace f4ncgb {

template<typename I>
size_t inline compute_hash(const std::array<I, 7>& values) noexcept {
  size_t seed = 0;
  for(I value : values)
    boost::hash_combine(seed, value);
  return seed;
}

// container to store when two monomials overlaps
template<typename I>
struct ambiguity {
  std::array<I, 7> values;
  size_t hash_;
  mutable I fi_;
  mutable I fj_;
  mutable I lm_;

  using mon_id = I;
  using poly_id = I;

  inline ambiguity() = default;

  inline ambiguity(I d, I i, I j, I ai, I ci, I aj, I cj)
    : values{ d, i, j, ai, ci, aj, cj }
    , hash_(compute_hash(values)) {}

  inline I degree() const noexcept { return values[0]; }
  inline I i() const noexcept { return values[1]; }
  inline I j() const noexcept { return values[2]; }
  inline I ai() const noexcept { return values[3]; }
  inline I ci() const noexcept { return values[4]; }
  inline I aj() const noexcept { return values[5]; }
  inline I cj() const noexcept { return values[6]; }
  inline I fi() const noexcept { return fi_; }
  inline I fj() const noexcept { return fj_; }
  inline I lm() const noexcept { return lm_; }

  // Equality operator for unordered_set
  inline bool operator==(const ambiguity& other) const noexcept {
    return (hash_ == other.hash_) and (values == other.values);
  }

  template<typename PS>
  void prepare(PS& poly,
               boost::unordered_map<mon_id, poly_id>& lm_to_poly) const {
    poly_id ii = lm_to_poly[i()];
    poly_id jj = lm_to_poly[j()];

    fi_ = poly.multiply_front_and_back(ai(), ii, ci());
    fj_ = poly.multiply_front_and_back(aj(), jj, cj());
    lm_ = poly.get_lm_id(values[7]);
  }
};

// Custom hash function
template<typename I>
struct ambiguity_hash {
  inline size_t operator()(const ambiguity<I>& a) const noexcept {
    return a.hash_;
  }
};

template<class T, std::size_t N, std::size_t M>
constexpr inline bool
startswith(std::span<T, N> data, std::span<T, M> prefix) noexcept {
  return data.size() >= prefix.size()
         && std::equal(prefix.begin(), prefix.end(), data.begin());
}

template<class T, std::size_t N, std::size_t M>
constexpr inline bool
endswith(std::span<T, N> data, std::span<T, M> suffix) noexcept {
  return data.size() >= suffix.size()
         && std::equal(data.end() - static_cast<long>(suffix.size()),
                       data.end(),
                       suffix.end() - static_cast<long>(suffix.size()));
}
}

#endif// AMBIGUITY_H
