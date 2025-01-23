#pragma once

#include <cstddef>
#include <cstdint>

// The operations in this file are inspired by:
//
//   - [1]: https://arxiv.org/pdf/2008.08654

namespace kommunopp {
[[nodiscard]] constexpr inline uint32_t
v_mod_2_31_1(uint32_t v) noexcept {
  // Algorithm 4 of [1]
  //
  // Mersenne Prime: 2^31-1
  // b = 31

  const uint32_t v_prime = v + 1;
  const uint32_t z = ((v_prime >> 31) + v_prime) >> 31;
  return (v + z) & 2'147'483'647;// 2^31-1
}

[[nodiscard]] constexpr inline uint32_t
v_mod_p(uint32_t v, uint32_t p) noexcept {
  // Algorithm 5 of [1]
  const uint64_t c = 2'147'483'648 - p;// 2^31 - p = c; b = 31
  const uint64_t b = 31;
  const uint64_t v_prime = (uint64_t)v + (uint64_t)c;
  uint64_t z = v_prime >> b;
  z = (z * c + v_prime) >> b;
  z = (z * c + v_prime) >> b;
  return (v - z * p);
}

[[nodiscard]] constexpr inline uint32_t
mult_mod_2_31_1(uint32_t a, uint32_t b, uint32_t x) noexcept {
  // Algorithm 4 of [1], merged with multiply and addition
  //
  // Mersenne Prime: 2^31-1
  // b = 31
  const uint64_t ax = (uint64_t)a * (uint64_t)x;
  const uint64_t ax_b_1 = ax + b + 1ul;
  const uint64_t z = ((ax_b_1 >> 31) + ax_b_1) >> 31;
  return (ax_b_1 + z - 1) & 2'147'483'647ul;// 2^31-1
}

[[nodiscard]] constexpr inline uint32_t
mult_mod_p(uint32_t a, uint32_t b, uint32_t x, uint32_t p) noexcept {
  // Algorithm 5 of [1] together with multiplication.
  return ((uint64_t)a * (uint64_t)x + (uint64_t)b) % p;
}

}
