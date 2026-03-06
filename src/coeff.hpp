#pragma once

#include <boost/multiprecision/gmp.hpp>

#include <flint/fmpz.h>
namespace f4ncgb {

inline void
fmpz_cleanup(fmpz* arr, size_t len) {
  for(size_t i = 0; i < len; i++)
    fmpz_clear(arr + i);
  delete[] arr;
}

struct coeff {
  fmpz_t value;

  coeff() { fmpz_init(value); }
  ~coeff() { fmpz_clear(value); }

  coeff(const coeff& other) {
    fmpz_init(value);
    fmpz_set(value, other.value);
  }

  coeff(coeff&& other) noexcept {
    fmpz_init(value);
    fmpz_swap(value, other.value);
  }

  // Leaves x in valid but unspecified state
  explicit coeff(fmpz_t&& x) {
    fmpz_init(value);
    fmpz_swap(value, x);
  }

  coeff(long x) {
    fmpz_init(value);
    fmpz_set_si(value, x);
  }

  explicit coeff(const mpz_t x) {
    fmpz_init(value);
    fmpz_set_mpz(value, x);
  }

  friend std::ostream& operator<<(std::ostream& os, const coeff& c) {
    char* s = fmpz_get_str(nullptr, 10, c.value);
    os << s;
    flint_free(s);
    return os;
  }

  bool operator==(const coeff& other) const {
    return fmpz_equal(value, other.value) != 0;
  }

  bool operator!=(const coeff& other) const { return !(*this == other); }

  coeff& operator=(const coeff& other) {
    if(this != &other)
      fmpz_set(value, other.value);
    return *this;
  }

  coeff& operator=(coeff&& other) noexcept {
    if(this != &other)
      fmpz_swap(value, other.value);
    return *this;
  }

  coeff& operator*=(const coeff& other) {
    fmpz_mul(value, value, other.value);
    return *this;
  }

  coeff& operator/=(long b) {
    fmpz_divexact_si(value, value, b);
    return *this;
  }

  friend coeff operator+(const coeff& a, const coeff& b) {
    coeff r;
    fmpz_add(r.value, a.value, b.value);
    return r;
  }

  friend coeff operator+(const coeff& a, ulong b) {
    coeff r;
    fmpz_add_ui(r.value, a.value, b);
    return r;
  }

  friend coeff operator*(const coeff& a, const coeff& b) {
    coeff r;
    fmpz_mul(r.value, a.value, b.value);
    return r;
  }

  friend coeff operator*(const coeff& a, ulong b) {
    coeff r;
    fmpz_mul_ui(r.value, a.value, b);
    return r;
  }

  friend coeff operator*(ulong a, const coeff& b) {
    coeff r;
    fmpz_mul_ui(r.value, b.value, a);
    return r;
  }

  friend coeff operator+(ulong a, const coeff& b) { return b + a; }

  friend bool operator<(const coeff& a, const coeff& b) {
    return fmpz_cmp(a.value, b.value) < 0;
  }

  int sign() const { return fmpz_sgn(value); }

  bool is_zero() const { return fmpz_is_zero(value); }

  bool is_unit() const { return fmpz_is_pm1(value); }

  coeff abs() const {
    coeff r;
    fmpz_abs(r.value, value);
    return r;
  }

  void lcm_inplace(coeff& b) { fmpz_lcm(value, value, b.value); }

  void lcm_inplace(long b) {
    coeff tmp(b);
    lcm_inplace(tmp);
  }
};

}
