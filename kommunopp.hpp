#pragma once

#include <algorithm>
#include <cstdint>
#include <limits>
#include <span>
#include <vector>

#include <boost/multiprecision/gmp.hpp>
#include <boost/unordered/unordered_flat_map.hpp>

namespace kommunopp {
namespace internal {
struct monomial_store_overrun_exception : public std::exception {
  virtual const char* what() const throw() {
    return "monomial_store is overful";
  }
};
struct monomial_length_overrun_exception : public std::exception {
  virtual const char* what() const throw() {
    return "tried to create a monomial that was too long";
  }
};

template<class I, class V>
class monomial_store {
  public:
  using monomial = std::span<V>;

  inline static bool monomial_equal(const monomial& a,
                                    const monomial& b) noexcept {
    return std::equal(a.begin(), a.end(), b.begin(), b.end());
  }

  struct monomial_equal_struct {
    bool operator()(const monomial& a, const monomial& b) const noexcept {
      return monomial_equal(a, b);
    }
  };

  monomial_store() = default;
  ~monomial_store() = default;

  inline V getlength(I idx) const noexcept {
    assert(idx > 0);
    V* len = monomials_.get() + idx;
    return *len;
  }
  inline const monomial operator[](I idx) const {
    assert(idx > 0);
    V len = getlength(idx);
    V* start = monomials_.get() + idx + 1;
    assert(len > 0);
    return monomial(start, len);
  }

  inline const I getidx(const monomial& m) {
    I idx = find(m);
    if(!idx)
      idx = insert(m);
    return idx;
  }
  inline const I getidx(std::vector<V> m) {
    monomial m_span(m.begin(), m.size());
    return getidx(m_span);
  }
  inline const monomial get(const monomial& m) {
    I idx = getidx(m);
    return (*this)[idx];
  }
  inline const monomial get(std::vector<V> m) {
    monomial m_span(m.begin(), m.size());
    return get(m_span);
  }

  inline I getproductidx(I a, I b) {
    std::pair<I, I> prod_pair{ a, b };

    {
      auto it = products_.find(prod_pair);
      if(it != products_.end()) {
        return it->second;
      }
    }

    size_t length_combined = getlength(a);
    length_combined += getlength(b);
    if(length_combined > std::numeric_limits<V>::max()) {
      throw monomial_length_overrun_exception();
    }
    V prod[length_combined];
    auto a_it = (*this)[a];
    auto b_it = (*this)[b];
    auto it = std::copy(a_it.begin(), a_it.end(), prod);
    std::copy(b_it.begin(), b_it.end(), it);
    I prod_idx = getidx(monomial(prod, length_combined));
    products_.insert(std::make_pair(prod_pair, prod_idx));
    return prod_idx;
  }
  inline const monomial getproduct(I a, I b) {
    return (*this)[getproductidx(a, b)];
  }

  private:
  I capacity_ = std::numeric_limits<I>::max();
  I size_ = 1;
  std::unique_ptr<V[]> monomials_
    = std::make_unique_for_overwrite<V[]>(capacity_);

  boost::unordered_flat_map<monomial,
                            I,
                            boost::hash<monomial>,
                            monomial_equal_struct>
    map_;

  boost::unordered_flat_map<std::pair<I, I>, I> products_;

  I find(const monomial& m) const {
    auto it = map_.find(m);
    if(it == map_.end())
      return 0;
    return it->second;
  }
  I insert(const monomial& m) {
    if(capacity_ - m.size() - 1 < size_) {
      throw monomial_store_overrun_exception();
    }
    auto idx = size_;
    *(monomials_.get() + size_++) = m.size();
    std::copy(m.begin(), m.end(), monomials_.get() + size_);
    map_.insert(std::pair((*this)[idx], idx));
    size_ += m.size();
    return idx;
  }
};

template<class I, class V>
class polynomial {
  public:
  using monomial_store = monomial_store<I, V>;
  using monomial = monomial_store::monomial;

  using rational = boost::multiprecision::mpq_rational;

  inline polynomial(monomial_store& store)
    : store_(store) {}

  inline I append(I idx, rational k = 1) {
    monomials_.push_back(idx);
    koefficients_.push_back(k);
    return idx;
  }
  inline I append(monomial v, rational k = 1) {
    I idx = store_.getidx(v);
    monomials_.push_back(idx);
    koefficients_.push_back(k);
    return idx;
  }
  inline I append(std::vector<V> v, rational k = 1) {
    I idx = store_.getidx(v);
    monomials_.push_back(idx);
    koefficients_.push_back(k);
    return idx;
  }

  inline void multiply_back(I idx) {
    for(I& m : monomials_) {
      m = store_.getproductidx(m, idx);
    }
  }
  inline void multiply_front(I idx) {
    for(I& m : monomials_) {
      m = store_.getproductidx(idx, m);
    }
  }

  inline I getidx(I idx) const {
    assert(idx < monomials_.size());
    return monomials_[idx];
  }
  inline rational getkoeff(I idx) const {
    assert(idx < koefficients_.size());
    return koefficients_[idx];
  }
  inline monomial get(I idx) const {
    assert(idx < monomials_.size());
    return store_[monomials_[idx]];
  }
  inline I size() const { return monomials_.size(); }

  const std::vector<I> monomials() const { return monomials_; }

  private:
  monomial_store& store_;
  // List of indices into the global monomial vector.
  std::vector<I> monomials_;
  std::vector<rational> koefficients_;
};
}

template<typename I = uint32_t, typename V = uint8_t>
struct impl {
  using idx = I;
  using var = V;
  using monomial_store = internal::monomial_store<I, V>;
  using polynomial = internal::polynomial<I, V>;
};
}
