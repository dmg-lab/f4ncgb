#pragma once

#include <algorithm>
#include <concepts>
#include <cstdint>
#include <limits>
#include <span>
#include <vector>

#include <boost/align/align_down.hpp>
#include <boost/align/align_up.hpp>
#include <boost/container/small_vector.hpp>
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

template<typename length_type = uint8_t>
struct metadata {
  length_type length = 0;
};

static_assert(sizeof(metadata<uint8_t>) == sizeof(uint8_t));
static_assert(alignof(metadata<uint8_t>) == alignof(uint8_t));

template<typename T>
concept metadata_concept = requires(T t) {
  { t.length } -> std::convertible_to<size_t>;
};

consteval bool
is_power_of_two(size_t n) {
  return (n & (n - 1)) == 0;
}

template<typename T>
concept value_concept = requires { is_power_of_two(sizeof(T)); };

/// A store holds some instances of structures like [metadata, V_1, ...
/// V_metadata.length].
template<metadata_concept M = metadata<>,
         value_concept V = uint8_t,
         typename I = uint32_t>
  requires(sizeof(M) >= sizeof(V))
class store {
  public:
  inline static bool V_equal(const std::span<const V>& a,
                             const std::span<const V>& b) noexcept {
    return std::equal(a.begin(), a.end(), b.begin(), b.end());
  }

  struct V_equality_struct {
    bool operator()(const std::span<const V>& a,
                    const std::span<const V>& b) const noexcept {
      return V_equal(a, b);
    }
  };

  protected:
  using self = store<M, V, I>;

  consteval static size_t capacity() {
    return std::min(std::numeric_limits<I>::max() * alignof(M),
                    std::numeric_limits<size_t>::max());
  }

  constexpr inline M& get_metadata_from_id(I id) {
    return *reinterpret_cast<M*>(std::assume_aligned<alignof(M)>(
      boost::alignment::align_up(pool_.get() + id, alignof(M))));
  }
  constexpr inline M& get_metadata_from_id(I id) const {
    return const_cast<self*>(this)->get_metadata_from_id(id);
  }
  constexpr inline V* get_value_from_id(I id) {
    void* ptr = boost::alignment::align_up(pool_.get() + id, alignof(M));
    return reinterpret_cast<V*>(
      std::assume_aligned<alignof(V)>(boost::alignment::align_up(
        static_cast<std::byte*>(ptr) + sizeof(M), alignof(V))));
  }
  constexpr inline const V* get_value_from_id(I id) const {
    return const_cast<self*>(this)->get_value_from_id(id);
  }

  I size_ = 0;
  M zero_metadata_;

  std::unique_ptr<std::byte[]> pool_
    = std::make_unique_for_overwrite<std::byte[]>(capacity());

  boost::unordered_flat_map<std::span<const V>,
                            I,
                            boost::hash<std::span<const V>>,
                            V_equality_struct>
    map_;

  I find(const std::span<const V>& v) const {
    if(v.size() == 0)
      return 0;

    const auto it = map_.find(v);
    if(it == map_.end())
      return 0;

    return it->second + 1;
  }

  /// Work on the current tip but do not commit anything. The tip can later be
  /// committed using insert_scratch.
  inline std::pair<M&, V*> new_scratch(I space = 32) {
    if(capacity() - space * sizeof(V) - sizeof(M) < size_) {
      throw monomial_store_overrun_exception();
    }

    void* ptr = boost::alignment::align_up(pool_.get() + size_, alignof(M));

    M* metadata = new(std::assume_aligned<alignof(M)>(ptr)) M;

    // Advance the pointer to the first possible location of V.
    ptr = boost::alignment::align_up(static_cast<std::byte*>(ptr) + sizeof(M),
                                     alignof(V));
    V* vv = reinterpret_cast<V*>(std::assume_aligned<alignof(V)>(ptr));

    return std::pair<M&, V*>(*metadata, vv);
  }

  I insert_scratch() {
    I id = size_;

    void* ptr = boost::alignment::align_up(pool_.get() + size_, alignof(M));

    M* metadata = reinterpret_cast<M*>(std::assume_aligned<alignof(M)>(ptr));
    assert(metadata->length > 0);

    std::byte* ptr_start = static_cast<std::byte*>(ptr);
    ptr = boost::alignment::align_up(static_cast<std::byte*>(ptr) + sizeof(M),
                                     alignof(V));

    std::byte* ptr_end
      = static_cast<std::byte*>(ptr) + metadata->length * sizeof(V);

    map_.insert(std::pair((*this)[id + 1], id));

    size_ += ptr_end - ptr_start;

    return id + 1;
  }

  I insert(const std::span<const V>& v) {
    auto [m, vv] = new_scratch(v.size());
    std::copy(v.begin(), v.end(), vv);
    m.length = v.size();
    return insert_scratch();
  }

  public:
  using length_type = decltype(M::length);
  using value_type = V;

  store() = default;
  ~store() = default;

  inline M& get_metadata(I id) noexcept {
    if(id == 0) {
      return zero_metadata_;
    }
    return get_metadata_from_id(id - 1);
  }
  inline const M& get_metadata(I id) const noexcept {
    return const_cast<self*>(this)->get_metadata(id);
  }

  inline length_type get_length(I id) const noexcept {
    return get_metadata(id).length;
  }

  inline const std::span<const V> operator[](I id) const noexcept {
    if(id == 0)
      return std::span<V>();

    size_t len = get_length(id);

    const V* start = get_value_from_id(id - 1);

    return std::span<const V>(start, len);
  }

  inline const I getid(const std::span<const V>& v) {
    I id = find(v);
    if(!id)
      id = insert(v);
    return id;
  }
  inline const I getid(std::vector<V> v) {
    std::span<const V> s(v.begin(), v.size());
    return getid(s);
  }

  inline const std::span<V> get(const std::span<V>& v) {
    I id = getid(v);
    return (*this)[id];
  }
  inline const std::span<V> get(std::vector<V> v) {
    std::span<V> s(v.begin(), v.size());
    return get(s);
  }
};

template<metadata_concept M = metadata<>,
         value_concept V = uint8_t,
         typename I = uint32_t>
class monomial_store : public store<M, V, I> {
  using base = store<M, V, I>;

  boost::unordered_flat_map<std::tuple<I, I, I>, I> products_;

  public:
  inline I get_product_id(I a, I b) {
    std::tuple<I, I, I> prod_tuple{ a, b, 0 };

    {
      auto it = products_.find(prod_tuple);
      if(it != products_.end()) {
        return it->second;
      }
    }

    const size_t length_combined
      = static_cast<size_t>(base::get_length(a)) + base::get_length(b);
    if(length_combined
       > std::numeric_limits<typename base::length_type>::max()) {
      throw monomial_length_overrun_exception();
    }
    auto [m, vv] = base::new_scratch(length_combined);
    auto a_it = (*this)[a];
    auto b_it = (*this)[b];
    auto it = std::copy(a_it.begin(), a_it.end(), vv);
    std::copy(b_it.begin(), b_it.end(), it);
    I prod_idx = base::find(std::span(vv, length_combined));
    if(!prod_idx) {
      prod_idx = base::insert_scratch();
    }
    products_.insert(std::make_pair(prod_tuple, prod_idx));
    return prod_idx;
  }

  inline I get_product_id(I a, I b, I c) {
    std::tuple<I, I, I> prod_tuple{ a, b, c };

    {
      auto it = products_.find(prod_tuple);
      if(it != products_.end()) {
        return it->second;
      }
    }

    const size_t length_combined = static_cast<size_t>(base::get_length(a))
                                   + base::get_length(b) + base::get_length(c);
    if(length_combined
       > std::numeric_limits<typename base::length_type>::max()) {
      throw monomial_length_overrun_exception();
    }
    auto [m, vv] = base::new_scratch(length_combined);
    auto a_it = (*this)[a];
    auto b_it = (*this)[b];
    auto c_it = (*this)[c];
    auto it = std::copy(a_it.begin(), a_it.end(), vv);
    std::copy(b_it.begin(), b_it.end(), it);
    std::copy(c_it.begin(), c_it.end(), it);
    I prod_idx = base::find(std::span(vv, length_combined));
    if(!prod_idx) {
      prod_idx = base::insert_scratch();
    }
    products_.insert(std::make_pair(prod_tuple, prod_idx));
    return prod_idx;
  }

  inline const std::span<V> get_product(I a, I b) {
    return (*this)[get_product_id(a, b)];
  }
};

template<metadata_concept M, value_concept V, typename I>
class polynomial {
  public:
  using monomial_store = monomial_store<M, V, I>;
  using monomial = monomial_store::value_type;

  using rational = boost::multiprecision::mpq_rational;

  inline polynomial(monomial_store& store)
    : store_(store) {}

  inline I append(I idx, rational k = 1) {
    monomials_.push_back(idx);
    koefficients_.push_back(k);
    return idx;
  }
  inline I append(monomial v, rational k = 1) {
    I idx = store_.getid(v);
    monomials_.push_back(idx);
    koefficients_.push_back(k);
    return idx;
  }
  inline I append(std::vector<V> v, rational k = 1) {
    I idx = store_.getid(v);
    monomials_.push_back(idx);
    koefficients_.push_back(k);
    return idx;
  }

  inline void multiply_back(I idx) {
    for(I& m : monomials_) {
      m = store_.get_product_id(m, idx);
    }
  }
  inline void multiply_front(I idx) {
    for(I& m : monomials_) {
      m = store_.get_product_id(idx, m);
    }
  }

  inline I getid(I id) const {
    assert(id < monomials_.size());
    return monomials_[id];
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

  inline const std::vector<I> monomials() const { return monomials_; }

  private:
  monomial_store& store_;
  // List of indices into the global monomial vector.
  std::vector<I> monomials_;
  std::vector<rational> koefficients_;
};
}

template<internal::metadata_concept M = internal::metadata<uint8_t>,
         internal::value_concept V = uint8_t,
         typename I = uint32_t>
struct impl {
  using var = V;
  using idx = I;
  using monomial_store = internal::monomial_store<M, V, I>;
  using polynomial = internal::polynomial<M, V, I>;
};
}
