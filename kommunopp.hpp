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
template<class B,
         metadata_concept M = metadata<>,
         value_concept V = uint8_t,
         typename I = uint32_t>
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

  friend B;

  protected:
  using self = store<B, M, V, I>;

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
  I inserted_count_ = 0;
  M zero_metadata_;
  bool scratch_metadata_created_ = false;

  std::unique_ptr<std::byte[]> pool_
    = std::make_unique_for_overwrite<std::byte[]>(capacity());

  /// Work on the current tip but do not commit anything. The tip can later be
  /// committed using insert_scratch.
  inline std::pair<M&, V*> new_scratch(I space = 32) {
    if(capacity() - space * sizeof(V) - sizeof(M) < size_) {
      throw monomial_store_overrun_exception();
    }

    void* ptr = boost::alignment::align_up(pool_.get() + size_, alignof(M));

    // Destruct old metadata at that position if it was created previously.
    if(scratch_metadata_created_)
      reinterpret_cast<M*>(std::assume_aligned<alignof(M)>(ptr))->~M();

    M* metadata = new(std::assume_aligned<alignof(M)>(ptr)) M;
    scratch_metadata_created_ = true;

    // Advance the pointer to the first possible location of V.
    ptr = boost::alignment::align_up(static_cast<std::byte*>(ptr) + sizeof(M),
                                     alignof(V));
    V* vv = reinterpret_cast<V*>(std::assume_aligned<alignof(V)>(ptr));

    return std::pair<M&, V*>(*metadata, vv);
  }

  inline I insert_scratch() {
    I id = size_;

    void* ptr = boost::alignment::align_up(pool_.get() + size_, alignof(M));

    M* metadata = reinterpret_cast<M*>(std::assume_aligned<alignof(M)>(ptr));
    assert(metadata->length > 0);

    std::byte* ptr_start = static_cast<std::byte*>(ptr);
    ptr = boost::alignment::align_up(static_cast<std::byte*>(ptr) + sizeof(M),
                                     alignof(V));

    std::byte* ptr_end
      = static_cast<std::byte*>(ptr) + metadata->length * sizeof(V);

    reinterpret_cast<B*>(this)->new_entry(id);

    size_ += ptr_end - ptr_start;
    ++inserted_count_;

    scratch_metadata_created_ = false;

    return id + 1;
  }

  inline I insert(const std::span<const V>& v) {
    auto [m, vv] = new_scratch(v.size());
    std::copy(v.begin(), v.end(), vv);
    m.length = v.size();
    return insert_scratch();
  }

  inline void new_entry(I id) { (void)id; }

  public:
  using length_type = decltype(M::length);
  using value_type = V;

  store() = default;
  ~store() {
    for(I i = 0; i < inserted_count_; ++i) {
      // Do not call the destructor of the 0 element, as this is special. Only
      // call higher ones.
      M& m = get_metadata(i + 1);
      m.~M();
    }

    // If there was some scratch space with metadata, destroy it too.
    if(scratch_metadata_created_)
      get_metadata(inserted_count_ + 1).~M();
  }

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
};

template<metadata_concept M = metadata<>,
         value_concept V = uint8_t,
         typename I = uint32_t>
class monomial_store : public store<monomial_store<M, V, I>, M, V, I> {
  using self = monomial_store<M, V, I>;
  using base = store<self, M, V, I>;
  friend base;

  boost::unordered_flat_map<std::span<const V>,
                            I,
                            boost::hash<std::span<const V>>,
                            typename base::V_equality_struct>
    map_;

  boost::unordered_flat_map<std::tuple<I, I, I>, I> products_;

  protected:
  inline void new_entry(I id) { map_.insert(std::pair((*this)[id + 1], id)); }

  I find(const std::span<const V>& v) const {
    if(v.size() == 0)
      return 0;

    const auto it = map_.find(v);
    if(it == map_.end())
      return 0;

    return it->second + 1;
  }

  public:
  inline const I getid(const std::span<const V>& v) {
    I id = find(v);
    if(!id)
      id = base::insert(v);
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
    m.length = length_combined;
    auto a_it = (*this)[a];
    auto b_it = (*this)[b];
    auto it = std::copy(a_it.begin(), a_it.end(), vv);
    std::copy(b_it.begin(), b_it.end(), it);
    I prod_idx = find(std::span(vv, length_combined));
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
    m.length = length_combined;
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

template<metadata_concept M, typename I>
struct polynomial_metadata : public M {
  I coefficients;
};

template<metadata_concept PM,
         metadata_concept MM,
         value_concept V,
         typename I,
         typename C = boost::multiprecision::gmp_rational>
class polynomial_store
  : public store<polynomial_store<PM, MM, V, I, C>,
                 polynomial_metadata<PM, I>,
                 I,
                 I> {
  public:
  using self = polynomial_store<PM, MM, V, I, C>;
  using base = store<polynomial_store<PM, MM, V, I, C>,
                     polynomial_metadata<PM, I>,
                     I,
                     I>;
  using monomial_store_ = monomial_store<MM, V, I>;
  using monomial = monomial_store_::value_type;
  using metadata = polynomial_metadata<PM, I>;

  using polynomial_vec = std::vector<std::pair<C, I>>;
  using polynomial_nested_vec = std::vector<std::pair<C, std::vector<V>>>;

  friend base;

  inline polynomial_store(monomial_store_& store)
    : base::store()
    , store_(store) {}

  ~polynomial_store() {
    for(I i = 0; i < cpool_size_; ++i) {
      C* c = get_coefficients_raw(i);
      // Explicitly destruct the coefficent again.
      c->~C();
    }
  }

  /// Initialize a new polynomial. Remember to initialize the coefficents!
  inline std::tuple<metadata&, I*, C*> add(I length = 32) {
    auto [m, vv] = this->new_scratch(length);
    return std::tuple<metadata&, I*, C*>(
      reinterpret_cast<metadata&>(m), vv, get_coefficients_raw(cpool_size_));
  }
  inline I commit() { return this->insert_scratch(); }

  inline I add_polynomial(polynomial_vec p) {
    auto [m, monomials, coefficients] = add(p.size());
    for(size_t i = 0; i < p.size(); ++i) {
      C* c = new(coefficients + i) C;
      *c = p[i].first;
      monomials[i] = p[i].second;
    }
    m.length = p.size();
    return commit();
  }
  inline I add_polynomial(polynomial_nested_vec p) {
    auto [m, monomials, coefficients] = add(p.size());
    for(size_t i = 0; i < p.size(); ++i) {
      C* c = new(coefficients + i) C;
      *c = p[i].first;
      monomials[i] = store_.getid(p[i].second);
    }
    m.length = p.size();
    return commit();
  }

  std::span<C> get_coefficients(I id) {
    if(id == 0)
      return std::span<C>();
    metadata& m = this->get_metadata(id);
    return std::span<C>(get_coefficients_raw(m.coefficients), m.length);
  }

  protected:
  monomial_store_& store_;
  std::unique_ptr<std::byte[]> cpool_
    = std::make_unique_for_overwrite<std::byte[]>(
      std::numeric_limits<I>::max());
  size_t cpool_size_ = 0;

  inline C* get_coefficients_raw(I id) {
    return reinterpret_cast<C*>(cpool_.get()
                                + static_cast<size_t>(id) * sizeof(C));
  }

  inline void new_entry(I id) {
    // Commit the new coefficients to the coefficient array.
    metadata& m = this->get_metadata(id + 1);
    m.coefficients = cpool_size_;
    cpool_size_ += m.length;
  }
};
}

template<internal::metadata_concept MM = internal::metadata<uint8_t>,
         internal::metadata_concept PM = internal::metadata<uint8_t>,
         internal::value_concept V = uint8_t,
         typename I = uint32_t,
         typename C = boost::multiprecision::gmp_rational>
struct impl {
  using var = V;
  using idx = I;
  using coefficient = C;
  using monomial_store = internal::monomial_store<MM, V, I>;
  using polynomial_store = internal::polynomial_store<PM, MM, V, I, C>;
};
}
