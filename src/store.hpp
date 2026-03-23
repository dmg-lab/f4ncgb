#pragma once

#include <algorithm>
#include <concepts>
#include <cstdint>
#include <cstdlib>
#include <format>
#include <fstream>
#include <iterator>
#include <limits>
#include <memory>
#include <new>
#include <span>
#include <utility>
#include <vector>

// The profiler shows boost::alignment to take a while, but this seems
// to be an artifact from waiting for memory to appear. It makes no
// difference if this is changed, so it will remain at this approch.
#include <boost/align/align_down.hpp>
#include <boost/align/align_up.hpp>

#include <boost/align/aligned_alloc.hpp>

#include <boost/container/small_vector.hpp>

#include <boost/container_hash/hash.hpp>
#include <boost/unordered/unordered_flat_map.hpp>
#include <boost/unordered/unordered_flat_set.hpp>

#include "coeff.hpp"
#include "debug.hpp"
#include "f4ncgb.hpp"
#include "profiling.hpp"
#include "signal_statistics.hpp"

#include "ambiguity.hpp"

#ifdef F4NCGB_ENABLE_STORE_TRACE
#include "store_tracer.hpp"
#endif

#ifdef __linux__
#include <sys/mman.h>
#endif

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
// Older GCC compilers warn about unknown warnings, which is really annoying.
#pragma GCC diagnostic ignored "-Wpragmas"
#pragma GCC diagnostic ignored "-Walloc-size"
#endif

namespace f4ncgb {

extern std::vector<coeff> input_denoms;
extern size_t proof_level;
extern bool interreduce;
extern bool reduce_mode;
extern coeff reduce_denom_c;
extern size_t nr_cols;

namespace internal {
struct monomial_store_overrun_exception : public std::exception {
  virtual const char* what() const throw() {
    return "monomial_store is overful";
  }
};
struct length_overrun_exception : public std::exception {
  size_t size;
  size_t capacity;
  std::string msg;
  length_overrun_exception(size_t size, size_t capacity)
    : size(size)
    , capacity(capacity)
    , msg(std::format("tried to create a length of {} but the highest possible "
                      "fitting length is {}",
                      size,
                      capacity)) {}
  virtual const char* what() const throw() { return msg.c_str(); }
};
struct scratch_insertion_with_zero_length_exception : public std::exception {
  size_t size;
  std::string msg;
  scratch_insertion_with_zero_length_exception(size_t size)
    : size(size)
    , msg(std::format(
        "Tried to call insert_scratch() on a store where the current metadata "
        "wasn't assigned a length, size of the store was {}",
        size)) {}
  virtual const char* what() const throw() { return msg.c_str(); }
};
struct coefficient_overrun_exception : public std::exception {
  size_t size;
  size_t capacity;
  std::string msg;
  coefficient_overrun_exception(size_t size, size_t capacity)
    : size(size)
    , capacity(capacity)
    , msg(std::format("tried to get a coefficient {} the capacity is only {}",
                      size,
                      capacity)) {}
  virtual const char* what() const throw() { return msg.c_str(); }
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
  protected:
  using self = store<B, M, V, I>;
  using V_span = std::span<const V>;
  using V_span_it = V_span::iterator;

  struct pool_deleter {
    void operator()(std::byte pool[]) { free(pool); }
  };

  public:
  __attribute__((always_inline)) inline static bool V_equal(
    const std::span<const V>& a,
    const std::span<const V>& b) noexcept {
    return std::equal(a.begin(), a.end(), b.begin(), b.end());
  }

  struct V_equality_struct {
    explicit V_equality_struct(self& s)
      : s(s) {}

    self& s;

    __attribute__((always_inline)) inline bool operator()(I a,
                                                          I b) const noexcept {
      const auto& a_span = s[a];
      const auto& b_span = s[b];
      return V_equal(a_span, b_span);
    }
  };

  struct V_hash_struct : public boost::hash<I> {
    explicit V_hash_struct(self& s)
      : s(s) {}

    self& s;

    typedef I argument_type;
    typedef std::size_t result_type;

    __attribute__((always_inline)) inline bool operator()(I a) const noexcept {
      const auto a_span = s[a];
      return boost::hash_value(a_span);
    }
  };

  struct V_span_equality_struct {
    __attribute__((always_inline)) inline bool operator()(
      const std::span<const V>& a,
      const std::span<const V>& b) const noexcept {
      return V_equal(a, b);
    }
  };

#ifdef F4NCGB_ENABLE_STORE_DUMP
  std::string dump_output_path_;
  void dump_to_binary() {
    msg("Dumping %d bytes to %s", size_, dump_output_path_.c_str());
    std::ofstream of(dump_output_path_, std::ios::binary | std::ios::out);
    for(size_t i = 0; i < size_; ++i) {
      of << static_cast<uint8_t>(pool_[i]);
    }
  }
#endif

  friend B;

  static_assert(
    alignof(M) >= alignof(V),
    "Alignment of M needs to be stricter, as elements are aligned this way.");

  protected:
  consteval static size_t capacity() {
#if defined(__has_feature)
#if __has_feature(address_sanitizer)
    size_t absolute_max = 0x10000000000;
#else
    size_t absolute_max = std::numeric_limits<size_t>::max();
#endif
#else
    size_t absolute_max = std::numeric_limits<size_t>::max();
#endif
    return std::min(
      static_cast<size_t>(std::numeric_limits<I>::max() * alignof(M)),
      absolute_max);
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
  I real_size_ = 0;
  I inserted_count_ = 0;
  M zero_metadata_;
  bool scratch_metadata_created_ = false;

  std::unique_ptr<std::byte[], pool_deleter> pool_;

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

  inline I scratch_id() const {
    I id = size_;

    return id + 1;
  }

  inline I insert_scratch() {
    I id = size_;

    void* ptr = boost::alignment::align_up(pool_.get() + size_, alignof(M));

    M* metadata = reinterpret_cast<M*>(std::assume_aligned<alignof(M)>(ptr));
    if(metadata->length == 0) {
      throw scratch_insertion_with_zero_length_exception(size_);
    }

    std::byte* ptr_start = static_cast<std::byte*>(ptr);
    ptr = boost::alignment::align_up(static_cast<std::byte*>(ptr) + sizeof(M),
                                     alignof(V));

    std::byte* ptr_end
      = static_cast<std::byte*>(ptr) + metadata->length * sizeof(V);
    std::byte* ptr_next = static_cast<std::byte*>(
      boost::alignment::align_up(ptr_end, alignof(M)));

    reinterpret_cast<B*>(this)->new_entry(id);

    size_ += ptr_next - ptr_start;
    ++inserted_count_;
    real_size_ = size_;
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
  using index_type = I;
  using value_type = V;

  size_t pool_capacity_ = 2097152 * capacity() / 2097152;

  store()
    : pool_(reinterpret_cast<std::byte*>(
        boost::alignment::aligned_alloc(2097152 /* 2^21, 2MB */,
                                        2097152 * (capacity() / 2097152)))) {

    size_t divisor = 1;
    while(!pool_.get()) {
      pool_.reset(reinterpret_cast<std::byte*>(
        std::malloc((2097152 / divisor) * (capacity() / 2097152))));
      pool_capacity_ = (2097152 / divisor) * (capacity() / 2097152);

      divisor *= 2;

      if(divisor > 32) {
        throw std::bad_alloc();
      }
    }
#ifdef __linux__
    madvise(pool_.get(), pool_capacity_, MADV_HUGEPAGE);
#endif
  }
  ~store() {
#ifdef F4NCGB_ENABLE_STORE_DUMP
    if(dump_output_path_ != "") {
      dump_to_binary();
    }
#endif

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

  inline void reset() { size_ = real_size_; }

  inline const std::span<const V> operator[](I id) const noexcept {
    if(id == 0)
      return std::span<V>();

    size_t len = get_length(id);

    const V* start = get_value_from_id(id - 1);

#ifdef F4NCGB_ENABLE_STORE_TRACE
    reinterpret_cast<B*>(const_cast<self*>(this))->tracer_.access(id);
#endif

    return std::span<const V>(start, len);
  }

  I next_pos(I pos) const {
    if(pos == 0)
      return 1;
    else
      --pos;
    // Increment the position to the next valid index.
    void* ptr = boost::alignment::align_up(pool_.get() + pos, alignof(M));

    assert(ptr == pool_.get() + pos);

    std::byte* ptr_start = static_cast<std::byte*>(ptr);
    ptr = boost::alignment::align_up(static_cast<std::byte*>(ptr) + sizeof(M),
                                     alignof(V));

    std::byte* ptr_end
      = static_cast<std::byte*>(ptr) + get_length(pos + 1) * sizeof(V);

    std::byte* ptr_next = static_cast<std::byte*>(
      boost::alignment::align_up(ptr_end, alignof(M)));

    return 1 + pos + ptr_next - ptr_start;
  }

  std::ostream& print(std::ostream& o) {
    o << "size: " << size_ << std::endl;
    for(auto pos : (*this)) {
      o << pos << ": =";
      for(auto v : (*this)[pos]) {
        o << " " << static_cast<int>(v);
      }
      o << ";" << std::endl;
    }
    return o;
  }

  struct pos_iterator {
    using iterator_category = std::forward_iterator_tag;
    using difference_type = std::ptrdiff_t;
    using value_type = I;
    using pointer = value_type*;
    using reference = const value_type&;

    pos_iterator(const self& s, I pos = 0)
      : s(s)
      , pos(pos) {}

    reference operator*() const { return pos; }
    pointer operator->() { return &pos; }
    pos_iterator& operator++() {
      pos = s.next_pos(pos);
      return *this;
    }
    pos_iterator operator++(int) {
      pos_iterator tmp = *this;
      ++(*this);
      return tmp;
    }
    friend bool operator==(const pos_iterator& a, const pos_iterator& b) {
      return &a.s == &b.s && a.pos == b.pos;
    };
    friend bool operator!=(const pos_iterator& a, const pos_iterator& b) {
      return &a.s != &b.s || a.pos != b.pos;
    };

    pos_iterator operator+(I count) const {
      pos_iterator it(*this);
      for(I i = 0; i < count; ++i) {
        ++it;
      }
      return it;
    }

    private:
    const self& s;
    I pos;
  };

  pos_iterator begin() const { return pos_iterator(*this, 0); }
  pos_iterator end() const { return pos_iterator(*this, size_ + 1); }
  size_t size() const { return (size_t)std::distance(begin(), end()); }

#ifdef F4NCGB_ENABLE_STORE_DUMP
  void set_binary_dump_path(const std::string& p) { dump_output_path_ = p; }
#endif
};

//================================================================
template<internal::metadata_concept M = internal::metadata<uint8_t>,
         internal::value_concept V = uint8_t,
         typename I = uint32_t,
         uint16_t Nblocks = 0>
class monomial_store;

template<metadata_concept M, value_concept V, typename I, uint16_t Nblocks>
class monomial_store : public store<monomial_store<M, V, I>, M, V, I> {
  using self = monomial_store<M, V, I>;
  using base = store<self, M, V, I>;
  friend base;

  using ambiguity_ = ambiguity<I>;
  using amb_hash = ambiguity_hash<I>;

#ifdef F4NCGB_ENABLE_STORE_TRACE
  store_tracer tracer_ = store_tracer("monomial_store");
#endif

#ifdef F4NCGB_USE_COMPACT_MONOMIAL_MAP
  using lookup_map
    = boost::unordered_flat_set<I,
                                typename base::V_hash_struct,
                                typename base::V_equality_struct>;

  lookup_map map_ = lookup_map(1000,
                               typename base::V_hash_struct(*this),
                               typename base::V_equality_struct(*this));
#else
  boost::unordered_flat_map<std::span<const V>,
                            I,
                            boost::hash<std::span<const V>>,
                            typename base::V_span_equality_struct>
    map_;
#endif

#ifdef F4NCGB_USE_MONOMIAL_PRODUCTS_MAP
  boost::unordered_flat_map<std::tuple<I, I, I>, I> products_;
#endif

  protected:
#ifdef F4NCGB_USE_COMPACT_MONOMIAL_MAP
  inline void new_entry(I id) { map_.insert(id + 1); }

  std::optional<I> findscratch() const {
    F4NCGB_PROFILE(gstats.store_find_calls++);

    I id = base::scratch_id();

    const auto it = map_.find(id);
    if(it == map_.end())
      return std::nullopt;

    F4NCGB_PROFILE(gstats.store_find_hits++);
    return *it;
  }

  std::optional<I> find(const std::span<const V>& v) {
    if(v.size() == 0)
      return 0;

    auto [m, vv] = this->new_scratch(v.size());
    m.length = v.size();
    std::copy(v.begin(), v.end(), vv);

    return findscratch();
  }
#else
  inline void new_entry(I id) { map_.insert(std::pair((*this)[id + 1], id)); }

  std::optional<I> find(const std::span<const V>& v) const {
    F4NCGB_PROFILE(gstats.store_find_calls++);
    if(v.size() == 0)
      return 0;

    const auto it = map_.find(v);
    if(it == map_.end())
      return std::nullopt;

    F4NCGB_PROFILE(gstats.store_find_hits++);
    return it->second + 1;
  }
#endif

  public:
  inline const I getid(const std::span<const V>& v) {
    auto id = find(v);
    if(!id) {
#ifdef F4NCGB_USE_COMPACT_MONOMIAL_MAP
      id = base::insert_scratch();
#else
      id = base::insert(v);
#endif
    }
    return *id;
  }
  inline const I getid(std::vector<V> v) {
    std::span<const V> s(v.begin(), v.size());
    return getid(s);
  }

  inline I get_product_id(I a, I b) {
    std::tuple<I, I, I> prod_tuple{ a, b, 0 };

#ifdef F4NCGB_USE_MONOMIAL_PRODUCTS_MAP
    {
      auto it = products_.find(prod_tuple);
      if(it != products_.end()) {
        return it->second;
      }
    }
#endif

    const size_t length_combined
      = static_cast<size_t>(base::get_length(a)) + base::get_length(b);
    if(length_combined
       > std::numeric_limits<typename base::length_type>::max()) {
      throw length_overrun_exception(
        length_combined,
        std::numeric_limits<typename base::length_type>::max());
    }
    auto [m, vv] = base::new_scratch(length_combined);
    m.length = length_combined;
    auto a_it = (*this)[a];
    auto b_it = (*this)[b];
    auto it = std::copy(a_it.begin(), a_it.end(), vv);
    std::copy(b_it.begin(), b_it.end(), it);
#ifdef F4NCGB_USE_COMPACT_MONOMIAL_MAP
    auto prod_idx = findscratch();
#else
    auto prod_idx = find(std::span(vv, length_combined));
#endif
    if(!prod_idx) {
      prod_idx = base::insert_scratch();
    }
#ifdef F4NCGB_USE_MONOMIAL_PRODUCTS_MAP
    products_.insert(std::make_pair(prod_tuple, *prod_idx));
#endif
    return *prod_idx;
  }

  inline I get_product_id(I a, I b, I c) {
    std::tuple<I, I, I> prod_tuple{ a, b, c };

#ifdef F4NCGB_USE_MONOMIAL_PRODUCTS_MAP
    F4NCGB_PROFILE(gstats.hashmap_calls++);

    {
      auto it = products_.find(prod_tuple);
      if(it != products_.end()) {
        F4NCGB_PROFILE(gstats.hashmap_hits++);
        return it->second;
      }
    }
#endif

    const size_t length_combined = static_cast<size_t>(base::get_length(a))
                                   + base::get_length(b) + base::get_length(c);
    if(length_combined
       > std::numeric_limits<typename base::length_type>::max()) {
      throw length_overrun_exception(
        length_combined,
        std::numeric_limits<typename base::length_type>::max());
    }
    auto [m, vv] = base::new_scratch(length_combined);
    m.length = length_combined;
    auto a_it = (*this)[a];
    auto b_it = (*this)[b];
    auto c_it = (*this)[c];
    auto it = std::copy(a_it.begin(), a_it.end(), vv);
    it = std::copy(b_it.begin(), b_it.end(), it);
    std::copy(c_it.begin(), c_it.end(), it);

#ifdef F4NCGB_USE_COMPACT_MONOMIAL_MAP
    auto prod_idx = findscratch();
#else
    auto prod_idx = find(std::span(vv, length_combined));
#endif

    if(!prod_idx) {
      prod_idx = base::insert_scratch();
    }
#ifdef F4NCGB_USE_MONOMIAL_PRODUCTS_MAP
    products_.insert(std::make_pair(prod_tuple, *prod_idx));
#endif
    return *prod_idx;
  }

  inline const std::span<V> get_product(I a, I b) {
    return (*this)[get_product_id(a, b)];
  }

  std::ostream& print_monomial(I i, std::ostream& o) const {
    auto m = (*this)[i];
    o << "(";
    for(const auto& v : m)
      o << static_cast<int>(v) << ", ";
    o << ")";
    return o;
  }

  inline bool is_divisible(const std::span<const V>& a,
                           const std::span<const V>& b) {
    return std::ranges::search(a, b).begin() != a.end();
  }
  inline bool is_divisible(I a, I b) {
    return is_divisible((*this)[a], (*this)[b]);
  }
  //-----------------------------------------------------------------
  // if v is  blocks[i] <= v <= blocks[i+1], then v is in block i
  std::vector<size_t> blocks;
  inline void set_blocks(std::vector<size_t> block_sizes) {
    blocks.resize(block_sizes.size() + 1);
    size_t n = 0;
    for(size_t i = 0; i < block_sizes.size(); i++) {
      blocks[i] = n;
      n += block_sizes[i];
    }
    blocks[blocks.size() - 1] = n;
  }

  template<bool block_order>
  inline bool cmp(I a, I b) {
    // this is a strict order
    if(a == b)
      return false;

    size_t la, lb;

    if constexpr(block_order) {
      auto a_it = (*this)[a];
      auto b_it = (*this)[b];

      // compare all blocks
      for(size_t i = 0; i < Nblocks; i++) {
        size_t n = Nblocks - i;
        la = (size_t)std::ranges::count_if(
          a_it,
          [l = blocks[n - 1], u = blocks[n]](V v) { return l < v && v <= u; });
        lb = (size_t)std::ranges::count_if(
          b_it,
          [l = blocks[n - 1], u = blocks[n]](V v) { return l < v && v <= u; });
        if(la != lb)
          return la < lb;
      }
    } else {
      // compare lengths
      la = base::get_length(a);
      lb = base::get_length(b);
      if(la != lb)
        return la < lb;
    }

    // compare monomials lexicographically
    auto a_it = (*this)[a];
    auto b_it = (*this)[b];

    int cmp = std::memcmp(a_it.data(), b_it.data(), a_it.size() * sizeof(V));
    return cmp < 0;
  }
  //-----------------------------------------------------------------
  std::ostream& print_ambiguity(const ambiguity_& a, std::ostream& o) const {
    o << "(" << a.degree() << ", ";
    print_monomial(a.ai(), o);
    o << ", ";
    print_monomial(a.ci(), o);
    o << ", ";
    print_monomial(a.aj(), o);
    o << ", ";
    print_monomial(a.cj(), o);
    o << ", " << static_cast<int>(a.i()) << ", " << static_cast<int>(a.j())
      << ")\n";
    return o;
  }
};

//================================================================

template<metadata_concept M, typename I>
struct polynomial_metadata : public M {
  I coefficients;
  uint16_t idx;
};

template<metadata_concept PM = internal::metadata<uint8_t>,
         metadata_concept MM = internal::metadata<uint8_t>,
         value_concept V = uint8_t,
         typename I = uint32_t,
         typename C = coeff,
         uint16_t Nblocks = 0>
class polynomial_store
  : public store<polynomial_store<PM, MM, V, I, C, Nblocks>,
                 polynomial_metadata<PM, I>,
                 I,
                 I> {

#ifdef F4NCGB_ENABLE_STORE_TRACE
  store_tracer tracer_ = store_tracer("polynomial_store");
#endif

  public:
  using self = polynomial_store<PM, MM, V, I, C, Nblocks>;
  using base = store<polynomial_store<PM, MM, V, I, C, Nblocks>,
                     polynomial_metadata<PM, I>,
                     I,
                     I>;
  using monomial_store_ = monomial_store<MM, V, I, Nblocks>;
  using monomial = monomial_store_::value_type;
  using metadata = polynomial_metadata<PM, I>;

  using polynomial_vec = std::vector<std::pair<C, I>>;
  using polynomial_nested_vec = std::vector<std::pair<C, std::vector<V>>>;

  friend base;

  constexpr static bool block_order = Nblocks > 1;

  constexpr static size_t capacity() {
#if defined(__has_feature)
#if __has_feature(address_sanitizer)
    size_t absolute_max = 0x10000000000;
#else
    size_t absolute_max = std::numeric_limits<size_t>::max();
#endif
#else
    size_t absolute_max = std::numeric_limits<size_t>::max();
#endif
    return std::min(
      static_cast<size_t>(static_cast<size_t>(std::numeric_limits<I>::max())
                          * sizeof(C)),
      absolute_max);
  }

  size_t cpool_capacity_ = 2097152 * capacity() / 2097152;

  inline polynomial_store(monomial_store_& store)
    : base::store()
    , store_(store)
    , cpool_(reinterpret_cast<std::byte*>(
        boost::alignment::aligned_alloc(2097152 /* 2^21, 2MB */,
                                        2097152 * (capacity() / 2097152)))) {
    cpool_capacity_ = 2097152 * (capacity() / 2097152);
    size_t divisor = 1;

    while(!cpool_.get()) {
      cpool_.reset(reinterpret_cast<std::byte*>(
        std::malloc((2097152 / divisor) * (capacity() / 2097152))));
      cpool_capacity_ = (2097152 / divisor) * (capacity() / 2097152);

      divisor *= 2;

      if(divisor > 32) {
        throw std::bad_alloc();
      }
    }
#ifdef __linux__
    madvise(cpool_.get(), cpool_capacity_, MADV_HUGEPAGE);
#endif
  }

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

  inline I add_polynomial(const polynomial_vec& p) {
    if(p.size() > std::numeric_limits<typename base::length_type>::max()) {
      throw length_overrun_exception(
        p.size(), std::numeric_limits<typename base::length_type>::max());
    }
    auto [m, monomials, coefficients] = add(p.size());
    for(size_t i = 0; i < p.size(); ++i) {
      C* c = new(coefficients + i) C;
      *c = p[i].first;
      monomials[i] = p[i].second;
    }
    m.length = p.size();
    return commit();
  }
  inline I add_polynomial(const polynomial_nested_vec& p) {
    auto [m, monomials, coefficients] = add(p.size());
    for(size_t i = 0; i < p.size(); ++i) {
      C* c = new(coefficients + i) C;
      *c = p[i].first;
      monomials[i] = store_.getid(p[i].second);
    }
    m.length = p.size();
    return commit();
  }

  inline std::span<C> get_coefficients(I id) {
    if(id == 0)
      return std::span<C>();
    metadata& m = this->get_metadata(id);
    return std::span<C>(get_coefficients_raw(m.coefficients), m.length);
  }

  inline std::span<const C> get_coefficients(I id) const {
    return const_cast<self*>(this)->get_coefficients(id);
  }

  // not needed
  // inline std::vector<I> get_monomial_ids(I id) {
  //   if(id == 0)
  //     return std::vector<I>(0);
  //   std::vector<I> res;
  //   res.reserve(this->get_length(id));
  //   for(size_t i = 0; i < this->get_length(id); i++)
  //     res.push_back((*this)[id][i]);
  //   return res;
  // }

  inline I get_monomial_id(I id, I idx) const {
    if(id == 0)
      return 0;
    assert(idx < this->get_metadata(id).length);
    return (*this)[id][idx];
  }

  inline void sort_polynomial(std::vector<std::pair<C, I>>& p) {
    std::stable_sort(p.begin(), p.end(), [this](const auto& a, const auto& b) {
      return this->store_.template cmp<block_order>(b.second, a.second);
    });
  }

  inline I get_lm_id(I id) const {
    // discuss if polynomials are sorted
    // increasing or decreasing -- decreasing
    // is probably better
    return get_monomial_id(id, 0);
  }

  inline void set_idx(I id, uint16_t idx) {
    assert(idx != 0);
    if(id == 0)
      return;
    metadata& m = this->get_metadata(id);
    m.idx = idx;
  }

  inline uint16_t get_idx(I id) {
    if(id == 0)
      return 0;
    metadata& m = this->get_metadata(id);
    return m.idx;
  }

  template<bool front, bool back>
  inline I multiply_front_or_back_or_both(I f, I p, I b) {
    assert(p != 0);

    size_t cpool_size_old = cpool_size_;
    size_t size_old = this->real_size_;
    metadata& p_metadata = this->get_metadata(p);
    auto [new_m, new_i, new_c] = add(p_metadata.length);
    new_m.length = p_metadata.length;
    for(I i = 0; i < p_metadata.length; ++i) {
      if constexpr(front && !back) {
        new_i[i] = store_.get_product_id(f, (*this)[p][i]);
      } else if constexpr(!front && back) {
        new_i[i] = store_.get_product_id((*this)[p][i], b);
      } else if constexpr(front && back) {
        new_i[i] = store_.get_product_id(f, (*this)[p][i], b);
      }
    }
    I g = commit();
    this->get_metadata(g).coefficients = p_metadata.coefficients;
    cpool_size_ = cpool_size_old;
    this->real_size_ = size_old;
    return g;
  }

  inline I multiply_front(I m, I p) {
    return multiply_front_or_back_or_both<true, false>(m, p, 0);
  }
  inline I multiply_back(I p, I m) {
    return multiply_front_or_back_or_both<false, true>(0, p, m);
  }
  inline I multiply_front_and_back(I f, I p, I b) {
    return multiply_front_or_back_or_both<true, true>(f, p, b);
  }

  const monomial_store_& get_monomial_store() const { return store_; }
  monomial_store_& get_monomial_store() { return store_; }

  std::ostream& print_polynomial(I i, std::ostream& o) const {
    auto c = this->get_coefficients(i);
    size_t j = 0;
    for(auto m : (*this)[i]) {
      o << c[j++];
      o << "*";
      store_.print_monomial(m, o);
    }
    return o;
  }

  protected:
  monomial_store_& store_;
  std::unique_ptr<std::byte[], typename base::pool_deleter> cpool_;
  size_t cpool_size_ = 0;

  inline C* get_coefficients_raw(I id) {
    if(static_cast<size_t>(id) * sizeof(C) + sizeof(C) > cpool_capacity_) {
      throw coefficient_overrun_exception(id * sizeof(C), cpool_capacity_);
    }
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

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif

class parser_context;
int
f4ncgb_main(parser_context& context,
            const std::string& output_name,
            std::function<void()>* stats_print_function = nullptr,
            bool leak_memory = false,
            bool print_read_problem = false,
            void* userdata = nullptr,
            f4ncgb_add_cb add_cb = nullptr,
            f4ncgb_end_poly_cb end_poly_cb = nullptr);
}
