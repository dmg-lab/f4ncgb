#ifndef F4NCGB_HPP
#define F4NCGB_HPP

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <gmp.h>

#ifdef __cplusplus
extern "C" {
#endif

const char*
f4ncgb_version();

typedef struct f4ncgb_handle f4ncgb_handle;

typedef enum f4ncgb_state {
  F4NCGB_STATE_READY,
  F4NCGB_STATE_SOLVING,
} f4ncgb_state;

typedef enum f4ncgb_result {
  F4NCGB_OK = 0,
  F4NCGB_ERROR,
  F4NCGB_ARGERROR,
  F4NCGB_UNKNOWN,
} f4ncgb_result;

void
f4ncgb_set_msg_printing(bool);

f4ncgb_handle*
f4ncgb_init();

void
f4ncgb_free(f4ncgb_handle*);

f4ncgb_state
f4ncgb_get_state(const f4ncgb_handle*);

// Add a monomial with some coefficient.
//
// This builds up the current polynomial. Polynomials are
// automatically sorted after they are added like this.
const char*
f4ncgb_add(f4ncgb_handle*,
           long numerator,
           long denominator,
           size_t varcount,
           uint32_t* vars);

const char*
f4ncgb_end_poly(f4ncgb_handle*);

const char*
f4ncgb_set_blocks(f4ncgb_handle*, uint32_t blockcount, uint32_t* blocklengths);

const char*
f4ncgb_set_nvars(f4ncgb_handle*, uint32_t);

const char*
f4ncgb_set_characteristic(f4ncgb_handle*, uint32_t);

const char*
f4ncgb_set_maxiter(f4ncgb_handle*, uint32_t);

const char*
f4ncgb_set_maxdeg(f4ncgb_handle*, uint32_t);

const char*
f4ncgb_set_threads(f4ncgb_handle*, uint32_t);

const char*
f4ncgb_set_tracer(f4ncgb_handle*, bool);

const char*
f4ncgb_set_expanded_proof(f4ncgb_handle*, bool);

const char*
f4ncgb_set_output_file(f4ncgb_handle*, const char*);

const char*
f4ncgb_set_proof_file(f4ncgb_handle*, const char*);

typedef void (*f4ncgb_add_cb)(void* userdata,
                              mpz_ptr numerator,
                              mpz_ptr denominator,
                              size_t varcount,
                              const uint32_t* vars);

typedef void (*f4ncgb_end_poly_cb)(void* userdata);

f4ncgb_result
f4ncgb_solve(f4ncgb_handle*,
             void* userdata,
             f4ncgb_add_cb add_cb,
             f4ncgb_end_poly_cb end_cb);

f4ncgb_result
f4ncgb_reduce(f4ncgb_handle*,
             void* userdata,
             f4ncgb_add_cb add_cb,
             f4ncgb_end_poly_cb end_cb);

#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
#include <functional>
#include <memory>
#include <span>
#include <vector>

namespace f4ncgb {
class Solver {
  struct Deleter {
    void operator()(f4ncgb_handle* h) { f4ncgb_free(h); }
  };
  std::unique_ptr<f4ncgb_handle, Deleter> handle_;

  public:
  Solver() { handle_.reset(f4ncgb_init()); }
  ~Solver() = default;

  f4ncgb_state state() const { return f4ncgb_get_state(handle_.get()); }

  void add(long numerator, long denominator, std::vector<uint32_t> vars) {
    const char* msg = f4ncgb_add(
      handle_.get(), numerator, denominator, vars.size(), vars.data());
    if(msg) {
      throw std::runtime_error(msg);
    }
  }

  void end_poly() {
    const char* msg = f4ncgb_end_poly(handle_.get());
    if(msg) {
      throw std::runtime_error(msg);
    }
  }

  void set_blocks(std::vector<uint32_t> b) {
    const char* msg = f4ncgb_set_blocks(handle_.get(), b.size(), b.data());
    if(msg) {
      throw std::runtime_error(msg);
    }
  }
  void set_nvars(uint32_t nvars) {
    const char* msg = f4ncgb_set_nvars(handle_.get(), nvars);
    if(msg) {
      throw std::runtime_error(msg);
    }
  }
  void set_characteristic(uint32_t characteristic) {
    const char* msg = f4ncgb_set_characteristic(handle_.get(), characteristic);
    if(msg) {
      throw std::runtime_error(msg);
    }
  }
  void set_maxiter(uint32_t maxiter) {
    const char* msg = f4ncgb_set_maxiter(handle_.get(), maxiter);
    if(msg) {
      throw std::runtime_error(msg);
    }
  }
  void set_maxdeg(uint32_t maxdeg) {
    const char* msg = f4ncgb_set_maxdeg(handle_.get(), maxdeg);
    if(msg) {
      throw std::runtime_error(msg);
    }
  }
  void set_threads(uint32_t threads) {
    const char* msg = f4ncgb_set_threads(handle_.get(), threads);
    if(msg) {
      throw std::runtime_error(msg);
    }
  }
  void set_tracer(bool tracer) {
    const char* msg = f4ncgb_set_tracer(handle_.get(), tracer);
    if(msg) {
      throw std::runtime_error(msg);
    }
  }
  void set_expanded_proof(bool expanded) {
    const char* msg = f4ncgb_set_expanded_proof(handle_.get(), expanded);
    if(msg) {
      throw std::runtime_error(msg);
    }
  }
  void set_output_file(const std::string& output_file) {
    const char* msg
      = f4ncgb_set_output_file(handle_.get(), output_file.c_str());
    if(msg) {
      throw std::runtime_error(msg);
    }
  }
  void set_proof_file(const std::string& proof_file) {
    const char* msg = f4ncgb_set_proof_file(handle_.get(), proof_file.c_str());
    if(msg) {
      throw std::runtime_error(msg);
    }
  }

  using add_cb
    = std::function<void(mpz_ptr, mpz_ptr, std::span<const uint32_t>)>;
  using end_poly_cb = std::function<void()>;

  f4ncgb_result solve(add_cb add, end_poly_cb end) {
    struct meta {
      add_cb a;
      end_poly_cb e;
    };
    meta m{ add, end };

    auto c_add_cb = [](void* userdata,
                       mpz_ptr numerator,
                       mpz_ptr denominator,
                       size_t varcount,
                       const uint32_t* vars) {
      meta& m = *static_cast<meta*>(userdata);
      return m.a(
        numerator, denominator, std::span<const uint32_t>(vars, varcount));
    };
    auto c_end_poly_cb = [](void* userdata) {
      meta& m = *static_cast<meta*>(userdata);
      return m.e();
    };

    return f4ncgb_solve(
      handle_.get(), static_cast<void*>(&m), c_add_cb, c_end_poly_cb);
  }

  f4ncgb_result reduce(add_cb add, end_poly_cb end) {
    struct meta {
      add_cb a;
      end_poly_cb e;
    };
    meta m{ add, end };

    auto c_add_cb = [](void* userdata,
                       mpz_ptr numerator,
                       mpz_ptr denominator,
                       size_t varcount,
                       const uint32_t* vars) {
      meta& m = *static_cast<meta*>(userdata);
      return m.a(
        numerator, denominator, std::span<const uint32_t>(vars, varcount));
    };
    auto c_end_poly_cb = [](void* userdata) {
      meta& m = *static_cast<meta*>(userdata);
      return m.e();
    };

    return f4ncgb_reduce(
      handle_.get(), static_cast<void*>(&m), c_add_cb, c_end_poly_cb);
  }

  using monomial = std::vector<uint32_t>;
  using polynomial = std::vector<std::tuple<mpz_ptr, mpz_ptr, monomial>>;

  std::pair<f4ncgb_result, std::vector<polynomial>> solve() {
    std::vector<polynomial> polys;
    bool create = true;
    f4ncgb_result res = solve(
      [&polys,
       &create](mpz_ptr num, mpz_ptr den, std::span<const uint32_t> vars) {
        if(create) {
          polys.emplace_back();
          create = false;
        }
        polys.back().emplace_back(
          num, den, std::vector<uint32_t>(vars.begin(), vars.end()));
      },
      [&create]() { create = true; });
    return std::make_pair(res, std::move(polys));
  }
};
}
#endif

#endif
