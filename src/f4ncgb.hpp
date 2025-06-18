#ifndef F4NCGB_HPP
#define F4NCGB_HPP

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

const char*
f4ncgb_version();

typedef struct f4ncgb_handle f4ncgb_handle;

typedef enum f4ncgb_state {
  F4NCGB_STATE_INITIAL,
  F4NCGB_STATE_ADD,
} f4ncgb_state;

f4ncgb_handle*
f4ncgb_init();

void
f4ncgb_free(f4ncgb_handle*);

f4ncgb_state
f4ncgb_get_state(const f4ncgb_handle*);

// Prepares the internal state from INITIAL to ADD. Uses all set parameters to
// decide on the implementation to use and prepares the solver for the addition
// of monomials and polynomials.
const char*
f4ncgb_prepare(f4ncgb_handle* h);

// Add a monomial with some coefficient.
//
// This builds up the current polynomial. In order to terminate a
// polynomial, issue a polynomial with varcount=0. Polynomials are
// automatically sorted after they are added like this.
const char*
f4ncgb_add(f4ncgb_handle*,
           long numerator,
           long denominator,
           size_t varcount,
           uint32_t* vars);

const char*
f4ncgb_set_nblocks(f4ncgb_handle*, uint32_t);

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
f4ncgb_set_output_file(f4ncgb_handle*, const char*);

const char*
f4ncgb_set_proof_file(f4ncgb_handle*, const char*);

int
f4ncgb_solve(f4ncgb_handle*);

#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
#include <memory>

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

  void set_nblocks(uint32_t nblocks) {
    const char* msg = f4ncgb_set_nblocks(handle_.get(), nblocks);
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
};
}
#endif

#endif
