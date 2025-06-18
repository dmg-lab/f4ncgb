#include "f4ncgb.hpp"
#include <climits>
#include <config.hpp>
#include <vector>

using namespace f4ncgb;

#define REQUIRE_HANDLE(h)     \
  if(!h) {                    \
    return "no handle given"; \
  }

typedef struct f4ncgb_handle {
  f4ncgb_state state = F4NCGB_STATE_INITIAL;
  uint32_t nblocks = 1;
  uint32_t nvars = 0;
  uint32_t maxiter = 10;
  uint32_t maxdeg = UINT_MAX;
  uint32_t threads = 1;
  uint32_t characteristic = 0;
  const char* output_file = nullptr;
  const char* proof_file = nullptr;

  std::vector<std::tuple<long, long, size_t>> polynomial;
} f4ncgb_handle;

extern "C" const char*
f4ncgb_version() {
  return F4NCGB_VERSION;
}

extern "C" f4ncgb_handle*
f4ncgb_init() {
  f4ncgb_handle* h = new f4ncgb_handle;
  return h;
}

extern "C" void
f4ncgb_free(f4ncgb_handle* h) {
  delete h;
}

extern "C" f4ncgb_state
f4ncgb_get_state(const f4ncgb_handle* h) {
  return h->state;
}

extern "C" const char*
f4ncgb_prepare(f4ncgb_handle* h) {
  REQUIRE_HANDLE(h);
  if(h->state != F4NCGB_STATE_INITIAL) {
    return "invalid state, must be in INITIAL";
  }

  if(h->nvars == 0) {
    return "nvars not set";
  }

  // TODO: Create f4ncgb algo. This needs some unbundling from the context in
  // f4.hpp.
}

extern "C" const char*
f4ncgb_add_monomial(f4ncgb_handle* h,
                    long numerator,
                    long denominator,
                    size_t varcount,
                    uint32_t* vars) {
  REQUIRE_HANDLE(h);
  if(h->state != F4NCGB_STATE_ADD) {
    return "invalid state, must be in ADD";
  }

  if(varcount) {
    // TODO: Add vector to f4ncgb algo.
  } else {
    // TODO: Add vector to f4ncgb algo.
  }
  return nullptr;
}

extern "C" const char*
f4ncgb_set_nblocks(f4ncgb_handle* h, uint32_t nblocks) {
  REQUIRE_HANDLE(h);
  if(h->state != F4NCGB_STATE_INITIAL) {
    return "invalid state, must be in INITIAL";
  }
  h->nblocks = nblocks;
  return nullptr;
}

extern "C" const char*
f4ncgb_set_nvars(f4ncgb_handle* h, uint32_t nvars) {
  REQUIRE_HANDLE(h);
  if(h->state != F4NCGB_STATE_INITIAL) {
    return "invalid state, must be in INITIAL";
  }
  h->nvars = nvars;
  return nullptr;
}

extern "C" const char*
f4ncgb_set_characteristic(f4ncgb_handle* h, uint32_t characteristic) {
  REQUIRE_HANDLE(h);
  if(h->state != F4NCGB_STATE_INITIAL) {
    return "invalid state, must be in INITIAL";
  }
  h->characteristic = characteristic;
  return nullptr;
}

extern "C" const char*
f4ncgb_set_maxiter(f4ncgb_handle* h, uint32_t maxiter) {
  REQUIRE_HANDLE(h);
  if(h->state != F4NCGB_STATE_INITIAL) {
    return "invalid state, must be in INITIAL";
  }
  h->maxiter = maxiter;
  return nullptr;
}

extern "C" const char*
f4ncgb_set_maxdeg(f4ncgb_handle* h, uint32_t maxdeg) {
  REQUIRE_HANDLE(h);
  if(h->state != F4NCGB_STATE_INITIAL) {
    return "invalid state, must be in INITIAL";
  }
  h->maxdeg = maxdeg;
  return nullptr;
}

extern "C" const char*
f4ncgb_set_threads(f4ncgb_handle* h, uint32_t threads) {
  REQUIRE_HANDLE(h);
  if(h->state != F4NCGB_STATE_INITIAL) {
    return "invalid state, must be in INITIAL";
  }
  h->threads = threads;
  return nullptr;
}

extern "C" const char*
f4ncgb_set_output_file(f4ncgb_handle* h, const char* output_file) {
  REQUIRE_HANDLE(h);
  if(h->state != F4NCGB_STATE_INITIAL) {
    return "invalid state, must be in INITIAL";
  }
  h->output_file = output_file;
  return nullptr;
}

extern "C" const char*
f4ncgb_set_proof_file(f4ncgb_handle* h, const char* proof_file) {
  REQUIRE_HANDLE(h);
  if(h->state != F4NCGB_STATE_INITIAL) {
    return "invalid state, must be in INITIAL";
  }
  h->proof_file = proof_file;
  return nullptr;
}

extern "C" int
f4ncgb_solve(f4ncgb_handle* h) {
  // TODO: Issue solve call
  return 0;
}
