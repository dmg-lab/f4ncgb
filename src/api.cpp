#include "f4ncgb.hpp"
#include "parser.hpp"
#include "signal_statistics.hpp"
#include <climits>
#include <config.hpp>
#include <numeric>
#include <vector>

#include <ulong_extras.h>

using namespace f4ncgb;

#define REQUIRE_HANDLE(h)     \
  if(!h) {                    \
    return "no handle given"; \
  }

using monomial = std::vector<uint32_t>;
using polynomial = std::vector<std::tuple<long, long, monomial>>;

typedef struct f4ncgb_handle {
  f4ncgb_state state = F4NCGB_STATE_READY;
  uint32_t maxiter = 10;
  uint32_t maxdeg = UINT_MAX;
  uint32_t threads = 1;
  uint32_t proof_level = 0;
  bool tracer = true;
  const char* output_file = "";
  const char* proof_file = "";

  parser_context ctx;

  polynomial current_polynomial;
  std::vector<polynomial> polynomials;
} f4ncgb_handle;

extern "C" void
f4ncgb_set_msg_printing(bool printing) {
  msg_printing = printing ? 1 : 0;
}

extern "C" const char*
f4ncgb_version() {
  return F4NCGB_VERSION;
}

extern "C" f4ncgb_handle*
f4ncgb_init() {
  f4ncgb_handle* h = new f4ncgb_handle;

  h->ctx.impl = [h](parse_add_cb add_cb,
                    void* add_cb_userdata,
                    parse_monomial_boundary_cb boundary_cb,
                    void* boundary_cb_userdata) {
    for(const auto& poly : h->polynomials) {
      for(const auto& [numerator, denominator, mono] : poly) {
        for(uint32_t v : mono) {
          add_cb(add_cb_userdata, v);
        }
        boundary_cb(boundary_cb_userdata,
                    numerator,
                    denominator,
                    numerator == denominator);
      }
      add_cb(add_cb_userdata, 0);
    }
    return std::nullopt;
  };

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
f4ncgb_add(f4ncgb_handle* h,
           long numerator,
           long denominator,
           size_t varcount,
           uint32_t* vars) {
  REQUIRE_HANDLE(h);
  if(h->state != F4NCGB_STATE_READY) {
    return "invalid state, must be in READY";
  }

  if(varcount) {
    h->current_polynomial.emplace_back(
      numerator, denominator, monomial(vars, vars + varcount));
  } else {
    h->current_polynomial.emplace_back(numerator, denominator, monomial());
  }
  return nullptr;
}

extern "C" const char*
f4ncgb_end_poly(f4ncgb_handle* h) {
  REQUIRE_HANDLE(h);
  if(h->state != F4NCGB_STATE_READY) {
    return "invalid state, must be in READY";
  }
  h->polynomials.emplace_back(h->current_polynomial);
  h->current_polynomial.clear();
  return nullptr;
}

extern "C" const char*
f4ncgb_set_blocks(f4ncgb_handle* h,
                  uint32_t blockcount,
                  uint32_t* blocklengths) {
  REQUIRE_HANDLE(h);
  if(h->state != F4NCGB_STATE_READY) {
    return "invalid state, must be in READY";
  }

  if(blockcount > F4NCGB_MAX_BLOCKS)
    return "more blocks than current compilation allows";

  h->ctx.blocks_.resize(blockcount);
  uint32_t v = 1;
  for(size_t i = 0; i < blockcount; ++i) {
    h->ctx.blocks_[i].resize(blocklengths[i]);
    auto begin = h->ctx.blocks_[i].begin();
    auto end = h->ctx.blocks_[i].end();
    std::iota(begin, end, v);
    v += blocklengths[i];
  }

  if(h->ctx.num_vars() == 0) {
    h->ctx.num_vars_ = v;
  }

  return nullptr;
}

extern "C" const char*
f4ncgb_set_nvars(f4ncgb_handle* h, uint32_t nvars) {
  REQUIRE_HANDLE(h);
  if(h->state != F4NCGB_STATE_READY) {
    return "invalid state, must be in READY";
  }
  h->ctx.num_vars_ = nvars;
  return nullptr;
}

extern "C" const char*
f4ncgb_set_characteristic(f4ncgb_handle* h, uint32_t characteristic) {
  REQUIRE_HANDLE(h);
  if(h->state != F4NCGB_STATE_READY) {
    return "invalid state, must be in READY";
  }

  if(characteristic > 2147483647l) {// 2^31 -1
    return "provided characteristic %lu is too large. Only p < 2^31 supported";
  } else if(characteristic != 0 and !n_is_prime(characteristic)) {
    return "provided nonzero characteristic is not prime";
  }

  h->ctx.characteristic_ = characteristic;
  return nullptr;
}

extern "C" const char*
f4ncgb_set_maxiter(f4ncgb_handle* h, uint32_t maxiter) {
  REQUIRE_HANDLE(h);
  if(h->state != F4NCGB_STATE_READY) {
    return "invalid state, must be in READY";
  }
  h->maxiter = maxiter;
  return nullptr;
}

extern "C" const char*
f4ncgb_set_maxdeg(f4ncgb_handle* h, uint32_t maxdeg) {
  REQUIRE_HANDLE(h);
  if(h->state != F4NCGB_STATE_READY) {
    return "invalid state, must be in READY";
  }
  h->maxdeg = maxdeg;
  return nullptr;
}

extern "C" const char*
f4ncgb_set_threads(f4ncgb_handle* h, uint32_t threads) {
  REQUIRE_HANDLE(h);
  if(h->state != F4NCGB_STATE_READY) {
    return "invalid state, must be in READY";
  }
  h->threads = threads;
  return nullptr;
}

extern "C" const char*
f4ncgb_set_output_file(f4ncgb_handle* h, const char* output_file) {
  REQUIRE_HANDLE(h);
  if(h->state != F4NCGB_STATE_READY) {
    return "invalid state, must be in READY";
  }
  h->output_file = output_file;
  return nullptr;
}

extern "C" const char*
f4ncgb_set_proof_file(f4ncgb_handle* h, const char* proof_file) {
  REQUIRE_HANDLE(h);
  if(h->state != F4NCGB_STATE_READY) {
    return "invalid state, must be in READY";
  }
  h->proof_file = proof_file;
  h->proof_level = 1;
  return nullptr;
}

extern "C" const char*
f4ncgb_set_expanded_proof(f4ncgb_handle* h, bool expanded) {
  REQUIRE_HANDLE(h);
  if(h->state != F4NCGB_STATE_READY) {
    return "invalid state, must be in READY";
  }
  if(expanded)
    h->proof_level = 2;
  return nullptr;
}

extern "C" const char*
f4ncgb_set_tracer(f4ncgb_handle* h, bool tracer) {
  REQUIRE_HANDLE(h);
  if(h->state != F4NCGB_STATE_READY) {
    return "invalid state, must be in READY";
  }
  h->tracer = tracer;
  return nullptr;
}

extern "C" f4ncgb_result
f4ncgb_solve(f4ncgb_handle* h,
             void* userdata,
             f4ncgb_add_cb add_cb,
             f4ncgb_end_poly_cb end_cb) {
  if(!h)
    return F4NCGB_ARGERROR;
  if(h->ctx.num_vars() > 0 && h->ctx.num_blocks() == 0) {
    h->ctx.blocks_.emplace_back();
    h->ctx.blocks_[0].resize(h->ctx.num_vars() + 1);
    auto begin = h->ctx.blocks_[0].begin();
    auto end = h->ctx.blocks_[0].end();
    std::iota(begin, end, 1);
  }
  h->ctx.num_blocks_ = h->ctx.blocks_.size();

  if(h->ctx.num_blocks() > F4NCGB_MAX_BLOCKS)
    return F4NCGB_ARGERROR;
  if(h->ctx.characteristic() > 2147483647l)// 2^31 -1
    return F4NCGB_ARGERROR;
  if(h->ctx.characteristic() != 0 and !n_is_prime(h->ctx.characteristic()))
    return F4NCGB_ARGERROR;
  if(!h->proof_file and h->proof_level > 0)
    return F4NCGB_ARGERROR;

  h->state = F4NCGB_STATE_SOLVING;

  int res = f4ncgb::f4ncgb_main(h->ctx,
                                h->ctx.num_blocks(),
                                h->ctx.num_vars(),
                                h->ctx.characteristic(),
                                h->maxiter,
                                h->maxdeg,
                                h->threads,
                                h->proof_level,
                                h->tracer,
                                h->output_file,
                                h->proof_file,
                                nullptr,
                                false, /* No leaking */
                                false, /* No problem printing */
                                userdata,
                                add_cb,
                                end_cb);

  h->state = F4NCGB_STATE_READY;

  switch(res) {
    case 0:
      return F4NCGB_OK;
    default:
      return F4NCGB_ERROR;
  }
}
