v0.3
====

## New Features

  - Introduce an API (C-based) in `f4ncgb.hpp` that can start basis
    calculations and return the computed basis.
  - Add a handler for SIGUSR1 that prints the basis that has already
    been computed up to this point.
  - Add support for shared library builds.
  - Add support for dumping the stores at destruction (enable with
    CMake-option `-DENABLE_STORE_DUMP=ON`).
    
## Changes

  - Both `proof_level` and `tracer` are no longer global symbols,
    instead they are given as parameters to a solve call.
