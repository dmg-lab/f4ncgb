#include <csignal>
#include <cstdio>
#include <cstdlib>

#include <sys/resource.h>
#include <sys/time.h>

#include "profiling.hpp"
#include "signal_statistics.hpp"

namespace kommunopp {
statistics gstats;

static size_t
maximum_resident_set_size() {
  struct rusage u;
  if(getrusage(RUSAGE_SELF, &u))
    return 0;
  return ((size_t)u.ru_maxrss) << 10;
}

/*------------------------------------------------------------------------*/

static double
process_time() {
  struct rusage u;
  if(getrusage(RUSAGE_SELF, &u))
    return 0;
  double res = u.ru_utime.tv_sec + 1e-6 * u.ru_utime.tv_usec;
  res += u.ru_stime.tv_sec + 1e-6 * u.ru_stime.tv_usec;
  return res;
}

static double
percent(unsigned a, unsigned b) {
  return b ? 100.0 * a / b : 0;
}

static double
percent_d(double a, double b) {
  return b ? 100.0 * a / b : 0;
}
/*------------------------------------------------------------------------*/

void
statistics::print() {
  auto total_time = process_time();
  statistics& s = gstats;
  msg("");
  msg("maximum resident set size:  %22.2f MB",
      maximum_resident_set_size() / static_cast<double>((1 << 20)));
  msg("monomial hashmap products:   %22d", gstats.hashmap_calls);
  msg("monomial hashmap prodhits: %22.2f %%",
      percent(s.hashmap_hits, s.hashmap_calls));
  msg("store find calls:            %22d", gstats.store_find_calls);
  msg("store find hits:           %22.2f %%",
      percent(s.store_find_hits, s.store_find_calls));
  msg("computing ambiguities:     %22.2f (%2.2f %%)",
      s.amb,
      percent_d(s.amb, total_time));
  msg("  computing overlaps:      %22.2f (%2.2f %%)",
      s.overlap,
      percent_d(s.overlap, total_time));
  msg("  computing inclusions:    %22.2f (%2.2f %%)",
      s.inclusion,
      percent_d(s.inclusion, total_time));
  msg("handling critical pairs:   %22.2f (%2.2f %%)",
      s.crit_pair,
      percent_d(s.crit_pair, total_time));
  msg("symbolic preprocessing:    %22.2f (%2.2f %%)",
      s.sym_pre,
      percent_d(s.sym_pre, total_time));
  msg("linear algebra:            %22.2f (%2.2f %%)",
      s.reduction,
      percent_d(s.reduction, total_time));
  msg("  Gauss elimination:       %22.2f (%2.2f %%)",
      s.rref,
      percent_d(s.rref, total_time));
  msg("  CRT:                     %22.2f (%2.2f %%)",
      s.crt,
      percent_d(s.crt, total_time));
  msg("  rat. reconstruction:     %22.2f (%2.2f %%)",
      s.ratrec,
      percent_d(s.ratrec, total_time));
  msg("construct new elements:    %22.2f (%2.2f %%)",
      s.new_elements,
      percent_d(s.new_elements, total_time));
  msg("other:                     %22.2f (%2.2f %%)",
      s.other,
      percent_d(s.other, total_time));
  msg("binary search:             %22.2f (%2.2f %%)",
      s.binarysearch,
      percent_d(s.binarysearch, total_time));
  msg("");
  msg("total process time:        %22.2f seconds", process_time());
}
}
