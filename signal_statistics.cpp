/*------------------------------------------------------------------------*/
/*! \file signal_statistics.cpp
    \brief used to handle signals, messages and statistics

  Part of FGLMultiLing
  Copyright(C) 2024 Daniela Kaufmann,TU Wien, Austria
*/
/*------------------------------------------------------------------------*/

#include "signal_statistics.hpp"

/*------------------------------------------------------------------------*/
// Global variable
void(*original_SIGINT_handler)(int);
void(*original_SIGSEGV_handler)(int);
void(*original_SIGABRT_handler)(int);
void(*original_SIGTERM_handler)(int);
/*------------------------------------------------------------------------*/
const char * bench = 0;
int van_mon_depth_count = 0;
int child_superset_count= 0;
int degree_three_polys_count= 0;
int degree_three_of_child_count= 0;
int grandchild_is_subset_of_child_count= 0;
int grandchild_is_superset_child_count= 0;
int two_subsets_count= 0;
int grand_children_are_equal_count=0;
int l_f_count=0;
int child_l_f_count=0;
int children_share_l_f_count = 0;
int totalcount = 0;
int equiv_gate_count = 0;
int lin_gb_count = 0;
int count_cocoa_calls = 0;
bool miter_inp = 0;
bool mult_inp = 0;
bool msolve = 0;

const char * signal_name(int sig) {
  switch (sig) {
    case SIGINT: return "SIGINT";
    case SIGSEGV: return "SIGSEGV";
    case SIGABRT: return "SIGABRT";
    case SIGTERM: return "SIGTERM";
    default: return "SIGUNKNOWN";
  }
}

/*------------------------------------------------------------------------*/

void init_all_signal_handers() {
  original_SIGINT_handler  = signal(SIGINT,  catch_signal);
  original_SIGSEGV_handler = signal(SIGSEGV, catch_signal);
  original_SIGABRT_handler = signal(SIGABRT, catch_signal);
  original_SIGTERM_handler = signal(SIGTERM, catch_signal);
}

/*------------------------------------------------------------------------*/

void catch_signal(int sig) {
  printf("c\nc caught signal '%s'(%d)\nc\n", signal_name(sig), sig);
  printf("c\nc raising signal '%s'(%d) again\n", signal_name(sig), sig);
  reset_all_signal_handlers();
  fflush(stdout);
  raise(sig);
}

/*------------------------------------------------------------------------*/

void reset_all_signal_handlers() {
  (void) signal(SIGINT, original_SIGINT_handler);
  (void) signal(SIGSEGV, original_SIGSEGV_handler);
  (void) signal(SIGABRT, original_SIGABRT_handler);
  (void) signal(SIGTERM, original_SIGTERM_handler);
}

/*------------------------------------------------------------------------*/
// Global variable
int verbose = 1;

/*------------------------------------------------------------------------*/

void msg(const char *fmt, ...) {
  va_list ap;
  fputs("[freegb] ", stdout);
  va_start(ap, fmt);
  vfprintf(stdout, fmt, ap);
  va_end(ap);
  fputc('\n', stdout);
  fflush(stdout);
}

/*------------------------------------------------------------------------*/

void die(int error_code, const char *fmt, ...) {
  fflush(stdout);
  va_list ap;
  fprintf(stderr, "*** [freegb] error code %i \n", error_code);
  fputs("*** [freegb] ", stderr);
  va_start(ap, fmt);
  vfprintf(stderr, fmt, ap);
  va_end(ap);
  fputc('\n', stderr);
  fflush(stderr);
  exit(error_code);
}

/*------------------------------------------------------------------------*/
// Global variables

long long hashmap_hits, hashmap_calls;
double amb_time, crit_pair_time, sym_pre_time, reduction_time;
double overlap_time, inclusion_time, crt_time, ratrec_time, rref_time;
double tt, other_time, overlap_map_time;

/*------------------------------------------------------------------------*/

size_t maximum_resident_set_size() {
  struct rusage u;
  if (getrusage(RUSAGE_SELF, &u)) return 0;
  return((size_t) u.ru_maxrss) << 10;
}

/*------------------------------------------------------------------------*/

double process_time() {
  struct rusage u;
  if (getrusage(RUSAGE_SELF, &u)) return 0;
  double res = u.ru_utime.tv_sec + 1e-6 * u.ru_utime.tv_usec;
  res += u.ru_stime.tv_sec + 1e-6 * u.ru_stime.tv_usec;
  return res;
}

static double percent(unsigned a, unsigned b) { return b ? 100.0*a/b : 0; }
static double percent_d(double a, double b) { return b ? 100.0*a/b : 0; }
/*------------------------------------------------------------------------*/

void print_statistics() {
  auto total_time = process_time();
  msg("");
  msg("maximum resident set size:  %22.2f MB",
  maximum_resident_set_size() / static_cast<double>((1<<20)));
  msg("monomial hashmap hitrate: %22.2f %%", percent(hashmap_hits,hashmap_calls)); 
  msg("computing ambiguities:    %22.2f (%2.2f %%)", amb_time, percent_d(amb_time,total_time));
  msg("  computing overlaps:     %22.2f (%2.2f %%)", overlap_time, percent_d(overlap_time,total_time));
  msg("  computing inclusions:   %22.2f (%2.2f %%)", inclusion_time, percent_d(inclusion_time,total_time));
  msg("handling critical pairs:  %22.2f (%2.2f %%)", crit_pair_time, percent_d(crit_pair_time,total_time));
  msg("symbolic preprocessing:   %22.2f (%2.2f %%)", sym_pre_time, percent_d(sym_pre_time,total_time));
  msg("linear algebra:           %22.2f (%2.2f %%)", reduction_time, percent_d(reduction_time,total_time));
  msg("  other:                  %22.2f (%2.2f %%)", other_time, percent_d(other_time,total_time));
  msg("  rref:                   %22.2f (%2.2f %%)", rref_time, percent_d(rref_time,total_time));
  msg("  CRT:                    %22.2f (%2.2f %%)", crt_time, percent_d(crt_time,total_time));
  msg("  rat. reconstruction:    %22.2f (%2.2f %%)", ratrec_time, percent_d(ratrec_time,total_time));
  msg("");
  msg("total process time:       %22.2f seconds",  process_time());
}
