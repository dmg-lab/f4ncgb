/*------------------------------------------------------------------------*/
/*! \file signal_statistics.cpp
    \brief used to handle signals, messages and statistics

  Part of FGLMultiLing
  Copyright(C) 2024 Daniela Kaufmann,TU Wien, Austria
*/
/*------------------------------------------------------------------------*/

#include <cstdio>
#include <cstdlib>

#include "signal_statistics.hpp"

/*------------------------------------------------------------------------*/
// Global variable
void(*original_SIGINT_handler)(int);
void(*original_SIGSEGV_handler)(int);
void(*original_SIGABRT_handler)(int);
void(*original_SIGTERM_handler)(int);
/*------------------------------------------------------------------------*/
const char * bench = 0;

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

