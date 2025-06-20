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
const char * bench = 0;

/*------------------------------------------------------------------------*/
// Global variable
int verbose = 1;
int msg_printing = 1;

/*------------------------------------------------------------------------*/

void
msg(const char* fmt, ...) {
  if(!msg_printing)
    return;

  va_list ap;
  fputs("[f4ncgb] ", stderr);
  va_start(ap, fmt);
  vfprintf(stderr, fmt, ap);
  va_end(ap);
  fputc('\n', stderr);
  fflush(stderr);
}

/*------------------------------------------------------------------------*/

void die(int error_code, const char *fmt, ...) {
  fflush(stdout);
  va_list ap;
  fprintf(stderr, "*** [f4ncgb] error code %i \n", error_code);
  fputs("*** [f4ncgb] ", stderr);
  va_start(ap, fmt);
  vfprintf(stderr, fmt, ap);
  va_end(ap);
  fputc('\n', stderr);
  fflush(stderr);
  exit(error_code);
}

