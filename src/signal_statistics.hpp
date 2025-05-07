#ifndef SIGNAL_STATISTICS_H_
#define SIGNAL_STATISTICS_H_

#include <cstdarg>
#include <cstddef>
#include <csignal>
#include <sys/time.h>
#include <sys/resource.h>

extern void(*original_SIGINT_handler)(int);
extern void(*original_SIGSEGV_handler)(int);
extern void(*original_SIGABRT_handler)(int);
extern void(*original_SIGTERM_handler)(int);

/*------------------------------------------------------------------------*/
// / Level of output verbosity, ranges from 0 to 4
extern int verbose;

/**
    Prints an error message to stderr and exits the program

    @param char* fmt message
    @param int error_code
*/
void die(int error_code, const char *fmt,...);

/**
    Prints a message to stdout

    @param char* fmt message
*/
void msg(const char *fmt, ...);
/*------------------------------------------------------------------------*/

/**
    Determines max used memory
*/
size_t maximum_resident_set_size();

/**
    Determines the used process time
*/
double process_time();

/**
    Print statistics of maximum memory and used process time depending on
    selected modus

    @param modus integer
*/
void print_statistics();

#endif  // SIGNAL_STATISTICS_H_
