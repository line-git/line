#ifndef GLOBAL_VARS_H_
#define GLOBAL_VARS_H_

// tolerance_digits limits the number of digits of the reationals
// used for the eta_path
extern int workprec;
extern int wp2;
// tolerance is used to decide whether something is zero
extern mpfr_t mpfr_tol;
// privatize global variables for parallelization
#pragma omp threadprivate(mpfr_tol)

// used for cache
extern char* filepath_cache;

// used for debug
extern int dbg;
#include <unistd.h>  // just for sleep function (debug)
extern double sleep_time;

// useful
extern FILE *dev_null_fptr;

// macros
#define MAX_PATH_LEN 1000
#define MAX_LEAF_LEN 256

// // used when finding numerical roots
// #define MAXM 100
// #define MR 8
// #define MT 10000
// #define MT_red 10
// #define MAXIT (MT*MR)
// #define MAXIT_red (MT_red*MR)
// extern double laguer_frac[MR+1];

// used when rationalizing with continued fractions
#define MAX_NUM_CONT_FRAC 100

#endif /* GLOBAL_VARS_H_ */

