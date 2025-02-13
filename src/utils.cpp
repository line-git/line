
#include <iostream>
#include <cstring>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
// #include <fstream>
// #include <sstream>
// #include <string>
#include <cmath>
// #include <complex>
#include "mpc.h"

using namespace std;

#include "tensor_utils.h"
#include "conversions.h"
#include "global_vars.h"
#include "malloc_defs.h"
#include "poly_frac.h"
#include "topology.h"

extern "C" {
  #include "cpoly.h"
  #include "algebra.h"
  #include "rel_err_mpc.h"
}


double clock_t_to_double(
  clock_t start, clock_t end
) {
  return ((double) end - start) / CLOCKS_PER_SEC;
}


double timespec_to_double(
  timespec start, timespec end
) {
  return end.tv_sec - start.tv_sec + (end.tv_nsec - start.tv_nsec)/1e9;
}


char* format_time(double seconds) {
  int minutes = (int)(seconds / 60);
  double remaining_seconds = fmod(seconds, 60);

  char* formatted_time = (char*) malloc(32*sizeof(char));
  snprintf(formatted_time, 32, "%6dm %09.6fs", minutes, remaining_seconds);

  return formatted_time;
}


void log_time_stats(
  double time_el_main,
  double *time_el_kira,
  double time_line, double time_el_line,
  double time_DE_preproc,
  double time_eps_loop, double time_el_eps_loop,
  double time_el_eps_iter, double time_el_DE, double time_el_prop,
  double time_el_regular_avg, double time_el_singular_avg,
  double time_el_normalize_avg,
  int eps_num, int nthreads, int opt_kira_parallel,
  int neta_values, int nsings
) {
  printf("\n------------------------------------------------------------\n");
  printf("\nTIME STATS\n");
  printf("- %d eps value(s), %d thread(s)", eps_num, nthreads);
  if (time_el_kira[0] > 0) {
  printf(", %d Kira thread(s)", opt_kira_parallel);
  }
  printf("\n");
  printf("- %d path points: %d regular, %d singular\n", neta_values, neta_values - nsings, nsings);

  printf("\n");
  if (time_el_kira[0] > 0) {
  printf("%-26s %s\n", "Kira (elapsed):", format_time(time_el_kira[0]+time_el_kira[1]));
  }
  if (nthreads > 1) {
  printf("%-26s %s, %s\n", "LINE (elapsed, cpu):", format_time(time_el_line), format_time(time_line));
  } else {
  printf("%-26s %s\n", "LINE:", format_time(time_el_line));
  }
  if (time_DE_preproc) {
  printf("%-26s %s\n", "  DE preproc:", format_time(time_DE_preproc));
  }
  if (nthreads > 1) {
  printf("%-26s %s, %s\n", "  eps loop (elapsed, cpu):", format_time(time_el_eps_loop), format_time(time_eps_loop));
  } else {
  printf("%-26s %s\n", "  eps loop:", format_time(time_eps_loop));
  }
  printf("%-26s %s\n", "    eps iteration:", format_time(time_el_eps_iter));
  if (time_el_DE > 0) {
  printf("%-26s %s\n", "      DE:", format_time(time_el_DE));
  }
  if (time_el_prop > 0) {
    printf("%-26s %s\n", "      propagation:", format_time(time_el_prop));
    printf("%-26s %s\n", "        regular:", format_time(time_el_regular_avg));
    printf("%-26s %s\n", "        singular:", format_time(time_el_singular_avg));
    printf("%-26s %s\n", "        normalize:", format_time(time_el_normalize_avg));
  }
  printf("\n");
  printf("%-26s %s\n", "total (elapsed):", format_time(time_el_main));

}


void generate_cache_filename(
  char **filename
) {
  time_t t;
  srand((unsigned) (time(&t) ^ getpid())); // combine time and PID for the seed

  int timestamp = (int)time(NULL);
  int random_component = rand() % 9000 + 1000; // random number between 1000 and 9999

  if (*filename) free(*filename);
  *filename = (char*) malloc(30 * sizeof(char));
  if (*filename == NULL) {
    perror("Failed to allocate memory for cache file name.");
    exit(EXIT_FAILURE);
  }

  snprintf(*filename, 30, "cache_%d_%d/", timestamp, random_component);
}


void print_result(
  FILE *resfptr, int precision,
  mpc_t **res, LI *MI, int dim, int order, int nloops
) {
  char res_format[20];
  snprintf(res_format, sizeof(res_format), "%%.%dRe", precision-1);
  int sign;
  for (int i=0; i<dim; i++) {
    fprintf(resfptr, "n. %d", i);
    if (MI) fprintf(resfptr, ", MI%s", MI[i].pows_str);
    fprintf(resfptr, "\n");
    for (int ep=0; ep<=order; ep++) {
      fprintf(resfptr, "eps^%d: ", ep - 2*nloops);
      if (ep - 2*nloops >= 0) fprintf(resfptr, " ");
      fprintf(resfptr, "(");

      // print real part
      sign = mpfr_cmp_ui(mpc_realref(res[i][ep]), 0);
      if (sign == 0) {
        fprintf(resfptr, "0");
      } else {
        if (sign > 0) fprintf(resfptr, "+");
        mpfr_fprintf(resfptr, res_format, mpc_realref(res[i][ep]));
        // mpfr_out_str(resfptr, 10, 0, mpc_realref(res[i][ep]), MPFR_RNDN);
      }

      fprintf(resfptr, " ");

      // print imag part
      sign = mpfr_cmp_ui(mpc_imagref(res[i][ep]), 0);
      if (sign == 0) {
        fprintf(resfptr, "0");
      } else {
        if (sign > 0) fprintf(resfptr, "+");
        mpfr_fprintf(resfptr, res_format, mpc_imagref(res[i][ep]));
        // mpfr_out_str(resfptr, 10, 0, mpc_imagref(res[i][ep]), MPFR_RNDN);
      }
      fprintf(resfptr, ")\n");
    }
    fprintf(resfptr, "\n");
  }
}


// void write_to_cache(const char* data) {
//   char* filename = generate_cache_filename();
//   FILE* cache_file = fopen(filename, "w");
//   if (cache_file == NULL) {
//     perror("Failed to open file");
//     free(filename);
//     exit(EXIT_FAILURE);
//   }

//   fprintf(cache_file, "%s", data);
//   fclose(cache_file);
//   printf("Data written to %s\n", filename);

//   free(filename); // Free the allocated memory
// }


int MIN(int a, int b) {
  return ((a) < (b) ? a : b);
}


inline int MAX(int a, int b) {
  return ((a) > (b) ? a : b);
}


int int_pow_int(int base, int exp) {
  int out = 1;
  for (int k=0; k<exp; k++) {
    out *= base;
  }
  return out;
}


int count_lines(char *filename) {
  FILE *fptr = fopen(filename, "r");
	int count = 0;
	char *line = NULL;
	size_t line_len = 0;
	while (getline(&line, &line_len, fptr) != -1) {
		if (line[0] != '\n' && line[0] != '\r') {
			count++;
		}
	}
  fclose(fptr);
  return count;
}


int count_lines_fptr(FILE *fptr) {
	int count = 0;
	char *line = NULL;
	size_t line_len = 0;
	while (getline(&line, &line_len, fptr) != -1) {
		if (line[0] != '\n' && line[0] != '\r') {
			count++;
		}
	}
  return count;
}


void join_path(
  // OUTPUT
  char **path,
  // INPUT
  char *parent, char *child
) {
  // cout << "parent: " << parent << endl;
  // cout << "child: " << child << endl;
  char *in1, *in2;
  int overwrt;
  // deal with overwrite
  if (*path == parent) {
    // cout << "overwrite 1" << endl;
    overwrt = 1;
    in1 = (char*) malloc((strlen(parent)+1)*sizeof(char));
    strcpy(in1, parent);
    in2 = child;
    (*path)[0] = '\0';
    // cout << "in1: " << in1 << endl;
    // cout << "in2: " << in2 << endl;
  } else if (*path == child) {
    // cout << "overwrite 2" << endl;
    overwrt = 2;
    in1 = parent;
    in2 = (char*) malloc((strlen(child)+1)*sizeof(char));
    strcpy(in2, child);
    (*path)[0] = '\0';
    // cout << "in1: " << in1 << endl;
    // cout << "in2: " << in2 << endl;
  } else {
    overwrt = 0;
    in1 = parent;
    in2 = child;
  }

  (*path) = (char*) realloc(*path, (strlen(in1)+strlen(in2)+1)*sizeof(char));  
  strcpy(*path, in1);
  strcat(*path, in2);
  // cout << "result: " << *path << endl;

  if (overwrt == 1) {
    free(in1);
  } else if (overwrt == 2) {
    free(in2);
  }
}


int mpfr_lessthan_tol(mpfr_t in) {
  if (mpfr_cmpabs(in, mpfr_tol) < 0) {
    return 1;
  } else {
    return 0;
  }
}


int mpfr_sign_within_tol(mpfr_t in) {
  if (mpfr_cmpabs(in, mpfr_tol) < 0) {
    return 0;
  } else {
    return mpfr_cmp_ui(in, 0);
  }
}


int mpc_lessthan_tol(mpc_t in) {
  if (mpfr_cmpabs(mpc_realref(in), mpfr_tol) < 0\
    && mpfr_cmpabs(mpc_imagref(in), mpfr_tol) < 0) {
    return 1;
  } else {
    return 0;
  }
}


int mpfr_equal_within_tol(
  mpfr_t in1, mpfr_t in2
) {
  // if (dbg) {
  //   cout << "in1 = "; mpfr_out_str(stdout, 10, 0, in1, MPFR_RNDN); cout << endl;
  //   cout << "in2 = "; mpfr_out_str(stdout, 10, 0, in2, MPFR_RNDN); cout << endl;
  //   cout << "mpfr_cmp = " << mpfr_cmp(in1, in2) << endl;
  // }
  if (mpfr_cmp(in1, in2) == 0) {
    return 1;
  }
  mpfr_t tmpfr;
  mpfr_init2(tmpfr, wp2);
  mpfr_sub(tmpfr, in1, in2, MPFR_RNDN);
  if (mpfr_cmpabs(in1, mpfr_tol) > 0) {
    mpfr_div(tmpfr, tmpfr, in1, MPFR_RNDN);
  }
  if (mpfr_cmpabs(tmpfr, mpfr_tol) < 0) {
    mpfr_clear(tmpfr);
    return 1;
  } else {
    mpfr_clear(tmpfr);
    return 0;
  }
}


int mpc_equal_within_tol(
  mpc_t in1, mpc_t in2
) {
  // if (dbg) {
  //   cout << "in1 = "; mpc_out_str(stdout, 10, 0, in1, MPFR_RNDN); cout << endl;
  //   cout << "in2 = "; mpc_out_str(stdout, 10, 0, in2, MPFR_RNDN); cout << endl;
  // }
  if (mpfr_equal_within_tol(mpc_realref(in1), mpc_realref(in2))) {
    if (mpfr_equal_within_tol(mpc_imagref(in1), mpc_imagref(in2))) {
      return 1;
    } else {
      return 0;
    }
  } else {
    return 0;
  }
}


int mpc_zero_p(
  mpc_t in
) {
  return (mpfr_zero_p(mpc_realref(in)) && mpfr_zero_p(mpc_imagref(in)));
}


void mpc_rk2_prune_abs_tol(
  mpc_t **tens,
  int dim1, int dim2
) {
  for (int i1=0; i1<dim1; i1++) {
    for (int i2=0; i2<dim2; i2++) {
      if (mpfr_lessthan_tol(mpc_realref(tens[i1][i2]))) {
        mpfr_set_ui(mpc_realref(tens[i1][i2]), 0, MPFR_RNDN);
      }
      if (mpfr_lessthan_tol(mpc_imagref(tens[i1][i2]))) {
        mpfr_set_ui(mpc_imagref(tens[i1][i2]), 0, MPFR_RNDN);
      }
    }
  }
}


int max_ints(int *arr, int dim) {
  int max = arr[0];
  for (int i=1; i<dim; i++) {
    if (arr[i] > max) {
      max = arr[i];
    }
  }
  return max;
}


void minmax0(double *min, double *max, double *arr, int MAX_SIZE) {
  // compute max and min for arrays of unique doubles, where
  // the min if found among non-vanishing doubles
  int i;

  *max = arr[0];
  *min = arr[0];
  if (*min == 0)
    *min = arr[1];
  for(i=1; i<MAX_SIZE; i++) {
    if(arr[i] > *max) {
      *max = arr[i];
    }
    if(arr[i] < *min & arr[i] != 0) {
      *min = arr[i];
    }
  }
}


void minmax0_mp(mpfr_t *min, mpfr_t *max, mpfr_t *arr, int MAX_SIZE) {
  // compute max and min for arrays of unique mpfr_t, where
  // the min is found among non-vanishing values
  int i;

  mpfr_set(*max, arr[0], MPFR_RNDN);
  mpfr_set(*min, arr[0], MPFR_RNDN);
  if (mpfr_lessthan_tol(*min))
    mpfr_set(*min, arr[1], MPFR_RNDN);
  for(i=1; i<MAX_SIZE; i++) {
    if(mpfr_cmp(arr[i], *max) > 0) {
      mpfr_set(*max, arr[i], MPFR_RNDN);
    }
    if(mpfr_cmp(arr[i], *min) < 0 && !mpfr_lessthan_tol(arr[i])) {
      mpfr_set(*min, arr[i], MPFR_RNDN);
    }
  }
}


double min1(double *arr, int MAX_SIZE) {
  // compute the min among non-vanishing elements of an arrays 
  // of unique doubles
  int i, j;
  double min = 0;
  for(i=0; i<MAX_SIZE; i++) {
    if (arr[i] != 0) {
      min = arr[i];
      break;
    }
  }
  if (min == 0) {
    fprintf(stderr, "non-zero minimum does not exit in input array\n");
    exit(EXIT_FAILURE);
  }

  for(j=i+1; j<MAX_SIZE; j++) {
    if(arr[j] < min & arr[j] != 0) {
      min = arr[j];
    }
  }
  return min;
}


void min1_mp(mpfr_t *min, mpfr_t *arr, int MAX_SIZE) {
  // compute the min among non-vanishing elements of an array 
  // of unique mpfr_t
  int i, j;
  mpfr_set_si(*min, 0, MPFR_RNDN);
  for(i=0; i<MAX_SIZE; i++) {
    if (!mpfr_lessthan_tol(arr[i])) {
      mpfr_set(*min, arr[i], MPFR_RNDN);
      break;
    }
  }
  if (mpfr_lessthan_tol(*min)) {
    fprintf(stderr, "non-zero minimum does not exist in input array\n");
    exit(EXIT_FAILURE);
  }

  for(j=i+1; j<MAX_SIZE; j++) {
    if(mpfr_cmp(arr[j], *min) < 0 && !mpfr_lessthan_tol(arr[j])) {
      mpfr_set(*min, arr[j], MPFR_RNDN);
    }
  }
}


void limit_denominator(
  // OUTPUT
  mpq_t out_rat,
  // INPUT
  mpq_t rat, mpz_t max_denominator
) {
  mpz_t num, den, p0, q0, p1, q1, n, d, k;

  mpz_init(num);
  mpz_init(den);
  mpz_init(p0);
  mpz_init(q0);
  mpz_init(p1);
  mpz_init(q1);
  mpz_init(n);
  mpz_init(d);
  mpz_init(k);
  
  mpq_get_num(num, rat);
  mpq_get_den(den, rat);

  if (mpz_cmp(den, max_denominator) <= 0) {
    mpq_set(out_rat, rat);
    mpz_clear(num);
    mpz_clear(den);
    mpz_clear(p0);
    mpz_clear(q0);
    mpz_clear(p1);
    mpz_clear(q1);
    mpz_clear(n);
    mpz_clear(d);
    mpz_clear(k);
    return;
  }

  mpz_set_ui(p0, 0); // p0 = 0
  mpz_set_ui(q0, 1); // q0 = 1
  mpz_set_ui(p1, 1); // p1 = 1
  mpz_set_ui(q1, 0); // q1 = 0
  mpz_set(n, num);
  mpz_set(d, den);

  mpz_t a; mpz_init(a);
  mpz_t q2; mpz_init(q2);
  mpz_t t; mpz_init(t);
  mpq_t bound1; mpq_init(bound1);
  mpq_t bound2; mpq_init(bound2);
  mpq_t tmp1; mpq_init(tmp1);
  mpq_t tmp2; mpq_init(tmp2);
  int num_it = 0;
  while (1) {
    num_it++;
    mpz_tdiv_q(a, n, d); // a = floor(n / d)
    mpz_mul(q2, a, q1);
    mpz_add(q2, q2, q0);
    if (mpz_cmp(q2, max_denominator) > 0) {
      break;
    }
    
    // aggiornamento delle variabili
    mpz_set(t, p0); // t = p0
    mpz_set(p0, p1); // p0 = p1
    mpz_set(q0, q1); // q0 = q1
    mpz_mul(p1, p1, a);
    mpz_add(p1, p1, t);
    mpz_set(q1, q2); // q1 = q2

    // scambia n e d
    mpz_set(t, n); // t = n
    mpz_set(n, d); // n = d
    mpz_mul(d, a, d);
    mpz_sub(d, t, d);
  }

  // calcola k
  mpz_sub(k, max_denominator, q0);
  mpz_tdiv_q(k, k, q1); // k = floor((max_denominator - q0) / q1)

  mpq_set_num(bound1, p0);
  mpz_addmul(mpq_numref(bound1), k, p1);
  mpq_set_den(bound1, q0);
  mpz_addmul(mpq_denref(bound1), k, q1);
  mpq_set_num(bound2, p1);
  mpq_set_den(bound2, q1);

  // cout << "num_it = " << num_it << endl;
  // cout << "bound1 = "; mpq_out_str(stdout, 10, bound1); cout << endl;
  // cout << "bound2 = "; mpq_out_str(stdout, 10, bound2); cout << endl;

  // confronta le differenze
  mpq_sub(tmp1, bound1, rat);
  mpq_sub(tmp2, bound2, rat);
  mpq_abs(tmp1, tmp1);
  mpq_abs(tmp2, tmp2);
  if (mpq_cmp(tmp2, tmp1) <= 0) {
    mpq_set(out_rat, bound2);
  } else {
    mpq_set(out_rat, bound1);
  }

  // libera la memoria
  mpz_clear(num);
  mpz_clear(den);
  mpz_clear(p0);
  mpz_clear(q0);
  mpz_clear(p1);
  mpz_clear(q1);
  mpz_clear(n);
  mpz_clear(d);
  mpz_clear(k);
  mpz_clear(a);
  mpz_clear(q2);
  mpz_clear(t);
  mpq_clear(bound1);
  mpq_clear(bound2);
  mpq_clear(tmp1);
  mpq_clear(tmp2);
}


void rationalize(mpq_t out_rat, double x, double tolerance) {
  mpq_t rat; mpq_init(rat);
  mpq_set_d(rat, x);

  mpz_t max_den; mpz_init(max_den);
  mpz_set_d(max_den, 1.0/tolerance);
  mpz_fdiv_q_ui(max_den, max_den, 1);
  // cout << "rat = "; mpq_out_str(stdout, 10, rat); cout << endl;
  // cout << "mpz max den = "; mpz_out_str(stdout, 10, max_den); cout << endl;  
  limit_denominator(out_rat, rat, max_den);
  // cout << "out_rat = "; mpq_out_str(stdout, 10, out_rat); cout << endl;
  if (x<0) {
    mpq_neg(out_rat, out_rat);
  }
}


void find_cf(double x, double eps) {
  long p[MAX_NUM_CONT_FRAC], q[MAX_NUM_CONT_FRAC], a[MAX_NUM_CONT_FRAC], len;
  int i;
  double approx, err, tmp = x;
  //The first two convergents are 0/1 and 1/0
  p[0] = 0; q[0] = 1;
  p[1] = 1; q[1] = 0;
  //The rest of the convergents (and continued fraction)
  for(i=2; i<MAX_NUM_CONT_FRAC; ++i) {
    printf("i = %d\n", i);
    a[i] = lrint(floor(tmp));
    p[i] = a[i]*p[i-1] + p[i-2];
    q[i] = a[i]*q[i-1] + q[i-2];
    printf("%ld:  %ld/%ld\n", a[i], p[i], q[i]);
    approx = ((double) p[i])/((double) q[i]);
    err = fabs(x - approx);
    printf("approx = %.16e\n", approx);
    printf("error = %.16e\n", err);
    if(err < eps) return;
    tmp = 1.0/(tmp - a[i]);
  }
}


void rationalize_cf(mpq_t rat, mpfr_t x, mpfr_t delta) {
  int sgn = mpfr_cmp_si(x, 0);
  if (sgn == 0) {
    mpq_set_ui(rat, 0, 1);
    return;
  }

  mpz_t p[MAX_NUM_CONT_FRAC], q[MAX_NUM_CONT_FRAC], a[MAX_NUM_CONT_FRAC];
  int i;
  mpfr_t approx, err, tmp;
  mpfr_init2(approx, wp2);
  mpfr_init2(err, wp2);
  mpfr_init2(tmp, wp2);

  mpfr_abs(tmp, x, MPFR_RNDN);
  if (sgn < 0) {
    mpfr_neg(tmp, x, MPFR_RNDN);
    mpfr_neg(x, x, MPFR_RNDN);
  } else {
    mpfr_set(tmp, x, MPFR_RNDN);
  }

  //The first two convergents are 0/1 and 1/0
  mpz_init(p[0]); mpz_set_ui(p[0], 0);
  mpz_init(q[0]); mpz_set_ui(q[0], 1);
  mpz_init(p[01]); mpz_set_ui(p[1], 1);
  mpz_init(q[1]); mpz_set_ui(q[1], 0);

  //The rest of the convergents (and continued fraction)
  for(i=2; i<MAX_NUM_CONT_FRAC; i++) {
    // printf("i = %d\n", i);
    // a[i] = lrint(floor(tmp));
    mpz_init(a[i]);
    mpfr_get_z(a[i], tmp, MPFR_RNDN);
    // p[i] = a[i]*p[i-1] + p[i-2];
    mpz_init(p[i]);
    mpz_mul(p[i], a[i], p[i-1]);
    mpz_add(p[i], p[i], p[i-2]);
    // q[i] = a[i]*q[i-1] + q[i-2];
    mpz_init(q[i]);
    mpz_mul(q[i], a[i], q[i-1]);
    mpz_add(q[i], q[i], q[i-2]);

    // cout << "frac = ";
    // mpz_out_str(stdout, 10, p[i]); cout << "/";
    // mpz_out_str(stdout, 10, q[i]); cout << endl;

    // approx = ((double) p[i])/((double) q[i]);
    mpfr_set_z(approx, p[i], MPFR_RNDN);
    mpfr_div_z(approx, approx, q[i], MPFR_RNDN);
    // cout << "approx = "; mpfr_out_str(stdout, 10, 0, approx, MPFR_RNDN); cout << endl;
    // err = fabs(x - approx);
    mpfr_sub(err, x, approx, MPFR_RNDN);
    mpfr_abs(err, err, MPFR_RNDN);
    // cout << "err = "; mpfr_out_str(stdout, 10, 0, err, MPFR_RNDN); cout << endl;
    
    if (mpfr_less_p(err, delta)) break;

    // tmp = 1.0/(tmp - a[i]);
    mpfr_sub_z(tmp, tmp, a[i], MPFR_RNDN);
    mpfr_ui_div(tmp, 1, tmp, MPFR_RNDN);
    
  }
  if (i == MAX_NUM_CONT_FRAC) {
    printf("reached max number of iterations in continued fractions algorithm\n");
    perror("reached max number of iterations in continued fractions algorithm");
    exit(1);
    i--;
  }

  mpz_set(mpq_numref(rat), p[i]);
  mpz_set(mpq_denref(rat), q[i]);
  mpq_canonicalize(rat);

  if (sgn < 0) {
    mpq_neg(rat, rat);
    mpfr_neg(x, x, MPFR_RNDN);
  }

  // FREE
  mpfr_clear(approx);
  mpfr_clear(err);
  mpfr_clear(tmp);
  
  for (int ip=0; ip<=i; ip++) {
    fflush(stdout);
    mpz_clear(p[ip]);
    mpz_clear(q[ip]);
    if (ip > 1) mpz_clear(a[ip]);
  }
}


void rationalize_cf_from_str(mpq_t rat, char *str) {
  int ndigits = 0, start_count = 0;
  for (int c=0; str[c] != '\0'; c++) {
    if (str[c] == 'e') break;
    if (start_count) {
      if (str[c] != '.') {
        ndigits++;
      }
    } else {
      if (str[c] != '.' && str[c] != '0' && str[c] != '-') {
        ndigits = 1;
        start_count = 1;
      }
    }
  }

  mpfr_t real; mpfr_init2(real, wp2);
  mpfr_t delta; mpfr_init2(delta, wp2);

  mpfr_set_str(real, str, 10, MPFR_RNDN);
  mp_exp_t exp10 = mpfr_get_exp(real) * log10(2);
  mpfr_set_ui(delta, 10, MPFR_RNDN);
  mpfr_pow_si(delta, delta, exp10 - ndigits + 1, MPFR_RNDN);

  // cout << "ndigits = " << ndigits << endl;
  // cout << "exp10 = " << exp10 << endl;
  // cout << "delta = 10^" << exp10 - ndigits + 1 << endl;
  // cout << "real = "; mpfr_out_str(stdout, 10, 0, real, MPFR_RNDN); cout << endl;
  // cout << "delta = "; mpfr_out_str(stdout, 10, 0, delta, MPFR_RNDN); cout << endl;

  rationalize_cf(rat, real, delta);

  mpfr_clear(real);
  mpfr_clear(delta);
}


void print_mpc(mpc_t *mpc_num) {
  mpc_out_str(stdout, 10, 0, *mpc_num, MPFR_RNDN);
}


void print_mpfr(mpfr_t *mpfr_num) {
  mpfr_out_str(stdout, 10, 0, *mpfr_num, MPFR_RNDN);
}


bool mpfr_is_diff_int(mpfr_t real1, mpfr_t real2) {
  mpfr_t diff, nearest_int;
  mpfr_init2(diff, wp2);
  mpfr_init2(nearest_int, wp2);
  mpfr_sub(diff, real1, real2, MPFR_RNDN);
  // round difference to the nearest integer
  mpfr_rint(nearest_int, diff, MPFR_RNDN);
  // subtract nearest integer
  mpfr_sub(diff, diff, nearest_int, MPFR_RNDN);
  // check whether what remain is compatible with zero within tolerance
  mpfr_abs(diff, diff, MPFR_RNDN);
  return mpfr_less_p(diff, mpfr_tol);
}


void mpc_swap(mpc_t *z1, mpc_t *z2) {
  mpc_t tmp;
  mpc_init3(tmp, wp2, wp2);
  mpc_set(tmp, *z1, MPFR_RNDN);
  mpc_set(*z1, *z2, MPFR_RNDN);
  mpc_set(*z2, tmp, MPFR_RNDN);
}


void mpfr_tol_enlarge(double wp2_rel_decr) {
  mpfr_t tmp;
  mpfr_init2(tmp, wp2);
  mpfr_set_ui(tmp, 2, MPFR_RNDN);
  mpfr_pow_si(mpfr_tol, tmp, -wp2_rel_decr*((int) wp2), MPFR_RNDN);
}


void mpfr_tol_enlarge_rel(double wp2_rel_decr) {
  mpfr_t pow; mpfr_init2(pow, wp2);
  mpfr_set_d(pow, wp2_rel_decr, MPFR_RNDN);
  mpfr_powr(mpfr_tol, mpfr_tol, pow, MPFR_RNDN);
}


int mpfr_log2_int(mpfr_t in) {
  mpfr_t log2_mpfr; mpfr_init2(log2_mpfr, wp2);
  mpfr_log2(log2_mpfr, in, MPFR_RNDN);
  int log2_int = mpfr_get_si(log2_mpfr, MPFR_RNDN);
  mpfr_clear(log2_mpfr);
  return log2_int;
}


size_t mpfr_get_memory_usage(mpfr_t x) {
	// sizeof(mpfr_t) restituisce la dimensione statica della struttura
	size_t static_size = sizeof(mpfr_t);

	// Calcola la memoria allocata dinamicamente per il significand
	size_t dynamic_size = (mpfr_get_prec(x) + 7) / 8; // byte necessari per la precisione

	return static_size + dynamic_size;
}


void interpolate_epsilon_orders(
  // OUTPUT
  mpc_t **sol_eps_ord,
  // INPUT
  mpc_t **sol_at_eps, char **eps_str,
  int dim, int eps_num, int loops,
  int precision,
  int *starting_ord
) {
  int wp_bin = (10*precision)/3;
  cout << "precision = " << precision << endl;
  cout << "wp_bin = " << wp_bin << endl;

  int pow_offset = 2*loops;
  int null_pow = 2*loops+1;
  if (starting_ord == NULL) {
    starting_ord = new int[dim];
    for (int m=0; m<dim; m++) {
      starting_ord[m] = 0;
    }
  }
  int *inv_computed = new int[null_pow];
  for (int k=0; k<null_pow; k++) {
    inv_computed[k] = 0;
  }
  mpc_t *mpc_eps_list;
  mpc_eps_list = new mpc_t[eps_num];
  init_rk1_mpc(mpc_eps_list, eps_num);
  mpc_t **eps_mat, **eps_mat_copy;
  malloc_rk2_tens(eps_mat, eps_num, eps_num + null_pow);
  init_rk2_mpc(eps_mat, eps_num, eps_num + null_pow);
  malloc_rk2_tens(eps_mat_copy, eps_num, eps_num);
  init_rk2_mpc(eps_mat_copy, eps_num, eps_num);
  mpc_t ***inv_eps_mat;
  malloc_rk3_tens(inv_eps_mat, null_pow, eps_num, eps_num);
  init_rk3_mpc(inv_eps_mat, null_pow, eps_num, eps_num);
  // ex ex_epv;
  for (int ep=0; ep<eps_num; ep++) {
    // ex_epv = (ex) (numeric) eps_list[ep];
    // gnc_to_mpc(&mpc_eps_list[ep], ex_epv*1.);
    // cout << "ep = " << ep << endl;
    // cout << "eps_str: " << eps_str[ep] << endl;
    mpc_set_str_rat(&mpc_eps_list[ep], eps_str[ep]);
    // cout << "eps = "; print_mpc(&mpc_eps_list[ep]); cout << endl;
    for (int k=0; k<eps_num+null_pow; k++) {
      mpc_pow_si(eps_mat[ep][k], mpc_eps_list[ep], -pow_offset+k, MPFR_RNDN);
    }
  }
  // mp_inverse(inv_eps_mat, eps_mat, eps_num);
  
  for (int m=0; m<dim; m++) {
    // cout << "m = " << m << endl;
    // cout << "starting order = " << starting_ord[m] << endl;
    if (inv_computed[starting_ord[m]] == 0) {
      for (int ep=0; ep<eps_num; ep++) {
        eps_mat[ep] += starting_ord[m];
      }
      copy_rk2_mpc(eps_mat_copy, eps_mat, eps_num, eps_num);
      // cout << "eps_mat:" << endl;
      // print_rk2_mpc(eps_mat_copy, eps_num, eps_num);
      mp_inverse(inv_eps_mat[starting_ord[m]], eps_mat_copy, eps_num);
      // cout << "inv_eps_mat:" << endl;
      // print_rk2_mpc(inv_eps_mat[starting_ord[m]], eps_num, eps_num);
      for (int ep=0; ep<eps_num; ep++) {
        eps_mat[ep] -= starting_ord[m];
      }
      inv_computed[starting_ord[m]] = 1;
    }
    // for (int k=0; k<null_pow; k++) {
    //   mpc_set_ui(sol_eps_ord[m][k], 0, MPFR_RNDN);
    // }
    mpc_rk2_mul_mpc_rk1(
      sol_eps_ord[m], inv_eps_mat[starting_ord[m]], sol_at_eps[m],
      eps_num, eps_num
    );
    // rel_err_mpc_rk2_mul_mpc_rk1(
    //   sol_eps_ord[m], inv_eps_mat[starting_ord[m]], sol_at_eps[m],
    //   eps_num, eps_num, wp_bin
    // );
    // cout << "orders:" << endl;
    // print_poly(sol_eps_ord[m], eps_num-1);
  }

  // check
  mpc_t mpc_check;
  mpc_init3(mpc_check, wp2, wp2);
  mpc_t mpc_diff;
  mpc_init3(mpc_diff, wp2, wp2);
  for (int m=0; m<dim; m++) {
    for (int ep=0; ep<eps_num; ep++) {
      // cout << "m, ep = " << m << ", " << ep << endl;
      // reconstruct series at selected point
      mpc_set(mpc_check, sol_eps_ord[m][eps_num-1], MPFR_RNDN);
      for (int epp=eps_num-1; epp>0; epp--) {
        mpc_fma(mpc_check, mpc_eps_list[ep], mpc_check, sol_eps_ord[m][epp-1], MPFR_RNDN);
      }
      mpc_pow_si(mpc_diff, mpc_eps_list[ep], -2*loops+starting_ord[m], MPFR_RNDN);
      mpc_mul(mpc_check, mpc_check, mpc_diff, MPFR_RNDN);
      // compare with pre-interpolation value
      mpc_sub(mpc_diff, mpc_check, sol_at_eps[m][ep], MPFR_RNDN);
      mpc_div(mpc_diff, mpc_diff, sol_at_eps[m][ep], MPFR_RNDN);
      // if (1) {
      if (mpc_lessthan_tol(mpc_diff)) {
        cout << "found difference at m, ep = " << m << ", " << ep << endl;
        cout << "pre-interpolation:" << endl;
        print_mpc(&sol_at_eps[m][ep]); cout << endl;
        cout << "post-interpolation:" << endl;
        print_mpc(&mpc_check); cout << endl;
        cout << "diff:" << endl;
        print_mpc(&mpc_diff); cout << endl;
      }      
    }
  }

  // FREE
  mpc_rk1_clear(mpc_eps_list, eps_num);
  delete[] mpc_eps_list;
  mpc_rk2_clear(eps_mat, eps_num, eps_num + null_pow);
  del_rk2_tens(eps_mat, eps_num);
  mpc_rk2_clear(eps_mat_copy, eps_num, eps_num);
  del_rk2_tens(eps_mat_copy, eps_num);
  mpc_rk3_clear(inv_eps_mat, null_pow, eps_num, eps_num);
  del_rk3_tens(inv_eps_mat, null_pow, eps_num);

}


void interpolate_epsilon_orders_prune(
  // OUTPUT
  mpc_t **sol_eps_ord,
  // INPUT
  mpc_t **sol_at_eps, char **eps_str,
  int dim, int eps_num, int loops,
  int precision, int order,
  int *starting_ord
) {
  // cout << "precision = " << precision << endl;
  // cout << "wp_bin = " << wp_bin << endl;

  int pow_offset = 2*loops;
  int null_pow = 2*loops+1;
  for (int m=0; m<dim; m++) {
    starting_ord[m] = 0;
  }
  int *inv_computed = new int[null_pow];
  for (int k=0; k<null_pow; k++) {
    inv_computed[k] = 0;
  }
  mpc_t *mpc_eps_list;
  mpc_eps_list = new mpc_t[eps_num];
  init_rk1_mpc(mpc_eps_list, eps_num);
  mpc_t **eps_mat, **eps_mat_copy;
  malloc_rk2_tens(eps_mat, eps_num, eps_num + null_pow);
  init_rk2_mpc(eps_mat, eps_num, eps_num + null_pow);
  malloc_rk2_tens(eps_mat_copy, eps_num, eps_num);
  init_rk2_mpc(eps_mat_copy, eps_num, eps_num);
  mpc_t ***inv_eps_mat;
  malloc_rk3_tens(inv_eps_mat, null_pow, eps_num, eps_num);
  init_rk3_mpc(inv_eps_mat, null_pow, eps_num, eps_num);
  mpc_t ***inv_eps_mat_red;
  malloc_rk3_tens(inv_eps_mat_red, null_pow, eps_num-1, eps_num-1);
  init_rk3_mpc(inv_eps_mat_red, null_pow, eps_num-1, eps_num-1);
  for (int ep=0; ep<eps_num; ep++) {
    // cout << "ep = " << ep << endl;
    // cout << "eps_str: " << eps_str[ep] << endl;
    mpc_set_str_rat(&mpc_eps_list[ep], eps_str[ep]);
    // cout << "eps = "; print_mpc(&mpc_eps_list[ep]); cout << endl;
    for (int k=0; k<eps_num+null_pow; k++) {
      mpc_pow_si(eps_mat[ep][k], mpc_eps_list[ep], -pow_offset+k, MPFR_RNDN);
    }
  }
  // mp_inverse(inv_eps_mat, eps_mat, eps_num);

  // get avarage base 2 exponent of epsilon
  int eps_exp = mpfr_get_exp(mpc_realref(mpc_eps_list[0]));
  for (int ep=1; ep<eps_num; ep++) {
    eps_exp += mpfr_get_exp(mpc_realref(mpc_eps_list[ep]));
  }
  eps_exp /= -eps_num;
  // cout << "epsilon avarage base 2 exponent = " << eps_exp << endl;
  // cout << "epsilon avarage base 10 exponent = " << eps_exp * log10(2) << endl;

  int wp_bin = log2(10)*precision + (order+1)*eps_exp;
  // cout << "lo wp_bin = " << wp_bin << endl;
  // cout << "lo wp_bin (base 10) = " << wp_bin * log10(2) << endl;
  
  // compute inverse for every offset
  for (int s=0; s<null_pow; s++) {
    // cout << "s = " << s << endl;
    copy_rk2_mpc(eps_mat_copy, eps_mat, eps_num, eps_num);
    // if (s == 0) {
    //   cout << "eps_mat:" << endl;
    //   print_rk2_mpc(eps_mat_copy, eps_num, eps_num);
    // }
    mp_inverse(inv_eps_mat[s], eps_mat_copy, eps_num);
    // cout << "inv_eps_mat:" << endl;
    // print_rk2_mpc(inv_eps_mat[s], eps_num, eps_num);
    copy_rk2_mpc(eps_mat_copy, eps_mat, eps_num-1, eps_num-1);
    mp_inverse(inv_eps_mat_red[s], eps_mat_copy, eps_num-1);
    // cout << "inv_eps_mat_red:" << endl;
    // print_rk2_mpc(inv_eps_mat_red[s], eps_num-1, eps_num-1);
    for (int ep=0; ep<eps_num; ep++) {
      eps_mat[ep]++;
    }
  }
  for (int ep=0; ep<eps_num; ep++) {
    eps_mat[ep] -= null_pow;
  }


  mpc_t lo; mpc_init3(lo, wp2, wp2);
  mpc_t *sol_eps_ord_red = new mpc_t[eps_num];
  init_rk1_mpc(sol_eps_ord_red, eps_num);
  for (int m=0; m<dim; m++) {
    // cout << "m = " << m << endl;    
    // INTERPOLATE AND FIND LEADING ORDER

    // compute epsilon orders with no shift
    mpc_rk2_mul_mpc_rk1(
      sol_eps_ord[m], inv_eps_mat[0], sol_at_eps[m],
      eps_num, eps_num
    );
    
    for (int s=0; s<null_pow-1; s++) {
      // cout << "s = " << s << endl;
      // compute lo with fewer epsilons
      mpc_set_ui(lo, 0, MPFR_RNDN);
      for (int k=0; k<eps_num-1; k++) {
        mpc_fma(lo, inv_eps_mat_red[s][0][k], sol_at_eps[m][k], lo, MPFR_RNDN);
      }

      // check for cancellation of lo
      // cout << "lo      = "; print_mpc(&sol_eps_ord[m][0]); cout << endl;
      // cout << "lo red  = "; print_mpc(&lo); cout << endl;
      mpc_neg(lo, lo, MPFR_RNDN);
      exp_rel_err_mpc_add(lo, sol_eps_ord[m][0], lo, MPFR_RNDN, wp_bin - eps_exp*s);
      // cout << "lo diff = "; print_mpc(&lo); cout << endl;
      
      if (
        (mpfr_zero_p(mpc_realref(lo)) && !mpfr_zero_p(mpc_realref(sol_eps_ord[m][0]))) || \
        (mpfr_zero_p(mpc_imagref(lo)) && !mpfr_zero_p(mpc_imagref(sol_eps_ord[m][0])))
      ) {
        // lo is stable
        break;
      }

      // shift leading power
      starting_ord[m]++;
      mpc_rk2_mul_mpc_rk1(
        sol_eps_ord[m], inv_eps_mat[s+1], sol_at_eps[m],
        eps_num, eps_num
      );
    }

    // PRUNE REAL AND IMAGINARY PARTS
    
    // compute orders with fewer epsilons
    mpc_rk2_mul_mpc_rk1(
      sol_eps_ord_red, inv_eps_mat_red[starting_ord[m]], sol_at_eps[m],
      eps_num-1, eps_num-1
    );

    // prune real and imaginary parts separately
    for (int k=0; k<eps_num-1; k++) {
      // check for cancellation
      mpc_neg(sol_eps_ord_red[k], sol_eps_ord_red[k], MPFR_RNDN);
      exp_rel_err_mpc_add(lo, sol_eps_ord[m][k], sol_eps_ord_red[k], MPFR_RNDN, wp_bin - eps_exp*(starting_ord[m]+k));
      if (!mpfr_zero_p(mpc_realref(lo))) {
        // prune real part
        mpfr_set_ui(mpc_realref(sol_eps_ord[m][k]), 0, MPFR_RNDN);
      }
      if (!mpfr_zero_p(mpc_imagref(lo))) {
        // prune imag part
        mpfr_set_ui(mpc_imagref(sol_eps_ord[m][k]), 0, MPFR_RNDN);
      }
    }

  }

  // relative prune of real and imaginary part
  for (int m=0; m<dim; m++) {
    mpc_rk1_prune_re_im(sol_eps_ord[m], precision*log2(10), eps_num);
  }

  // check
  mpc_t mpc_check;
  mpc_init3(mpc_check, wp2, wp2);
  mpc_t mpc_diff;
  mpc_init3(mpc_diff, wp2, wp2);
  for (int m=0; m<dim; m++) {
    for (int ep=0; ep<eps_num; ep++) {
      // cout << "m, ep = " << m << ", " << ep << endl;
      // reconstruct series at selected point
      mpc_set(mpc_check, sol_eps_ord[m][eps_num-1], MPFR_RNDN);
      for (int epp=eps_num-1; epp>0; epp--) {
        mpc_fma(mpc_check, mpc_eps_list[ep], mpc_check, sol_eps_ord[m][epp-1], MPFR_RNDN);
      }
      mpc_pow_si(mpc_diff, mpc_eps_list[ep], -2*loops+starting_ord[m], MPFR_RNDN);
      mpc_mul(mpc_check, mpc_check, mpc_diff, MPFR_RNDN);
      // compare with pre-interpolation value
      mpc_sub(mpc_diff, mpc_check, sol_at_eps[m][ep], MPFR_RNDN);
      mpc_div(mpc_diff, mpc_diff, sol_at_eps[m][ep], MPFR_RNDN);
      // if (1) {
      if (mpc_lessthan_tol(mpc_diff)) {
        cout << "found difference at m, ep = " << m << ", " << ep << endl;
        cout << "pre-interpolation:" << endl;
        print_mpc(&sol_at_eps[m][ep]); cout << endl;
        cout << "post-interpolation:" << endl;
        print_mpc(&mpc_check); cout << endl;
        cout << "diff:" << endl;
        print_mpc(&mpc_diff); cout << endl;
      }      
    }
  }

  // FREE
  mpc_rk1_clear(mpc_eps_list, eps_num);
  delete[] mpc_eps_list;
  mpc_rk2_clear(eps_mat, eps_num, eps_num + null_pow);
  del_rk2_tens(eps_mat, eps_num);
  mpc_rk2_clear(eps_mat_copy, eps_num, eps_num);
  del_rk2_tens(eps_mat_copy, eps_num);
  mpc_rk3_clear(inv_eps_mat, null_pow, eps_num, eps_num);
  del_rk3_tens(inv_eps_mat, null_pow, eps_num);
  mpc_rk3_clear(inv_eps_mat_red, null_pow, eps_num-1, eps_num-1);
  del_rk3_tens(inv_eps_mat_red, null_pow, eps_num-1);
  
}

