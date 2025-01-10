#ifndef UTILS_H
#define UTILS_H

int MIN(int a, int b);

inline int MAX(int a, int b);

int int_pow_int(int base, int exp);

int count_lines(char *filename);

int count_lines_fptr(FILE *fptr);

void join_path(
  // OUTPUT
  char **path,
  // INPUT
  char *parent, char *name
);

int mpfr_lessthan_tol(mpfr_t in);

int mpfr_sign_within_tol(mpfr_t in);

int mpc_lessthan_tol(mpc_t in);

int mpfr_equal_within_tol(
  mpfr_t in1, mpfr_t in2
);

int mpc_equal_within_tol(
  mpc_t in1, mpc_t in2
);

int mpc_zero_p(
  mpc_t in
);

void mpc_rk2_prune_abs_tol(
  mpc_t **tens,
  int dim1, int dim2
);

int max_ints(int *arr, int dim);

void minmax0(double *min, double *max, double *arr, int MAX_SIZE);

double min1(double *arr, int MAX_SIZE);

int length_file(ifstream &fmat);

void rationalize(mpq_t out_rat, double x, double tolerance);

void limit_denominator(
  // OUTPUT
  mpq_t out_rat,
  // INPUT
  mpq_t rat, mpz_t max_denominator
);

void rationalize_cf_from_str(mpq_t rat, char *str);

void print_mpc(mpc_t *mpc_num);

bool mpfr_is_diff_int(mpfr_t real1, mpfr_t real2);

void mpc_swap(mpc_t *z1, mpc_t *z2);

void mpfr_tol_enlarge(double wp2_rel_decr);

void mpfr_tol_enlarge_rel(double wp2_rel_decr);

int mpfr_log2_int(mpfr_t log2_mpfr);

size_t mpfr_get_memory_usage(mpfr_t x);

void interpolate_epsilon_orders(
  // OUTPUT
  mpc_t **sol_eps_ord,
  // INPUT
  mpc_t **sol_at_eps, char **eps_str,
  int dim, int eps_num, int loops,
  int precision,
  int *starting_ord
);

void interpolate_epsilon_orders_prune(
  // OUTPUT
  mpc_t **sol_eps_ord,
  // INPUT
  mpc_t **sol_at_eps, char **eps_str,
  int dim, int eps_num, int loops,
  int precision, int order,
  int *starting_ord
);

#endif
