#ifndef TENSOR_UTILS_H
#define TENSOR_UTILS_H

//////
// PRINT
//////
void print_rk1_mpc(mpc_t *tens, int dim);

void print_rk2_mpc(mpc_t **tens, int dim1, int dim2);

void print_rk3_mpc(mpc_t ***tens, int dim1, int dim2, int dim3);

//////
// INITIALIZE WITH MULTIPLE PRECISION
//////
void init_rk1_mpfr(mpfr_t *tens, int dim1);

void init_rk1_mpc(mpc_t *tens, int dim1);

void init_rk1_mpc_wp2(int wp2, mpc_t *tens, int dim1);

void init_rk2_mpfr(mpfr_t **tens, int dim1, int dim2);

void init_rk2_mpc(mpc_t **tens, int dim1, int dim2);

void init_rk3_mpc(mpc_t ***tens, int dim1, int dim2, int dim3);

void init_rk4_mpc(mpc_t ****tens, int dim1, int dim2, int dim3, int dim4);

//////
// CLEAR
//////
void mpz_rk1_clear(
  mpz_t *tens, int dim
);

void mpq_rk1_clear(
  mpq_t *tens, int dim
);

void mpfr_rk1_clear(
  mpfr_t *tens, int dim
);

void mpc_rk1_clear(
  mpc_t *tens, int dim
);

void mpc_rk2_clear(
  mpc_t **tens, int dim1, int dim2
);

void mpc_rk3_clear(
  mpc_t ***tens, int dim1, int dim2, int dim3
);

void mpc_rk4_clear(
  mpc_t ****tens, int dim1, int dim2, int dim3, int dim4
);

//////
// COPY
//////
void copy_rk1_int(int *copy, int *tens, int dim1);

void copy_rk1_mpc(mpc_t *copy, mpc_t *tens, int dim1);

void copy_rk1_mpfr(mpfr_t *copy, mpfr_t *tens, int dim1);

// void copy_rk2_tens(ex **tens_copy, ex **tens, int dim1, int dim2);

void copy_rk2_mpc(mpc_t **copy, mpc_t **tens, int dim1, int dim2);

//////
// COMPARE
//////
void int_rk1_compare(int *tens1, int *tens2, int dim);

void int_rk1_compare_perm(
  int *perm, int *tens1, int *tens2, int dim
);

void mpc_rk1_compare(mpc_t *tens1, mpc_t *tens2, int dim);

void mpc_rk1_compare_double(mpc_t *tens1, mpc_t *tens2, int dim);

void mpc_rk1_compare_perm(
	// OUTPUT
	int *perm,
	// INPUT
	int dim, mpc_t *tens1, mpc_t *tens2
);

int mpc_rk1_prop(
  // OUTPUT
  mpc_t *cf,
  // INPUT
  mpc_t *v1, mpc_t *v2, int dim,
  int *is_zero_v1, int *is_zero_v2
);

void mpc_rk2_compare(mpc_t **tens1, mpc_t **tens2, int dim1, int dim2);

void poly_frac_rk2_compare(
  struct poly_frac **tens1,
  struct poly_frac **tens2,
  int dim1, int dim2,
  mpc_t *roots
);

//////
// COMPARE
//////
void int_rk1_append_zero(
  // IN-OUT
  int *dim, int **vec,
  // INPUT
  int ap_dim
);

void int_rk1_append(
  // IN-OUT
  int *dim, int **vec,
  // INPUT
  int ap_dim, int *ap_vec
);

void mpc_rk1_append(
  // IN-OUT
  int *dim, mpc_t **vec,
  // INPUT
  int ap_dim, mpc_t *ap_vec
);

void mpfr_rk1_append(
  // IN-OUT
  int *dim, mpfr_t **vec,
  // INPUT
  int ap_dim, mpfr_t *ap_vec
);

//////
// TRANSFORM
//////
void mpfr_rk2f_switch_rows(
  // IN-OUT
  mpfr_t *tens,
  // INPUT
  int idx1, int idx2, int dim
);

void mpfr_rk2f_switch_cols(
  // IN-OUT
  mpfr_t *tens,
  // INPUT
  int idx1, int idx2, int dim1, int dim2
);

void mpc_rk2_switch_rows(
  // IN-OUT
  mpc_t **tens,
  int idx1, int idx2, int dim
);

void mpc_rk2_switch_cols(
  // IN-OUT
  mpc_t **tens,
  int idx1, int idx2, int dim
);

//////
// TENSOR ALGEBRA
//////
void mpc_rk2_mul_mpc_rk1(
  // OUTPUT
  mpc_t *out_tens,
  // INPUT
  mpc_t **tens1, mpc_t *tens2,
  int dim1, int dim2
);

void mpc_rk2_mul_mpc_rk2(
  // OUTPUT
  mpc_t **out_tens,
  // INPUT
  mpc_t **tens1, mpc_t **tens2,
  int dim1, int dim, int dim2
);

void mpc_rk2_mul_mpc_rk2_slice(
  // OUTPUT
  mpc_t *out_tens,
  // INPUT
  mpc_t **tens1, mpc_t **tens2,
  int dim1, int dim2, int slice
);

//////
// OPERATIONS
//////
void int_rk1_to_str(
  // OUTPUT
  char **out,
  // INPUT
  int dim, int *vec
);

int int_rk1_count_postivie(
  // INPUT
  int *vec, int dim
);

int int_rk1_sum_postivie(
  // INPUT
  int *vec, int dim
);

int int_rk1_sum_negative_abs(
  // INPUT
  int *vec, int dim
);

int int_rk1_is_positive_to_decimal(
  int *vec, int dim
);

int int_rk1_is_zero_to_decimal(
  int *vec, int dim
);

#endif