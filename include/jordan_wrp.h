#ifndef JORDAN_WRP_H
#define JORDAN_WRP_H

void jordan_triv_1(
  // OUTPUT
  mpc_t **jor, mpc_t **t_mat, int **ein_pos, int *num_ein,
  // INPUT
  mpc_t *ac, int N, mpc_t *prop, mpc_t *eig, int first_non_zero,
  int print
);

void jordan_triv_2(
  // OUTPUT
  mpc_t **jor, mpc_t **t_mat, int **ein_pos, int *num_ein,
  // INPUT
  mpc_t *ac, int N, mpc_t *prop, int first_non_zero, int last_non_zero_row,
  int *is_zero,
  int print
);

void jordan_wrp(
  // OUTPUT
  mpc_t **jor, mpc_t **t_mat, int **ein_pos, int *num_ein,
  // INPUT
  mpfr_t *ar, mpfr_t *ai, int N, int nm, int debug
);

#endif