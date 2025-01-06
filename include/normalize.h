#ifndef NORMALIZE_H
#define NORMALIZE_H

void RLeeAlg1(
  // OUTPUT
  int *k0, int *dimS, int *&S, mpc_t **Delta,
  // INPUT
  mpc_t **L0, int num_ein, int r, int *chain_lengths
);

void find_projector(
  // OUTPUT
  mpc_t **proj,
  // INPUT
  mpfr_t *A0r, mpfr_t *A0i, mpc_t **A1, int dim
);

void find_eigen_duplicates(
  // IN-OUT
  int* eig_grid,
  //OUTPUT
  mpc_t *eigenvalues, int dim
);

void decompose_eigenvalues(
  // OUTPUT
  int *eig_int_diff,
  // INPUT
  int dim, mpc_t *eigenvalues, int *eig_grid,
  int *eq_class, int *class_min
);

void group_eigenvalues(
  // OUTPUT
  int *&out_class_min, int *num_classes,
  // INPUT
  int dim, mpc_t *eigenvalues, int *eig_grid,
  int *eq_class, int *eig_int_diff
);

#endif