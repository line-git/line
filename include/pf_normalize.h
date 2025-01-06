#ifndef PF_NORMALIZE_H
#define PF_NORMALIZE_H

void pf_to_Fuchsian(
  // OUTPUT
  struct poly_frac **bal,
  struct poly_frac **inv_bal,
  // IN-OUT
  struct poly_frac **matb,
  // INPUT
  int b_len,
  mpc_t *roots, int nroots
);

void print_eigenvalues(
  int dim, int num_classes, int *eig_grid, int *eq_class,
  mpc_t *eigenvalues
);

void pf_NormalizeDiagonal(
  // OUTPUT
  struct poly_frac **tmat, struct poly_frac **inv_tmat,
  mpc_t **eig_list, int *num_classes, int *eq_class, int *eig_grid,
  // IN-OUT
  struct poly_frac **mat_ep,
  // INPUT
  int dim, int nblocks, int **prof, int **sb_grid,
  mpc_t *roots, int nroots,
  FILE *terminal
);

void pf_apply_sub_diagonal_p_reduction(
  // IN-OUT
  struct poly_frac **mat,
  // INPUT
  int p_rank, int b, int sb, int nblocks, int **prof,
  mpc_t *sb_sys_sol, mpc_t *roots, int nroots
);

void pf_to_Fuchsian_global(
  // IN-OUT
  struct poly_frac **mat_ep,
  struct poly_frac **tmat, struct poly_frac **inv_tmat,
  // INPUT
  int dim, int nblocks, int **prof, int **sb_grid,
  mpc_t *roots, int nroots,
  FILE *terminal
);

void pf_NormalizeMat(
  // OUTPUT
  struct poly_frac **tmat, struct poly_frac **inv_tmat,
  int *num_classes, int *eq_class, mpc_t **eig_list, int *eig_grid, 
  // IN-OUT
  struct poly_frac **mat_ep,
  // INPUT
  int dim, int nblocks, int **prof, int **sb_grid,
  mpc_t *roots, int nroots,
  FILE *terminal
);

void pf_check_NormalizeMat(
  struct poly_frac **tmat, struct poly_frac **inv_tmat,
  int num_classes, int *eq_class, mpc_t *eig_list, int *eig_grid,
  struct poly_frac **mat_ep,
  int dim, int nblocks, int **prof, int **sb_grid,
  mpc_t *roots, int nroots,
  FILE *terminal
);

#endif