#ifndef KIRA_H
#define KIRA_H

void kira_to_DE_pf(
  // OUTPUT
  poly_frac ***pfmat,
  int *zero_label, int *nroots, mpc_t **roots, mpfr_t **tols,
  // INPUT
  int nroots_branch, mpc_t *roots_branch, mpfr_t *tols_branch,
  int ninvs, char **symbols, int *is_mass,
  poly_frac *pspf,
  int *skip_inv, char ***ep_kin,
  int dim_eta, char ****mats_str,
  int nbranches, int *branch_deg,
  int eps_num, char **eps_str,
  LI* MI_eta, char *dir_amflow,
  FILE *terminal
);

void call_kira(
  // OUTPUT
  LI **MI_eta, int *dim_eta, int **MI_idx, int *dim,
  // INPUT
  int redo, int opt_kira_parallel,
  char **kin, char *dir_parent, char *dir_amflow,
  FILE *terminal
);

#endif