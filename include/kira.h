#ifndef KIRA_H
#define KIRA_H

typedef struct lnode_pf {
  int sign;
  int info;  // 1: pf is one, -1: pf is zero, 2: other
  int den_is_one;  // 1: den is one
  int num_nterms;  // number of terms in the numerator
  int *num_pows;  // power of terms in the numerator
  lnode *num_terms;  // term trees in the numerator
  int den_nterms;  // number of terms in the denominator
  int *den_pows;  // power of terms in the denominator
  lnode *den_terms;  // trees in the denominator
} lnode_pf;

void lnode_pf_build(lnode_pf *ndpf);

void lnode_pf_free(lnode_pf *ndpf);

void process_kira_IBPs(
  // OUTPUT
  mpz_t ******coeffs_num_den, int ******pows_num_den,
  int *****nterms_num_den, lnode_pf ***der_ndpf,
  char ****kira_str,
  // INPUT
  int ninvs, char **symbols, int *is_mass,
  int dim_eta, int max_num_contr,
  LI* MI_eta, char *dir_amflow,
  FILE *terminal
);

void kira_IBPs_to_DE_pf(
  // OUTPUT
  poly_frac ***pfmat,
  int *zero_label, int *nroots, mpc_t **roots, mpfr_t **tols,
  // INPUT
  int ep, int nroots_branch, mpc_t *roots_branch, mpfr_t *tols_branch,
  int ninvs, char **symbols, int *is_mass,
  poly_frac *pspf,
  int *skip_inv, char ***ep_kin,
  int dim_eta, char ****mats_str,
  int nbranches, int *branch_deg,
  char **eps_str,
  LI* MI_eta, char *dir_amflow,
  int max_num_contr,
  mpz_t ******coeffs_num_den, int ******pows_num_den,
  int *****nterms_num_den, lnode_pf ***der_ndpf,
  char ****kira_str,
  FILE *terminal
);

// void kira_to_DE_pf(
//   // OUTPUT
//   poly_frac ***pfmat,
//   int *zero_label, int *nroots, mpc_t **roots, mpfr_t **tols,
//   // INPUT
//   int nroots_branch, mpc_t *roots_branch, mpfr_t *tols_branch,
//   int ninvs, char **symbols, int *is_mass,
//   poly_frac *pspf,
//   int *skip_inv, char ***ep_kin,
//   int dim_eta, char ****mats_str,
//   int nbranches, int *branch_deg,
//   int eps_num, char **eps_str,
//   LI* MI_eta, char *dir_amflow,
//   FILE *terminal
// );

void call_kira(
  // OUTPUT
  LI **MI_eta, int *dim_eta, int **MI_idx, int *dim,
  // INPUT
  int redo, int opt_kira_parallel,
  char **kin, char *dir_parent, char *dir_amflow,
  FILE *terminal
);

#endif