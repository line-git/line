#ifndef SYSTEM_ANALYZER_H
#define SYSTEM_ANALYZER_H

void generate_sub_diag_grid(
  // OUTPUT
  int ***subblock,
  // INPUT
  int **profile, int nblocks, struct poly_frac **pfmat
);

void poly_frac_rk2_find_profile(
  // OUTPUT
  int ***profile, int ***subblock, int *nblocks,
  // INPUT
  struct poly_frac **pfmat, int g_dim
);

void update_new_roots(
  // OUTPUT
  int *label, int *num_out_roots, mpc_t *out_roots, mpfr_t *out_tols,
  // INPUT
  int *num_in_roots, mpc_t *in_roots, mpfr_t *in_tols,
  mpc_t *root, mpfr_t *tol
);

void linearize_masses_nd(
	// IN-OUT
	struct lnode *nd,
	// INPUT
	int *is_mass
);

void process_branch_points(
  // OUTPUT
  int *num_branch_roots, mpc_t *branch_roots, mpfr_t *branch_tols,
  int **branch_sing_lab,
  int *branch_deg, mpc_t **branch_poly, 
  // IN-OUT
  int *nroots, mpc_t **roots, mpfr_t **tols,
  // INPUT
  int nbranches, char **branch_cut_str,
  int ninvs, char **symbols, int *is_mass,
  int *skip_inv, char ***ep_kin, struct poly_frac* pspf
);

void linearize_mass_DE(
  // IN-OUT
  struct lnode *nd,
  // INPUT
  int sym_idx
);

void wrt_cmp_DE(
  // OUTPUT
  int **perm,
  // INPUT
  int *zero_label, int *nroots, mpc_t **roots,
  poly_frac ***pfmat, int eps_num, int dim,
  double wp2_rel_decr,
  char *file_ext, char *filepath_matrix, char *filepath_roots,
  int opt_write
);

void generate_poly_frac_DE(
  // OUTPUT
  poly_frac ***pfmat,
  int *zero_label, int *nroots, mpc_t **roots, mpfr_t **tols,
  // INPUT
  int nroots_branch, mpc_t *roots_branch, mpfr_t *tols_branch,
  int ninvs, char **symbols, int *is_mass,
  poly_frac *pspf,
  int *skip_inv, char ***ep_kin,
  int dim, char ****mats_str,
  int nbranches, int *branch_deg,
  int eps_num, char **eps_str,
  FILE *terminal
);

void build_block_indices(int *&b_idx, int *&sb_idx, int *b_len, int *sb_len, int b, int** prof, int **sb_grid);

template <typename type>
void select_block(type **&block, int b, int sb, type **mat, int** prof);

void generate_constant_term(
  mpc_t **constant_term,
  mpc_t ***sblock_num, mpc_t ***sblock_den, mpc_t **solutions,
  int b_len, int sb_len, int *sb_idx,
  int num_max_deg, int den_max_deg, int eta_ord,
  int log_sb_idx_dim = 0, int *log_sb_idx = NULL);

void sys_block_info(
  // OUTPUT
  mpc_t *&lcm, mpc_t ***&block, mpc_t ****&const_term,
  // IN-OUT
  bool *flag_lcm, int *lcm_deg, mpc_t *&lcm_orig, int *&lcm_roots,
  int *block_max_deg, int *num_max_deg, int *den_max_deg,
  // INPUT
  int b_len, int sb_len,
  int *b_idx, int *sb_idx,
  struct poly_frac **lcm_mat_ep, mpc_t ****mpc_lcm_mat_ep,
  int npoles, mpc_t *poles,
  mpc_t eta_shift,
  int eta_ord,
  mpc_t ****solutions, int num_prev_eig=-1, int *prev_eig=NULL,
  int *prev_log_len=NULL, int *sol_log_len=NULL, int **log_prof=NULL
);

#endif