#ifndef EQ_SOLVER_H
#define EQ_SOLVER_H

void solve_block(mpc_t **solutions,
  mpc_t *lcm, mpc_t ***block, mpc_t **const_term,
  int *b_idx, int b_len, int sb_len,
  int eta_ord, int lcm_deg, int block_max_deg
);

void zero_solve_block(mpc_t *****solutions,
  mpc_t *lcm, mpc_t ***block,
  int *b_idx, int b_len, int sb_len,
  int eta_ord, int lcm_deg, int block_max_deg,
  mpc_t *eig_list, int *cum_eig, int num_cum_eig,
  int *block_log_len, int *sol_log_len, int num_classes,
  bool particular=false, mpc_t ****const_term = NULL, int *rhs_log_len = NULL
  // ,int match_in_zero = 0, int lam0 = 0, mpc_t **boundary = NULL
);

void match_next_eta(
  mpc_t *eta_values, int et,
  mpc_t **solutions, int dim, int eta_ord
);

void propagate_regular(
  // IN-OUT
  mpc_t ***solutions,
  // INPUT
  int neta_values, mpc_t *eta_values,
  int dim, struct poly_frac **mat_ep,
  int npoles, mpc_t *poles,
  int nblocks, int **prof, int **sb_grid, int eta_ord
);

void get_block_log_len(
  // OUTPUT
  int **block_log_prof, int *block_log_len, int *rhs_log_len,
  // INPUT
  int b_len, int num_classes, int offset, int *eig_grid, int *eq_class
);

int count_log_len(mpfr_t mpfr_tol, mpc_t ***tens, int dim1, int dim2, int dim3, int *sb_idx);

void evaluate_series_around_zero(
  mpc_t ***output,
  mpc_t *mpc_eta,
  mpc_t *****solutions,
  // int block_max_deg,
  mpc_t *eig_list, int num_eig, int *eig_labels, int *log_prof, int dim1, int dim2, int eta_ord, int num_classes, int lam0,
  bool rk4 = false, mpc_t ****rk4_sol = NULL,
  int *match_constr_dim = NULL, int **match_constr_idx = NULL, bool particular = true,
  int analytic_continuation = 0
);

void solve_log_constraints(
  // IN-OUT
  mpc_t *****psol,
  // INPUT
  int eta_ord, int b_len, int offset, int num_cum_eig, int *cum_eig, int *eq_class,
  int *block_log_len, int *sol_log_len, int *rhs_log_len, mpc_t *****hsol
);

void match_sol_around_zero(
  // OUTPUT
  mpc_t ****sol,
  // INPUT
  int eta_ord, int b_len, int sb_len, int offset, int num_cum_eig, int *cum_eig,
  int num_classes, int *eq_class, mpc_t *eig_list,
  int *block_log_len, int *sol_log_len, int *rhs_log_len,
  mpc_t *****hsol, mpc_t *****psol,
  mpc_t ***solutions, mpc_t *mpc_eta, int lam0,
  int analytic_cont,
  mpc_t ****prev_sol, struct poly_frac **tmat, int nroots, mpc_t *roots,
  int **bound_behav, int **mi_eig, int *mi_eig_num, int **log_prof,
  struct poly_frac **pfmat
);

void pf_limit_in_zero(
  // OUTPUT
  mpc_t *final_sol,
  // INTPUT
  mpc_t ****sol, struct poly_frac **tmat,
  int dim, int lam0, int tp_rank,
  mpc_t *roots, int nroots,
  int rad_exp
);

void pf_limit_in_zero_block(
  // OUTPUT
  mpc_t *final_sol,
  // IN-OUT
  struct poly_frac *pf_sol,
  // INTPUT
  mpc_t ****sol, struct poly_frac **tmat,
  int offset, int b_len, int lam0, int tp_rank,
  mpc_t *roots, int nroots,
  int rad_exp
);

void solve_zero(
  // IN-OUTPUT
  mpc_t ***solutions,
  // INPUT
  int dim, struct poly_frac **mat_ep,
  struct poly_frac **tmat, struct poly_frac **inv_tmat,
  int npoles, mpc_t *poles,
  mpc_t *matching_point,
  int nblocks, int **prof, int **sb_grid, int eta_ord,
  int num_classes, mpc_t *eig_list, int *eq_class, int *eig_grid,
  mpc_t *roots, int nroots,
  int cross, mpc_t *eta_target, int analytic_cont,
  int *is_mass, int *skip_inv, int ninvs, mpc_t *PS_ini, mpc_t *PS_fin, char *eps_str,
  int try_analytic,
  int **bound_behav, int **mi_eig, int *mi_eig_num,
  FILE *logfptr, FILE *terminal
);

void center_around_pole(
  // OUTPUT
  struct poly_frac **sh_pfmat, mpc_t *sh_roots,
  // INPUT
  struct poly_frac **pfmat, mpc_t *roots, int label,
  int nroots, int dim1, int dim2
);

void set_analytic_cont_sign(
  // OUTPUT
  int *analytic_cont,
  // INPUT
  mpc_t *m1, mpc_t *m2,
  int branch_deg, mpc_t *branch_poly
);

void propagate_along_path(
  // OUTPUT
  mpc_t **sol_at_eps,
  // IN-OUT
  mpc_t ***solutions,
  // INPUT
  int dim, int eta_ord,
  int ep, char *eps_str,
  int neta_values, mpc_t *path, int *path_tags, int nsings,
  struct poly_frac **pfmat, int nroots, mpc_t *roots, int *sing_lab,
  int nblocks, int **prof, int **sb_grid,
  int nbranches, int *branch_deg, mpc_t **branch_poly, int **branch_sing_lab,
  int ninvs, mpc_t *PS_ini, mpc_t *PS_fin, char **symbols,
  int *is_mass, int *skip_inv,
  int **bound_behav, int **mi_eig, int *mi_eig_num,
  FILE *logfptr, FILE *terminal
);

void propagate_eps(
  // OUTPUT
  mpc_t **sol_at_eps,
  // INPUT
  int ep, int exit_sing, int nloops,
  poly_frac ***pfmat,
  int *zero_label, int *nroots, mpc_t **roots,
  int *neta_values, mpc_t **path, int **path_tags, int *nsings, int **sing_lab,
  mpc_t ***solutions, int dim, int eta_ord,
  int eps_num, char **eps_str,
  int ninvs, mpc_t *PS_ini, mpc_t *PS_fin, char **symbols,
  int *is_mass, int *skip_inv,
  int nbranches, int *branch_deg, mpc_t **branch_poly, int **branch_sing_lab,
  // int **bound_behav, int **mi_eig, int *mi_eig_num,
  int gen_bound, char *filepath_bound, mpc_t **bound,
  char *filepath_bound_build, char *filepath_bound_behav,
  // char *filepath_matrix, char *filepath_roots, // char *filepath_branch_sing_lab,
  char *filepath_path, char *filepath_path_tags, char *filepath_sol, char* dir_partial,
  char *file_ext, FILE *logfptr, int opt_write, int opt_checkpoint,
  FILE *terminal
);

#endif