#ifndef POLY_FRAC_H
#define POLY_FRAC_H

void poly_eval_sym_zero_by_roots_wo_0(
  // OUTPUT
  mpc_t *value,
  // INPUT,
  mpc_t *roots, int *mults, int nroots
);

void poly_coeff1_by_roots_wo_0(
  // OUPUT
  mpc_t *value,
  // INPUT,
  mpc_t *roots, int *mults, int nroots
);

void poly_eval_by_roots(
  // OUTPUT
  mpc_t *out,
  // INPUT
  mpc_t *roots, int *mults, int nroots, mpc_t *val
);

void poly_eval(
  // OUTPUT
  mpc_t *out,
  // INPUT
  mpc_t *pol, int deg, mpc_t *val
);

void mul_vpoly(
  mpc_t *out_pol, mpc_t *pol1, mpc_t *pol2,
  int deg1, int deg2,
  int vdeg1, int vdeg2,
  int ldeg1, int ldeg2
);

void poly_mul_root_new_out(
  mpc_t *out, mpc_t *coeffs, int deg, mpc_t root
);

typedef struct poly_frac {
  int num_deg;
  int num_vdeg;
  int den_deg;
  int nroots;
  int *mults;
  mpc_t *coeffs;
} poly_frac;

void poly_frac_build(
  struct poly_frac *pf
);

void poly_frac_rk1_build(
  struct poly_frac *pf_tens,
  int dim1
);

void poly_frac_rk2_build(
  struct poly_frac **pf_tens,
  int dim1, int dim2
);

void poly_frac_rk3_build(
  struct poly_frac ***pf_tens,
  int dim1, int dim2, int dim3
);

void poly_frac_free(
  struct poly_frac *pf
);

void poly_frac_rk1_free(
  struct poly_frac *pf_tens, int dim1
);

void poly_frac_rk2_free(
  struct poly_frac **pf_tens, int dim1, int dim2
);

void poly_frac_rk3_free(
  struct poly_frac ***pf_tens, int dim1, int dim2, int dim3
);

void poly_frac_alloc_init_coeffs(
  struct poly_frac *pf
);

void del_rk1_poly_frac(
  struct poly_frac *pftens,
  int dim  
);

void del_rk2_poly_frac(
  struct poly_frac **pftens,
  int dim1, int dim2
);

void poly_frac_print(
  struct poly_frac *pf
);

void poly_frac_rk2_print(
  struct poly_frac **pf_tens,
  int dim1, int dim2
);

void poly_frac_print_to_math(
  struct poly_frac *pf,
  mpc_t *roots
);

void poly_frac_rk2_print_to_math(
  struct poly_frac **pf_tens,
  int dim1, int dim2,
  mpc_t *roots
);

void poly_frac_set_info(
  struct poly_frac *pf,
  int num_deg, int num_vdeg, int den_deg
);

void poly_frac_set_info_pf(
  struct poly_frac *pfout,
  struct poly_frac *pfin
);

void poly_frac_set_zero(
    struct poly_frac *pf
);

void poly_frac_set_ui(
  struct poly_frac *pf,
  int ui_num,
  int nroots
);

void poly_frac_set_si(
  struct poly_frac *pf,
  int si_num,
  int nroots
);

void poly_frac_set_ui_wp2(
  int wp2,
  struct poly_frac *pf,
  int ui_num,
  int nroots
);

void poly_frac_set_mpc(
  // OUTPUT
  struct poly_frac *pf,
  // INPUT
  mpc_t *mpcin,
  int nroots
);

void poly_frac_rk2_set_mpc_rk2(
  struct poly_frac **pf_tens,
  mpc_t **mpc_tens,
  int dim1, int dim2,
  int nroots
);

void poly_frac_set_pf(
  struct poly_frac *pfout,
  struct poly_frac *pfin
);

void poly_frac_rk2_set_pf_rk2(
  struct poly_frac **pf_tens_out,
  struct poly_frac **pf_tens_in,
  int dim1, int dim2
);

void poly_frac_set_coeffs(
  struct poly_frac *pf,
  mpc_t *coeffs, int deg
);

void poly_frac_set_coeffs_wp2(
  struct poly_frac *pf,
  mpc_t *coeffs, int deg,
  int wp2
);

void poly_frac_rk2_set_id(
  struct poly_frac **pf_tens,
  int dim,
  int nroots
);

int poly_frac_pow_behav(
  struct poly_frac *pf
);

void poly_frac_normal(
  struct poly_frac *pf
);

void poly_frac_prune(
  struct poly_frac *pf
);

void poly_frac_prune_tol(
  struct poly_frac *pf, double wp2_rel_decr
);

int mpc_get_exp2(
  mpc_t in
);

void poly_frac_prune_rel_tol(
  struct poly_frac *pf, int wp_bin
);

void mpc_rk1_prune_rel_tol_real_max(
  mpc_t *coeffs, int dim, int wp_bin, mpc_t *val
);

void poly_frac_rk2_prune_rel_tol_real_max(
  struct poly_frac **pf, int wp_bin,
  int dim1, int dim2
);

void poly_frac_rk2_trim_zero_p(
  struct poly_frac **pf,
  int dim1, int dim2
);

void poly_frac_rk1_prune_rel_tol(
  struct poly_frac *pf, int wp_bin,
  int dim1
);

void poly_frac_rk2_normal(
  struct poly_frac **pf,
  int dim1, int dim2
);

void poly_frac_prune_radius(
  struct poly_frac *pf, int wp_bin, int rad_exp
);

void mpc_rk1_prune_radius(
  mpc_t *coeffs, int dim, int wp_bin, mpc_t *val
);

void poly_frac_rk1_prune_radius(
  struct poly_frac *pf, int wp_bin, int rad_exp,
  int dim1
);

void poly_frac_rk2_prune_radius(
  struct poly_frac **pf, int wp_bin, int wp_rad_exp,
  int dim1, int dim2
);

void poly_frac_rk2_prune(
  struct poly_frac **pf,
  int dim1, int dim2, double wp2_rel_decr
);

void poly_frac_rk2_prune_rel_tol(
  struct poly_frac **pf, int wp_bin,
  int dim1, int dim2
);

void poly_frac_shift(
  struct poly_frac *pfout,
  struct poly_frac *pfin,
  mpc_t *shift,
  int swap_root_lab
);

void poly_frac_rk2_shift(
  struct poly_frac **pfout,
  struct poly_frac **pfin,
  mpc_t *shift,
  int swap_root_lab,
  int dim1, int dim2
);

void poly_frac_perm_roots(
  // IN-OUT
  struct poly_frac *pf,
  // INPUT
  int *perm
);

void poly_frac_rk2_perm_roots(
  // IN-OUT
  struct poly_frac **pf,
  // INPUT
  int *perm, int dim1, int dim2
);

void poly_frac_rk2_infty(
  // OUTPUT
  struct poly_frac **pfout, mpc_t *out_roots, int *out_zero_label,
  // INPUT
  struct poly_frac **pfin, mpc_t *in_roots,
  int dim1, int dim2, int nroots
);

void poly_frac_neg(
  struct poly_frac *pf
);

void poly_frac_rk2_neg(
  struct poly_frac **pfmat, int dim1, int dim2
);

void poly_frac_mul_ui(
  struct poly_frac *pfout,
  struct poly_frac *pfin,
  int ui
);

void poly_frac_mul_si(
  struct poly_frac *pfout,
  struct poly_frac *pfin,
  int si
);

void poly_frac_mul_mpc(
  struct poly_frac *pfout,
  struct poly_frac *pfin,
  mpc_t *mpcin
);

void poly_frac_div_mpc(
  struct poly_frac *pfout,
  struct poly_frac *pfin,
  mpc_t *mpcin
);

void poly_frac_mul_pf(
  struct poly_frac *pfout,
  struct poly_frac *pfin1,
  struct poly_frac *pfin2
);

void poly_frac_mul_sym_pow(
  struct poly_frac *pfout,
  struct poly_frac *pfin,
  int pow
);

void poly_frac_rk2_mul_sym_pow(
  struct poly_frac **pfmat, int pow,
  int dim1, int dim2
);

void poly_frac_mul_root(
  struct poly_frac *pf,
  mpc_t *root, int mult
);

void poly_frac_rk2_mul_root(
  struct poly_frac **pfmat,
  mpc_t *root, int mult,
  int dim1, int dim2
);

void poly_frac_den_to_poly(
  // OUTPUT
  mpc_t *coeffs,
  // INPUT
  struct poly_frac *pfin, mpc_t *roots
);

void poly_frac_den_to_num(
  // OUPUT
  struct poly_frac *pfout,
  // INPUT
  struct poly_frac *pfin,
  mpc_t *roots,
  int wp2
);

void poly_frac_pow_ui(
  struct poly_frac *pfout,
  struct poly_frac *pfin,
  int powui
);

void poly_frac_pow_si(
  // OUTPUT
  struct poly_frac *pfout,
  mpc_t **roots, mpfr_t **tols,
  // INPUT
  struct poly_frac *pfin,
  int pow,
  mpc_t *in_roots
);

void poly_frac_add_pf(
  struct poly_frac *pfout,
  struct poly_frac *pfin1,
  struct poly_frac *pfin2,
  mpc_t *roots,
  mpc_t *flop_mpc
);

void poly_frac_rk2_add_pf_rk2(
  struct poly_frac **pf_tens_out,
  struct poly_frac **pf_tens_in1,
  struct poly_frac **pf_tens_in2,
  mpc_t *roots, mpc_t *flop_mpc,
  int dim1, int dim2
);

int poly_frac_eval_sym_zero(
  // OUTPUT
  mpc_t *value,
  // INPUT
  struct poly_frac *pf, mpc_t *roots
);


void poly_frac_rk1_eval_sym_zero(
  // OUTPUT
  mpc_t *mpc_tens,
  // INPUT
  struct poly_frac *pf_tens, mpc_t *roots, int dim
);


void poly_frac_eval_value(
  // OUTPUT
  mpc_t *out,
  // INPUT
  struct poly_frac *pf, mpc_t *roots, mpc_t *val
);

void poly_frac_rk2_eval_value(
  // OUTPUT
  mpc_t **mpc_tens,
  // INPUT
  struct poly_frac **pf_tens, mpc_t *roots, mpc_t *val,
  int dim1, int dim2
);

void poly_frac_extract_LO(
  // OUTPUT
  mpfr_t *out_re, mpfr_t *out_im,
  // INPUT
  struct poly_frac *pf, mpc_t *roots, int pow
);

void poly_frac_mat_extract_LO(
  // OUTPUT
  mpfr_t *LOre, mpfr_t *LOim,
  // INPUT
  struct poly_frac  **pfmat, mpc_t *roots, int pow,
  int dim1, int dim2
);

void poly_frac_mat_extract_LOc(
  // OUTPUT
  mpc_t **LO,
  // INPUT
  struct poly_frac  **pfmat, mpc_t *roots, int pow,
  int dim1, int dim2
);

void poly_frac_extract_NLO(
  // OUTPUT
  mpc_t *NLO,
  // INPUT
  mpc_t *LO, struct poly_frac *pf, mpc_t *roots, int pow
);

void poly_frac_mat_extract_NLOc(
  // OUTPUT
  mpc_t **NLO,
  // INPUT
  mpfr_t *LOre, mpfr_t *LOim,
  struct poly_frac  **pfmat, mpc_t *roots, int pow,
  int dim1, int dim2
);

void poly_frac_extract_order(
  // OUTPUT
  mpc_t *out,
  // INPUT
  struct poly_frac *pf, mpc_t *roots, int nroots, int order
);

void poly_frac_derivative(
  struct poly_frac *pfout,
  struct poly_frac *pfin,
  mpc_t *roots
);

void mpc_rk2_mul_poly_frac_rk2(
  struct poly_frac **pf_tens_out,
  mpc_t **mpc_tens,
  struct poly_frac **pf_tens_in,
  int dim,
  mpc_t *roots
);

void poly_frac_rk2_mul_mpc_rk2(
  struct poly_frac **pf_tens_out,
  struct poly_frac **pf_tens_in,
  mpc_t **mpc_tens,
  int dim,
  mpc_t *roots
);

void poly_frac_rk2_mul_pf_rk1(
  struct poly_frac *pf_tens_out,
  struct poly_frac **pf_tens_in1,
  struct poly_frac *pf_tens_in2,
  int dim1, int dim2,
  mpc_t *roots
);

void poly_frac_rk2_mul_pf_rk2(
  struct poly_frac **pf_tens_out,
  struct poly_frac **pf_tens_in1,
  struct poly_frac **pf_tens_in2,
  int dim1, int dim, int dim2,
  mpc_t *roots
);

void pf_find_profile(
  // OUTPUT
  int **&profile, int **&subblock, int &nblocks,
  // INPUT
  struct poly_frac **mat, int g_dim
);

int poly_frac_Poinc_rank(
  struct poly_frac **pfmat, int dim1, int dim2
);

int poly_frac_Poinc_rank_sh(
  struct poly_frac **pfmat, int dim1, int dim2, int sh
);

int poly_frac_k2_num_deg_max(
  struct poly_frac **pfmat, int dim1, int dim2
);

void poly_frac_k2_num_exp2_avg(
  // OUTPUT
  double *num_exp2_avg,
  // INPUT
  struct poly_frac **pfmat, int dim1, int dim2,
  int Poinc_rank, int num_deg_max
);

void poly_frac_roots_append_zero(
  // IN-OUT
  struct poly_frac *pf,
  // INPUT
  int nroots
);

void poly_frac_rk2_slice(
  struct poly_frac ***sub_tens,
  struct poly_frac ***tens,
  int dim1_start, int dim1_end, int dim2_start
);

void poly_frac_rk2_unslice(
  struct poly_frac ***tens,
  int dim1_start, int dim1_end, int dim2_start
);

void poly_frac_rk2_subtens(
  struct poly_frac ***sub_tens,
  struct poly_frac ***tens,
  int dim1_start, int dim1_end, int dim2_start
);

int poly_frac_cmp_pf(
  struct poly_frac *pf1,
  struct poly_frac *pf2
);

void root_prof_to_poly_frac_den(
  struct poly_frac *pf,
  int *root_prof, int deg,
  int nroots //, int zero_label
);

void root_prof_rk2_to_poly_frac_den_rk2(
  struct poly_frac **pf_tens,
  int ***root_prof,
  int nroots, // int zero_label,
  int dim1, int dim2
);

void poly_frac_rk2_to_file(
  char *file_name,
  struct poly_frac **pfmat, int dim1, int dim2
);

void poly_frac_rk2_from_file(
  char *file_name,
  struct poly_frac **pfmat, int dim1, int dim2
);

//////
// LOCAL PRECISION
//////
void poly_frac_set_pf_wp2(
	int wp2,
  struct poly_frac *pfout,
  struct poly_frac *pfin
);

void poly_frac_add_pf_wp2(
	int wp2,
  struct poly_frac *pfout,
  struct poly_frac *pfin1,
  struct poly_frac *pfin2,
  mpc_t *roots,
  mpc_t *flop_mpc
);

void poly_frac_mul_pf_wp2(
	int wp2,
  struct poly_frac *pfout,
  struct poly_frac *pfin1,
  struct poly_frac *pfin2
);

void poly_frac_pow_ui_wp2(
	int wp2,
  struct poly_frac *pfout,
  struct poly_frac *pfin,
  int pow
);

void poly_frac_set_mpc_wp2(
	int wp2,
  // OUTPUT
  struct poly_frac *pf,
  // INPUT
  mpc_t *mpcin,
  int nroots
);

void rel_err_poly_frac_add_pf(
  struct poly_frac *pfout,
  struct poly_frac *pfin1,
  struct poly_frac *pfin2,
  mpc_t *roots,
  mpc_t *flop_mpc,
  int wp_bin
);

void rel_err_poly_frac_mul_pf(
  struct poly_frac *pfout,
  struct poly_frac *pfin1,
  struct poly_frac *pfin2,
  int wp_bin
);

void rel_err_poly_frac_rk2_mul_pf_rk1(
  struct poly_frac *pf_tens_out,
  struct poly_frac **pf_tens_in1,
  struct poly_frac *pf_tens_in2,
  int dim1, int dim2,
  mpc_t *roots,
  int wp_bin
);

void rel_err_poly_frac_rk2_mul_pf_rk2(
  struct poly_frac **pf_tens_out,
  struct poly_frac **pf_tens_in1,
  struct poly_frac **pf_tens_in2,
  int dim1, int dim, int dim2,
  mpc_t *roots,
  int wp_bin
);

void rel_err_poly_frac_rk2_mul_mpc_rk2(
  struct poly_frac **pf_tens_out,
  struct poly_frac **pf_tens_in,
  mpc_t **mpc_tens,
  int dim,
  mpc_t *roots,
  int wp_bin
);

void rel_err_mpc_rk2_mul_poly_frac_rk2(
  struct poly_frac **pf_tens_out,
  mpc_t **mpc_tens,
  struct poly_frac **pf_tens_in,
  int dim,
  mpc_t *roots,
  int wp_bin
);

void rel_err_poly_frac_rk2_add_pf_rk2(
  struct poly_frac **pf_tens_out,
  struct poly_frac **pf_tens_in1,
  struct poly_frac **pf_tens_in2,
  mpc_t *roots, mpc_t *flop_mpc,
  int dim1, int dim2,
  int wp_bin
);

void rel_err_poly_frac_extract_LO(
  // OUTPUT
  mpfr_t *out_re, mpfr_t *out_im,
  // INPUT
  struct poly_frac *pf, mpc_t *roots, int pow,
  int wp_bin
);

void rel_err_poly_frac_mat_extract_LO(
  // OUTPUT
  mpfr_t *LOre, mpfr_t *LOim,
  // INPUT
  struct poly_frac  **pfmat, mpc_t *roots, int pow,
  int dim1, int dim2,
  int wp_bin
);

void rel_err_poly_frac_extract_NLO(
  // OUTPUT
  mpc_t *NLO,
  // INPUT
  mpc_t *LO, struct poly_frac *pf, mpc_t *roots, int pow,
  int wp_bin
);

void rel_err_poly_frac_mat_extract_NLOc(
  // OUTPUT
  mpc_t **NLO,
  // INPUT
  mpfr_t *LOre, mpfr_t *LOim,
  struct poly_frac  **pfmat, mpc_t *roots, int pow,
  int dim1, int dim2,
  int wp_bin
);

void rel_err_poly_frac_extract_order(
  // OUTPUT
  mpc_t *out,
  // INPUT
  struct poly_frac *pf, mpc_t *roots, int nroots, int order,
  int wp_bin
);

void rel_err_poly_frac_rk2_infty(
  // OUTPUT
  struct poly_frac **pfout, mpc_t *out_roots, int *out_zero_label,
  // INPUT
  struct poly_frac **pfin, mpc_t *in_roots,
  int dim1, int dim2, int nroots,
  int wp_bin
);

void rel_err_poly_frac_rk2_shift(
  struct poly_frac **pfout,
  struct poly_frac **pfin,
  mpc_t *shift,
  int swap_root_lab,
  int dim1, int dim2,
  int wp_bin
);

void rel_err_center_around_pole(
  // OUTPUT
  struct poly_frac **sh_pfmat, mpc_t *sh_roots,
  // INPUT
  struct poly_frac **pfmat, mpc_t *roots, int label,
  int nroots, int dim1, int dim2,
  int wp_bin
);

void mpc_rk1_prune_re_im(
  // IN-OUT
  mpc_t *vec,
  // INPUT
  int wp_bin, int dim
);

void rel_err_mpc_rk2_mul_mpc_rk1(
  // OUTPUT
  mpc_t *out_tens,
  // INPUT
  mpc_t **tens1, mpc_t *tens2,
  int dim1, int dim2,
  int wp_bin
);

#endif