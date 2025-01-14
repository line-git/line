#ifndef SETUP_H
#define SETUP_H

#include <complex>
typedef complex<double> complex_d;

#define MAX_VALUE_LEN 2000
#define MAX_POW_DIGITS 5
#define MAX_SYM_LEN 50

void read_mat_to_str(
  // OUTPUT
  char ****matrix, int *dim,
  // INPUT
  char *filepath
);


void compute_setup(
  // OUTPUT
  int *workprec, int *xorder, int *number, char ***eps_list,
  // INPUT
  int order, int precision , int loops
);

void mpfr_tol_set();

void mpfr_tol_set_wp(int workprec);

void info_to_log(
  // INPUT
  FILE *logfptr,
  int loops, int order, int precision, int eta_ord, int workprec, int wp2
);

void parse_list(
  // OUTPUT
  char ***out, int *dim,
  // INPUT
  char *list
);

int sym_to_idx(
  char *key, char **symbol, int ninvs
);

int end_at_char(
  // IN-OUT
  char *line,
  // INPUT
  char ch
);

int read_param_str(
  // OUTPUT
  char **value,
  // INPUT
  char *key, FILE *fptr
);

int read_param_int(
  // OUTPUT
  int *value,
  // INPUT
  char *key, FILE *fptr
);

int read_list_str(
  // OUTPUT
  char ***list,
  // INPUT
  int *dim, char *key, FILE *fptr
);

int read_PS_point(
  // OUTPUT
  char ***kin,
  // INPUT
  int *ninvs, char **symbols, char *key, FILE *fptr
);

int read_bound_factors(
  // OUTPUT
  char ****out, int *nfac,
  // INPUT
  int nfac_max, int nparam, char ***symbols, char *key, FILE *fptr
);

void load_parameters(
  // OUTPUT
  int *loops, int *order, int *precision,
  char **prev_run,
  char ***kin,
  // INPUT
  int *ninvs, char ***symbols,
  char *filename
);

void param_string(
  // OUTPUT
  char **out,
  // INPUT
  int num, char **keys, int *values
);

void point_string(
	// OUTPUT
	char **out,
	// INPUT
	// int order, int precision,
	int ninvs, char** symbols, char** PS
);

void target_string(
	// OUTPUT
	char **out,
	// INPUT
	int nparam, char **param_keys, int *param_values,
	int ninvs, char** symbols, char** PS
);

char* hash_to_string(unsigned char *hash);

void calculate_hash(
	// OUTPUT
	char **out,
	// INPUT
	char *in
);

int file_exists(const char *path);

int directory_exists(const char *path);

int make_dir(char* path);

int remove_dir(const char* path);

void log_param(
  char *key, int value,
  char *filepath, const char* mode
);

void log_params(
  int num, char **keys, int *values,
  char *filepath, const char* mode
);

void log_point(
  int ninvs, char **symbols, char **kin,
  char *key, char *filepath, const char* mode
);

void prepare_prev_run(
  // OUTPUT
  char **filepath_prev_run_PS, char **filepath_prev_run_sol,
  // INPUT
  char *dir_parent, char *prev_run
);

void load_symbols(
	// OUTPUT
	int *ninvs, char ***symbols,
	// IN-OUT
	char *filepath_vars
);

void detect_masses(
  // OUTPUT
  int *nmass, int **mass, int **is_mass,
  // INPUT
  int ninvs, char **symbols
);

void load_DE_matrices(
  // OUTPUT
  char ****mats_str, int *dim,
  // INPUT
  int ninvs, char **symbols, int *skip_inv,
  char *filepath_mats, char *file_ext
);

void str_rk1_to_mpc_rk1 (
  // OUT
  mpc_t *out,
  // INPUT
  char **kin, int ninvs
);

void detect_active_invs(
  // OUTPUT
  int *skip_inv,
  // INPUT
  int ninvs, mpc_t *PS_ini, mpc_t *PS_fin
);

void load_PS(
  // OUT
  mpc_t *PS_ini, mpc_t *PS_fin,
  int *skip_inv, char ***kin,
  // INPUT
  int ninvs, char *filepath_PS
);

void load_PS_prev_run(
  // OUT
  mpc_t *PS_ini, mpc_t *PS_fin,
  int *skip_inv, char ***kin,
  // INPUT
  int ninvs, char *filepath_PS
);

void linearize_masses_mpc(
  // IN-OUT
  mpc_t *PS_ini, mpc_t *PS_fin,
  // INPUT
  int nmass, int *mass
);

void linearize_masses_str(
  // OUTPUT
  char ***lin_symbols,
  // INPUT
  int ninvs, char **symbols, int *is_mass
);

void load_branch_expr_str(
	// OUTPUT
	int *nexpr, char ***expr,
	// IN-OUT
	int linearize_mass, char *filepath
);

void bound_behav_from_file(
  char *filepath,
  int ***bound_behav, int ***mi_eig, int **mi_eig_num, int *num_region_classes,
  int dim
);

void bound_behav_from_file_old(
  char *file_name,
  int ***pow_behav, int ***mi_eig, int **mi_eig_num, int *num_classes,
  int num_mi
);

void generate_coordinates(
  // OUTPUT
  double *pts,
  int *npts,
  // INPUT
  complex_d *poles, int npoles,
  int RunRadius, int len_max, int mall
);

void generate_regular_points(
  // OUTPUT
  double *pts, int *npts,
  // INPUT
  complex_d *poles, int npoles,
  int RunRadius, int len_max,
  complex_d pt1, complex_d pt2
);

void get_path_PS(
  // OUTPUT
  mpc_t **path, int **path_tags, int *neta_vals, int **sing_lab, int *nsings,
  // INPUT
  mpc_t *roots, int nroots, int zero_label
);

void get_path_PS_mp(
  // OUTPUT
  mpc_t **path, int **path_tags, int *neta_vals, int **sing_lab, int *nsings,
  // INPUT
  mpc_t *roots, int nroots, int zero_label
);

void get_path(mpc_t **eta_values, int* neta_vals, complex_d *poles, int npoles);

void get_path_PS_infty(
  // OUTPUT
  mpc_t **path, int **path_tags, int *neta_vals, int **sing_lab, int *nsings,
  // INPUT
  mpc_t *roots, int nroots, int zero_label
);

void get_path_PS_infty_mp(
  // OUTPUT
  mpc_t **path, int **path_tags, int *neta_vals, int **sing_lab, int *nsings,
  // INPUT
  mpc_t *roots, int nroots, int zero_label
);

void path_to_math(
  mpc_t *path, int *path_tags, int neta_values,
  mpc_t *roots, int nroots
);

void wrt_cmp_path(
  int ep, mpc_t **path, int **path_tags, int eps_num, int *neta_values,
  int *nsings, int **sing_lab, int **perm,
  char *file_ext, char *filepath_path, char *filepath_path_tags, char* filepath_sing_lab,
  int opt_write
);

//////
// NESTED LISTS
//////
typedef struct nlist {
  int nitems;
  nlist *item;
  char *str;
} nlist;

void nlist_self_parse(nlist *nl);

int nlist_read_file(
  // OUTPUT
  nlist *nl,
  // INPUT
  char *key, FILE *fptr
);

void nlist_free(nlist *nl);

void nlist_rk1_free(nlist *nl, int dim);

#endif