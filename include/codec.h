#ifndef CODEC_H
#define CODEC_H

void decode_tree(
	// OUTPUT
	mpc_t *out,
	// INPUT
	mpc_t *in,
	char **s_tree,
	char *sep
);

void decode_tree_mpq(
	// OUTPUT
	mpq_t *out,
	// INPUT
	mpq_t *in,
	char **s_tree,
	char *sep
);

void decode_tree_poly(
	// OUTPUT
	struct cpoly *out,
	// INPUT
	struct cpoly *in,
	char **s_tree,
	char *sep
);

void poly_mpq_find_roots(
	// OUTPUT
	int *found_mul_roots,
	int *vdeg, mpc_t **roots, mpfr_t **tols,
	// IN-OUT
	// int *in_nroots, mpc_t **in_roots, mpfr_t **in_tols,
	// INPUT
	mpq_t *coeffs_mpq, int deg,
	int wp2, mpfr_t mpfr_tol, int incr_prec, int curr_mul,
	int defl_nroots, int *defl_mults, mpc_t *defl_roots
);

void generate_PS_pf(
	// OUTPUT
	struct poly_frac *out,
	// INPUT
	char ***kin, int skip_inv, int is_mass, int s, int nroots,
	// int ninvs, ***char kin, int *skip_inv, int *is_mass
	int wp2
);

void decode_tree_pf(
	// OUTPUT
	struct poly_frac *out,
	// IN-OUT
	int *nroots, mpc_t **roots, mpfr_t **tols,
	// INPUT
	struct poly_frac *in,
	char ***kin, int *skip_inv, int *is_mass,
	char **s_tree, char *sep,
	int wp2, mpfr_t mpfr_tol, int incr_prec
);

void lnode_to_pf(
	// OUTPUT
	struct poly_frac *out,
	// IN-OUT
	int *nroots, mpc_t **roots, mpfr_t **tols,
	// INPUT
	struct poly_frac *in,
	char ***kin, int *skip_inv, int *is_mass,
	struct lnode *nd,
	int wp2, mpfr_t mpfr_tol, int incr_prec
);

#endif