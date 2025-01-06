#ifndef MP_ROOTS_POLY_H
#define MP_ROOTS_POLY_H

void mp_laguer(mpc_t *a, int m, mpc_t *x, int *its, int wp2, mpfr_t EPSS, mpfr_t *err);

void mp_zroots(mpc_t *a, int m, mpc_t *roots, int polish, int wp2, mpfr_t mpfr_tol, mpfr_t *err);

int test_root(
	mpc_t *a, int m, mpc_t *x, int wp2, mpfr_t EPSS, mpfr_t *err,
	mpc_t *b, mpc_t *d, mpc_t *f, mpfr_t *abx, mpfr_t *tmpfr
);

void deflation(mpc_t *ad, int j, mpc_t *x, int wp2);

int mp_laguer_impr(mpc_t *a, int m, mpc_t *x, int *its, int wp2, mpfr_t EPSS, mpfr_t *err, int mul);

void mp_zroots_impr(
	// OUTPUT
	mpc_t *roots, mpfr_t *err,
	// INPUT
	mpc_t *a, int m, mpc_t *dens,
	int polish, int wp2, mpfr_t mpfr_tol, int curr_mul
);

void mp_zroots_impr_out_mul(
	// OUTPUT
	mpc_t *roots, mpfr_t *err, int *out_mul,
	// INPUT
	mpc_t *a, int m, mpc_t *dens,
	int polish, int wp2, mpfr_t mpfr_tol, int curr_mul
);

#endif