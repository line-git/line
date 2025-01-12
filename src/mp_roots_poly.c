#include <stdio.h>
#include <stdlib.h>
// #include <math.h>
#include "mpc.h"
#include "global_vars.h"
#define MAXM 100
#define MR 8
#define MT 10000
#define MT_red 10
#define MAXIT (MT*MR)
#define MAXIT_red (MT_red*MR)
double laguer_frac[MR+1] = {0.0,0.5,0.25,0.75,0.13,0.38,0.62,0.88,1.0};

void mp_laguer(mpc_t *a, int m, mpc_t *x, int *its, int wp2, mpfr_t EPSS, mpfr_t *err)
{
	// printf("I'm in mp_laguer with wp2 = %d\n", wp2);
	int iter, j;
	mpfr_t abx, abp, abm;
	mpfr_init2(abx, wp2);
	mpfr_init2(abp, wp2);
	mpfr_init2(abm, wp2);
	mpc_t dx, x1, b, d, f, g, h, sq, gp, gm, g2;
	mpc_init3(dx, wp2, wp2);
	mpc_init3(x1, wp2, wp2);
	mpc_init3(b, wp2, wp2);
	mpc_init3(d, wp2, wp2);
	mpc_init3(f, wp2, wp2);
	mpc_init3(g, wp2, wp2);
	mpc_init3(h, wp2, wp2);
	mpc_init3(sq, wp2, wp2);
	mpc_init3(gp, wp2, wp2);
	mpc_init3(gm, wp2, wp2);
	mpc_init3(g2, wp2, wp2);
	// mpfr_t laguer_frac[MR+1];
	// for (int i=0; i<=MR; i++) {
	// 	mpfr_init2(laguer_frac[i], wp2);
	// }
	// mpfr_set_d(laguer_frac[0], 0.0, MPFR_RNDN);
	// mpfr_set_d(laguer_frac[1], 0.5, MPFR_RNDN);
	// mpfr_set_d(laguer_frac[2], 0.25, MPFR_RNDN);
	// mpfr_set_d(laguer_frac[3], 0.75, MPFR_RNDN);
	// mpfr_set_d(laguer_frac[4], 0.13, MPFR_RNDN);
	// mpfr_set_d(laguer_frac[5], 0.38, MPFR_RNDN);
	// mpfr_set_d(laguer_frac[6], 0.62, MPFR_RNDN);
	// mpfr_set_d(laguer_frac[7], 0.88, MPFR_RNDN);
	// mpfr_set_d(laguer_frac[8], 1.0, MPFR_RNDN);
	// static double laguer_frac[MR+1] = {0.0,0.5,0.25,0.75,0.13,0.38,0.62,0.88,1.0};

	mpfr_t tmpfr;
	mpfr_init2(tmpfr, wp2);
	mpc_t tmpc;
	mpc_init3(tmpc, wp2, wp2);
	for (iter = 1; iter <= MAXIT; iter++){
		// printf("iter = :%d\n", iter);
		*its = iter;
		mpc_set(b, a[m], MPFR_RNDN);
		mpc_abs(*err, b, MPFR_RNDN);
		mpc_set_d_d(d, 0, 0, MPFR_RNDN);
		mpc_set_d_d(f, 0, 0, MPFR_RNDN);
		mpc_abs(abx, *x, MPFR_RNDN);
		for (j=m-1; j>=0; j--) {
			mpc_fma(f, *x, f, d, MPFR_RNDN);
			mpc_fma(d, *x, d, b, MPFR_RNDN);
			mpc_fma(b, *x, b, a[j], MPFR_RNDN);
			mpc_abs(tmpfr, b, MPFR_RNDN);
			mpfr_fma(*err, abx, *err, tmpfr, MPFR_RNDN);
		}
		mpfr_mul(*err, *err, EPSS, MPFR_RNDN);
		// printf("cabs(b) = ");
		// mpfr_out_str(stdout, 10, 0, tmpfr, MPFR_RNDN);
		// printf("\n");
		// printf("*err = ");
		// mpfr_out_str(stdout, 10, 0, *err, MPFR_RNDN);
		// printf("\n");
		if (mpfr_lessequal_p(tmpfr, *err)) {
			// we are on the root
			// printf("root:\n");
			// mpc_out_str(stdout, 10, 0, *x, MPFR_RNDN); printf("\n");
			// printf("polynomial at root:\n");
			// mpfr_out_str(stdout, 10, 0, tmpfr, MPFR_RNDN); printf("\n");
			// printf("roundoff error:\n");
			// mpfr_out_str(stdout, 10, 0, *err, MPFR_RNDN); printf("\n");
			goto goto_return;
			// return;
		}
		mpc_div(g, d, b, MPFR_RNDN);
		mpc_sqr(g2, g, MPFR_RNDN);
		mpc_mul_ui(tmpc, f, 2, MPFR_RNDN);
		mpc_div(tmpc, tmpc, b, MPFR_RNDN);
		mpc_sub(h, g2, tmpc, MPFR_RNDN);
		mpc_mul_ui(tmpc, h, m, MPFR_RNDN);
		mpc_sub(tmpc, tmpc, g2, MPFR_RNDN);
		mpc_mul_ui(tmpc, tmpc, m-1, MPFR_RNDN);
		mpc_sqrt(sq, tmpc, MPFR_RNDN);
		mpc_add(gp, g, sq, MPFR_RNDN);
		mpc_sub(gm, g, sq, MPFR_RNDN);
		mpc_abs(abp, gp, MPFR_RNDN);
		mpc_abs(abm, gm, MPFR_RNDN);
		if (mpfr_lessequal_p(abp, abm)) {
			mpc_set(gp, gm, MPFR_RNDN);
		}
		mpfr_max(tmpfr, abp, abm, MPFR_RNDN);
		if (mpfr_sgn(tmpfr) > 0) {
			mpc_pow_si(dx, gp, -1, MPFR_RNDN);
			mpc_mul_ui(dx, dx, m, MPFR_RNDN);
		} else {
			mpc_set_ui_ui(tmpc, 0, 1, MPFR_RNDN);
			mpc_mul_ui(tmpc, tmpc, iter, MPFR_RNDN);
			mpc_exp(dx, tmpc, MPFR_RNDN);
			mpfr_add_ui(tmpfr, abx, 1, MPFR_RNDN);
			mpc_mul_fr(dx, dx, tmpfr, MPFR_RNDN);
		}
		
		mpc_sub(x1, *x, dx, MPFR_RNDN);
		mpc_sub(tmpc, *x, x1, MPFR_RNDN);
		if (mpfr_equal_p(mpc_realref(*x), mpc_realref(x1)) && mpfr_equal_p(mpc_imagref(*x), mpc_imagref(x1))) {
			goto goto_return;
			return;
		}
		
		if (iter % MT) {
			mpc_set(*x, x1, MPFR_RNDN);
		} else {
			// mpc_mul_fr(tmpc, dx, laguer_frac[iter/MT], MPFR_RNDN);
			mpfr_mul_d(mpc_realref(tmpc), mpc_realref(dx), laguer_frac[iter/MT], MPFR_RNDN);
			mpfr_mul_d(mpc_imagref(tmpc), mpc_imagref(dx), laguer_frac[iter/MT], MPFR_RNDN);
			mpc_sub(*x, *x, tmpc, MPFR_RNDN);
		} 
	}
	perror("too many iterationc in mp_laguer");
	exit(1);

	goto_return:
	// FREE
	mpfr_clear(abx);
	mpfr_clear(abp);
	mpfr_clear(abm);
	mpfr_clear(tmpfr);
	mpc_clear(dx);
	mpc_clear(x1);
	mpc_clear(b);
	mpc_clear(d);
	mpc_clear(f);
	mpc_clear(g);
	mpc_clear(h);
	mpc_clear(sq);
	mpc_clear(gp);
	mpc_clear(gm);
	mpc_clear(g2);
	mpc_clear(tmpc);
	return;
}


void mp_zroots(mpc_t *a, int m, mpc_t *roots, int polish, int wp2, mpfr_t mpfr_tol, mpfr_t *err)
{
	mpfr_t EPSS;
	mpfr_init2(EPSS, wp2);
	mpfr_mul_ui(EPSS, mpfr_tol, 10, MPFR_RNDN);
	mpfr_t EPS;
	mpfr_init2(EPS, wp2);
	mpfr_mul_ui(EPS, EPSS, 20, MPFR_RNDN);

	mpfr_t tmpfr;
	mpfr_init2(tmpfr, wp2);
	mpfr_t prev_err;
	mpfr_init2(prev_err, wp2);

	int i, its, j, jj;
	mpc_t x, b, c;
	mpc_init3(x, wp2, wp2);
	mpc_init3(b, wp2, wp2);
	mpc_init3(c, wp2, wp2);
	mpc_t ad[m+1]; // = (mpc_t*) malloc((m+1)*sizeof(mpc_t));
	for (j=0; j<=m; j++){
		mpc_init3(ad[j], wp2, wp2);
		mpc_set(ad[j], a[j], MPFR_RNDN);
	}

	mpfr_t abs_re, abs_im;
	mpfr_init2(abs_re, wp2);
	mpfr_init2(abs_im, wp2);
	mpfr_set_d(prev_err, 0, MPFR_RNDN);
	for (j=m; j>=1; j--) {
		// printf("\n");
		// printf("j = %d\n", j);
		// printf("polyomial ad\n");
		// for (int k=0; k<=j; k++) {
		// 	mpc_out_str(stdout, 10, 0, ad[k], MPFR_RNDN);
		// 	printf("\n");
		// }
		// printf("\n");
		mpc_set_d_d(x, 0, 0, MPFR_RNDN);
		mp_laguer(ad, j, &x, &its, wp2, EPSS, &err[j]);
		mpfr_add(err[j], err[j], prev_err, MPFR_RNDN);
		mpfr_set(prev_err, err[j], MPFR_RNDN);
		// printf("output of lauguer:\n");
		// printf("x = ");
		// mpc_out_str(stdout, 10, 0, x, MPFR_RNDN);
		// printf("\n");
		// printf("\n");
		mpfr_abs(abs_re, mpc_realref(x), MPFR_RNDN);
		mpfr_abs(abs_im, mpc_imagref(x), MPFR_RNDN);
		mpfr_mul(abs_re, abs_re, EPS, MPFR_RNDN);
		mpfr_mul_ui(abs_re, abs_re, 2, MPFR_RNDN);
		if (mpfr_lessequal_p(abs_im, abs_re)) {
			mpc_real(tmpfr, x, MPFR_RNDN);
			mpc_set_fr(x, tmpfr, MPFR_RNDN);
			// mpfr_set_d(mpc_imagref(x), 0, MPFR_RNDN);
		}
		mpc_set(roots[j], x, MPFR_RNDN);
		// mpc_out_str(stdout, 10, 0, roots[j], MPFR_RNDN); printf("\n");
		// forward deflation
		mpc_set(b, ad[j], MPFR_RNDN);
		for (jj=j-1; jj>=0; jj--) {
			mpc_set(c, ad[jj], MPFR_RNDN);
			mpc_set(ad[jj], b, MPFR_RNDN);
			mpc_fma(b, x, b, c, MPFR_RNDN);
			// mpc_out_str(stdout, 10, 0, ad[jj], MPFR_RNDN); printf("\n");
		}
	}

	if (polish) {
		mpfr_set_d(prev_err, 0, MPFR_RNDN);
		for (j=1; j<=m; j++) {
			mp_laguer(a, m, &roots[j], &its, wp2, EPSS, &err[j]);
			mpfr_add(err[j], err[j], prev_err, MPFR_RNDN);
			mpfr_set(prev_err, err[j], MPFR_RNDN);
		}
	}
	
	for (j=2; j<=m; j++) {
		mpc_set(x, roots[j], MPFR_RNDN);
		for (i=j-1; i>=1; i--) {
			if (mpfr_lessequal_p(mpc_realref(roots[i]), mpc_realref(x))) {
				break;
			}
			mpc_set(roots[i+1], roots[i], MPFR_RNDN);
		}
		mpc_set(roots[i+1], x, MPFR_RNDN);
	}
}


int test_root(
	mpc_t *a, int m, mpc_t *x, int wp2, mpfr_t EPSS, mpfr_t *err,
	mpc_t *b, mpc_t *d, mpc_t *f, mpfr_t *abx, mpfr_t *tmpfr
) {
	int j;

	mpc_set(*b, a[m], MPFR_RNDN);
	mpc_abs(*err, *b, MPFR_RNDN);
	mpc_set_d_d(*d, 0, 0, MPFR_RNDN);
	mpc_set_d_d(*f, 0, 0, MPFR_RNDN);
	mpc_abs(*abx, *x, MPFR_RNDN);
	for (j=m-1; j>=0; j--) {
		mpc_fma(*f, *x, *f, *d, MPFR_RNDN);
		mpc_fma(*d, *x, *d, *b, MPFR_RNDN);
		mpc_fma(*b, *x, *b, a[j], MPFR_RNDN);
		mpc_abs(*tmpfr, *b, MPFR_RNDN);
		mpfr_fma(*err, *abx, *err, *tmpfr, MPFR_RNDN);
	}
	mpfr_mul(*err, *err, EPSS, MPFR_RNDN);
	if (dbg) {
	printf("err = "); mpfr_out_str(stdout, 10, 0, *err, MPFR_RNDN); printf("\n");
	printf("tmpfr = "); mpfr_out_str(stdout, 10, 0, *tmpfr, MPFR_RNDN); printf("\n");
	}

	if (mpfr_lessequal_p(*tmpfr, *err)) {
		// we are on the root
		return 1;
	}
	return 0;
}


void deflation(mpc_t *ad, int j, mpc_t *x, int wp2) {
	int jj;

	// // copy input polynomial
	// for (jj=j; jj>=0; jj--) {
	// 	mpc_set(ad[jj], in[jj], MPFR_RNDN);
	// }

	mpc_t b, c;
	mpc_init3(b, wp2, wp2);
	mpc_init3(c, wp2, wp2);
	mpc_set(b, ad[j], MPFR_RNDN);
	for (jj=j-1; jj>=0; jj--) {
		mpc_set(c, ad[jj], MPFR_RNDN);
		mpc_set(ad[jj], b, MPFR_RNDN);
		mpc_fma(b, *x, b, c, MPFR_RNDN);
		// mpc_out_str(stdout, 10, 0, ad[jj], MPFR_RNDN); printf("\n");
	}

	// FREE
	mpc_clear(b);
	mpc_clear(c);
}


int mp_laguer_impr(mpc_t *a, int m, mpc_t *x, int *its, int wp2, mpfr_t EPSS, mpfr_t *err, int mul)
{
	int iter, j;
	mpfr_t abx, abp, abm;
	mpfr_init2(abx, wp2);
	mpfr_init2(abp, wp2);
	mpfr_init2(abm, wp2);
	mpc_t dx, x1, b, d, f, g, h, sq, gp, gm, g2;
	mpc_init3(dx, wp2, wp2);
	mpc_init3(x1, wp2, wp2);
	mpc_init3(b, wp2, wp2);
	mpc_init3(d, wp2, wp2);
	mpc_init3(f, wp2, wp2);
	mpc_init3(g, wp2, wp2);
	mpc_init3(h, wp2, wp2);
	mpc_init3(sq, wp2, wp2);
	mpc_init3(gp, wp2, wp2);
	mpc_init3(gm, wp2, wp2);
	mpc_init3(g2, wp2, wp2);
	static mpfr_t frac[MR+1];
	for (int i=0; i<=MR; i++) {
		mpfr_init2(frac[i], wp2);
	}
	mpfr_set_d(frac[0], 0.0, MPFR_RNDN);
	mpfr_set_d(frac[1], 0.5, MPFR_RNDN);
	mpfr_set_d(frac[2], 0.25, MPFR_RNDN);
	mpfr_set_d(frac[3], 0.75, MPFR_RNDN);
	mpfr_set_d(frac[4], 0.13, MPFR_RNDN);
	mpfr_set_d(frac[5], 0.38, MPFR_RNDN);
	mpfr_set_d(frac[6], 0.62, MPFR_RNDN);
	mpfr_set_d(frac[7], 0.88, MPFR_RNDN);
	mpfr_set_d(frac[8], 1.0, MPFR_RNDN);
	// static double frac[MR+1] = {0.0,0.5,0.25,0.75,0.13,0.38,0.62,0.88,1.0};

	mpfr_t tmpfr;
	mpfr_init2(tmpfr, wp2);
	mpc_t tmpc;
	mpc_init3(tmpc, wp2, wp2);
	mpc_t *defl;
	defl = (mpc_t*) malloc((m+1)*sizeof(mpc_t));
	int mm, converged;
	for (int k=0; k<=m; k++) {
		mpc_init3(defl[k], wp2, wp2);
	}
	for (iter = 1; iter <= MAXIT_red; iter++){
		*its = iter;

		// // perform enhanced test
		for (int k=0; k<=m; k++) {
			mpc_set(defl[k], a[k], MPFR_RNDN);
		}

		for (mm=0; mm<mul; mm++) {
			// printf("mm = %d\n", mm);
			converged = test_root(
				defl, m-mm, x, wp2, EPSS, err,
				&b, &d, &f, &abx, &tmpfr
			);
			if (!converged) {
				// printf("roundoff error:\n");
				// mpfr_out_str(stdout, 10, 0, *err, MPFR_RNDN); printf("\n");
				// printf("poly at root:\n");
				// mpfr_out_str(stdout, 10, 0, tmpfr, MPFR_RNDN); printf("\n");
				break;
			}
			deflation(defl, m-mm, x, wp2);
		}
		if (mm == mul) {
			// printf("root:\n");
			// mpc_out_str(stdout, 10, 0, *x, MPFR_RNDN); printf("\n");
			// printf("polynomial at root:\n");
			// mpfr_out_str(stdout, 10, 0, tmpfr, MPFR_RNDN); printf("\n");
			// printf("roundoff error:\n");
			// mpfr_out_str(stdout, 10, 0, *err, MPFR_RNDN); printf("\n");
			return 1;
		}
		
		// converged = test_root(
		// 	defl, m, x, wp2, EPSS, err,
		// 	&b, &d, &f, &abx, &tmpfr
		// );
		// if (converged) {
		// 	return 1;
		// }

		mpc_div(g, d, b, MPFR_RNDN);
		mpc_sqr(g2, g, MPFR_RNDN);
		mpc_mul_ui(tmpc, f, 2, MPFR_RNDN);
		mpc_div(tmpc, tmpc, b, MPFR_RNDN);
		mpc_sub(h, g2, tmpc, MPFR_RNDN);
		mpc_mul_ui(tmpc, h, m, MPFR_RNDN);
		mpc_sub(tmpc, tmpc, g2, MPFR_RNDN);
		mpc_mul_ui(tmpc, tmpc, m-1, MPFR_RNDN);
		mpc_sqrt(sq, tmpc, MPFR_RNDN);
		mpc_add(gp, g, sq, MPFR_RNDN);
		mpc_sub(gm, g, sq, MPFR_RNDN);
		mpc_abs(abp, gp, MPFR_RNDN);
		mpc_abs(abm, gm, MPFR_RNDN);
		if (mpfr_lessequal_p(abp, abm)) {
			mpc_set(gp, gm, MPFR_RNDN);
		}
		mpfr_max(tmpfr, abp, abm, MPFR_RNDN);
		if (mpfr_sgn(tmpfr) > 0) {
			mpc_pow_si(dx, gp, -1, MPFR_RNDN);
			mpc_mul_ui(dx, dx, m, MPFR_RNDN);
		} else {
			mpc_set_ui_ui(tmpc, 0, 1, MPFR_RNDN);
			mpc_mul_ui(tmpc, tmpc, iter, MPFR_RNDN);
			mpc_exp(dx, tmpc, MPFR_RNDN);
			mpfr_add_ui(tmpfr, abx, 1, MPFR_RNDN);
			mpc_mul_fr(dx, dx, tmpfr, MPFR_RNDN);
		}
		
		mpc_sub(x1, *x, dx, MPFR_RNDN);
		mpc_sub(tmpc, *x, x1, MPFR_RNDN);
		if (mpfr_equal_p(mpc_realref(*x), mpc_realref(x1)) && mpfr_equal_p(mpc_imagref(*x), mpc_imagref(x1))) {
			mpfr_set(*err, tmpfr, MPFR_RNDN);
		}
		
		if (iter % MT_red) {
			mpc_set(*x, x1, MPFR_RNDN);
		} else {
			mpc_mul_fr(tmpc, dx, frac[iter/MT_red], MPFR_RNDN);
			mpc_sub(*x, *x, tmpc, MPFR_RNDN);
		} 
	}
	perror("too many iterationc in mp_laguer");
	return 0;
}


void mp_zroots_impr2(mpc_t *a, int m, mpc_t *roots, int polish, int wp2, mpfr_t mpfr_tol, mpfr_t *err)
{
	// mpfr_set_ui(mpfr_one, 1, MPFR_RNDN);
	mpfr_t EPSS;
	mpfr_init2(EPSS, wp2);
	mpfr_mul_ui(EPSS, mpfr_tol, 10, MPFR_RNDN);
	mpfr_t EPS;
	mpfr_init2(EPS, wp2);
	mpfr_mul_ui(EPS, EPSS, 20, MPFR_RNDN);

	mpfr_t EPSS_large;
	mpfr_init2(EPSS_large, wp2);
	mpfr_set_ui(EPSS_large, 10, MPFR_RNDN);
	mpfr_pow_si(EPSS_large, EPSS_large, wp2*3/10/10, MPFR_RNDN);
	mpfr_mul(EPSS_large, EPSS, EPSS_large, MPFR_RNDN);
	// mpfr_out_str(stdout, 10, 0, EPSS_large, MPFR_RNDN); printf("\n");

	mpfr_t tmpfr;
	mpfr_init2(tmpfr, wp2);
	mpfr_t prev_err;
	mpfr_init2(prev_err, wp2);

	int i, its, j, jj;
	mpc_t x, b, c;
	mpc_init3(x, wp2, wp2);
	mpc_init3(b, wp2, wp2);
	mpc_init3(c, wp2, wp2);
	mpc_t ad[m+1]; // = (mpc_t*) malloc((m+1)*sizeof(mpc_t));
	for (j=0; j<=m; j++){
		mpc_init3(ad[j], wp2, wp2);
		mpc_set(ad[j], a[j], MPFR_RNDN);
	}

	mpfr_t abs_re, abs_im;
	mpfr_init2(abs_re, wp2);
	mpfr_init2(abs_im, wp2);
	mpfr_set_d(prev_err, 0, MPFR_RNDN);
	int mul, converged;
	mpc_t *defl = (mpc_t*) malloc((m+1)*sizeof(mpc_t));	
	for (int k=0; k<=m; k++) {
		mpc_init3(defl[k], wp2, wp2);
	}
	mpfr_t err_tmp, tmpfr_tmp, abx_tmp, EPSS_red, EPSS_red_tmp;
	mpfr_init2(err_tmp, wp2);
	mpfr_init2(tmpfr_tmp, wp2);
	mpfr_init2(abx_tmp, wp2);
	mpfr_init2(EPSS_red, wp2);
	mpfr_init2(EPSS_red_tmp, wp2);
	mpc_t b_tmp, d_tmp, f_tmp;
	mpc_init3(b_tmp, wp2, wp2);
	mpc_init3(d_tmp, wp2, wp2);
	mpc_init3(f_tmp, wp2, wp2);
	for (j=m; j>=1; j--) {
		// printf("\n");
		// printf("j = %d\n", j);
		// printf("polyomial ad\n");
		// for (int k=0; k<=j; k++) {
		// 	mpc_out_str(stdout, 10, 0, ad[k], MPFR_RNDN);
		// 	printf("\n");
		// }
		// printf("\n");

		// if (j == 5) {
		// 	for (int k=0; k<=j; k++) {
		// 		mpc_out_str(stdout, 10, 0, ad[k], MPFR_RNDN); printf("\n");
		// 	}
		// }
		mpc_set_d_d(x, 0, 0, MPFR_RNDN);
		mp_laguer(ad, j, &x, &its, wp2, EPSS, &err[j]);
		// printf("root:\n");
		// mpc_out_str(stdout, 10, 0, x, MPFR_RNDN); printf("\n");

		for (int k=0; k<=j; k++) {
			mpc_set(defl[k], ad[k], MPFR_RNDN);
		}

		// DETECT MULTIPLICITY
		mul=1;
		while (1) {
			if (mul==j) {
				break;
			}
			// printf("try mul = %d\n", mul);
			deflation(defl, j-mul+1, &x, wp2);
			// if (j == 5) {
			// 	for (int k=0; k<j; k++) {
			// 		mpc_out_str(stdout, 10, 0, defl[k], MPFR_RNDN); printf("\n");
			// 	}
			// }
			mpfr_rootn_ui(EPSS_red, EPSS, mul+1, MPFR_RNDN);
			mpfr_mul_d(EPSS_red, EPSS_red, 1e10, MPFR_RNDN);
			// printf("reduced precision:\n");
			// mpfr_out_str(stdout, 10, 0, EPSS_red, MPFR_RNDN); printf("\n");
			converged = test_root(
				defl, j-mul, &x, wp2, EPSS_red, &err_tmp,
				&b_tmp, &d_tmp, &f_tmp, &abx_tmp, &tmpfr_tmp
			);
			// printf("converged = %d\n", converged);
			// printf("reduced roundoff error:\n");
			// mpfr_out_str(stdout, 10, 0, err_tmp, MPFR_RNDN); printf("\n");
			// printf("poly at root:\n");
			// mpfr_out_str(stdout, 10, 0, tmpfr_tmp, MPFR_RNDN); printf("\n");
			if (!converged) {
				break;
			}
			mul++;
		}
		// if (j == 5) {
		// 	mul = 2;
		// }
		// find roots with multiplicity higher than one
		if (mul > 1) {
			// printf("look for root with multiplicity %d\n", mul);
			// mpc_set_d_d(x, 0, 0, MPFR_RNDN);
			converged = mp_laguer_impr(ad, j, &x, &its, wp2, EPSS_large, &err[j], mul);
			// printf("n. iterations = %d\n", its);
		}

		// // for (mul=j; mul>j-1; mul--) {
		// for (mul=2; mul>0; mul--) {
		// 	printf("try mul = %d\n", mul);
		// 	mpc_set_d_d(x, 0, 0, MPFR_RNDN);
		// 	converged = mp_laguer_impr(ad, j, &x, &its, wp2, EPSS, &err[j], mul);
		// 	printf("n. iterations = %d\n", its);
		// 	if (converged) {
		// 		break;
		// 	}
		// }
		
		mpfr_add(err[j], err[j], prev_err, MPFR_RNDN);
		mpfr_set(prev_err, err[j], MPFR_RNDN);
		// printf("output of lauguer:\n");
		// printf("x = ");
		// mpc_out_str(stdout, 10, 0, x, MPFR_RNDN);
		// printf("\n");
		// printf("\n");
		mpfr_abs(abs_re, mpc_realref(x), MPFR_RNDN);
		mpfr_abs(abs_im, mpc_imagref(x), MPFR_RNDN);
		mpfr_mul(abs_re, abs_re, EPS, MPFR_RNDN);
		mpfr_mul_ui(abs_re, abs_re, 2, MPFR_RNDN);
		if (mpfr_lessequal_p(abs_im, abs_re)) {
			mpc_real(tmpfr, x, MPFR_RNDN);
			mpc_set_fr(x, tmpfr, MPFR_RNDN);
			// mpfr_set_d(mpc_imagref(x), 0, MPFR_RNDN);
		}

		for (int k=j; k>j-mul; k--) {
			mpc_set(roots[k], x, MPFR_RNDN);
			// mpc_out_str(stdout, 10, 0, roots[j], MPFR_RNDN); printf("\n");
			// forward deflation
			mpc_set(b, ad[k], MPFR_RNDN);
			for (jj=k-1; jj>=0; jj--) {
				mpc_set(c, ad[jj], MPFR_RNDN);
				mpc_set(ad[jj], b, MPFR_RNDN);
				mpc_fma(b, x, b, c, MPFR_RNDN);
				// mpc_out_str(stdout, 10, 0, ad[jj], MPFR_RNDN); printf("\n");
			}
		}
		j = j - mul + 1;
	}

	if (polish) {
		mpfr_set_d(prev_err, 0, MPFR_RNDN);
		for (j=1; j<=m; j++) {
			mp_laguer(a, m, &roots[j], &its, wp2, EPSS, &err[j]);
			mpfr_add(err[j], err[j], prev_err, MPFR_RNDN);
			mpfr_set(prev_err, err[j], MPFR_RNDN);
		}
	}
	
	for (j=2; j<=m; j++) {
		mpc_set(x, roots[j], MPFR_RNDN);
		mpfr_set(tmpfr, err[j], MPFR_RNDN);
		for (i=j-1; i>=1; i--) {
			if (mpfr_lessequal_p(mpc_realref(roots[i]), mpc_realref(x))) {
				break;
			}
			mpc_set(roots[i+1], roots[i], MPFR_RNDN);
			mpfr_set(err[i+1], err[i], MPFR_RNDN);
		}
		mpc_set(roots[i+1], x, MPFR_RNDN);
		mpfr_set(err[i+1], tmpfr, MPFR_RNDN);
	}

	// empirical
	mpfr_set_ui(tmpfr, 10, MPFR_RNDN);
	mpfr_pow_si(tmpfr, tmpfr, wp2*3/10/10, MPFR_RNDN);
	for (j=1; j<=m; j++) {
		mpfr_mul(err[j], err[j], tmpfr, MPFR_RNDN);
	}
}


int mp_laguer_impr2(mpc_t *a, int m, mpc_t *x, int *its, int wp2, mpfr_t EPSS, mpfr_t *err, int mul)
{
	int iter, j;
	mpfr_t abx, abp, abm;
	mpfr_init2(abx, wp2);
	mpfr_init2(abp, wp2);
	mpfr_init2(abm, wp2);
	mpc_t dx, x1, b, d, f, g, h, sq, gp, gm, g2;
	mpc_init3(dx, wp2, wp2);
	mpc_init3(x1, wp2, wp2);
	mpc_init3(b, wp2, wp2);
	mpc_init3(d, wp2, wp2);
	mpc_init3(f, wp2, wp2);
	mpc_init3(g, wp2, wp2);
	mpc_init3(h, wp2, wp2);
	mpc_init3(sq, wp2, wp2);
	mpc_init3(gp, wp2, wp2);
	mpc_init3(gm, wp2, wp2);
	mpc_init3(g2, wp2, wp2);
	static mpfr_t frac[MR+1];
	for (int i=0; i<=MR; i++) {
		mpfr_init2(frac[i], wp2);
	}
	mpfr_set_d(frac[0], 0.0, MPFR_RNDN);
	mpfr_set_d(frac[1], 0.5, MPFR_RNDN);
	mpfr_set_d(frac[2], 0.25, MPFR_RNDN);
	mpfr_set_d(frac[3], 0.75, MPFR_RNDN);
	mpfr_set_d(frac[4], 0.13, MPFR_RNDN);
	mpfr_set_d(frac[5], 0.38, MPFR_RNDN);
	mpfr_set_d(frac[6], 0.62, MPFR_RNDN);
	mpfr_set_d(frac[7], 0.88, MPFR_RNDN);
	mpfr_set_d(frac[8], 1.0, MPFR_RNDN);
	// static double frac[MR+1] = {0.0,0.5,0.25,0.75,0.13,0.38,0.62,0.88,1.0};

	mpfr_t tmpfr;
	mpfr_init2(tmpfr, wp2);
	mpc_t tmpc;
	mpc_init3(tmpc, wp2, wp2);
	mpc_t *defl;
	defl = (mpc_t*) malloc((m+1)*sizeof(mpc_t));
	int mm, converged;
	for (int k=0; k<=m; k++) {
		mpc_init3(defl[k], wp2, wp2);
	}
	for (iter = 1; iter <= MAXIT; iter++){
		*its = iter;

		// // perform enhanced test
		for (int k=0; k<=m; k++) {
			mpc_set(defl[k], a[k], MPFR_RNDN);
		}

		for (mm=0; mm<mul; mm++) {
			printf("mm = %d\n", mm);
			converged = test_root(
				defl, m-mm, x, wp2, EPSS, err,
				&b, &d, &f, &abx, &tmpfr
			);
			if (!converged) {
				printf("roundoff error:\n");
				mpfr_out_str(stdout, 10, 0, *err, MPFR_RNDN); printf("\n");
				printf("poly at root:\n");
				mpfr_out_str(stdout, 10, 0, tmpfr, MPFR_RNDN); printf("\n");
				break;
			}
			deflation(defl, m-mm, x, wp2);
		}
		if (mm == mul) {
			printf("root:\n");
			mpc_out_str(stdout, 10, 0, *x, MPFR_RNDN); printf("\n");
			printf("polynomial at root:\n");
			mpfr_out_str(stdout, 10, 0, tmpfr, MPFR_RNDN); printf("\n");
			printf("roundoff error:\n");
			mpfr_out_str(stdout, 10, 0, *err, MPFR_RNDN); printf("\n");
			return 1;
		}
		
		// converged = test_root(
		// 	defl, m, x, wp2, EPSS, err,
		// 	&b, &d, &f, &abx, &tmpfr
		// );
		// if (converged) {
		// 	return 1;
		// }

		mpc_div(g, d, b, MPFR_RNDN);
		mpc_sqr(g2, g, MPFR_RNDN);
		mpc_mul_ui(tmpc, f, 2, MPFR_RNDN);
		mpc_div(tmpc, tmpc, b, MPFR_RNDN);
		mpc_sub(h, g2, tmpc, MPFR_RNDN);
		mpc_mul_ui(tmpc, h, m, MPFR_RNDN);
		mpc_sub(tmpc, tmpc, g2, MPFR_RNDN);
		mpc_mul_ui(tmpc, tmpc, m-1, MPFR_RNDN);
		mpc_sqrt(sq, tmpc, MPFR_RNDN);
		mpc_add(gp, g, sq, MPFR_RNDN);
		mpc_sub(gm, g, sq, MPFR_RNDN);
		mpc_abs(abp, gp, MPFR_RNDN);
		mpc_abs(abm, gm, MPFR_RNDN);
		if (mpfr_lessequal_p(abp, abm)) {
			mpc_set(gp, gm, MPFR_RNDN);
		}
		mpfr_max(tmpfr, abp, abm, MPFR_RNDN);
		if (mpfr_sgn(tmpfr) > 0) {
			mpc_pow_si(dx, gp, -1, MPFR_RNDN);
			mpc_mul_ui(dx, dx, m, MPFR_RNDN);
		} else {
			mpc_set_ui_ui(tmpc, 0, 1, MPFR_RNDN);
			mpc_mul_ui(tmpc, tmpc, iter, MPFR_RNDN);
			mpc_exp(dx, tmpc, MPFR_RNDN);
			mpfr_add_ui(tmpfr, abx, 1, MPFR_RNDN);
			mpc_mul_fr(dx, dx, tmpfr, MPFR_RNDN);
		}
		
		mpc_sub(x1, *x, dx, MPFR_RNDN);
		mpc_sub(tmpc, *x, x1, MPFR_RNDN);
		if (mpfr_equal_p(mpc_realref(*x), mpc_realref(x1)) && mpfr_equal_p(mpc_imagref(*x), mpc_imagref(x1))) {
			mpfr_set(*err, tmpfr, MPFR_RNDN);
		}
		
		if (iter % MT_red) {
			mpc_set(*x, x1, MPFR_RNDN);
		} else {
			mpc_mul_fr(tmpc, dx, frac[iter/MT_red], MPFR_RNDN);
			mpc_sub(*x, *x, tmpc, MPFR_RNDN);
		} 
	}
	perror("too many iterationc in mp_laguer");
	return 0;
}


void mp_zroots_impr(
	// OUTPUT
	mpc_t *roots, mpfr_t *err,
	// INPUT
	mpc_t *a, int m, mpc_t *dens,
	int polish, int wp2, mpfr_t mpfr_tol, int curr_mul
) {
	// increased precision variables
	int p_wp2;
	mpc_t *p_a = NULL, *p_dens = NULL, *p_roots = NULL;
	mpfr_t *p_tols = NULL;
	mpfr_t p_mpfr_tol;
	mpfr_init(p_mpfr_tol);
	mpc_t p_x;
	mpc_init3(p_x, wp2, wp2);
	mpc_t tmpc;
	mpc_init3(tmpc, wp2, wp2);
	mpfr_t p_err, p_EPSS;
	mpfr_init(p_err);
	mpfr_init(p_EPSS);
	mpfr_t tmp_mpfr;
	mpfr_init(tmp_mpfr);
	

	// mpfr_set_ui(mpfr_one, 1, MPFR_RNDN);
	mpfr_t EPSS;
	mpfr_init2(EPSS, wp2);
	mpfr_mul_ui(EPSS, mpfr_tol, 10, MPFR_RNDN);
	mpfr_t EPS;
	mpfr_init2(EPS, wp2);
	mpfr_mul_ui(EPS, EPSS, 20, MPFR_RNDN);

	mpfr_t EPSS_large;
	mpfr_init2(EPSS_large, wp2);
	mpfr_set_ui(EPSS_large, 10, MPFR_RNDN);
	mpfr_pow_si(EPSS_large, EPSS_large, wp2*3/10/20, MPFR_RNDN);
	mpfr_mul(EPSS_large, EPSS, EPSS_large, MPFR_RNDN);

	mpfr_t tmpfr;
	mpfr_init2(tmpfr, wp2);
	mpfr_t prev_err;
	mpfr_init2(prev_err, wp2);

	int i, its, j, jj, k;
	mpc_t x, b, c;
	mpc_init3(x, wp2, wp2);
	mpc_init3(b, wp2, wp2);
	mpc_init3(c, wp2, wp2);
	mpc_t ad[m+1]; // = (mpc_t*) malloc((m+1)*sizeof(mpc_t));
	for (j=0; j<=m; j++){
		mpc_init3(ad[j], wp2, wp2);
		mpc_set(ad[j], a[j], MPFR_RNDN);
	}

	mpfr_t abs_re, abs_im;
	mpfr_init2(abs_re, wp2);
	mpfr_init2(abs_im, wp2);
	mpfr_set_d(prev_err, 0, MPFR_RNDN);
	int mul, converged;
	mpc_t *defl = (mpc_t*) malloc((m+1)*sizeof(mpc_t));	
	for (k=0; k<=m; k++) {
		mpc_init3(defl[k], wp2, wp2);
	}
	mpfr_t err_tmp, tmpfr_tmp, abx_tmp, EPSS_red, EPSS_red_tmp;
	mpfr_init2(err_tmp, wp2);
	mpfr_init2(tmpfr_tmp, wp2);
	mpfr_init2(abx_tmp, wp2);
	mpfr_init2(EPSS_red, wp2);
	mpfr_init2(EPSS_red_tmp, wp2);
	mpc_t b_tmp, d_tmp, f_tmp;
	mpc_init3(b_tmp, wp2, wp2);
	mpc_init3(d_tmp, wp2, wp2);
	mpc_init3(f_tmp, wp2, wp2);

	// divide by denominators if necessary
	if (dens) {
		// printf("divide by denominators\n");
		for (k=0; k<=m; k++) {
			mpc_div(ad[k], ad[k], dens[k], MPFR_RNDN);
		}
	}

	for (j=m; j>=1; j--) {
		// printf("\n");
		// printf("j = %d\n", j);
		// printf("polyomial ad\n");
		// for (int k=0; k<=j; k++) {
		// 	mpc_out_str(stdout, 10, 0, ad[k], MPFR_RNDN);
		// 	printf("\n");
		// }
		// printf("\n");

		// if (j == 5) {
		// 	for (int k=0; k<=j; k++) {
		// 		mpc_out_str(stdout, 10, 0, ad[k], MPFR_RNDN); printf("\n");
		// 	}
		// }
		mpc_set_d_d(x, 0, 0, MPFR_RNDN);
		mp_laguer(ad, j, &x, &its, wp2, EPSS, &err[j]);
		// printf("root:\n");
		// mpc_out_str(stdout, 10, 0, x, MPFR_RNDN); printf("\n");

		for (k=0; k<=j; k++) {
			mpc_set(defl[k], ad[k], MPFR_RNDN);
		}

		// DETECT MULTIPLICITY
		mul=1;
		while (1) {
			if (mul==j) {
				break;
			}
			// printf("try mul = %d\n", mul);
			deflation(defl, j-mul+1, &x, wp2);
			// if (j == 5) {
			// 	for (int k=0; k<j; k++) {
			// 		mpc_out_str(stdout, 10, 0, defl[k], MPFR_RNDN); printf("\n");
			// 	}
			// }
			mpfr_rootn_ui(EPSS_red, EPSS, mul+1, MPFR_RNDN);
			mpfr_mul_d(EPSS_red, EPSS_red, 1e10, MPFR_RNDN);
			// printf("reduced precision:\n");
			// mpfr_out_str(stdout, 10, 0, EPSS_red, MPFR_RNDN); printf("\n");
			converged = test_root(
				defl, j-mul, &x, wp2, EPSS_red, &err_tmp,
				&b_tmp, &d_tmp, &f_tmp, &abx_tmp, &tmpfr_tmp
			);
			// printf("converged = %d\n", converged);
			// printf("reduced roundoff error:\n");
			// mpfr_out_str(stdout, 10, 0, err_tmp, MPFR_RNDN); printf("\n");
			// printf("poly at root:\n");
			// mpfr_out_str(stdout, 10, 0, tmpfr_tmp, MPFR_RNDN); printf("\n");
			if (!converged) {
				break;
			}
			mul++;
		}
		// if (j == 5) {
		// 	mul = 2;
		// }
		// find roots with multiplicity higher than one
		mpfr_rootn_ui(EPSS_red, EPSS, mul, MPFR_RNDN);
		if (mul > curr_mul) {
			// printf("look for root with multiplicity %d\n", mul);
			// temporarily increase precision
			p_wp2 = mul/curr_mul * wp2;
			if (p_a) {
				free(p_a);
			}
			p_a = (mpc_t*) malloc((m+1)*sizeof(mpc_t));
			for (k=0; k<=m; k++) {
				// printf("set increased precision\n");
				mpc_init3(p_a[k], p_wp2, p_wp2);
				// printf("set coefficient\n");
				mpc_set(p_a[k], a[k], MPFR_RNDN);
				// mpc_out_str(stdout, 10, 0, p_a[k], MPFR_RNDN); printf("\n");
			}
			if (p_dens) {
				free(p_dens);
				p_dens = (mpc_t*) malloc((m+1)*sizeof(mpc_t));
				for (k=0; k<=m; k++) {
					// printf("set increased precision\n");
					mpc_init3(p_dens[k], p_wp2, p_wp2);
					// printf("set coefficient\n");
					mpc_set(p_dens[k], dens[k], MPFR_RNDN);
					// mpc_out_str(stdout, 10, 0, p_a[k], MPFR_RNDN); printf("\n");
				}
			}
			if (p_roots) {
				free(p_roots);
			}
			p_roots = (mpc_t*) malloc((m+1)*sizeof(mpc_t));
			if (p_tols) {
				free(p_tols);
			}
			p_tols = (mpfr_t*) malloc((m+1)*sizeof(mpfr_t));
			for (k=0; k<=m; k++) {
				mpc_init3(p_roots[k], p_wp2, p_wp2);
				mpfr_init2(p_tols[k], p_wp2);
			}
			mpfr_set_prec(p_mpfr_tol, p_wp2);
			mpfr_rootn_ui(p_mpfr_tol, mpfr_tol, curr_mul, MPFR_RNDN);
			mpfr_pow_si(p_mpfr_tol, p_mpfr_tol, mul, MPFR_RNDN);
			// mpc_set_prec(p_x, p_wp2);
			// mpc_set_d_d(p_x, 0, 0, MPFR_RNDN);
			// mpc_out_str(stdout, 10, 0, p_x, MPFR_RNDN); printf("\n");
			// mpfr_set_prec(p_err, p_wp2);
			// mpfr_set_prec(p_EPSS, p_wp2);
			// mpfr_mul_ui(p_EPSS, mpfr_tol, 10, MPFR_RNDN);
			// mpfr_out_str(stdout, 10, 0, p_EPSS, MPFR_RNDN); printf("\n");
			// mpfr_pow_ui(p_EPSS, p_EPSS, mul, MPFR_RNDN);
			// mpfr_out_str(stdout, 10, 0, p_EPSS, MPFR_RNDN); printf("\n");
			// mp_laguer(p_a, j, &p_x, &its, p_wp2, p_EPSS, &p_err);
			// printf("n. iterations = %d\n", its);
			// printf("root:\n");
			// mpc_out_str(stdout, 10, 0, p_x, MPFR_RNDN); printf("\n");
			mp_zroots_impr(
				p_roots, p_tols,
				p_a, m, dens,
				1, p_wp2, p_mpfr_tol, mul);
			mpfr_set_prec(tmp_mpfr, p_wp2);
			mpc_set_prec(tmpc, p_wp2);
			// printf("augmented roots:\n");
			for (k=0; k<=m; k++) {
				mpc_set(roots[k], p_roots[k], MPFR_RNDN);
				mpc_sub(tmpc, roots[k], x, MPFR_RNDN);
				mpc_abs(tmp_mpfr, tmpc, MPFR_RNDN);
				// mpfr_sub(tmp_mpfr, tmp_mpfr, mpfr_tol, MPFR_RNDN);
				// mpfr_abs(tmp_mpfr, tmp_mpfr, MPFR_RNDN);
				if (mpfr_less_p(tmp_mpfr, EPSS_red)) {
					mpfr_rootn_ui(err[k], p_tols[k], mul, MPFR_RNDN);
				} else {
					mpfr_set(err[k], p_tols[k], MPFR_RNDN);
				}
			}
			return;
		}

		// // for (mul=j; mul>j-1; mul--) {
		// for (mul=2; mul>0; mul--) {
		// 	printf("try mul = %d\n", mul);
		// 	mpc_set_d_d(x, 0, 0, MPFR_RNDN);
		// 	converged = mp_laguer_impr(ad, j, &x, &its, wp2, EPSS, &err[j], mul);
		// 	printf("n. iterations = %d\n", its);
		// 	if (converged) {
		// 		break;
		// 	}
		// }
		
		mpfr_add(err[j], err[j], prev_err, MPFR_RNDN);
		mpfr_set(prev_err, err[j], MPFR_RNDN);
		// printf("output of lauguer:\n");
		// printf("x = ");
		// mpc_out_str(stdout, 10, 0, x, MPFR_RNDN);
		// printf("\n");
		// printf("\n");
		mpfr_abs(abs_re, mpc_realref(x), MPFR_RNDN);
		mpfr_abs(abs_im, mpc_imagref(x), MPFR_RNDN);
		mpfr_mul(abs_re, abs_re, EPS, MPFR_RNDN);
		mpfr_mul_ui(abs_re, abs_re, 2, MPFR_RNDN);
		if (mpfr_lessequal_p(abs_im, abs_re)) {
			mpc_real(tmpfr, x, MPFR_RNDN);
			mpc_set_fr(x, tmpfr, MPFR_RNDN);
			// mpfr_set_d(mpc_imagref(x), 0, MPFR_RNDN);
		}

		for (int k=j; k>j-mul; k--) {
			mpc_set(roots[k], x, MPFR_RNDN);
			// mpc_out_str(stdout, 10, 0, roots[j], MPFR_RNDN); printf("\n");
			// forward deflation
			mpc_set(b, ad[k], MPFR_RNDN);
			for (jj=k-1; jj>=0; jj--) {
				mpc_set(c, ad[jj], MPFR_RNDN);
				mpc_set(ad[jj], b, MPFR_RNDN);
				mpc_fma(b, x, b, c, MPFR_RNDN);
				// mpc_out_str(stdout, 10, 0, ad[jj], MPFR_RNDN); printf("\n");
			}
		}
		j = j - mul + 1;
	}

	if (polish) {
		mpfr_set_d(prev_err, 0, MPFR_RNDN);
		for (j=1; j<=m; j++) {
			mp_laguer(a, m, &roots[j], &its, wp2, EPSS, &err[j]);
			mpfr_add(err[j], err[j], prev_err, MPFR_RNDN);
			mpfr_set(prev_err, err[j], MPFR_RNDN);
		}
	}
	
	for (j=2; j<=m; j++) {
		mpc_set(x, roots[j], MPFR_RNDN);
		mpfr_set(tmpfr, err[j], MPFR_RNDN);
		for (i=j-1; i>=1; i--) {
			if (mpfr_lessequal_p(mpc_realref(roots[i]), mpc_realref(x))) {
				break;
			}
			mpc_set(roots[i+1], roots[i], MPFR_RNDN);
			mpfr_set(err[i+1], err[i], MPFR_RNDN);
		}
		mpc_set(roots[i+1], x, MPFR_RNDN);
		mpfr_set(err[i+1], tmpfr, MPFR_RNDN);
	}

	// empirical
	mpfr_set_ui(tmpfr, 10, MPFR_RNDN);
	mpfr_pow_si(tmpfr, tmpfr, wp2*3/10/10, MPFR_RNDN);
	for (j=1; j<=m; j++) {
		mpfr_mul(err[j], err[j], tmpfr, MPFR_RNDN);
	}
}


void mp_zroots_impr_out_mul(
	// OUTPUT
	mpc_t *roots, mpfr_t *err, int *out_mul,
	// INPUT
	mpc_t *a, int m, mpc_t *dens,
	int polish, int wp2, mpfr_t mpfr_tol, int curr_mul
) {	
	// mpfr_set_ui(mpfr_one, 1, MPFR_RNDN);
	mpfr_t EPSS;
	mpfr_init2(EPSS, wp2);
	mpfr_mul_ui(EPSS, mpfr_tol, 10, MPFR_RNDN);
	mpfr_t EPS;
	mpfr_init2(EPS, wp2);
	mpfr_mul_ui(EPS, EPSS, 20, MPFR_RNDN);

	mpfr_t EPSS_rel_increase;
	mpfr_init2(EPSS_rel_increase, wp2);
	mpfr_set_ui(EPSS_rel_increase, 10, MPFR_RNDN);
	mpfr_pow_si(EPSS_rel_increase, EPSS_rel_increase, wp2*3/10/20, MPFR_RNDN);

	mpfr_t tmpfr;
	mpfr_init2(tmpfr, wp2);
	mpfr_t prev_err;
	mpfr_init2(prev_err, wp2);

	int i, its, j, jj, k;
	mpc_t x, b, c;
	mpc_init3(x, wp2, wp2);
	mpc_init3(b, wp2, wp2);
	mpc_init3(c, wp2, wp2);
	mpc_t ad[m+1]; // = (mpc_t*) malloc((m+1)*sizeof(mpc_t));
	for (j=0; j<=m; j++) {
		mpc_init3(ad[j], wp2, wp2);
		mpc_set(ad[j], a[j], MPFR_RNDN);
	}

	mpfr_t abs_re, abs_im;
	mpfr_init2(abs_re, wp2);
	mpfr_init2(abs_im, wp2);
	mpfr_set_d(prev_err, 0, MPFR_RNDN);
	int mul, converged;
	mpc_t *defl = (mpc_t*) malloc((m+1)*sizeof(mpc_t));	
	for (k=0; k<=m; k++) {
		mpc_init3(defl[k], wp2, wp2);
	}
	mpfr_t err_tmp, tmpfr_tmp, abx_tmp, EPSS_red, EPSS_red_tmp;
	mpfr_init2(err_tmp, wp2);
	mpfr_init2(tmpfr_tmp, wp2);
	mpfr_init2(abx_tmp, wp2);
	mpfr_init2(EPSS_red, wp2);
	mpfr_init2(EPSS_red_tmp, wp2);
	mpc_t b_tmp, d_tmp, f_tmp;
	mpc_init3(b_tmp, wp2, wp2);
	mpc_init3(d_tmp, wp2, wp2);
	mpc_init3(f_tmp, wp2, wp2);

	// divide by denominators if necessary
	if (dens) {
		// printf("divide by denominators\n");
		for (k=0; k<=m; k++) {
			mpc_div(ad[k], ad[k], dens[k], MPFR_RNDN);
		}
	}

	for (j=m; j>=1; j--) {
		// printf("\n");
		// printf("j = %d\n", j);
		// printf("polyomial ad\n");
		// for (int k=0; k<=j; k++) {
		// 	mpc_out_str(stdout, 10, 0, ad[k], MPFR_RNDN);
		// 	printf("\n");
		// }
		// printf("\n");

		// if (j == 5) {
		// 	for (int k=0; k<=j; k++) {
		// 		mpc_out_str(stdout, 10, 0, ad[k], MPFR_RNDN); printf("\n");
		// 	}
		// }
		mpc_set_d_d(x, 0, 0, MPFR_RNDN);
		mp_laguer(ad, j, &x, &its, wp2, EPSS, &err[j]);
		// printf("root:\n");
		// mpc_out_str(stdout, 10, 0, x, MPFR_RNDN); printf("\n");
		if (polish) {
			// // store unpolished root
			// mpc_set(unpol_root, x, MPFR_RNDN);
			// mpfr_set(unpol_tol, err[j], MPFR_RNDN);
			
			// find polished root
			mp_laguer(a, m, &x, &its, wp2, EPSS, &err[j]);
			// printf("polished root:\n");
			// mpc_out_str(stdout, 10, 0, x, MPFR_RNDN); printf("\n");
			// printf("polished tol:\n");
			// mpfr_out_str(stdout, 10, 0, err[j], MPFR_RNDN); printf("\n");
		}

		for (k=0; k<=j; k++) {
			mpc_set(defl[k], ad[k], MPFR_RNDN);
		}

		// DETECT MULTIPLICITY
		mul=1;
		while (1) {
			if (mul==j) {
				break;
			}
			// printf("try mul = %d\n", mul);
			deflation(defl, j-mul+1, &x, wp2);
			// if (j == 5) {
			// 	for (int k=0; k<j; k++) {
			// 		mpc_out_str(stdout, 10, 0, defl[k], MPFR_RNDN); printf("\n");
			// 	}
			// }
			mpfr_rootn_ui(EPSS_red, EPSS, mul+1, MPFR_RNDN);
			mpfr_mul(EPSS_red, EPSS_red, EPSS_rel_increase, MPFR_RNDN);
			// printf("reduced precision:\n");
			// mpfr_out_str(stdout, 10, 0, EPSS_red, MPFR_RNDN); printf("\n");
			converged = test_root(
				defl, j-mul, &x, wp2, EPSS_red, &err_tmp,
				&b_tmp, &d_tmp, &f_tmp, &abx_tmp, &tmpfr_tmp
			);
			// printf("converged = %d\n", converged);
			// printf("reduced roundoff error:\n");
			// mpfr_out_str(stdout, 10, 0, err_tmp, MPFR_RNDN); printf("\n");
			// printf("poly at root:\n");
			// mpfr_out_str(stdout, 10, 0, tmpfr_tmp, MPFR_RNDN); printf("\n");
			if (!converged) {
				break;
			}
			mul++;
		}
		// printf("detected mul = %d\n", mul);
		*out_mul = mul;
		if (mul > curr_mul) {
			*out_mul = mul;
			// mpc_init3(roots[0], wp2, wp2);
			mpc_set(roots[0], x, MPFR_RNDN);
			goto goto_return;
			return;
		}

		// // for (mul=j; mul>j-1; mul--) {
		// for (mul=2; mul>0; mul--) {
		// 	printf("try mul = %d\n", mul);
		// 	mpc_set_d_d(x, 0, 0, MPFR_RNDN);
		// 	converged = mp_laguer_impr(ad, j, &x, &its, wp2, EPSS, &err[j], mul);
		// 	printf("n. iterations = %d\n", its);
		// 	if (converged) {
		// 		break;
		// 	}
		// }
		
		// lower precision of multiple roots
		if (mul > 1) {
			mpfr_rootn_ui(err[j], err[j], mul, MPFR_RNDN);
		}
		mpfr_add(err[j], err[j], prev_err, MPFR_RNDN);
		mpfr_set(prev_err, err[j], MPFR_RNDN);

		// // check that the polished root is compatible with the unpolished one
		// mpfr_sub(unpol_root, unpol_root, x, MPFR_RNDN);
		// if(mpfr_cmpabs(unpol_root, unpol_tol) > 0) {
		// 	fprintf(stderr, "polished root is incompatible with unpolished one\n");
		// 	exit(1);
		// };

		// printf("cumulated tol:\n");
		// mpfr_out_str(stdout, 10, 0, err[j], MPFR_RNDN); printf("\n");
		// printf("output of lauguer:\n");
		// printf("x = ");
		// mpc_out_str(stdout, 10, 0, x, MPFR_RNDN);
		// printf("\n");
		// printf("\n");
		mpfr_abs(abs_re, mpc_realref(x), MPFR_RNDN);
		mpfr_abs(abs_im, mpc_imagref(x), MPFR_RNDN);
		mpfr_mul(abs_re, abs_re, EPS, MPFR_RNDN);
		mpfr_mul_ui(abs_re, abs_re, 2, MPFR_RNDN);
		if (mpfr_lessequal_p(abs_im, abs_re)) {
			mpc_real(tmpfr, x, MPFR_RNDN);
			mpc_set_fr(x, tmpfr, MPFR_RNDN);
			// mpfr_set_d(mpc_imagref(x), 0, MPFR_RNDN);
		}

		for (int k=j; k>j-mul; k--) {
			mpc_set(roots[k], x, MPFR_RNDN);
			mpfr_set(err[k], err[j], MPFR_RNDN);
			// mpc_out_str(stdout, 10, 0, roots[j], MPFR_RNDN); printf("\n");
			// forward deflation
			mpc_set(b, ad[k], MPFR_RNDN);
			for (jj=k-1; jj>=0; jj--) {
				mpc_set(c, ad[jj], MPFR_RNDN);
				mpc_set(ad[jj], b, MPFR_RNDN);
				mpc_fma(b, x, b, c, MPFR_RNDN);
				// mpc_out_str(stdout, 10, 0, ad[jj], MPFR_RNDN); printf("\n");
			}
		}
		j = j - mul + 1;
	}

	// if (polish) {
	// 	mpfr_set_d(prev_err, 0, MPFR_RNDN);
	// 	for (j=1; j<=m; j++) {
	// 		mp_laguer(a, m, &roots[j], &its, wp2, EPSS, &err[j]);
	// 		printf("polished root:\n");
	// 		mpc_out_str(stdout, 10, 0, roots[j], MPFR_RNDN); printf("\n");
	// 		mpfr_add(err[j], err[j], prev_err, MPFR_RNDN);
	// 		mpfr_set(prev_err, err[j], MPFR_RNDN);
	// 	}
	// }
	
	// sort roots
	for (j=2; j<=m; j++) {
		mpc_set(x, roots[j], MPFR_RNDN);
		mpfr_set(tmpfr, err[j], MPFR_RNDN);
		for (i=j-1; i>=1; i--) {
			if (mpfr_lessequal_p(mpc_realref(roots[i]), mpc_realref(x))) {
				break;
			}
			mpc_set(roots[i+1], roots[i], MPFR_RNDN);
			mpfr_set(err[i+1], err[i], MPFR_RNDN);
		}
		mpc_set(roots[i+1], x, MPFR_RNDN);
		mpfr_set(err[i+1], tmpfr, MPFR_RNDN);
	}

	// empirical
	mpfr_set_ui(tmpfr, 10, MPFR_RNDN);
	mpfr_pow_si(tmpfr, tmpfr, wp2*3/10/10, MPFR_RNDN);
	// printf("final tols:\n");
	for (j=1; j<=m; j++) {
		mpfr_mul(err[j], err[j], tmpfr, MPFR_RNDN);
		// mpfr_out_str(stdout, 10, 0, err[j], MPFR_RNDN); printf("\n");
	}

	goto_return:
	// FREE
	for (j=0; j<=m; j++) {
		mpc_clear(ad[j]);
		mpc_clear(defl[j]);
	}
	free(defl);

	mpfr_clear(EPSS);
	mpfr_clear(EPS);
	mpfr_clear(EPSS_rel_increase);
	mpfr_clear(tmpfr);
	mpfr_clear(prev_err);

	mpc_clear(x);
	mpc_clear(b);
	mpc_clear(c);

	mpfr_clear(abs_re);
	mpfr_clear(abs_im);

	mpfr_clear(err_tmp);
	mpfr_clear(tmpfr_tmp);
	mpfr_clear(abx_tmp);
	mpfr_clear(EPSS_red);
	mpfr_clear(EPSS_red_tmp);
	mpc_clear(b_tmp);
	mpc_clear(d_tmp);
	mpc_clear(f_tmp);

}

