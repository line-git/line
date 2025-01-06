#include <stdlib.h>
#include <stdio.h>
// #include <math.h>
#include "mpc.h"
#include "global_vars.h"


void mp_ludcmp(mpc_t **a, int n, int *indx, float *d)
{
	int i, imax, j, k;
	mpc_t sum;
	mpc_init3(sum, wp2, wp2);
	mpfr_t temp, big, dum;
	mpfr_init2(temp, wp2);
	mpfr_init2(big, wp2);
	mpfr_init2(dum, wp2);
	mpfr_t *vv = (mpfr_t*) malloc(n*sizeof(mpfr_t));
	--vv;
	for (i = 1; i <= n; i++) {
		mpfr_init2(vv[i], wp2);
	}

	*d = 1.0;
	for (i = 1; i <= n; i++) {
		mpfr_set_d(big, 0, MPFR_RNDN);
		for (j = 1; j <= n; j++) {
			mpc_abs(temp, a[i][j], MPFR_RNDN);
			if (mpfr_greater_p(temp, big)) {
				mpfr_set(big, temp, MPFR_RNDN);
			}
		}
		if (mpfr_zero_p(big)) {
			printf("Singular matrix in routine mp_ludcmp\n");
			return;
		}
		mpfr_pow_si(vv[i], big, -1, MPFR_RNDN);
	}
	
	
	for (j = 1; j <= n; j++) {
		for (i = 1; i < j; i++) {
			mpc_set(sum, a[i][j], MPFR_RNDN);
			for (k = 1; k < i; k++) {
				mpc_neg(sum, sum, MPFR_RNDN);
				mpc_fma(sum, a[i][k], a[k][j], sum, MPFR_RNDN);
				mpc_neg(sum, sum, MPFR_RNDN);
			}
			mpc_set(a[i][j], sum, MPFR_RNDN);
		}
		mpfr_set_d(big, 0, MPFR_RNDN);
		for (i = j; i <= n; i++) {
			mpc_set(sum, a[i][j], MPFR_RNDN);
			for (k = 1; k < j; k++) {
				mpc_neg(sum, sum, MPFR_RNDN);
				mpc_fma(sum, a[i][k], a[k][j], sum, MPFR_RNDN);
				mpc_neg(sum, sum, MPFR_RNDN);
			}
			mpc_set(a[i][j], sum, MPFR_RNDN);
			mpc_abs(temp, sum, MPFR_RNDN);
			mpfr_mul(dum, vv[i], temp, MPFR_RNDN);
			if (mpfr_greaterequal_p(dum, big)) {
				mpfr_set(big, dum, MPFR_RNDN);
				imax = i;
			}
		}
		if (j != imax) {
			for (k = 1; k <= n; k++){
				mpc_swap(a[imax][k], a[j][k]);
			}
			*d = -(*d);

			mpfr_set(vv[imax], vv[j], MPFR_RNDN);
		}
		indx[j] = imax;
		if (j != n) {
			for (i = j + 1; i <= n; i++) {
				mpc_div(a[i][j], a[i][j], a[j][j], MPFR_RNDN);
			}
		}
		// free_vector(vv, 1, n);	
	}
}


void mp_lubksb(mpc_t **a, int n, int *indx, mpc_t *b)
{	
	int i, ii = 0, ip, j;
	mpc_t sum, dum;
	mpc_init3(sum, wp2, wp2);
	mpc_init3(dum, wp2, wp2);

	for (i = 1; i <= n; i++) {
		ip = indx[i];
		mpc_set(sum, b[ip], MPFR_RNDN);
		mpc_set(b[ip], b[i], MPFR_RNDN);

		if (ii) {
			for (j = ii; j <= i - 1; j++) {
				mpc_neg(sum, sum, MPFR_RNDN);
				mpc_fma(sum, a[i][j], b[j], sum, MPFR_RNDN);
				mpc_neg(sum, sum, MPFR_RNDN);
			}
		}
		else if (mpfr_zero_p((mpc_realref(sum))) || mpfr_zero_p((mpc_imagref(sum)))) {
			ii = i;
		}
		mpc_set(b[i], sum, MPFR_RNDN);
	}

	for (i = n; i >= 1; i--) {
		mpc_set(sum, b[i], MPFR_RNDN);
		for (j = i + 1; j <= n; j++) {
			mpc_neg(sum, sum, MPFR_RNDN);
			mpc_fma(sum, a[i][j], b[j], sum, MPFR_RNDN);
			mpc_neg(sum, sum, MPFR_RNDN);
			}
		mpc_div(b[i], sum, a[i][i], MPFR_RNDN);
	}
}

//////
//
// pointers are accepted with the 0, ..., n-1 convention
// and shifted so that numerical recipes routines can work
// with the 1, ..., n convention.
//
// The ouput is returned with 0, ..., n-1 convetion (the same
// as the input).
//////
void mp_inverse(mpc_t **a_inv, mpc_t **a, int n)
{
	// shift pointers so that indices goes from 1 to n
	for (int i=0; i<n; i++) {
		a[i]--;
	}
	a--;

	int *indx = (int*) malloc(n*sizeof(int));
	indx--;
	float d;

	// find the LU decomposition and store it into a
	mp_ludcmp(a, n, indx, &d);
	
	mpc_t* col = (mpc_t*) malloc(n*sizeof(mpc_t));
	col--;
	for (int i=1; i<=n; i++) {
		mpc_init3(col[i], wp2, wp2);
	}

	// compute the inverse by solving for identity columns
	for (int j=1; j<=n; j++) {
		for (int i=1; i<=n; i++) {
			if (i != j) {
				mpc_set_d_d(col[i], 0, 0, MPFR_RNDN);
			}
		}
		mpc_set_d_d(col[j], 1, 0, MPFR_RNDN);

		// solve linear system using col as constant term and store the result in col
		mp_lubksb(a, n, indx, col);

		// col is the j-th column of the inverse
		for (int i=0; i<n; i++) {
			mpc_set(a_inv[i][j-1], col[i+1], MPFR_RNDN);
		}
	}

	// shift back row pointers
	// (because row pointers would otherwise remain shifted in the main)
	a++;
	for (int i=0; i<n; i++) {
		a[i]++;
	}
}
