#include <stdio.h>
// #include <stdlib.h>
// #include <string.h>
// #include <errno.h>
#include "mpc.h"
// #include "rel_err_mpc.h"
// #include "cpoly.h"
#include "global_vars.h"


void print_poly(mpc_t *pol, int deg) {
    int i;
    for (i=0; i<=deg; i++) {
        printf("pow = %d:  ", i);
        mpc_out_str(stdout, 10, 0, pol[i], MPFR_RNDN);
        printf("\n");
    }
}


int ipow(int x, int power) {
    int p;
    if (x == 1)
        return 1;
    else if (x == -1) {
        if (power % 2 == 0)
            return 1;
        else if (power % 2 == 1)
            return -1;
    }
    else { 
        for (p = 0; p<power-1; p++)
            x *= x;
    }
    return x;
}


void copy_poly(mpc_t *out, mpc_t *in, int deg) {
    for (int k=0; k<=deg; k++){
        mpc_set(out[k], in[k], MPFR_RNDN);
    }
}


void init_poly(mpc_t *pol, int deg) {
    for (int k=0; k<=deg; k++){
        mpc_init3(pol[k], wp2, wp2);
    }
}


void set_null_poly(mpc_t *pol, int deg) {
    for (int k=0; k<=deg; k++){
        mpc_init3(pol[k], wp2, wp2);
        mpc_set_d(pol[k], 0, MPFR_RNDN);
    }
}


void add_poly(mpc_t *out_pol, mpc_t *pol1, mpc_t *pol2, int deg1, int deg2) {
    int k=0;
    int k1, k2;
    int pol1_longest = 1;
    if (deg1 <= deg2) {
        k1 = deg1;
        k2 = deg2;
        pol1_longest = 0;
    } else {
        k1 = deg2;
        k2 = deg1;
    }
    for (k=0; k<=k1; k++) {
        // mpc_init3(out_pol[k], wp2, wp2);
        mpc_add(out_pol[k], pol1[k], pol2[k], MPFR_RNDN);
    }
    if (pol1_longest) {
        for (k=k1+1; k<=k2; k++) {
            // mpc_init3(out_pol[k], wp2, wp2);
            mpc_set(out_pol[k], pol1[k], MPFR_RNDN);
        }
    } else {
        for (k=k1+1; k<=k2; k++) {
            // mpc_init3(out_pol[k], wp2, wp2);
            mpc_set(out_pol[k], pol2[k], MPFR_RNDN);
        }
    }
}


void mul_poly(mpc_t *out_pol, mpc_t *pol1, mpc_t *pol2, int deg1, int deg2) {
    /*
    Perform multiplication of two polynomials:

        pol1 = \sum_{i=0}^{deg1} a_i * x^i
        pol2 = \sum_{j=0}^{deg2} b_j * x^j

        pol1*pol2 = \sum_{k=0}^{deg1+deg2} \sum_{i=\max{0, k-deg2}}^{k} a_i * b_{k-i}

        Note that deg1 has to be greater than deg2.
    */
    int i, k, i_min=0;
    mpc_t tmp;
    mpc_init3(tmp, wp2, wp2);
    for (k=0; k<=deg1; k++) {
        // mpc_init3(out_pol[k], wp2, wp2);
        mpc_set_d(out_pol[k], 0, MPFR_RNDN);
        if (k > deg2)
            i_min = k - deg2;
        for (i=i_min; i<=k; i++) {
            mpc_mul(tmp, pol1[i], pol2[k-i], MPFR_RNDN);
            mpc_add(out_pol[k], out_pol[k], tmp, MPFR_RNDN);
            // mpc_fma(out_pol[k], pol1[i], pol2[k-i], out_pol[k], MPFR_RNDN);
        }
    }
}


void mul_poly_mpq(mpq_t *out_pol, mpq_t *pol1, mpq_t *pol2, int deg1, int deg2) {
    /*
    Perform multiplication of two polynomials:

        pol1 = \sum_{i=0}^{deg1} a_i * x^i
        pol2 = \sum_{j=0}^{deg2} b_j * x^j

        pol1*pol2 = \sum_{k=0}^{deg1+deg2} \sum_{i=\max{0, k-deg2}}^{k} a_i * b_{k-i}

        Note that deg1 has to be greater than deg2.
    */

    // printf("pol1:\n");
    // for (int k=0; k<=deg1; k++) {
    //     mpq_out_str(stdout, 10, pol1[k]); printf("\n");
    // }

    int i, k, i_min=0;
    mpq_t tmp;
    mpq_init(tmp);
    for (k=0; k<=deg1; k++) {
        // mpc_init3(out_pol[k], wp2, wp2);
        mpq_set_ui(out_pol[k], 0, 1);
        if (k > deg2)
            i_min = k - deg2;
        for (i=i_min; i<=k; i++) {
            // printf("i = %d\n", i);
            // mpq_out_str(stdout, 10, pol1[i]); printf("\n");
            // mpq_out_str(stdout, 10, pol2[k-i]); printf("\n");
            mpq_mul(tmp, pol1[i], pol2[k-i]);
            // mpq_out_str(stdout, 10, tmp); printf("\n");
            // mpq_out_str(stdout, 10, out_pol[k]); printf("\n");
            mpq_add(out_pol[k], out_pol[k], tmp);
            // mpq_out_str(stdout, 10, out_pol[k]); printf("\n");
        }
    }
}


void mul_poly_mpq_tot(mpq_t *out_pol, mpq_t *pol1, mpq_t *pol2, int deg1, int deg2) {    
    int i, k, i_min = 0, i_max = deg1;
    mpq_t tmp;
    mpq_init(tmp);
    for (k=0; k<=deg1; k++) {
        // mpc_init3(out_pol[k], wp2, wp2);
        mpq_set_ui(out_pol[k], 0, 1);
        if (k > deg2)
            i_min = k - deg2;
        for (i=i_min; i<=k; i++) {
            // printf("i = %d\n", i);
            // mpq_out_str(stdout, 10, pol1[i]); printf("\n");
            // mpq_out_str(stdout, 10, pol2[k-i]); printf("\n");
            mpq_mul(tmp, pol1[i], pol2[k-i]);
            // mpq_out_str(stdout, 10, tmp); printf("\n");
            // mpq_out_str(stdout, 10, out_pol[k]); printf("\n");
            mpq_add(out_pol[k], out_pol[k], tmp);
            // mpq_out_str(stdout, 10, out_pol[k]); printf("\n");
        }
    }
    for (k=deg1+1; k<=deg1+deg2; k++) {
        mpq_set_ui(out_pol[k], 0, 1);
        i_min = k - deg2;
        for (i=i_min; i<=deg1; i++) {
            mpq_mul(tmp, pol1[i], pol2[k-i]);
            mpq_add(out_pol[k], out_pol[k], tmp);
        }
    }
}


void mul_poly_mpq_wrp(mpq_t *out_pol, mpq_t *pol1, mpq_t *pol2, int deg1, int deg2) {
    if (deg1 < deg2) {
        // printf("scambio argomenti\n");
        mul_poly_mpq(out_pol, pol2, pol1, deg2, deg1);
    } else {
        mul_poly_mpq(out_pol, pol1, pol2, deg1, deg2);
    }
}


void mul_poly_mpq_tot_wrp(mpq_t *out_pol, mpq_t *pol1, mpq_t *pol2, int deg1, int deg2) {
    if (deg1 < deg2) {
        // printf("scambio argomenti\n");
        mul_poly_mpq_tot(out_pol, pol2, pol1, deg2, deg1);
    } else {
        mul_poly_mpq_tot(out_pol, pol1, pol2, deg1, deg2);
    }
}


void div_poly(mpc_t *out_pol, mpc_t *pol1, mpc_t *pol2, int deg1, int deg2) {
    /*
    Perform division of two polynomials:

        pol1 = \sum_{i=0}^{deg1} a_i * x^i
        pol2 = \sum_{j=0}^{deg2} b_j * x^j

        Note that pol2[0] has to be non zero and deg1 has to be lower than out_deg.
    */
    mpc_t tmp;
    mpc_init3(tmp, wp2, wp2);
    int k;
    // if deg2 = 0, the division is trivial and exact
    if (deg2 == 0) {
        for (k=0; k<=deg1; k++) {
            mpc_div(tmp, pol1[k], pol2[0], MPFR_RNDN);
            mpc_set(out_pol[k], tmp, MPFR_RNDD);
        }
        return;        
    }

    int i, i_min=0;
    // mpc_init3(out_pol[0], wp2, wp2);
    // printf("start division\n");
    mpc_div(out_pol[0], pol1[0], pol2[0], MPFR_RNDN);
    // printf("done k=0\n");
    for (k=1; k<=deg1; k++) {
        // printf("k = %d\n", k);
        // mpc_init3(out_pol[k], wp2, wp2);
        mpc_set(out_pol[k], pol1[k], MPFR_RNDN);
        if (k > deg2)
            i_min = k - deg2;
        for (i=i_min; i<k; i++) {
            mpc_mul(tmp, out_pol[i], pol2[k-i], MPFR_RNDN);
            mpc_sub(out_pol[k], out_pol[k], tmp, MPFR_RNDN);
        }
	    mpc_div(out_pol[k], out_pol[k], pol2[0], MPFR_RNDN);
    }
}


int div_poly_mpq(mpq_t *out_pol, mpq_t *pol1, mpq_t *pol2, int deg1, int deg2) {
    /*
    Perform division of two polynomials:

        pol1 = \sum_{i=0}^{deg1} a_i * x^i
        pol2 = \sum_{j=0}^{deg2} b_j * x^j

        Note that pol2[0] has to be non zero and deg1 has to be lower than out_det.
    */
    mpq_t tmp;
    mpq_init(tmp);
    int k;
    // if deg2 = 0, the division is trivial and exact
    if (deg2 == 0) {
        for (k=0; k<=deg1; k++) {
            mpq_div(out_pol[k], pol1[k], pol2[0]);
        }
        return 0;
    }

    int i, i_min=0;
    mpq_div(out_pol[0], pol1[0], pol2[0]);
    for (k=1; k<=deg1; k++) {
        mpq_set(out_pol[k], pol1[k]);
        if (k > deg2)
            i_min = k - deg2;
        for (i=i_min; i<k; i++) {
            mpq_mul(tmp, out_pol[i], pol2[k-i]);
            mpq_sub(out_pol[k], out_pol[k], tmp);
        }
	    mpq_div(out_pol[k], out_pol[k], pol2[0]);
    }

    return mpq_cmp_ui(out_pol[deg1-deg2+1], 0, 1);
}


void mp_pcshft(mpc_t a, mpc_t b, mpc_t *d, int n)
{
	int k,j;
	mpc_t fac;
    mpc_init3(fac, wp2, wp2);
    mpc_t cnst;
    mpc_init3(cnst, wp2, wp2);

	mpc_t tmp;
    mpc_init3(tmp, wp2, wp2);

    mpc_set_d(cnst, 2, MPFR_RNDN);
    mpc_sub(tmp, b, a, MPFR_RNDN);
    mpc_div(cnst, cnst, tmp, MPFR_RNDN);
    mpc_set(fac, cnst, MPFR_RNDN);

	for (j=1; j<n; j++) {
		mpc_mul(d[j], d[j], fac, MPFR_RNDN);
		mpc_mul(fac, fac, cnst, MPFR_RNDN);
	}
    mpc_add(cnst, b, a, MPFR_RNDN);
    mpc_div_ui(cnst, cnst, 2, MPFR_RNDN);
	for (j=0; j<=n-2; j++) {
		for (k=n-2; k>=j; k--) {
            mpc_neg(d[k], d[k], MPFR_RNDN);
            mpc_fma(d[k], cnst, d[k+1], d[k], MPFR_RNDN);
            mpc_neg(d[k], d[k], MPFR_RNDN);
        }
    }
}


void poly_shift(mpc_t cnst, mpc_t *d, int n)
{
	int k,j;
	for (j=0; j<n; j++) {
		for (k=n-1; k>=j; k--) {
            mpc_fma(d[k], cnst, d[k+1], d[k], MPFR_RNDN);
        }
    }
}


void poly_mul_root(mpc_t *coeffs, int deg, mpc_t root) {
    /*
    Perform the product between a polynomial

        \sum_{i=0}^{deg} a_i * x^i

    (whose coefficients are given as input) and the polynomial

        x - root

    */
    mpc_set(coeffs[deg+1], coeffs[deg], MPFR_RNDN);
    for (int k=deg; k>0; k--) {
        mpc_neg(coeffs[k], coeffs[k], MPFR_RNDN);
        mpc_fma(coeffs[k], coeffs[k], root, coeffs[k-1], MPFR_RNDN);
    }
    mpc_mul(coeffs[0], coeffs[0], root, MPFR_RNDN);
    mpc_neg(coeffs[0], coeffs[0], MPFR_RNDN);
}


void poly_mul_root_cp(mpc_t *out, mpc_t *coeffs, int deg, mpc_t root) {
    /*
    Perform the product between a polynomial

        \sum_{i=0}^{deg} a_i * x^i

    (whose coefficients are given as input) and the polynomial

        x - root

    */

    mpc_set(out[deg+1], coeffs[deg], MPFR_RNDN);
    for (int k=deg; k>0; k--) {
        // ERROR: coeff should not change sign
        // mpc_neg(out[k], coeffs[k], MPFR_RNDN);
        // mpc_fma(out[k], out[k], root, coeffs[k-1], MPFR_RNDN);
        mpc_neg(coeffs[k], coeffs[k], MPFR_RNDN);
        mpc_fma(out[k], coeffs[k], root, coeffs[k-1], MPFR_RNDN);
    }
    mpc_mul(out[0], coeffs[0], root, MPFR_RNDN);
    mpc_neg(out[0], out[0], MPFR_RNDN);
}
