#include <stdio.h>
// #include <math.h>
#include <stdbool.h>
#include <stdlib.h>

#include "mpc.h"
// #include "rel_err_mpc.h"
#include "global_vars.h"


static int c__0 = 0;

#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))


//////
// TRANSLATED SUBROUTINES
//////

int cbal_(int *nm, int *n, mpfr_t *ar, mpfr_t *ai, int *low, int *igh, mpfr_t *scale)
{
	
	mpfr_t mpfr_zero, mpfr_one;
	mpfr_init2(mpfr_zero, wp2);
	mpfr_set_d(mpfr_zero, (double) 0, MPFR_RNDN);
	mpfr_init2(mpfr_one, wp2);
	mpfr_set_d(mpfr_one, (double) 1, MPFR_RNDN);
    /* System generated locals */
    int ar_dim1, ar_offset, ai_dim1, ai_offset, i__1, i__2;
    mpfr_t r__1;
    mpfr_init2(r__1, wp2);
    mpfr_t r__2;
    mpfr_init2(r__2, wp2);

    /* Local variables */
    static mpfr_t c__;
    mpfr_init2(c__, wp2);
    static mpfr_t f;
    mpfr_init2(f, wp2);
    static mpfr_t g;
    mpfr_init2(g, wp2);
    static int i__, j, k, l, m;
    static mpfr_t r__;
    mpfr_init2(r__, wp2);
    static mpfr_t s;
    mpfr_init2(s, wp2);
    static mpfr_t b2;
    mpfr_init2(b2, wp2);
    static int jj, iexc;
    static mpfr_t radix;
    mpfr_init2(radix, wp2);
    static bool noconv;

	mpfr_t cmp0;
	mpfr_init2(cmp0, wp2);
	mpfr_t cmp1;
	mpfr_init2(cmp1, wp2);

	/*     REAL ABS */

	/*     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE */
	/*     CBALANCE, WHICH IS A COMPLEX VERSION OF BALANCE, */
	/*     NUM. MATH. 13, 293-304(1969) BY PARLETT AND REINSCH. */
	/*     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 315-326(1971). */

	/*     THIS SUBROUTINE BALANCES A COMPLEX MATRIX AND ISOLATES */
	/*     EIGENVALUES WHENEVER POSSIBLE. */

	/*     ON INPUT- */

	/*       NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL */
	/*         ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM */
	/*         DIMENSION STATEMENT, */

	/*       N IS THE ORDER OF THE MATRIX, */

	/*       AR AND AI CONTAINS THE REAL AND IMAGINARY PARTS, */
	/*         RESPECTIVELY, OF THE COMPLEX MATRIX TO BE BALANCED. */

	/*     ON OUTPUT- */

	/*       AR AND AI CONTAIN THE REAL AND IMAGINARY PARTS, */
	/*         RESPECTIVELY, OF THE BALANCED MATRIX, */

	/*         LOW AND IGH ARE TWO intS SUCH THAT AR(I,J) AND AI(I,J) */
	/*          ARE EQUAL TO ZERO IF */
	/*           (1) I IS GREATER THAN J AND */
	/*           (2) J=1,...,LOW-1 OR I-IGH+1,...,N, */

	/*         SCALE CONTAINS INFORMATION DETERMINING THE */
	/*           PERMUTATIONS AND SCALING FACTORS USED. */

	/*     SUPPOSE THAT THE PRINCIPAL SUBMATRIX IN ROWS LOW THROUGH IGH */
	/*     HAS BEEN BALANCED, THAT P(J) DENOTES THE INDEX INTERCHANGED */
	/*     WITH J DURING THE PERMUTATION STEP, AND THAT THE ELEMENTS */
	/*     OF THE DIAGONAL MATRIX USED ARE DENOTED BY D(I,J). THEN */
	/*         SCALE(J) = P(J),   FOR J = 1,...,LOW-1 */
	/*                  = D(J,J)      J = LOW,...,IGH */
	/*                  = P(J)        J = IGH+1,...,N. */
	/*     THE ORDER IN WHICH THE INTERCHANGES ARE MADE IS N TO IGH+1, */
	/*     THEN 1 TO LOW-1. */

	/*     NOTE THAT I IS RETURNED FOR IGH IF IGH IS ZERO FORMALLY. */

	/*     THE ALGOL PROCEDURE EXC CONTAINED IN CBALANCE APPEARS IN */
	/*     CBAL IN LINE. (NOTE THAT THE ALGOL ROLES OF IDENTIFIERS */
	/*     K,L HAVE BFEN REVERSED.) */

	/*     ARITHMETIC IS REAL THROUGHOUT. */

	/*     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW, */
	/*     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY */

	/*     ------------------------------------------------------------------ */

	/*     ********** RADIX IS A MACHINE DEPENDENT PARAMETER SPECIFYING */
	/*                THE BASE OF BHE MACHINE FLOATING POINT REPRESENTATION. */

	/*                *********** */
	
	/* Parameter adjustments */
    --scale;
	ai_dim1 = *nm;
	ai_offset = 1 + ai_dim1;
	ai = ai - ai_offset;
	ar_dim1 = *nm;
	ar_offset = 1 + ar_dim1;
	ar = ar - ar_offset;

    /* Function Body */
mpfr_set_d(radix, 2e0, MPFR_RNDN);

mpfr_mul(b2, radix, radix, MPFR_RNDN);
k = 1;
l = *n;
    goto L100;
/*     ********** IN-LINE PROCEDURE FOR ROW AND */
/*                COLUMN EXCHANGE ********** */
L20:;
// printf("L20:;\n");
mpfr_set_si(scale[m], j, MPFR_RNDN);
    if(j == m) {
	goto L50;
    }

    if(1 > l) {
	goto L10000;
    }
    for (int i__ = 1; i__ <= l; ++i__) {
mpfr_set(f, ar[i__ + j * ar_dim1], MPFR_RNDN);
mpfr_set(ar[i__ + j * ar_dim1], ar[i__ + m * ar_dim1], MPFR_RNDN);
mpfr_set(ar[i__ + m * ar_dim1], f, MPFR_RNDN);
mpfr_set(f, ai[i__ + j * ai_dim1], MPFR_RNDN);
mpfr_set(ai[i__ + j * ai_dim1], ai[i__ + m * ai_dim1], MPFR_RNDN);
mpfr_set(ai[i__ + m * ai_dim1], f, MPFR_RNDN);
/* L30: */
    }
L10000:;
// printf("L10000:;\n");

    if(k > *n) {
	goto L10010;
    }
    for (int i__ = k; i__ <= *n; ++i__) {
mpfr_set(f, ar[j + i__ * ar_dim1], MPFR_RNDN);
mpfr_set(ar[j + i__ * ar_dim1], ar[m + i__ * ar_dim1], MPFR_RNDN);
mpfr_set(ar[m + i__ * ar_dim1], f, MPFR_RNDN);
mpfr_set(f, ai[j + i__ * ai_dim1], MPFR_RNDN);
mpfr_set(ai[j + i__ * ai_dim1], ai[m + i__ * ai_dim1], MPFR_RNDN);
mpfr_set(ai[m + i__ * ai_dim1], f, MPFR_RNDN);
/* L40: */
    }
L10010:;
// printf("L10010:;\n");

L50:;
// printf("L50:;\n");
    if(iexc <= 1) {
	goto L80;
    }
    goto L130;
/*     ********** SEARCH FOR ROWS ISOLATING AN EIGENVALUE */
/*                AND PUSH THEM DOWN ********** */
L80:;
// printf("L80:;\n");
    if(l == 1) {
	goto L280;
    }
    --l;
/*     ********** FOR J=L STEP -1 UNTIL 1 DO -- ********** */
L100:;
// printf("L100:;\n");
    if(1 > l) {
	goto L10020;
    }
    for (jj = 1; jj <= l; ++jj) {
j = l + 1;
j = j - jj;

	if(1 > l) {
	    goto L10030;
	}
	for (int i__ = 1; i__ <= l; ++i__) {
	    if(i__ == j) {
		goto L110;
	    }
	    if(!mpfr_equal_p(ar[j + i__ * ar_dim1], mpfr_zero) || !mpfr_equal_p(ai[j + i__ * ai_dim1], mpfr_zero)) {
		goto L120;
	    }
L110:;
// printf("L110:;\n");
	    ;
	}
L10030:;
// printf("L10030:;\n");

m = l;
iexc = 1;
	goto L20;
L120:;
// printf("L120:;\n");
	;
    }
L10020:;
// printf("L10020:;\n");

    goto L140;
/*     ********** SEARCH FOR COLUMNS ISOLATING AN EIGENVALUE */
/*                AND PUSH THEM LEFT ********** */
L130:;
// printf("L130:;\n");
    ++k;

L140:;
// printf("L140:;\n");
    if(k > l) {
	goto L10040;
    }
    for (j = k; j <= l; ++j) {

	if(k > l) {
	    goto L10050;
	}
	for (int i__ = k; i__ <= l; ++i__) {
	    if(i__ == j) {
		goto L150;
	    }
	    if(!mpfr_equal_p(ar[i__ + j * ar_dim1], mpfr_zero) || !mpfr_equal_p(ai[i__ + j * ai_dim1], mpfr_zero)) {
		goto L170;
	    }
L150:;
// printf("L150:;\n");
	    ;
	}
L10050:;
// printf("L10050:;\n");

m = k;
iexc = 2;
	goto L20;
L170:;
// printf("L170:;\n");
	;
    }
L10040:;
// printf("L10040:;\n");
/*     ********** NOW BALANCE THE SUBMATRIX IN ROWS K TO L ********** */
    if(k > l) {
	goto L10060;
    }
    for (int i__ = k; i__ <= l; ++i__) {
/* L180: */
mpfr_set_d(scale[i__], 1e0, MPFR_RNDN);
    }
L10060:;
// printf("L10060:;\n");
/*     ********** ITERATIVE LOOP FOR NORM REDUCTION ********** */
L190:;
// printf("L190:;\n");
noconv = false;

    if(k > l) {
	goto L10070;
    }
    for (int i__ = k; i__ <= l; ++i__) {
mpfr_set_d(c__, 0e0, MPFR_RNDN);
mpfr_set_d(r__, 0e0, MPFR_RNDN);

	if(k > l) {
	    goto L10080;
	}
	for (j = k; j <= l; ++j) {
	    if(j == i__) {
		goto L200;
	    }
mpfr_abs(r__1, ar[j + i__ * ar_dim1], MPFR_RNDN);
mpfr_abs(r__2, ai[j + i__ * ai_dim1], MPFR_RNDN);
mpfr_add(c__, c__, r__1, MPFR_RNDN);
mpfr_add(c__, c__, r__2, MPFR_RNDN);
mpfr_abs(r__1, ar[i__ + j * ar_dim1], MPFR_RNDN);
mpfr_add(r__, r__, r__1, MPFR_RNDN);
mpfr_add(r__, r__, r__2, MPFR_RNDN);
L200:;
// printf("L200:;\n");
	    ;
	}
L10080:;
// printf("L10080:;\n");

mpfr_div(g, r__, radix, MPFR_RNDN);
mpfr_set_d(f, 1, MPFR_RNDN);
mpfr_add(s, c__, r__, MPFR_RNDN);
L210:;
// printf("L210:;\n");
	if(mpfr_greaterequal_p(c__, g)) {
	    goto L220;
	}
mpfr_mul(f, f, radix, MPFR_RNDN);
mpfr_mul(c__, c__, b2, MPFR_RNDN);
	goto L210;
L220:;
// printf("L220:;\n");
// printf("L220\n");
mpfr_mul(g, r__, radix, MPFR_RNDN);
L230:;
// printf("L230:;\n");
// printf("g230\n");
	if(mpfr_less_p(c__, g)) {
	    goto L240;
	}
mpfr_div(f, f, radix, MPFR_RNDN);
mpfr_div(c__, c__, b2, MPFR_RNDN);
	goto L230;
/*     ********** NOW BALANCE ********** */
L240:;
// printf("L240:;\n");
mpfr_add(cmp0, c__, r__, MPFR_RNDN);
mpfr_div(cmp0, cmp0, f, MPFR_RNDN);
mpfr_mul_d(cmp1, s, 95, MPFR_RNDN);
mpfr_div_d(cmp1, cmp1, 100, MPFR_RNDN);
	if(mpfr_greaterequal_p(cmp0, cmp1)) {
	    goto L270;
	}
mpfr_pow_si(g, f, -1, MPFR_RNDN);
mpfr_mul(scale[i__], scale[i__], f, MPFR_RNDN);
noconv = true;

	if(k > *n) {
	    goto L10090;
	}
	for (j = k; j <= *n; ++j) {
mpfr_mul(ar[i__ + j * ar_dim1], ar[i__ + j * ar_dim1], g, MPFR_RNDN);
mpfr_mul(ai[i__ + j * ai_dim1], ai[i__ + j * ai_dim1], g, MPFR_RNDN);
/* L250: */
	}
L10090:;
// printf("L10090:;\n");

	if(1 > l) {
	    goto L10100;
	}
	for (j = 1; j <= l; ++j) {
mpfr_mul(ar[j + i__ * ar_dim1], ar[j + i__ * ar_dim1], f, MPFR_RNDN);
mpfr_mul(ai[j + i__ * ai_dim1], ai[j + i__ * ai_dim1], f, MPFR_RNDN);
/* L260: */
	}
L10100:;
// printf("L10100:;\n");

L270:;
// printf("L270:;\n");
	;
    }
L10070:;
// printf("L10070:;\n");

        if (noconv) {
	goto L190;
    }

L280:;
// printf("L280:;\n");
*low = k;
*igh = l;
    return 0;
/*     ********** LAST CARD OF CBAL ********** */
} /* cbal_ */


int comhes_(int *nm, int *n, int *low, int *igh, mpfr_t *ar, mpfr_t *ai, int *int__)
{
	
	mpfr_t mpfr_zero, mpfr_one;
	mpfr_init2(mpfr_zero, wp2);
	mpfr_set_d(mpfr_zero, (double) 0, MPFR_RNDN);
	mpfr_init2(mpfr_one, wp2);
	mpfr_set_d(mpfr_one, (double) 1, MPFR_RNDN);
	mpfr_t tmp0;
	mpfr_init2(tmp0, wp2);
	mpfr_t tmp1;
	mpfr_init2(tmp1, wp2);
    /* System generated locals */
    int ar_dim1, ar_offset, ai_dim1, ai_offset, i__1, i__2, i__3;
    mpfr_t r__1;
    mpfr_init2(r__1, wp2);
    mpfr_t r__2;
    mpfr_init2(r__2, wp2);
    mpc_t q__1;
    mpc_init3(q__1, wp2, wp2);
    static mpc_t *equiv_1 = NULL;
    if (!equiv_1) {
    	equiv_1 = (mpc_t*) malloc(1*sizeof(mpc_t));
    	for (int i=0; i<1; i++) {
    		mpc_init3(equiv_1[i], wp2, wp2);
    	}
    }
    	static mpc_t *equiv_3 = NULL;
    	if (!equiv_3) {
    		equiv_3 = (mpc_t*) malloc(1*sizeof(mpc_t));
    		for (int i=0; i<1; i++) {
    			mpc_init3(equiv_3[i], wp2, wp2);
    		}
    	}

    /* Builtin functions */
    

    /* Local variables */
    static int i__, j, m;
    static int la;
    static int mm1, kp1, mp1;

	mpc_t *x = equiv_1;
	mpc_t *y = equiv_3;
	mpfr_t *t1 = (mpfr_t*) equiv_1;
	mpfr_t *t2 = (mpfr_t*) equiv_3;
	mpfr_t *xi = (mpfr_t*) equiv_1 + 1;
	mpfr_t *yi = (mpfr_t*) equiv_3 + 1;
	mpfr_t *xr = (mpfr_t*) equiv_1;
	mpfr_t *yr = (mpfr_t*) equiv_3;


/*     REAL ABS */

/*     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE COMHES, */
/*     NUM. MATH. 12, 349-368(1968) BY MARTIN AND WILKINSON. */
/*     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 339-358(1971). */

/*     GIVEN A COMPLEX GENERAL MATRIX, THIS SUBROUTINE */
/*     REDUCES A SUBMATRIX SITUATED IN ROWS AND COLUMNS */
/*     LOW THROUGH IGH TO UPPER HESSENBERG FORM BY */
/*     STABILIZED ELEMENTARY SIMILARITY TRANSFORMATIONS. */

/*     ON INPUT- */

/*        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL */
/*           ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM */
/*           DIMENSION STATEMENT, */

/*        N IS THE ORDER OF THE MATRIX, */

/*        LOW AND IGH ARE intS DETERMINED BY THE BALANCING */
/*           SUBROUTINE CBAL. IF CBAL HAS NOT BEEN USED, */
/*           SET LOW=1, IGH=N, */

/*        AR AND AI CONTAIN THE REAL AND IMAGINARY PARTS, */
/*           RESPECTIVELY, OF THE COMPLEX INPUT MATRIX. */

/*     ON OUTPUT- */

/*     AR AND AI CONTAIN THE REAL AND IMAGINARY PARTS, */
/*           RESPECTIVELY, OF THE HESSENBERG MATRIX. THE */
/*           MULTIPLIERS WHICH WERE USED IN THE REDUCTION */
/*           ARE STORED IN THE REMAINING TRIANGLES UNDER THE */
/*           HESSENBERG MATRIX, */

/*        INT CONTAINS INFORMATION ON THE ROWS AND COLUMNS */
/*           INTERCHANGED IN THE REDUCTION. */
/*           ONLY ELEMENTS LOW THROUGH IGH ARE USED. */

/*     ARITHMETIC IS REAL EXCEPT FOR THE REPLACEMENT OF THE ALGOL */
/*     PROCEDURE CDIV BY COMPLEX DIVISION. */

/*     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW, */
/*     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY */

/*     ------------------------------------------------------------------ */

    /* Parameter adjustments */
ai_dim1 = *nm;
ai_offset = 1 + ai_dim1;
ai = ai - ai_offset;
ar_dim1 = *nm;
ar_offset = 1 + ar_dim1;
ar = ar - ar_offset;
    --int__;

    /* Function Body */
la = *igh - 1;
kp1 = *low + 1;
    if(la < kp1) {
	goto L200;
    }

    if(kp1 > la) {
	goto L10000;
    }
    for (m = kp1; m <= la; ++m) {
mm1 = m - 1;
mpfr_set_d(*xr, 0, MPFR_RNDN);
mpfr_set_d(*xi, 0, MPFR_RNDN);
i__ = m;

	if(m > *igh) {
	    goto L10010;
	}
	for (j = m; j <= *igh; ++j) {
mpfr_abs(r__1, ar[j + mm1 * ar_dim1], MPFR_RNDN);
mpfr_abs(r__2, ai[j + mm1 * ai_dim1], MPFR_RNDN);
mpfr_add(r__1, r__1, r__2, MPFR_RNDN);
mpfr_abs(tmp0, *xr, MPFR_RNDN);
mpfr_abs(r__2, *xi, MPFR_RNDN);
mpfr_add(r__2, tmp0, r__2, MPFR_RNDN);
	    if(mpfr_lessequal_p(r__1, r__2)) {
		goto L100;
	    }
mpfr_set(*xr, ar[j + mm1 * ar_dim1], MPFR_RNDN);
mpfr_set(*xi, ai[j + mm1 * ai_dim1], MPFR_RNDN);
i__ = j;
L100:;
// printf("L100:;\n");
	    ;
	}
L10010:;
// printf("L10010:;\n");

int__[m] = i__;
	if(i__ == m) {
	    goto L130;
	}
/*     ********** INTERCHANGE ROWS AND COLUMNS OF AR AND AI ********** */
	if(mm1 > *n) {
	    goto L10020;
	}
	for (j = mm1; j <= *n; ++j) {
mpfr_set(*yr, ar[i__ + j * ar_dim1], MPFR_RNDN);
mpfr_set(ar[i__ + j * ar_dim1], ar[m + j * ar_dim1], MPFR_RNDN);
mpfr_set(ar[m + j * ar_dim1], *yr, MPFR_RNDN);
mpfr_set(*yi, ai[i__ + j * ai_dim1], MPFR_RNDN);
mpfr_set(ai[i__ + j * ai_dim1], ai[m + j * ai_dim1], MPFR_RNDN);
mpfr_set(ai[m + j * ai_dim1], *yi, MPFR_RNDN);
/* L110: */
	}
L10020:;
// printf("L10020:;\n");

	if(1 > *igh) {
	    goto L10030;
	}
	for (j = 1; j <= *igh; ++j) {
mpfr_set(*yr, ar[j + i__ * ar_dim1], MPFR_RNDN);
mpfr_set(ar[j + i__ * ar_dim1], ar[j + m * ar_dim1], MPFR_RNDN);
mpfr_set(ar[j + m * ar_dim1], *yr, MPFR_RNDN);
mpfr_set(*yi, ai[j + i__ * ai_dim1], MPFR_RNDN);
mpfr_set(ai[j + i__ * ai_dim1], ai[j + m * ai_dim1], MPFR_RNDN);
mpfr_set(ai[j + m * ai_dim1], *yi, MPFR_RNDN);
/* L120: */
	}
L10030:;
// printf("L10030:;\n");
/*     **********END INTERCHANGE ********** */
L130:;
// printf("L130:;\n");
	if(mpfr_equal_p(*xr, mpfr_zero) && mpfr_equal_p(*xi, mpfr_zero)) {
	    goto L180;
	}
mp1 = m + 1;

	if(mp1 > *igh) {
	    goto L10040;
	}
	for (i__ = mp1; i__ <= *igh; ++i__) {
mpfr_set(*yr, ar[i__ + mm1 * ar_dim1], MPFR_RNDN);
mpfr_set(*yi, ai[i__ + mm1 * ai_dim1], MPFR_RNDN);
	    if(mpfr_equal_p(*yr, mpfr_zero) && mpfr_equal_p(*yi, mpfr_zero)) {
		goto L160;
	    }
	    mpc_div(*y, *y, *x, MPFR_RNDN);
mpfr_set(ar[i__ + mm1 * ar_dim1], *yr, MPFR_RNDN);
mpfr_set(ai[i__ + mm1 * ai_dim1], *yi, MPFR_RNDN);
	    if(m > *n) {
		goto L10050;
	    }
	    for (j = m; j <= *n; ++j) {
mpfr_mul(tmp0, *yr, ar[m + j * ar_dim1], MPFR_RNDN);
mpfr_mul(tmp1, *yi, ai[m + j * ai_dim1], MPFR_RNDN);
mpfr_sub(ar[i__ + j * ar_dim1], ar[i__ + j * ar_dim1], tmp0, MPFR_RNDN);
mpfr_add(ar[i__ + j * ar_dim1], ar[i__ + j * ar_dim1], tmp1, MPFR_RNDN);
mpfr_mul(tmp0, *yr, ai[m + j * ai_dim1], MPFR_RNDN);
mpfr_mul(tmp1, *yi, ar[m + j * ar_dim1], MPFR_RNDN);
mpfr_sub(ai[i__ + j * ai_dim1], ai[i__ + j * ai_dim1], tmp0, MPFR_RNDN);
mpfr_sub(ai[i__ + j * ai_dim1], ai[i__ + j * ai_dim1], tmp1, MPFR_RNDN);
/* L140: */
	    }
L10050:;
// printf("L10050:;\n");

	    if(1 > *igh) {
		goto L10060;
	    }
i__3 = *igh;
	    for (j = 1; j <= i__3; ++j) {
mpfr_mul(tmp0, *yr, ar[j + i__ *ar_dim1], MPFR_RNDN);
mpfr_mul(tmp1, *yi, ai[j + i__ * ai_dim1], MPFR_RNDN);
mpfr_add(ar[j + m * ar_dim1], ar[j + m * ar_dim1], tmp0, MPFR_RNDN);
mpfr_sub(ar[j + m * ar_dim1], ar[j + m * ar_dim1], tmp1, MPFR_RNDN);
mpfr_mul(tmp0, *yr, ai[j + i__ *ai_dim1], MPFR_RNDN);
mpfr_mul(tmp1, *yi, ar[j + i__ * ar_dim1], MPFR_RNDN);
mpfr_add(ai[j + m * ai_dim1], ai[j + m * ai_dim1], tmp0, MPFR_RNDN);
mpfr_add(ai[j + m * ai_dim1], ai[j + m * ai_dim1], tmp1, MPFR_RNDN);
/* L150: */
	    }
L10060:;
// printf("L10060:;\n");

L160:;
// printf("L160:;\n");
	    ;
	}
L10040:;
// printf("L10040:;\n");

L180:;
// printf("L180:;\n");
	;
    }
L10000:;
// printf("L10000:;\n");

L200:;
// printf("L200:;\n");
    return 0;
/*     ********** LAST CARD OF COMHES ********** */
} /* comhes_ */


int comlr2_(int *nm, int *n, int *low, int *igh, int *int__, mpfr_t *hr, mpfr_t *hi, mpfr_t *wr, mpfr_t *wi, mpfr_t *zr, mpfr_t *zi, int *ierr, int *nobak)
{	
	// printf("low = %d, igh = %d\n", *low, *igh);
	//getchar();
	mpfr_t mpfr_zero, mpfr_one;
	mpfr_init2(mpfr_zero, wp2);
	mpfr_set_d(mpfr_zero, 0, MPFR_RNDN);
	mpfr_init2(mpfr_one, wp2);
	mpfr_set_d(mpfr_one, 1, MPFR_RNDN);
	mpfr_t tmp0;
	mpfr_init2(tmp0, wp2);
	mpfr_t tmp1;
	mpfr_init2(tmp1, wp2);
	mpfr_t cmp0;
	mpfr_init2(cmp0, wp2);
	mpfr_t cmp1;
	mpfr_init2(cmp1, wp2);
	mpfr_t cmp2;
	mpfr_init2(cmp2, wp2);
	mpfr_t cmp3;
	mpfr_init2(cmp3, wp2);

    /* System generated locals */
    int hr_dim1, hr_offset, hi_dim1, hi_offset, zr_dim1, zr_offset, zi_dim1, zi_offset, i__1, i__2, i__3;
    mpfr_t r__1;
    mpfr_init2(r__1, wp2);
    mpfr_t r__2;
    mpfr_init2(r__2, wp2);
    mpfr_t r__3;
    mpfr_init2(r__3, wp2);
    mpfr_t r__4;
    mpfr_init2(r__4, wp2);
    mpfr_t r__5;
    mpfr_init2(r__5, wp2);
    mpfr_t r__6;
    mpfr_init2(r__6, wp2);
    mpc_t q__1;
    mpc_init3(q__1, wp2, wp2);
    mpc_t q__2;
    mpc_init3(q__2, wp2, wp2);
    static mpc_t *equiv_1 = NULL;
    if (!equiv_1) {
    	equiv_1 = (mpc_t*) malloc(1*sizeof(mpc_t));
    	for (int i=0; i<1; i++) {
    		mpc_init3(equiv_1[i], wp2, wp2);
    	}
    }
    	static mpc_t *equiv_3 = NULL;
    	if (!equiv_3) {
    		equiv_3 = (mpc_t*) malloc(1*sizeof(mpc_t));
    		for (int i=0; i<1; i++) {
    			mpc_init3(equiv_3[i], wp2, wp2);
    		}
    	}
    		static mpc_t *equiv_5 = NULL;
    		if (!equiv_5) {
    			equiv_5 = (mpc_t*) malloc(1*sizeof(mpc_t));
    			for (int i=0; i<1; i++) {
    				mpc_init3(equiv_5[i], wp2, wp2);
    			}
    		}

    /* Builtin functions */

    /* Local variables */
    static int i__, j, k, l, m;
    static int ii, en, jj, ll, mm, nn;
    static mpfr_t si;
    mpfr_init2(si, wp2);
    static mpfr_t ti;
    mpfr_init2(ti, wp2);
    static mpfr_t sr;
    mpfr_init2(sr, wp2);
    static mpfr_t tr;
    mpfr_init2(tr, wp2);
    static int im1, ip1, mp1, its;
    static int enm1, iend;
    static mpfr_t norm;
    mpfr_init2(norm, wp2);
    static mpfr_t machep;
    mpfr_init2(machep, wp2);

	mpc_t *x = equiv_1;
	mpc_t *y = equiv_3;
	mpc_t *z__ = equiv_5;
	mpfr_t *t1 = (mpfr_t*) equiv_1;
	mpfr_t *t2 = (mpfr_t*) equiv_3;
	mpfr_t *t3 = (mpfr_t*) equiv_5;
	mpfr_t *xi = (mpfr_t*) equiv_1 + 1;
	mpfr_t *yi = (mpfr_t*) equiv_3 + 1;
	mpfr_t *xr = (mpfr_t*) equiv_1;
	mpfr_t *yr = (mpfr_t*) equiv_3;
	mpfr_t *zzi = (mpfr_t*) equiv_5 + 1;
	mpfr_t *zzr = (mpfr_t*) equiv_5;


/*     REAL ABS */
/*     int MIN0 */
/*     COMPLEX CSQRT,CMPLX */

/*     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE COMLR2, */
/*     NUM. MATH. 16, 181-204(1970) BY PETERS AND WILKINSON. */
/*     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 372-395(1971). */

/*     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS */
/*     OF A COMPLEX UPPER HESSENBERG MATRIX BY THE MODIFIED LR */
/*     METHOD. THE EIGENVECTORS OF A COMPLEX GENERAL MATRIX */
/*     CAN ALSO BE FOUND IF COMHES HAS BEEN USED TO REDUCE */
/*     THIS GENERAL MATRIX TO HESSENBERG FORM. */

/*     ON INPUT- */

/*        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL */
/*           ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM */
/*           DIMENSION STATEMENT, */

/*        N IS THE ORDER OF THE MATRIX, */

/*     NOBAK IS AN EXTRA  PARAMETER WHICH  TELLS IF THE USER WANTS */
/*     1     THE EIGENVECTORS  OR  JUST  THE  TRANSFORMATIONMATRICES (=0) */

/*        LOW AND IGH ARE intS DETRMINED BY THE BALANCING */
/*           SUBROUTINE CBAL. IF CBAL HAS NOT BEEN USED, */
/*           SET LOW=1, IGH=N, */

/*        INT CONTAINS INFORMATION ON THE ROWS AND COLUMNS INTERCHANGED */
/*           IN THE REDUCTION BY COMHES, IF PERFORMED. ONLY ELEMENTS */
/*           LOW THROUGH IGH ARE USED. IF THE EIGENVECTORS OF THE HESSEN- */
/*           BERG MATRIX ARE DESIRED, SET INT(J)=J FOR THESE ELEMENTS, */

/*        HR AND HI CONTAINS THE REAL AND IMAGINARY PARTS, */
/*           RESPECTIVELY, OF THE COMPLEX UPPER HESSENBERG MATRIX. */
/*           THEIR LOWER TRIANGLES BELOW THE SUBDIAGONAL CONTAINS THE */
/*           MULTIPLIERS WHICH WERE USED IN THE REDUCTION BY COMHES, */
/*           IF PERFORMED. IF THE EIGENVECTORS OF THE HESSENBERG */
/*           MATRIX ARE DESIRED, THESE ELEMENTS MUST BE SET TO ZERO. */

/*     ON OUTPUT- */

/*        THE UPPER HESSENBERG PORTIONS OF HR AND HI HAVE BEEN */
/*           DESTROYED, */

/*        WR AND WI CONTAINS THE REAL AND IMAGINARY PARTS, */
/*           RESPECTIVELY, OF THE EIGENVALUES. IF AN ERROR */
/*           EXIT IS MADE, THE EIGENVALUES SHOULD BE CORRECT */
/*           FOR INDICES IERR+1,...,N, */

/*        ZR AND ZI CONTAINS THE REAL AND IMAGINARY PARTS, */
/*           RESPECTIVELY, OF THE EIGENVECTORS. THE EIGENVECTORS */
/*           ARE UNNORMALIZED. IF AN ERROR EXIT IS MADE, NONE OF */
/*           THE EIGENVECTORS HAS BEEN FOUND, */

/*        IERR IS SET TO */
/*           ZERO     FOR NORMAL RETURN, */
/*           J         IF THE J-TH EIGENVALUE HAS NOT BEEN */
/*                     DETERMINED AFTER 30 ITERATIONS. */

/*     ARITHMETIC IS REAL EXCEPT FOR THE REPLACEMENT OF THE ALGOL */
/*     PROCEDURE CDIV BY COMPLEX DIVISION AND USE OF THE SUBROUTINES */
/*     CSQRT AND CMPLX IN COMPUTING COMPLEX SQUARE ROOTS. */

/*     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW, */
/*     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY */

/*     ------------------------------------------------------------------ */
/*     *********** MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING */
/*                 THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC. */

/*                 ********** */
    /* Parameter adjustments */
zi_dim1 = *nm;
zi_offset = 1 + zi_dim1;
zi = zi - zi_offset;
zr_dim1 = *nm;
zr_offset = 1 + zr_dim1;
zr = zr - zr_offset;
    --wi;
    --wr;
hi_dim1 = *nm;
hi_offset = 1 + hi_dim1;
hi = hi - hi_offset;
hr_dim1 = *nm;
hr_offset = 1 + hr_dim1;
hr = hr - hr_offset;
    --int__;

    /* Function Body */
// mpfr_set_d(machep, 10, MPFR_RNDN);
// mpfr_pow_si(machep, machep, -15, MPFR_RNDN);
// // mpfr_mul_d(machep, machep, 7.1054273576010019e0, MPFR_RNDN);
// mpfr_mul_ui(machep, machep, 7, MPFR_RNDN);
    mpfr_set_d(machep, 2, MPFR_RNDN);
	mpfr_pow_si(machep, machep, (int) (-0.95 * ((double) wp2)), MPFR_RNDN);
	// printf("machep = ");
	// mpfr_out_str(stdout, 10, 0, machep, MPFR_RNDN);
	//getchar();
	


*ierr = 0;
/*     ********** INITIALIZE EIGENVECTOR MATRIX ********** */
    if(1 > *n) {
	goto L10000;
    }
    for (i__ = 1; i__ <= *n; ++i__) {

	if(1 > *n) {
	    goto L10000;
	}
	for (j = 1; j <= *n; ++j) {
mpfr_set_d(zr[i__ + j * zr_dim1], 0, MPFR_RNDN);
mpfr_set_d(zi[i__ + j * zi_dim1], 0, MPFR_RNDN);
	    if(i__ == j) {
mpfr_set_d(zr[i__ + j * zr_dim1], 1, MPFR_RNDN);
	    }
/* L100: */
	}
    }
L10000:;
// printf("L10000:;\n");
/*     ********** FORM THE MATRIX OF ACCUMULATED TRANSFORMATIONS */
/*                FROM THE INFORMATION LEFT BY COMHES ********** */
iend = *igh - *low - 1;
    if(iend <= 0) {
	goto L180;
    }
/*     ********** FOR I=IGH-1 STEP -1 UNTIL LOW+1 DO -- ********** */
    if(1 > iend) {
	goto L10010;
    }
    for (ii = 1; ii <= iend; ++ii) {
i__ = *igh - ii;
ip1 = i__ + 1;

	if(ip1 > *igh) {
	    goto L10020;
	}
	for (k = ip1; k <= *igh; ++k) {
mpfr_set(zr[k + i__ * zr_dim1], hr[k + (i__ - 1) * hr_dim1], MPFR_RNDN);
mpfr_set(zi[k + i__ * zi_dim1], hi[k + (i__ - 1) * hi_dim1], MPFR_RNDN);
/* L120: */
	}
L10020:;
// printf("L10020:;\n");

j = int__[i__];
	if(i__ == j) {
	    goto L160;
	}

	if(i__ > *igh) {
	    goto L10030;
	}
	for (k = i__; k <= *igh; ++k) {
mpfr_set(zr[i__ + k * zr_dim1], zr[j + k * zr_dim1], MPFR_RNDN);
mpfr_set(zi[i__ + k * zi_dim1], zi[j + k * zi_dim1], MPFR_RNDN);
mpfr_set_d(zr[j + k * zr_dim1], 0, MPFR_RNDN);
mpfr_set_d(zi[j + k * zi_dim1], 0, MPFR_RNDN);
/* L140: */
	}
L10030:;
// printf("L10030:;\n");

mpfr_set_d(zr[j + i__ * zr_dim1], 1, MPFR_RNDN);
L160:;
// printf("L160:;\n");
	;
    }
L10010:;
// printf("L10010:;\n");
/*     ********** STORE ROOTS ISOLATED BY CBAL ********** */
L180:;
// printf("L180:;\n");
    if(1 > *n) {
	goto L10040;
    }
    for (i__ = 1; i__ <= *n; ++i__) {
	if(i__ >= *low && i__ <= *igh) {
	    goto L200;
	}
mpfr_set(wr[i__], hr[i__ + i__ * hr_dim1], MPFR_RNDN);
mpfr_set(wi[i__], hi[i__ + i__ * hi_dim1], MPFR_RNDN);
// printf("wr[%d] = ", i__);
// mpfr_out_str(stdout, 10, 0, wr[i__], MPFR_RNDN);
// printf("\n");

L200:;
// printf("L200:;\n");
	;
    }
L10040:;
// printf("L10040:;\n");

en = *igh;
mpfr_set_d(tr, 0, MPFR_RNDN);
mpfr_set_d(ti, 0, MPFR_RNDN);
/*     ********** SEARCH FOR NEXT EIGENVALUE ********** */
L220:;
// printf("L220:;\n");
    if(en < *low) {
	goto L680;
    }
its = 0;
enm1 = en - 1;
/*     ********** LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT */
/*                FOR L=EN STEP -1 UNTIL LOW DO -- ********** */
L240:;
// printf("L240:;\n");
    if(*low > en) {
	goto L10050;
    }
    for (ll = *low; ll <= en; ++ll) {
l = en + *low - ll;
	if(l == *low) {
	    goto L300;
	}
mpfr_abs(r__1, hr[l + (l - 1) * hr_dim1], MPFR_RNDN);
mpfr_abs(r__2, hi[l + (l - 1) * hi_dim1], MPFR_RNDN);
mpfr_abs(r__3, hr[l - 1 + (l - 1) * hr_dim1], MPFR_RNDN);
mpfr_abs(r__4, hi[l - 1 + (l - 1) * hi_dim1], MPFR_RNDN);
mpfr_abs(r__5, hr[l + l * hr_dim1], MPFR_RNDN);
mpfr_abs(r__6, hi[l + l * hi_dim1], MPFR_RNDN);
mpfr_add(r__1, r__1, r__2, MPFR_RNDN);
mpfr_add(r__3, r__3, r__4, MPFR_RNDN);
mpfr_add(r__3, r__3, r__5, MPFR_RNDN);
mpfr_add(r__3, r__3, r__6, MPFR_RNDN);
mpfr_mul(r__3, machep, r__3, MPFR_RNDN);
	if(mpfr_lessequal_p(r__1, r__3)) {
	    goto L300;
	}
/* L260: */
    }
L10050:;
// printf("L10050:;\n");
/*     ********** FORM SHIFT *********** */
L300:;
// printf("L300:;\n");
    if(l == en) {
	goto L660;
    }
    // if(its == 30) {
    if(its == 300000) {
	goto L1000;
    }
    // if(its == 10 || its == 20) {
    if(its == 100000 || its == 200000) {
	goto L320;
    }
mpfr_set(sr, hr[en + en * hr_dim1], MPFR_RNDN);
mpfr_set(si, hi[en + en * hi_dim1], MPFR_RNDN);
mpfr_mul(tmp0, hr[enm1 + en * hr_dim1], hr[en + enm1 * hr_dim1], MPFR_RNDN);
mpfr_mul(tmp1, hi[enm1 + en * hi_dim1], hi[en + enm1 * hi_dim1], MPFR_RNDN);
mpfr_sub(*xr, tmp0, tmp1, MPFR_RNDN);
mpfr_mul(tmp0, hr[enm1 + en * hr_dim1], hi[en + enm1 * hi_dim1], MPFR_RNDN);
mpfr_mul(tmp1, hi[enm1 + en * hi_dim1], hr[en + enm1 * hr_dim1], MPFR_RNDN);
mpfr_add(*xi, tmp0, tmp1, MPFR_RNDN);
    if(mpfr_equal_p(*xr, mpfr_zero) && mpfr_equal_p(*xi, mpfr_zero)) {
	goto L340;
    }
mpfr_sub(tmp0, hr[enm1 + enm1 * hr_dim1], sr, MPFR_RNDN);
mpfr_div_ui(*yr, tmp0, 2, MPFR_RNDN);
mpfr_sub(tmp0, hi[enm1 + enm1 * hi_dim1], si, MPFR_RNDN);
mpfr_div_ui(*yi, tmp0, 2, MPFR_RNDN);
mpfr_mul(tmp0, *yr, *yr, MPFR_RNDN);
mpfr_mul(tmp1, *yi, *yi, MPFR_RNDN);
mpfr_sub(tmp0, tmp0, tmp1, MPFR_RNDN);
mpfr_add(tmp0, tmp0, *xr, MPFR_RNDN);
mpfr_mul_ui(tmp1, *yr, 2, MPFR_RNDN);
mpfr_mul(tmp1, tmp1, *yi, MPFR_RNDN);
mpfr_add(tmp1, tmp1, *xi, MPFR_RNDN);
mpc_set_fr_fr(*z__, tmp0, tmp1, MPFR_RNDN);
	mpc_sqrt(*z__, *z__, MPFR_RNDN);
mpfr_mul(tmp0, *yr, *zzr, MPFR_RNDN);
mpfr_mul(tmp1, *yi, *zzi, MPFR_RNDN);
mpfr_add(cmp0, tmp0, tmp1, MPFR_RNDN);
    if(mpfr_less_p(cmp0, mpfr_zero)) {
		mpc_neg(*z__, *z__, MPFR_RNDN);
	}
	mpc_add(q__2, *y, *z__, MPFR_RNDN);
    mpc_div(*x, *x, q__2, MPFR_RNDN);
mpfr_sub(sr, sr, *xr, MPFR_RNDN);
mpfr_sub(si, si, *xi, MPFR_RNDN);
    goto L340;
/*     ********** FORM EXCEPTIONAL SHIFT ********** */
L320:;
// printf("L320:;\n");
mpfr_abs(r__1, hr[en + enm1 * hr_dim1], MPFR_RNDN);
mpfr_abs(sr, hr[enm1 + (en - 2) * hr_dim1], MPFR_RNDN);
mpfr_add(sr, r__1, sr, MPFR_RNDN);
mpfr_abs(r__1, hi[en + enm1 * hi_dim1], MPFR_RNDN);
mpfr_abs(si, hi[enm1 + (en - 2) * hi_dim1], MPFR_RNDN);
mpfr_add(si, r__1, si, MPFR_RNDN);

L340:;
// printf("L340:;\n");
    if(*low > en) {
	goto L10060;
    }
    for (i__ = *low; i__ <= en; ++i__) {
mpfr_sub(hr[i__ + i__ * hr_dim1], hr[i__ + i__ * hr_dim1], sr, MPFR_RNDN);
mpfr_sub(hi[i__ + i__ * hi_dim1], hi[i__ + i__ * hi_dim1], si, MPFR_RNDN);
/* L360: */
    }
L10060:;
// printf("L10060:;\n");

mpfr_add(tr, tr, sr, MPFR_RNDN);
mpfr_add(ti, ti, si, MPFR_RNDN);
    ++its;
	// printf("its = %d\n", its);
/*     ********** LOOK FOR TWO CONSECUTIVE SMALL */
/*                SUB-DIAGONAL ELEMENTS ********** */
mpfr_abs(r__1, hr[enm1 + enm1 * hr_dim1], MPFR_RNDN);
mpfr_abs(*xr, hi[enm1 + enm1 * hi_dim1], MPFR_RNDN);
mpfr_add(*xr, r__1, *xr, MPFR_RNDN);
mpfr_abs(r__1, hr[en + enm1 * hr_dim1], MPFR_RNDN);
mpfr_abs(*yr, hi[en + enm1 * hi_dim1], MPFR_RNDN);
mpfr_add(*yr, r__1, *yr, MPFR_RNDN);
mpfr_abs(r__1, hr[en + en * hr_dim1], MPFR_RNDN);
mpfr_abs(*zzr, hi[en + en * hi_dim1], MPFR_RNDN);
mpfr_add(*zzr, r__1, *zzr, MPFR_RNDN);
	
/*     ********** FOR M=EN-1 STEP -1 UNTIL L DO -- ********** */
    if(l > enm1) {
	goto L10070;
    }
    for (mm = l; mm <= enm1; ++mm) {
m = enm1 + l - mm;
	if(m == l) {
	    goto L420;
	}
mpfr_set(*yi, *yr, MPFR_RNDN);
mpfr_abs(r__1, hr[m + (m - 1) * hr_dim1], MPFR_RNDN);
mpfr_abs(*yr, hi[m + (m - 1) * hi_dim1], MPFR_RNDN);
mpfr_add(*yr, r__1, *yr, MPFR_RNDN);
mpfr_set(*xi, *zzr, MPFR_RNDN);
mpfr_set(*zzr, *xr, MPFR_RNDN);
mpfr_abs(r__1, hr[m - 1 + (m - 1) * hr_dim1], MPFR_RNDN);
mpfr_abs(*xr, hi[m - 1 + (m - 1) * hi_dim1], MPFR_RNDN);
mpfr_add(*xr, r__1, *xr, MPFR_RNDN);
mpfr_add(cmp0, *zzr, *xr, MPFR_RNDN);
mpfr_add(cmp0, cmp0, *xi, MPFR_RNDN);
mpfr_mul(tmp0, cmp0, machep, MPFR_RNDN);
mpfr_mul(tmp0, tmp0, *zzr, MPFR_RNDN);
mpfr_div(cmp0, tmp0, *yi, MPFR_RNDN);
	if(mpfr_lessequal_p(*yr, cmp0)) {
	    goto L420;
	}
/* L380: */
    }
L10070:;
// printf("L10070:;\n");
/*     ********** TRIANGULAR DECOMPOSITION H=L* R ********** */
L420:;
// printf("L420:;\n");
mp1 = m + 1;

    if(mp1 > en) {
	goto L10080;
    }
    for (i__ = mp1; i__ <= en; ++i__) {
im1 = i__ - 1;
mpfr_set(*xr, hr[im1 + im1 * hr_dim1], MPFR_RNDN);
mpfr_set(*xi, hi[im1 + im1 * hi_dim1], MPFR_RNDN);
mpfr_set(*yr, hr[i__ + im1 * hr_dim1], MPFR_RNDN);
mpfr_set(*yi, hi[i__ + im1 * hi_dim1], MPFR_RNDN);
mpfr_abs(cmp0, *xr, MPFR_RNDN);
mpfr_abs(cmp1, *xi, MPFR_RNDN);
mpfr_add(cmp0, cmp0, cmp1, MPFR_RNDN);
mpfr_abs(cmp2, *yr, MPFR_RNDN);
mpfr_abs(cmp3, *yi, MPFR_RNDN);
mpfr_add(cmp2, cmp2, cmp3, MPFR_RNDN);
	if(mpfr_greaterequal_p(cmp0, cmp2)) {
	    goto L460;
	}
/*     ********** INTERCHANGE ROWS OF HR AND HI ********** */
	if(im1 > *n) {
	    goto L10090;
	}
	for (j = im1; j <= *n; ++j) {
mpfr_set(*zzr, hr[im1 + j * hr_dim1], MPFR_RNDN);
mpfr_set(hr[im1 + j * hr_dim1], hr[i__ + j * hr_dim1], MPFR_RNDN);
mpfr_set(hr[i__ + j * hr_dim1], *zzr, MPFR_RNDN);
mpfr_set(*zzi, hi[im1 + j * hi_dim1], MPFR_RNDN);
mpfr_set(hi[im1 + j * hi_dim1], hi[i__ + j * hi_dim1], MPFR_RNDN);
mpfr_set(hi[i__ + j * hi_dim1], *zzi, MPFR_RNDN);
/* L440: */
	}
L10090:;
// printf("L10090:;\n");

	mpc_div(*z__, *x, *y, MPFR_RNDN);
mpfr_set_d(wr[i__], 1, MPFR_RNDN);
// printf("setto wr[%d] = ", i__);
// mpfr_out_str(stdout, 10, 0, wr[i__], MPFR_RNDN);
// printf("\n");
// 
	goto L480;
L460:;
// printf("L460:;\n");
	mpc_div(*z__, *y, *x, MPFR_RNDN);

mpfr_set_d(wr[i__], -1, MPFR_RNDN);
// printf("setto wr[%d] = ", i__);
// mpfr_out_str(stdout, 10, 0, wr[i__], MPFR_RNDN);
// printf("\n");
// 
L480:;
// printf("L480:;\n");
mpfr_set(hr[i__ + im1 * hr_dim1], *zzr, MPFR_RNDN);
mpfr_set(hi[i__ + im1 * hi_dim1], *zzi, MPFR_RNDN);

	if(i__ > *n) {
	    goto L10100;
	}

	for (j = i__; j <= *n; ++j) {
mpfr_mul(tmp0, *zzr, hr[im1 + j * hr_dim1], MPFR_RNDN);
mpfr_mul(tmp1, *zzi, hi[im1 + j * hi_dim1], MPFR_RNDN);
mpfr_sub(hr[i__ + j * hr_dim1], hr[i__ + j * hr_dim1], tmp0, MPFR_RNDN);
mpfr_add(hr[i__ + j * hr_dim1], hr[i__ + j * hr_dim1], tmp1, MPFR_RNDN);
mpfr_mul(tmp0, *zzr, hi[im1 + j * hi_dim1], MPFR_RNDN);
mpfr_mul(tmp1, *zzi, hr[im1 + j * hr_dim1], MPFR_RNDN);
mpfr_sub(hi[i__ + j * hi_dim1], hi[i__ + j * hi_dim1], tmp0, MPFR_RNDN);
mpfr_sub(hi[i__ + j * hi_dim1], hi[i__ + j * hi_dim1], tmp1, MPFR_RNDN);
/* L500: */
	}
L10100:;
// printf("L10100:;\n");

/* L520: */
	;
    }
L10080:;
// printf("L10080:;\n");
/*     ********** COMPOSITION R*L=H ********** */
    if(mp1 > en) {
	goto L10110;
    }
    for (j = mp1; j <= en; ++j) {
mpfr_set(*xr, hr[j + (j - 1) * hr_dim1], MPFR_RNDN);
mpfr_set(*xi, hi[j + (j - 1) * hi_dim1], MPFR_RNDN);
mpfr_set_d(hr[j + (j - 1) * hr_dim1], 0e0, MPFR_RNDN);
mpfr_set_d(hi[j + (j - 1) * hi_dim1], 0e0, MPFR_RNDN);
/*     ********** INTERCHANGE COLUMNS OF HR, HI, ZR, AND ZI, */
/*                IF NECESSARY ********** */
	if(mpfr_lessequal_p(wr[j], mpfr_zero)) {
	    goto L580;
	}

	if(1 > j) {
	    goto L10120;
	}
	for (i__ = 1; i__ <= j; ++i__) {
mpfr_set(*zzr, hr[i__ + (j - 1) * hr_dim1], MPFR_RNDN);
mpfr_set(hr[i__ + (j - 1) * hr_dim1], hr[i__ + j * hr_dim1], MPFR_RNDN);
mpfr_set(hr[i__ + j * hr_dim1], *zzr, MPFR_RNDN);
mpfr_set(*zzi, hi[i__ + (j - 1) * hi_dim1], MPFR_RNDN);
mpfr_set(hi[i__ + (j - 1) * hi_dim1], hi[i__ + j * hi_dim1], MPFR_RNDN);
mpfr_set(hi[i__ + j * hi_dim1], *zzi, MPFR_RNDN);
/* L540: */
	}
L10120:;
// printf("L10120:;\n");

	if(*low > *igh) {
	    goto L10130;
	}
	for (i__ = *low; i__ <= *igh; ++i__) {
mpfr_set(*zzr, zr[i__ + (j - 1) * zr_dim1], MPFR_RNDN);
mpfr_set(zr[i__ + (j - 1) * zr_dim1], zr[i__ + j * zr_dim1], MPFR_RNDN);
mpfr_set(zr[i__ + j * zr_dim1], *zzr, MPFR_RNDN);
mpfr_set(*zzi, zi[i__ + (j - 1) * zi_dim1], MPFR_RNDN);
mpfr_set(zi[i__ + (j - 1) * zi_dim1], zi[i__ + j * zi_dim1], MPFR_RNDN);
mpfr_set(zi[i__ + j * zi_dim1], *zzi, MPFR_RNDN);
/* L560: */
	}
L10130:;
// printf("L10130:;\n");

L580:;
// printf("L580:;\n");
	if(1 > j) {
	    goto L10140;
	}
	for (i__ = 1; i__ <= j; ++i__) {
mpfr_mul(tmp0, *xr, hr[i__ + j * hr_dim1], MPFR_RNDN);
mpfr_mul(tmp1, *xi, hi[i__ + j * hi_dim1], MPFR_RNDN);
mpfr_add(hr[i__ + (j - 1) * hr_dim1], hr[i__ + (j - 1) * hr_dim1], tmp0, MPFR_RNDN);
mpfr_sub(hr[i__ + (j - 1) * hr_dim1], hr[i__ + (j - 1) * hr_dim1], tmp1, MPFR_RNDN);
mpfr_mul(tmp0, *xr, hi[i__ + j * hi_dim1], MPFR_RNDN);
mpfr_mul(tmp1, *xi, hr[i__ + j * hr_dim1], MPFR_RNDN);
mpfr_add(hi[i__ + (j - 1) * hi_dim1], hi[i__ + (j - 1) * hi_dim1], tmp0, MPFR_RNDN);
mpfr_add(hi[i__ + (j - 1) * hi_dim1], hi[i__ + (j - 1) * hi_dim1], tmp1, MPFR_RNDN);
/* L600: */
	}
L10140:;
// printf("L10140:;\n");
/*     ********** ACCUMULATE TRANSFORMATIONS ********** */
	if(*low > *igh) {
	    goto L10150;
	}
	for (i__ = *low; i__ <= *igh; ++i__) {
mpfr_mul(tmp0, *xr, zr[i__ + j * zr_dim1], MPFR_RNDN);
mpfr_mul(tmp1, *xi, zi[i__ + j * zi_dim1], MPFR_RNDN);
mpfr_add(zr[i__ + (j - 1) * zr_dim1], zr[i__ + (j - 1) * zr_dim1], tmp0, MPFR_RNDN);
mpfr_sub(zr[i__ + (j - 1) * zr_dim1], zr[i__ + (j - 1) * zr_dim1], tmp1, MPFR_RNDN);
mpfr_mul(tmp0, *xr, zi[i__ + j * zi_dim1], MPFR_RNDN);
mpfr_mul(tmp1, *xi, zr[i__ + j * zr_dim1], MPFR_RNDN);
mpfr_add(zi[i__ + (j - 1) * zi_dim1], zi[i__ + (j - 1) * zi_dim1], tmp0, MPFR_RNDN);
mpfr_add(zi[i__ + (j - 1) * zi_dim1], zi[i__ + (j - 1) * zi_dim1], tmp1, MPFR_RNDN);
/* L620: */
	}
L10150:;
// printf("L10150:;\n");

/* L640: */
	;
    }
L10110:;
// printf("L10110:;\n");

    goto L240;
/*     ********** A ROOT FOUND ********** */
L660:;
// printf("L660:;\n");
mpfr_add(hr[en + en * hr_dim1], hr[en + en * hr_dim1], tr, MPFR_RNDN);
mpfr_set(wr[en], hr[en + en * hr_dim1], MPFR_RNDN);
// printf("wr[%d] = ", en);
// mpfr_out_str(stdout, 10, 0, wr[en], MPFR_RNDN);
// printf("\n");

mpfr_add(hi[en + en * hi_dim1], hi[en + en * hi_dim1], ti, MPFR_RNDN);
mpfr_set(wi[en], hi[en + en * hi_dim1], MPFR_RNDN);
en = enm1;
    goto L220;
/*     ********** ALL ROOTS FOUND. BACKSUBSTITUTE TO FIND */
/*                VECTORS OF UPPER TRIANGULAR FORM ********** */
L680:;
// printf("L680:;\n");
    if(*n == 1) {
	goto L1001;
    }
    if(*nobak == 0) {
	goto L1001;
    }
mpfr_set_d(norm, 0, MPFR_RNDN);

    if(1 > *n) {
	goto L10160;
    }
    for (i__ = 1; i__ <= *n; ++i__) {

	if(i__ > *n) {
	    goto L10160;
	}
	for (j = i__; j <= *n; ++j) {
mpfr_abs(r__1, hr[i__ + j * hr_dim1], MPFR_RNDN);
mpfr_abs(r__2, hi[i__ + j * hi_dim1], MPFR_RNDN);
mpfr_add(norm, norm, r__1, MPFR_RNDN);
mpfr_add(norm, norm, r__2, MPFR_RNDN);
/* L720: */
	}
    }
L10160:;
// printf("L10160:;\n");
/*     ********** FOR EN=N STEP -1 UNTIL 2 DO -- ********** */
    if(2 > *n) {
	goto L10170;
    }
    for (nn = 2; nn <= *n; ++nn) {
en = *n + 2 - nn;
mpfr_set(*xr, wr[en], MPFR_RNDN);
mpfr_set(*xi, wi[en], MPFR_RNDN);
enm1 = en - 1;
/*     ********** FOR I=EN-1 STEP -1 UNTIL 1 DO -- ********** */
	if(1 > enm1) {
	    goto L10180;
	}
	for (ii = 1; ii <= enm1; ++ii) {
i__ = en - ii;
mpfr_set(*zzr, hr[i__ + en * hr_dim1], MPFR_RNDN);
mpfr_set(*zzi, hi[i__ + en * hi_dim1], MPFR_RNDN);
	    if(i__ == enm1) {
		goto L760;
	    }
ip1 = i__ + 1;

	    if(ip1 > enm1) {
		goto L10190;
		}
	    for (j = ip1; j <= enm1; ++j) {
mpfr_mul(tmp0, hr[i__ + j * hr_dim1], hr[j + en * hr_dim1], MPFR_RNDN);
mpfr_mul(tmp1, hi[i__ + j * hi_dim1], hi[j + en * hi_dim1], MPFR_RNDN);
mpfr_add(*zzr, *zzr, tmp0, MPFR_RNDN);
mpfr_sub(*zzr, *zzr, tmp1, MPFR_RNDN);
mpfr_mul(tmp0, hr[i__ + j * hr_dim1], hi[j + en * hi_dim1], MPFR_RNDN);
mpfr_mul(tmp1, hi[i__ + j * hi_dim1], hr[j + en * hr_dim1], MPFR_RNDN);
mpfr_add(*zzi, *zzi, tmp0, MPFR_RNDN);
mpfr_add(*zzi, *zzi, tmp1, MPFR_RNDN);
/* L740: */
	    }
L10190:;
// printf("L10190:;\n");

L760:;
// printf("L760:;\n");
mpfr_sub(*yr, *xr, wr[i__], MPFR_RNDN);
mpfr_sub(*yi, *xi, wi[i__], MPFR_RNDN);
	    if(mpfr_equal_p(*yr, mpfr_zero) && mpfr_equal_p(*yi, mpfr_zero)) {
mpfr_mul(*yr, machep, norm, MPFR_RNDN);
	    }
mpc_div(*z__, *z__, *y, MPFR_RNDN);
mpfr_set(hr[i__ + en * hr_dim1], t3[0], MPFR_RNDN);
mpfr_set(hi[i__ + en * hi_dim1], t3[1], MPFR_RNDN);
/* L780: */
	}
L10180:;
// printf("L10180:;\n");

/* L800: */
	;
    }
L10170:;
// printf("L10170:;\n");
/*     ********** END BACKSUBSTITUTION ********** */
enm1 = *n - 1;
/*     ********** VECTORS OF ISOLATED ROOTS ********** */
    if(1 > enm1) {
	goto L10200;
    }
    for (i__ = 1; i__ <= enm1; ++i__) {
	if(i__ >= *low && i__ <= *igh) {
	    goto L840;
	}
ip1 = i__ + 1;

	if(ip1 > *n) {
	    goto L10210;
	}
	for (j = ip1; j <= *n; ++j) {
mpfr_set(zr[i__ + j * zr_dim1], hr[i__ + j * hr_dim1], MPFR_RNDN);
mpfr_set(zi[i__ + j * zi_dim1], hi[i__ + j * hi_dim1], MPFR_RNDN);
/* L820: */
	}
L10210:;
// printf("L10210:;\n");

L840:;
// printf("L840:;\n");
	;
    }
L10200:;
// printf("L10200:;\n");
/*     ********** MULTIPLY BY TRANSFORMATION MATRIX TO GIVE */
/*                VECTORS OF ORIGINAL FULL MATRIX. */
/*                FOR J=N STEP -1 UNTIL LOW+1 DO -- ********** */
    if(*low > enm1) {
	goto L10220;
    }
    for (jj = *low; jj <= enm1; ++jj) {
j = *n + *low - jj;
/* Computing MIN */
m = min(j - 1,*igh);

	if(*low > *igh) {
	    goto L10220;
	}
	for (i__ = *low; i__ <= *igh; ++i__) {
mpfr_set(*zzr, zr[i__ + j * zr_dim1], MPFR_RNDN);
mpfr_set(*zzi, zi[i__ + j * zi_dim1], MPFR_RNDN);

	    if(*low > m) {
		goto L10230;
	    }
	    for (k = *low; k <= m; ++k) {
mpfr_mul(tmp0, zr[i__ + k * zr_dim1], hr[k + j * hr_dim1], MPFR_RNDN);
mpfr_mul(tmp1, zi[i__ + k * zi_dim1], hi[k + j * hi_dim1], MPFR_RNDN);
mpfr_add(*zzr, *zzr, tmp0, MPFR_RNDN);
mpfr_sub(*zzr, *zzr, tmp1, MPFR_RNDN);
mpfr_mul(tmp0, zr[i__ + k * zr_dim1], hi[k + j * hi_dim1], MPFR_RNDN);
mpfr_mul(tmp1, zi[i__ + k * zi_dim1], hr[k + j * hr_dim1], MPFR_RNDN);
mpfr_add(*zzi, *zzi, tmp0, MPFR_RNDN);
mpfr_add(*zzi, *zzi, tmp1, MPFR_RNDN);
/* L860: */
	    }
L10230:;
// printf("L10230:;\n");

mpfr_set(zr[i__ + j * zr_dim1], *zzr, MPFR_RNDN);
mpfr_set(zi[i__ + j * zi_dim1], *zzi, MPFR_RNDN);
/* L880: */
	}
    }
L10220:;
// printf("L10220:;\n");

    goto L1001;
/*     ********** SET ERROR -- NO CONVERGENCE TO AN */
/*                EIGENVALUE AFTER 30 ITERATIONS ********** */

L1000:;
// printf("L1000:;\n");
*ierr = en;
L1001:;
// printf("L1001:;\n");
// printf("WR:\n");
for (int i=1; i<=*n; i++) {
	// mpfr_out_str(stdout, 10, 0, wr[i], MPFR_RNDN);
	// printf("\n");
}
//getchar();
    return 0;
/*     ********** LAST CARD OF COMLR2 ********** */
} /* comlr2_ */


int cbabk2_(int *nm, int *n, int *low, int *igh, mpfr_t *scale, int *m, mpfr_t *zr, mpfr_t *zi)
{
	
	mpfr_t mpfr_zero, mpfr_one;
	mpfr_init2(mpfr_zero, wp2);
	mpfr_set_d(mpfr_zero, (double) 0, MPFR_RNDN);
	mpfr_init2(mpfr_one, wp2);
	mpfr_set_d(mpfr_one, (double) 1, MPFR_RNDN);
    /* System generated locals */
    int zr_dim1, zr_offset, zi_dim1, zi_offset, i__1, i__2;

    /* Local variables */
    static int i__, j, k;
    static mpfr_t s;
    mpfr_init2(s, wp2);
    static int ii;



/*     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE */
/*     CBARK2, WHICH IS A COMPLEX VERSION OF BALBAK, */
/*     NUM. MATH. 13, 293-304(1969) BY PARLETT AND REINSCH. */
/*     HANDBOOK FOR AUTO. COM., VOL.II-LINEAR ALGEBRA, 315-326(1971). */

/*     THIS SUBROUTINE FORMS THE EIGENVECTORS OF A COMPLEX GENERAL */
/*     MATRIX BY BACK TRANSFORMING THOSE OF THE CORRESPONDING */
/*     BALANCED MATRIX DETRMINED BY CBAL. */

/*     ON INPUT- */

/*        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL */
/*           ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM */
/*           DIMENSION STATEMENT, */

/*        N IS THE ORDER OF THE MATRIX, */

/*        LOW AND IGH ARE intS DETRMINED BY CBAL, */

/*        SCALE CONTAINS INFORMATION DETERMINING THE PERMUTATIONS */
/*           AND SCALING FACTORS USED BY CBAL, */

/*        M IS THE NUMBER OF EIGENVECTORS OT BE BACK TRANSFORMED, */

/*     ZR AND ZI CONTAIN THE REAL AND IMAGINALY PARTS, */
/*           RESPECTIVELY, OF THE EIGENVECTORS TO BE */
/*           BACK TRANSFORMED IN THEIR FIRST M COLUMNS. */

/*     ON OUTPUT- */

/*        ZR AND ZI CONTAIN THE REAL AND IMAGINARY PARTS, */
/*           RESPECTIVELY, OF THE TRANSFORMED EIGENVECTORS */
/*           IN THEIR FIRST M COLUMNS. */

/*     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW, */
/*     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONALLABORATORY */

/*     ------------------------------------------------------------------ */

    /* Parameter adjustments */
    --scale;
zi_dim1 = *nm;
zi_offset = 1 + zi_dim1;
zi = zi - zi_offset;
zr_dim1 = *nm;
zr_offset = 1 + zr_dim1;
zr = zr - zr_offset;

    /* Function Body */
    if(*igh == *low) {
	goto L120;
    }

    if(*low > *igh) {
	goto L10000;
    }
    for (i__ = *low; i__ <= *igh; ++i__) {
mpfr_set(s, scale[i__], MPFR_RNDN);
/*     ********** LEFT HAND EIGENVECTORS ARE BACK TRANSFORMED */
/*                IF THE FOREGOING STATEMENT IS REPLACED BY */
/*                S=1.0/SCALE(I). ********** */
	if(1 > *m) {
	    goto L10010;
	}
	for (j = 1; j <= *m; ++j) {
mpfr_mul(zr[i__ + j * zr_dim1], zr[i__ + j * zr_dim1], s, MPFR_RNDN);
mpfr_mul(zi[i__ + j * zi_dim1], zi[i__ + j * zi_dim1], s, MPFR_RNDN);
/* L100: */
	}
L10010:;
// printf("L10010:;\n");

/* L110: */
	;
    }
L10000:;
// printf("L10000:;\n");
/*     ********** FOR I=LOW-1 STEP -1 UNTIL 1, */
/*                IGH+1 STEP 1 UNTIL N DO -- ********** */
L120:;
// printf("L120:;\n");
    if(1 > *n) {
	goto L10020;
    }
    for (ii = 1; ii <= *n; ++ii) {
i__ = ii;
	if(i__ >= *low && i__ <= *igh) {
	    goto L140;
	}
	if(i__ < *low) {
i__ = *low - ii;
	}
k = mpfr_get_si(scale[i__], MPFR_RNDD);
	if(k == i__) {
	    goto L140;
	}

	if(1 > *m) {
	    goto L10030;
	}
	for (j = 1; j <= *m; ++j) {
mpfr_set(s, zr[i__ + j * zr_dim1], MPFR_RNDN);
mpfr_set(zr[i__ + j * zr_dim1], zr[k + j * zr_dim1], MPFR_RNDN);
mpfr_set(zr[k + j * zr_dim1], s, MPFR_RNDN);
mpfr_set(s, zi[i__ + j * zi_dim1], MPFR_RNDN);
mpfr_set(zi[i__ + j * zi_dim1], zi[k + j * zi_dim1], MPFR_RNDN);
mpfr_set(zi[k + j * zi_dim1], s, MPFR_RNDN);
/* L130: */
	}
L10030:;
// printf("L10030:;\n");

L140:;
// printf("L140:;\n");
	;
    }
L10020:;
// printf("L10020:;\n");

    return 0;
/*     ********** LAST CARD OF CBABK2 ********** */
} /* cbabk2_ */


int decide_(int *nm, int *n, int *ndel, mpfr_t *cshtr, mpfr_t *cshti, int *nblock, mpfr_t *hr, mpfr_t *hi)
{
	
	mpfr_t mpfr_zero, mpfr_one;
	mpfr_init2(mpfr_zero, wp2);
	mpfr_set_d(mpfr_zero, (double) 0, MPFR_RNDN);
	mpfr_init2(mpfr_one, wp2);
	mpfr_set_d(mpfr_one, (double) 1, MPFR_RNDN);

    /* System generated locals */
    int hr_dim1, hr_offset, hi_dim1, hi_offset, i__1;

    /* Builtin functions */
    

    /* Local variables */
    static int k, mult;



/*   ************************************************************* */
/*   * */
/*   *  DECIDE IS A USER WRITTEN SUBROUTINE WHICH MAKES IT POSSIBLE */
/*   *  FOR THE USER TO CHANGE THE GROUPING OF THE NUMERICAL MUL- */
/*   *  TIPLE EIGENVALUES,AND/OR CHANGE THE VALUES OF THE MULTIPLE */
/*   *  EIGENVALUES. */
/*   *  THE FORMAL PARAMETER LIST */
/*   * */
/*   *  NM       THE ROW DIMENSION OF THE TWO-DIMENSIONAL ARRAY */
/*   *           PARAMETERS(HR,HI)AS DECLARED IN THE CALLING */
/*   *           PROGRAM DIMENSION STATEMENT */
/*   *  N        THE ORDER OF THE ORIGINAL MATRIX */
/*   *  NBLOCK   THE NUMBER OF BLOCKS I.E THE NUMBER OF NUMERICAL */
/*   *           MULTIPLE EIGENVALUES */
/*   *  CSHTR,   THE REAL AND IMAGINARY PARTS,RESPECTIVELY,OF THE */
/*   *  CSHTI    NUMERICAL MULTIPLE EIGENVALUES */
/*   *  NDEL     INDICATES THE MULTIPLICITY OF THE NUMERICAL */
/*   *           MULTIPLE EIGENVALUES.NDEL(I+1)-NDEL(I)=THE */
/*   *           MULTIPLICITY OF EIGENVALUE I FOR I=1,..,NBLOCK */
/*   *  HR,HI    THE REAL AND IMAGINARY PARTS,RESPECTIVELY, OF */
/*   *           THE UPPER TRIANGULAR MATRIX RESULTING FROM */
/*   *           STEPS 1-3 OF THE ALGORITHM */
/*   * */
/*   *  NOTE -- THE USER MUST NOT CHANGE HR,HI */
/*   * */
/*   ****************************************************************** */
    /* Parameter adjustments */
hi_dim1 = *nm;
hi_offset = 1 + hi_dim1;
hi = hi - hi_offset;
hr_dim1 = *nm;
hr_offset = 1 + hr_dim1;
hr = hr - hr_offset;
    --ndel;
    --cshtr;
    --cshti;


    for (k = 1; k <= *nblock; ++k) {
mult = ndel[k + 1] - ndel[k];
    }
    return 0;
} /* decide_ */


int csvd1_(mpc_t *a, int *mmax, int *nmax, int *m, int *n, int *p, int *nu, int *nv, mpfr_t *s, mpc_t *u, mpc_t *v, int N, int reallc)
{
	// for (int i=0; i<*m; i++) {
	// 	for (int j=0; j<*n; j++) {
	// 		mpc_out_str(stdout, 10, 0, a[i+j*(*m)], MPFR_RNDN); printf("\n");
	// 	}
	// }
	mpfr_t mpfr_zero, mpfr_one;
	mpfr_init2(mpfr_zero, wp2);
	mpfr_set_d(mpfr_zero, 0, MPFR_RNDN);
	mpfr_init2(mpfr_one, wp2);
	mpfr_set_d(mpfr_one, 1, MPFR_RNDN);
	mpfr_t tmp0;
	mpfr_init2(tmp0, wp2);
	mpfr_t tmp1;
	mpfr_init2(tmp1, wp2);
	mpfr_t cmp0;
	mpfr_init2(cmp0, wp2);
    /* Initialized data */

    static mpfr_t eta;
    mpfr_init2(eta, wp2);
	mpfr_set_d(eta, 2, MPFR_RNDN);
    // mpfr_set_d(eta, 10, MPFR_RNDN);
	// mpfr_pow_si(eta, eta, -160, MPFR_RNDN);
	// mpfr_pow_si(eta, eta, (int) (-0.95 * ((double) wp2)), MPFR_RNDN);
	mpfr_pow_si(eta, eta, (int) (-0.95 * ((double) wp2)), MPFR_RNDN);

    /* System generated locals */
    int a_dim1, a_offset, u_dim1, u_offset, v_dim1, v_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    mpfr_t r__1;
    mpfr_init2(r__1, wp2);
    mpfr_t r__2;
    mpfr_init2(r__2, wp2);
    mpc_t q__1;
    mpc_init3(q__1, wp2, wp2);
    mpc_t q__2;
    mpc_init3(q__2, wp2, wp2);
    mpc_t q__3;
    mpc_init3(q__3, wp2, wp2);

    /* Builtin functions */
    

    /* Local variables */
    static mpfr_t *b = NULL;
		if (reallc) {
			// printf("reallc option activated\n");
			if (b) {
				// printf("free b\n");
				free(b);
				b = NULL;
			}
		}
		if(!b) {
			// printf("allocating b with N+1 = %d\n", N+1);
			b = (mpfr_t*) malloc((N+1) * sizeof(mpfr_t));
			for (int i=0; i<N+1; i++) {
				mpfr_init2(b[i], wp2);
			}
		}
		static mpfr_t *c__ = NULL;
		if (reallc) {
			// printf("reallc option activated\n");
			if (c__) {
				// printf("free c__\n");
				free(c__);
				c__ = NULL;
			}
		}
		if(!c__) {
			// printf("allocating c__ with N+1 = %d\n", N+1);
			c__ = (mpfr_t*) malloc((N+1) * sizeof(mpfr_t));
			for (int i=0; i<N+1; i++) {
				mpfr_init2(c__[i], wp2);
			}
		}
		static mpfr_t *t = NULL;
		if (reallc) {
			// printf("reallc option activated\n");
			if (t) {
				// printf("free t\n");
				free(t);
				t = NULL;
			}
		}
		if(!t) {
			// printf("allocating t with N+1 = %d\n", N+1);
			t = (mpfr_t*) malloc((N+1) * sizeof(mpfr_t));
			for (int i=0; i<N+1; i++) {
				mpfr_init2(t[i], wp2);
			}
		}

    static mpfr_t f;
    mpfr_init2(f, wp2);
    static mpfr_t g;
    mpfr_init2(g, wp2);
    static mpfr_t h__;
    mpfr_init2(h__, wp2);
    static int i__, j, k, l;
    static mpc_t q;
    mpc_init3(q, wp2, wp2);
    static mpc_t r__;
    mpc_init3(r__, wp2, wp2);
    static mpfr_t w;
    mpfr_init2(w, wp2);
    static mpfr_t x;
    mpfr_init2(x, wp2);
    static mpfr_t y;
    mpfr_init2(y, wp2);
    static mpfr_t z__;
    mpfr_init2(z__, wp2);
    static int k1, l1, n1, kk;
    static mpfr_t cs;
    mpfr_init2(cs, wp2);
    static int ll, np;
    static mpfr_t sn;
    mpfr_init2(sn, wp2);
    static mpfr_t eps;
    mpfr_init2(eps, wp2);
    static mpfr_t tol;
    mpfr_init2(tol, wp2);

	/* *********************************************************************** */


	/*     CSDV  COMPUTES THE SINGULAR VALUES TO A COMPLEX  M BY N MATRIX A */

	/*     INPUT0 */
	/*           A   M BY N COMPLEX MATRIX */
	/*           MMAX  MAXIMUM ROW-DIMENSION */
	/*           NMAX  MAXIMUM COLUMN-DIMENSION */
	/*           M     ACTUAL  ROW-DIMENSION */
	/*           N     ACTUAL  COLUMN-DIMENSION */
	/*           NU  0,N OR M , THE NUMBER OF COLUMNS OF U */
	/*           NV  0 OR N   ,THE NUMBER OF COLUMNS OF V */

	/*     OUTPUT0 */
	/*           S   THE COMPUTED SORTED SINGULAR VALUES */
	/*           U   M BY M  UNITARY MATRIX */
	/*           V   N BY N  UNITARY MATRIX */


	/*     CSVD CAN ALSO BE USED TO */

	/*       - FIND THE LEAST  SQUARES SOLUTION OF MINIMAL LENGTH */
	/*       - SOLVE A  HOMOGENEOUS SYSTEM OF LINEAR EQUATIONS */
	/*     THIS      IMPLEMENTATION OF THE ALGORITH 358(CACM) IS MADE */
	/*     BY BO KAGSTROM,AVD.F.INFORMATIONSBEH.,901 87 UME% UNIVERSITET */


	/* *********************************************************************** */
	/*   ******************************************************** */
	/*   * */
	/*   *  THE CONSTANTS USED FOR ETA AND TOL ARE MACHINE- */
	/*   *  DEPENDENT .ETA IS THE RELATIVE MACHINE PRECISION */
	/*   *  TOL IS THE SMALLEST NORMALIZED POSITIVE NUMBER */
	/*   *  DIVIDED BY ETA */
	/*   * */
	/*   *********************************************************** */
    
	/* Parameter adjustments */
v_dim1 = *nmax;
v_offset = 1 + v_dim1;
v = v - v_offset;
u_dim1 = *mmax;
u_offset = 1 + u_dim1;
u = u - u_offset;
    --s;
a_dim1 = *mmax;
a_offset = 1 + a_dim1;
a = a - a_offset;

// if (*mmax == *m && *nmax == *n) {
// 	printf("a_offset = %d\n", a_offset);
// 	for (int i=1; i<=*m; i++) {
// 		for (int j=1; j<=*n; j++) {
// 			mpc_out_str(stdout, 10, 0, a[i+j*(*m)], MPFR_RNDN); printf("\n");
// 		}
// 	}
// }

    /* Function Body */
// mpfr_set_d(tol, 10, MPFR_RNDN);
// mpfr_pow_si(tol, tol, -50, MPFR_RNDN);
mpfr_set_d(tol, 2, MPFR_RNDN);
mpfr_pow_si(tol, tol, mpfr_get_emin() + ((int) wp2), MPFR_RNDN);
// printf("TOLERANCE:\n");
// mpfr_out_str(stdout, 10, 0, tol, MPFR_RNDN);
// printf("\n");
//getchar();
np = *n + *p;
n1 = *n + 1;

	/* HOUSEHOLDER REDUCTION */
mpfr_set_d(c__[0], 0e0, MPFR_RNDN);
k = 1;
L10:;
// printf("L10:;\n");
k1 = k + 1;

	/*    ELIMINATION OF A(I,K), I=K+1,...,M */
mpfr_set_d(z__, 0e0, MPFR_RNDN);
    if(k > *m) {
	goto L10000;
    }
    for (i__ = k; i__ <= *m; ++i__) {
	mpc_norm(cmp0, a[i__ + k * a_dim1], MPFR_RNDN);
mpfr_add(z__, z__, cmp0, MPFR_RNDN);
    }
L10000:;
// printf("L10000:;\n");
mpfr_set_d(b[k - 1], 0e0, MPFR_RNDN);
    if(mpfr_lessequal_p(z__, tol)) {
	goto L70;
    }
mpfr_sqrt(z__, z__, MPFR_RNDN);
mpfr_set(b[k - 1], z__, MPFR_RNDN);
mpc_abs(w, a[k + k * a_dim1], MPFR_RNDN);
mpc_set_d_d(q, 1, 0, MPFR_RNDN);
// printf("w = ");
// mpfr_out_str(stdout, 10, 0, w, MPFR_RNDN);
// printf("\n");
    if(!mpfr_equal_p(w, mpfr_zero)) {
mpc_div_fr(q, a[k + k * a_dim1], w, MPFR_RNDN);
// printf("q");
// mpc_out_str(stdout, 10, 0, q, MPFR_RNDN);
// printf("\n");
    }
mpfr_add(tmp0, z__, w, MPFR_RNDN);
mpc_mul_fr(a[k + k * a_dim1], q, tmp0, MPFR_RNDN);
    if(k == np) {
	goto L70;
    }
    if(k1 > np) {
	goto L10010;
    }
    for (j = k1; j <= np; ++j) {
mpc_set_d_d(q, 0e0, 0e0, MPFR_RNDN);
	if(k > *m) {
	    goto L10020;
	}
	for (i__ = k; i__ <= *m; ++i__) {
/* L30: */
mpc_conj(q__1, a[i__ + k * a_dim1], MPFR_RNDN);
mpc_mul(q__1, q__1, a[i__ + j * a_dim1], MPFR_RNDN);
mpc_add(q, q, q__1, MPFR_RNDN);
	}
L10020:;
// printf("L10020:;\n");
mpfr_add(tmp0, z__, w, MPFR_RNDN);
mpfr_mul(tmp0, tmp0, z__, MPFR_RNDN);
mpc_div_fr(q, q, tmp0, MPFR_RNDN);
// printf("q\n");
// mpc_out_str(stdout, 10, 0, q, MPFR_RNDN);
// printf("\n");

	if(k > *m) {
	    goto L10030;
	}
	for (i__ = k; i__ <= *m; ++i__) {
	/* L40: */
	mpc_mul(q__1, q, a[i__ + k * a_dim1], MPFR_RNDN);
	mpc_sub(a[i__ + j * a_dim1], a[i__ + j * a_dim1], q__1, MPFR_RNDN);
	}
L10030:;
// printf("L10030:;\n");
	/* L50: */
	;
    }
L10010:;
// printf("L10010:;\n");

	/*    PHASE TRANSFORMATION */
    mpc_conj(q, a[k + k * a_dim1], MPFR_RNDN);
	mpc_neg(q, q, MPFR_RNDN);
mpc_abs(tmp0, a[k + k * a_dim1], MPFR_RNDN);
	mpc_div_fr(q, q, tmp0, MPFR_RNDN);
// printf("q\n");
// mpc_out_str(stdout, 10, 0, q, MPFR_RNDN);
// printf("\n");

    if(k1 > np) {
	goto L10040;
    }
    for (j = k1; j <= np; ++j) {
	/* L60: */
	mpc_mul(a[k + j * a_dim1], a[k + j * a_dim1], q, MPFR_RNDN);
	}
L10040:;
// printf("L10040:;\n");

	/*    ELIMINATION OF A(K,J), J=K+2,...,N */
L70:;
// printf("L70:;\n");
    if(k == *n) {
	goto L140;
    }
mpfr_set_d(z__, 0e0, MPFR_RNDN);
    if(k1 > *n) {
	goto L10050;
    }
    for (j = k1; j <= *n; ++j) {
	mpc_norm(cmp0, a[k + j * a_dim1], MPFR_RNDN);
mpfr_add(z__, z__, cmp0, MPFR_RNDN);
    }
L10050:;
// printf("L10050:;\n");
// printf("k1 = %d\n", k1);
// printf("z__ = ");
// mpfr_out_str(stdout, 10, 0, z__, MPFR_RNDN);
// printf("\n");
// if (k1 == 5 && *n == 6) {
// 	// printf("c__[%d] = ", k1-1);
// 	mpfr_out_str(stdout, 10, 0, c__[k1-1], MPFR_RNDN);
// 	// printf("\n");
// }
mpfr_set_d(c__[k1 - 1], 0e0, MPFR_RNDN);
// printf("c__[%d] = ", k1-1);
// mpfr_out_str(stdout, 10, 0, c__[k1-1], MPFR_RNDN);
// printf("\n");
// printf("c__[%d] = ", k1-1);
// mpfr_out_str(stdout, 10, 0, c__[k1-1], MPFR_RNDN);
// printf("\n");
    if(mpfr_lessequal_p(z__, tol)) {
		// printf("entrato nell'if\n");
	goto L130;
    }
mpfr_sqrt(z__, z__, MPFR_RNDN);
// printf("c__[%d] = ", k1-1);
// mpfr_out_str(stdout, 10, 0, c__[k1-1], MPFR_RNDN);
// printf("\n");
// printf("z__ = ");
// mpfr_out_str(stdout, 10, 0, z__, MPFR_RNDN);
// printf("\n");
mpfr_set(c__[k1 - 1], z__, MPFR_RNDN);
// printf("c__[%d] = ", k1-1);
// mpfr_out_str(stdout, 10, 0, c__[k1-1], MPFR_RNDN);
// printf("\n");
// if (k1 == 5 && *n == 6)
// 	getchar();
// printf("a_[%d, %d] = ", k, k1);
// mpc_out_str(stdout, 10, 0, a[k + k1 * a_dim1], MPFR_RNDN);
// printf("\n");
mpc_abs(w, a[k + k1 * a_dim1], MPFR_RNDN);
mpc_set_d_d(q, 1e0, 0e0, MPFR_RNDN);
// printf("w = ");
// mpfr_out_str(stdout, 10, 0, w, MPFR_RNDN);
// printf("\n");
    if(!mpfr_equal_p(w, mpfr_zero)) {
mpc_div_fr(q, a[k + k1 * a_dim1], w, MPFR_RNDN);
// printf("q\n");
// mpc_out_str(stdout, 10, 0, q, MPFR_RNDN);
// printf("\n");

    }
mpfr_add(tmp0, z__, w, MPFR_RNDN);
mpc_mul_fr(a[k + k1 * a_dim1], q, tmp0, MPFR_RNDN);
    if(k1 > *m) {
	goto L10060;
    }
    for (i__ = k1; i__ <= *m; ++i__) {
mpc_set_d_d(q, 0e0, 0e0, MPFR_RNDN);
	if(k1 > *n) {
	    goto L10070;
	}
	for (j = k1; j <= *n; ++j) {
/* L90: */
mpc_conj(q__1, a[k + j * a_dim1], MPFR_RNDN);
mpc_mul(q__1, q__1, a[i__ + j * a_dim1], MPFR_RNDN);
mpc_add(q, q, q__1, MPFR_RNDN);
	}
L10070:;
// printf("L10070:;\n");
mpfr_add(tmp0, z__, w, MPFR_RNDN);
mpfr_mul(tmp0, tmp0, z__, MPFR_RNDN);
mpc_div_fr(q, q, tmp0, MPFR_RNDN);
// printf("q\n");
// mpc_out_str(stdout, 10, 0, q, MPFR_RNDN);
// printf("\n");


	if(k1 > *n) {
	    goto L10080;
	}
	for (j = k1; j <= *n; ++j) {
/* L100: */
mpc_mul(q__1, q, a[k + j * a_dim1], MPFR_RNDN);
mpc_sub(a[i__ + j * a_dim1], a[i__ + j * a_dim1], q__1, MPFR_RNDN);
	}
L10080:;
// printf("L10080:;\n");
/* L110: */
	;
    }
L10060:;
// printf("L10060:;\n");

	/*    PHASE TRANSFORMATION */
    mpc_conj(q, a[k + k1 * a_dim1], MPFR_RNDN);
	mpc_neg(q, q, MPFR_RNDN);
mpc_abs(tmp0, a[k + k1 * a_dim1], MPFR_RNDN);
	mpc_div_fr(q, q, tmp0, MPFR_RNDN);
// printf("q\n");
// mpc_out_str(stdout, 10, 0, q, MPFR_RNDN);
// printf("\n");

    if(k1 > *m) {
	goto L10090;
    }
    for (i__ = k1; i__ <= *m; ++i__) {
/* L120: */
	mpc_mul(a[i__ + k1 * a_dim1], a[i__ + k1 * a_dim1], q, MPFR_RNDN);
    }
L10090:;
// printf("L10090:;\n");
L130:;
// printf("L130:;\n");
k = k1;
    goto L10;

	/* TOLERANCE FOR NEGLIGIBLE ELEMENTS */
L140:;
// printf("L140:;\n");
mpfr_set_d(eps, 0, MPFR_RNDN);
    if(1 > *n) {
	goto L10100;
    }
    for (k = 1; k <= *n; ++k) {
mpfr_set(s[k], b[k - 1], MPFR_RNDN);
mpfr_set(t[k - 1], c__[k - 1], MPFR_RNDN);
/* L150: */
	/* Computing MAX */
mpfr_add(tmp0, s[k], t[k - 1], MPFR_RNDN);
mpfr_max(eps, eps, tmp0, MPFR_RNDN);
    }
L10100:;
// printf("L10100:;\n");
mpfr_mul(eps, eps, eta, MPFR_RNDN);

	/* INITIALIZATION OF U AND V */
    if(*nu == 0) {
	goto L180;
    }
    if(1 > *nu) {
	goto L10110;
    }
    for (j = 1; j <= *nu; ++j) {
	if(1 > *m) {
	    goto L10120;
	}
	for (i__ = 1; i__ <= *m; ++i__) {
/* L160: */
mpc_set_d_d(u[i__ + j * u_dim1], 0e0, 0e0, MPFR_RNDN);
	}
L10120:;
// printf("L10120:;\n");
/* L170: */
mpc_set_d_d(u[j + j * u_dim1], 1e0, 0e0, MPFR_RNDN);
    }
L10110:;
// printf("L10110:;\n");
L180:;
// printf("L180:;\n");
    if(*nv == 0) {
	goto L210;
    }
    if(1 > *nv) {
	goto L10130;
    }
    for (j = 1; j <= *nv; ++j) {
	if(1 > *n) {
	    goto L10140;
	}
	for (i__ = 1; i__ <= *n; ++i__) {
/* L190: */
mpc_set_d_d(v[i__ + j * v_dim1], 0e0, 0e0, MPFR_RNDN);
	}
L10140:;
// printf("L10140:;\n");
/* L200: */
mpc_set_d_d(v[j + j * v_dim1], 1e0, 0e0, MPFR_RNDN);
    }
L10130:;
// printf("L10130:;\n");

	/* QR DIAGONALIZATION */
L210:;
// printf("L210:;\n");
    if(1 > *n) {
	goto L10150;
    }
    for (kk = 1; kk <= *n; ++kk) {
k = n1 - kk;

	/*    TEST FOR SPLIT */
L220:;
// printf("L220:;\n");
	if(1 > k) {
	    goto L10160;
	}
	for (ll = 1; ll <= k; ++ll) {
l = k + 1 - ll;
mpfr_abs(tmp0, t[l - 1], MPFR_RNDN);
	    if(mpfr_lessequal_p(tmp0, eps)) {
		goto L290;
	    }
mpfr_abs(tmp0, s[l - 1], MPFR_RNDN);
	    if(mpfr_lessequal_p(tmp0, eps)) {
		goto L240;
	    }
/* L230: */
	}
L10160:;
// printf("L10160:;\n");

	/*    CANCELLATION */
L240:;
// printf("L240:;\n");
mpfr_set_d(cs, 0e0, MPFR_RNDN);
mpfr_set_d(sn, 1e0, MPFR_RNDN);
l1 = l - 1;
	if(l > k) {
	    goto L10170;
	}
	for (i__ = l; i__ <= k; ++i__) {
mpfr_mul(f, sn, t[i__ - 1], MPFR_RNDN);
mpfr_abs(tmp0, f, MPFR_RNDN);
mpfr_mul(t[i__ - 1], cs, t[i__ - 1], MPFR_RNDN);
	    if(mpfr_lessequal_p(tmp0, eps)) {
		goto L290;
	    }
mpfr_set(h__, s[i__], MPFR_RNDN);
// printf("h__ 1\n");
// mpfr_out_str(stdout, 10, 0, h__, MPFR_RNDN);
// printf("\n");

mpfr_mul(tmp0, f, f, MPFR_RNDN);
mpfr_mul(tmp1, h__, h__, MPFR_RNDN);
mpfr_add(w, tmp0, tmp1, MPFR_RNDN);
mpfr_sqrt(w, w, MPFR_RNDN);
mpfr_set(s[i__], w, MPFR_RNDN);

// FT modification
if(mpfr_equal_p(w, mpfr_zero)) {
	mpfr_set_d(w,1,MPFR_RNDN);
	mpfr_set_d(cs,1,MPFR_RNDN);
	mpfr_set_d(sn,0,MPFR_RNDN);
 } else {
	mpfr_div(cs, h__, w, MPFR_RNDN);
	mpfr_div(sn, f, w, MPFR_RNDN);
 }

mpfr_neg(sn, sn, MPFR_RNDN);
	    if(*nu == 0) {
		goto L260;
	    }
	    if(1 > *m) {
		goto L10180;
	    }
	    for (j = 1; j <= *m; ++j) {
mpfr_set(x, mpc_realref(u[j + l1 * u_dim1]), MPFR_RNDN);
mpfr_set(y, mpc_realref(u[j + i__ * u_dim1]), MPFR_RNDN);
mpfr_mul(tmp0, x, cs, MPFR_RNDN);
mpfr_mul(tmp1, y, sn, MPFR_RNDN);
mpfr_add(tmp0, tmp0, tmp1, MPFR_RNDN);
mpc_set_d_d(u[j + l1 * u_dim1], 0e0, 0e0, MPFR_RNDN);
mpc_set_fr(u[j + l1 * u_dim1], tmp0, MPFR_RNDN);
/* L250: */
mpfr_mul(tmp0, y, cs, MPFR_RNDN);
mpfr_mul(tmp1, x, sn, MPFR_RNDN);
mpfr_sub(tmp0, tmp0, tmp1, MPFR_RNDN);
mpc_set_d_d(u[j + i__ * u_dim1], 0e0, 0e0, MPFR_RNDN);
mpc_set_fr(u[j + i__ * u_dim1], tmp0, MPFR_RNDN);
	    }
L10180:;
// printf("L10180:;\n");
L260:;
// printf("L260:;\n");
  if(np == *n) {
    goto L280;
  }
  if(n1 > np) {
    goto L10190;
  }
  for (j = n1; j <= np; ++j) {
    mpc_set(q,a[l1 + j * a_dim1], MPFR_RNDN);
    mpc_set(r__,a[i__ + j * a_dim1], MPFR_RNDN);
    mpc_mul_fr(q__2,q,cs, MPFR_RNDN);
    mpc_mul_fr(q__3,r__,sn, MPFR_RNDN);
    mpc_add(a[l1 + j * a_dim1],q__2,q__3, MPFR_RNDN);
/* L270: */
    mpc_mul_fr(q__2,r__,cs, MPFR_RNDN);
    mpc_mul_fr(q__3,q,sn, MPFR_RNDN);
    mpc_sub(a[i__ + j * a_dim1],q__2,q__3, MPFR_RNDN);
	    }
L10190:;
// printf("L10190:;\n");
L280:;
// printf("L280:;\n");
	    ;
	}
L10170:;
// printf("L10170:;\n");

	/*    TEST FOR CONVERGENCE */
L290:;
// printf("L290:;\n");
mpfr_set(w, s[k], MPFR_RNDN);
	if(l == k) {
	    goto L360;
	}

	/*    ORIGIN SHIFT */
mpfr_set(x, s[l], MPFR_RNDN);
mpfr_set(y, s[k - 1], MPFR_RNDN);
mpfr_set(g, t[k - 2], MPFR_RNDN);
mpfr_set(h__, t[k - 1], MPFR_RNDN);
// printf("h__ 2\n");
// mpfr_out_str(stdout, 10, 0, h__, MPFR_RNDN);
// printf("\n");

mpfr_mul(tmp0, y, y, MPFR_RNDN);
mpfr_mul(tmp1, w, w, MPFR_RNDN);
mpfr_sub(f, tmp0, tmp1, MPFR_RNDN);
mpfr_mul(tmp0, g, g, MPFR_RNDN);
mpfr_mul(tmp1, h__, h__, MPFR_RNDN);
mpfr_add(f, f, tmp0, MPFR_RNDN);
mpfr_sub(f, f, tmp1, MPFR_RNDN);
mpfr_div(tmp0, f, h__, MPFR_RNDN);
// printf("tmp0\n");
// mpfr_out_str(stdout, 10, 0, tmp0, MPFR_RNDN);
// printf("\n");

mpfr_div_d(tmp0, tmp0, 2e0, MPFR_RNDN);
// printf("tmp0\n");
// mpfr_out_str(stdout, 10, 0, tmp0, MPFR_RNDN);
// printf("\n");

mpfr_div(f, tmp0, y, MPFR_RNDN);
// printf("f\n");
// mpfr_out_str(stdout, 10, 0, f, MPFR_RNDN);
// printf("\n");

mpfr_mul(cmp0, f, f, MPFR_RNDN);
mpfr_add_d(cmp0, cmp0, 1e0, MPFR_RNDN);
mpfr_sqrt(g, cmp0, MPFR_RNDN);
// printf("f and g before neg\n");
// mpfr_out_str(stdout, 10, 0, f, MPFR_RNDN);
// printf("\n");

// mpfr_out_str(stdout, 10, 0, g, MPFR_RNDN);
// printf("\n");

	if(mpfr_less_p(f, mpfr_zero)) {
// mpfr_set(g, g, MPFR_RNDN);
mpfr_neg(g, g, MPFR_RNDN);
	}
mpfr_add(f, f, g, MPFR_RNDN);
mpfr_div(f, y, f, MPFR_RNDN);
// printf("f\n");
// mpfr_out_str(stdout, 10, 0, f, MPFR_RNDN);
// printf("\n");

mpfr_sub(f, f, h__, MPFR_RNDN);
mpfr_mul(f, f, h__, MPFR_RNDN);
mpfr_mul(tmp0, x, x, MPFR_RNDN);
mpfr_mul(tmp1, w, w, MPFR_RNDN);
mpfr_add(f, f, tmp0, MPFR_RNDN);
mpfr_sub(f, f, tmp1, MPFR_RNDN);
mpfr_div(f, f, x, MPFR_RNDN);
// printf("f\n");
// mpfr_out_str(stdout, 10, 0, f, MPFR_RNDN);
// printf("\n");



	/*    QR STEP */
mpfr_set_d(cs, 1e0, MPFR_RNDN);
mpfr_set_d(sn, 1e0, MPFR_RNDN);
l1 = l + 1;
	if(l1 > k) {
	    goto L10200;
	}
	for (i__ = l1; i__ <= k; ++i__) {
mpfr_set(g, t[i__ - 1], MPFR_RNDN);
mpfr_set(y, s[i__], MPFR_RNDN);
mpfr_mul(h__, sn, g, MPFR_RNDN);
// printf("h__ 3\n");
// mpfr_out_str(stdout, 10, 0, h__, MPFR_RNDN);
// printf("\n");

mpfr_mul(g, cs, g, MPFR_RNDN);
mpfr_mul(tmp0, h__, h__, MPFR_RNDN);
mpfr_mul(tmp1, f, f, MPFR_RNDN);
mpfr_add(w, tmp0, tmp1, MPFR_RNDN);
mpfr_sqrt(w, w, MPFR_RNDN);
mpfr_set(t[i__ - 2], w, MPFR_RNDN);

// FT modification
if(mpfr_equal_p(w, mpfr_zero)) {
	mpfr_set_d(w,1,MPFR_RNDN);
	mpfr_set_d(cs,1,MPFR_RNDN);
	mpfr_set_d(sn,0,MPFR_RNDN);
 } else {
	mpfr_div(cs, f, w, MPFR_RNDN);
	mpfr_div(sn, h__, w, MPFR_RNDN);
 }


mpfr_mul(tmp0, x, cs, MPFR_RNDN);
mpfr_mul(tmp1, g, sn, MPFR_RNDN);
mpfr_add(f, tmp0, tmp1, MPFR_RNDN);
mpfr_mul(tmp0, g, cs, MPFR_RNDN);
mpfr_mul(tmp1, x, sn, MPFR_RNDN);
mpfr_sub(g, tmp0, tmp1, MPFR_RNDN);
mpfr_mul(h__, y, sn, MPFR_RNDN);
// printf("h__ 4\n");
// mpfr_out_str(stdout, 10, 0, h__, MPFR_RNDN);
// printf("\n");

mpfr_mul(y, y, cs, MPFR_RNDN);
	    if(*nv == 0) {
		goto L310;
	    }
	    if(1 > *n) {
		goto L10210;
	    }
	    for (j = 1; j <= *n; ++j) {
mpfr_set(x, mpc_realref(v[j +  (i__ - 1) * v_dim1]), MPFR_RNDN);
mpfr_set(w, mpc_realref(v[j + i__ * v_dim1]), MPFR_RNDN);
mpfr_mul(tmp0, x, cs, MPFR_RNDN);
mpfr_mul(tmp1, w, sn, MPFR_RNDN);
mpfr_add(tmp0, tmp0, tmp1, MPFR_RNDN);
mpc_set_d_d(v[j +  (i__ - 1) * v_dim1], 0e0, 0e0, MPFR_RNDN);
mpc_set_fr(v[j +  (i__ - 1) * v_dim1], tmp0, MPFR_RNDN);
/* L300: */
mpfr_mul(tmp0, w, cs, MPFR_RNDN);
mpfr_mul(tmp1, x, sn, MPFR_RNDN);
mpfr_sub(tmp0, tmp0, tmp1, MPFR_RNDN);
mpc_set_d_d(v[j +  i__ * v_dim1], 0e0, 0e0, MPFR_RNDN);
mpc_set_fr(v[j +  i__ * v_dim1], tmp0, MPFR_RNDN);
	    }
L10210:;
// printf("L10210:;\n");
L310:;
// printf("L310:;\n");

// printf("h__\n");
// mpfr_out_str(stdout, 10, 0, h__, MPFR_RNDN);
// printf("\n");

// printf("f\n");
// mpfr_out_str(stdout, 10, 0, f, MPFR_RNDN);
// printf("\n");

mpfr_mul(tmp0, h__, h__, MPFR_RNDN);
mpfr_mul(tmp1, f, f, MPFR_RNDN);
mpfr_add(w, tmp0, tmp1, MPFR_RNDN);
mpfr_sqrt(w, w, MPFR_RNDN);
mpfr_set(s[i__ - 1], w, MPFR_RNDN);

// FT modification
if(mpfr_equal_p(w, mpfr_zero)) {
	mpfr_set_d(w,1,MPFR_RNDN);
	mpfr_set_d(cs,1,MPFR_RNDN);
	mpfr_set_d(sn,0,MPFR_RNDN);
} else {
	mpfr_div(cs, f, w, MPFR_RNDN);
	mpfr_div(sn, h__, w, MPFR_RNDN);
 }

mpfr_mul(tmp0, cs, g, MPFR_RNDN);
mpfr_mul(tmp1, sn, y, MPFR_RNDN);
mpfr_add(f, tmp0, tmp1, MPFR_RNDN);
mpfr_mul(tmp0, cs, y, MPFR_RNDN);
mpfr_mul(tmp1, sn, g, MPFR_RNDN);
mpfr_sub(x, tmp0, tmp1, MPFR_RNDN);
	    if(*nu == 0) {
		goto L330;
	    }
	    if(1 > *m) {
		goto L10220;
	    }
	    for (j = 1; j <= *m; ++j) {
mpfr_set(y, mpc_realref(u[j +  (i__ - 1) * u_dim1]), MPFR_RNDN);
mpfr_set(w, mpc_realref(u[j + i__ * u_dim1]), MPFR_RNDN);
mpfr_mul(tmp0, y, cs, MPFR_RNDN);
mpfr_mul(tmp1, w, sn, MPFR_RNDN);
mpfr_add(tmp0, tmp0, tmp1, MPFR_RNDN);
mpc_set_d_d(u[j +  (i__ - 1) * u_dim1], 0e0, 0e0, MPFR_RNDN);
mpc_set_fr(u[j +  (i__ - 1) * u_dim1], tmp0, MPFR_RNDN);
/* L320: */
mpfr_mul(tmp0, w, cs, MPFR_RNDN);
mpfr_mul(tmp1, y, sn, MPFR_RNDN);
mpfr_sub(tmp0, tmp0, tmp1, MPFR_RNDN);
mpc_set_d_d(u[j +  i__ * u_dim1], 0e0, 0e0, MPFR_RNDN);
mpc_set_fr(u[j +  i__ * u_dim1], tmp0, MPFR_RNDN);
	    }
L10220:;
// printf("L10220:;\n");
L330:;
// printf("L330:;\n");
	    if(*n == np) {
		goto L350;
	    }
	    if(n1 > np) {
		goto L10230;
	    }
	    for (j = n1; j <= np; ++j) {
mpc_set(q, a[i__ - 1 + j * a_dim1], MPFR_RNDN);
mpc_set(r__, a[i__ + j * a_dim1], MPFR_RNDN);
mpc_mul_fr(a[i__ - 1 + j * a_dim1], q, cs, MPFR_RNDN);
mpc_mul_fr(q__1, r__, sn, MPFR_RNDN);
mpc_add(a[i__ - 1 + j * a_dim1], a[i__ - 1 + j * a_dim1], q__1, MPFR_RNDN);
/* L340: */
mpc_mul_fr(a[i__ + j * a_dim1], r__, cs, MPFR_RNDN);
mpc_mul_fr(q__1, q, sn, MPFR_RNDN);
mpc_sub(a[i__ + j * a_dim1], a[i__ + j * a_dim1], q__1, MPFR_RNDN);
	    }
L10230:;
// printf("L10230:;\n");
L350:;
// printf("L350:;\n");
	    ;
	}
L10200:;
// printf("L10200:;\n");
mpfr_set_d(t[l - 1], 0e0, MPFR_RNDN);
mpfr_set(t[k - 1], f, MPFR_RNDN);
mpfr_set(s[k], x, MPFR_RNDN);
	goto L220;

	/*  CONVERGENCE */
L360:;
// printf("L360:;\n");
	if(mpfr_greaterequal_p(w, mpfr_zero)) {
	    goto L380;
	}
mpfr_set(s[k], w, MPFR_RNDN);
mpfr_neg(s[k], s[k], MPFR_RNDN);
	if(*nv == 0) {
	    goto L380;
	}
	if(1 > *n) {
	    goto L10240;
	}
	for (j = 1; j <= *n; ++j) {
/* L370: */
	  mpc_neg(v[j + k * v_dim1], v[j + k * v_dim1], MPFR_RNDN);
	}
L10240:;
// printf("L10240:;\n");
L380:;
// printf("L380:;\n");
	;
    }
L10150:;
// printf("L10150:;\n");

	/* SORT SINGULAR VALUES */
    if(1 > *n) {
	goto L10250;
    }
    for (k = 1; k <= *n; ++k) {
mpfr_set(g, s[k], MPFR_RNDN);
j = k;
	if(k > *n) {
	    goto L10260;
	}
	for (i__ = k; i__ <= *n; ++i__) {
	    if(mpfr_greaterequal_p(s[i__], g)) {
		goto L390;
	    }
mpfr_set(g, s[i__], MPFR_RNDN);
j = i__;
L390:;
// printf("L390:;\n");
	    ;
	}
L10260:;
// printf("L10260:;\n");
	if(j == k) {
	    goto L450;
	}
// printf("j, k = %d, %d\n", j, k);
mpfr_set(s[j], s[k], MPFR_RNDN);
mpfr_set(s[k], g, MPFR_RNDN);
	if(*nv == 0) {
	    goto L410;
	}
	if(1 > *n) {
	    goto L10270;
	}
	for (i__ = 1; i__ <= *n; ++i__) {
mpc_set(q, v[i__ + j * v_dim1], MPFR_RNDN);
mpc_set(v[i__ + j * v_dim1], v[i__ + k * v_dim1], MPFR_RNDN);
/* L400: */
mpc_set(v[i__ + k * v_dim1], q, MPFR_RNDN);
	}
L10270:;
// printf("L10270:;\n");
L410:;
// printf("L410:;\n");
	if(*nu == 0) {
	    goto L430;
	}
	if(1 > *m) {
	    goto L10280;
	}
	for (i__ = 1; i__ <= *m; ++i__) {
// printf("i__, j, u_dim1 = %d, %d, %d\n", i__, j, u_dim1);
mpc_set(q, u[i__ + j * u_dim1], MPFR_RNDN);
mpc_set(u[i__ + j * u_dim1], u[i__ + k * u_dim1], MPFR_RNDN);
/* L420: */
mpc_set(u[i__ + k * u_dim1], q, MPFR_RNDN);
	}
L10280:;
// printf("L10280:;\n");
L430:;
// printf("L430:;\n");
	if(*n == np) {
	    goto L450;
	}
	if(n1 > np) {
	    goto L10290;
	}
	for (i__ = n1; i__ <= np; ++i__) {
mpc_set(q, a[j + i__ * a_dim1], MPFR_RNDN);
mpc_set(a[j + i__ * a_dim1], a[k + i__ * a_dim1], MPFR_RNDN);
/* L440: */
mpc_set(a[k + i__ * a_dim1], q, MPFR_RNDN);
	}
L10290:;
// printf("L10290:;\n");
L450:;
// printf("L450:;\n");
	;
    }
L10250:;
// printf("L10250:;\n");

	/* BACK TRANSFORMATION */
    if(*nu == 0) {
	goto L510;
    }
    if(1 > *n) {
	goto L10300;
    }
    for (kk = 1; kk <= *n; ++kk) {
k = n1 - kk;
	if(mpfr_equal_p(b[k - 1], mpfr_zero)) {
	    goto L500;
	}
	mpc_neg(q, a[k + k * a_dim1], MPFR_RNDN);
mpc_abs(r__1, a[k + k * a_dim1], MPFR_RNDN);
	mpc_div_fr(q,q,r__1, MPFR_RNDN);
// printf("q\n");
// mpc_out_str(stdout, 10, 0, q, MPFR_RNDN);
// printf("\n");

	if(1 > *nu) {
	    goto L10310;
	}
	for (j = 1; j <= *nu; ++j) {
/* L460: */
	  mpc_mul(u[k + j * u_dim1],q,u[k + j * u_dim1], MPFR_RNDN);
	}
L10310:;
// printf("L10310:;\n");
	if(1 > *nu) {
	    goto L10320;
	}
	for (j = 1; j <= *nu; ++j) {
	  mpc_set_d_d(q,0e0,0e0, MPFR_RNDN);
	    if(k > *m) {
		goto L10330;
	    }
	    for (i__ = k; i__ <= *m; ++i__) {
/* L470: */
	      mpc_conj(q__3, a[i__ + k * a_dim1], MPFR_RNDN);
	      mpc_mul(q__3,q__3,u[i__ + j * u_dim1], MPFR_RNDN);
	      mpc_add(q,q,q__3, MPFR_RNDN);
	    }
L10330:;
// printf("L10330:;\n");
mpc_abs(tmp0, a[k + k * a_dim1], MPFR_RNDN);
mpfr_mul(r__1, tmp0, b[k - 1], MPFR_RNDN);
	    mpc_div_fr(q,q,r__1, MPFR_RNDN);
// printf("q\n");
// mpc_out_str(stdout, 10, 0, q, MPFR_RNDN);
// printf("\n");

	    if(k > *m) {
		goto L10340;
	    }
	    for (i__ = k; i__ <= *m; ++i__) {
/* L480: */
	      mpc_mul(q__2,q,a[i__ + k * a_dim1], MPFR_RNDN);
	      mpc_sub(u[i__ + j * u_dim1],u[i__ + j * u_dim1],q__2, MPFR_RNDN);
	    }
L10340:;
// printf("L10340:;\n");
/* L490: */
	    ;
	}
L10320:;
// printf("L10320:;\n");
L500:;
// printf("L500:;\n");
	;
    }
L10300:;
// printf("L10300:;\n");
L510:;
// printf("L510:;\n");
    if(*nv == 0) {
	goto L570;
    }
    if(*n < 2) {
	goto L570;
    }
    if(2 > *n) {
	goto L10350;
    }
    for (kk = 2; kk <= *n; ++kk) {
k = n1 - kk;
k1 = k + 1;
	if(mpfr_equal_p(c__[k1 - 1], mpfr_zero)) {
	    goto L560;
	}
	mpc_conj(q, a[k + k1 * a_dim1], MPFR_RNDN);
	mpc_neg(q, q, MPFR_RNDN);
mpc_abs(r__1, a[k + k1 * a_dim1], MPFR_RNDN);
	mpc_div_fr(q,q,r__1, MPFR_RNDN);
// printf("q\n");
// mpc_out_str(stdout, 10, 0, q, MPFR_RNDN);
// printf("\n");

	if(1 > *nv) {
	    goto L10360;
	}
	for (j = 1; j <= *nv; ++j) {
/* L520: */
	  mpc_mul(v[k1 + j * v_dim1],q,v[k1 + j * v_dim1], MPFR_RNDN);
	}
L10360:;
// printf("L10360:;\n");
	if(1 > *nv) {
	    goto L10370;
	}
	for (j = 1; j <= *nv; ++j) {
	  mpc_set_d_d(q,0e0,0e0, MPFR_RNDN);
	    if(k1 > *n) {
		goto L10380;
	    }
	    for (i__ = k1; i__ <= *n; ++i__) {
/* L530: */
	      mpc_mul(q__2,a[k + i__ * a_dim1],v[i__ + j * v_dim1], MPFR_RNDN);
	      mpc_add(q,q,q__2, MPFR_RNDN);
	    }
L10380:;
// printf("L10380:;\n");
mpc_abs(r__1, a[k + k1 * a_dim1], MPFR_RNDN);
mpfr_mul(r__1, r__1, c__[k1 - 1], MPFR_RNDN);
	    mpc_div_fr(q,q,r__1, MPFR_RNDN);
// printf("q\n");
// mpc_out_str(stdout, 10, 0, q, MPFR_RNDN);
// printf("\n");

	    if(k1 > *n) {
		goto L10390;
	    }
	    for (i__ = k1; i__ <= *n; ++i__) {
/* L540: */
	      mpc_conj(q__3, a[k + i__ * a_dim1], MPFR_RNDN);
	      mpc_mul(q__3,q,q__3, MPFR_RNDN);
	      mpc_sub(v[i__ + j * v_dim1],v[i__ + j * v_dim1],q__3, MPFR_RNDN);
	    }
L10390:;
// printf("L10390:;\n");
/* L550: */
	    ;
	}
L10370:;
// printf("L10370:;\n");
L560:;
// printf("L560:;\n");
	;
    }
L10350:;
// printf("L10350:;\n");
L570:;
// printf("L570:;\n");
    return 0;
} /* csvd1_ */


int rdefl_(int *nm, int *n, int *kb, int *ndef, mpfr_t *hr, mpfr_t *hi, mpfr_t *zr, mpfr_t *zi, mpfr_t *told, int *d__, mpfr_t *c1, mpfr_t *c2, mpfr_t *dl, int last_call, int reallc)
{	
	
	mpfr_t mpfr_zero, mpfr_one;
	mpfr_init2(mpfr_zero, wp2);
	mpfr_set_d(mpfr_zero, 0, MPFR_RNDN);
	mpfr_init2(mpfr_one, wp2);
	mpfr_set_d(mpfr_one, 1, MPFR_RNDN);
    /* Format strings */

    /* System generated locals */
    int hr_dim1, hr_offset, hi_dim1, hi_offset, zr_dim1, zr_offset, zi_dim1, zi_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    mpfr_t r__1;
    mpfr_init2(r__1, wp2);
    mpc_t q__1;
    mpc_init3(q__1, wp2, wp2);
    mpc_t q__2;
    mpc_init3(q__2, wp2, wp2);
    mpc_t q__3;
    mpc_init3(q__3, wp2, wp2);
    mpc_t q__4;
    mpc_init3(q__4, wp2, wp2);

	int nnm = *n * *nm;
	int N = *n;

	// printf("nm = %d\n", *nm);
	// printf("n = %d\n", *n);
	// printf("kb = %d\n", *kb);
	// printf("ndef = %d\n", *ndef);

    /* Builtin functions */
    
    /* Local variables */
	/* was [20][20] */
    static mpc_t *h__ = NULL;
	static mpfr_t *s = NULL;
	static mpc_t *u = NULL;
	static mpc_t *v = NULL;

	// if (first_call) {
	// 	h__ = NULL;
	// 	s = NULL;
	// 	u = NULL;
	// 	v = NULL;
	// }
	
	if (reallc) {
		// printf("reallc option activated\n");
		if (h__) {
			// printf("free h__\n");
			free(h__);
			h__ = NULL;
		}
	}
	if(!h__) {
		// printf("allocating h__ with nnm = %d\n", nnm);
		h__ = (mpc_t*) malloc(nnm * sizeof(mpc_t));
		for (int i=0; i<nnm; i++) {
			mpc_init3(h__[i], wp2, wp2);
		}
	}
	 
    static int i__, j, k;
    
	if (reallc) {
		// printf("reallc option activated\n");
		if (s) {
			// printf("free s\n");
			free(s);
			s = NULL;
		}
	}
	if(!s) {
		// printf("allocating s with N = %d\n", N);
		s = (mpfr_t*) malloc(N * sizeof(mpfr_t));
		for (int i=0; i<N; i++) {
			mpfr_init2(s[i], wp2);
		}
	}
	/* was [20][20] */
    
	if (reallc) {
		// printf("reallc option activated\n");
		if (u) {
			// printf("free u\n");
			free(u);
			u = NULL;
		}
	}
	if(!u) {
		// printf("allocating u with nnm = %d\n", nnm);
		u = (mpc_t*) malloc(nnm * sizeof(mpc_t));
		for (int i=0; i<nnm; i++) {
			mpc_init3(u[i], wp2, wp2);
		}
	}
    
    
	if (reallc) {
		// printf("reallc option activated\n");
		if (v) {
			// printf("free v\n");
			free(v);
			v = NULL;
		}
	}
	if(!v) {
		// printf("allocating v with nnm = %d\n", nnm);
		v = (mpc_t*) malloc(nnm * sizeof(mpc_t));
		for (int i=0; i<nnm; i++) {
			mpc_init3(v[i], wp2, wp2);
		}
	}

    static mpfr_t t1;
    mpfr_init2(t1, wp2);
    static mpfr_t t2;
    mpfr_init2(t2, wp2);
    static int id;
    static mpc_t hij;
    mpc_init3(hij, wp2, wp2);
    static mpc_t zij;
    mpc_init3(zij, wp2, wp2);
    static int mmax;

    /* Fortran I/O blocks */



	/*   ******************************************************************** */
	/*   *    THE  ROUTINE  MAKES  A  CSVD- DECOMPOSITION OF A BLOCK TO */
	/*   *    DETERMINE THE RANK OF THE BLOCK */
	/*   * */
	/*   *    THE FORMAL PARAMETER LIST */
	/*   * */
	/*   *    ON INPUT */
	/*   * */
	/*   *    NM        THE ROW DIMENSION OF THE TWO-DIMENSIONAL ARRAY */
	/*   *              PARAMETERS(HR,HI,ZR,ZI) AS DECLARED IN THE */
	/*   *              CALLING PROGRAM DIMENSION STATEMENT */
	/*   *    N         THE ORDER OF THE ORIGINAL MATRIX */
	/*   *    KB        THE DIMENSION OF THE BLOCK STARTING AT THE */
	/*   *              ADDRESS HR,HI */
	/*   *    HR,HI     THE ADDRESS OF THE REAL AND IMAGINARY PARTS, */
	/*   *              RESPECTIVELY,OF THE ACTUAL BLOCK S (1,1)- */
	/*   *              ELEMENT.THE BLOCK IS AFFECTED BY THE PROCESS */
	/*   *    ZR,ZI     THE REAL AND IMAGINARY PARTS,RESPECTIVELY, */
	/*   *              OF THE TRANSFORMATION MATRIX(WILL BE ACCUMULATED */
	/*   *              DURING THE RANK DETERMINATION PROCESS) */
	/*   *    C1,C2     TOLERANCE GAP IS C2/C1 IN THE DEFLATION PROCESS */
	/*   *              (THE PARAMETERS ARE GIVEN VALUES BY JNF,C1=1.0 */
	/*   *               AND C2=1000.0) */
	/*   *    TOLD      TOLERANCE PARAMETER USED IN THE DETERMINATION */
	/*   *              OF THE RANK OF THE ACTUAL BLOCK */
	/*   *              (TOLD=TOL*THE FOBENIUS NORM OF THE BALANCED */
	/*   *              MATRIX FROM CBAL,THE PARAMETER IS COMPUTED BY JNF) */
	/*   * */
	/*   *    ON OUTPUT */
	/*   * */
	/*   *    NDEF      THE NUMBER OF DEFLATION STEPS I.E THE NUMBER OF */
	/*   *              SINGULAR VALUES ,INTERPRETED AS ZEROES */
	/*   *    DL        THE SQUARE OF DELETED SINGULAR VALUES */
	/*   * */
	/*   *    THE DIMENSIONS OF H,U,V AND S HAVE TO BE CHANGED FOR */
	/*   *    MATRICES HR,HI WITH DIMENSION GREATER THEN 20 */
	/*   ******************************************************************** */

    /* Parameter adjustments */
zi_dim1 = *nm;
zi_offset = 1 + zi_dim1;
zi = zi - zi_offset;
zr_dim1 = *nm;
zr_offset = 1 + zr_dim1;
zr = zr - zr_offset;
hi_dim1 = *nm;
hi_offset = 1 + hi_dim1;
hi = hi - hi_offset;
hr_dim1 = *nm;
hr_offset = 1 + hr_dim1;
hr = hr - hr_offset;

    /* Function Body */
mmax = *nm;
int reallc_csvd = 1;
// printf("INPUT CSVD\nmatrix h__\n");
    for (i__ = 1; i__ <= *kb; ++i__) {
	for (j = 1; j <= *kb; ++j) {
	  // printf("mm[[%d,%d]]=(",i__,j);
	  // mpfr_out_str(stdout, 10, 0, hr[i__ + j * hr_dim1], MPFR_RNDN);
	  // printf(", sign");
	  // mpfr_out_str(stdout, 10, 0, hi[i__ + j * hr_dim1], MPFR_RNDN);
	  // printf("*I);\n");

	  // printf("(i__, j) = (%d, %d)\n", i__, j);
	  // printf("index = %d\n", i__ + j * *nm - *nm - 1);
		// fflush(stdout);
	  mpc_set_fr_fr(h__[i__ + j * *nm - *nm - 1],hr[i__ + j * hr_dim1],hi[i__ + j * hi_dim1], MPFR_RNDN);
	  mpc_set_fr_fr(h__[i__ + j * *nm - *nm - 1],hr[i__ + j * hr_dim1],hi[i__ + j * hi_dim1], MPFR_RNDN);
	//   mpc_out_str(stdout, 10, 0, h__[i__ + j * *nm - *nm - 1], MPFR_RNDN);
	//   printf("\n");
	  /* L1: */
	}
    }
	/*   ******  CSVD1 IS THE CACM ALGORITHM 358 WITH ONE CHANCE, */
	/*           NAMELY,THE SINGULAR VALUES ARE SORTED IN INCREASING */
	/*           ORDER  (INSTEAD OF DECREASING ORDER). */
	// printf("kb = %d\n", *kb);
	// getchar();

		csvd1_(h__, &mmax, &mmax, kb, kb, &c__0, kb, kb, s, u, v, *n, reallc_csvd);
		reallc_csvd = 0;
    // printf("vector s\n");
	// printf("{");
    // for (int i=0; i<*kb; i++) {
	// 	mpfr_out_str(stdout, 10, 0, s[i], MPFR_RNDN);
	// 	printf("\n");
    // }
    // printf("}\n");
	//getchar();
	/*   * */
	/*   *    PRINT THE SINGULAR VALUES TO GET THE THEORETICAL */
	/*   *    COUPLING ELEMENT */

	/*   *****  DETERMINE THE NUMBER OF DEFLATION STEPS NDEF= NUMBER OF */
	/*          SINGULAR VALUES EQUAL TO ZERO                             *** */
mpfr_mul(t1, *c1, *told, MPFR_RNDN);
mpfr_mul(t2, *c2, *told, MPFR_RNDN);
    for (i__ = 1; i__ <= *kb; ++i__) {
mpfr_abs(r__1, s[i__ - 1], MPFR_RNDN);
	// printf("t1 = "); mpfr_out_str(stdout, 10, 0, t1, MPFR_RNDN); printf("\n");
	// printf("r__1 = "); mpfr_out_str(stdout, 10, 0, r__1, MPFR_RNDN); printf("\n");
	// printf("s[i__-1]= "); mpfr_out_str(stdout, 10, 0, s[i__-1], MPFR_RNDN); printf("\n");
	if(mpfr_greaterequal_p(r__1, t1)) {
		// printf("L4\n");
	    goto L4;
	}
/* L2: */
    }
*ndef = *kb;
    for (i__ = 1; i__ <= *ndef; ++i__) {
/* L3: */
	/* Computing 2nd power */
mpfr_mul(r__1, s[i__ - 1], s[i__ - 1], MPFR_RNDN);
mpfr_add(*dl, *dl, r__1, MPFR_RNDN);
    }
	// free memory if this is the last call of rdefl
	if (last_call) {
		// printf("last call: free static vectors h__, s, u, v\n");
		free(h__);
		free(s);
		free(u);
		free(v);
		h__ = NULL;
		s = NULL;
		u = NULL;
		v = NULL;
	}
    return 0;
L4:;
// printf("L4:;\n");
mpfr_abs(r__1, s[i__ - 1], MPFR_RNDN);
    if(mpfr_greater_p(r__1, t2)) {
	goto L5;
    }
    if(i__ >= *kb) {
	goto L5;
    }
    ++i__;
    goto L4;
L5:;
// printf("L5:;\n");
*ndef = i__ - 1;
    if(*ndef == 0) {
*ndef = 1;
    }
    for (i__ = 1; i__ <= *ndef; ++i__) {
	/* Computing 2nd power */
mpfr_mul(r__1, s[i__ - 1], s[i__ - 1], MPFR_RNDN);
mpfr_add(*dl, *dl, r__1, MPFR_RNDN);
/* L6: */
    }

	/*   *****   OVERWRITE  HR,HI WITH  V(TRANS,CONJG)*U* DIAG(S(I))      *** */
    for (i__ = 1; i__ <= *kb; ++i__) {
	for (j = 1; j <= *kb; ++j) {
	  mpc_set_d_d(hij, 0e0, 0e0, MPFR_RNDN);
	    for (k = 1; k <= *kb; ++k) {
/* L7: */
		mpc_conj(q__3, v[k + i__ * *nm - *nm - 1], MPFR_RNDN);
		mpc_mul(q__3,q__3,u[k + j * *nm - *nm - 1], MPFR_RNDN);
		mpc_add(hij,hij,q__3, MPFR_RNDN);
	    }
	    mpc_mul_fr(hij,hij,s[j - 1], MPFR_RNDN);
mpfr_set(hr[i__ + j * hr_dim1], mpc_realref(hij), MPFR_RNDN);
/* L8: */
mpfr_set(hi[i__ + j * hi_dim1], mpc_imagref(hij), MPFR_RNDN);
	}
    }

	/*   *****  MULTIPLY COLUMNS NF,...,NF+ID WITH V(TRANS,CONJG) FROM LEFT * */
    if(1 > *d__) {
	goto L10060;
    }
    for (j = 1; j <= *d__; ++j) {
id = j - *d__;
	for (i__ = 1; i__ <= *kb; ++i__) {
	  mpc_set_d_d(hij, 0e0, 0e0, MPFR_RNDN);
	  for (k = 1; k <= *kb; ++k) {
/* L9: */
		mpc_conj(q__3, v[k + i__ * *nm - *nm - 1], MPFR_RNDN);
		mpc_set_fr_fr(q__4,hr[k + id * hr_dim1],hi[k + id * hi_dim1], MPFR_RNDN);
		mpc_mul(q__3,q__3,q__4, MPFR_RNDN);
		mpc_add(hij,hij,q__3, MPFR_RNDN);
	    }
/* L10: */
mpc_set(h__[i__ * *nm - *nm], hij, MPFR_RNDN);
	}
	for (i__ = 1; i__ <= *kb; ++i__) {
mpfr_set(hr[i__ + id * hr_dim1], mpc_realref(h__[i__ * *nm - *nm]), MPFR_RNDN);
/* L11: */
mpfr_set(hi[i__ + id * hi_dim1], mpc_imagref(h__[i__ * *nm - *nm]), MPFR_RNDN);
	}
/* L12: */
    }
L10060:;
// printf("L10060:;\n");

	/*   *****  FINISH THE SIMILARITY TRANSFORMATION AND ACCUMULATE */
	/*          TRANSFORMATIONS                                           *** */
    for (i__ = 1; i__ <= *n; ++i__) {
id = i__ - *d__;
      for (j = 1; j <= *kb; ++j) {
	mpc_set_d_d(hij,0e0,0e0, MPFR_RNDN);
	mpc_set_d_d(zij,0e0,0e0, MPFR_RNDN);
	for (k = 1; k <= *kb; ++k) {
	  if(id > 0) {
	    goto L13;
	  }
	  mpc_set_fr_fr(q__3,hr[id + k * hr_dim1],hi[id + k * hr_dim1], MPFR_RNDN);
	  mpc_mul(q__3,q__3,v[k + j * *nm - *nm - 1], MPFR_RNDN);
	  mpc_add(hij,hij,q__3, MPFR_RNDN);
L13:;
// printf("L13:;\n");
/* L14: */
	  mpc_set_fr_fr(q__3,zr[i__ + k * zr_dim1],zi[i__ + k * zi_dim1], MPFR_RNDN);
	  mpc_mul(q__3,q__3,v[k + j * *nm - *nm - 1], MPFR_RNDN);
	  mpc_add(zij,zij,q__3, MPFR_RNDN);
	}
	if(id > 0) {
	  goto L15;
	    }
	mpc_set(h__[j * *nm - *nm],hij, MPFR_RNDN);
L15:;
// printf("L15:;\n");
	mpc_set(u[j * *nm - *nm],zij, MPFR_RNDN);
/* L16: */
      }
      for (j = 1; j <= *kb; ++j) {
	if(id > 0) {
	  goto L17;
	}
mpfr_set(hr[id + j * hr_dim1], mpc_realref(h__[j * *nm - *nm]), MPFR_RNDN);
mpfr_set(hi[id + j * hi_dim1], mpc_imagref(h__[j * *nm - *nm]), MPFR_RNDN);
L17:;
// printf("L17:;\n");
mpfr_set(zr[i__ + j * zr_dim1], mpc_realref(u[j * *nm - *nm]), MPFR_RNDN);
mpfr_set(zi[i__ + j * zi_dim1], mpc_imagref(u[j * *nm - *nm]), MPFR_RNDN);
/* L18: */
	}
/* L19: */
    }

	// free memory if this is the last call of rdefl
	if (last_call) {
		free(h__);
		free(s);
		free(u);
		free(v);
		h__ = NULL;
		s = NULL;
		u = NULL;
		v = NULL;
	}
    return 0;
} /* rdefl_ */

// #TBD: with a 2x2 (exactly) null matrix as input, the output has a NaN in the i=1,j=1 element.
// a temporary fix has been implemented, with a check on the input to detect null input,
// but further investigation is needed.
void jordan(mpc_t **jor, mpc_t **t_mat, int **ein_pos, int *num_ein, mpfr_t *ar, mpfr_t *ai, int N, int nm, int debug)
{
	// check if input is null
	int found_non_zero = 0;
	for (int i=0; i<N; i++) {
		for (int j=0; j<nm; j++) {
			if (mpfr_zero_p(ar[i+j*nm]) != 1 || mpfr_zero_p(ai[i+j*nm]) != 1) {
				found_non_zero = 1;
				break;
			}
		}
		if (found_non_zero) {
			break;
		}
	}
	if (!found_non_zero) {
		// null input detected

		// set jordan form to zero, transformation matrix to identity,
		// number of eigenvalues to row dimension and
		// eigenvalue positions to {1, ..., row dimension}
		*num_ein = N;
		(*ein_pos) = malloc(N*sizeof(int));
		for (int i=0; i<N; i++) {
			(*ein_pos)[i] = i;
			for (int j=0; j<nm; j++) {
				mpc_set_ui(jor[i][j], 0, MPFR_RNDN);
				if (i == j) {
					mpc_set_ui(t_mat[i][j], 1, MPFR_RNDN);
				} else {
					mpc_set_ui(t_mat[i][j], 0, MPFR_RNDN);
				}
			}
		}
		return;
	}
	// end check

	mpfr_t mpfr_zero, mpfr_one;
	mpfr_init2(mpfr_zero, wp2);
	mpfr_set_d(mpfr_zero, (double) 0, MPFR_RNDN);
	mpfr_init2(mpfr_one, wp2);
	mpfr_set_d(mpfr_one, (double) 1, MPFR_RNDN);
	mpfr_t tmp0;
	mpfr_init2(tmp0, wp2);
	mpfr_t tmp1;
	mpfr_init2(tmp1, wp2);
	mpfr_t tmp2;
	mpfr_init2(tmp2, wp2);
	mpfr_t tmp3;
	mpfr_init2(tmp3, wp2);
	mpfr_t cmp0;
	mpfr_init2(cmp0, wp2);

    mpfr_t *rad = (mpfr_t*) malloc((N+1)*sizeof(mpfr_t));
    for (int i=0; i<N+1; i++) {
    	mpfr_init2(rad[i], wp2);
    }
    mpfr_t *radr = (mpfr_t*) malloc((N+1)*sizeof(mpfr_t));
    for (int i=0; i<N+1; i++) {
    	mpfr_init2(radr[i], wp2);
    }

    mpfr_t *sm = (mpfr_t*) malloc((N+1)*sizeof(mpfr_t));
    for (int i=0; i<N+1; i++) {
    	mpfr_init2(sm[i], wp2);
    }
    mpfr_t *dele = (mpfr_t*) malloc((N+1)*sizeof(mpfr_t));
    for (int i=0; i<N+1; i++) {
    	mpfr_init2(dele[i], wp2);
    }
    mpfr_t *supd = (mpfr_t*) malloc((N+1)*sizeof(mpfr_t));
    for (int i=0; i<N+1; i++) {
    	mpfr_init2(supd[i], wp2);
    }
    mpfr_t *evr = (mpfr_t*) malloc((N+1)*sizeof(mpfr_t));
    for (int i=0; i<N+1; i++) {
    	mpfr_init2(evr[i], wp2);
    }
    mpfr_t *evi = (mpfr_t*) malloc((N+1)*sizeof(mpfr_t));
    for (int i=0; i<N+1; i++) {
    	mpfr_init2(evi[i], wp2);
    }
    int *ndel = (int*) malloc((N+1)*sizeof(int));
    int *ndefl = (int*) malloc((N+1)*sizeof(int));
    int *ndb = (int*) malloc((N+1)*sizeof(int));
    int *nxt = (int*) malloc((N+1)*sizeof(int));

    mpfr_t *zr = (mpfr_t*) malloc((N+1)*(nm+1)*sizeof(mpfr_t));
    for (int i=0; i<(N+1)*(nm+1); i++) {
    	mpfr_init2(zr[i], wp2);
    }
    mpfr_t *zi = (mpfr_t*) malloc((N+1)*(nm+1)*sizeof(mpfr_t));
    for (int i=0; i<(N+1)*(nm+1); i++) {
    	mpfr_init2(zi[i], wp2);
    }
    

/************************************************/
    
    int low = 0, igh = 0;

    static int nblock;
    int hr_dim1, hr_offset, hi_dim1, hi_offset, zr_dim1, zr_offset, zi_dim1, zi_offset, i__1, i__2, i__3, i__4, i__5, i__6, i__7, i__8, i__9, i__10, i__11;
    mpfr_t r__1;
    mpfr_init2(r__1, wp2);
    mpfr_t r__2;
    mpfr_init2(r__2, wp2);
    mpfr_t r__3;
    mpfr_init2(r__3, wp2);
    mpfr_t r__4;
    mpfr_init2(r__4, wp2);
    // static int i__, j, k, l;
    static mpfr_t fn;
    mpfr_init2(fn, wp2);
    static mpfr_t told;
    mpfr_init2(told, wp2);
    int int__[N];
    int ierr;
    int j, j1, j2;
    int k2, kprim;
    mpfr_t min__;
    mpfr_init2(min__, wp2);
    mpfr_t sk;
    mpfr_init2(sk, wp2);
    mpfr_t ksi;
    mpfr_init2(ksi, wp2);
    mpfr_t myr;
    mpfr_init2(myr, wp2);
    mpfr_t myi;
    mpfr_init2(myi, wp2);
    mpfr_t nyr;
    mpfr_init2(nyr, wp2);
    mpfr_t nyi;
    mpfr_init2(nyi, wp2);
    mpfr_t sl;
    mpfr_init2(sl, wp2);
    mpfr_t hilj;
    mpfr_init2(hilj, wp2);
    mpfr_t hrlj;
    mpfr_init2(hrlj, wp2);
    mpfr_t hilj1;
    mpfr_init2(hilj1, wp2);
    mpfr_t hrlj1;
    mpfr_init2(hrlj1, wp2);
    mpc_t q__1;
    mpc_init3(q__1, wp2, wp2);
    mpc_t q__2;
    mpc_init3(q__2, wp2, wp2);
    mpc_t q__3;
    mpc_init3(q__3, wp2, wp2);
    mpc_t q__4;
    mpc_init3(q__4, wp2, wp2);
    mpc_t q__5;
    mpc_init3(q__5, wp2, wp2);
    mpc_t q__6;
    mpc_init3(q__6, wp2, wp2);
    mpc_t q__7;
    mpc_init3(q__7, wp2, wp2);

    // tolerance
    /*   *  PC1        E20.5  TOLERANCE PARAMETER,IS USED TO COMPUTE */
    /*   *                    EINF=PC1*MACEP. */
    /*   *                    EINF CORRESPONDS TO PERTURBATIONS OF */
    /*   *                    THE INPUT MATRIX.EINF IS USED IN THE */
    /*   *                    GROUPING OF NUMERICAL MULTIPLE EIGEN- */
    /*   *                    VALUES AND IS A PARAMETER TO JNF. */
    /*   *                    PC1 IS PROBLEM-DEPENDENT BUT MACHINE- */
    /*   *                    INDEPENDENT. */
    /*   *  PC2        E20.5  TOLERANCE PARAMETER,IS USED TO COMPUTE */
    /*   *                    TOL=PC2*MACHEP. */
    /*   *                    TOL IS USED IN THE CONSTRUCTION OF NIL- */
    /*   *                    POTENT MATRICES I.E IN THE DETERMINATION */
    /*   *                    OF THE STRUCTURE OF THE INPUT MATRIX AND */
    /*   *                    IS A PARAMETER TO JNF. */
    /*   *                    PC2 IS PROBLEM-DEPENDENT BUT MACHINE- */
    /*   *                    INDEPENDENT. */

    /*   *  MAIN      MACHEP    RELATIVE MACHINE PRECISION (F.P.) */
    static mpfr_t pc1;
    mpfr_init2(pc1, wp2);
    static mpfr_t pc2;
    mpfr_init2(pc2, wp2);
// mpfr_set_d(pc1, 10, MPFR_RNDN);
// mpfr_pow_si(pc1, pc1, 2, MPFR_RNDN);
	// mpfr_set_ui(pc1, 1, MPFR_RNDN);
	mpfr_set_d(pc1, 2, MPFR_RNDN);
	// mpfr_pow_si(pc1, pc1, 0.05*((int) wp2), MPFR_RNDN);
	// mpfr_pow_si(pc1, pc1, 0.1*((int) wp2), MPFR_RNDN);
	// mpfr_pow_si(pc1, pc1, 0.2*((int) wp2), MPFR_RNDN);
	// mpfr_pow_si(pc1, pc1, 0.495*((int) wp2), MPFR_RNDN);
	mpfr_pow_si(pc1, pc1, 0.55*((int) wp2), MPFR_RNDN);
	// mpfr_pow_si(pc1, pc1, 0.80*((int) wp2), MPFR_RNDN);

	/* 1e20f * sqrt(10); */
	// original value was 1e6, used 1e8 for a while,
	// then switched to 30% of the digits
	// mpfr_set_d(pc2, 10, MPFR_RNDN);
	// mpfr_pow_si(pc2, pc2, 8, MPFR_RNDN);
	// mpfr_pow_si(pc2, pc2, 20, MPFR_RNDN);
	mpfr_set_d(pc2, 2, MPFR_RNDN);
	// mpfr_pow_si(pc2, pc2, 0.20*((int) wp2), MPFR_RNDN);
	// mpfr_pow_si(pc2, pc2, 0.27*((int) wp2), MPFR_RNDN);
	mpfr_pow_si(pc2, pc2, 0.30*((int) wp2), MPFR_RNDN);
	// mpfr_pow_si(pc2, pc2, 0.50*((int) wp2), MPFR_RNDN);
	// mpfr_pow_si(pc2, pc2, 0.90*((int) wp2), MPFR_RNDN);

	/* pc1; */
    static mpfr_t machep;
    mpfr_init2(machep, wp2);
    mpfr_set_d(machep, 2, MPFR_RNDN);
	mpfr_pow_si(machep, machep, - ((int) wp2), MPFR_RNDN);
	// printf("external machep = ");
	// mpfr_out_str(stdout, 10, 0, machep, MPFR_RNDN);
	// printf("\n");
	//getchar();

	// mpfr_set_d(machep, 10, MPFR_RNDN);
	// mpfr_pow_si(machep, machep, -14, MPFR_RNDN);
    mpfr_t ein;
    mpfr_init2(ein, wp2);
    mpfr_t tol;
    mpfr_init2(tol, wp2);
    mpfr_t einf;
    mpfr_init2(einf, wp2);
mpfr_mul(ein, pc1, machep, MPFR_RNDN);
mpfr_mul(tol, pc2, machep, MPFR_RNDN);


	int nf = 1;

    int ne, nd, mb, nd2, kp1;
    mpfr_t mpr;
    mpfr_init2(mpr, wp2);
    mpfr_t mpi;
    mpfr_init2(mpi, wp2);
    mpfr_t rgmax;
    mpfr_init2(rgmax, wp2);
    mpfr_t cabs1;
    mpfr_init2(cabs1, wp2);
    mpfr_t rigmin;
    mpfr_init2(rigmin, wp2);
    mpfr_t radi;
    mpfr_init2(radi, wp2);
    mpfr_t shtr;
    mpfr_init2(shtr, wp2);
    mpfr_t shti;
    mpfr_init2(shti, wp2);
	static mpfr_t rsum;
	mpfr_init2(rsum, wp2);
    mpfr_t *dt = (mpfr_t*) malloc(N*sizeof(mpfr_t));
    for (int i=0; i<N; i++) {
    	mpfr_init2(dt[i], wp2);
    }
    mpfr_t *cshtr = (mpfr_t*) malloc(N*sizeof(mpfr_t));
    for (int i=0; i<N; i++) {
    	mpfr_init2(cshtr[i], wp2);
    }
    mpfr_t *cshti = (mpfr_t*) malloc(N*sizeof(mpfr_t));
    for (int i=0; i<N; i++) {
    	mpfr_init2(cshti[i], wp2);
    }

	int ns, ibl, ibl2, k;    
	int l1, i__, k5, i1, k1, k3, nss;
    mpc_t m;
    mpc_init3(m, wp2, wp2);
    mpc_t x;
    mpc_init3(x, wp2, wp2);
    mpfr_t xr;
    mpfr_init2(xr, wp2);
    mpfr_t xi;
    mpfr_init2(xi, wp2);
    int l;
    mpfr_t c1;
    mpfr_init2(c1, wp2);
    mpfr_t c2;
    mpfr_init2(c2, wp2);
    int id, nfs, kb, ndef;
    
    int idef;

   	int nkf, nf1, j3, nrs, nrf, nks, ip, kp, im, ip1;
    mpfr_t hrik;
    mpfr_init2(hrik, wp2);
    mpfr_t hiik;
    mpfr_init2(hiik, wp2);
    mpc_t hri;
    mpc_init3(hri, wp2, wp2);
    mpc_t zri;
    mpc_init3(zri, wp2, wp2);
    mpfr_t b;
    mpfr_init2(b, wp2);
    mpfr_t t1;
    mpfr_init2(t1, wp2);
    mpfr_t t2;
    mpfr_init2(t2, wp2);

	static int noback;


    // shift pointers
    --sm;
    --dele;
    --ndb;
    --ndefl;
    --ndel;
    --nxt;
    --supd;
    --evi;
    --evr;
	zi_dim1 = nm;
	zi_offset = 1 + zi_dim1;
	zi = zi - zi_offset;
	zr_dim1 = nm;
	zr_offset = 1 + zr_dim1;
	zr = zr - zr_offset;
	hi_dim1 = nm;
	hi_offset = 1 + hi_dim1;
	ai = ai - hi_offset;
	hr_dim1 = nm;
	hr_offset = 1 + hr_dim1;
	ar = ar - hr_offset;

	if (debug < 0 || debug >= 5) {
		printf("----- INPUT:\n");
		printf("matrix ar\n");
		for (int i=1; i<N+1; i++) {
			for (int j=1; j<N+1; j++) {
				mpfr_out_str(stdout, 10, 0, ar[i + hr_dim1*j], MPFR_RNDN);
				printf("\n");
			}
			printf("\n");
		}

		printf("matrix ai\n");
		for (int i=1; i<N+1; i++) {
			for (int j=1; j<N+1; j++) {
				mpfr_out_str(stdout, 10, 0, ai[i + hr_dim1*j], MPFR_RNDN);
				printf("\n");
			}
			printf("\n");
		}
		// getchar();
	}

    /***************/
    /*     STEP 1. */
    /***************/

    /*   ******************************************************************** */
    /*   *    TRANSFORM THE N BY N  COMPLEX MATRIX (HR,HI) INTO TRIANGULAR */
    /*   *    FORM BY USING THE EISPACK-ROUTINES (REF.  )                  * */
    /*   *    CBAL,COMMHES,COMLR2 AND CBABK2. */
    /*   *    NOTE0 */
    /*   *    A NEW PARAMETER (NOBACK) HAS BEEN ADDED TO COMLR2. */
    /*   *    NOBACK= 0  IMPLIES THAT THE EIGENVECTORS ARE NOT COMPUTED I.E */
    /*   *    THERE IS NO BACKSUBSTITUTION ,AND THE USER GETS THE TRANSFOR- */
    /*   *    MATIONS MADE DURING THE LR-PROCESS IN THE MATRICES ZR,ZI. */
    /*   *    THIS CHANGE CONSIST OF ONE EXTRA STATEMENT IN THE CODE OF */
    /*   *    COMLR20 */
    /*   *        IF(NOBACK.EQ. 0) GO TO 1001 */
    /*   * */
    /*   *    WHICH IS PLACED DIRECTLY  AFTER THE STATEMENT0 */
    /*   *    680   IF(N.EQ.1) GO TO 1001 */
    /*   * */
    /*   ******************************************************************** */


    /*   ******     BALANCE THE MATRIX                                    *** */
	cbal_(&nm, &N, &ar[hr_offset], &ai[hi_offset], &low, &igh, rad);
	
	if (debug < 0 || debug >= 10) {	
		printf("----- CBAL OUTPUT\n");
		printf("matrix ar\n");
		for (int i=1; i<N+1; i++) {
			for (int j=1; j<N+1; j++) {
				mpfr_out_str(stdout, 10, 0, ar[i + hr_dim1*j], MPFR_RNDN);
				printf("\n");
			}
			printf("\n");
		}

		printf("matrix ai\n");
		for (int i=1; i<N+1; i++) {
			for (int j=1; j<N+1; j++) {
				mpfr_out_str(stdout, 10, 0, ai[i + hr_dim1*j], MPFR_RNDN);
				printf("\n");
			}
			printf("\n");
		}
		// getchar();
	}

    /*   *****     COMPUTE EINF=EIN*FN,TOLD=TOL*FN WHERE */
    /*             FN=THE FROBENIUS NORM OF THE BALANCED MATRIX HR,HI  **** */
    /*                               MATRIX AR,AI)                    ***** */
    
	mpfr_set_d(fn, 0e0, MPFR_RNDN);
    for (int i__ = 1; i__ <= N; ++i__) {
        for (int j = 1; j <= N; ++j) {
			mpfr_mul(tmp0, ar[i__ + j * hr_dim1], ar[i__ + j * hr_dim1], MPFR_RNDN);
			mpfr_add(fn, fn, tmp0, MPFR_RNDN);
			mpfr_mul(tmp0, ai[i__ + j * hi_dim1], ai[i__ + j * hi_dim1], MPFR_RNDN);
			mpfr_add(fn, fn, tmp0, MPFR_RNDN);
        }
    }
	mpfr_sqrt(fn, fn, MPFR_RNDN);
	mpfr_mul(einf, fn, ein, MPFR_RNDN);
	mpfr_mul(told, fn, tol, MPFR_RNDN);

    /*   ******     TRANSFORM TO UPPER HESSENBERG FORM                    *** */
    comhes_(&nm, &N, &low, &igh, &ar[hr_offset], &ai[hi_offset], int__);

	if (debug < 0 || debug >= 20) {	
		printf("----- COMHES OUTPUT\n");
		printf("matrix ar\n");
		for (int i=1; i<N+1; i++) {
			for (int j=1; j<N+1; j++) {
				mpfr_out_str(stdout, 10, 0, ar[i + hr_dim1*j], MPFR_RNDN);
				printf("\n");
			}
			printf("\n");
		}
		// getchar();
	}

    /*   ******     TRANSFORM TO UPPER  TRIANGULAR FORM                   *** */
    comlr2_(&nm, &N, &low, &igh, int__,
        &ar[hr_offset], &ai[hi_offset],
        &evr[1], &evi[1], &zr[zr_offset], &zi[zi_offset],
        &ierr, &noback);
    
	if (debug < 0 || debug >= 30) {	
		printf("----- COMLR2 OUTPUT\n");
		printf("matrix ar\n");
		for (int i=1; i<N+1; i++) {
			for (int j=1; j<N+1; j++) {
				mpfr_out_str(stdout, 10, 0, ar[i + hr_dim1*j], MPFR_RNDN);
				printf("\n");
			}
			printf("\n");
		}
		printf("}\n");
		printf("vector evr\n");
		for (int i=1; i<N+1; i++) {
			mpfr_out_str(stdout, 10, 0, evr[i], MPFR_RNDN);
			printf("\n");
		}
		printf("}\n");
		// getchar();
	}

    /*   *****   TRANSFORM THE VECTORS  SO THAT ZR,ZI CONTAIN THE */
    /*           TRANSFORMATIONS WHICH TRANSFORMED HR,HI TO UPPER */
    /*           TRIANGULAR  FORM                                     ***** */
    cbabk2_(&nm, &N, &low, &igh, rad, &N, &zr[zr_offset], &zi[zi_offset]);

	if (debug < 0 || debug >= 40) {	
		printf("----- CBABK2 OUTPUT\n");
		printf("matrix zr\n");
		for (int i=1; i<N+1; i++) {
			for (int j=1; j<N+1; j++) {
				mpfr_out_str(stdout, 10, 0, zr[i + hr_dim1*j], MPFR_RNDN);
				printf("\n");
			}
			printf("\n");
		}
		// getchar();
	}

    /***************/
    /*     STEP 2. */
    /***************/

    /*   ******************************************************************** */
    /*   *    SORT THE DIAGONAL ELEMENTS OF THE TRIANGULAR MATRIX SO THAT */
    /*   *    THE NUMERICAL MULTIPLE EIGENVALUES APPEAR IN ADJACENT POSITIONS */
    /*   ******************************************************************** */

    for (int i__ = 2; i__ <= N; ++i__) {
		k = N - i__ + 2;

		/*   ******     K=N,N-1,...,2  (INDEX FOR ACTUAL EIGENVALUE)      ******* */

		/* Computing 2nd power */
		mpfr_sub(r__1, ar[k + k * hr_dim1], ar[k - 1 + (k - 1) * hr_dim1], MPFR_RNDN);
		/* Computing 2nd power */
		mpfr_sub(r__2, ai[k + k * hi_dim1], ai[k - 1 + (k - 1) * hi_dim1], MPFR_RNDN);
		mpfr_mul(r__1, r__1, r__1, MPFR_RNDN);
		mpfr_mul(r__2, r__2, r__2, MPFR_RNDN);
		mpfr_add(r__1, r__1, r__2, MPFR_RNDN);
		mpfr_sqrt(min__, r__1, MPFR_RNDN);
		kprim = k;
		k2 = k - 2;
		if(1 > k2) {
			goto L10020;
		}
		for (int j1 = 1; j1 <= k2; ++j1) {
			j = k2 - j1 + 1;
			/* Computing 2nd power */
			mpfr_sub(r__1, ar[j + j * hr_dim1], ar[k + k * hr_dim1], MPFR_RNDN);
			/* Computing 2nd power */
			mpfr_sub(r__2, ai[j + j * hi_dim1], ai[k + k * hi_dim1], MPFR_RNDN);
			mpfr_mul(r__1, r__1, r__1, MPFR_RNDN);
			mpfr_mul(r__2, r__2, r__2, MPFR_RNDN);
			mpfr_add(r__1, r__1, r__2, MPFR_RNDN);
			mpfr_sqrt(sk, r__1, MPFR_RNDN);
			if(mpfr_greaterequal_p(sk, min__)) {
				goto L2;
			}
			kprim = j;
			mpfr_set(min__, sk, MPFR_RNDN);
L2:;
// printf("L2:;\n");
;
		}
L10020:;
// printf("L10020:;\n");
		if(kprim == k) {
			goto L7;
		}

		if(kprim > k2) {
			goto L7;
		}
		for (int j = kprim; j <= k2; ++j) {
			/*   ******************************************************************** */
			/*   *    CHANGE PLACES OF THE DIAGONALELEMENTS (J,J) AND(J+1,J+1) WITH */
			/*   *    A UNITARY TRANSFORMATION */
			/*   ******************************************************************** */

			/* Computing 2nd power */
			mpfr_sub(r__1, ar[j + 1 + (j + 1) * hr_dim1], ar[j + j * hr_dim1], MPFR_RNDN);
			/* Computing 2nd power */
			mpfr_sub(r__2, ai[j + 1 + (j + 1) * hi_dim1], ai[j + j * hi_dim1], MPFR_RNDN);
			/* Computing 2nd power */
			mpfr_set(r__3, ar[j + (j + 1) * hr_dim1], MPFR_RNDN);
			/* Computing 2nd power */
			mpfr_set(r__4, ai[j + (j + 1) * hi_dim1], MPFR_RNDN);
			mpfr_mul(r__1, r__1, r__1, MPFR_RNDN);
			mpfr_mul(r__2, r__2, r__2, MPFR_RNDN);
			mpfr_mul(r__3, r__3, r__3, MPFR_RNDN);
			mpfr_mul(r__4, r__4, r__4, MPFR_RNDN);
			mpfr_add(r__1, r__1, r__2, MPFR_RNDN);
			mpfr_add(r__1, r__1, r__3, MPFR_RNDN);
			mpfr_add(r__1, r__1, r__4, MPFR_RNDN);
			mpfr_sqrt(ksi, r__1, MPFR_RNDN);
			mpfr_div(myr, ar[j + (j + 1) * hr_dim1], ksi, MPFR_RNDN);
			mpfr_div(myi, ai[j + (j + 1) * hi_dim1], ksi, MPFR_RNDN);
			mpfr_sub(nyr, ar[j + 1 + (j + 1) * hr_dim1], ar[j + j * hr_dim1], MPFR_RNDN);
			mpfr_div(nyr, nyr, ksi, MPFR_RNDN);
			mpfr_sub(nyi, ai[j + 1 + (j + 1) * hi_dim1], ai[j + j * hi_dim1], MPFR_RNDN);
			mpfr_div(nyi, nyi, ksi, MPFR_RNDN);
			mpfr_set(sl, ar[j + j * hr_dim1], MPFR_RNDN);
			mpfr_set(ar[j + j * hr_dim1], ar[j + 1 + (j + 1) * hr_dim1], MPFR_RNDN);
			mpfr_set(ar[j + 1 + (j + 1) * hr_dim1], sl, MPFR_RNDN);
			mpfr_set(sl, ai[j + j * hi_dim1], MPFR_RNDN);
			mpfr_set(ai[j + j * hi_dim1], ai[j + 1 + (j + 1) * hi_dim1], MPFR_RNDN);
			mpfr_set(ai[j + 1 + (j + 1) * hi_dim1], sl, MPFR_RNDN);

			mpfr_set(ai[j + (j + 1) * hi_dim1], ai[j + (j + 1) * hi_dim1], MPFR_RNDN);
			mpfr_neg(ai[j + (j + 1) * hi_dim1], ai[j + (j + 1) * hi_dim1], MPFR_RNDN);

			j1 = j - 1;
			if(1 > j1) {
				goto L10040;
			}
			for (int l = 1; l <= j1; ++l) {
				mpfr_set(hrlj, ar[l + j * hr_dim1], MPFR_RNDN);
				mpfr_set(hilj, ai[l + j * hi_dim1], MPFR_RNDN);
				mpfr_set(hrlj1, ar[l + (j + 1) * hr_dim1], MPFR_RNDN);
				mpfr_set(hilj1, ai[l + (j + 1) * hi_dim1], MPFR_RNDN);
				mpfr_mul(tmp0, myr, hrlj, MPFR_RNDN);
				mpfr_mul(tmp1, myi, hilj, MPFR_RNDN);
				mpfr_mul(tmp2, nyr, hrlj1, MPFR_RNDN);
				mpfr_mul(tmp3, nyi, hilj1, MPFR_RNDN);
				mpfr_sub(ar[l + j * hr_dim1], tmp0, tmp1, MPFR_RNDN);
				mpfr_add(ar[l + j * hr_dim1], ar[l + j * hr_dim1], tmp2, MPFR_RNDN);
				mpfr_sub(ar[l + j * hr_dim1], ar[l + j * hr_dim1], tmp3, MPFR_RNDN);
				mpfr_mul(tmp0, myr, hilj, MPFR_RNDN);
				mpfr_mul(tmp1, myi, hrlj, MPFR_RNDN);
				mpfr_mul(tmp2, nyr, hilj1, MPFR_RNDN);
				mpfr_mul(tmp3, nyi, hrlj1, MPFR_RNDN);
				mpfr_add(ai[l + j * hi_dim1], tmp0, tmp1, MPFR_RNDN);
				mpfr_add(ai[l + j * hi_dim1], ai[l + j * hi_dim1], tmp2, MPFR_RNDN);
				mpfr_add(ai[l + j * hi_dim1], ai[l + j * hi_dim1], tmp3, MPFR_RNDN);
				mpfr_mul(tmp0, nyr, hrlj, MPFR_RNDN);
				mpfr_mul(tmp1, nyi, hilj, MPFR_RNDN);
				mpfr_mul(tmp2, myr, hrlj1, MPFR_RNDN);
				mpfr_mul(tmp3, myi, hilj1, MPFR_RNDN);
				mpfr_set(ar[l + (j + 1) * hr_dim1], tmp0, MPFR_RNDN);
				mpfr_neg(ar[l + (j + 1) * hr_dim1], ar[l + (j + 1) * hr_dim1], MPFR_RNDN);
				mpfr_sub(ar[l + (j + 1) * hr_dim1], ar[l + (j + 1) * hr_dim1], tmp1, MPFR_RNDN);
				mpfr_add(ar[l + (j + 1) * hr_dim1], ar[l + (j + 1) * hr_dim1], tmp2, MPFR_RNDN);
				mpfr_add(ar[l + (j + 1) * hr_dim1], ar[l + (j + 1) * hr_dim1], tmp3, MPFR_RNDN);
				mpfr_mul(tmp0, nyr, hilj, MPFR_RNDN);
				mpfr_mul(tmp1, nyi, hrlj, MPFR_RNDN);
				mpfr_mul(tmp2, myr, hilj1, MPFR_RNDN);
				mpfr_mul(tmp3, myi, hrlj1, MPFR_RNDN);
				mpfr_set(ai[l + (j + 1) * hi_dim1], tmp0, MPFR_RNDN);
				mpfr_neg(ai[l + (j + 1) * hi_dim1], ai[l + (j + 1) * hi_dim1], MPFR_RNDN);
				mpfr_add(ai[l + (j + 1) * hi_dim1], ai[l + (j + 1) * hi_dim1], tmp1, MPFR_RNDN);
				mpfr_add(ai[l + (j + 1) * hi_dim1], ai[l + (j + 1) * hi_dim1], tmp2, MPFR_RNDN);
				mpfr_sub(ai[l + (j + 1) * hi_dim1], ai[l + (j + 1) * hi_dim1], tmp3, MPFR_RNDN);
			/* L3: */
			}
L10040:;
// printf("L10040:;\n");

			/*   ******     ACCUMULATE TRANSFORMATIONS                        ******* */
			for (int l = 1; l <= N; ++l) {
				mpfr_set(hrlj, zr[l + j * zr_dim1], MPFR_RNDN);
				mpfr_set(hilj, zi[l + j * zi_dim1], MPFR_RNDN);
				mpfr_set(hrlj1, zr[l + (j + 1) * zr_dim1], MPFR_RNDN);
				mpfr_set(hilj1, zi[l + (j + 1) * zi_dim1], MPFR_RNDN);
				mpfr_mul(tmp0, myr, hrlj, MPFR_RNDN);
				mpfr_mul(tmp1, myi, hilj, MPFR_RNDN);
				mpfr_mul(tmp2, nyr, hrlj1, MPFR_RNDN);
				mpfr_mul(tmp3, nyi, hilj1, MPFR_RNDN);
				mpfr_sub(zr[l + j * zr_dim1], tmp0, tmp1, MPFR_RNDN);
				mpfr_add(zr[l + j * zr_dim1], zr[l + j * zr_dim1], tmp2, MPFR_RNDN);
				mpfr_sub(zr[l + j * zr_dim1], zr[l + j * zr_dim1], tmp3, MPFR_RNDN);
				mpfr_mul(tmp0, myr, hilj, MPFR_RNDN);
				mpfr_mul(tmp1, myi, hrlj, MPFR_RNDN);
				mpfr_mul(tmp2, nyr, hilj1, MPFR_RNDN);
				mpfr_mul(tmp3, nyi, hrlj1, MPFR_RNDN);
				mpfr_add(zi[l + j * zi_dim1], tmp0, tmp1, MPFR_RNDN);
				mpfr_add(zi[l + j * zi_dim1], zi[l + j * zi_dim1], tmp2, MPFR_RNDN);
				mpfr_add(zi[l + j * zi_dim1], zi[l + j * zi_dim1], tmp3, MPFR_RNDN);
				mpfr_mul(tmp0, nyr, hrlj, MPFR_RNDN);
				mpfr_mul(tmp1, nyi, hilj, MPFR_RNDN);
				mpfr_mul(tmp2, myr, hrlj1, MPFR_RNDN);
				mpfr_mul(tmp3, myi, hilj1, MPFR_RNDN);
				mpfr_set(zr[l + (j + 1) * zr_dim1], tmp0, MPFR_RNDN);
				mpfr_neg(zr[l + (j + 1) * zr_dim1], zr[l + (j + 1) * zr_dim1], MPFR_RNDN);
				mpfr_sub(zr[l + (j + 1) * zr_dim1], zr[l + (j + 1) * zr_dim1], tmp1, MPFR_RNDN);
				mpfr_add(zr[l + (j + 1) * zr_dim1], zr[l + (j + 1) * zr_dim1], tmp2, MPFR_RNDN);
				mpfr_add(zr[l + (j + 1) * zr_dim1], zr[l + (j + 1) * zr_dim1], tmp3, MPFR_RNDN);
				mpfr_mul(tmp0, nyr, hilj, MPFR_RNDN);
				mpfr_mul(tmp1, nyi, hrlj, MPFR_RNDN);
				mpfr_mul(tmp2, myr, hilj1, MPFR_RNDN);
				mpfr_mul(tmp3, myi, hrlj1, MPFR_RNDN);
				mpfr_set(zi[l + (j + 1) * zi_dim1], tmp0, MPFR_RNDN);
				mpfr_neg(zi[l + (j + 1) * zi_dim1], zi[l + (j + 1) * zi_dim1], MPFR_RNDN);
				mpfr_add(zi[l + (j + 1) * zi_dim1], zi[l + (j + 1) * zi_dim1], tmp1, MPFR_RNDN);
				mpfr_add(zi[l + (j + 1) * zi_dim1], zi[l + (j + 1) * zi_dim1], tmp2, MPFR_RNDN);
				mpfr_sub(zi[l + (j + 1) * zi_dim1], zi[l + (j + 1) * zi_dim1], tmp3, MPFR_RNDN);
			/* L4: */
			}

			j2 = j + 2;
	    	if(j2 > N) {
				goto L6;
	    	}
			for (int l = j2; l <= N; ++l) {
				mpfr_set(hrlj, ar[j + l * hr_dim1], MPFR_RNDN);
				mpfr_set(hilj, ai[j + l * hi_dim1], MPFR_RNDN);
				mpfr_set(hrlj1, ar[j + 1 + l * hr_dim1], MPFR_RNDN);
				mpfr_set(hilj1, ai[j + 1 + l * hi_dim1], MPFR_RNDN);
				mpfr_mul(tmp0, myr, hrlj, MPFR_RNDN);
				mpfr_mul(tmp1, myi, hilj, MPFR_RNDN);
				mpfr_mul(tmp2, nyr, hrlj1, MPFR_RNDN);
				mpfr_mul(tmp3, nyi, hilj1, MPFR_RNDN);
				mpfr_add(ar[j + l * hr_dim1], tmp0, tmp1, MPFR_RNDN);
				mpfr_add(ar[j + l * hr_dim1], ar[j + l * hr_dim1], tmp2, MPFR_RNDN);
				mpfr_add(ar[j + l * hr_dim1], ar[j + l * hr_dim1], tmp3, MPFR_RNDN);
				mpfr_mul(tmp0, myr, hilj, MPFR_RNDN);
				mpfr_mul(tmp1, myi, hrlj, MPFR_RNDN);
				mpfr_mul(tmp2, nyr, hilj1, MPFR_RNDN);
				mpfr_mul(tmp3, nyi, hrlj1, MPFR_RNDN);
				mpfr_sub(ai[j + l * hi_dim1], tmp0, tmp1, MPFR_RNDN);
				mpfr_add(ai[j + l * hi_dim1], ai[j + l * hi_dim1], tmp2, MPFR_RNDN);
				mpfr_sub(ai[j + l * hi_dim1], ai[j + l * hi_dim1], tmp3, MPFR_RNDN);
				mpfr_mul(tmp0, nyr, hrlj, MPFR_RNDN);
				mpfr_mul(tmp1, nyi, hilj, MPFR_RNDN);
				mpfr_mul(tmp2, myr, hrlj1, MPFR_RNDN);
				mpfr_mul(tmp3, myi, hilj1, MPFR_RNDN);
				mpfr_set(ar[j + 1 + l * hr_dim1], tmp0, MPFR_RNDN);
				mpfr_neg(ar[j + 1 + l * hr_dim1], ar[j + 1 + l * hr_dim1], MPFR_RNDN);
				mpfr_add(ar[j + 1 + l * hr_dim1], ar[j + 1 + l * hr_dim1], tmp1, MPFR_RNDN);
				mpfr_add(ar[j + 1 + l * hr_dim1], ar[j + 1 + l * hr_dim1], tmp2, MPFR_RNDN);
				mpfr_sub(ar[j + 1 + l * hr_dim1], ar[j + 1 + l * hr_dim1], tmp3, MPFR_RNDN);
				mpfr_mul(tmp0, nyr, hilj, MPFR_RNDN);
				mpfr_mul(tmp1, nyi, hrlj, MPFR_RNDN);
				mpfr_mul(tmp2, myr, hilj1, MPFR_RNDN);
				mpfr_mul(tmp3, myi, hrlj1, MPFR_RNDN);
				mpfr_set(ai[j + 1 + l * hi_dim1], tmp0, MPFR_RNDN);
				mpfr_neg(ai[j + 1 + l * hi_dim1], ai[j + 1 + l * hi_dim1], MPFR_RNDN);
				mpfr_sub(ai[j + 1 + l * hi_dim1], ai[j + 1 + l * hi_dim1], tmp1, MPFR_RNDN);
				mpfr_add(ai[j + 1 + l * hi_dim1], ai[j + 1 + l * hi_dim1], tmp2, MPFR_RNDN);
				mpfr_add(ai[j + 1 + l * hi_dim1], ai[j + 1 + l * hi_dim1], tmp3, MPFR_RNDN);
/* L5: */
			}
L6:;
// printf("L6:;\n");
;
		}

		/*   ******     CONTINUE THE SORTING PROCEDURE                    ******* */
L7:;
// printf("L7:;\n");
;
    }
	i__2 = nm;
    for (int i__ = 1; i__ <= i__2; ++i__) {
		ndel[i__] = 0;
		ndefl[i__] = 0;
		nxt[i__] = 0;
		ndb[i__] = 0;
		int__[i__ - 1] = 0;
		mpfr_set_d(supd[i__], 0e0, MPFR_RNDN);
/* L8: */
    }

	ns = N;
	ibl = 1;

L9:;
// printf("L9:;\n");

	ne = ns;
	nd = ne;
	ndel[ibl] = ne;
	if(ns < 1) {
		goto L21;
	}
L10:;
// printf("L10:;\n");
	mpfr_set_d(mpr, 0e0, MPFR_RNDN);
	mpfr_set_d(mpi, 0e0, MPFR_RNDN);
    if(nd > ne) {
		goto L10080;
    }
    for (int k = nd; k <= ne; ++k) {
		mpfr_add(mpr, mpr, ar[k + k * hr_dim1], MPFR_RNDN);
		mpfr_add(mpi, mpi, ai[k + k * hi_dim1], MPFR_RNDN);
/* L11: */
    }
L10080:;
// printf("L10080:;\n");
	mb = ne - nd + 1;
	mpfr_div_si(mpr, mpr, mb, MPFR_RNDN);
	mpfr_div_si(mpi, mpi, mb, MPFR_RNDN);
    if(nd <= nf) {
		goto L20;
    }
	mpfr_set_d(rgmax, 0e0, MPFR_RNDN);
    if(nd > ne) {
		goto L10090;
    }
    for (int k = nd; k <= ne; ++k) {
		/* Computing 2nd power */
		mpfr_sub(r__1, ar[k + k * hr_dim1], mpr, MPFR_RNDN);
		/* Computing 2nd power */
		mpfr_sub(r__2, ai[k + k * hi_dim1], mpi, MPFR_RNDN);
		mpfr_mul(r__1, r__1, r__1, MPFR_RNDN);
		mpfr_mul(r__2, r__2, r__2, MPFR_RNDN);
		mpfr_add(r__1, r__1, r__2, MPFR_RNDN);
		mpfr_sqrt(cabs1, r__1, MPFR_RNDN);
		/* L12: */
		mpfr_max(rgmax, rgmax,cabs1, MPFR_RNDN);
    }
L10090:;
// printf("L10090:;\n");
	nd2 = nd - 1;
	/* Computing 2nd power */
	mpfr_sub(r__1, ar[nd2 + nd2 * hr_dim1], mpr, MPFR_RNDN);
	/* Computing 2nd power */
	mpfr_sub(r__2, ai[nd2 + nd2 * hi_dim1], mpi, MPFR_RNDN);
	mpfr_mul(r__1, r__1, r__1, MPFR_RNDN);
	mpfr_mul(r__2, r__2, r__2, MPFR_RNDN);
	mpfr_add(r__1, r__1, r__2, MPFR_RNDN);
	mpfr_sqrt(rigmin, r__1, MPFR_RNDN);
	nd2 = nd - 2;
    if(nf > nd2) {
		goto L10100;
    }
    for (int k = nf; k <= nd2; ++k) {
		/* Computing 2nd power */
		mpfr_sub(r__1, ar[k + k * hr_dim1], mpr, MPFR_RNDN);
		/* Computing 2nd power */
		mpfr_sub(r__2, ai[k + k * hi_dim1], mpi, MPFR_RNDN);
		mpfr_mul(r__1, r__1, r__1, MPFR_RNDN);
		mpfr_mul(r__2, r__2, r__2, MPFR_RNDN);
		mpfr_add(r__1, r__1, r__2, MPFR_RNDN);
		mpfr_sqrt(cabs1, r__1, MPFR_RNDN);
/* L13: */
		mpfr_min(rigmin, rigmin,cabs1, MPFR_RNDN);
    }
L10100:;
// printf("L10100:;\n");
	nd2 = ne + 1;
    if(nd2 > N) {
		goto L10110;
    }
    for (int k = nd2; k <= N; ++k) {
		/* Computing 2nd power */
		mpfr_sub(r__1, ar[k + k * hr_dim1], mpr, MPFR_RNDN);
		/* Computing 2nd power */
		mpfr_sub(r__2, ai[k + k * hi_dim1], mpi, MPFR_RNDN);
		mpfr_mul(r__1, r__1, r__1, MPFR_RNDN);
		mpfr_mul(r__2, r__2, r__2, MPFR_RNDN);
		mpfr_add(r__1, r__1, r__2, MPFR_RNDN);
		mpfr_sqrt(cabs1, r__1, MPFR_RNDN);
/* L14: */
		mpfr_min(rigmin, rigmin,cabs1, MPFR_RNDN);
    }
L10110:;
// printf("L10110:;\n");
	mpfr_add(radi, rgmax, rigmin, MPFR_RNDN);
	mpfr_div_d(radi, radi, 2e0, MPFR_RNDN);
    for (int k = 1; k <= N; ++k) {
		/* Computing 2nd power */
		mpfr_sub(r__1, ar[k + k * hr_dim1], mpr, MPFR_RNDN);
		/* Computing 2nd power */
		mpfr_sub(r__2, ai[k + k * hi_dim1], mpi, MPFR_RNDN);
		mpfr_mul(r__1, r__1, r__1, MPFR_RNDN);
		mpfr_mul(r__2, r__2, r__2, MPFR_RNDN);
		mpfr_add(r__1, r__1, r__2, MPFR_RNDN);
		mpfr_sqrt(cabs1, r__1, MPFR_RNDN);
/* L15: */
		mpfr_sub(cabs1, cabs1, radi, MPFR_RNDN);
		mpfr_abs(rad[k - 1], cabs1, MPFR_RNDN);
    }

	/*   ******************************************************************** */
	/*   *    ABOVE WE HAVE  FIXED THE CIRCLE WITH THE RADIUS RGMAX AND */
	/*   *    DRAWN CIRCLES  AROUND EVERY EIGENVALUE WHICH TOUCH AT THE */
	/*   *    HALF DISTRANCE (=(RGMAX + RIGMIN)/2) */
	/*   *    THESE RADIUS ARE STORED IN THE VECTOR RAD(K) */
	/*   *    NOW COMPUTE  DT(K)= THE DIAGONAL SIMILARITY TRANSFORMATION */
	/*   *    AND INVESTIGATE IF WE MANAGED TO ISOLATE THE LAST M EIGEN- */
	/*   *    VALUES */
	/*   ******************************************************************** */

	k = N;
	L16:;
// printf("	L16:;\n");
	mpfr_set_d(rsum, 0e0, MPFR_RNDN);
    if(k == N) {
		goto L18;
    }
	kp1 = k + 1;
	/*     RSUM=   SUM ( T(K,I)/D(I))  I=K+1,N */
    if(kp1 > N) {
		goto L10130;
    }
    for (int i__ = kp1; i__ <= N; ++i__) {
/* L17: */
		/* Computing 2nd power */
		mpfr_set(r__1, ar[k + i__ * hr_dim1], MPFR_RNDN);
		/* Computing 2nd power */
		mpfr_set(r__2, ai[k + i__ * hi_dim1], MPFR_RNDN);
		mpfr_mul(r__1, r__1, r__1, MPFR_RNDN);
		mpfr_mul(r__2, r__2, r__2, MPFR_RNDN);
		mpfr_add(r__1, r__1, r__2, MPFR_RNDN);
		mpfr_sqrt(r__1, r__1, MPFR_RNDN);
		mpfr_div(r__1, r__1, dt[i__ - 1], MPFR_RNDN);
		mpfr_add(rsum, rsum, r__1, MPFR_RNDN);
    }
L10130:;
// printf("L10130:;\n");
L18:;
// printf("L18:;\n");
	mpfr_add(dt[k - 1], einf, rsum, MPFR_RNDN);
	mpfr_div(dt[k - 1], rad[k - 1], dt[k - 1], MPFR_RNDN);
    if(mpfr_less_p(dt[k - 1], mpfr_one)) {
		goto L19;
    }
    --k;
    if(k >= nf) {
		goto L16;
    }
    goto L20;
L19:;
// printf("L19:;\n");
    --nd;
    goto L10;
L20:;
// printf("L20:;\n");

	/*   ******     STORE THE NUMERICAL MULTIPLE EIGENVALUE IN CSHTR(IBL) */
	/*              , CSHTI(IBL)                                      ******* */
	mpfr_set(cshtr[ibl - 1], mpr, MPFR_RNDN);
	mpfr_set(cshti[ibl - 1], mpi, MPFR_RNDN);
    ++ibl;
	ns = nd - 1;
    goto L9;
L21:;
// printf("L21:;\n");
	nblock = ibl - 1;

	ibl2 = ibl / 2;
	j = nblock + 1;
    if(1 > ibl2) {
		goto L10140;
    }
    for (int i__ = 1; i__ <= ibl2; ++i__) {
		k = ndel[j];
		ndel[j] = ndel[i__];
		ndel[i__] = k;
		if(i__ == j - 1) {
			goto L22;
		}
		mpfr_set(shtr, cshtr[j - 2], MPFR_RNDN);
		mpfr_set(shti, cshti[j - 2], MPFR_RNDN);
		mpfr_set(cshtr[j - 2], cshtr[i__ - 1], MPFR_RNDN);
		mpfr_set(cshti[j - 2], cshti[i__ - 1], MPFR_RNDN);
		mpfr_set(cshtr[i__ - 1], shtr, MPFR_RNDN);
		mpfr_set(cshti[i__ - 1], shti, MPFR_RNDN);
L22:;
// printf("L22:;\n");
		--j;
    }
L10140:;
// printf("L10140:;\n");
ndel[1] = 0;

    /*   ******************************************************************** */
    /*   *    CALL THE USER WRITTEN ROUTINE DECIDE ,WHICH POSSIBLY */
    /*   *    CHANGES THE GROUPING  AND/OR THE VALUES OF THE EIGENVALUES */
    /*   ******************************************************************** */
    decide_(&nm, &N, &ndel[1], cshtr, cshti, &nblock, &ar[hr_offset], &ai[
	    hi_offset]);

	
    // printf("vector ndel\n");
	// printf("{");
    // for (int i=1; i<N+1; i++) {
    //     printf("%d, ", ndel[i]);
    // }
    // printf("}\n");
	// getchar();

    /***************/
    /*     STEP 4. */
    /***************/

    /*   ******************************************************************** */
    /*   *    ELIMINATE THE ELEMENTS ABOVE THE MAIN DIAGONAL AND OUTSIDE */
    /*   *    THE BLOCKS CORRESPONDING TO NUMERICAL MULTIPLE EIGENVALUES */
    /*   *    THE PROCESS IS CARRIED OUT IN ROWS WORKING UPWARDS,AND IN */
    /*   *    COLUMNS FROM LEFT TO RIHGT. */
    /*   *    COMPUTE  THE BLOCKS M(K,J) AND M(K,J)-INVERSE(SEE THE ALGO- */
    /*   *    RITHM) IN ORDER TO ESTIMATE THE SPECTRAL PROJECTORS SM(K) */
    /*   *    FOR EIGENVALUE  K, K=1,2,...,NBLOCK */
    /*   *    THE M(K,J)-INVERSE BLOCKS ARE STORED IN THE CORRESPONDING */
    /*   *    ELIMINATED BLOCKS OF HR,HI. */
    /*   ******************************************************************** */

    for (int i__ = 1; i__ <= nblock; ++i__) {
		mpfr_set_d(sm[i__], 0e0, MPFR_RNDN);
/* L23: */
		mpfr_set_d(dt[i__ - 1], 0e0, MPFR_RNDN);
    }
	l1 = nblock + 1;
    for (j2 = 1; j2 <= nblock; ++j2) {
		ibl = nblock - j2 + 1;
		ns = ndel[ibl];
		if(ns < 1) {
			goto L33;
		}
		nf = ndel[ibl - 1] + 1;

		if(nf > ns) {
			goto L10170;
		}
		for (int i2 = nf; i2 <= ns; ++i2) {
			k5 = 0;
			i__ = ns - i2 + nf;
			k2 = ns + 1;

			if(k2 > N) {
				goto L10190;
			}
			for (k = k2; k <= N; ++k) {
				mpfr_set_d(radr[k - 1], 0e0, MPFR_RNDN);
/* L24: */
				mpfr_set_d(rad[k - 1], 0e0, MPFR_RNDN);
	    	}
			for (k = k2; k <= N; ++k) {
				mpc_set_fr_fr(m, ar[i__ + k * hr_dim1], ai[i__ + k * hi_dim1], MPFR_RNDN);
				mpfr_sub(r__1, ar[k + k * hr_dim1], ar[i__ + i__ * hr_dim1], MPFR_RNDN);
				mpfr_sub(r__2, ai[k + k * hi_dim1], ai[i__ + i__ * hi_dim1], MPFR_RNDN);
				mpc_set_fr_fr(q__1, r__1, r__2, MPFR_RNDN);
				mpc_div(m, m, q__1, MPFR_RNDN);
				i1 = i__ - 1;
				if(1 > i1) {
					goto L10200;
				}
				for (j = 1; j <= i1; ++j) {
					mpc_set_fr_fr(x, ar[j + i__ * hr_dim1], ai[j + i__ * hi_dim1], MPFR_RNDN);
					mpc_mul(x, x, m, MPFR_RNDN);
					mpc_set_fr_fr(q__1, ar[j + k * hr_dim1], ai[j + k * hi_dim1], MPFR_RNDN);
					mpc_add(x, x, q__1, MPFR_RNDN);
					mpfr_set(ar[j + k * hr_dim1], mpc_realref(x), MPFR_RNDN);
/* L25: */
					mpfr_set(ai[j + k * hi_dim1], mpc_imagref(x), MPFR_RNDN);
				}
L10200:;
// printf("L10200:;\n");
				k1 = k + 1;
				if(k <= ndel[l1]) {
					goto L26;
				}
				++l1;
				++k5;
L26:;
// printf("L26:;\n");
				nss = ndel[l1];
				if(k1 > nss) {
					goto L10210;
				}
				for (j = k1; j <= nss; ++j) {
					mpc_set_fr_fr(x, ar[k + j * hr_dim1], ai[k + j * hi_dim1], MPFR_RNDN);
					mpc_mul(x, x, m, MPFR_RNDN);
					mpc_set_fr_fr(q__1, ar[i__ + j * hr_dim1], ai[i__ + j * hi_dim1], MPFR_RNDN);
					mpc_sub(x, q__1, x, MPFR_RNDN);
					mpfr_set(ar[i__ + j * hr_dim1], mpc_realref(x), MPFR_RNDN);
/* L27: */
					mpfr_set(ai[i__ + j * hi_dim1], mpc_imagref(x), MPFR_RNDN);
				}
L10210:;
// printf("L10210:;\n");

				/*   ******     ACCUMULATE  TRANSFORMATIONS                       ******* */
				for (j = 1; j <= N; ++j) {
					mpc_set_fr_fr(x, zr[j + i__ * zr_dim1], zi[j + i__ * zi_dim1], MPFR_RNDN);
					mpc_mul(x, x, m, MPFR_RNDN);
					mpc_set_fr_fr(q__1, zr[j + k * zr_dim1], zi[j + k * zi_dim1], MPFR_RNDN);
					mpc_add(x, x, q__1, MPFR_RNDN);
					// printf("x:\n");
					// mpc_out_str(stdout, 10, 0, x, MPFR_RNDN);
					// printf("\n");
					// mpfr_out_str(stdout, 10, 0, mpc_realref(x), MPFR_RNDN);
					// printf("\n");
					// mpfr_out_str(stdout, 10, 0, mpc_imagref(x), MPFR_RNDN);
					// printf("\n");
					mpfr_set(zr[j + k * zr_dim1], mpc_realref(x), MPFR_RNDN);
/* L28: */
					mpfr_set(zi[j + k * zi_dim1], mpc_imagref(x), MPFR_RNDN);
				}
				mpfr_sub(ar[i__ + k * hr_dim1], radr[k - 1], mpc_realref(m), MPFR_RNDN);
				mpfr_sub(ai[i__ + k * hi_dim1], rad[k - 1], mpc_imagref(m), MPFR_RNDN);
				k3 = nss + 1;

				if(k3 > N) {
					goto L10230;
				}
				for (int l = k3; l <= N; ++l) {
				mpc_set_fr_fr(x, ar[k + l * hr_dim1], ai[k + l * hi_dim1], MPFR_RNDN);
				mpc_mul(x, x, m, MPFR_RNDN);
				mpc_neg(x, x, MPFR_RNDN);
				mpfr_add(radr[l - 1], radr[l - 1], mpc_realref(x), MPFR_RNDN);
				mpfr_add(rad[l - 1], rad[l - 1], mpc_imagref(x), MPFR_RNDN);
/* L29: */
				}
L10230:;
// printf("L10230:;\n");

				mpc_norm(tmp0, m, MPFR_RNDN);
				mpfr_add(sm[l1 - 1], sm[l1 - 1], tmp0, MPFR_RNDN);
				mpfr_mul(tmp0, ar[i__ + k * hr_dim1], ar[i__ + k * hr_dim1], MPFR_RNDN);
				mpfr_mul(tmp1, ai[i__ + k * hi_dim1], ai[i__ + k * hi_dim1], MPFR_RNDN);
				mpfr_add(dt[ibl - 2], tmp0, tmp1, MPFR_RNDN);
/* L30: */
	    	}
L10190:;
// printf("L10190:;\n");

			l1 = l1 - k5;
/* L31: */
		}
L10170:;
// printf("L10170:;\n");
		--l1;
/* L32: */
    }
L33:;
// printf("L33:;\n");
	for (int i__ = 1; i__ <= nblock; ++i__) {
/* L34: */
		mpfr_sqrt(sm[i__], sm[i__], MPFR_RNDN);
		mpfr_add_d(sm[i__], sm[i__], 1e0, MPFR_RNDN);
		mpfr_sqrt(cmp0, dt[i__ - 1], MPFR_RNDN);
		mpfr_add_d(cmp0, cmp0, 1e0, MPFR_RNDN);
		mpfr_mul(sm[i__], sm[i__], cmp0, MPFR_RNDN);
	}


	/*   *****   ZERO THE ELEMENTS BEYOND THE MAIN DIAGONAL. */
	/*           THIS IS MADE BECAUSE THE ALGORITHMS WHICH COMPUTE */
	/*           THE NILPOTENT MATRIX OF A DIAGONAL BLOCK, */
	/*           WORK ON THE WHOLE BLOCK,WHICH IS ASSUMED TO */
	/*           BE UPPER TRIANGULAR.                                 ***** */

	for (int i__ = 1; i__ <= N; ++i__) {
		for (j = 1; j <= N; ++j) {
			if(i__ <= j) {
				goto L38;
			}
			mpfr_set_d(ar[i__ + j * hr_dim1], -0e0, MPFR_RNDN);
			mpfr_set_d(ai[i__ + j * hi_dim1], -0e0, MPFR_RNDN);
L38:;
// printf("L38:;\n");
;
		}
	}

	if (debug < 0 || debug >= 50) {	
		printf("----- STEP 4 OUTPUT\n");
		printf("matrix ar\n");
		for (int i=1; i<N+1; i++) {
			for (int j=1; j<N+1; j++) {
				mpfr_out_str(stdout, 10, 0, ar[i + hr_dim1*j], MPFR_RNDN);
				printf("\n");
			}
			printf("\n");
		}

		printf("matrix ai\n");
		for (int i=1; i<N+1; i++) {
			for (int j=1; j<N+1; j++) {
				mpfr_out_str(stdout, 10, 0, ai[i + hr_dim1*j], MPFR_RNDN);
				printf("\n");
			}
			printf("\n");
		}

		// printf("matrix zr\n");
		for (int i=1; i<N+1; i++) {
			for (int j=1; j<N+1; j++) {
				// mpfr_out_str(stdout, 10, 0, zr[i + hr_dim1*j], MPFR_RNDN);
				// printf("\n");
			}
			// printf("\n");
		}

		// printf("matrix zi\n");
		for (int i=1; i<N+1; i++) {
			for (int j=1; j<N+1; j++) {
				// mpfr_out_str(stdout, 10, 0, zi[i + hr_dim1*j], MPFR_RNDN);
				// printf("\n");
			}
			// printf("\n");
		}

    // printf("vector evr\n");
    for (int i=1; i<N+1; i++) {
		// mpfr_out_str(stdout, 10, 0, evr[i], MPFR_RNDN);
		// printf("\n");
    }
	// printf("}\n");
    
	// printf("vector evi\n");
    for (int i=1; i<N+1; i++) {
		// mpfr_out_str(stdout, 10, 0, evi[i], MPFR_RNDN);
		// printf("\n");
    }
	// printf("}\n");

		//getchar();
	}

    /***************/
    /*     STEP 5. */
    /***************/

    /*   ****************************************************************** */
    /*   *    ORTHONORMALIZE THE VECTORS CORRESPONDING TO A NUMERICAL */
    /*   *    MULTIPLE EIGENVALUE (1,..,NBLOCK) BY THE MODIFIED GRAMM */
    /*   *    SCHMIDT  PROCESS. */
    /*   *    THE VECTORS ARE UPDATED AND STORED IN ZR,ZI */
    /*   ****************************************************************** */


	ns = 0;
    for (ibl = 1; ibl <= nblock; ++ibl) {
		nf = ns + 1;
		// printf("nblock = %d\n", nblock);
		// printf("ibl = %d\n", ibl);
		// printf("ndel = %d\n", ndel[ibl+1]);
		ns = ndel[ibl + 1];
		if(nf > ns) {
			goto L49;
		}
		for (k = nf; k <= ns; ++k) {
			mpfr_set_d(rsum, 0e0, MPFR_RNDN);
	    	for (i__ = 1; i__ <= N; ++i__) {
/* L39: */
				mpfr_mul(tmp0, zr[i__ + k * zr_dim1], zr[i__ + k * zr_dim1], MPFR_RNDN);
				mpfr_mul(tmp1, zi[i__ + k * zi_dim1], zi[i__ + k * zi_dim1], MPFR_RNDN);
				mpfr_add(rsum, rsum, tmp0, MPFR_RNDN);
				mpfr_add(rsum, rsum, tmp1, MPFR_RNDN);
	    	}
			mpfr_sqrt(rsum, rsum, MPFR_RNDN);
			mpfr_set(evr[1], rsum, MPFR_RNDN);
			mpfr_set_d(evi[1], 0e0, MPFR_RNDN);
	    	for (i__ = 1; i__ <= N; ++i__) {
				mpfr_div(zr[i__ + k * zr_dim1], zr[i__ + k * zr_dim1], rsum, MPFR_RNDN);
/* L40: */
				mpfr_div(zi[i__ + k * zi_dim1], zi[i__ + k * zi_dim1], rsum, MPFR_RNDN);
	    	}
			k1 = k + 1;
			if(k1 > ns) {
				goto L10300;
			}
	    	for (j = k1; j <= ns; ++j) {
				mpfr_set_d(xr, 0e0, MPFR_RNDN);
				mpfr_set_d(xi, 0e0, MPFR_RNDN);
				for (i__ = 1; i__ <= N; ++i__) {
					mpfr_mul(tmp0, zr[i__ + k * zr_dim1], zr[i__ + j * zr_dim1], MPFR_RNDN);
					mpfr_mul(tmp1, zi[i__ + k * zi_dim1], zi[i__ + j * zi_dim1], MPFR_RNDN);
					mpfr_add(xr, xr, tmp0, MPFR_RNDN);
					mpfr_add(xr, xr, tmp1, MPFR_RNDN);
/* L41: */
					mpfr_mul(tmp0, zr[i__ + k * zr_dim1], zi[i__ + j * zi_dim1], MPFR_RNDN);
					mpfr_mul(tmp1, zi[i__ + k * zi_dim1], zr[i__ + j * zr_dim1], MPFR_RNDN);
					mpfr_add(xi, xi, tmp0, MPFR_RNDN);
					mpfr_sub(xi, xi, tmp1, MPFR_RNDN);
				}
				l = j - k + 1;
				mpfr_set(evr[l], xr, MPFR_RNDN);
				mpfr_set(evi[l], xi, MPFR_RNDN);
				for (i__ = 1; i__ <= N; ++i__) {
					mpfr_mul(tmp0, evr[l], zr[i__ + k * zr_dim1], MPFR_RNDN);
					mpfr_mul(tmp1, evi[l], zi[i__ + k * zi_dim1], MPFR_RNDN);
					mpfr_sub(zr[i__ + j * zr_dim1], zr[i__ + j * zr_dim1], tmp0, MPFR_RNDN);
					mpfr_add(zr[i__ + j * zr_dim1], zr[i__ + j * zr_dim1], tmp1, MPFR_RNDN);
					mpfr_mul(tmp0, evr[l], zi[i__ + k * zi_dim1], MPFR_RNDN);
					mpfr_mul(tmp1, evi[l], zr[i__ + k * zr_dim1], MPFR_RNDN);
					mpfr_sub(zi[i__ + j * zi_dim1], zi[i__ + j * zi_dim1], tmp0, MPFR_RNDN);
					mpfr_sub(zi[i__ + j * zi_dim1], zi[i__ + j * zi_dim1], tmp1, MPFR_RNDN);
/* L42: */
				}
/* L43: */
	    	}
L10300:;
// printf("L10300:;\n");
	    	if(k > ns) {
				goto L10330;
	    	}
	    	for (i__ = k; i__ <= ns; ++i__) {
				mpfr_set_d(xr, 0e0, MPFR_RNDN);
				mpfr_set_d(xi, 0e0, MPFR_RNDN);
				for (j = k; j <= i__; ++j) {
					l = j - k + 1;
					mpfr_mul(tmp0, evr[l], ar[j + i__ * hr_dim1], MPFR_RNDN);
					mpfr_mul(tmp1, evi[l], ai[j + i__ * hi_dim1], MPFR_RNDN);
					mpfr_add(xr, xr, tmp0, MPFR_RNDN);
					mpfr_sub(xr, xr, tmp1, MPFR_RNDN);
/* L44: */
					mpfr_mul(tmp0, evr[l], ai[j + i__ * hi_dim1], MPFR_RNDN);
					mpfr_mul(tmp1, evi[l], ar[j + i__ * hr_dim1], MPFR_RNDN);
					mpfr_add(xi, xi, tmp0, MPFR_RNDN);
					mpfr_add(xi, xi, tmp1, MPFR_RNDN);
				}
				mpfr_set(ar[k + i__ * hr_dim1], xr, MPFR_RNDN);
				mpfr_set(ai[k + i__ * hi_dim1], xi, MPFR_RNDN);
/* L45: */
	    	}
L10330:;
// printf("L10330:;\n");
	    	if(nf > k) {
				goto L48;
	    	}
	    	for (i__ = nf; i__ <= k; ++i__) {
				if(k1 > ns) {
		    		goto L10360;
				}
				for (j = k1; j <= ns; ++j) {
					l = j - k + 1;
					mpfr_mul(tmp0, evi[l], ai[i__ + k * hi_dim1], MPFR_RNDN);
					mpfr_mul(tmp1, evr[l], ar[i__ + k * hr_dim1], MPFR_RNDN);
					mpfr_sub(cmp0, tmp0, tmp1, MPFR_RNDN);
					mpfr_div(cmp0, cmp0, rsum, MPFR_RNDN);
					mpfr_add(ar[i__ + j * hr_dim1], cmp0, ar[i__ + j * hr_dim1], MPFR_RNDN);
					mpfr_mul(tmp0, evr[l], ai[i__ + k * hi_dim1], MPFR_RNDN);
					mpfr_mul(tmp1, evi[l], ar[i__ + k * hr_dim1], MPFR_RNDN);
					mpfr_add(cmp0, tmp0, tmp1, MPFR_RNDN);
					mpfr_div(cmp0, cmp0, rsum, MPFR_RNDN);
					mpfr_sub(ai[i__ + j * hi_dim1], ai[i__ + j * hi_dim1], cmp0, MPFR_RNDN);
/* L46: */
				}
L10360:;
// printf("L10360:;\n");
				mpfr_div(ar[i__ + k * hr_dim1], ar[i__ + k * hr_dim1], rsum, MPFR_RNDN);
				mpfr_div(ai[i__ + k * hi_dim1], ai[i__ + k * hi_dim1], rsum, MPFR_RNDN);
/* L47: */
	   		}
L48:;
// printf("L48:;\n");
;
		}
L49:;
// printf("L49:;\n");
;
    }

	if (debug < 0 || debug >= 60) {	
		printf("----- STEP 5 OUTPUT\n");
		printf("matrix ar\n");
		for (int i=1; i<N+1; i++) {
			for (int j=1; j<N+1; j++) {
				mpfr_out_str(stdout, 10, 0, ar[i + hr_dim1*j], MPFR_RNDN);
				printf("\n");
			}
			printf("\n");
		}

		printf("matrix ai\n");
		for (int i=1; i<N+1; i++) {
			for (int j=1; j<N+1; j++) {
				mpfr_out_str(stdout, 10, 0, ai[i + hr_dim1*j], MPFR_RNDN);
				printf("\n");
			}
			printf("\n");
		}

		// printf("matrix zr\n");
		for (int i=1; i<N+1; i++) {
			for (int j=1; j<N+1; j++) {
				// mpfr_out_str(stdout, 10, 0, zr[i + hr_dim1*j], MPFR_RNDN);
				// printf("\n");
			}
			// printf("\n");
		}

		// printf("matrix zi\n");
		for (int i=1; i<N+1; i++) {
			for (int j=1; j<N+1; j++) {
				// mpfr_out_str(stdout, 10, 0, zi[i + hr_dim1*j], MPFR_RNDN);
				// printf("\n");
			}
			// printf("\n");
		}
		//getchar();
	}


    /***************/
    /*     STEP 6. */
    /***************/

    /*   ****************************************************************** */
    /*   *    TO EACH DIAGONAL BLOCK,CORRESPONDING TO A NUMERICAL */
    /*   *    MULTIPLE EIGENVALUE, FIND  A NILPOTENT MATRIX(I.E */
    /*   *    DETERMINE THE STRUCTURE OF THE BLOCK) */
    /*   *    USE  RQ-DECOMPOSITIONS OR CSVD-DECOMPOSITIONS(SEE THE */
    /*   *    ALGORITHM). */
    /*   *    NOTE THE STRUCTURE IN NDEFL. */
    /*   ****************************************************************** */
	idef = 1;
	ndefl[1] = 1;

	/*   *****  TOLERANCE GAP IS C2/C1 IN THE DEFLATION PROCESS      ***** */
	mpfr_set_d(c1, 1e0, MPFR_RNDN);
	mpfr_set_d(c2, 1e6, MPFR_RNDN);
	// mpfr_set_d(c1, 1e0, MPFR_RNDN);
	// mpfr_set_d(c2, 1e10, MPFR_RNDN);
	ns = 0;
	int last_call, reallc = 1;
	// printf("nblock = %d\n", nblock);
    for (ibl = 1; ibl <= nblock; ++ibl) {
		id = 0;
		nf = ns + 1;
		
		// printf("ibl = %d\n", ibl);
		// printf("ndel = %d\n", ndel[ibl+1]);
		ns = ndel[ibl + 1];
		nfs = nf;
		mpfr_set_d(dele[ibl], 0, MPFR_RNDN);
		ndb[ibl] = idef;

		// printf("ibl = %d\n", ibl);
		// printf("nf = %d\n", nf);
		// printf("ns = %d\n", ns);



		mpfr_set(shtr, cshtr[ibl - 1], MPFR_RNDN);
		mpfr_set(shti, cshti[ibl - 1], MPFR_RNDN);
		for (i__ = nfs; i__ <= ns; ++i__) {
			mpfr_sub(ar[i__ + i__ * hr_dim1], ar[i__ + i__ * hr_dim1], shtr, MPFR_RNDN);
/* L50: */
			mpfr_sub(ai[i__ + i__ * hi_dim1], ai[i__ + i__ * hi_dim1], shti, MPFR_RNDN);
		}

		ndef = 1;
		if(nf == ns) {
	    	goto L10385;
		}
L51:;
// printf("L51:;\n");
		if(nf > ns) {
			goto L52;
		}

		kb = ns - nf + 1;
		if (debug < 0 || debug >= 65) {	
			printf("----- RDEFL INPUT\n");
			printf("matrix ar\n");
			for (int i=1; i<N+1; i++) {
				for (int j=1; j<N+1; j++) {
					mpfr_out_str(stdout, 10, 0, ar[i + hr_dim1*j], MPFR_RNDN);
					printf("\n");
				}
				printf("\n");
			}
		}
		if (ibl == nblock) {
			last_call = 1;
		} else if (ibl < nblock) {
			last_call = 0;
		}
		rdefl_(&nm, &N, &kb, &ndef, &ar[nf + nf * hr_dim1], &ai[nf + nf * 
			hi_dim1], &zr[nf * zr_dim1 + 1], &zi[nf * zi_dim1 + 1], &told,
			&id, &c1, &c2, &dele[ibl], last_call, reallc);
		reallc = 0;
		if (debug < 0 || debug >= 65) {	
			printf("----- RDEFL OUTPUT\n");
			printf("matrix ar\n");
			for (int i=1; i<N+1; i++) {
				for (int j=1; j<N+1; j++) {
					mpfr_out_str(stdout, 10, 0, ar[i + hr_dim1*j], MPFR_RNDN);
					printf("\n");
				}
				printf("\n");
			}
			// getchar();
		}
L10385:;
// printf("L10385:;\n");
		id = id + ndef;
		nf = nf + ndef;
		++idef;
		ndefl[idef] = nf;
		goto L51;
L52:;
// printf("L52:;\n");
		for (i__ = nfs; i__ <= ns; ++i__) {
			mpfr_add(ar[i__ + i__ * hr_dim1], ar[i__ + i__ * hr_dim1], shtr, MPFR_RNDN);
/* L53: */
			mpfr_add(ai[i__ + i__ * hi_dim1], ai[i__ + i__ * hi_dim1], shti, MPFR_RNDN);
		}

		mpfr_sqrt(dele[ibl], dele[ibl], MPFR_RNDN);

/* L54: */
    }


	if (debug < 0 || debug >= 70) {	
		printf("----- STEP 6 OUTPUT\n");
		printf("matrix ar\n");
		for (int i=1; i<N+1; i++) {
			for (int j=1; j<N+1; j++) {
				mpfr_out_str(stdout, 10, 0, ar[i + hr_dim1*j], MPFR_RNDN);
				printf("\n");
			}
			printf("\n");
		}

		printf("matrix ai\n");
		for (int i=1; i<N+1; i++) {
			for (int j=1; j<N+1; j++) {
				mpfr_out_str(stdout, 10, 0, ai[i + hr_dim1*j], MPFR_RNDN);
				printf("\n");
			}
			printf("\n");
		}

		printf("matrix zr\n");
		for (int i=1; i<N+1; i++) {
			for (int j=1; j<N+1; j++) {
				mpfr_out_str(stdout, 10, 0, zr[i + hr_dim1*j], MPFR_RNDN);
				printf("\n");
			}
			printf("\n");
		}

		printf("matrix zi\n");
		for (int i=1; i<N+1; i++) {
			for (int j=1; j<N+1; j++) {
				mpfr_out_str(stdout, 10, 0, zi[i + hr_dim1*j], MPFR_RNDN);
				printf("\n");
			}
			printf("\n");
		}
		// getchar();
	}

    /***************/
    /*     STEP 7. */
    /***************/
    
    /*   ****************************************************************** */
    /*   *    TO EACH DIAGONAL BLOCK,CORRESPONDING TO A NUMERICAL MULTIPLE */
    /*   *    EIGENVALUE,COMPUTE THE COUPLING ELEMENTS OF THE PRINCIPAL */
    /*   *    CHAINS. */
    /*   *    THIS IS MADE BY STABILIZED ELIMINATION MATRICES. */
    /*   *    (SEE THE ALGORITHM) */
    /*   ****************************************************************** */
// printf("nblock = %d\n", nblock);
// printf("ndel before STEP 7:\n");
// for (i1 = 1; i1 <= nblock; ++i1) {
// printf("ndel[%d] = %d\n", i1, ndel[i1]);
// printf("ndb[%d] = %d\n", i1, ndb[i1]);
// printf("ndb[%d] = %d\n", i1, ndb[i1]);
// }
	ns = 0;
    for (i1 = 1; i1 <= nblock; ++i1) {
		nf = ns + 1;
		// printf("nblock = %d\n", nblock);
		// printf("i1 = %d\n", i1);
		// printf("ndel = %d\n", ndel[i1+1]);
		ns = ndel[i1 + 1];
		if(nf >= ns) {
			goto L78;
		}
		mpfr_set(shtr, cshtr[i1 - 1], MPFR_RNDN);
		mpfr_set(shti, cshti[i1 - 1], MPFR_RNDN);
		for (i__ = nf; i__ <= ns; ++i__) {
			mpfr_sub(ar[i__ + i__ * hr_dim1], ar[i__ + i__ * hr_dim1], shtr, MPFR_RNDN);
/* L55: */
			mpfr_sub(ai[i__ + i__ * hi_dim1], ai[i__ + i__ * hi_dim1], shti, MPFR_RNDN);
		}
		nkf = ns + 1;
		j1 = ndb[i1 + 1];
		if(j1 > 0) {
			goto L56;
		}
		j1 = idef;
L56:;
// printf("L56:;\n");
		nf1 = ndb[i1];
		// printf("nf1 = %d\n", nf1);
		// printf("j1 prima = %d\n", j1);
		j1 = j1 - nf1 - 1;
		// printf("j1 dopo = %d\n", j1);

		/*   ******  J1= THE NUMBER OF T(K,K+1)- BLOCKS(SEE THE ALGORITHM) **** */
		if(1 > j1) {
			goto L10420;
		}
		for (j2 = 1; j2 <= j1; ++j2) {
			// printf("#####################################################\n");
			// printf("#####################################################\n");
			// printf("#####################################################\n");
			j3 = j1 - j2 + nf1 + 1;
			nrs = ndefl[j3] - 1;
			nrf = ndefl[j3 - 1];
			// printf("j2 = %d\n", j2);
			// printf("nrs = %d\n", nrs);
			// printf("nrf = %d\n", nrf);
			nks = nkf - 1;
			nkf = nrs + 1;
			// printf("nks = %d\n", nks);
			// printf("nkf = %d\n", nkf);
			ip = nrs;
			// printf("ip = %d\n", ip);
        	/*   *****  PASS 1  IN  THE  ELIMINATION PROCESS                    ***** */

			if(nkf > nks) {
				goto L10430;
			}
			for (k = nkf; k <= nks; ++k) {
				// printf("-----------------------------------------\n");
				// printf("k = %d\n", k);
				// printf("nks = %d\n", nks);
				// printf("nkf = %d\n", nkf);
				// printf("ip = %d\n", ip);
				/*   *****  DETERMINE THE MAXIMUM OF ABS(HR(I,KP)) +ABS(HI(I,KP)) */
				/*          WHERE  I=NRF,NRF+1,...,IP-1 */

				kp = nks - k + nkf;
				im = ip;
				ip1 = ip - 1;
				mpfr_set(xr, ar[ip + kp * hr_dim1], MPFR_RNDN);
				mpfr_set(xi, ai[ip + kp * hi_dim1], MPFR_RNDN);
				if(nrf > ip1) {
					goto L10440;
				}
				for (i__ = nrf; i__ <= ip1; ++i__) {
					mpfr_abs(r__1, ar[i__ + kp * hr_dim1], MPFR_RNDN);
					mpfr_abs(r__2, ai[i__ + kp * hi_dim1], MPFR_RNDN);
					mpfr_add(r__1, r__1, r__2, MPFR_RNDN);
					mpfr_abs(tmp0, xr, MPFR_RNDN);
					mpfr_abs(tmp1, xi, MPFR_RNDN);
					mpfr_add(r__2, tmp0, tmp1, MPFR_RNDN);
					if(mpfr_lessequal_p(r__1, r__2)) {
						goto L57;
					}

					mpfr_set(xr, ar[i__ + kp * hr_dim1], MPFR_RNDN);
					mpfr_set(xi, ai[i__ + kp * hi_dim1], MPFR_RNDN);
					im = i__;
L57:;
// printf("L57:;\n");
;
				}
L10440:;
// printf("L10440:;\n");
				if(im == ip) {
					goto L61;
				}

				/*   *****  INTERCHANGE THE ROWS  IP AND IM OF HR AND HI             **** */

				if(nkf > kp) {
					goto L10450;
				}
				for (j = nkf; j <= kp; ++j) {
					mpfr_set(xr, ar[ip + j * hr_dim1], MPFR_RNDN);
					mpfr_set(ar[ip + j * hr_dim1], ar[im + j * hr_dim1], MPFR_RNDN);
					mpfr_set(ar[im + j * hr_dim1], xr, MPFR_RNDN);
					mpfr_set(xi, ai[ip + j * hi_dim1], MPFR_RNDN);
					mpfr_set(ai[ip + j * hi_dim1], ai[im + j * hi_dim1], MPFR_RNDN);
/* L58: */
					mpfr_set(ai[im + j * hi_dim1], xi, MPFR_RNDN);
				}
L10450:;
// printf("L10450:;\n");

        		/*   *****  INTERCHANGE  THE COLUMNS IP AND IM OF HR,JI AND ZR,ZI    **** */
				for (j = 1; j <= N; ++j) {

					if(j < nf || j >= nrf) {
						goto L59;
		    		}
					mpfr_set(xr, ar[j + ip * hr_dim1], MPFR_RNDN);
					mpfr_set(ar[j + ip * hr_dim1], ar[j + im * hr_dim1], MPFR_RNDN);
					mpfr_set(ar[j + im * hr_dim1], xr, MPFR_RNDN);
					mpfr_set(xi, ai[j + ip * hi_dim1], MPFR_RNDN);
					mpfr_set(ai[j + ip * hi_dim1], ai[j + im * hi_dim1], MPFR_RNDN);
					mpfr_set(ai[j + im * hi_dim1], xi, MPFR_RNDN);
L59:;
// printf("L59:;\n");
					mpfr_set(xr, zr[j + ip * zr_dim1], MPFR_RNDN);
					mpfr_set(zr[j + ip * zr_dim1], zr[j + im * zr_dim1], MPFR_RNDN);
					mpfr_set(zr[j + im * zr_dim1], xr, MPFR_RNDN);
					mpfr_set(xi, zi[j + ip * zi_dim1], MPFR_RNDN);
					mpfr_set(zi[j + ip * zi_dim1], zi[j + im * zi_dim1], MPFR_RNDN);
/* L60: */
					mpfr_set(zi[j + im * zi_dim1], xi, MPFR_RNDN);
				}

L61:;
// printf("L61:;\n");
				/*   *****  ELIMINATE  THE ELEMENTS  HR,HI(J,KP)  J=NF,NF+1,...,IP-1 */
				/*          PIVOTE ELEMENT IS  HR,HI(IP,KP)                         ***** */
				kp1 = kp - 1;
				mpfr_set(hrik, ar[ip + kp * hr_dim1], MPFR_RNDN);
				mpfr_set(hiik, ai[ip + kp * hi_dim1], MPFR_RNDN);

				if(nkf > kp1) {
					goto L10470;
				}
				for (l = nkf; l <= kp1; ++l) {
					if(nf > ip1) {
						goto L63;
					}
					for (j = nf; j <= ip1; ++j) {
						mpc_set_fr_fr(m, ar[j + l * hr_dim1], ai[j + l * hi_dim1], MPFR_RNDN);
						mpc_set_fr_fr(q__2, ar[j + kp * hr_dim1], ai[j + kp * hi_dim1], MPFR_RNDN);
						mpc_set_fr_fr(q__3, hrik,hiik, MPFR_RNDN);
						mpc_div(q__2,q__2,q__3, MPFR_RNDN);
						mpc_set_fr_fr(q__3, ar[ip + l * hr_dim1], ai[ip + l * hi_dim1], MPFR_RNDN);
						mpc_mul(q__2,q__2,q__3, MPFR_RNDN);
						mpc_sub(m,m,q__2, MPFR_RNDN);
						mpfr_set(ar[j + l * hr_dim1], mpc_realref(m), MPFR_RNDN);
/* L62: */
						mpfr_set(ai[j + l * hi_dim1], mpc_imagref(m), MPFR_RNDN);
		    		}
L63:;
// printf("L63:;\n");
;
				}
L10470:;
// printf("L10470:;\n");

				/*   *****  FINISH THE SIMILARITY TRANSFORMATION  AND ACCUMULATE */
				/*               TRANSFORMATIONS                                    ***** */
				for (j = 1; j <= N; ++j) {
					// printf("entrato nel for\n");
					// printf("(ip, j) = (%d, %d)\n", ip, j);
					// printf("index = %d\n", j + ip * hr_dim1);
					// printf("ar\n");
					// mpfr_out_str(stdout, 10, 0, ar[j + ip * hr_dim1], MPFR_RNDN);
					// printf("\n");
					// printf("ai\n");
					// mpfr_out_str(stdout, 10, 0, ai[j + ip * hr_dim1], MPFR_RNDN);
					// printf("\n");
					mpc_set_fr_fr(hri,ar[j + ip * hr_dim1], ai[j + ip * hi_dim1], MPFR_RNDN);
					// printf("fatta operazione 1\n");
					mpc_set_fr_fr(zri,zr[j + ip * hr_dim1], zi[j + ip * hi_dim1], MPFR_RNDN);
					// printf("assegnati valori a hri e zri\n");
				
					if(nf > ip1) {
						goto L10500;
					}
					for (l = nf; l <= ip1; ++l) {
						if(j < nf || j >= nrf) {
							goto L64;
						}
						// printf("entrato nel secondo for\n");
						mpc_set_fr_fr(q__2,ar[l + kp * hr_dim1], ai[l + kp * hi_dim1], MPFR_RNDN);
						// printf("fatta operazione 1\n");
						mpc_set_fr_fr(q__3, hrik,hiik, MPFR_RNDN);
						// printf("fatta operazione 2\n");
						mpc_div(q__2, q__2, q__3, MPFR_RNDN);
						// printf("fatta operazione 3\n");
						mpc_set_fr_fr(q__3,ar[j + l * hr_dim1], ai[j + l * hi_dim1], MPFR_RNDN);
						// printf("fatta operazione 4\n");
						mpc_mul(q__4,q__2,q__3, MPFR_RNDN);
						// printf("fatta operazione 5\n");
						mpc_add(hri,hri,q__4, MPFR_RNDN);
						// printf("fatta operazione 6\n");
L64:;
// printf("L64:;\n");
						mpc_set_fr_fr(q__2,ar[l + kp * hr_dim1], ai[l + kp * hi_dim1], MPFR_RNDN);
						mpc_set_fr_fr(q__3, hrik,hiik, MPFR_RNDN);
						mpc_div(q__2, q__2, q__3, MPFR_RNDN);
						mpc_set_fr_fr(q__3,zr[j + l * zr_dim1], zi[j + l * zi_dim1], MPFR_RNDN);
						mpc_mul(q__4,q__2,q__3, MPFR_RNDN);
						mpc_add(zri,zri,q__4, MPFR_RNDN);
/* L65: */
   					}   
L10500:;
// printf("L10500:;\n");
   					if(j < nf || j >= nrf) {
     					goto L66;
   					}
					mpfr_set(ar[j + ip * hr_dim1], mpc_realref(hri), MPFR_RNDN);
					mpfr_set(ai[j + ip * hi_dim1], mpc_imagref(hri), MPFR_RNDN);
L66:;
// printf("L66:;\n");
					mpfr_set(zr[j + ip * zr_dim1], mpc_realref(zri), MPFR_RNDN);
					mpfr_set(zi[j + ip * zi_dim1], mpc_imagref(zri), MPFR_RNDN);
/* L67: */
				}
				 --ip;
/* L69: */
	    	}
L10430:;
// printf("L10430:;\n");
    
			/*   *****   END  PASS  1                                           ***** */

			/*   *****  PASS 2  IN THE ELIMINATION - PROCESS                     **** */
			if(nkf > nks) {
				goto L76;
			}
	    	for (kp = nkf; kp <= nks; ++kp) {
				++ip;

        		/*   *****  ELIMINATE THE ELEMTS HR,HI(J,KP),J=IP+1,...,NRS          **** */
				ip1 = ip + 1;

				/*   *****  FINISH THE SIMILARITY TRANSFORMATION AND  ACCUMULATE */
				/*          TRANSFORMATIONS                                          **** */
				mpfr_set(hrik, ar[ip + kp * hr_dim1], MPFR_RNDN);
				mpfr_set(hiik, ai[ip + kp * hi_dim1], MPFR_RNDN);

        		/* Computing 2nd power */
				mpfr_set(r__1, hrik, MPFR_RNDN);
        		/* Computing 2nd power */
				mpfr_set(r__2, hiik, MPFR_RNDN);
				mpfr_mul(r__1, r__1, r__1, MPFR_RNDN);
				mpfr_mul(r__2, r__2, r__2, MPFR_RNDN);
				mpfr_add(r__1, r__1, r__2, MPFR_RNDN);
				mpfr_sqrt(b, r__1, MPFR_RNDN);
				mpfr_div(t1, hrik, b, MPFR_RNDN);
				mpfr_div(tmp0, hiik, b, MPFR_RNDN);
				mpfr_set(t2, tmp0, MPFR_RNDN);
				mpfr_neg(t2, t2, MPFR_RNDN);
 				for (j = 1; j <= N; ++j) {
					mpc_set_fr_fr(hri,ar[j + ip * hr_dim1], ai[j + ip * hi_dim1], MPFR_RNDN);
					mpc_set_fr_fr(zri,zr[j + ip * zr_dim1], zi[j + ip * zi_dim1], MPFR_RNDN);
					if(ip1 > nrs) {
						goto L10530;
					}
					for (l = ip1; l <= nrs; ++l) {
						if(j < nf || j >= nrf) {
							goto L70;
						}
						mpc_set_fr_fr(q__4,ar[l + kp * hr_dim1], ai[l + kp * hi_dim1], MPFR_RNDN);
						mpc_set_fr_fr(q__5, hrik,hiik, MPFR_RNDN);
						mpc_div(q__4, q__4, q__5, MPFR_RNDN);
						mpc_set_fr_fr(q__5,ar[j + l * hr_dim1], ai[j + l * hi_dim1], MPFR_RNDN);
						mpc_mul(q__4,q__4,q__5, MPFR_RNDN);
						mpc_add(hri,hri,q__4, MPFR_RNDN);
L70:;
// printf("L70:;\n");
						mpc_set_fr_fr(q__4,ar[l + kp * hr_dim1], ai[l + kp * hi_dim1], MPFR_RNDN);
						mpc_set_fr_fr(q__5, hrik,hiik, MPFR_RNDN);
						mpc_div(q__4, q__4, q__5, MPFR_RNDN);
						mpc_set_fr_fr(q__5,zr[j + l * zr_dim1], zi[j + l * zi_dim1], MPFR_RNDN);
						mpc_mul(q__4,q__4,q__5, MPFR_RNDN);
						mpc_add(zri,zri,q__4, MPFR_RNDN);
/* L71: */
		    		}
L10530:;
// printf("L10530:;\n");
					if(j < nf || j >= nrf) {
						goto L72;
					}
					mpc_set_fr_fr(q__2,t1,t2, MPFR_RNDN);
					mpc_div(hri, hri, q__2, MPFR_RNDN);
					mpfr_set(ar[j + ip * hr_dim1], mpc_realref(hri), MPFR_RNDN);
					mpfr_set(ai[j + ip * hi_dim1], mpc_imagref(hri), MPFR_RNDN);
L72:;
// printf("L72:;\n");
					mpc_set_fr_fr(q__2,t1,t2, MPFR_RNDN);
					mpc_div(zri, zri, q__2, MPFR_RNDN);
					mpfr_set(zr[j + ip * zr_dim1], mpc_realref(zri), MPFR_RNDN);
					mpfr_set(zi[j + ip * zi_dim1], mpc_imagref(zri), MPFR_RNDN);
/* L73: */
				}

        		/*     MAKE THE (IP,KP)- ELEMENT REAL */
				mpfr_mul(tmp0, hrik, t1, MPFR_RNDN);
				mpfr_mul(tmp1, hiik, t2, MPFR_RNDN);
				mpfr_sub(ar[ip + kp * hr_dim1], tmp0, tmp1, MPFR_RNDN);
				mpfr_mul(tmp0, hrik, t2, MPFR_RNDN);
				mpfr_mul(tmp1, hiik, t1, MPFR_RNDN);
				mpfr_add(ai[ip + kp * hi_dim1], tmp0, tmp1, MPFR_RNDN);
				nxt[ip] = kp;
				mpfr_set(supd[ip], ar[ip + kp * hr_dim1], MPFR_RNDN);
/* L75: */
	    	}
        	/*   ****   END PASS  2                                              **** */

        	/*   ***   A  NEW  T(K,K+1) - BLOCK                                  **** */
L76:;
// printf("L76:;\n");
;
		}
L10420:;
// printf("L10420:;\n");
		for (i__ = nf; i__ <= ns; ++i__) {
			mpfr_add(ar[i__ + i__ * hr_dim1], ar[i__ + i__ * hr_dim1], shtr, MPFR_RNDN);
			mpfr_add(ai[i__ + i__ * hi_dim1], ai[i__ + i__ * hi_dim1], shti, MPFR_RNDN);
/* L77: */
		}

    	/*   ****   A  NEW  EIGENVALUE                                       **** */
L78:;
// printf("L78:;\n");
;
    }
L10400:;
// printf("L10400:;\n");
	i__1 = N;
    for (i__ = 1; i__ <= i__1; ++i__) {
		mpfr_set(evr[i__], ar[i__ + i__ * hr_dim1], MPFR_RNDN);
/* L79: */
		mpfr_set(evi[i__], ai[i__ + i__ * hi_dim1], MPFR_RNDN);
    }



	/* CONSTRUCTION OF THE JORDAN NORMAL FORM */
	// find index of eingenvalues (i.e. index not appearing into nxt)
	
	// initialize to zero elements outsite diagonal and super-diagonal
	for (int i=0; i<N; i++) {
		for (int j=0; j<N; j++) {
			if (j != i && j!= i+1) {
				mpc_set_d_d(jor[i][j], 0, 0, MPFR_RNDN);
			}
		}
	}

    int *zeros = (int*) malloc(N*sizeof(int));
    int *inv = (int*) malloc(N*sizeof(int));
    int *evec = (int*) malloc(N*sizeof(int));
    int *perm = (int*) malloc(N*sizeof(int));

    for (int i=0; i<N; i++) {
		zeros[i]=0;
		inv[i]=0;
		evec[i]=0;
    }
      
    
    k = 0;
    for (int j=1; j<=N; j++) {
		int i = j-1;
		if (nxt[i+1] == 0) {
			zeros[k] = j;
			k++;
		} else {
			inv[nxt[i+1]-1] = j;
      	}
    }
    
    k = 0;
    for (int j=1; j<=N; j++) {
    	int i = j-1;
    	if (inv[i] == 0) {
			evec[k] = j;
			k++;
    	}
    }
	// number of pure eingenvecotors
	*num_ein = k;
    
        
	// building the permutation
    mpfr_t *factors = (mpfr_t*) malloc(N*sizeof(mpfr_t));
    for (int i=0; i<N; i++) {
    	mpfr_init2(factors[i], wp2);
    }

    int restart1 = 1;
    int restart2 = 0;
    int h=-1;
    for (int j=1; j<=N; j++) {
		int i = j-1;
		if (evec[i] == 0) {
			break;
		}
		h++;
		perm[h] = evec[i];
		mpc_set_fr_fr(jor[perm[h]-1][perm[h]-1], evr[h+1], evi[h+1], MPFR_RNDN);
		mpfr_set_d(factors[h], 1, MPFR_RNDN);
		if (h>0) {
			mpc_set_d_d(jor[h-1][h], 0, 0, MPFR_RNDN);
		}
		restart1 = 0;
		restart2 = 0;
		while (1)  {
			if (nxt[perm[h]] == 0) {
				restart1 = 1;
				break;
			} else {
				h++;
				perm[h]=nxt[perm[h-1]];
				mpc_set_fr_fr(jor[perm[h]-1][perm[h]-1], evr[h+1], evi[h+1], MPFR_RNDN);
				mpfr_mul(factors[h], factors[h-1], supd[perm[h-1]], MPFR_RNDN);
				if (h>0) {
					mpc_set_d_d(jor[h-1][h], 1, 0, MPFR_RNDN);
				}
				restart2 = 1;
			}
			if (restart1 == 0 && restart2 == 0) {
				break;
			}
		}
    }

	// applying permutation and rescaling
    for (int i=1; i<N+1; i++) {
        for (int j=1; j<N+1; j++) {
			mpc_set_fr_fr(t_mat[i-1][j-1], zr[i + hr_dim1*perm[j-1]], zi[i + hi_dim1*perm[j-1]], MPFR_RNDN);
			mpc_div_fr(t_mat[i-1][j-1], t_mat[i-1][j-1], factors[j-1], MPFR_RNDN);
        }
    }

	// get indices how the eingenvectors
	k = 0;
	*ein_pos = (int*) malloc(*num_ein*sizeof(int));
	for (int i=0; i<*num_ein; i++) {
		// printf("i = %d\n", i);
		(*ein_pos)[i] = perm[evec[i] - 1] - 1;
	}


	if (debug < 0 || debug >= 80) {
		/***********/
		/* RESULTS */
		/***********/
		printf("----- RESULTS\n");
		printf("matrix ar\n");
		for (int i=1; i<N+1; i++) {
			for (int j=1; j<N+1; j++) {
				mpfr_out_str(stdout, 10, 0, ar[i + hr_dim1*j], MPFR_RNDN);
				printf("\n");
			}
			printf("\n");
		}

		printf("matrix ai\n");
		for (int i=1; i<N+1; i++) {
			for (int j=1; j<N+1; j++) {
				mpfr_out_str(stdout, 10, 0, ai[i + hr_dim1*j], MPFR_RNDN);
				printf("\n");
			}
			printf("\n");
		}

		printf("matrix zr\n");
		for (int i=1; i<N+1; i++) {
			for (int j=1; j<N+1; j++) {
				mpfr_out_str(stdout, 10, 0, zr[i + hr_dim1*j], MPFR_RNDN);
				printf("\n");
			}
			printf("\n");
		}

		printf("matrix zi\n");
		for (int i=1; i<N+1; i++) {
			for (int j=1; j<N+1; j++) {
				mpfr_out_str(stdout, 10, 0, zi[i + hr_dim1*j], MPFR_RNDN);
				printf("\n");
			}
			printf("\n");
		}

		printf("vector nxt\n");
		printf("{");
		for (int i=1; i<N+1; i++) {
			printf("%d, ", nxt[i]);
		}
		printf("}\n");

		printf("vector supd\n");
		for (int i=1; i<N+1; i++) {
			mpfr_out_str(stdout, 10, 0, supd[i], MPFR_RNDN);
			printf("\n");
		}
		printf("}\n");

		printf("vector evr\n");
		for (int i=1; i<N+1; i++) {
			mpfr_out_str(stdout, 10, 0, evr[i], MPFR_RNDN);
			printf("\n");
		}
		printf("}\n");
		
		printf("vector evi\n");
		for (int i=1; i<N+1; i++) {
			mpfr_out_str(stdout, 10, 0, evi[i], MPFR_RNDN);
			printf("\n");
		}
		printf("}\n");

		printf("vector ndefl\n");
		printf("{");
		for (int i=1; i<N+1; i++) {
			printf("%d, ", ndefl[i]);
		}
		printf("}\n");

		printf("vector ndel\n");
		printf("{");
		for (int i=1; i<N+1; i++) {
			printf("%d, ", ndel[i]);
		}
		printf("}\n");

		printf("vector ndb\n");
		printf("{");
		for (int i=1; i<N+1; i++) {
			printf("%d, ", ndb[i]);
		}
		printf("}\n");

		printf("zr={");
		for (int i=1; i<N+1; i++) {
			printf("{");
			for (int j=1; j<N+1; j++) {
				mpfr_out_str(stdout, 10, 0, zr[i + hr_dim1*j], MPFR_RNDN);
				printf(", ");
			}
			printf("}, ");
		}
		printf("}\n");


		printf("zi={");
		for (int i=1; i<N+1; i++) {
			printf("{");
			for (int j=1; j<N+1; j++) {
				mpfr_out_str(stdout, 10, 0, zi[i + hi_dim1*j], MPFR_RNDN);
				printf(", ");
			}
			printf("}, ");
		}
		printf("}\n");

		printf("vector factors\n");
		for (int i=0; i<N; i++) {
		mpfr_out_str(stdout, 10, 0, factors[i], MPFR_RNDN);
		printf("\n");
		}
		printf("}\n");

		printf("perm:\n");
		for (int j=1; j<=N; j++) {
			printf("%d ", perm[j-1]);
		}
		printf("\n");

		printf("num eig = %d\n", *num_ein);
		printf("evec: \n");
		for (int i=0; i<*num_ein; i++) {
			printf("%d, ", evec[i]);
		}
		printf("\n");

		printf("ein_pos: \n");
		for (int i=0; i<*num_ein; i++) {
			printf("%d, ", (*ein_pos)[i]);
		}
		printf("\n");
	}
}
