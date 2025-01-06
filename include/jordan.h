#ifndef JORDAN_H
#define JORDAN_H

int cbal_(int *nm, int *n, mpfr_t *ar, mpfr_t *ai, int *low, int *igh, mpfr_t *scale);

int comhes_(int *nm, int *n, int *low, int *igh, mpfr_t *ar, mpfr_t *ai, int *int__);

int comlr2_(int *nm, int *n, int *low, int *igh, int *int__, mpfr_t *hr, mpfr_t *hi, mpfr_t *wr, mpfr_t *wi, mpfr_t *zr, mpfr_t *zi, int *ierr, int *nobak);

int cbabk2_(int *nm, int *n, int *low, int *igh, mpfr_t *scale, int *m, mpfr_t *zr, mpfr_t *zi);

int decide_(int *nm, int *n, int *ndel, mpfr_t *cshtr, mpfr_t *cshti, int *nblock, mpfr_t *hr, mpfr_t *hi);

int csvd1_(mpc_t *a, int *mmax, int *nmax, int *m, int *n, int *p, int *nu, int *nv, mpfr_t *s, mpc_t *u, mpc_t *v, int N, int reallc);

int rdefl_(int *nm, int *n, int *kb, int *ndef, mpfr_t *hr, mpfr_t *hi, mpfr_t *zr, mpfr_t *zi, mpfr_t *told, int *d__, mpfr_t *c1, mpfr_t *c2, mpfr_t *dl, int last_call, int reallc);

void jordan(mpc_t **jor, mpc_t **t_mat, int **ein_pos, int *num_ein, mpfr_t *ar, mpfr_t *ai, int N, int nm, int debug);

#endif