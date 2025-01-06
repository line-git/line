#ifndef ALGEBRA_H
#define ALGEBRA_H

void mp_ludcmp(mpc_t **a, int n, int *indx, float *d);

void mp_lubksb(mpc_t **a, int n, int *indx, mpc_t *b);

void mp_inverse(mpc_t **a_inv, mpc_t **a, int n);

#endif