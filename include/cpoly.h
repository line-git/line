#ifndef CPOLY_H
#define CPOLY_H

void print_poly(mpc_t *, int);

int ipow(int , int);

void copy_poly(mpc_t *out, mpc_t *in, int deg);

void init_poly(mpc_t *, int);

void set_null_poly(mpc_t *, int);

void add_poly(mpc_t *, mpc_t *, mpc_t *, int, int);

void mul_poly(mpc_t *, mpc_t *, mpc_t *, int, int);

void mul_poly_mpq(mpq_t *out_pol, mpq_t *pol1, mpq_t *pol2, int deg1, int deg2);

void mul_poly_mpq_wrp(mpq_t *out_pol, mpq_t *pol1, mpq_t *pol2, int deg1, int deg2);

void mul_poly_mpq_tot_wrp(mpq_t *out_pol, mpq_t *pol1, mpq_t *pol2, int deg1, int deg2);

void div_poly(mpc_t *, mpc_t *, mpc_t *, int, int);

int div_poly_mpq(mpq_t *out_pol, mpq_t *pol1, mpq_t *pol2, int deg1, int deg2);

void poly_shift(mpc_t cnst, mpc_t *d, int n);

// void poly_shift_cp(mpc_t *out, mpc_t cnst, mpc_t *d, int n);

void poly_mul_root(mpc_t *coeffs, int deg, mpc_t root);

void poly_mul_root_cp(mpc_t *out, mpc_t *coeffs, int deg, mpc_t root);

#endif