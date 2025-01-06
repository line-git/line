#ifndef REL_ERR_MPC_H
#define REL_ERR_MPC_H

void exp_rel_err_mpc_add(
  mpc_t out, mpc_t in1, mpc_t in2, mpfr_rnd_t round, int wp_bin
);

void exp_rel_err_mpc_mul(
  mpc_t out, mpc_t in1, mpc_t in2, mpfr_rnd_t round, int wp_bin
);

void exp_rel_err_mpc_div(
  mpc_t out, mpc_t in1, mpc_t in2, mpfr_rnd_t round, int wp_bin
);

void rel_err_mpfr_add(
  mpfr_t out, mpfr_t in1, mpfr_t in2, mpfr_rnd_t round 
);

void rel_err_mpfr_sub(
  mpfr_t out, mpfr_t in1, mpfr_t in2, mpfr_rnd_t round 
);

void rel_err_mpc_add(
  mpc_t out, mpc_t in1, mpc_t in2, mpfr_rnd_t round 
);

void rel_err_mpc_sub(
  mpc_t out, mpc_t in1, mpc_t in2, mpfr_rnd_t round 
);

void rel_err_mpfr_fma(
  mpfr_t out, mpfr_t in0, mpfr_t in1, mpfr_t in2, mpfr_rnd_t round
);

void rel_err_mpfr_fmma(
  mpfr_t out, mpfr_t in0, mpfr_t in1, mpfr_t in2, mpfr_t in3, mpfr_rnd_t round
);

void rel_err_mpc_fma(
  mpc_t out, mpc_t in0, mpc_t in1, mpc_t in2, mpfr_rnd_t round
);

#endif