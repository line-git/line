#include <iostream>
#include <cstring>
// #include <fstream>
// #include <sstream>
// #include <string>
// #include <cmath>
// #include <complex>
#include "mpc.h"

using namespace std;

#include "global_vars.h"
#include "utils.h"
#include "setup.h"
#include "tensor_utils.h"
#include "conversions.h"
#include "malloc_defs.h"


void mpc_parg(mpfr_t *out, mpc_t *in) {
  mpc_arg(*out, *in, MPFR_RNDN);
  if (mpfr_sign_within_tol(*out) < 0) {
    mpfr_t mpfr_2pi;
    mpfr_init2(mpfr_2pi, wp2);
    mpfr_const_pi(mpfr_2pi, MPFR_RNDN);
    mpfr_mul_ui(mpfr_2pi, mpfr_2pi, 2, MPFR_RNDN);
    mpfr_add(*out, *out, mpfr_2pi, MPFR_RNDN);
  }
}


void mpc_plog(mpc_t *out, mpc_t *in) {
  mpfr_t tmpc;
  mpfr_init2(tmpc, wp2);
  mpc_parg(&mpc_imagref(*out), in);
  mpc_parg(&tmpc, in);
  mpc_abs(mpc_realref(*out), *in, MPFR_RNDN);
  mpfr_log(mpc_realref(*out), mpc_realref(*out), MPFR_RNDN);
}


void hypergeom(
  // OUTPUT
  mpc_t *out,
  // INPUT
  int na, int nb, mpc_t *a, mpc_t *b, mpc_t z
) {
  mpc_t tmp, term;
  mpc_init3(tmp, wp2, wp2);
  mpc_init3(term, wp2, wp2);

  mpfr_t tmpfr;
  mpfr_init2(tmpfr, wp2);

  int k, j = 1;
  mpc_set_ui(term, 1, MPFR_RNDN);
  mpc_set_ui(*out, 1, MPFR_RNDN);
  while(1) {
    // update term
    for (k=0; k<na; k++) {
      mpc_mul(term, term, a[k], MPFR_RNDN);
    }
    for (k=0; k<nb; k++) {
      mpc_div(term, term, b[k], MPFR_RNDN);
    }
    mpc_mul(term, term, z, MPFR_RNDN);
    mpc_div_ui(term, term, j, MPFR_RNDN);

    // cumulate partial sum
    mpc_add(*out, *out, term, MPFR_RNDN);

    // update arguments
    for (k=0; k<na; k++) {
      mpc_add_ui(a[k], a[k], 1, MPFR_RNDN);
    }
    for (k=0; k<nb; k++) {
      mpc_add_ui(b[k], b[k], 1, MPFR_RNDN);
    }
    j++;
    
    // estimate error
    mpfr_mul(tmpfr, mpfr_tol, mpc_realref(*out), MPFR_RNDN);
    if (mpfr_cmpabs(mpc_realref(term), tmpfr) < 0) {
      if (mpfr_zero_p(mpc_imagref(*out))) {
        break;
      }
      mpfr_mul(tmpfr, mpfr_tol, mpc_imagref(*out), MPFR_RNDN);
      if (mpfr_cmpabs(mpc_imagref(term), tmpfr) < 0) {
        break;
      }
    }

    if (j == 10000) {
      perror("hypergeometric function not converging fast enough");
      // exit(1);
      mpc_set_ui(*out, 0, MPFR_RNDN);
      return;
    }
  }

  // cout << "num terms = " << j << endl;
}


void formulary(mpc_t *out, int loop, int id, mpc_t eps, mpc_t *scale, int *intargs) {
  int print = 0;
  if (loop == 1) {
    mpc_t f_eps;
    mpc_init3(f_eps, wp2, wp2);
    mpfr_t gamma_f_eps;
    mpfr_init2(gamma_f_eps, wp2);

    switch(id) {
      case -1:
        // tapdole of mass m2 with n-1 dots
        //
        // (-1)^n m2^(2 - n - eps)*Gamma[-2 + n + eps]/Gamma[n];
        //
        mpc_set(f_eps, eps, MPFR_RNDN);  // eps
        mpc_add_si(f_eps, f_eps, -2+intargs[0], MPFR_RNDN);  // eps+2-n

        mpfr_gamma(gamma_f_eps, mpc_realref(f_eps), MPFR_RNDN);  // Gamma[eps+2-n]
        mpc_set_fr(*out, gamma_f_eps, MPFR_RNDN);

        // 1/Gamma[n]
        mpc_set_ui(f_eps, intargs[0], MPFR_RNDN);  // n
        mpfr_gamma(gamma_f_eps, mpc_realref(f_eps), MPFR_RNDN);  // Gamma[n]
        mpc_div_fr(*out, *out, gamma_f_eps, MPFR_RNDN);

        // E^((2-n-eps)*pLog[m2])
        mpc_neg(f_eps, eps, MPFR_RNDN);  // -eps
        mpc_add_si(f_eps, f_eps, 2-intargs[0], MPFR_RNDN);  // 2-n-eps
        mpc_pow(f_eps, *scale, f_eps, MPFR_RNDN);
        mpc_mul(*out, *out, f_eps, MPFR_RNDN);

        // (-1)^n
        if (intargs[0]%2 == 1) {
          mpc_neg(*out, *out, MPFR_RNDN);
        }
        break;

      case 0:
        // tapdole of mass m2
        //
        // - m2^(1 - eps)*Gamma[-1 + eps];
        //
        // cout << "compute tadpole" << endl;
        mpc_t epsm1;
        mpc_init3(epsm1, wp2, wp2);
        mpc_set(epsm1, eps, MPFR_RNDN);
        mpc_sub_ui(epsm1, epsm1, 1, MPFR_RNDN);

        // mpfr_t gamma_eps;
        // mpfr_init2(gamma_eps, wp2);
        // mpfr_gamma(gamma_eps, mpc_realref(eps), MPFR_RNDN);

        mpfr_t gamma_epsm1;
        mpfr_init2(gamma_epsm1, wp2);
        mpfr_gamma(gamma_epsm1, mpc_realref(epsm1), MPFR_RNDN);
        // mpc_set_fr(*out, gamma_eps, MPFR_RNDN);
        // mpc_mul_fr(*out, *out, gamma_epsm1, MPFR_RNDN);
        mpc_set_fr(*out, gamma_epsm1, MPFR_RNDN);

        mpc_neg(*out, *out, MPFR_RNDN);

        // E^((1-eps)*pLog[m2])
        mpc_neg(epsm1, epsm1, MPFR_RNDN);  // 1-eps
        mpc_pow(epsm1, *scale, epsm1, MPFR_RNDN);
        mpc_mul(*out, *out, epsm1, MPFR_RNDN);

        // mpc_t tmpc;
        // mpc_init3(tmpc, wp2, wp2);
        // mpc_plog(&tmpc, scale);  // pLog[m2]
        // mpc_mul(tmpc, tmpc, epsm1, MPFR_RNDN);  // (1-eps)*pLog[m2]
        // mpc_exp(tmpc, tmpc, MPFR_RNDN);  // Exp[(1-eps)*pLog[m2]]
        // mpc_mul(*out, *out, tmpc, MPFR_RNDN);
        break;

      case 1:
        // massless bubble with squared momentum s
        //
        // (Gamma[1 - eps]^2 Gamma[eps])/Gamma[2 - 2 eps]*E^(I*Pi*eps)*E^(-eps*pLog[s]);
        //
        // cout << "compute bubble" << endl;

        // Gamma[eps]
        mpfr_gamma(gamma_f_eps, mpc_realref(eps), MPFR_RNDN);  // Gamma[eps]
        mpc_set_fr(*out, gamma_f_eps, MPFR_RNDN);
        // cout << "gamma = "; print_mpc(out); cout << endl;

        // Gamma[1-eps]^2
        mpc_neg(f_eps, eps, MPFR_RNDN);  // -eps
        mpc_add_ui(f_eps, f_eps, 1, MPFR_RNDN);  // 1-eps
        mpfr_gamma(gamma_f_eps, mpc_realref(f_eps), MPFR_RNDN);  // Gamma[1-eps]
        mpfr_sqr(gamma_f_eps, gamma_f_eps, MPFR_RNDN);  // Gamma[1-eps]^2
        // cout << "gamma = "; mpfr_out_str(stdout, 10, 0, gamma_f_eps, MPFR_RNDN); cout << endl;
        mpc_mul_fr(*out, *out, gamma_f_eps, MPFR_RNDN);
        // cout << "gamma = "; print_mpc(out); cout << endl;

        // 1/Gamma[2-2*eps]
        mpc_mul_ui(f_eps, f_eps, 2, MPFR_RNDN);  // 2-2*eps
        mpfr_gamma(gamma_f_eps, mpc_realref(f_eps), MPFR_RNDN);  // Gamma[2-2*eps]
        mpc_div_fr(*out, *out, gamma_f_eps, MPFR_RNDN);
        // cout << "gamma product = "; print_mpc(out); cout << endl;

        // E^(I*Pi*eps)
        mpfr_set_ui(mpc_realref(f_eps), 0, MPFR_RNDN);
        mpfr_const_pi(mpc_imagref(f_eps), MPFR_RNDN);
        mpc_mul(f_eps, f_eps, eps, MPFR_RNDN);
        mpc_exp(f_eps, f_eps, MPFR_RNDN);
        mpc_mul(*out, *out, f_eps, MPFR_RNDN);

        // E^(-eps*pLog[sf])
        mpc_plog(&f_eps, scale);
        mpc_mul(f_eps, f_eps, eps, MPFR_RNDN);
        mpc_neg(f_eps, f_eps, MPFR_RNDN);
        mpc_exp(f_eps, f_eps, MPFR_RNDN);
        mpc_mul(*out, *out, f_eps, MPFR_RNDN);
        break;

    }
  } else if (loop == 2) {
    mpc_t tmpc;
    mpc_init3(tmpc, wp2, wp2);
    mpc_t f_eps;
    mpc_init3(f_eps, wp2, wp2);
    mpfr_t f_epsr;
    mpfr_init2(f_epsr, wp2);
    mpfr_t gamma_f_eps;
    mpfr_init2(gamma_f_eps, wp2);

    if (print) cout << "called formulary with id = " << id << endl;
    switch(id) {
      case 0: {
        // Gamma[1-eps]^2
        mpc_neg(f_eps, eps, MPFR_RNDN);  // -eps
        mpc_add_ui(f_eps, f_eps, 1, MPFR_RNDN);  // 1-eps
        mpfr_gamma(gamma_f_eps, mpc_realref(f_eps), MPFR_RNDN);  // Gamma[1-eps]
        mpfr_pow_ui(gamma_f_eps, gamma_f_eps, 3, MPFR_RNDN);  // Gamma[1-eps]^3
        // cout << "gamma = "; mpfr_out_str(stdout, 10, 0, gamma_f_eps, MPFR_RNDN); cout << endl;
        mpc_set_fr(*out, gamma_f_eps, MPFR_RNDN);
        // cout << "gamma = "; print_mpc(out); cout << endl;

        // 1/Gamma[3-3*eps]
        mpc_mul_ui(f_eps, f_eps, 3, MPFR_RNDN);  // 3-3*eps
        mpfr_gamma(gamma_f_eps, mpc_realref(f_eps), MPFR_RNDN);  // Gamma[3-3*eps]
        mpc_div_fr(*out, *out, gamma_f_eps, MPFR_RNDN);

        mpc_mul_ui(f_eps, eps, 2, MPFR_RNDN);  // 2*eps
        mpc_sub_ui(f_eps, f_eps, 1, MPFR_RNDN);  // 2*eps-1
        mpfr_gamma(gamma_f_eps, mpc_realref(f_eps), MPFR_RNDN);  // Gamma[2*eps-1]
        mpc_mul_fr(*out, *out, gamma_f_eps, MPFR_RNDN);

        // E^(I*Pi*(2*eps-1))
        mpfr_set_ui(mpc_realref(tmpc), 0, MPFR_RNDN);
        mpfr_const_pi(mpc_imagref(tmpc), MPFR_RNDN);
        mpc_mul(tmpc, tmpc, f_eps, MPFR_RNDN);
        mpc_exp(tmpc, tmpc, MPFR_RNDN);
        mpc_mul(*out, *out, tmpc, MPFR_RNDN);

        // E^((1-2*eps)*pLog[sf])
        mpc_plog(&tmpc, scale);
        mpc_neg(f_eps, f_eps, MPFR_RNDN);  // 1-2*eps
        mpc_mul(tmpc, tmpc, f_eps, MPFR_RNDN);
        mpc_exp(tmpc, tmpc, MPFR_RNDN);
        mpc_mul(*out, *out, tmpc, MPFR_RNDN);

        mpc_neg(*out, *out, MPFR_RNDN);

        break;
      }
      case 1: {
        // Gamma[-eps]
        mpc_neg(f_eps, eps, MPFR_RNDN);  // -eps
        mpfr_gamma(gamma_f_eps, mpc_realref(f_eps), MPFR_RNDN);  // Gamma[-eps]
        mpc_set_fr(*out, gamma_f_eps, MPFR_RNDN);

        // Gamma[2-3*eps]
        mpc_mul_ui(f_eps, f_eps, 3, MPFR_RNDN);  // -3*eps
        mpc_add_ui(f_eps, f_eps, 2, MPFR_RNDN);  // 2-3*eps
        mpfr_gamma(gamma_f_eps, mpc_realref(f_eps), MPFR_RNDN);  // Gamma[2-3*eps]
        mpc_div_fr(*out, *out, gamma_f_eps, MPFR_RNDN);

        // Pi
        mpfr_const_pi(mpc_realref(f_eps), MPFR_RNDN);  // Pi
        mpfr_set_ui(mpc_imagref(f_eps), 0, MPFR_RNDN);  // Pi + I*0
        mpc_mul(*out, *out, f_eps, MPFR_RNDN);

        // eps*Pi
        mpc_mul(f_eps, f_eps, eps, MPFR_RNDN);  // eps*Pi
        mpc_mul(*out, *out, f_eps, MPFR_RNDN);

        // Csc[eps*Pi]^2        
        mpfr_csc(gamma_f_eps, mpc_realref(f_eps), MPFR_RNDN);  // Csc[eps*Pi]
        mpfr_sqr(gamma_f_eps, gamma_f_eps, MPFR_RNDN);  // Csc[eps*Pi]^2
        mpc_mul_fr(*out, *out, gamma_f_eps, MPFR_RNDN);

        // Sec[eps*Pi]
        mpfr_sec(gamma_f_eps, mpc_realref(f_eps), MPFR_RNDN);  // Sec[eps*Pi]
        mpc_mul_fr(*out, *out, gamma_f_eps, MPFR_RNDN);

        // 1/(2*eps-1)/2
        mpc_mul_ui(f_eps, eps, 2, MPFR_RNDN);  // 2*eps
        mpc_sub_ui(f_eps, f_eps, 1, MPFR_RNDN);  // 2*eps-1
        mpc_mul_ui(f_eps, f_eps, 2, MPFR_RNDN);  // (2*eps-1)*2
        mpc_div(*out, *out, f_eps, MPFR_RNDN);

        // E^(2*I*Pi*eps)
        mpfr_set_ui(mpc_realref(f_eps), 0, MPFR_RNDN);  // I*Pi
        mpfr_const_pi(mpc_imagref(f_eps), MPFR_RNDN);   // 0 + I*Pi
        mpc_mul(f_eps, f_eps, eps, MPFR_RNDN);  // I*Pi*eps
        mpc_mul_ui(f_eps, f_eps, 2, MPFR_RNDN);  // 2*I*Pi*eps
        mpc_exp(f_eps, f_eps, MPFR_RNDN);  // E^(2*I*Pi*eps)
        mpc_mul(*out, *out, f_eps, MPFR_RNDN);

        // E^(-2*eps*pLog[sf])
        mpc_plog(&f_eps, scale);
        mpc_mul(f_eps, f_eps, eps, MPFR_RNDN);
        mpc_mul_si(f_eps, f_eps, -2, MPFR_RNDN);
        mpc_exp(f_eps, f_eps, MPFR_RNDN);
        mpc_mul(*out, *out, f_eps, MPFR_RNDN);

        break;
      }
      case 2: {
        mpc_set_str(
          *out,
          // "(3.1684445908121786003796188869681207469262054567966336188473972321905767552446233224778214280384741205391846215310141631793166119945813680783406630695853741439471614362681273652610578621515303334130420e7 0)",
          "(9.48492441003640825260962197129821386821935925179745358655703451954439836015317916332124781308352224437848919330244005516616880421890377423138789483293981556e13 -1.00947753147919615613456932874428726657221737865114652664327271646661672479748071744091520831603494308707860478643749428735527122265757502131690243917754981e-86)",
          10,
          MPFR_RNDN
          );
        break;
      }
      case 3: {
        mpc_set_str(
          *out,
          "(1.0556885073537226376204154581819559565162485520340765226718881571635022449216909215333020217081031872224999407879008893705082364943347751840882546266282524561005532543000156205486419738105082371684575e7 0)",
          10,
          MPFR_RNDN
          );
        break;
      }
      case 4: {
        if (print) {
        cout << "intargs[0] = " << intargs[0] << endl;
        cout << "intargs[1] = " << intargs[1] << endl;
        cout << "intargs[2] = " << intargs[2] << endl;
        }
        mpc_t tmpc;
        mpc_init3(tmpc, wp2, wp2);

        mpfr_t f_epsr1, f_epsr2, f_epsr3;
        mpfr_init2(f_epsr1, wp2);
        mpfr_init2(f_epsr2, wp2);
        mpfr_init2(f_epsr3, wp2);

        mpc_t z;
        mpc_init3(z, wp2, wp2);
        mpc_set_ui(z, 1, MPFR_RNDN);
        mpc_div_ui(z, z, 4, MPFR_RNDN);

        mpc_t *a1 = new mpc_t[4];
        mpc_t *a2 = new mpc_t[4];
        mpc_t *b1 = new mpc_t[3];
        mpc_t *b2 = new mpc_t[3];
        init_rk1_mpc(a1, 4);
        init_rk1_mpc(a2, 4);
        init_rk1_mpc(b1, 3);
        init_rk1_mpc(b2, 3);

        mpc_set_si(a1[1], intargs[0], MPFR_RNDN);
        mpc_set_si(a1[2], intargs[1], MPFR_RNDN);
        mpc_set_si(b1[0], intargs[0]+intargs[1], MPFR_RNDN);
        mpc_add_si(b1[1], b1[0], 1, MPFR_RNDN);
        mpc_div_ui(b1[0], b1[0], 2, MPFR_RNDN);
        mpc_div_ui(b1[1], b1[1], 2, MPFR_RNDN);
        mpc_set_si(a2[0], intargs[2], MPFR_RNDN);

        //////
        // GAMMA FACTORS
        //////
        mpfr_t gamma_fact1;
        mpfr_init2(gamma_fact1, wp2);
        mpfr_t gamma_fact2;
        mpfr_init2(gamma_fact2, wp2);

        // 1/Gamma[2-eps]
        mpfr_set_ui(f_epsr, 2, MPFR_RNDN);  // 2
        mpfr_sub(f_epsr, f_epsr, mpc_realref(eps), MPFR_RNDN);  // 2-eps
        mpc_set_fr(a1[0], f_epsr, MPFR_RNDN);
        mpfr_gamma(gamma_fact2, f_epsr, MPFR_RNDN);
        
        mpfr_neg(f_epsr2, f_epsr, MPFR_RNDN);  // eps-2
        mpfr_add_si(f_epsr3, f_epsr2, intargs[0]+intargs[1], MPFR_RNDN);  // eps-2+n1+n2
        mpc_set_fr(a1[3], f_epsr2, MPFR_RNDN);
        mpc_add_si(a1[3], a1[3], intargs[0]+intargs[1], MPFR_RNDN);
        mpfr_mul_ui(f_epsr2, f_epsr2, 2, MPFR_RNDN);  // 2*eps-4
        
        // Gamma[2-eps-n3]
        mpfr_add_si(f_epsr, f_epsr, -intargs[2], MPFR_RNDN);  // 2-eps-n3
        mpc_set_fr(b1[2], f_epsr, MPFR_RNDN);
        mpc_add_ui(b1[2], b1[2], 1, MPFR_RNDN);  // 3-eps-n3
        mpc_set_fr(b2[0], f_epsr, MPFR_RNDN);  
        mpc_neg(b2[0], b2[0], MPFR_RNDN);  // -2+eps+n3
        mpc_add_ui(b2[0], b2[0], 1, MPFR_RNDN);  // -1+eps+n3
        mpfr_gamma(gamma_f_eps, f_epsr, MPFR_RNDN);
        mpfr_div(gamma_fact2, gamma_f_eps, gamma_fact2, MPFR_RNDN);

        // Gamma[-2+eps+n3]
        mpfr_neg(f_epsr1, f_epsr, MPFR_RNDN);  // -2+eps+n3
        mpfr_gamma(gamma_fact1, f_epsr1, MPFR_RNDN);

        // Gamma[-2+eps+n3+n1]
        mpfr_add_si(f_epsr, f_epsr1, intargs[0], MPFR_RNDN);  // -2+eps+n3+n1
        mpc_set_fr(a2[1], f_epsr, MPFR_RNDN);
        mpfr_gamma(gamma_f_eps, f_epsr, MPFR_RNDN);
        mpfr_mul(gamma_fact2, gamma_fact2, gamma_f_eps, MPFR_RNDN);

        // Gamma[-2+eps+n3+n2]
        mpfr_add_si(f_epsr, f_epsr1, intargs[1], MPFR_RNDN);  // -2+eps+n3+n2
        mpc_set_fr(a2[2], f_epsr, MPFR_RNDN);
        mpfr_gamma(gamma_f_eps, f_epsr, MPFR_RNDN);
        mpfr_mul(gamma_fact2, gamma_fact2, gamma_f_eps, MPFR_RNDN);

        // Gamma[-2+eps+n1+n2]
        mpfr_gamma(gamma_f_eps, f_epsr3, MPFR_RNDN);
        mpfr_mul(gamma_fact1, gamma_fact1, gamma_f_eps, MPFR_RNDN);

        // Gamma[-4+2*eps+n1+n2+n3]
        mpfr_add_si(f_epsr2, f_epsr2, intargs[0]+intargs[1]+intargs[2], MPFR_RNDN);  // -4+2*eps+n1+n2+n3
        mpc_set_fr(a2[3], f_epsr2, MPFR_RNDN);
        mpfr_gamma(gamma_f_eps, f_epsr2, MPFR_RNDN);
        mpfr_mul(gamma_fact2, gamma_fact2, gamma_f_eps, MPFR_RNDN);

        // 1/Gamma[-4+2*eps+n1+n2+2*n3]
        mpfr_add_si(f_epsr2, f_epsr2, intargs[2], MPFR_RNDN);
        mpc_set_fr(b2[1], f_epsr2, MPFR_RNDN);
        mpc_add_ui(b2[2], b2[1], 1, MPFR_RNDN);
        mpc_div_ui(b2[1], b2[1], 2, MPFR_RNDN);
        mpc_div_ui(b2[2], b2[2], 2, MPFR_RNDN);
        mpfr_gamma(gamma_f_eps, f_epsr2, MPFR_RNDN);
        mpfr_div(gamma_fact2, gamma_fact2, gamma_f_eps, MPFR_RNDN);

        // 1/Gamma[n1]
        mpfr_set_si(f_epsr, intargs[0], MPFR_RNDN);
        mpfr_gamma(gamma_f_eps, f_epsr, MPFR_RNDN);
        mpfr_div(gamma_fact2, gamma_fact2, gamma_f_eps, MPFR_RNDN);

        // 1/Gamma[n2]
        mpfr_set_si(f_epsr, intargs[1], MPFR_RNDN);
        mpfr_gamma(gamma_f_eps, f_epsr, MPFR_RNDN);
        mpfr_div(gamma_fact2, gamma_fact2, gamma_f_eps, MPFR_RNDN);

        // 1/Gamma[n1+n2]
        mpfr_set_si(f_epsr, intargs[0]+intargs[1], MPFR_RNDN);
        mpfr_gamma(gamma_f_eps, f_epsr, MPFR_RNDN);
        mpfr_div(gamma_fact1, gamma_fact1, gamma_f_eps, MPFR_RNDN);

        // 1/Gamma[n3]
        mpfr_set_si(f_epsr, intargs[2], MPFR_RNDN);
        mpfr_gamma(gamma_f_eps, f_epsr, MPFR_RNDN);
        mpfr_div(gamma_fact1, gamma_fact1, gamma_f_eps, MPFR_RNDN);

        //////
        // HYPERGEOMETRIC FUNCTIONS 2F1
        //////
        // cout << "gamma_fact1:" << endl;
        // mpfr_out_str(stdout, 10, 0, gamma_fact1, MPFR_RNDN); cout << endl;
        // cout << "gamma_fact2:" << endl;
        // mpfr_out_str(stdout, 10, 0, gamma_fact2, MPFR_RNDN); cout << endl;

        // cout << "hg1 args:" << endl;
        // cout << "a:" << endl;
        // print_poly(a1, 3);
        // cout << "b:" << endl;
        // print_poly(b1, 2);
        // cout << "z:" << endl;
        // print_mpc(&z); cout << endl;

        // cout << "hg2 args:" << endl;
        // cout << "a:" << endl;
        // print_poly(a2, 3);
        // cout << "b:" << endl;
        // print_poly(b2, 2);
        // cout << "z:" << endl;
        // print_mpc(&z); cout << endl;

        hypergeom(out, 4, 3, a1, b1, z);
        // cout << "hg1 = "; print_mpc(out); cout << endl;
        mpc_mul_fr(*out, *out, gamma_fact1, MPFR_RNDN);

        hypergeom(&tmpc, 4, 3, a2, b2, z);
        // cout << "hg2 = "; print_mpc(&tmpc); cout << endl;
        mpc_mul_fr(tmpc, tmpc, gamma_fact2, MPFR_RNDN);
        mpc_add(*out, *out, tmpc, MPFR_RNDN);

        if ((intargs[0]+intargs[1]+intargs[2])%2 == 1) {
          mpc_neg(*out, *out, MPFR_RNDN);
        }

        // E^((4-n1-n2-n3-2eps)*pLog[m2])
        mpc_set(f_eps, eps, MPFR_RNDN);  // eps
        mpc_mul_si(f_eps, f_eps, -2, MPFR_RNDN);  // -2*eps
        mpc_add_si(f_eps, f_eps, 4-intargs[0]-intargs[1]-intargs[2], MPFR_RNDN);  // -2*eps+4-n1-n2-n3
        mpc_pow(f_eps, *scale, f_eps, MPFR_RNDN);
        if (print) {cout << "scale factor = "; print_mpc(&f_eps); cout << endl;}
        mpc_mul(*out, *out, f_eps, MPFR_RNDN);

        break;
      }
      case 5: {
        // find negative integer
        int neg_int, m, n;
        if (intargs[0] < 0) {
          neg_int = intargs[0];
          m = intargs[1];
          n = intargs[2];
        } else if (intargs[1] < 0) {
          m = intargs[0];
          neg_int = intargs[1];
          n = intargs[2];
        } else if (intargs[2] < 0) {
          m = intargs[0];
          n = intargs[1];
          neg_int = intargs[2];
        }
        
        if (print) cout << "m = " << m << endl;
        if (print) cout << "n = " << n << endl;
        if (print) cout << "neg_int = " << neg_int << endl;
        if (neg_int == -1) {
          // Gamma[-3+m+eps]
          mpc_add_si(f_eps, eps, -3+m, MPFR_RNDN);  // -3+m+eps
          mpfr_gamma(gamma_f_eps, mpc_realref(f_eps), MPFR_RNDN);  // Gamma[-3+m+eps]
          mpc_set_fr(*out, gamma_f_eps, MPFR_RNDN);
          
          // Gamma[-3+n+eps]
          mpc_add_si(f_eps, eps, -3+n, MPFR_RNDN);  // -3+n+eps
          mpfr_gamma(gamma_f_eps, mpc_realref(f_eps), MPFR_RNDN);  // Gamma[-3+n+eps]
          mpc_mul_fr(*out, *out, gamma_f_eps, MPFR_RNDN);

          // 1/Gamma[m]
          mpc_set_ui(f_eps, m, MPFR_RNDN);  // m
          mpfr_gamma(gamma_f_eps, mpc_realref(f_eps), MPFR_RNDN);  // Gamma[m]
          mpc_div_fr(*out, *out, gamma_f_eps, MPFR_RNDN);

          // 1/Gamma[n]
          mpc_set_ui(f_eps, n, MPFR_RNDN);  // n
          mpfr_gamma(gamma_f_eps, mpc_realref(f_eps), MPFR_RNDN);  // Gamma[n]
          mpc_div_fr(*out, *out, gamma_f_eps, MPFR_RNDN);

          // 3+(-4+eps)*eps+m+n-m*n
          mpc_set(f_eps, eps, MPFR_RNDN);  // eps
          mpc_add_si(f_eps, f_eps, -4, MPFR_RNDN);  // -4+eps
          mpc_mul(f_eps, f_eps, eps, MPFR_RNDN);  // (-4+eps)*eps
          mpc_add_si(f_eps, f_eps, 3+m+n-m*n, MPFR_RNDN);  // (-4+eps)*eps+3+m+n-m*n
          mpc_mul(*out, *out, f_eps, MPFR_RNDN);

          // (-1)^(m+n)
          if ((m+n)%2 == 1) {
            mpc_neg(*out, *out, MPFR_RNDN);
          }
          
          // E^((5-m-n-2eps)*pLog[m2])
          mpc_set(f_eps, eps, MPFR_RNDN);  // eps
          mpc_mul_si(f_eps, f_eps, -2, MPFR_RNDN);  // -2*eps
          mpc_add_si(f_eps, f_eps, 5-m-n, MPFR_RNDN);  // -2*eps+5-m-n
          mpc_pow(f_eps, *scale, f_eps, MPFR_RNDN);
          if (print) {cout << "scale factor = "; print_mpc(&f_eps); cout << endl;}
          mpc_mul(*out, *out, f_eps, MPFR_RNDN);

        } else if (neg_int == -2) {
          // perror("hypergeometric with negative index -2");
          // mpc_set_ui(*out, 0, MPFR_RNDN);

          // Gamma[eps+m-2]
          mpc_add_si(f_eps, eps, m-2, MPFR_RNDN);  // eps+m-2
          mpfr_gamma(gamma_f_eps, mpc_realref(f_eps), MPFR_RNDN);  // Gamma[eps+m-2]
          mpc_set_fr(*out, gamma_f_eps, MPFR_RNDN);

          // 1/(eps+m-3)
          mpc_add_si(f_eps, f_eps, -1, MPFR_RNDN);  // eps+m-3
          mpc_div(*out, *out, f_eps, MPFR_RNDN);

          // 1/(eps+m-4)
          mpc_add_si(f_eps, f_eps, -1, MPFR_RNDN);  // eps+m-4
          mpc_div(*out, *out, f_eps, MPFR_RNDN);

          // Gamma[eps+n-4]
          mpc_add_si(f_eps, eps, n-4, MPFR_RNDN);  // eps+n-4
          mpfr_gamma(gamma_f_eps, mpc_realref(f_eps), MPFR_RNDN);  // Gamma[eps+n-4]
          mpc_mul_fr(*out, *out, gamma_f_eps, MPFR_RNDN);

          // 1/Gamma[m]
          mpc_set_ui(f_eps, m, MPFR_RNDN);  // m
          mpfr_gamma(gamma_f_eps, mpc_realref(f_eps), MPFR_RNDN);  // Gamma[m]
          mpc_div_fr(*out, *out, gamma_f_eps, MPFR_RNDN);

          // 1/Gamma[n]
          mpc_set_ui(f_eps, n, MPFR_RNDN);  // n
          mpfr_gamma(gamma_f_eps, mpc_realref(f_eps), MPFR_RNDN);  // Gamma[n]
          mpc_div_fr(*out, *out, gamma_f_eps, MPFR_RNDN);

          // (-1)^(m+n)
          if ((m+n)%2 == 1) {
            mpc_neg(*out, *out, MPFR_RNDN);
          }

          // 96+m^2*(-2+n)*(-1+n)+2*(-7+n)*n+m*(-14+(5-3*n)*n) +
          // 2*(-58+2*m+2*n+3*m*n)*eps +
          // (55-2*m*n)*eps^2 +
          // -12*eps^3 +
          // eps^4
          mpc_add_si(f_eps, eps, -12, MPFR_RNDN);  // -12+eps
          mpc_mul(f_eps, f_eps, eps, MPFR_RNDN);  // -12*eps+eps^2
          mpc_add_si(f_eps, f_eps, 55-2*m*n, MPFR_RNDN);  // (55-2*m*n)-12*eps+eps^2
          mpc_mul(f_eps, f_eps, eps, MPFR_RNDN);  // (55-2*m*n)*eps-12*eps^2+eps^3
          mpc_add_si(f_eps, f_eps, 2*(-58+2*m+2*n+3*m*n), MPFR_RNDN);  // 2*(-58*+2*m+2*n+3*m*n)+(55-2*m*n)*eps-12*eps^2+eps^3
          mpc_mul(f_eps, f_eps, eps, MPFR_RNDN);  // 2*(-58*+2*m+2*n+3*m*n)*eps+(55-2*m*n)*eps^2-12*eps^3+eps^4
          mpc_add_si(f_eps, f_eps, 96+m*m*(-2+n)*(-1+n)+2*(-7+n)*n+m*(-14+(5-3*n)*n), MPFR_RNDN);  // 96+m^2*(-2+n)*(-1+n)+2(-7+n)*n+m*(-14+(5-3*n)*n)+2*(-58*+2*m+2*n+3*m*n)*eps+(55-2*m*n)*eps^2-12*eps^3+eps^4
          mpc_mul(*out, *out, f_eps, MPFR_RNDN);

          // E^((6-m-n-2eps)*pLog[m2])
          mpc_set(f_eps, eps, MPFR_RNDN);  // eps
          mpc_mul_si(f_eps, f_eps, -2, MPFR_RNDN);  // -2*eps
          mpc_add_si(f_eps, f_eps, 6-m-n, MPFR_RNDN);  // -6*eps+2-m-n
          mpc_pow(f_eps, *scale, f_eps, MPFR_RNDN);
          if (print) {cout << "scale factor = "; print_mpc(&f_eps); cout << endl;}
          mpc_mul(*out, *out, f_eps, MPFR_RNDN);

          break;
        } else if (neg_int == -3) {
          // Gamma[eps+m-2]
          mpc_add_si(f_eps, eps, m-5, MPFR_RNDN);
          mpfr_gamma(gamma_f_eps, mpc_realref(f_eps), MPFR_RNDN);
          mpc_set_fr(*out, gamma_f_eps, MPFR_RNDN);

          // Gamma[eps+n-5]
          mpc_add_si(f_eps, eps, n-5, MPFR_RNDN);
          mpfr_gamma(gamma_f_eps, mpc_realref(f_eps), MPFR_RNDN);
          mpc_mul_fr(*out, *out, gamma_f_eps, MPFR_RNDN);

          // 1/Gamma[m]
          mpc_set_ui(f_eps, m, MPFR_RNDN);
          mpfr_gamma(gamma_f_eps, mpc_realref(f_eps), MPFR_RNDN);
          mpc_div_fr(*out, *out, gamma_f_eps, MPFR_RNDN);

          // 1/Gamma[n]
          mpc_set_ui(f_eps, n, MPFR_RNDN);
          mpfr_gamma(gamma_f_eps, mpc_realref(f_eps), MPFR_RNDN);
          mpc_div_fr(*out, *out, gamma_f_eps, MPFR_RNDN);

          // (-1)^(m+n+1)
          if ((m+n) % 2 == 0) {
            mpc_neg(*out, *out, MPFR_RNDN);
          }

          // polynomial expression
          mpc_neg(f_eps, eps, MPFR_RNDN);
          mpc_add_si(f_eps, f_eps, 24, MPFR_RNDN);
          mpc_mul(f_eps, f_eps, eps, MPFR_RNDN);
          mpc_add_si(f_eps, f_eps, -238 + 3*n+ 3*m*(1 + n), MPFR_RNDN);
          mpc_mul(f_eps, f_eps, eps, MPFR_RNDN);
          mpc_add_si(f_eps, f_eps, -6*(-206 + 7*n + m*(7 + 5*n)), MPFR_RNDN);
          mpc_mul(f_eps, f_eps, eps, MPFR_RNDN);
          mpc_add_si(f_eps, f_eps, 5*(-701 + 39*n) + 3*m*(65 + n*(40 + m + n - m*n)), MPFR_RNDN);
          mpc_mul(f_eps, f_eps, eps, MPFR_RNDN);
          mpc_add_si(f_eps, f_eps, 6*(846 + m*m*(-1 + n)*(3 + n) + m*(1 + n)*(-53 + 2*n) - n*(53 + 3*n)), MPFR_RNDN);
          mpc_mul(f_eps, f_eps, eps, MPFR_RNDN);
          mpc_add_si(f_eps, f_eps, m*m*m*(-3 + n)*(-2 + n)*(-1 + n) - 6*m*m*(-1 + n)*(12 + (-5 + n)*n) - 6*(480 + (-13 + n)*n*(1 + n)) + m*(78 + n*(445 + n*(-102 + 11*n))), MPFR_RNDN);
          mpc_mul(*out, *out, f_eps, MPFR_RNDN);

          // E^((7-m-n-2eps)*pLog[m2])
          mpc_set(f_eps, eps, MPFR_RNDN);  // eps
          mpc_mul_si(f_eps, f_eps, -2, MPFR_RNDN);  // -2*eps
          mpc_add_si(f_eps, f_eps, 7-m-n, MPFR_RNDN);  // -2*eps+7-m-n
          mpc_pow(f_eps, *scale, f_eps, MPFR_RNDN);
          if (print) {cout << "scale factor = "; print_mpc(&f_eps); cout << endl;}
          mpc_mul(*out, *out, f_eps, MPFR_RNDN);

          break;
        } else {
          perror("hypergeometric with negative index lower than -3");
          exit(1);
        }

        break;
      }
      case 6: {
        if (intargs[0] < 0 || intargs[1] < 0 ||intargs[2] < 0) {
          // at least one current is negative
          if (intargs[0] < 0 && intargs[1] <  0 && intargs[2] < 0) {
            // all currents are negative
            perror("all of the loop currents are negative in boundary generation");
            exit(1);
          } else if (intargs[0]*intargs[1]*intargs[2] > 0) {
            // two currents are negative
            perror("two loop currents are negative in boundary generation");
            exit(1);
          } else if (intargs[0]*intargs[1]*intargs[2] < 0) {
            // only one current is negative
            formulary(out, 2, 5, eps, scale, intargs);
            return;
          }
        } else if (intargs[0] == 0 && intargs[1] == 0 && intargs[2] == 0) {
          // null boundary
          mpc_set_ui(*out, 0, MPFR_RNDN);
          return;
        } else if (intargs[0] == 0 || intargs[1] == 0 || intargs[2] == 0) {
          // at least one current is zero
          mpc_set_ui(*out, 1, MPFR_RNDN);
          if (intargs[0] != 0) {
            formulary(&tmpc, 1, -1, eps, scale, intargs);
            mpc_mul(*out, *out, tmpc, MPFR_RNDN);
          }
          if (intargs[1] != 0) {
            formulary(&tmpc, 1, -1, eps, scale, intargs+1);
            mpc_mul(*out, *out, tmpc, MPFR_RNDN);
          }
          if (intargs[2] != 0) {
            formulary(&tmpc, 1, -1, eps, scale, intargs+2);
            mpc_mul(*out, *out, tmpc, MPFR_RNDN);
          }
          return;
        } else {
          formulary(out, 2, 4, eps, scale, intargs);
          return;
        }
      }
      case 7: {
        // Sum[Binomial[n1,k1]*(-1)^(n1-k1)*(l1^2+1)^k1/(l1^2+1)^d1, {k1, 0, n1}] * (n1->n2,d1->d2,l1->l2) * (n1->n3,d1->d3,l1->l1+l2)
        int n1 = -intargs[3];
        int n2 = -intargs[4];
        int n3 = -intargs[5];
        int d1 = intargs[0];
        int d2 = intargs[1];
        int d3 = intargs[2];
        int fac1, fac2, fac3;
        int *inner_intargs = new int[3];
        if (print) {
        cout << "n1, n2, n3 = " << n1 << "," << n2 << "," << n3 << "," << endl;
        cout << "d1, d2, d3 = " << d1 << "," << d2 << "," << d3 << "," << endl;
        }
        mpc_set_ui(*out, 0, MPFR_RNDN);
        fac1 = 1;
        mpc_t inner_scale; mpc_init3(inner_scale, wp2, wp2);
        mpc_set_ui(inner_scale, 1, MPFR_RNDN);
        for (int k1=0; k1<=n1; k1++) {
          fac2 = 1;
          for (int k2=0; k2<=n2; k2++) {
            fac3 = 1;
            for (int k3=0; k3<=n3; k3++) {
              inner_intargs[0] = d1-k1;
              inner_intargs[1] = d2-k2;
              inner_intargs[2] = d3-k3;
              if (print) {
              cout << "k1, k2, k3 = " << k1 << "," << k2 << "," << k3 << "," << endl;
              cout << "fac1, fac2, fac3 = " << fac1 << "," << fac2 << "," << fac3 << "," << endl;
              cout << "d1-k1, d2-k2, d3-k3 = " << inner_intargs[0] << "," << inner_intargs[1] << "," << inner_intargs[2] << "," << endl;
              }
              formulary(&tmpc, 2, 6, eps, &inner_scale, inner_intargs);
              if (print) {
              cout << "term = "; print_mpc(&tmpc); cout << endl;
              }
              mpc_mul_ui(tmpc, tmpc, fac1*fac2*fac3, MPFR_RNDN);
              // if (n1-k1+n2-k2+n3-k3 % 2 == 1) {
              //   if (print) cout << "sign = -1" << endl;
              //   mpc_neg(tmpc, tmpc, MPFR_RNDN);
              // }
              mpc_add(*out, *out, tmpc, MPFR_RNDN);
              fac3 *= (n3-k3);
              fac3 /= k3+1;
            }
            fac2 *= (n2-k2);
            fac2 /= k2+1;
          }
          fac1 *= (n1-k1);
          fac1 /= k1+1;
        }

        // scale factor
        // E^((4-d1-d2-d3+n1+n2+n3-2eps)*pLog[m2])
        mpc_set(f_eps, eps, MPFR_RNDN);  // eps
        mpc_mul_si(f_eps, f_eps, -2, MPFR_RNDN);  // -2*eps
        mpc_add_si(f_eps, f_eps, 4-d1-d2-d3+n1+n2+n3, MPFR_RNDN);  // -2*eps+4-d1-d2-d3+n1+n2+n3
        mpc_pow(f_eps, *scale, f_eps, MPFR_RNDN);
        if (print) {cout << "scale factor = "; print_mpc(&f_eps); cout << endl;}
        if (print) {cout << "before scale factor = "; print_mpc(out); cout << endl;}
        mpc_mul(*out, *out, f_eps, MPFR_RNDN);
        if (print) {cout << "after scale factor = "; print_mpc(out); cout << endl;}
        
      }
    }
  }
}


void dispatcher(
  // OUTPUT
  mpc_t *bound,
  // INPUT
  // char *id_str, mpc_t eps,  int *is_mass, 
  char *eps_str, char ***kin, char **symbols, int ninvs,
  int dim, int *nadd, int **nfac, int ***loop, int ***id, char ****sym, char ****point, int ****intargs
) {
  int print = 0;
  mpc_t eps;
  mpc_init3(eps, wp2, wp2);
  mpc_set_str_rat(&eps, eps_str);
  mpc_t *PS_ini = new mpc_t[ninvs];
  init_rk1_mpc(PS_ini, ninvs);
  mpc_t *PS_fin = new mpc_t[ninvs];
  init_rk1_mpc(PS_fin, ninvs);
  str_rk1_to_mpc_rk1(PS_ini, kin[0], ninvs);
  str_rk1_to_mpc_rk1(PS_fin, kin[1], ninvs);
  mpc_t *addend = new mpc_t[2];
  init_rk1_mpc(addend, 2);
  
  for (int m=0; m<dim; m++) {
    if (print) cout << "m = " << m << endl;

    mpc_set_ui(bound[m], 0, MPFR_RNDN);
    for (int a=0; a<nadd[m]; a++) {
      if (print) cout << "a = " << a << endl;
      if (nfac[m][a] == 0) {
        mpc_set_ui(addend[a], 0, MPFR_RNDN);
        continue;
      }

      mpc_set_ui(addend[a], 1, MPFR_RNDN);
      mpc_t factor;
      mpc_init3(factor, wp2, wp2);
      for (int s, k=0; k<nfac[m][a]; k++) {
        if (print) cout << "k = " << k << endl;
        if (print) cout << sym[m][a][k] << endl;
        // get index of symbol
        s = sym_to_idx(sym[m][a][k], symbols+1, ninvs);
        if (print) cout << "s = " << s << endl;
        if (print) cout << "id = " << id[m][a][k] << endl;
        // get factor from formulary
        if (strcmp(point[m][a][k], "ini") == 0) {
          // cout << "call ini" << endl;
          formulary(&factor, loop[m][a][k], id[m][a][k], eps, &PS_ini[s], intargs[m][a][k]);
        } else if (strcmp(point[m][a][k], "fin") == 0) {
          // cout << "call fin" << endl;
          formulary(&factor, loop[m][a][k], id[m][a][k], eps, &PS_fin[s], intargs[m][a][k]);
        }
        if (print) {cout << "factor = "; print_mpc(&factor); cout << endl;}
        // cumulate product
        mpc_mul(addend[a], addend[a], factor, MPFR_RNDN);
        if (print) {cout << "cum prod = "; print_mpc(&addend[a]); cout << endl;}
      }
    }
    
    // cumulate
    if (nadd[m] == 0) {
      mpc_set_ui(bound[m], 0, MPFR_RNDN);
    } else if (nadd[m] == 1) {
      mpc_set(bound[m], addend[0], MPFR_RNDN);
    } else if (nadd[m] == 2) {
      mpc_sub(bound[m], addend[0], addend[1], MPFR_RNDN);  // #hardcoded: always subtract
    }

    if (print) {cout << "cum bound = "; print_mpc(&bound[m]); cout << endl;}
    
  }

}


void generate_boundaries(
  // OUTPUT
  mpc_t **bound,
  // INPUT
  char *filepath, int eps_num, int dim, int loops,
  char **eps_str, char ***kin, char **symbols, int ninvs
) {
  int print = 0;
  int *nadd = new int[dim];
  int **nfac; // = new int[dim];
  malloc_rk2_tens(nfac, dim, 2);
  int ***loop, ***id, ****intargs;
  malloc_rk3_tens(loop, dim, 2, loops);
  malloc_rk3_tens(id, dim, 2, loops);
  malloc_rk3_tens(intargs, dim, 2, loops);
  char ****sym, ****point;
  malloc_rk3_tens(sym, dim, 2, loops);
  malloc_rk3_tens(point, dim, 2, loops);
  for (int m=0; m<dim; m++) {
    for (int a=0; a<2; a++) {
      for (int k=0; k<loops; k++) {
        sym[m][a][k] = (char*) malloc(10*sizeof(char));
        point[m][a][k] = (char*) malloc(10*sizeof(char));
      }
    }
  }

  int nparam = 5;
  char **param_symbols = new char*[nparam];
  for (int s=0; s<nparam; s++) {
    param_symbols[s] = (char*) malloc(100*sizeof(char));
  }
  strcpy(param_symbols[0], "loop");
  strcpy(param_symbols[1], "id");
  strcpy(param_symbols[2], "sym");
  strcpy(param_symbols[3], "point");
  strcpy(param_symbols[4], "intargs");

  char ***factors;
  malloc_rk2_tens(factors, loops, nparam);

  nlist *BC = new nlist[dim];

  FILE *bound_fptr = fopen(filepath, "r");

  char MI_key[13];
  int param_idx;
  for (int m=0; m<dim; m++) {
    if (print) cout << "m = " << m << endl;
    snprintf(MI_key, sizeof(MI_key), "MI%d", m);
    // read_bound_factors(
    //   &factors, &nfac[m],
    //   loops, nparam, &param_symbols, MI_key, bound_fptr
    // );

    nlist_read_file(&BC[m], MI_key, bound_fptr);
    if (print) cout << BC[m].str << endl;
    
    nadd[m] = BC[m].nitems;

    char **intargs_str = NULL;
    int nintargs;

    for (int a=0; a<BC[m].nitems; a++) {
      if (print) cout << "a = " << a << endl;
      nfac[m][a] = BC[m].item[a].nitems;
      for (int k=0; k<BC[m].item[a].nitems; k++) {
        if (print) cout << "k = " << k << ": ";
        for (int j=0; j<BC[m].item[a].item[k].nitems; j++) {
          if (print) cout << "j = " << j << endl;
          if (print) cout << BC[m].item[a].item[k].item[j].item[0].str << endl;
          param_idx = sym_to_idx(BC[m].item[a].item[k].item[j].item[0].str, param_symbols, nparam);
          switch(param_idx) {
            case 0:
              loop[m][a][k] = atoi(BC[m].item[a].item[k].item[j].item[1].str);
              break;
            case 1:
              id[m][a][k] = atoi(BC[m].item[a].item[k].item[j].item[1].str);
              break;
            case 2:
              strcpy(sym[m][a][k], BC[m].item[a].item[k].item[j].item[1].str);
              break;
            case 3:
              strcpy(point[m][a][k], BC[m].item[a].item[k].item[j].item[1].str);
              break;
            case 4:
              nintargs = BC[m].item[a].item[k].item[j].item[1].nitems;
              intargs[m][a][k] = new int[nintargs];
              for (int i=0; i<nintargs; i++) {
                intargs[m][a][k][i] = atoi(BC[m].item[a].item[k].item[j].item[1].item[i].str);
              }
            break;
          }
        }
        
        if (print) {
        cout << loop[m][a][k] << ", ";
        cout << id[m][a][k] << ", ";
        cout << sym[m][a][k] << ", ";
        cout << point[m][a][k] << ", ";
        cout << endl;
        }

      }
      if (print) cout << "nfac = " << nfac[m][a] << endl;
    }
    
  }

  for (int ep=0; ep<eps_num; ep++) {
    dispatcher(
      bound[ep],
      eps_str[ep], kin, symbols, ninvs,
      dim, nadd, nfac, loop, id, sym, point, intargs
    );
  }

  // FREE
  nlist_rk1_free(BC, dim);
  delete[] BC;
}

