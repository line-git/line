#include <stdio.h>
#include <stdlib.h>
#include "mpc.h"
#include "global_vars.h"


//////
// EXPONENT BASED
//////
void exp_rel_err_mpfr_add(
  mpfr_t out, mpfr_t in1, mpfr_t in2, mpfr_rnd_t round, int wp_bin
) {
  int print = 0;
  if (print) {printf("in1 = "); mpfr_out_str(stdout, 10, 0, in1, MPFR_RNDN); printf("\n");}
  if (print) {printf("in2 = "); mpfr_out_str(stdout, 10, 0, in2, MPFR_RNDN); printf("\n");}
  // check for zero in the input
  if (mpfr_zero_p(in1)) {
    if (mpfr_zero_p(in2)) {
      mpfr_set_ui(out, 0, MPFR_RNDN);
      return;
    } else {
      mpfr_set(out, in2, round);
      return;
    }
  } else if (mpfr_zero_p(in2)) {
    if (mpfr_zero_p(in1)) {
      mpfr_set_ui(out, 0, MPFR_RNDN);
      return;
    } else {
      mpfr_set(out, in1, round);
      return;
    }
  }

  // check sign
  if (mpfr_sgn(in1) == mpfr_sgn(in2)) {
    mpfr_add(out, in1, in2, round);
    if (print) {printf("out = "); mpfr_out_str(stdout, 10, 0, out, MPFR_RNDN); printf("\n");}
    return;
  }

  mp_exp_t exp1 = mpfr_get_exp(in1);
  mp_exp_t exp2 = mpfr_get_exp(in2);
  if (print) printf("exp1 = %ld\n", exp1);
  if (print) printf("exp2 = %ld\n", exp2);
  if (exp1 > exp2) {
    if (exp1 - exp2 > 1) {
      mpfr_add(out, in1, in2, round);
      if (print) {printf("out = "); mpfr_out_str(stdout, 10, 0, out, MPFR_RNDN); printf("\n");}
      return;
    }
  } else {
    if (exp2 - exp1 > 1) {
      mpfr_add(out, in1, in2, round);
      if (print) {printf("out = "); mpfr_out_str(stdout, 10, 0, out, MPFR_RNDN); printf("\n");}
      return;
    } else {
      exp1 = exp2;  // store maximum in exp1
    }
  }

  mpfr_add(out, in1, in2, round);
  if (print) {printf("sum = "); mpfr_out_str(stdout, 10, 0, out, MPFR_RNDN); printf("\n");}
  exp2 = mpfr_get_exp(out);  // store result exponent in exp2
  if (print) printf("exp max = %ld\n", exp1);
  if (print) printf("exp sum = %ld\n", exp2);

  // check for cancellation
  if (exp1-exp2 > wp_bin) {
    if (print) printf("exp_rel_err_mpfr_add: set to zero\n");
    mpfr_set_ui(out, 0, round);
    if (print) {printf("out = "); mpfr_out_str(stdout, 10, 0, out, MPFR_RNDN); printf("\n");}
  }
  
  // if (exp2 < -140*10/3 && !mpfr_zero_p(out)) {
  //   printf("exp max = %ld\n", exp1);
  //   printf("exp sum = %ld\n", exp2);
  //   fflush(stdout);
  //   perror("exp_rel_err_mpc_add: cancellation missed");
  //   exit(1);
  // }

}


void exp_rel_err_mpc_add(
  mpc_t out, mpc_t in1, mpc_t in2, mpfr_rnd_t round, int wp_bin
) {
  exp_rel_err_mpfr_add(mpc_realref(out), mpc_realref(in1), mpc_realref(in2), round, wp_bin);
  exp_rel_err_mpfr_add(mpc_imagref(out), mpc_imagref(in1), mpc_imagref(in2), round, wp_bin);
}


void exp_rel_err_mpc_mul(
  mpc_t out, mpc_t in1, mpc_t in2, mpfr_rnd_t round, int wp_bin
) {
  int print = 0;
  if (print) {
    printf("rel_err_mpc_mul:\n");
    printf("in1 = "); mpc_out_str(stdout, 10, 0, in1, MPFR_RNDN); printf("\n");
    printf("in1 = "); mpc_out_str(stdout, 10, 0, in2, MPFR_RNDN); printf("\n");
  }
  

  mpfr_ptr re1 = mpc_realref(in1);
  mpfr_ptr im1 = mpc_imagref(in1);
  mpfr_ptr re2 = mpc_realref(in2);
  mpfr_ptr im2 = mpc_imagref(in2);
  // the order is important in case of overwriting
  if (mpfr_zero_p(im1)) {
    // in1 is real
    if (mpfr_zero_p(im2)) {
      // in2 is real
      // output is real
      mpfr_mul(mpc_realref(out), re1, re2, round);
      mpfr_set_ui(mpc_imagref(out), 0, round);
      if (print) {
        printf("real real\n");
        printf("out = "); mpc_out_str(stdout, 10, 0, out, round); printf("\n");
      }
      return;
    } else if (mpfr_zero_p(re2)) {
      // in2 is imaginary
      // output is imaginary      
      mpfr_mul(mpc_imagref(out), re1, im2, round);
      mpfr_set_ui(mpc_realref(out), 0, round);
      if (print) {
        printf("real imag\n");
        printf("out = "); mpc_out_str(stdout, 10, 0, out, round); printf("\n");
      }
      return;
    } else {
      // in2 is complex
      mpfr_mul(mpc_imagref(out), re1, im2, MPFR_RNDN);
      mpfr_mul(mpc_realref(out), re1, re2, MPFR_RNDN);
      return;
    }
  } else if (mpfr_zero_p(re1)) {
    // in1 is imaginary
    if (mpfr_zero_p(im2)) {
      // in2 is real
      // output is imaginary
      mpfr_mul(mpc_imagref(out), im1, re2, round);
      mpfr_set_ui(mpc_realref(out), 0, round);
      if (print) {
        printf("imag real\n");
        printf("out = "); mpc_out_str(stdout, 10, 0, out, round); printf("\n");
      }
      return;
    } else if (mpfr_zero_p(re2)) {
      // in2 is imaginary
      // output is real
      mpfr_mul(mpc_realref(out), im1, im2, round);
      mpfr_neg(mpc_realref(out), mpc_realref(out), round);
      mpfr_set_ui(mpc_imagref(out), 0, round);
      if (print) {
        printf("imag imag\n");
        printf("out = "); mpc_out_str(stdout, 10, 0, out, round); printf("\n");
      }
      return;
    } else {
      // in2 is complex
      // a copy is needed to prevent overwriting, otherwsie call mpc_mul
      mpc_mul(out, in1, in2, MPFR_RNDN);
      // mpfr_mul(mpc_realref(out), im1, im2, MPFR_RNDN);
      // mpfr_neg(mpc_realref(out), mpc_realref(out), MPFR_RNDN);
      // mpfr_mul(mpc_imagref(out), im1, re2, MPFR_RNDN);
      return;
    }
  } else {
    // in1 is complex
    if (mpfr_zero_p(im2)) {
      // in2 is real
      mpfr_mul(mpc_imagref(out), im1, re2, MPFR_RNDN);
      mpfr_mul(mpc_realref(out), re1, re2, MPFR_RNDN);
      return;
    } else if (mpfr_zero_p(re2)) {
      // in2 is imaginary
      // a copy is needed to prevent overwriting, otherwsie call mpc_mul
      mpc_mul(out, in1, in2, MPFR_RNDN);
      // mpfr_mul(mpc_realref(out), im1, im2, MPFR_RNDN);
      // mpfr_neg(mpc_realref(out), mpc_realref(out), MPFR_RNDN);
      // mpfr_mul(mpc_imagref(out), re1, im2, MPFR_RNDN);
      return;      
    }
    // in2 is complex
  }

  // save exponents
  mp_exp_t exp_re1 = mpfr_get_exp(re1)+mpfr_get_exp(re2);
  mp_exp_t exp_re2 = mpfr_get_exp(im1)+mpfr_get_exp(im2);
  mp_exp_t exp_im1 = mpfr_get_exp(re1)+mpfr_get_exp(im2);
  mp_exp_t exp_im2 = mpfr_get_exp(im1)+mpfr_get_exp(re2);
  if (print || dbg) {
  printf("exp_rel_err_mul:\n");
  printf("exp_re1 = %ld\n", exp_re1);
  printf("exp_re2 = %ld\n", exp_re2);
  printf("exp_im1 = %ld\n", exp_im1);
  printf("exp_im2 = %ld\n", exp_im2);
  }

  // save signs
  int sign_re = mpfr_sgn(re1)*mpfr_sgn(re2) == mpfr_sgn(im1)*mpfr_sgn(im2) ? 1 : 0;
  int sign_im = mpfr_sgn(re1)*mpfr_sgn(im2) != mpfr_sgn(im1)*mpfr_sgn(re2) ? 1 : 0;
  if (print || dbg) {
  printf("sign_re = %d\n", sign_re);
  printf("sign_im = %d\n", sign_im);
  }

  mpc_mul(out, in1, in2, round);
  if (print || dbg) {
    printf("exp_prod_re = %ld\n", mpfr_get_exp(mpc_realref(out)));
    printf("exp_prod_im = %ld\n", mpfr_get_exp(mpc_imagref(out)));
  }

  //////
  // REAL PART
  //////
  // check sign
  if (sign_re) {
    // cancellation is possible
    if (exp_re1 > exp_re2) {
      if (exp_re1 == exp_re2 + 1) {
        if (exp_re1 - mpfr_get_exp(mpc_realref(out)) > wp_bin) {
          // cancellation
          mpfr_set_ui(mpc_realref(out), 0, round);
        }
      }
    } else {
      if (exp_re2 <= exp_re1 + 1) {
        if (exp_re2 - mpfr_get_exp(mpc_realref(out)) > wp_bin) {
          // cancellation
          mpfr_set_ui(mpc_realref(out), 0, round);
        }
      }
    }
  }

  //////
  // IMAGINARY PART
  //////
  // check sign
  if (sign_im) {
    // cancellation is possible
    if (exp_im1 > exp_im2) {
      if (exp_im1 == exp_im2 + 1) {
        if (exp_im1 - mpfr_get_exp(mpc_imagref(out)) > wp_bin) {
          // cancellation
          mpfr_set_ui(mpc_imagref(out), 0, round);
        }
      }
    } else {
      if (exp_im2 <= exp_im1 + 1) {
        if (exp_im2 - mpfr_get_exp(mpc_imagref(out)) > wp_bin) {
          // cancellation
          mpfr_set_ui(mpc_imagref(out), 0, round);
        }
      }
    }
  }

  if (print) {printf("out = "); mpc_out_str(stdout, 10, 0, out, round); printf("\n");}

  // if (mpfr_get_exp(mpc_realref(out)) < -140*10/3  && !mpfr_zero_p(mpc_realref(out))) {
  //   printf("exp_re1 = %ld\n", exp_re1);
  //   printf("exp_re2 = %ld\n", exp_re2);
  //   printf("out = "); mpc_out_str(stdout, 10, 0, out, round); printf("\n");
  //   fflush(stdout);
  //   perror("exp_rel_err_mpc_mul: cancellation missed");
  //   exit(1);
  // }

  // if (mpfr_get_exp(mpc_imagref(out)) < -140*10/3  && !mpfr_zero_p(mpc_imagref(out))) {
  //   printf("exp_im1 = %ld\n", exp_im1);
  //   printf("exp_im2 = %ld\n", exp_im2);
  //   printf("out = "); mpc_out_str(stdout, 10, 0, out, round); printf("\n");
  //   fflush(stdout);
  //   perror("exp_rel_err_mpc_mul: cancellation missed");
  //   exit(1);
  // }
  
}


void exp_rel_err_mpc_div(
  mpc_t out, mpc_t in1, mpc_t in2, mpfr_rnd_t round, int wp_bin
) {
  int print = 0;
  if (print) {
    printf("rel_err_mpc_mul:\n");
    printf("in1 = "); mpc_out_str(stdout, 10, 0, in1, MPFR_RNDN); printf("\n");
    printf("in1 = "); mpc_out_str(stdout, 10, 0, in2, MPFR_RNDN); printf("\n");
  }
  

  mpfr_ptr re1 = mpc_realref(in1);
  mpfr_ptr im1 = mpc_imagref(in1);
  mpfr_ptr re2 = mpc_realref(in2);
  mpfr_ptr im2 = mpc_imagref(in2);
  if (mpfr_zero_p(im1)) {
    // in1 is real
    if (mpfr_zero_p(im2)) {
      // in2 is real
      // output is real
      mpfr_div(mpc_realref(out), re1, re2, round);
      mpfr_set_ui(mpc_imagref(out), 0, round);
      if (print) {
        printf("real real\n");
        printf("out = "); mpc_out_str(stdout, 10, 0, out, round); printf("\n");
      }
      return;
    } else if (mpfr_zero_p(re2)) {
      // in2 is imaginary
      // output is imaginary      
      mpfr_div(mpc_imagref(out), re1, im2, round);
      mpfr_neg(mpc_imagref(out), mpc_imagref(out), MPFR_RNDN);
      mpfr_set_ui(mpc_realref(out), 0, round);
      if (print) {
        printf("real imag\n");
        printf("out = "); mpc_out_str(stdout, 10, 0, out, round); printf("\n");
      }
      return;
    } else {
      // in2 is complex
      mpc_div(out, in1, in2, MPFR_RNDN);
      return;
    }
  } else if (mpfr_zero_p(re1)) {
    // in1 is imaginary
    if (mpfr_zero_p(im2)) {
      // in2 is real
      // output is imaginary
      mpfr_div(mpc_imagref(out), im1, re2, round);
      mpfr_set_ui(mpc_realref(out), 0, round);
      if (print) {
        printf("imag real\n");
        printf("out = "); mpc_out_str(stdout, 10, 0, out, round); printf("\n");
      }
      return;
    } else if (mpfr_zero_p(re2)) {
      // in2 is imaginary
      // output is real
      mpfr_div(mpc_realref(out), im1, im2, round);
      mpfr_set_ui(mpc_imagref(out), 0, round);
      if (print) {
        printf("imag imag\n");
        printf("out = "); mpc_out_str(stdout, 10, 0, out, round); printf("\n");
      }
      return;
    } else {
      // in2 is complex
      mpc_div(out, in1, in2, MPFR_RNDN);
      return;      
    }
  } else {
      // in1 is complex
    if (mpfr_zero_p(im2)) {
      // in2 is real
      mpfr_div(mpc_imagref(out), im1, re2, MPFR_RNDN);
      mpfr_div(mpc_realref(out), re1, re2, MPFR_RNDN);
      return;
    } else if (mpfr_zero_p(re2)) {
      // in2 is imaginary
      // a copy is needed to prevent overwriting, otherwsie call mpc_mul
      mpc_div(out, in1, in2, MPFR_RNDN);
      return;      
    }
    // in2 is complex
  }

  // save exponents
  mp_exp_t exp_mod2 = mpfr_get_exp(re2) > mpfr_get_exp(im2) ? mpfr_get_exp(re2) : mpfr_get_exp(im2);
  mp_exp_t exp_re1 = mpfr_get_exp(re1)+mpfr_get_exp(re2) - exp_mod2;
  mp_exp_t exp_re2 = mpfr_get_exp(im1)+mpfr_get_exp(im2) - exp_mod2;
  mp_exp_t exp_im1 = mpfr_get_exp(re1)+mpfr_get_exp(im2) - exp_mod2;
  mp_exp_t exp_im2 = mpfr_get_exp(im1)+mpfr_get_exp(re2) - exp_mod2;
  
  if (print) {
  printf("exp_rel_err_mul:\n");
  printf("exp_re1 = %ld\n", exp_re1);
  printf("exp_re2 = %ld\n", exp_re2);
  printf("exp_im1 = %ld\n", exp_im1);
  printf("exp_im2 = %ld\n", exp_im2);
  }

  // save signs
  int sign_re = mpfr_sgn(re1)*mpfr_sgn(re2) != mpfr_sgn(im1)*mpfr_sgn(im2) ? 1 : 0;
  int sign_im = mpfr_sgn(re1)*mpfr_sgn(im2) == mpfr_sgn(im1)*mpfr_sgn(re2) ? 1 : 0;
  if (print) {
  printf("sign_re = %d\n", sign_re);
  printf("sign_im = %d\n", sign_im);
  }

  mpc_div(out, in1, in2, round);

  //////
  // REAL PART
  //////
  // check sign
  if (sign_re) {
    // cancellation is possible
    if (exp_re1 > exp_re2) {
      if (exp_re1 == exp_re2 + 1) {
        if (exp_re1 - mpfr_get_exp(mpc_realref(out)) > wp_bin) {
          // cancellation
          mpfr_set_ui(mpc_realref(out), 0, round);
        }
      }
    } else {
      if (exp_re2 <= exp_re1 + 1) {
        if (exp_re2 - mpfr_get_exp(mpc_realref(out)) > wp_bin) {
          // cancellation
          mpfr_set_ui(mpc_realref(out), 0, round);
        }
      }
    }
  }

  //////
  // IMAGINARY PART
  //////
  // check sign
  if (sign_im) {
    // cancellation is possible
    if (exp_im1 > exp_im2) {
      if (exp_im1 == exp_im2 + 1) {
        if (exp_im1 - mpfr_get_exp(mpc_imagref(out)) > wp_bin) {
          // cancellation
          mpfr_set_ui(mpc_imagref(out), 0, round);
        }
      }
    } else {
      if (exp_im2 <= exp_im1 + 1) {
        if (exp_im2 - mpfr_get_exp(mpc_imagref(out)) > wp_bin) {
          // cancellation
          mpfr_set_ui(mpc_imagref(out), 0, round);
        }
      }
    }
  }

  if (print) {printf("out = "); mpc_out_str(stdout, 10, 0, out, round); printf("\n");}

  // if (mpfr_get_exp(mpc_realref(out)) < -140*10/3  && !mpfr_zero_p(mpc_realref(out))) {
  //   printf("exp_re1 = %ld\n", exp_re1);
  //   printf("exp_re2 = %ld\n", exp_re2);
  //   printf("out = "); mpc_out_str(stdout, 10, 0, out, round); printf("\n");
  //   fflush(stdout);
  //   perror("exp_rel_err_div_mul: cancellation missed");
  //   exit(1);
  // }

  // if (mpfr_get_exp(mpc_imagref(out)) < -140*10/3  && !mpfr_zero_p(mpc_imagref(out))) {
  //   printf("exp_im1 = %ld\n", exp_im1);
  //   printf("exp_im2 = %ld\n", exp_im2);
  //   printf("out = "); mpc_out_str(stdout, 10, 0, out, round); printf("\n");
  //   fflush(stdout);
  //   perror("exp_rel_err_mpc_div: cancellation missed");
  //   exit(1);
  // }
  
}
//
//////


int mpfr_sign_within_tol(mpfr_t in) {
  if (mpfr_cmpabs(in, mpfr_tol) < 0) {
    return 0;
  } else {
    return mpfr_cmp_ui(in, 0);
  }
}


void rel_err_mpfr_add(
  mpfr_t out, mpfr_t in1, mpfr_t in2, mpfr_rnd_t round 
) {
  // check signs
  if (mpfr_sgn(in1) == mpfr_sgn(in2)) {
    mpfr_add(out, in1, in2, round);
  } else {
    mpfr_t sum;
    mpfr_init2(sum, wp2);
    mpfr_t tmpfr;
    mpfr_init2(tmpfr, wp2);
    mpfr_add(sum, in1, in2, round);
    // printf("sum = "); mpfr_out_str(stdout, 10, 0, sum, MPFR_RNDN); printf("\n");

    if (mpfr_sign_within_tol(in1) != 0) {
      mpfr_div(tmpfr, sum, in1, round);
    } else if (mpfr_sign_within_tol(in2) != 0) {
      mpfr_div(tmpfr, sum, in2, round);
    } else {
      mpfr_set(tmpfr, sum, round);
    }
    // printf("tmpfr = "); mpfr_out_str(stdout, 10, 0, tmpfr, MPFR_RNDN); printf("\n");
    // printf("mpfr_tol = "); mpfr_out_str(stdout, 10, 0, mpfr_tol, MPFR_RNDN); printf("\n");

    if (mpfr_sign_within_tol(tmpfr) == 0) {
      mpfr_set_ui(out, 0, round);
    } else {
      mpfr_set(out, sum, round);
    }

    // FREE
    mpfr_clear(sum);
    mpfr_clear(tmpfr);
  }

}


void rel_err_mpfr_sub(
  mpfr_t out, mpfr_t in1, mpfr_t in2, mpfr_rnd_t round 
) {
  // check signs
  if (mpfr_sgn(in1) != mpfr_sgn(in2)) {
    mpfr_add(out, in1, in2, round);
  } else {
    mpfr_t sum;
    mpfr_init2(sum, wp2);
    mpfr_t tmpfr;
    mpfr_init2(tmpfr, wp2);
    mpfr_sub(sum, in1, in2, round);
    // cout << "sum = "; print_mpfr(&sum); cout << endl;

    if (mpfr_sign_within_tol(in1) != 0) {
      mpfr_div(tmpfr, sum, in1, round);
    } else if (mpfr_sign_within_tol(in2) != 0) {
      mpfr_div(tmpfr, sum, in2, round);
    } else {
      mpfr_set(tmpfr, sum, round);
    }
    // cout << "tmpfr = "; print_mpfr(&tmpfr); cout << endl;

    if (mpfr_sign_within_tol(tmpfr) == 0) {
      mpfr_set_ui(out, 0, round);
    } else {
      mpfr_set(out, sum, round);
    }
  }

}


void rel_err_mpc_add(
  mpc_t out, mpc_t in1, mpc_t in2, mpfr_rnd_t round 
) {
  rel_err_mpfr_add(mpc_realref(out), mpc_realref(in1), mpc_realref(in2), round);
  rel_err_mpfr_add(mpc_imagref(out), mpc_imagref(in1), mpc_imagref(in2), round);
}


void rel_err_mpc_sub(
  mpc_t out, mpc_t in1, mpc_t in2, mpfr_rnd_t round 
) {
  rel_err_mpfr_sub(mpc_realref(out), mpc_realref(in1), mpc_realref(in2), round);
  rel_err_mpfr_sub(mpc_imagref(out), mpc_imagref(in1), mpc_imagref(in2), round);
}


void rel_err_mpfr_fma(
  mpfr_t out, mpfr_t z0, mpfr_t z1, mpfr_t z2, mpfr_rnd_t round
) {
  __mpfr_struct *in0, *in1, *in2;
  if (out == z0) {
    in0 = (__mpfr_struct*) malloc(sizeof(__mpfr_struct));
    mpfr_init2(in0, wp2);
    mpfr_set(in0, z0, MPFR_RNDN);
    if (out == z1) {
      in1 = (__mpfr_struct*) malloc(sizeof(__mpfr_struct));
      mpfr_init2(in1, wp2);
      mpfr_set(in1, z1, MPFR_RNDN);
      if (out == z2) {
        in2 = (__mpfr_struct*) malloc(sizeof(__mpfr_struct));
        mpfr_init2(in2, wp2);
        mpfr_set(in2, z2, MPFR_RNDN);
      } else {
        in2 = z2;
      }
    } else {
      in1 = z1;
    }
    if (out == z2) {
      in2 = (__mpfr_struct*) malloc(sizeof(__mpfr_struct));
      mpfr_init2(in2, wp2);
      mpfr_set(in2, z2, MPFR_RNDN);
    } else {
      in2 = z2;
    }
  } else if (out == z1) {
    in1 = (__mpfr_struct*) malloc(sizeof(__mpfr_struct));
    mpfr_init2(in1, wp2);
    mpfr_set(in1, z1, MPFR_RNDN);
    if (out == z0) {
      in0 = (__mpfr_struct*) malloc(sizeof(__mpfr_struct));
      mpfr_init2(in0, wp2);
      mpfr_set(in0, z0, MPFR_RNDN);
      if (out == z2) {
        in2 = (__mpfr_struct*) malloc(sizeof(__mpfr_struct));
        mpfr_init2(in2, wp2);
        mpfr_set(in2, z2, MPFR_RNDN);
      } else {
        in2 = z2;
      }
    } else {
      in0 = z0;
    }
    if (out == z2) {
      in2 = (__mpfr_struct*) malloc(sizeof(__mpfr_struct));
      mpfr_init2(in2, wp2);
      mpfr_set(in2, z2, MPFR_RNDN);
    } else {
      in2 = z2;
    }
  } else if (out == z2) {
    in2 = (__mpfr_struct*) malloc(sizeof(__mpfr_struct));
    mpfr_init2(in2, wp2);
    mpfr_set(in2, z2, MPFR_RNDN);
    if (out == z1) {
      in1 = (__mpfr_struct*) malloc(sizeof(__mpfr_struct));
      mpfr_init2(in1, wp2);
      mpfr_set(in1, z1, MPFR_RNDN);
      if (out == z0) {
        in0 = (__mpfr_struct*) malloc(sizeof(__mpfr_struct));
        mpfr_init2(in0, wp2);
        mpfr_set(in0, z0, MPFR_RNDN);
      } else {
        in0 = z0;
      }
    } else {
      in1 = z1;
    }
    if (out == z0) {
      in0 = (__mpfr_struct*) malloc(sizeof(__mpfr_struct));
      mpfr_init2(in0, wp2);
      mpfr_set(in0, z0, MPFR_RNDN);
    } else {
      in0 = z0;
    }
  } else {
    in0 = z0;
    in1 = z1;
    in2 = z2;
  }

  //
  // in0*in1 + in2
  //
  mpfr_t tmpfr;
  mpfr_init2(tmpfr, wp2);

  if (mpfr_sign_within_tol(in2) != 0) {
    mpfr_fma(out, in0, in1, in2, round);
    mpfr_div(tmpfr, out, in2, round);
  } else {
    // out = in0*in1
    mpfr_mul(out, in0, in1, round);
    return;
  }
  
  if (mpfr_sign_within_tol(tmpfr) == 0) {
    // out = 0
    mpfr_set_ui(out, 0, round);
    return;
  } else {
    return;
  }

}


void rel_err_mpfr_fmma(
  mpfr_t out, mpfr_t in0, mpfr_t in1, mpfr_t in2, mpfr_t in3, mpfr_rnd_t round
) {
  //
  // in0*in1 + in2*in3
  //
  mpfr_t tmpfr;
  mpfr_init2(tmpfr, wp2);

  if (mpfr_sign_within_tol(in0) != 0 && mpfr_sign_within_tol(in1) != 0) {
    mpfr_fmma(out, in0, in1, in2, in3, round);
    mpfr_div(tmpfr, out, in0, round);
    mpfr_div(tmpfr, out, in1, round);
  } else {
    // out = in2*in3
    mpfr_mul(out, in2, in3, round);
    return;
  }

  if (mpfr_sign_within_tol(tmpfr) == 0) {
    // out = 0
    mpfr_set_ui(out, 0, round);
    return;
  } else {
    return;
  }

}


void rel_err_mpfr_dflop(
  mpfr_t out,
  mpfr_t in0, mpfr_t in1, mpfr_t in2, mpfr_t in3, mpfr_t in4,
  mpfr_rnd_t round
) {
  //
  // in0*in1 + in2*in3 + i4
  //
  int sign, sign1, sign2, sign3;
  sign1 = mpfr_sgn(in0)*mpfr_sgn(in1);  // in0*in1
  sign2 = mpfr_sgn(in2)*mpfr_sgn(in3);  // in2*in3
  sign3 = mpfr_sgn(in4);  // in4
  sign = sign1 + sign2 + sign3;
  // cout << "sign1 = " << sign1 << endl;
  // cout << "sign2 = " << sign2 << endl;
  // cout << "sign3 = " << sign3 << endl;
  // cout << "sign = " << sign << endl;
  if (sign == 1 || sign == -1) {
    mpfr_t tmpfr;
    mpfr_init2(tmpfr, wp2);
    if (sign1 == sign2) {
      if (sign1 == 0) {
        // out = in4
        mpfr_set(out, in4, MPFR_RNDN);
        return;
      }

      // in0*in1 same sign as in2*in3
      int exp1 = mpfr_get_exp(in0) + mpfr_get_exp(in1);  // in0*in1
      int exp2 = mpfr_get_exp(in2) + mpfr_get_exp(in3);  // in2*in3
      if (exp1 >= exp2) {
        rel_err_mpfr_fma(tmpfr, in0, in1, in4, round);  // in0*in1 + in4
        rel_err_mpfr_fma(out, in2, in3, tmpfr, round);  // (in0*in1 + in4) + in2*in3
      } else {
        rel_err_mpfr_fma(tmpfr, in2, in3, in4, round);  // in2*in3 + in4
        rel_err_mpfr_fma(out, in0, in1, tmpfr, round);  // (in2*in3 + in4) + in0*in1
      }

    } else if (sign1 == sign3) {
      // cout << "sign1 == sing3" << endl;
      if (sign1 == 0) {
        // out = in2*in3
        mpfr_mul(out, in2, in3, MPFR_RNDN);
        return;
      }

      // in0*in1 same sign as in4
      int exp1 = mpfr_get_exp(in0) + mpfr_get_exp(in1);  // in0*in1
      int exp3 = mpfr_get_exp(in4);  // in4
      if (exp1 >= exp3) {        
        rel_err_mpfr_fmma(out, in0, in1, in2, in3, round);  // in0*in1 + in2*in3
        rel_err_mpfr_add(out, out, in4, round);
      } else {
        rel_err_mpfr_fma(tmpfr, in2, in3, in4, round);  // in4 + in2*in3
        rel_err_mpfr_fma(out, in0, in1, tmpfr, round);  // (in4 + in2*in3) + in0*in1
      }

    } else if (sign2 == sign3) {
      // cout << "sign2 == sing3" << endl;
      if (sign2 == 0) {
        // out = in0*in1
        mpfr_mul(out, in0, in1, MPFR_RNDN);
        return;
      }

      // in2*in3 same sign as in4
      int exp2 = mpfr_get_exp(in2) + mpfr_get_exp(in3);  // in2*in3
      int exp3 = mpfr_get_exp(in4);  // in4
      if (exp2 >= exp3) {        
        rel_err_mpfr_fmma(out, in2, in3, in0, in1, round);  // in2*in3 + in0*in1
        rel_err_mpfr_add(out, out, in4, round);
      } else {
        rel_err_mpfr_fma(tmpfr, in0, in1, in4, round);  // in4 + in0*in1
        rel_err_mpfr_fma(out, in2, in3, tmpfr, round);  // (in4 + in0*in1) + in2*in3
      }

    }
  } else if (sign == 0) {
    mpfr_t tmpfr;
    mpfr_init2(tmpfr, wp2);
    if (sign1 == 0) {
      if (sign2 == 0) {
        // out = 0
        mpfr_set_ui(out, 0, round);
        return;
      }

      // out = in2*in3 + in4
      rel_err_mpfr_fma(out, in2, in3, in4, round);

    } else if (sign2 == 0) {
      if (sign1 == 0) {
        // out = 0
        mpfr_set_ui(out, 0, round);
        return;
      }

      // in0*in1 + in4
      rel_err_mpfr_fma(out, in0, in1, in4, round);

    } else if (sign3 == 0) {
      if (sign1 == 0) {
        // out = 0
        mpfr_set_ui(out, 0, round);
        return;
      }

      // in0*in1 + in2*in3
      rel_err_mpfr_fmma(out, in0, in1, in2, in3, round);

    }

  } else {
    // cout << "no check needed" << endl;
    mpfr_fmms(
      out,
      in0, in1, in2, in3,
      round
    );
    mpfr_add(out, out, in4, round);

  }

}


void rel_err_mpc_fma(
  mpc_t out, mpc_t z0, mpc_t z1, mpc_t z2, mpfr_rnd_t round
) {
  __mpc_struct *in0, *in1, *in2;
  int overwrite;
  if (out == z0) {
    in0 = (__mpc_struct*) malloc(sizeof(__mpc_struct));
    mpc_init3(in0, wp2, wp2);
    mpc_set(in0, z0, MPFR_RNDN);
    if (out == z1) {
      in1 = (__mpc_struct*) malloc(sizeof(__mpc_struct));
      mpc_init3(in1, wp2, wp2);
      mpc_set(in1, z1, MPFR_RNDN);
      if (out == z2) {
        in2 = (__mpc_struct*) malloc(sizeof(__mpc_struct));
        mpc_init3(in2, wp2, wp2);
        mpc_set(in2, z2, MPFR_RNDN);
        overwrite = 7;  // 111
      } else {
        in2 = z2;
        overwrite = 6;  // 110
      }
    } else {
      in1 = z1;
    }
    if (out == z2) {
      in2 = (__mpc_struct*) malloc(sizeof(__mpc_struct));
      mpc_init3(in2, wp2, wp2);
      mpc_set(in2, z2, MPFR_RNDN);
      overwrite = 5; // 101
    } else {
      in2 = z2;
      overwrite = 4; // 100
    }
  } else if (out == z1) {
    in1 = (__mpc_struct*) malloc(sizeof(__mpc_struct));
    mpc_init3(in1, wp2, wp2);
    mpc_set(in1, z1, MPFR_RNDN);
    if (out == z0) {
      in0 = (__mpc_struct*) malloc(sizeof(__mpc_struct));
      mpc_init3(in0, wp2, wp2);
      mpc_set(in0, z0, MPFR_RNDN);
      if (out == z2) {
        in2 = (__mpc_struct*) malloc(sizeof(__mpc_struct));
        mpc_init3(in2, wp2, wp2);
        mpc_set(in2, z2, MPFR_RNDN);
        overwrite = 7;  // 111
      } else {
        in2 = z2;
        overwrite = 6; // 110
      }
    } else {
      in0 = z0;
    }
    if (out == z2) {
      in2 = (__mpc_struct*) malloc(sizeof(__mpc_struct));
      mpc_init3(in2, wp2, wp2);
      mpc_set(in2, z2, MPFR_RNDN);
      overwrite = 3;  // 011
    } else {
      in2 = z2;
      overwrite = 2;  // 010
    }
  } else if (out == z2) {
    in2 = (__mpc_struct*) malloc(sizeof(__mpc_struct));
    mpc_init3(in2, wp2, wp2);
    mpc_set(in2, z2, MPFR_RNDN);
    if (out == z1) {
      in1 = (__mpc_struct*) malloc(sizeof(__mpc_struct));
      mpc_init3(in1, wp2, wp2);
      mpc_set(in1, z1, MPFR_RNDN);
      if (out == z0) {
        in0 = (__mpc_struct*) malloc(sizeof(__mpc_struct));
        mpc_init3(in0, wp2, wp2);
        mpc_set(in0, z0, MPFR_RNDN);
        overwrite = 7;  // 111
      } else {
        in0 = z0;
        overwrite = 3;  // 011
      }
    } else {
      in1 = z1;
    }
    if (out == z0) {
      in0 = (__mpc_struct*) malloc(sizeof(__mpc_struct));
      mpc_init3(in0, wp2, wp2);
      mpc_set(in0, z0, MPFR_RNDN);
      overwrite = 5;  // 101
    } else {
      in0 = z0;
      overwrite = 1;  // 001
    }
    // cout << "z0 = "; mpc_out_str(stdout, 10, 0, z0, MPFR_RNDN); cout << endl;
    // cout << "in0 = "; mpc_out_str(stdout, 10, 0, in0, MPFR_RNDN); cout << endl;
    // cout << "z1 = "; mpc_out_str(stdout, 10, 0, z1, MPFR_RNDN); cout << endl;
    // cout << "in1 = "; mpc_out_str(stdout, 10, 0, in1, MPFR_RNDN); cout << endl;
    // cout << "z2 = "; mpc_out_str(stdout, 10, 0, z2, MPFR_RNDN); cout << endl;
    // cout << "in2 = "; mpc_out_str(stdout, 10, 0, in2, MPFR_RNDN); cout << endl;
  } else {
    in0 = z0;
    in1 = z1;
    in2 = z2;
    overwrite = 0;  // 000
  }

  // REAL PART
  //
  // x_out = x0*x1 - y0*y1 + x2
  //
  mpfr_neg(mpc_imagref(in0), mpc_imagref(in0), round);
  rel_err_mpfr_dflop(
    mpc_realref(out),
    mpc_realref(in0), mpc_realref(in1), mpc_imagref(in0), mpc_imagref(in1), mpc_realref(in2),
    round
  );
  mpfr_neg(mpc_imagref(in0), mpc_imagref(in0), round);

  // IMAGINARY PART
  //
  // y_out = x0*y1 + y0*x1 + y2
  //
  rel_err_mpfr_dflop(
    mpc_imagref(out),
    mpc_realref(in0), mpc_imagref(in1), mpc_imagref(in0), mpc_realref(in1), mpc_imagref(in2),
    round
  );

  switch (overwrite) {
    case 0:  // 000
    break;
    case 1:  // 001
    mpc_clear(in2); free(in2);
    break;
    case 2:  // 010
    mpc_clear(in1); free(in1);
    break;
    case 3:  // 011
    mpc_clear(in1); free(in1);
    mpc_clear(in2); free(in2);
    break;
    case 4:  // 100
    mpc_clear(in0); free(in0);
    break;
    case 5:  // 101
    mpc_clear(in0); free(in0);
    mpc_clear(in2); free(in2);
    break;
    case 6:  // 110
    mpc_clear(in0); free(in0);
    mpc_clear(in1); free(in1);
    break;
    case 7:  // 111
    mpc_clear(in0); free(in0);
    mpc_clear(in1); free(in1);
    mpc_clear(in2); free(in2);
    break;
  }

  
}

