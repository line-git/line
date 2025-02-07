#include <iostream>
// #include <complex>
#include "mpc.h"

using namespace std;

#include "global_vars.h"
#include "malloc_defs.h"
#include "utils.h"
#include "setup.h"
#include "tensor_utils.h"
#include "conversions.h"
#include "system_analyzer.h"
#include "poly_frac.h"
#include "pf_normalize.h"
// #include "bc_gen.h"
extern "C" {
  #include "algebra.h"
  #include "cpoly.h"
  // #include "in_out.h"
  #include "jordan.h"
  #include "rel_err_mpc.h"
  #include "in_out.h"
}

// needed for csvd
static int c__0 = 0;


int factorial(int n) {
  int fact;
  if (n == 1 || n == 0) {
    return 1;
  } else if (n < 0) {
    perror("trying to compute the factorial of a negative number.");
    exit(1);
  } else {
    return fact = n*factorial(n-1);
  }
}


void mpc_to_math(
  string *output,
  mpc_t input
) {
  char *mantissa_real = NULL, *mantissa_imag = NULL;
  string sign_real = "", sign_imag = "", id_mat = "", mpc_str;
  long exponent_real, exponent_imag;
  int bor, boi;

  if (mpfr_nan_p(mpc_realref(input)) && mpfr_nan_p(mpc_imagref(input))) {
    *output = "0";
    return;
  }

  mantissa_real = mpfr_get_str(NULL, &exponent_real, 10, 0, mpc_realref(input), MPFR_RNDN);
  mantissa_imag = mpfr_get_str(NULL, &exponent_imag, 10, 0, mpc_imagref(input), MPFR_RNDN);

  sign_real = "";
  bor = 0;
  if (mpfr_sgn(mpc_realref(input)) <= 0) {
    sign_real = "-";
    bor = 1;
  }
  sign_imag = "";
  boi = 0;
  if (mpfr_sgn(mpc_imagref(input)) <= 0) {
    sign_imag = "-";
    boi = 1;
  }
  *output = sign_real+"0."+((string) (mantissa_real+bor))+"*^"+to_string(exponent_real)+\
    " + I*("+sign_imag+"0."+((string) (mantissa_imag+boi))+"*^"+to_string(exponent_imag)+")";
  
  mpfr_free_str(mantissa_real);
  mpfr_free_str(mantissa_imag);
}


void mpc_rk2_to_math(
  string *output,
  mpc_t **input, int dim1, int dim2
) {
  string tmp_str;
  *output = "{\n";
  for (int i1=0; i1<dim1; i1++) {
    *output += "{\n";
    for (int i2=0; i2<dim2; i2++) {
      mpc_to_math(&tmp_str, input[i1][i2]);
      *output += tmp_str;
      if (i2 < dim2-1) {
        *output += ",";
      }
      *output += "\n";
    }
    *output +="}";
    if (i1 < dim1-1) {
      *output += ",";
    }
    *output += "\n";
  }
  *output +="};\n";
}


void solve_block(mpc_t **solutions,
  mpc_t *lcm, mpc_t ***block, mpc_t **const_term,
  int *b_idx, int b_len, int sb_len,
  int eta_ord, int lcm_deg, int block_max_deg
) {
  //////
  // CHECK IF INPUT BLOCK IS ZERO
  /////
  // #2BD: check also sub-diagonal block (se non e' zero allora c'è il termine
  // costante e quindi la derivata non fa zero, quindi c'e' il termine ordine 1)
  // (probabilmente non e' possibile perche' i MI sono tali perche' irriducibili
  // a MI di sottotopologie)
  int found_non_zero = 0;
  mpfr_t mpfr_zero;
  mpfr_init2(mpfr_zero, wp2);
  mpfr_set_ui(mpfr_zero, 0, MPFR_RNDN);
  for (int i=0; i<b_len; i++) {
    for (int j=0; j<b_len; j++) {
      for (int k=0; k<=block_max_deg; k++) {
        if(!mpc_lessthan_tol(block[i][j][k])) {
        // if(mpfr_equal_p(mpc_realref(block[i][j][k]), mpfr_zero) && mpfr_equal_p(mpc_imagref(block[i][j][k]), mpfr_zero)) {
          // cout << "non-zero: i, j = " << i << ", " << j << endl;
          found_non_zero = 1;
          break;
        }
      }
      if (found_non_zero) {
        break;
      }
    }
    if (found_non_zero) {
      break;
    }
  }
  if (!found_non_zero) {
    for (int i=0; i<b_len; i++) {
      for (int k=1; k<=eta_ord; k++) {
        mpc_set_d(solutions[b_idx[i]][k], 0, MPFR_RNDN);
      }
    }
    // cout << "hit new return" << endl;
    return;
  }

  int max_j;
  mpc_t tmp_rec;
  mpc_init3(tmp_rec, wp2, wp2);
  mpc_t inv;
  mpc_init3(inv, wp2, wp2);
  for (int i=eta_ord; i>0; i--) {
    // cout << "i = " << i << endl;
    // iteration to construct solutions[:][eta_ord+1-i] in terms of
    // the lower order soulutions[:][eta_ord+1-j] with j>i

    //////
    // LHS MOVED TO RHS
    //////
    // initialize b_len-dim vector to zero
    for (int m=0; m<b_len; m++) {
      // mpc_init3(solutions[b_idx[m]][eta_ord+1-i], wp2, wp2);
      mpc_set_d(solutions[b_idx[m]][eta_ord+1-i], 0, MPFR_RNDN);
    }
    // cout << "initialized to zero solutions " << endl;
    // temporary variable needed for intermediate multiplications
    // cout << "initialized tmp_rec" << endl;

    // BLOCK CONTRIBUTION
    max_j = min(i+block_max_deg+1, eta_ord+1);
    // cout << "max_j = " << max_j << endl;
    for (int j=i+1; j<=max_j; j++) {
      // cout << "j = " << j << endl;
      // row by column product (block times solutions of previous iterations)
      // (b_len x b_len) times b_len
      // block[:][:][j-i+1] * solutions[:, eta_ord+1-j]
      for (int m=0; m<b_len; m++) {
        // fill m-th component
        for (int n=0; n<b_len; n++) {
          // cout << "n = " << n << endl;
          // print_mpc(&solutions[0][eta_ord+1-j]);
          // cout << endl;
          // print_mpc(&block[m][n][j-i-1]);
          // cout << endl;
          // print_mpc(&solutions[b_idx[n]][eta_ord+1-j]);
          // cout << endl;
          // cout << endl;
          mpc_mul(tmp_rec, block[m][n][j-i-1], solutions[b_idx[n]][eta_ord+1-j], MPFR_RNDN);
          // print_mpc(&tmp_rec);
          // cout << endl;
          // print_mpc(&solutions[b_idx[m]][eta_ord+1-i]);
          // cout << endl;
          mpc_add(solutions[b_idx[m]][eta_ord+1-i], solutions[b_idx[m]][eta_ord+1-i], tmp_rec, MPFR_RNDN);
          // print_mpc(&solutions[b_idx[m]][eta_ord+1-i]);
          // cout << endl;
        }
      }
    }
    // cout << "done block contribution" << endl;

    // LCM CONTRIBUTION
    max_j = min(i+lcm_deg, eta_ord+1);
    // cout << "j_max = " << max_j << endl;
    for (int j=i+1; j<=max_j; j++) {
      // cout << "j = " << j << endl;
      // cout << eta_ord+1-j << endl;
      // scalar times b_len-dim vector: (eta_ord+1-j) * lcm[j-i] * solutions[:][eta_ord+1-j]
      for (int m=0; m<b_len; m++) {
        mpc_mul(tmp_rec, lcm[j-i], solutions[b_idx[m]][eta_ord+1-j], MPFR_RNDN);
        mpc_mul_si(tmp_rec, tmp_rec, eta_ord+1-j, MPFR_RNDN);
        mpc_sub(solutions[b_idx[m]][eta_ord+1-i], solutions[b_idx[m]][eta_ord+1-i], tmp_rec, MPFR_RNDN);
      }
    }
    // cout << "done lcm contribution" << endl;
    // cout << endl;



    //////
    // CONSTANT TERM
    //////
    if (sb_len > 0) {
      // cout << "building RHS" << endl;
      // sum of b_len-dim vectors: solutions[:][eta_ord+1-i] + const_term[:][eta_ord-i]
      for (int m=0; m<b_len; m++) {
        // cout << "m = " << m << endl;
        mpc_add(solutions[b_idx[m]][eta_ord+1-i], solutions[b_idx[m]][eta_ord+1-i], const_term[m][eta_ord-i], MPFR_RNDN);
      }
    }
    // cout << "built RHS" << endl;

    //////
    // DIVISION BY SCALAR
    //////
    // cout << eta_ord+1-i << endl;
    mpc_mul_si(inv, lcm[0], eta_ord+1-i, MPFR_RNDN);
    // print_mpc(&inv);
    // cout << endl;
    // b_len-dim vector divided by scalar: solutions[:][eta_ord+1-j] / inv
    for (int m=0; m<b_len; m++) {
      // print_mpc(&solutions[b_idx[m]][eta_ord+1-i]);
      // cout << endl;
      mpc_div(solutions[b_idx[m]][eta_ord+1-i], solutions[b_idx[m]][eta_ord+1-i], inv, MPFR_RNDN);
      // print_mpc(&solutions[b_idx[m]][eta_ord+1-i]);
      // cout << endl;
    }
  }

  // FREE
  mpfr_clear(mpfr_zero);
  mpc_clear(tmp_rec);
  mpc_clear(inv);
}


void zero_solve_block(mpc_t *****solutions,
  mpc_t *lcm, mpc_t ***block,
  int *b_idx, int b_len, int sb_len,
  int eta_ord, int lcm_deg, int block_max_deg,
  mpc_t *eig_list, int *cum_eig, int num_cum_eig,
  int *block_log_len, int *sol_log_len, int num_classes,
  bool particular=false, mpc_t ****const_term = NULL, int *rhs_log_len = NULL
  // ,int match_in_zero = 0, int lam0 = 0, mpc_t **boundary = NULL
) {
  int print = 0;
  // lcm_deg++;
  // lcm--;
  int dim2 = b_len, k_start;
  if (particular) {
    dim2 = 1;
  }
  mpc_t *ptr, mpc_tmp;
  mpc_init3(mpc_tmp, wp2, wp2);
  mpc_t **mpc_mat;
  malloc_rk2_tens(mpc_mat, b_len, dim2);
  init_rk2_mpc(mpc_mat, b_len, dim2);
  mpc_t **mpc_mat1;
  malloc_rk2_tens(mpc_mat1, b_len, b_len);
  init_rk2_mpc(mpc_mat1, b_len, b_len);
  mpc_t **mpc_mat2;
  malloc_rk2_tens(mpc_mat2, b_len, b_len);
  init_rk2_mpc(mpc_mat2, b_len, b_len);
  int max_pow, min_pow, tmp_pow, tmp_pow2; 
  bool lcm_greater;
  if (lcm_deg+1 - 1 > block_max_deg) {
    max_pow = lcm_deg+1 - 1;
    min_pow = block_max_deg;
    lcm_greater = true;
  } else {
    min_pow = lcm_deg+1 - 1;
    max_pow = block_max_deg;
    lcm_greater = false;
  }
  if (print) {
  cout << "lcm_deg+1 = " << lcm_deg+1 << endl;
  cout << "block deg = "<< block_max_deg << endl;
  cout << "min pow = " << min_pow << endl;
  cout << "max pow = " << max_pow << endl;
  cout << "lcm_greter = " << lcm_greater << endl;
  }
  // getchar();

  //////
  // SET CONSTRAINTS
  //////
  //
  // - homogeneous: set zero-order coefficients to identity matrix
  // - particular: set zero-order coefficients to zero
  //
  //////
  if (!particular) {
    for (int lam, n=0; n<num_cum_eig; n++) {
      lam = cum_eig[n];
      if (block_log_len[lam] == 0) {
        continue;
      }
      // set to identity matrix
      for (int i=0; i<b_len; i++) {
        for (int j=0; j<dim2; j++) {
          if (i == j) {
            mpc_set_d(solutions[n][0][i][j][0], 1, MPFR_RNDN);
          } else {
            mpc_set_d(solutions[n][0][i][j][0], 0, MPFR_RNDN);
          }
        }
      }
    }
  // } else if (particular && match_in_zero == 0) {
  } else if (particular) {
    // cout << "entro in match_in_zero == 0" << endl;
    for (int lam, n=0; n<num_cum_eig; n++) {
      lam = cum_eig[n];
      if (sol_log_len[lam] == 0) {
        continue;
      }
      // set to zero
      for (int i=0; i<b_len; i++) {
        for (int j=0; j<dim2; j++) {
          mpc_set_d(solutions[n][0][i][j][0], 0, MPFR_RNDN);
        }
      }
    }
  // } else if (particular && match_in_zero == 1) {
  //   cout << "entro in match_in_zero == 1" << endl;
  //   for (int lam, n=0; n<num_cum_eig; n++) {
  //     lam = cum_eig[n];
  //     if (sol_log_len[lam] == 0) {
  //       continue;
  //     }
  //     cout << "lam = " << lam << endl;
  //     if (lam == lam0) {
  //       for (int i=0; i<b_len; i++) {
  //         for (int j=0; j<dim2; j++) {
  //           print_mpc(&boundary[i][0]); cout << endl;
  //           mpc_set(solutions[n][0][i][j][0], boundary[i][0], MPFR_RNDN);
  //         }
  //       }
  //     } else {
  //       // set to zero
  //       for (int i=0; i<b_len; i++) {
  //         for (int j=0; j<dim2; j++) {
  //           mpc_set_d(solutions[n][0][i][j][0], 0, MPFR_RNDN);
  //         }
  //       }
  //     }
  //   }
  }


  // #2BD: construct the matrix block - lam*G and remove checks about the diagonal,
  // since right now this operation is reapeted more than once (every l)

  //////
  // SOLVE BY RECURRENCE
  //////
  // cout << "solve by recurrence" << endl;
  int found_non_zero = 0;
  mpfr_t mpfr_zero;
  mpfr_init2(mpfr_zero, wp2);
  mpfr_set_ui(mpfr_zero, 0, MPFR_RNDN);
  for (int i=0; i<b_len; i++) {
    for (int j=0; j<b_len; j++) {
      for (int k=0; k<=block_max_deg; k++) {
        // if(mpfr_equal_p(mpc_realref(block[i][j][k]), mpfr_zero) && mpfr_equal_p(mpc_imagref(block[i][j][k]), mpfr_zero)) {
        if(!mpc_lessthan_tol(block[i][j][k])) {
          found_non_zero = 1;
          break;
        }
      }
      if (found_non_zero) {
        break;
      }
    }
    if (found_non_zero) {
      break;
    }
  }
  if (!found_non_zero && !particular) {
    // cout << "block is zero" << endl;
    // for (int i=0; i<b_len; i++) {
    //   for (int j=0; j<b_len; j++) {
    //     cout << "i, j = " << i << ", " << j << endl;
    //     print_poly(block[i][j], block_max_deg);
    //   }
    // }
    for (int lam, n=0; n<num_cum_eig; n++) {
      lam = cum_eig[n];
      if (block_log_len[lam] == 0) {
        continue;
      }

      for (int i=0; i<b_len; i++) {
        for (int j=0; j<dim2; j++) {
          mpc_set_d(solutions[n][sol_log_len[lam]-1][i][j][0], 0, MPFR_RNDN);
        }
      }
      k_start = 1;
      if (block_log_len[lam] == 0) {
        if (particular) {
          k_start = 0;
        }
      }
      // cout << "k_start = " << k_start << endl;
      for (int l=sol_log_len[lam]-2; l>=1; l--) {
        // cout << "l = " << l << endl;
        for (int k=0; k<=eta_ord; k++) {
          // cout << "k = " << k << endl;
          for (int i=0; i<b_len; i++) {
            // cout << "i = " << i << endl;
            for (int j=0; j<dim2; j++) {
              mpc_set_d(solutions[n][l][i][j][k], 0, MPFR_RNDN);
            }
          }
        }
      }

      for (int k=k_start; k<=eta_ord; k++) {
        // cout << "k = " << k << endl;
        for (int i=0; i<b_len; i++) {
          // cout << "i = " << i << endl;
          for (int j=0; j<dim2; j++) {
            mpc_set_d(solutions[n][0][i][j][k], 0, MPFR_RNDN);
          }
        }
      }
    }
    // cout << "hit new return 2" << endl;
    goto goto_return;
    return;
  }
  for (int lam, n=0; n<num_cum_eig; n++) {
    lam = cum_eig[n];
    // cout << "lam = " << lam << endl;
    // cout << "sol_log_len = " << sol_log_len[lam] << endl;
    // if (particular) {
    //   cout << "rhs_log_len = " << rhs_log_len[lam] << endl;
    // }

    //////
    // ZERO-ORDER LOG COEFFICIENTS
    //////
    // #2BD: do less operations separting solutions for which lambda belongs to that column index
    // from the solutions associated with the other lambda:
    // in the first case the max log power is that of the block
    // and, therefore, it is useless to keep going with matrix multiplication,
    // we already now that after a certain point we get all zeros
    // cout << "fill zero-order logs" << endl;
    if (block_log_len[lam] > 0) {
    // if (1) {
      // cout << "sol log len = " << sol_log_len[lam] << endl;
      for (int l=0; l<sol_log_len[lam]-1; l++) {
        // cout << "l = " << l << endl;
        // multiply the diagonal block of the L0 times the solution at previous log power
        for (int i=0; i<b_len; i++) {
          // cout << "i = " << i << endl;
          for (int j=0; j<dim2; j++) {
            // cout << "j = " << j << endl;
            mpc_set_d(solutions[n][l+1][i][j][0], 0, MPFR_RNDN);
            // consider the contribution from the diagonal
            mpc_neg(mpc_tmp, lcm[-1+1], MPFR_RNDN);
            // cout << "lam = "; print_mpc(&eig_list[lam]); cout << endl;
            mpc_fma(mpc_tmp, mpc_tmp, eig_list[lam], block[i][i][0], MPFR_RNDN);
            // cout << "ecco:" << endl;
            // print_mpc(&block[i][i][0]); cout << endl;
            // print_mpc(&lcm[-1+1]); cout << endl;
            // print_mpc(&eig_list[lam]); cout << endl;
            // print_mpc(&mpc_tmp); cout << endl;
            // print_mpc(&mpc_tmp);
            // cout << endl;
            // print_mpc(&solutions[n][l][i][j][0]);
            // cout << endl;
            mpc_fma(solutions[n][l+1][i][j][0], mpc_tmp, solutions[n][l][i][j][0], solutions[n][l+1][i][j][0], MPFR_RNDN);

            // #2BD: si può eliminare mettendo semplicemente a=i+1 e un if che controlla se i = b_len-1
            for (int a=i+1; a<b_len; a++) {
              mpc_fma(solutions[n][l+1][i][j][0], block[i][a][0], solutions[n][l][a][j][0], solutions[n][l+1][i][j][0], MPFR_RNDN);
            }
          }

          // CONSTANT TERM
          if (particular && sb_len > 0 && rhs_log_len[lam] > l) {
            // cout << "constant term" << endl;
            // print_mpc(&const_term[n][l][i][0]);
            // cout << endl;
            mpc_add(solutions[n][l+1][i][0][0], solutions[n][l+1][i][0][0], const_term[n][l][i][0], MPFR_RNDN);
          }
        }

        // divide by the constant term of the LCM
        for (int i=0; i<b_len; i++) {
          for (int j=0; j<dim2; j++) {
            mpc_div(solutions[n][l+1][i][j][0], solutions[n][l+1][i][j][0], lcm[-1+1], MPFR_RNDN);
          }
        }
      }
    // }
    } else {
      for (int l=1; l<sol_log_len[lam]; l++) {
        for (int i=0; i<b_len; i++) {
          for (int j=0; j<dim2; j++) {
            mpc_set_d(solutions[n][l][i][j][0], 0, MPFR_RNDN);
          }
        }
      }
    }

    //////
    // HIGHER-ORDER ETA COEFFICIENTS
    //////
    // cout << "fill eta powers" << endl;
    k_start = 1;
    if (block_log_len[lam] == 0) {
      if (particular) {
        k_start = 0;
      } else {
        continue;
      }
    }
    // k_start = 1;
    int i_print = 1, j_print = 1;
    if (print) cout << "k_start = " << k_start << endl;
    for (int l=sol_log_len[lam]-2; l>=0; l--) {
      if (print) cout << "-------- l = " << l << endl;
      for (int k=k_start; k<=eta_ord; k++) {
        if (print && k<5) cout << "---- k = " << k  << endl;
        for (int i=0; i<b_len; i++) {
          // cout << "i = " << i << endl;
          for (int j=0; j<dim2; j++) {
            // cout << "j = " << j << endl;
            mpc_set_d(solutions[n][l][i][j][k], 0, MPFR_RNDN);
            mpc_set_d(mpc_mat[i][j], 0, MPFR_RNDN);
            // ptr = &solutions[lam][l][k][i][j];
            // multiply by the constant term of the LCM
            // print_mpc(&solutions[lam][l+1][k][i][j]);
            // print_mpc(&lcm[0]);
            // if (l == sol_log_len[lam]-2 && k != 0) {
            //   // cout << "1st way" << endl;
            //   mpc_set_d(mpc_mat[i][j], 0, MPFR_RNDN);
            // } else {
            if (l != sol_log_len[lam]-2 || k == 0) {
              // cout << "2nd way" << endl;
              // print_mpc(&solutions[n][l+1][i][j][k]);
              // cout << endl;
              if (print && k<5 && i == i_print && j == j_print) cout << "contribution from l+1, k" << endl;
              mpc_mul(mpc_mat[i][j], solutions[n][l+1][i][j][k], lcm[-1+1], MPFR_RNDN);
              // print_mpc(&mpc_mat[i][j]); cout << endl;
            }
            if (print && k<5 && i == i_print && j == j_print) {cout << "cum = "; print_mpc(&mpc_mat[i][j]); cout << endl;}
            // cumulate LCM times higher log and lower eta-orders
            tmp_pow = MIN(k, lcm_deg+1-1);
            if (print && k<5 && i == i_print && j == j_print) cout << "cumulate LCM" << endl;
            if (print && k<5 && i == i_print && j == j_print) cout << "tmp_pow = " << tmp_pow << endl;
            for (int kp=1; kp<=tmp_pow; kp++) {
              if (print && k<5 && i == i_print && j == j_print) cout << "kp = " << kp << endl;
              if (l == sol_log_len[lam]-2 && k-kp>0) {
                // cout << "continue" << endl;
                continue;
              }
              if (print && k<5 && i == i_print && j == j_print) {
                cout << "picked" << endl;
                cout << "sol = "; print_mpc(&solutions[n][l+1][i][j][k-kp]); cout << endl;
                cout << "lcm = "; print_mpc(&lcm[-1+kp+1]); cout << endl;               
              }

              // if (print && k<5 && i == i_print && j == j_print) {print_mpc(&solutions[n][l+1][i][j][k-kp]); cout << endl;}
              // print_mpc(&lcm[-1+kp+1]);
              // cout << endl;
              mpc_fma(mpc_mat[i][j], lcm[-1+kp+1], solutions[n][l+1][i][j][k-kp], mpc_mat[i][j], MPFR_RNDN);
            }
            if (print && k<5 && i == i_print && j == j_print) {cout << "cum = "; print_mpc(&mpc_mat[i][j]); cout << endl;}

            if (print && k<5 && i == i_print && j == j_print) cout << "cumulate block contribution" << endl;
            mpc_neg(mpc_mat[i][j], mpc_mat[i][j], MPFR_RNDN);
            tmp_pow = MIN(k, block_max_deg);
            if (print && k<5 && i == i_print && j == j_print) cout << "tmp_pow = " << tmp_pow << endl;
            for (int kp=1; kp<=tmp_pow; kp++) {
              if (print && k<5 && i == i_print && j == j_print) cout << "kp = " << kp << endl;

              // // contribution from the diagonal of the block
              // // print_mpc(&solutions[n][l][i][j][k-kp]);
              // // cout << endl;
              // if (kp <= lcm_deg+1-1) {
              //   mpc_add_ui(mpc_tmp, eig_list[lam], k-kp, MPFR_RNDN);
              //   mpc_neg(mpc_tmp, mpc_tmp, MPFR_RNDN);
              //   mpc_fma(mpc_tmp, mpc_tmp, lcm[-1+kp+1], block[i][i][kp], MPFR_RNDN);
              //   mpc_fma(mpc_mat[i][j], mpc_tmp, solutions[n][l][i][j][k-kp], mpc_mat[i][j], MPFR_RNDN);
              // } else {
              //   mpc_fma(mpc_mat[i][j], block[i][i][kp], solutions[n][l][i][j][k-kp], mpc_mat[i][j], MPFR_RNDN);
              // }
              // for (int a=i+1; a<b_len; a++) {
              //   // cout << "a = " << a << endl;
              //   // print_mpc(&solutions[n][l][a][j][k-kp]);
              //   // cout << endl;
              //   mpc_fma(mpc_mat[i][j], block[i][a][kp], solutions[n][l][a][j][k-kp], mpc_mat[i][j], MPFR_RNDN);
              // }
              
              for (int a=0; a<b_len; a++) {
                if (print && k<5 && i == i_print && j == j_print) cout << "a = " << a << endl;
                // print_mpc(&solutions[n][l][a][j][k-kp]);
                // cout << endl;
                if (i == a && kp <= lcm_deg+1-1) {
                  mpc_add_ui(mpc_tmp, eig_list[lam], k-kp, MPFR_RNDN);
                  mpc_neg(mpc_tmp, mpc_tmp, MPFR_RNDN);
                  mpc_fma(mpc_tmp, mpc_tmp, lcm[-1+kp+1], block[i][i][kp], MPFR_RNDN);
                  if (print && k<5 && i == i_print && j == j_print) {
                    cout << "using lcm on the diagonal" << endl;
                    cout << "sol = "; print_mpc(&solutions[n][l][i][j][k-kp]); cout << endl;
                    cout << "block = "; print_mpc(&block[i][i][kp]); cout << endl;
                    cout << "block tilde = "; print_mpc(&mpc_tmp); cout << endl;
                  }
                  mpc_fma(mpc_mat[i][j], mpc_tmp, solutions[n][l][i][j][k-kp], mpc_mat[i][j], MPFR_RNDN);
                  // print_mpc(&mpc_tmp);
                  // cout << endl;
                  // print_mpc(&solutions[n][l][i][j][k-kp]);
                  // cout << endl;
                } else {
                  if (print && k<5 && i == i_print && j == j_print) {
                    cout << "not using lcm" << endl;
                    cout << "sol = "; print_mpc(&solutions[n][l][a][j][k-kp]); cout << endl;
                    cout << "block = "; print_mpc(&block[i][a][kp]); cout << endl;
                  }
                  mpc_fma(mpc_mat[i][j], block[i][a][kp], solutions[n][l][a][j][k-kp], mpc_mat[i][j], MPFR_RNDN);  
                }
                if (print && k<5 && i == i_print && j == j_print) {cout << "cum = "; print_mpc(&mpc_mat[i][j]); cout << endl;}
              }
            }
            mpc_neg(mpc_mat[i][j], mpc_mat[i][j], MPFR_RNDN);
            if (print && k<5 && i == i_print && j == j_print) {cout << "cum = "; print_mpc(&mpc_mat[i][j]); cout << endl;}

            if (lcm_deg+1-1 > tmp_pow) {
              tmp_pow2 = MIN(k, lcm_deg+1-1);
              if (print && k<5 && i == i_print && j == j_print) cout << "cumulate block-lcm contribution" << endl;
              if (print && k<5 && i == i_print && j == j_print) cout << "tmp_pow2 = " << tmp_pow2 << endl;
              for (int kp=tmp_pow+1; kp<=tmp_pow2; kp++) {
                if (print && k<5 && i == i_print && j == j_print) cout << "kp = " << kp << endl;
                mpc_add_ui(mpc_tmp, eig_list[lam], k-kp, MPFR_RNDN);
                mpc_mul(mpc_tmp, mpc_tmp, lcm[-1+kp+1], MPFR_RNDN);
                // cout << "lcm: "; print_mpc(&lcm[-1+kp+1]); cout << endl;
                // cout << "lam + ui: "; print_mpc(&mpc_tmp); cout << endl;
                mpc_fma(mpc_mat[i][j], mpc_tmp, solutions[n][l][i][j][k-kp], mpc_mat[i][j], MPFR_RNDN);
              }
            }
            if (print && k<5 && i == i_print && j == j_print) {cout << "cum = "; print_mpc(&mpc_mat[i][j]); cout << endl;}

            // cout << "SUM OVER PREVIOUS ETA ORDERS" << endl;
            // mpc_neg(mpc_mat[i][j], mpc_mat[i][j], MPFR_RNDN);
            // tmp_pow = MIN(k, min_pow);
            // cout << "tmp_pow = " << tmp_pow << endl;
            // cout << "kp < min" << endl;
            // for (int kp=1; kp<=tmp_pow; kp++) {
            //   cout << "kp = " << kp << endl;
            //   // contribution from the diagonal of the block
            //   mpc_add_ui(mpc_tmp, eig_list[lam], k-kp, MPFR_RNDN);
            //   mpc_neg(mpc_tmp, mpc_tmp, MPFR_RNDN);
            //   mpc_fma(mpc_tmp, mpc_tmp, lcm[-1+kp+1], block[i][i][kp], MPFR_RNDN);
            //   mpc_fma(mpc_mat[i][j], mpc_tmp, solutions[n][l][i][j][k-kp], mpc_mat[i][j], MPFR_RNDN);
            //   for (int a=i+1; a<b_len; a++) {
            //     // cout << "a = " << a << endl;
            //     mpc_fma(mpc_mat[i][j], block[i][a][kp], solutions[n][l][a][j][k-kp], mpc_mat[i][j], MPFR_RNDN);
            //   }
            // }
            // tmp_pow = MIN(k, max_pow);
            // cout << "tmp_pow = " << tmp_pow << endl;
            // cout << "kp > min" << endl;
            // if (lcm_greater) {
            //   cout << "lcm greater" << endl;
            //   for (int kp=min_pow+1; kp<=tmp_pow; kp++) {
            //     cout << "kp = " << kp << endl;
            //     mpc_add_ui(mpc_tmp, eig_list[lam], k-kp, MPFR_RNDN);
            //     mpc_mul(mpc_tmp, mpc_tmp, lcm[-1+kp+1], MPFR_RNDN);
            //     mpc_neg(mpc_tmp, mpc_tmp, MPFR_RNDN);
            //     mpc_fma(mpc_mat[i][j], mpc_tmp, solutions[n][l][i][j][k-kp], mpc_mat[i][j], MPFR_RNDN);
            //   }
            // } else {
            //   cout << "lcm smaller" << endl;
            //   for (int kp=min_pow+1; kp<=tmp_pow; kp++) {
            //     cout << "kp = " << kp << endl;
            //     for (int a=i; a<b_len; a++) {
            //       mpc_fma(mpc_mat[i][j], block[i][a][kp], solutions[n][l][a][j][k-kp], mpc_mat[i][j], MPFR_RNDN);
            //     }
            //   }
            // }
            // mpc_neg(mpc_mat[i][j], mpc_mat[i][j], MPFR_RNDN);
          }
          // CONSTANT TERM
          // if (particular) {
          //   cout << "CONSTANT TERM" << endl;
          //   cout << "l = " << l << endl;
          //   cout << "rhs log len = " << rhs_log_len[lam] << endl;
          //   cout << "sb_len = " << sb_len << endl;
          //   cout << "particular = " << particular << endl;
          // }
          if (particular && sb_len > 0 && rhs_log_len[lam] > l) {
            // if (k == 0 && l == 0) {
            //   print_mpc(&const_term[n][l][i][k]);
            //   cout << endl;
            // }
            mpc_sub(mpc_mat[i][0], mpc_mat[i][0], const_term[n][l][i][k], MPFR_RNDN);
          }
        }

        // INVERSE OF THE LO
        // #2BD: move computation of the inverse ouside of loop over l
        //// compute the matrix to be inverted
        // copy_rk2_mpc(mpc_mat1, block[0], b_len, b_len);
        ////// modify the diagonal
        for (int i=0; i<b_len; i++) {
          for (int j=0; j<b_len; j++) {
            if (i == j) {
              // cout << "ecco:" << endl;
              // print_mpc(&block[i][i][0]); cout << endl;
              // print_mpc(&lcm[-1+1]); cout << endl;
              // print_mpc(&eig_list[lam]); cout << endl;
              mpc_add_ui(mpc_tmp, eig_list[lam], k, MPFR_RNDN);
              mpc_neg(mpc_tmp, mpc_tmp, MPFR_RNDN);
              mpc_fma(mpc_mat1[i][i], mpc_tmp, lcm[-1+1], block[i][i][0], MPFR_RNDN);
            } else {
              mpc_set(mpc_mat1[i][j], block[i][j][0], MPFR_RNDN);
            }
          }
        }
        //// compute the inverse
        if (print && k<5) {
          cout << "LO to invert" << endl;
          cout << "particular = " << particular << endl;
          cout << "block log len = " << block_log_len[lam] << endl;
          print_rk2_mpc(mpc_mat1, b_len, b_len);
        }
        mp_inverse(mpc_mat2, mpc_mat1, b_len);
        if (print && k<5) {
          cout << "B0 inverse:" << endl;
          print_rk2_mpc(mpc_mat2, b_len, b_len);
        }

        //// multiply the inverse
        // copy_rk2_mpc(mpc_mat, solutions[lam][l][k], b_len, b_len);
        for (int i=0; i<b_len; i++) {
          for (int j=0; j<dim2; j++) {
            mpc_set_d(solutions[n][l][i][j][k], 0,  MPFR_RNDN);
            for (int a=0; a<b_len; a++) {
              mpc_fma(solutions[n][l][i][j][k],
                mpc_mat2[i][a], mpc_mat[a][j], solutions[n][l][i][j][k],
                MPFR_RNDN
              );
            }
          }
        }

        if (print && k<5) {
          cout << "solution:" << endl;
          for (int i=0; i<b_len; i++) {
            for (int j=0; j<dim2; j++) {
              cout << "i, j = " << i << ", " << j << ": ";
              print_mpc(&solutions[n][l][i][j][k]); cout << endl;
            }
            cout << endl;
          }
        }
      }
    }
  }

  goto_return:
  // FREE
  mpc_clear(mpc_tmp);
  mpc_rk2_clear(mpc_mat, b_len, dim2);
  del_rk2_tens(mpc_mat, b_len);
  mpc_rk2_clear(mpc_mat1, b_len, b_len);
  del_rk2_tens(mpc_mat1, b_len);
  mpc_rk2_clear(mpc_mat2, b_len, b_len);
  del_rk2_tens(mpc_mat2, b_len);
  
  return;
}


void match_next_eta(
  mpc_t *eta_values, int et,
  mpc_t **solutions, int dim, int eta_ord
) {
  mpc_t eta_diff;
  mpc_init3(eta_diff, wp2, wp2);
  // cln_to_mpc(&eta_diff, eta_values[et+1] - eta_values[et]);
  mpc_sub(eta_diff, eta_values[et+1], eta_values[et], MPFR_RNDN);
  mpc_t eta_power;
  mpc_init3(eta_power, wp2, wp2);
  mpc_set_ui(eta_power, 1, MPFR_RNDN);
  mpc_t tmp_match;
  mpc_init3(tmp_match, wp2, wp2);
  for (int j=1; j<eta_ord; j++) {
    mpc_mul(eta_power, eta_power, eta_diff, MPFR_RNDN);
    for (int i=0; i<dim; i++) {
      // #2BD: substitute with flop and eliminate tmp_match
      mpc_mul(tmp_match, solutions[i][j], eta_power, MPFR_RNDN);
      mpc_add(solutions[i][0], solutions[i][0], tmp_match, MPFR_RNDN);
    }
  }
}


void propagate_regular(
  // IN-OUT
  mpc_t ***solutions,
  // INPUT
  int neta_values, mpc_t *eta_values,
  int dim, struct poly_frac **mat_ep,
  int npoles, mpc_t *poles,
  int nblocks, int **prof, int **sb_grid, int eta_ord
) {
  int print = 0;
  int prt_ord = 10;

  //////
  // START ETA VALUES
  //////
  // numeric etav;
  mpc_t eta_shift;
  mpc_init3(eta_shift, wp2, wp2);
  int b_len, sb_len, *b_idx, *sb_idx;
  int *lcm_deg, *block_max_deg, *num_max_deg, *den_max_deg;
  lcm_deg = new int[nblocks];
  block_max_deg = new int[nblocks];
  num_max_deg = new int[nblocks];
  den_max_deg = new int[nblocks];
  bool flag_lcm = true;
  // ex *gnc_lcm;
  // gnc_lcm = new ex[nblocks];
  mpc_t **lcm, ***block, ****const_term;
  lcm = new mpc_t*[nblocks];
  mpc_t **lcm_orig;
  lcm_orig = new mpc_t*[nblocks];
  int **lcm_roots;
  lcm_roots = new int*[nblocks];
  const_term = new mpc_t***[1];
  const_term[0] = new mpc_t**[1];
  mpc_t ****mpc_lcm_mat_ep;
  malloc_rk3_tens(mpc_lcm_mat_ep, dim, dim, 2);

  // copy mat_ep
  struct poly_frac **lcm_mat_ep;
  malloc_rk2_tens(lcm_mat_ep, dim, dim);
  poly_frac_rk2_build(lcm_mat_ep, dim, dim);
  poly_frac_rk2_set_pf_rk2(lcm_mat_ep, mat_ep, dim, dim);
  // poly_frac_rk2_mul_sym_pow(lcm_mat_ep, 1+p_rank, dim, dim);

  for (int et=0; et<neta_values-1; et++) {
  // for (int et=0; et<1; et++) {
    // etav = (numeric) eta_values[et].real() + I * (numeric) eta_values[et].imag();
    if (print) {
      cout << "eta_" << et << " = "; print_mpc(&eta_values[et]); cout << endl;
    }
    //////
    // BLOCK LOOP
    //////
    for (int b=0; b<nblocks; b++) {
      // build block indices
      build_block_indices(b_idx, sb_idx, &b_len, &sb_len, b, prof, sb_grid);
      if (print) {
      cout << endl;
      cout << "-----------------" << endl;
      cout << "b = " << b << endl;
      cout << "b_len = " << b_len << endl;
      cout << "sb_len = " << sb_len << endl;

      cout << endl;
      cout << "input poly_frac:" << endl;
      for (int i=0; i<b_len; i++) {
        for (int j=0; j<b_len; j++) {
          cout << "i, j = " << i << ", " << j << endl;
          poly_frac_print(&mat_ep[prof[b][0]+i][prof[b][0]+j]);
        }
      }

      cout << endl;
      cout << "boundary:"  << endl;
      for (int i=0; i<b_len; i++) {
        print_mpc(&(*solutions)[prof[b][0]+i][0]); cout << endl;
      }
      }

      // coefficients of LCM, block and constant term
      sys_block_info(
        lcm[b], block, const_term,
        &flag_lcm, &lcm_deg[b], lcm_orig[b], lcm_roots[b],
        &block_max_deg[b], &num_max_deg[b], &den_max_deg[b],
        b_len, sb_len,
        b_idx, sb_idx,
        lcm_mat_ep, mpc_lcm_mat_ep,
        npoles, poles,
        eta_values[et],
        eta_ord,
        &solutions
      );
      
      if (print) {
      cout << endl;
      cout << "block LCM:" << endl;
      print_poly(lcm[b], lcm_deg[b]);
      cout << endl;
      cout << "block numerator:" << endl;
      print_rk3_mpc(block, b_len, b_len, block_max_deg[b]+1);
      cout << endl;
      if (sb_len > 0) {
        cout << endl; cout << "constant term: " << endl;
        print_rk2_mpc(const_term[0][0], b_len, prt_ord-1);
      }
      }

      solve_block(*solutions,
        lcm[b], block, const_term[0][0],
        b_idx, b_len, sb_len,
        eta_ord, lcm_deg[b], block_max_deg[b]);
      
      if (print) {
      cout << endl;
      cout << "block solution:"  << endl;
      for (int i=0; i<b_len; i++) {
        cout << "i = " << i << endl;
        print_poly((*solutions)[prof[b][0]+i], 5);
      }
      }
      
      // FREE
      mpc_rk1_clear(lcm[b], lcm_deg[b]+1);
      delete[] lcm[b];
      if (sb_len>0) {
        mpc_rk2_clear(const_term[0][0], b_len, eta_ord+1);
        del_rk2_tens(const_term[0][0], b_len);
      }
      delete[] b_idx;
      delete[] sb_idx;
      // delete[] lcm;
      // del_rk3_tens(block, b_len, b_len);
      for (int i=0; i<b_len; i++) {
        for (int j=0; j<b_len; j++) {
          mpc_rk1_clear(block[i][j], block_max_deg[b]+1);
          delete[] block[i][j];
        }
        delete[] block[i];
      }
      delete[] block;
    }
    flag_lcm = false;
    // print_poly(solutions[15], 5);
    // print_poly(solutions[15]+eta_ord-3, 3);

    // MATCH WITH NEXT ETA VALUE
    match_next_eta(eta_values, et, *solutions, dim, eta_ord);

    // cout << endl << "MI computed in eta_" << et+1 << ":" << endl;
    // for (int i=0; i<dim; i++) {
    //   cout << "master n." << i << ": ";
    //   print_mpc(&(*solutions)[i][0]);
    //   cout << endl;
    // }
  }

  // FREE
  for (int b=0; b<nblocks; b++) {
    build_block_indices(b_idx, sb_idx, &b_len, &sb_len, b, prof, sb_grid);
    for (int i=0; i<b_len; i++) {
      for (int j=0; j<b_len; j++) {
        mpc_rk1_clear(mpc_lcm_mat_ep[b_idx[i]][b_idx[j]][0], block_max_deg[b]+1);
        delete[] mpc_lcm_mat_ep[b_idx[i]][b_idx[j]][0];
      }
      for (int j=0; j<sb_len; j++) {
        mpc_rk1_clear(mpc_lcm_mat_ep[b_idx[i]][sb_idx[j]][0], num_max_deg[b]+1);
        delete[] mpc_lcm_mat_ep[b_idx[i]][sb_idx[j]][0];
        mpc_rk1_clear(mpc_lcm_mat_ep[b_idx[i]][sb_idx[j]][1], den_max_deg[b]+1);
        delete[] mpc_lcm_mat_ep[b_idx[i]][sb_idx[j]][1];
      }
    }
  }
  del_rk3_tens(mpc_lcm_mat_ep, dim, dim);

  poly_frac_rk2_free(lcm_mat_ep, dim, dim);
  del_rk2_tens(lcm_mat_ep, dim);
  
  for (int b=0; b<nblocks; b++) {
    mpc_rk1_clear(lcm_orig[b], lcm_deg[b]+1);
    delete[] lcm_orig[b];
    delete[] lcm_roots[b];
  }
  delete[] lcm;
  delete[] lcm_orig;
  delete[] lcm_roots;
  del_rk2_tens(const_term, 1);

  // print solutions at last eta value
  // etav = (numeric) eta_values[neta_values-1].real() + I * (numeric) eta_values[neta_values-1].imag();
  // cout << endl << "MI computed in eta_" << neta_values-1 << ": " << endl;
  // for (int i=0; i<dim; i++) {
  //   cout << "master n." << i << ": ";
  //   print_mpc(&(*solutions)[i][0]);
  //   cout << endl;
  // }
  //
  //////

}


void get_block_log_len(
  // OUTPUT
  int **block_log_prof, int *block_log_len, int *rhs_log_len,
  // INPUT
  int b_len, int num_classes, int offset, int *eig_grid, int *eq_class
) {
  /*
  Compute the length of every Jordan chain in the diagonal block.
  Assign to every eigenvalue a log-length equal to L+1,
  L being the max length among the Jordan chains associated with
  that eigenvalue.
  Build also a log_profile, where each column gets the length of
  its Jordan chain plus one.
  */

  int chain_len, curr_eig;

  // init block log profile
  for (int n=0; n<num_classes; n++) {
    block_log_prof[n] = new int[b_len];
    for (int i=0; i<b_len; i++) {
      block_log_prof[n][i] = 0;
    }
  } 

  // get log-lenghts of current block
  for (int i=0; i<num_classes; i++) {
    block_log_len[i] = 0;
    rhs_log_len[i] = 0;
  }
  curr_eig = eq_class[offset];
  chain_len = 1;
  // cout << "curr_eig = " << curr_eig << endl;
  for (int i=1; i<b_len; i++) {
    if (eig_grid[offset+i] != 0) {
      // close the previous chain and register its length
      // (if the new length is greater than the old one for the same eigenvalue)
      if (block_log_len[curr_eig] < chain_len + 1) {
        block_log_len[curr_eig] = chain_len + 1;
      }
      // assign the length of the chain to every index of the chain for the profile
      for (int ip=i-1; ip>=i-chain_len; ip--) {
        block_log_prof[curr_eig][ip] = chain_len + 1;
      }
      // select new eigenvalues as current eigenvalue and start a new chain
      curr_eig = eq_class[offset+i];
      chain_len = 1;
    } else if (eig_grid[offset+i] == 0) {
      chain_len++;
    }
  }
  // close last chain
  if (block_log_len[curr_eig] < chain_len + 1) {
    block_log_len[curr_eig] = chain_len + 1;
  }
  for (int ip=b_len-1; ip>=b_len-chain_len; ip--) {
    block_log_prof[curr_eig][ip] = chain_len + 1;
  }
  // cout << "block log-lenghts" << endl;
  // for (int i=0; i<num_classes; i++) {
  //   cout << "eq class = " << i << ", length = " << block_log_len[i] << endl;
  // }
  // cout << endl;
  // cout << "block log-profile" << endl;
  // for (int lam=0; lam<num_classes; lam++) {
  //   cout << "lam = " << lam << endl;
  //   for (int i=0; i<b_len; i++) {
  //     cout << "i = " << i << ", length = " << block_log_prof[lam][i] << endl;
  //   }
  // }
  // cout << endl;

}


int count_log_len(mpfr_t mpfr_tol, mpc_t ***tens, int dim1, int dim2, int dim3, int *sb_idx) {
  mpfr_t tmp;
  mpfr_init2(tmp, wp2);
  bool found_non_zero;
  int i1;
  // cout << "tol: "; mpfr_out_str(stdout, 10, 0, mpfr_tol, MPFR_RNDN); cout << endl;
  for (i1=0; i1<dim1; i1++) {
    found_non_zero=false;
    for (int i3=0; i3<dim3; i3++) {
      for (int i2=0; i2<dim2; i2++) {
        // check whether tens[i2][i3] is zero
        mpfr_abs(tmp, mpc_realref(tens[i1][sb_idx[i2]][i3]), MPFR_RNDN);
        // mpfr_div_d(tmp, tmp, 1e50, MPFR_RNDN);
        if (mpfr_greater_p(tmp, mpfr_tol)) {
          // cout << "found non-zero real: "; mpfr_out_str(stdout, 10, 0, tmp, MPFR_RNDN); cout << endl;
          // cout << "i1, i2, i3 = " << i1 << ", " << i2 << ", " << i3 << ", " << endl;
          found_non_zero = true;
          break;
        } else {
          mpfr_abs(tmp, mpc_imagref(tens[i1][sb_idx[i2]][i3]), MPFR_RNDN);
          // mpfr_div_d(tmp, tmp, 1e50, MPFR_RNDN);
          if (mpfr_greater_p(tmp, mpfr_tol)) {
            // cout << "found non-zero imag: "; mpfr_out_str(stdout, 10, 0, tmp, MPFR_RNDN); cout << endl;
            // cout << "i1, i2, i3 = " << i1 << ", " << i2 << ", " << i3 << ", " << endl;
            found_non_zero = true;
            break;
          }
        }
      }
      if (found_non_zero) {
       break;
      }
    }
    if (!found_non_zero) {
      break;
    }
  }
  return i1;
}


void evaluate_series_around_zero(
  mpc_t ***output,
  mpc_t *mpc_eta,
  mpc_t *****solutions,
  // int block_max_deg,
  mpc_t *eig_list, int num_eig, int *eig_labels, int *log_prof, int dim1, int dim2, int eta_ord, int num_classes, int lam0,
  bool rk4 = false, mpc_t ****rk4_sol = NULL,
  int *match_constr_dim = NULL, int **match_constr_idx = NULL, bool particular = true,
  int analytic_continuation = 0
) {
  int print=0;
  // ANALYTIC CONTINUATION CONTROL SWITCH
  if (mpfr_sign_within_tol(mpc_realref(*mpc_eta)) > 0 ||\
    mpfr_sign_within_tol(mpc_imagref(*mpc_eta)) != 0) {
    analytic_continuation = 0;
  }
  if (print) cout << "analytic continuation = " << analytic_continuation << endl;
  // #2BD: riscrivere partendo dall'ultima potenza di eta e di log e facendo solo flop,
  // senza cumulare potenze
  // mpc_t mpc_eta;
  // mpc_init3(mpc_eta, wp2, wp2);
  // cln_to_mpc(&mpc_eta, *eta_val);
  mpc_t log_eta;
  mpc_init3(log_eta, wp2, wp2);
  mpc_t tmpc;
  mpc_init3(tmpc, wp2, wp2);
  // analytic continuation
  mpfr_t mpfr_zero;
  mpfr_init2(mpfr_zero, wp2);
  mpfr_set_ui(mpfr_zero, 0, MPFR_RNDN);
  if (analytic_continuation) {
    // cout << "ANALYTIC CONTINUATION: eta = ";
    // print_mpc(mpc_eta);
    // cout << endl;
    mpc_neg(*mpc_eta, *mpc_eta, MPFR_RNDN);
    mpc_log(log_eta, *mpc_eta, MPFR_RNDN);
    mpc_set_ui(tmpc, 0, MPFR_RNDN);
    mpfr_const_pi(mpc_imagref(tmpc), MPFR_RNDN);
    if (analytic_continuation == -1) {
      mpc_neg(tmpc, tmpc, MPFR_RNDN);
    }
    mpc_add(log_eta, log_eta, tmpc, MPFR_RNDN);
    // mpc_sub(log_eta, log_eta, tmpc, MPFR_RNDN);
    mpc_neg(*mpc_eta, *mpc_eta, MPFR_RNDN);
    // cout << "log = ";
    // print_mpc(&log_eta);
    // cout << endl;
  } else {
    mpc_log(log_eta, *mpc_eta, MPFR_RNDN);
  }
  if (print) {cout << "log(eta) = "; print_mpc(&log_eta); cout << endl;}
  // mpc_set_d(log_eta, 0, MPFR_RNDN);
  mpc_t log_pow, eta_pow;
  mpc_init3(log_pow, wp2, wp2);
  mpc_set(log_pow, log_eta, MPFR_RNDN);
  mpc_init3(eta_pow, wp2, wp2);
  mpc_set_d(eta_pow, 1, MPFR_RNDN);
  mpc_t eta_lam;
  mpc_init3(eta_lam, wp2, wp2);

  if (particular) {
    match_constr_dim = new int[num_classes];
    match_constr_idx = new int*[num_classes];
    for (int lam, n=0; n<num_eig; n++) {
      lam = eig_labels[n];
      match_constr_dim[lam] = 1;
      match_constr_idx[lam] = new int[1];
      match_constr_idx[lam][0] = 0;
    }
  }

  // DISTRIBUITE LOG POWERS
  if (print) cout << "DISTRIBUITE LOG POWERS" << endl;
  int *max_log = new int[num_eig];
  mpc_t *****sol_times_log = new mpc_t****[num_eig];
  mpc_t *****in_tens;
  for (int lam, n=0; n<num_eig; n++) {
    lam = eig_labels[n];
    // cout << "lam = " << lam << endl;
    // if (particular && rk4) {
    //   log_prof[lam]--;
    // }
    // cout << "log_prof = " << log_prof[lam] << endl;
    max_log[n] = 0;
    for (int i=0; i<dim1; i++) {
      if (log_prof[lam] > max_log[n]) {
        max_log[n] = log_prof[lam];
      }
    }
  }

  // point input tensor to the correct variable
  // (if rk4 is true it adds a dummy column dimension)
  if (particular && rk4) {
    in_tens = new mpc_t****[num_eig];
    // malloc_rk4_tens(in_tens, num_eig, max_log, dim1, dim2);
    for (int lam, n=0; n<num_eig; n++) {
      lam = eig_labels[n];
      in_tens[n] = new mpc_t***[max_log[n]];
      for (int l=0; l<max_log[n]; l++) {
        // malloc_rk3_tens(in_tens[n][l], dim1, dim2, eta_ord+1);
        malloc_rk2_tens(in_tens[n][l], dim1, dim2);
        for (int i=0; i<dim1; i++) {
          for (int j=0; j<dim2; j++) {
            in_tens[n][l][i][j] = rk4_sol[lam][l][i];
          }
        }
      }
    }
  } else {
    in_tens = solutions;
  }

  for (int lam, n=0; n<num_eig; n++) {
    sol_times_log[n] = new mpc_t***[max_log[n]];
    for (int l=0; l<max_log[n]; l++) {
      sol_times_log[n][l] = new mpc_t**[dim1];
      for (int i=0; i<dim1; i++) {
        sol_times_log[n][l][i] = new mpc_t*[dim2];
        for (int j=0; j<dim2; j++) {
          sol_times_log[n][l][i][j] = new mpc_t[eta_ord+1];
          for (int k=0; k<=eta_ord; k++) {
            mpc_init3(sol_times_log[n][l][i][j][k], wp2, wp2);
            if (l == 0) {
              mpc_set(sol_times_log[n][l][i][j][k], in_tens[n][l][i][j][k], MPFR_RNDN);
            }
          }
        }
      }
    }
  }

  int l = 1;
  bool found_log = false;
  while(1) {
    if (print) cout << "l = " << l << endl;
    found_log = false;
    for (int lam, n=0; n<num_eig; n++) {
      lam = eig_labels[n];
      if (print) cout << "lam = " << lam << endl;
      if (print) cout << "log_prof = " << log_prof[lam] << endl;
      for (int i=0; i<dim1; i++) {
        if (print) cout << "i = " << i << endl;        
        if (log_prof[lam] - 1 > l) {
          for (int j=0; j<match_constr_dim[lam]; j++) {
            // cout << "j = " << j << endl;
            // cout << "match_constr_idx[lam][j] = " << match_constr_idx[lam][j] << endl;
            for (int k=0; k<=eta_ord; k++) {
              // cout << "k = " << k << endl;
              //   cout << "l = " << l << endl;
              //   cout << "coeff of max log" << endl;
              //   print_mpc(&in_tens[n][l][i][j][k]);
              //   cout << endl;
              // print_mpc(&log_pow);
              // cout << endl;
              mpc_mul(sol_times_log[n][l][i][match_constr_idx[lam][j]][k], in_tens[n][l][i][match_constr_idx[lam][j]][k], log_pow, MPFR_RNDN);
              // print_mpc(&sol_times_log[n][l][i][match_constr_idx[lam][j]][k]);
              // cout << endl;
            }
            if (print) {
              cout << "LO: ";
              print_mpc(&in_tens[n][l][i][match_constr_idx[lam][j]][0]);
              cout << endl;
            }
            found_log = true;
          }
          if (print) cout << "found log" << endl;
        }
      }
    }
    if (!found_log) {
      if (print) cout << "log not found" << endl;
      break;
    } else {
      mpc_mul(log_pow, log_pow, log_eta, MPFR_RNDN);
      l++;
      mpc_div_ui(log_pow, log_pow, l, MPFR_RNDN);
    }
  }

  // SUM OVER ETA POWERS
  if (print) cout << "SUM OVER ETA POWERS" << endl;
  //// initialize to zero the output
  for (int n=0; n<num_eig; n++) {
    // cout << "n = " << n << endl;
    for (int i=0; i<dim1; i++) {
      for (int j=0; j<dim2; j++) {
        mpc_set_d(output[n][i][j], 0, MPFR_RNDN);
      }
    }
  }
  for (int k=0; k<=eta_ord; k++) {
    // cout << "k = " << k << endl;
    for (int lam, n=0; n<num_eig; n++) {
      lam = eig_labels[n];
      if (print && k == 0) cout << "lam = " << lam << endl;
      for (int i=0; i<dim1; i++) {
        // cout << "i = " << i << endl;
        // cout << "log_prof = " << log_prof[lam] << endl;
        for (int lll=0; lll<log_prof[lam]-1; lll++) {
          if (print && k == 0 && i == 0) cout << "l = " << lll << endl;
          for (int j=0; j<match_constr_dim[lam]; j++) {
            // cout << "match_constr_idx[lam][j] = " << match_constr_idx[lam][j] << endl;
            // print_mpc(&sol_times_log[n][l][i][match_constr_idx[lam][j]][k]);
            // cout << endl;
            mpc_fma(output[n][i][match_constr_idx[lam][j]],
              sol_times_log[n][lll][i][match_constr_idx[lam][j]][k],
              eta_pow,
              output[n][i][match_constr_idx[lam][j]],
              MPFR_RNDN);
          }
        }
      }
    }
    mpc_mul(eta_pow, eta_pow, *mpc_eta, MPFR_RNDN);
  }

  // MULTIPLY BY ETA^LAM
  for (int lam, n=0; n<num_eig; n++) {
    lam = eig_labels[n];
    // cout << "lam = " << lam << endl;
    // if (lam == lam0) {
    //   mpc_set_ui(eta_lam, 1, MPFR_RNDN);
    // } else {
    if (analytic_continuation) {
      // analytic continuation
      //// eta^lam = exp{lam*log(eta)}
      if (print) {
      cout << "ANALYTIC CONTINUATION: eta = "; print_mpc(mpc_eta); cout << endl;
      }
      mpc_mul(eta_lam, eig_list[lam], log_eta, MPFR_RNDN);
      mpc_exp(eta_lam, eta_lam, MPFR_RNDN);
    } else {
      mpc_pow(eta_lam, *mpc_eta, eig_list[lam], MPFR_RNDN);
    }
    if (print) {
    cout << "eta_lam = "; print_mpc(&eta_lam); cout << endl;
    }
    for (int i=0; i<dim1; i++) {
      if (log_prof[lam] == 0) {
        continue;
      }
      for (int j=0; j<match_constr_dim[lam]; j++) {
        // cout << "output = "; print_mpc(&output[n][i][match_constr_idx[lam][j]]); cout << endl;
        mpc_mul(output[n][i][match_constr_idx[lam][j]],
          output[n][i][match_constr_idx[lam][j]],
          eta_lam,
          MPFR_RNDN);
        // cout << "output after = "; print_mpc(&output[n][i][match_constr_idx[lam][j]]); cout << endl;
      }
    }
  }

  // FREE
  for (int lam, n=0; n<num_eig; n++) {
    for (int l=0; l<max_log[n]; l++) {
      for (int i=0; i<dim1; i++) {
        for (int j=0; j<dim2; j++) {
          for (int k=0; k<=eta_ord; k++) {
            mpc_clear(sol_times_log[n][l][i][j][k]);
          }
          delete[] sol_times_log[n][l][i][j];
        }
        delete[] sol_times_log[n][l][i];
      }
      delete[] sol_times_log[n][l];
    }
    delete[] sol_times_log[n];
  }

  if (particular && rk4) {
    for (int lam, n=0; n<num_eig; n++) {
      lam = eig_labels[n];
      for (int l=0; l<max_log[n]; l++) {
        del_rk2_tens(in_tens[n][l], dim1);
      }
      delete[] in_tens[n];
    }
    delete[] in_tens;
  }

  if (particular) {
    for (int lam, n=0; n<num_eig; n++) {
      lam = eig_labels[n];
      delete[] match_constr_idx[lam];
    }
    delete[] match_constr_dim;
    delete[] match_constr_idx;
  }
}


void solve_log_constraints(
  // IN-OUT
  mpc_t *****psol,
  // INPUT
  int eta_ord, int b_len, int offset, int num_cum_eig, int *cum_eig, int *eq_class,
  int *block_log_len, int *sol_log_len, int *rhs_log_len, mpc_t *****hsol
) {
  int print = 0;
  int log_constr_dim = 0;
  int *log_constr_idx = new int[b_len];
  mpc_t *log_constr_coeffs = new mpc_t[b_len];
  init_rk1_mpc(log_constr_coeffs, b_len);
  mpc_t **log_constr_mat,  **log_constr_inv_mat;
  mpc_t tmpc; mpc_init3(tmpc, wp2, wp2);

  int wp_bin = - mpfr_log2_int(mpfr_tol);

  // assemble matrices by considering every block eigenvalue
  // but only in the positions not associated with that eiganvalue
  for (int lam, n=0; n<num_cum_eig; n++) {
    lam = cum_eig[n];
    if (block_log_len[lam] == 0) {
      continue;
    }
    if (print) cout << "lam = " << lam << endl;
    // first count dimensions for allocation and select indices
    log_constr_dim = 0;
    for (int j=0; j<b_len; j++) {
      // cout << "j = " << j << endl;
      // cout << "eq_class = " << eq_class[offset+j] << endl;
      if (eq_class[offset+j] != lam) {
        // cout << "picked" << endl;
        log_constr_idx[log_constr_dim] = j;
        log_constr_dim++;
      }
    }
    malloc_rk2_tens(log_constr_mat, log_constr_dim, log_constr_dim);
    init_rk2_mpc(log_constr_mat, log_constr_dim, log_constr_dim);
    malloc_rk2_tens(log_constr_inv_mat, log_constr_dim, log_constr_dim);
    init_rk2_mpc(log_constr_inv_mat, log_constr_dim, log_constr_dim);
    // cout << "log_constr_dim = " << log_constr_dim << endl;
    if (log_constr_dim == 0) {
      return;
    }
    
    // fill matrix
    // cout << "fill matrix" << endl;
    for (int j=0; j<log_constr_dim; j++) {
      for (int i=0; i<log_constr_dim; i++) {
        // cout << "n = " << n << endl;
        // cout << "log pow = " << sol_log_len[lam] << endl;
        // cout << "i = " << i << endl;
        // cout << "j = " << j << endl;
        // cout << "log i = " << log_constr_idx[i] << endl;
        // cout << "log j = " << log_constr_idx[j] << endl;
        // print_mpc(&hsol[n][sol_log_len[lam]-1][log_constr_idx[i]][log_constr_idx[j]][0]);
        // cout << endl;
        mpc_set(log_constr_mat[i][j], hsol[n][sol_log_len[lam]-1][log_constr_idx[i]][log_constr_idx[j]][0], MPFR_RNDN);
        mpc_neg(log_constr_mat[i][j], log_constr_mat[i][j], MPFR_RNDN);
      }
    }
    
    // solve the system
    //// invert the matrix
    if (print) {
    cout << "solve the system" << endl;
    cout << "linear system matrix:" << endl;
    print_rk2_mpc(log_constr_mat, log_constr_dim, log_constr_dim);
    }
    mp_inverse(log_constr_inv_mat, log_constr_mat, log_constr_dim);
    if (print) {
    cout << "constant term:" << endl;
    for (int j=0; j<log_constr_dim; j++) {
      print_mpc(&psol[n][sol_log_len[lam]-1][log_constr_idx[j]][0][0]); cout << endl;
    }
    }
    //// multiply the inverse times the constant term of the system
    for (int i=0; i<log_constr_dim; i++) {
      mpc_set_d(log_constr_coeffs[i], 0, MPFR_RNDN);
      int j;
      for (j=0; j<log_constr_dim; j++) {
        // mpc_fma(
        //   log_constr_coeffs[i],
        //   log_constr_inv_mat[i][j],
        //   psol[n][sol_log_len[lam]-1][log_constr_idx[j]][0][0],
        //   log_constr_coeffs[i],
        //   MPFR_RNDN
        // );

        exp_rel_err_mpc_mul(tmpc, log_constr_inv_mat[i][j], psol[n][sol_log_len[lam]-1][log_constr_idx[j]][0][0], MPFR_RNDN, wp_bin);
        exp_rel_err_mpc_add(log_constr_coeffs[i], log_constr_coeffs[i], tmpc, MPFR_RNDN, wp_bin);
      }
      
      // // rel err prune on last contribution
      // mpc_mul(tmpc, log_constr_inv_mat[i][j], psol[n][sol_log_len[lam]-1][log_constr_idx[j]][0][0], MPFR_RNDN);
      // rel_err_mpc_add(log_constr_coeffs[i], log_constr_coeffs[i], tmpc, MPFR_RNDN);
    }

    if (print) {
    cout << "solution coefficients:" << endl;
    print_poly(log_constr_coeffs, log_constr_dim-1);
    cout << endl;
    }

    // update particular solution
    if (print) cout << "update particular solution" << endl;
    int higher_log = 1;
    for (int l=sol_log_len[lam]-1; l>=0; l--) {
      int k_start = 0;
      // ZERO ORDER
      if (higher_log && l>=rhs_log_len[lam]) {
        k_start = 1;
        int col_is_zero = 1;
        for (int i=0; i<b_len; i++) {
          int j=0;
          for (j=0; j<log_constr_dim; j++) {
            // mpc_fma(
            //   psol[n][l][i][0][0],
            //   log_constr_coeffs[j],
            //   hsol[n][l][i][log_constr_idx[j]][0],
            //   psol[n][l][i][0][0],
            //   MPFR_RNDN
            // );

            exp_rel_err_mpc_mul(tmpc, log_constr_coeffs[j], hsol[n][l][i][log_constr_idx[j]][0], MPFR_RNDN, wp_bin);
            exp_rel_err_mpc_add(psol[n][l][i][0][0], psol[n][l][i][0][0], tmpc, MPFR_RNDN, wp_bin);
          }

          // rel err prune on last contribution
          // cout << "j = " << j << endl;
          // cout << "log_constr_dim = " << log_constr_dim << endl;
          // cout << "log_constr_idx[j] = " << log_constr_idx[j] << endl;
          // print_mpc(&log_constr_coeffs[j]); cout << endl;
          // print_mpc(&hsol[n][l][i][log_constr_idx[j]][0]); cout << endl;
          // mpc_mul(tmpc, log_constr_coeffs[j], hsol[n][l][i][log_constr_idx[j]][0], MPFR_RNDN);
          // rel_err_mpc_add(psol[n][l][i][0][0], psol[n][l][i][0][0], tmpc, MPFR_RNDN);

          if (!mpc_zero_p(psol[n][l][i][0][0])) {
            col_is_zero = 0;
          }
        }

        if (col_is_zero) {
          for (int i=0; i<b_len; i++) {
            for (int k=1; k<=eta_ord; k++) {
              mpc_set_ui(psol[n][l][i][0][k], 0, MPFR_RNDN);
            }
          }
          continue;
        } else {
          higher_log = 0;
        }
          
      }

      for (int k=k_start; k<=eta_ord; k++) {
        for (int i=0; i<b_len; i++) {
          // mpc_neg(psol[n][l][log_constr_idx[i]][0][k], psol[n][l][log_constr_idx[i]][0][k], MPFR_RNDN);
          for (int j=0; j<log_constr_dim; j++) {
            mpc_fma(
              psol[n][l][i][0][k],
              log_constr_coeffs[j],
              hsol[n][l][i][log_constr_idx[j]][k],
              psol[n][l][i][0][k],
              MPFR_RNDN
            );
          }
          // mpc_neg(psol[n][l][log_constr_idx[i]][0][k], psol[n][l][log_constr_idx[i]][0][k], MPFR_RNDN);
        }
      }
    }
    mpc_rk2_clear(log_constr_mat, log_constr_dim, log_constr_dim);
    del_rk2_tens(log_constr_mat, log_constr_dim);
    mpc_rk2_clear(log_constr_inv_mat, log_constr_dim, log_constr_dim);
    del_rk2_tens(log_constr_inv_mat, log_constr_dim);
  }
  delete[] log_constr_idx;
  delete[] log_constr_coeffs;

  if (print) cout << "solved" << endl;
  // if (1) {
  //   cout << "enter to check new particular solution" << endl;
  //   getchar();
  //   for (int lam, n=0; n<num_cum_eig; n++) {
  //     lam = cum_eig[n];
  //     if (sol_log_len[lam] == 0) {
  //       continue;
  //     }
  //     cout << "lam = " << lam << endl;
  //     for (int l=0; l<sol_log_len[lam]; l++) {
  //       cout << "l = " << l << endl;
  //       print_rk3_mpc(psol[n][l], b_len, 1, 5+1);
  //     }
  //   }
  //   getchar();
  // }

  // FREE
  mpc_clear(tmpc);
}


void build_csvd_mat(
  // OUTPUT
  mpc_t *csvd_mat,
  // INPUT
  int row_dim, int col_dim, mpc_t **constr, mpc_t *const_term,
  int b_len,
  // int num_cum_eig, int *cum_eig,
  // int *block_log_len, int *is_in_bound,
  // int *match_constr_dim, int **match_constr_idx,
  int *coeff_is_null
) {
  // int count_col = 0;
  int count = 0;
  for (int j=0; j<b_len; j++) {
    if (coeff_is_null[j]) {continue;}
    for (int i=0; i<row_dim; i++) {      
      // cout << "i, j, count_col = " << i << ", " << j << ", " << count_col << endl;
      mpc_set(csvd_mat[count++], constr[i][j], MPFR_RNDN);
      // if (i == 0 && j == 0) {
      //   mpc_set_ui(csvd_mat[i+row_dim*count_col], 1, MPFR_RNDN);
      // } else if (i == 1 && j == 0) {
      //   mpc_set_ui(csvd_mat[i+row_dim*count_col], 1, MPFR_RNDN);
      // } else {
      //   mpc_set_ui(csvd_mat[i+row_dim*count_col], 0, MPFR_RNDN);
      // }
      // count_col++;
    }
  }
  // attach constant term to the bottom
  for (int i=0; i<row_dim; i++) {
    mpc_set(csvd_mat[i+row_dim*col_dim], const_term[i], MPFR_RNDN);
  }
}


int compute_rank(
  mpc_t *mat, int row_dim, int col_dim
) {
  int print=0;
  int min_dim = min(row_dim, col_dim);
  if (print) cout << "min dim = " << min_dim << endl;
  mpfr_t told, tmp;
  mpfr_init2(told, wp2);
  mpfr_init2(tmp, wp2);

  // mpc_t *u = new mpc_t[row_dim*row_dim];
  // init_rk1_mpc(u, row_dim*row_dim);
  // mpc_t *v= new mpc_t[col_dim*col_dim];
  // init_rk1_mpc(v, col_dim*col_dim);
  // mpfr_t *s = new mpfr_t[row_dim*col_dim];
  // for (int i=0; i<row_dim*col_dim; i++) {
  //   mpfr_init2(s[i], wp2);
  // }
  mpc_t *u = new mpc_t[100];
  init_rk1_mpc(u, 100);
  mpc_t *v= new mpc_t[100];
  init_rk1_mpc(v, 100);
  mpfr_t *s = new mpfr_t[100];
  for (int i=0; i<100; i++) {
    mpfr_init2(s[i], wp2);
  }

  // CSVD
  csvd1_(mat, &row_dim, &col_dim, &row_dim, &col_dim, &c__0, &row_dim, &col_dim, s, u, v, 100, 1);

  if (print) {
  cout << "singular values" << endl;
  for (int i=0; i<min_dim; i++) {
    mpfr_out_str(stdout, 10, 0, s[i], MPFR_RNDN); cout << endl;
  }
  }

  // determine rank
  int rank = 0;
  mpfr_mul_ui(told, mpfr_tol, row_dim*col_dim, MPFR_RNDN);
  for (int i=0; i<min_dim; i++) {
    if (!mpfr_lessthan_tol(s[i])) {
      rank++;
    }
  }

  return rank;

}


void match_sol_around_zero(
  // OUTPUT
  mpc_t ****sol,
  // INPUT
  int eta_ord, int b_len, int sb_len, int offset, int num_cum_eig, int *cum_eig,
  int num_classes, int *eq_class, mpc_t *eig_list,
  int *block_log_len, int *sol_log_len, int *rhs_log_len,
  mpc_t *****hsol, mpc_t *****psol,
  mpc_t ***solutions, mpc_t *mpc_eta, int lam0,
  int analytic_cont,
  mpc_t ****prev_sol, struct poly_frac **tmat, int nroots, mpc_t *roots,
  int **bound_behav, int **mi_eig, int *mi_eig_num, int **log_prof,
  struct poly_frac **pfmat
) {
  int print = 0;
  mpc_t ***eval_hsol, ***eval_psol;
  int *match_constr_dim = new int[num_classes];
  int **match_constr_idx = new int*[num_classes];
  mpc_t **match_constr_mat, **match_constr_inv_mat, *match_constr_coeffs;

  malloc_rk2_tens(match_constr_mat, b_len, b_len);
  init_rk2_mpc(match_constr_mat, b_len, b_len);
  malloc_rk2_tens(match_constr_inv_mat, b_len, b_len);
  init_rk2_mpc(match_constr_inv_mat, b_len, b_len);
  match_constr_coeffs = new mpc_t[b_len];
  init_rk1_mpc(match_constr_coeffs, b_len);

  int wp_bin = - mpfr_log2_int(mpfr_tol);
  // wp_bin *= 0.90;  // #hard-coded

  // count and select indices for every eigenvalue
  // cout << "count and select indices for every eigenvalue" << endl;
  for (int lam, n=0; n<num_cum_eig; n++) {
    lam = cum_eig[n];
    // cout << "lam = " << lam << endl;
    if (block_log_len[lam] == 0) {
      continue;
    }
    match_constr_idx[lam] = new int[b_len];

    // select column indices associated with the current eigenvalue
    match_constr_dim[lam] = 0;
    for (int j=0; j<b_len; j++) {
      // cout << "j = " << j << endl;
      // cout << "eq_class = " << eq_class[offset+j] << endl;
      if (eq_class[offset+j] == lam) {
        // cout << "picked" << endl;
        match_constr_idx[lam][match_constr_dim[lam]] = j;
        match_constr_dim[lam]++;
      }
    }
  }

  // initialize solution to zero
  for (int lam, n=0; n<num_cum_eig; n++) {
    lam = cum_eig[n];
    for (int l=0; l<sol_log_len[lam]-1; l++) {
      for (int k=0; k<=eta_ord; k++) {
        for (int i=0; i<b_len; i++) {
          mpc_set_d(sol[lam][l][offset+i][k], 0, MPFR_RNDN);
        }
      }
    }
  }
  
  // evaluate homogeneous solution
  int match_in_zero;
  if (mpfr_zero_p(mpc_realref(*mpc_eta)) && mpfr_zero_p(mpc_imagref(*mpc_eta))) {
    match_in_zero = 1;
  } else {
    match_in_zero = 0;
  }
  if (print) cout << "match in zero = " << match_in_zero << endl;

  int ldeg, tmp_pow_behav, ***pow_behav, **pow_behav_min, **const_term_pow_behav;
  int **DE_const_term_pow_behav, **DE_hom_pow_behav;
  malloc_rk2_tens(const_term_pow_behav, b_len, num_cum_eig);
  malloc_rk2_tens(DE_const_term_pow_behav, num_classes, b_len);
  malloc_rk2_tens(DE_hom_pow_behav, num_classes, b_len);
  
  if (match_in_zero == 0) {
    malloc_rk3_tens(eval_hsol, num_cum_eig, b_len, b_len);
    init_rk3_mpc(eval_hsol, num_cum_eig, b_len, b_len);
    if (print) cout << "enter to HOMO SERIES EVALUATED" << endl;
    evaluate_series_around_zero(eval_hsol, mpc_eta, hsol,
      eig_list, num_cum_eig, cum_eig, block_log_len, b_len, b_len, eta_ord, num_classes, lam0,
      false, NULL,
      match_constr_dim, match_constr_idx, false,
      analytic_cont
    );
    if (print) {
    cout << "sum of the homogeneous:" << endl;
    print_rk3_mpc(eval_hsol, num_cum_eig, b_len, b_len);
    }

    // evaluate particular solution
    if (sb_len > 0) {
      malloc_rk3_tens(eval_psol, num_cum_eig, b_len, 1);
      init_rk3_mpc(eval_psol, num_cum_eig, b_len, 1);
      if (print) cout << "enter to PAR SERIES EVALUATED" << endl;
      // getchar();
      // cln_to_mpc(&mpc_eta, *eta_val);
      evaluate_series_around_zero(eval_psol, mpc_eta, psol,
        eig_list, num_cum_eig, cum_eig, sol_log_len, b_len, 1, eta_ord, num_classes, lam0,
        false, NULL,
        NULL, NULL, true,
        analytic_cont
      );
      if (print) {
      cout << "sum of the particular:" << endl;
      print_rk3_mpc(eval_psol, num_cum_eig, b_len, 1);
      }
    }
    
    // fill the matrix of the matching linear system
    for (int lam, n=0; n<num_cum_eig; n++) {
      lam = cum_eig[n];
      // cout << "lam = " << lam << endl;
      if (block_log_len[lam] == 0) {
        continue;
      }
      for (int j=0; j<match_constr_dim[lam]; j++) {
        // fill the entire column
        for (int i=0; i<b_len; i++) {
          mpc_set(match_constr_mat[i][match_constr_idx[lam][j]], eval_hsol[n][i][match_constr_idx[lam][j]], MPFR_RNDN);
        }
      }
    }
    //// build constant term
    // cout << "build constant term" << endl;
    if (sb_len > 0) {
      for (int i=0; i<b_len; i++) {
        for (int lamp, np=0; np<num_cum_eig; np++) {
          lamp = cum_eig[np];
          if (sol_log_len[lamp] == 0) {
            continue;
          }
          // cout << "lamp = " << lamp << endl;
          // cout << "boundary = "; print_mpc(&(*solutions)[offset+i][0]); cout << endl;
          // cout << "particul = "; print_mpc(&eval_psol[np][i][0]); cout << endl;
          mpc_neg(eval_psol[np][i][0], eval_psol[np][i][0], MPFR_RNDN);
          exp_rel_err_mpc_add(
            (*solutions)[offset+i][0],
            (*solutions)[offset+i][0], eval_psol[np][i][0],
            MPFR_RNDN,
            wp_bin
            );
          mpc_neg(eval_psol[np][i][0], eval_psol[np][i][0], MPFR_RNDN);
          // cout << "result = "; print_mpc(&(*solutions)[offset+i][0]); cout << endl;
        }
      }
    }
    mpc_rk3_clear(eval_hsol, num_cum_eig, b_len, b_len);
    del_rk3_tens(eval_hsol, num_cum_eig, b_len);
    if (sb_len) {
      mpc_rk3_clear(eval_psol, num_cum_eig, b_len, 1);
      del_rk3_tens(eval_psol, num_cum_eig, b_len);
    }
  } else {
    if (print) {
    cout << endl; cout << "BOUNDARY BEHAVIOUR:" << endl;
    for (int i=0; i<b_len; i++) {
      cout << "i = " << i << ": {";
      for (int lam, n=0; n<mi_eig_num[offset+i]; n++) {
        lam = mi_eig[offset+i][n];
        cout << "lam" << lam << ": ";
        cout << bound_behav[offset+i][lam] << ", ";
      }
      cout << "}" << endl;
    }
    }

    //////
    // PREPARATION
    //////
    int *col_idx;
    col_idx = new int[b_len];
    for (int jc=0; jc<b_len; jc++) {
      col_idx[jc] = jc;
    }
    int *mi_constr_idx;
    int count_mi_idx;
    mpc_t *free_constr_const_term, *ext_constr_const_term;
    int free_constr;
    pow_behav = new int**[num_classes];
    for (int lam, n=0; n<num_cum_eig; n++) {
      lam = cum_eig[n];
      if (block_log_len[lam] == 0) {
        continue;
      }
      malloc_rk2_tens(pow_behav[lam], b_len, match_constr_dim[lam]);
    }
    malloc_rk2_tens(pow_behav_min, num_classes, b_len);
    malloc_rk2_tens(const_term_pow_behav, num_classes, b_len);
    mpc_t tmpc;
    mpc_init3(tmpc, wp2, wp2);
    struct poly_frac *pf_sol = new struct poly_frac[offset+b_len];
    poly_frac_rk1_build(pf_sol, offset+b_len);
    struct poly_frac *tmp_pf = new struct poly_frac[b_len];
    poly_frac_rk1_build(tmp_pf, b_len);
    // #mem
    struct poly_frac ***tmat_times_hsol;
    tmat_times_hsol = new struct poly_frac**[num_classes];
    for (int lam, n=0; n<num_cum_eig; n++) {
      lam = cum_eig[n];
      if (block_log_len[lam] == 0) {
        continue;
      }
      malloc_rk2_tens(tmat_times_hsol[lam], b_len, match_constr_dim[lam]);
      poly_frac_rk2_build(tmat_times_hsol[lam], b_len, match_constr_dim[lam]);
    }
    struct poly_frac ***mat_times_tmat_times_hsol;// = new struct poly_frac**[b_len];
    mat_times_tmat_times_hsol = new struct poly_frac**[num_classes];
    for (int lam, n=0; n<num_cum_eig; n++) {
      lam = cum_eig[n];
      if (block_log_len[lam] == 0) {
        continue;
      }
      malloc_rk2_tens(mat_times_tmat_times_hsol[lam], b_len, match_constr_dim[lam]);
      poly_frac_rk2_build(mat_times_tmat_times_hsol[lam], b_len, match_constr_dim[lam]);
    }
    struct poly_frac **pf_const_term;
    pf_const_term = new struct poly_frac*[num_classes];
    for (int lam, n=0; n<num_cum_eig; n++) {
      lam = cum_eig[n];
      pf_const_term[lam] = new struct poly_frac[b_len];
      for (int i=0; i<b_len; i++) {
        poly_frac_build(&pf_const_term[lam][i]);
      }
    }
    struct poly_frac **tmat_times_psol;
    tmat_times_psol = new struct poly_frac*[num_classes];
    for (int lam, n=0; n<num_cum_eig; n++) {
      lam = cum_eig[n];
      tmat_times_psol[lam] = new struct poly_frac[b_len];
      for (int i=0; i<b_len; i++) {
        poly_frac_build(&tmat_times_psol[lam][i]);
      }
    }
    struct poly_frac **mat_times_tmat_times_psol;// = new struct poly_frac*[b_len];
    mat_times_tmat_times_psol = new struct poly_frac*[num_classes];
    for (int lam, n=0; n<num_cum_eig; n++) {
      lam = cum_eig[n];
      mat_times_tmat_times_psol[lam] = new struct poly_frac[b_len];
      for (int i=0; i<b_len; i++) {
        poly_frac_build(&mat_times_tmat_times_psol[lam][i]);
      }
    }
    struct poly_frac **mat_times_prevsol;// = new struct poly_frac*[b_len];
    mat_times_prevsol = new struct poly_frac*[num_classes];
    for (int lam, n=0; n<num_cum_eig; n++) {
      lam = cum_eig[n];
      mat_times_prevsol[lam] = new struct poly_frac[b_len];
      for (int i=0; i<b_len; i++) {
        poly_frac_build(&mat_times_prevsol[lam][i]);
      }
    }
    struct poly_frac **DE_const_term ;//= new struct poly_frac*[b_len];
    DE_const_term = new struct poly_frac*[num_classes];
    for (int lam, n=0; n<num_cum_eig; n++) {
      lam = cum_eig[n];
      DE_const_term[lam] = new struct poly_frac[b_len];
      for (int i=0; i<b_len; i++) {
        poly_frac_build(&DE_const_term[lam][i]);
      }
    }

    struct poly_frac pf_tmp;
    poly_frac_build(&pf_tmp);

    int wp_bin = - mpfr_log2_int(mpfr_tol);
    //
    //////

    // contribution for the constant term from previous blocks
    // cout << "GENERATE CONSTANT TERM" << endl;
    tmat += offset;

    int tp_rank = poly_frac_Poinc_rank(tmat, b_len, offset+b_len);
    if (print) {
    cout << "tp rank = " << tp_rank << endl;
    cout << "tmat:" << endl;
    poly_frac_rk2_print(tmat, b_len, offset+b_len);
    }
    if (offset > 0) {
      for (int lam, n=0; n<num_cum_eig; n++) {
        lam = cum_eig[n];
        if (sol_log_len[lam] == 0) {
          continue;
        }
        for (int i=0; i<offset; i++) {
          // cout << "lam, i = " << lam << ", " << i << endl;
          // cout << "log prof = " << log_prof[lam][i] << endl;
          // poly_frac_build(&pf_sol[i]);
          if (log_prof[lam][i] == 0) {
            // cout << "set pf to zero" << endl;
            poly_frac_set_zero(&pf_sol[i]);
          } else {
            // cout << "set pf coeffs from poly" << endl;
            // print_poly(prev_sol[lam][0][i], tp_rank+1);
            poly_frac_set_coeffs(&pf_sol[i], prev_sol[lam][0][i], tp_rank+1);
            pf_sol[i].den_deg = 0;
            pf_sol[i].nroots = nroots;
            pf_sol[i].mults = new int[nroots];
            for (int k=0; k<nroots; k++) {
              pf_sol[i].mults[k] = 0;
            }
          }
          // poly_frac_print(&pf_sol[i]);
        }
        // cout << "mul tmat" << endl;
        // poly_frac_rk2_print(tmat, b_len, offset);
        poly_frac_rk2_mul_pf_rk1(pf_const_term[lam], tmat, pf_sol, b_len, offset, roots);
        poly_frac_rk1_prune_rel_tol(pf_const_term[lam], wp_bin, b_len);

        if (print) {
        cout << "tmat times prev sol:" << endl;
        for (int i=0; i<b_len; i++) {
          cout << "lam, i = " << lam << ", " << i << endl;
          // poly_frac_normal(&pf_const_term[lam][i]);
          poly_frac_print(&pf_const_term[lam][i]);
          // const_term_pow_behav[lam][i] = poly_frac_pow_behav(&pf_const_term[lam][i]);
          // cout << "behav = " << const_term_pow_behav[lam][i] << endl;
        }
        }
      }
    }

    // shift tmat to make it point to the diagonal block
    for (int i=0; i<b_len; i++) {
      tmat[i] += offset;
    }

    if (sb_len>0) {
      if (print) cout << "diag tmat times particular" << endl;
      for (int lamp, np=0; np<num_cum_eig; np++) {
        lamp = cum_eig[np];
        // cout << "lamp = " << lamp << endl;
        if (sol_log_len[lamp] == 0) {
          continue;
        }

        for (int i=0; i<b_len; i++) {
          if (print) cout << "lamp, i = " << lamp << ", " << i << endl;
          // poly_frac_build(&pf_sol[i]);        
          // cout << "set pf coeffs from poly" << endl;
          // print_poly(prev_sol[lam][0][i], tp_rank+1);
          poly_frac_set_coeffs(&pf_sol[i], psol[np][0][i][0], tp_rank+1);
          pf_sol[i].den_deg = 0;
          pf_sol[i].nroots = nroots;
          pf_sol[i].mults = new int[nroots];
          for (int k=0; k<nroots; k++) {
            pf_sol[i].mults[k] = 0;
          }
          if (print) poly_frac_print(&pf_sol[i]);
        }
        if (print) cout << "mul tmat" << endl;
        // poly_frac_rk2_print(tmat, b_len, offset);
        // poly_frac_rk2_mul_pf_rk1(tmat_times_psol, tmat, pf_sol, b_len, b_len, roots);

        for (int i=0; i<b_len; i++) {
          if (print) cout << "lamp, i = " << lamp << ", " << i << endl;
          poly_frac_mul_pf(
            &tmat_times_psol[lamp][i], &tmat[i][0], &pf_sol[0]
          );
          for (int k=1; k<b_len; k++) {
            poly_frac_mul_pf(
              &pf_tmp, &tmat[i][k], &pf_sol[k]
            );
            poly_frac_add_pf(
              &tmat_times_psol[lamp][i], &tmat_times_psol[lamp][i], &pf_tmp,
              roots,
              NULL
            );
          }
          poly_frac_add_pf(
            &pf_const_term[lamp][i], &pf_const_term[lamp][i], &tmat_times_psol[lamp][i],
            roots,
            NULL
          );
          poly_frac_prune_rel_tol(&pf_const_term[lamp][i], wp_bin);
          if (print) poly_frac_print(&pf_const_term[lamp][i]);
        }
      }
    }

    if (sb_len > 0 || offset > 0 ) {
      // POWER BEHAVIOUR OF THE CONSTANT TERM
      if (print) cout << "pow behav of the constant term:" << endl;
      for (int lam, n=0; n<num_cum_eig; n++) {
        lam = cum_eig[n];
        if (sol_log_len[lam] == 0) {
          continue;
        }
        
        for (int i=0; i<b_len; i++) {
          const_term_pow_behav[lam][i] = poly_frac_pow_behav(&pf_const_term[lam][i]);
          if (print) {
          cout << "lam, i = " << lam << ", " << i << ", pow behav = " << const_term_pow_behav[lam][i] << endl;
          poly_frac_print(&pf_const_term[lam][i]);
          }
        }
      }
    } else {
      for (int lam, n=0; n<num_cum_eig; n++) {
        lam = cum_eig[n];
        if (sol_log_len[lam] == 0) {
          continue;
        }
        
        for (int i=0; i<b_len; i++) {
          const_term_pow_behav[lam][i] = 0;
        }
      }
    }

    // BUILD HOMOGENEOUS SOLUTION POLY FRAC
    if (print) cout << "BUILD HOMOGENEOUS SOLUTION POLY FRAC" << endl;
    // initialize min pow behav
    if (print) cout << "pow behav of the homogeneous:" << endl;
    for (int lam, n=0; n<num_cum_eig; n++) {
      lam = cum_eig[n];
      if (block_log_len[lam] == 0) {
        continue;
      }
      for (int i=0; i<b_len; i++) {
        pow_behav_min[lam][i] = tp_rank+1;
      }
    }
    // compute poly_frac and get pow behav, then update min pow behav
    for (int lam, n=0; n<num_cum_eig; n++) {
      lam = cum_eig[n];
      if (block_log_len[lam] == 0) {
        continue;
      }
      for (int jc, j=0; j<match_constr_dim[lam]; j++) {
        jc = match_constr_idx[lam][j];
        for (int i=0; i<b_len; i++) {
          // cout << "i = " << i << endl;
          // poly_frac_build(&pf_sol[offset+i]);
          // cout << "hsol poly:" << endl;
          // print_poly(hsol[n][0][i][match_constr_idx[lam][j]], tp_rank+1);
          // cout << "set coeffs..." << endl;
          poly_frac_set_coeffs(&pf_sol[offset+i], hsol[n][0][i][match_constr_idx[lam][j]], tp_rank+1);
          pf_sol[offset+i].den_deg = 0;
          pf_sol[offset+i].nroots = nroots;
          pf_sol[offset+i].mults = new int[nroots];
          // cout << "set mults..." << endl;
          for (int k=0; k<nroots; k++) {
            pf_sol[offset+i].mults[k] = 0;
          }
          // cout << "poly_frac:" << endl;
          // poly_frac_print(&pf_sol[offset+i]);
        }

        // cout << "mul tmat..." << endl;
        // poly_frac_rk2_mul_pf_rk1(tmat_times_hsol[lam], tmat, pf_sol+offset, b_len, b_len, roots);
        for (int i=0; i<b_len; i++) {
          if (print) cout << "lam, i, match_constr_idx = " << lam << ", " << i << ", " << match_constr_idx[lam][j] << ", ";
          // cout << endl;
          poly_frac_mul_pf(
            &tmat_times_hsol[lam][i][j], &tmat[i][0], &pf_sol[offset]
          );
          for (int k=1; k<b_len; k++) {
            poly_frac_mul_pf(
              &pf_tmp, &tmat[i][k], &pf_sol[offset+k]
            );
            poly_frac_add_pf(
              &tmat_times_hsol[lam][i][j], &tmat_times_hsol[lam][i][j], &pf_tmp,
              roots,
              NULL
            );
          }
          poly_frac_prune_rel_tol(&tmat_times_hsol[lam][i][j], wp_bin);
          pow_behav[lam][i][j] = poly_frac_pow_behav(&tmat_times_hsol[lam][i][j]);
          if (print) {
          cout << "pow behav = " << pow_behav[lam][i][j] << endl;
          poly_frac_print(&tmat_times_hsol[lam][i][j]);
          }
          if (pow_behav[lam][i][j] < pow_behav_min[lam][i]) {
            pow_behav_min[lam][i] = pow_behav[lam][i][j];
          }
          // cout << "poly_frac[" << i << "][" << j << "]:" << endl;
          // poly_frac_print(&tmat_times_hsol[lam][i][j]);
        }
      }
    }

    for (int lam, n=0; n<num_cum_eig; n++) {
      lam = cum_eig[n];
      if (block_log_len[lam] == 0) {
        continue;
      }
      // for (int i=0; i<b_len; i++) {
      //   cout << "lam, i = " << lam << ", " << i << ", pow behav min = " << pow_behav_min[lam][i] << endl;
      // }
    }

    // #DEPRECATED
    // it was used to predict pow_behav in order without computing poly_frac
    // // DECIDE TYPE OF MATCHING
    // cout << "DECIDE TYPE OF MATCHING" << endl;
    // for (int i=0; i<b_len; i++) {
    //   for (int lam, n=0; n<num_cum_eig; n++) {
    //     lam = cum_eig[n];
    //     if (block_log_len[lam] == 0) {
    //       continue;
    //     }
    //     pow_behav_min[lam][i] = tp_rank+1;
    //     for (int j=0; j<match_constr_dim[lam]; j++) {
    //       cout << "lam, i, j = " << lam << ", " << i << ", " << j << ", match_constr_idx = " << match_constr_idx[lam][j] << endl;
    //       // cout << "i = " << i << endl;
    //       // get ldeg of the solution
    //       for (int k=0; k<b_len; k++) {
    //         ldeg = 0;
    //         for (int ld=0; ld < tp_rank+1; ld++) {
    //           if (mpc_lessthan_tol(hsol[n][0][k][match_constr_idx[lam][j]][ld])) {
    //             ldeg++;
    //           } else {
    //             break;
    //           }
    //         }
    //         // cout << "ldeg = " << ldeg << endl;

    //         // get power behaviour
    //         tmp_pow_behav = ldeg + poly_frac_pow_behav(&tmat[i][k]);
    //         if (tmp_pow_behav < ) {}
    //         cout << "pow behav = " << pow_behav[lam][i][j] << endl;
    //         if (pow_behav[lam][i][j] < pow_behav_min[lam][i]) {
    //           pow_behav_min[lam][i] = pow_behav[lam][i][j];
    //         }
    //       }
    //     }
    //     cout << "pow behav min = " << pow_behav_min[lam][i] << endl;
    //   }
    // }
    
    // particular
    // int *part_pow_behav;
    // if (sb_len>0) {
    //   cout << "particular" << endl;
    //   for (int lamp, np=0; np<num_cum_eig; np++) {
    //     lamp = cum_eig[np];
    //     cout << "lamp = " << lamp << endl;
    //     if (sol_log_len[lamp] == 0) {
    //       continue;
    //     }
    //     part_pow_behav[lamp] = tp_rank+1;
    //     for (int i=0; i<b_len; i++) {
    //       // get ldeg of the solution
    //       ldeg = 0;
    //       for (int ld=0; ld < tp_rank+1; ld++) {
    //         if (mpc_lessthan_tol(psol[np][0][i][0][ld])) {
    //           ldeg++;
    //         } else {
    //           break;
    //         }
    //       }

    //       // get power behaviour
    //       for (int k=0; k<b_len; k++) {
    //         tmp_pow_behav = ldeg + tmat[i][k].num_deg - tmat[i][k].num_vdeg - tmat[i][k].mults[0];
    //         if (tmp_pow_behav < part_pow_behav[lamp]) {
    //           part_pow_behav[lamp] = tmp_pow_behav;
    //         }
    //       }
    //     }
    //   }
    // } 

    if (print) cout << endl << "DETECT NULL COEFFICIENTS" << endl;
    int *coeff_is_null = new int[b_len];  // array that is one when a column gets a null coefficient
    for (int j=0; j<b_len; j++) {
      coeff_is_null[j] = 0;
    }
    int num_constr_max = 0, num_null_coeff = 0;
    int *is_in_bound = new int[num_classes];
    // malloc_rk2_tens(is_in_bound, b_len, num_classes);
    for (int lam, n=0; n<num_cum_eig; n++) {
      lam = cum_eig[n];
      if (print) cout << "lam = " << lam << endl;
      if (block_log_len[lam] == 0) {
        continue;
      }
      if (print) cout << "is in block" << endl;

      // check whether eigenvalue is present in the boundary
      is_in_bound[lam] = 0;
      for (int i=0; i<b_len; i++) {
        for (int lamb, nb=0; nb<mi_eig_num[offset+i]; nb++) {
          lamb = mi_eig[offset+i][nb];
          if (lamb == lam) {
            // cout << "i, lam = " << i << ", " << lam << " found in boundary" << endl;
            is_in_bound[lam] = 1;
            break;
          }
        }
        if (is_in_bound[lam]) {
          break;
        }
      }

      if (is_in_bound[lam]) {
        for (int i=0; i<b_len; i++) {
          if (print) {
          cout << "lam, i = " << lam << ", " << i << endl;
          cout << "bound behav = " << bound_behav[offset+i][lam] << endl;
          cout << "min pow behav = " << pow_behav_min[lam][i] << endl;
          }
          num_constr_max += bound_behav[offset+i][lam] - pow_behav_min[lam][i] + 1;
        }
      } else {
        num_null_coeff += match_constr_dim[lam];
        for (int jc, j=0; j<match_constr_dim[lam]; j++) {
          jc = match_constr_idx[lam][j];
          coeff_is_null[jc] = 1;
        }
      }
    }
    int num_constr = b_len - num_null_coeff;
    if (print) {
    cout << "max num constraints = " << num_constr_max << endl;
    cout << "num null coeff = " << num_null_coeff << endl;
    for (int j=0; j<b_len; j++) {
      cout << coeff_is_null[j] << ", ";
    }
    cout << endl;
    cout << "num constr needed = " << num_constr << endl;
    }

    if (num_constr == 0) {
      if (print) cout << "skip to linear combination" << endl;
      // shift back
      for (int i=0; i<b_len; i++) {
        tmat[i] -= offset;
      }
      tmat -= offset;
      goto goto_lin_comb;
    }

    mpc_t **constr, **constr_inv, *constr_const_term;
    malloc_rk2_tens(constr, num_constr_max, b_len);
    init_rk2_mpc(constr, num_constr_max, b_len);
    constr_const_term = new mpc_t[b_len];
    init_rk1_mpc(constr_const_term, b_len);

    mpc_t **ext_constr;
    malloc_rk2_tens(ext_constr, b_len, b_len);
    init_rk2_mpc(ext_constr, b_len, b_len);
    ext_constr_const_term = new mpc_t[b_len];
    init_rk1_mpc(ext_constr_const_term, b_len);

    // int *match_constr_dim = new int[num_classes];
    int **constr_idx;
    malloc_rk2_tens(constr_idx, num_classes, b_len);
    // mpc_t **match_constr_mat, **match_constr_inv_mat, *match_constr_coeffs;

    // count and select indices for every eigenvalue
    int count;
    for (int lam, n=0; n<num_cum_eig; n++) {
      lam = cum_eig[n];
      // cout << "lam = " << lam << endl;
      if (block_log_len[lam] == 0) {
        continue;
      }

      // select column indices associated with the current eigenvalue
      count = 0;
      for (int j=0; j<b_len; j++) {
        // cout << "j = " << j << endl;
        // cout << "eq_class = " << eq_class[offset+j] << endl;
        if (eq_class[offset+j] == lam) {
          // cout << "picked" << endl;
          constr_idx[lam][count] = j;
        }
      }
    }

    if (print) cout << "DETECT FREE CONSTRAINTS" << endl;
    free_constr = 0;

    // // #hard-coded for the 2L-box
    // if (offset == 6) {
    //   mpc_set_str(
    //     match_constr_coeffs[0],
    //     (const char*)"(-1.744113385330068957560639044431000777372935544224743942542009046206541688046885798163793237827964842923234640357906751061987850843440195465327521119355869e15 -3.406958109513536047789675364584630841127871347849646315211023721098627020308516820876605743125001374661590306601249493712015905847096925950901368927226e11)",
    //     10,
    //     MPFR_RNDN
    //   );
    //   cout << "hard-coded coeff:" << endl; print_mpc(&match_constr_coeffs[0]); cout << endl;
    //   goto goto_hardcoded_2L_box;
    // }

    int **free_constr_idx;
    malloc_rk2_tens(free_constr_idx, num_constr_max, b_len);

    if (print) cout << endl << "first detect additional null coefficients" << endl;
    for (int i=0; i<b_len; i++) {
      for (int lam, n=0; n<num_cum_eig; n++) {
        lam = cum_eig[n];
        if (block_log_len[lam] == 0) {
          continue;
        }
        if (print) {
        cout << "i, lam = " << i << ", " << lam << endl;
        cout << "is in bound = " << is_in_bound[lam] << endl;
        }
        if (is_in_bound[lam] == 0) {
          continue;
        }
        if (print) {
        cout << "pow_behav_min = " << pow_behav_min[lam][i] << endl;
        cout << "bound_behav = " << bound_behav[offset+i][lam] << endl;
        cout << "const_term_pow_behav = " << const_term_pow_behav[lam][i] << endl;
        }
        // if (offset == 3) {exit(0);}
        if (
          pow_behav_min[lam][i] < bound_behav[offset+i][lam] &&\
          pow_behav_min[lam][i] < const_term_pow_behav[lam][i]
        ) {
          int num_non_zero = 0;
          for (int jc, j=0; j<match_constr_dim[lam]; j++) {
            jc = match_constr_idx[lam][j];
            if (coeff_is_null[jc] == 1) {continue;}
            if (poly_frac_pow_behav(&tmat_times_hsol[lam][i][j]) == pow_behav_min[lam][i]) {
              num_non_zero++;
            }
          }
          if (num_non_zero == 1) {
            if (print) cout << "found new null coeff" << endl;
            num_null_coeff++;
            for (int jc, j=0; j<match_constr_dim[lam]; j++) {
              jc = match_constr_idx[lam][j];
              if (poly_frac_pow_behav(&tmat_times_hsol[lam][i][j]) == pow_behav_min[lam][i]) {
                coeff_is_null[jc] = 1;
                break;
              }
            }
          }
        }
      }
    }

    // detect null homogeneous columns
    if (print) {cout << endl; cout << "detect null homogeneous columns" << endl;}
    int null_col;
    for (int lam, n=0; n<num_cum_eig; n++) {
      lam = cum_eig[n];
      if (block_log_len[lam] == 0) {
        continue;
      }

      for (int jc, j=0; j<match_constr_dim[lam]; j++) {
        jc = match_constr_idx[lam][j];
        if (print) cout << "lam, j, jc = " << lam << ", " << j << ", " << jc << endl;
        null_col = 1;
        for (int i=0; i<b_len; i++) {
          if (tmat_times_hsol[lam][i][j].num_vdeg != -1) {
            null_col = 0;
            break;
          }
        }

        if (null_col) {
          if (print) cout << "found new null coeff" << endl;
          num_null_coeff++;
          coeff_is_null[jc] = 1;
        }
      }
    }

    num_constr = b_len - num_null_coeff;
    if (print) {
    cout << "max num constraints = " << num_constr_max << endl;
    cout << "num null coeff = " << num_null_coeff << endl;
    for (int j=0; j<b_len; j++) {
      cout << coeff_is_null[j] << ", ";
    }
    cout << endl;
    cout << "num constr needed = " << num_constr << endl;
    }

    if (num_constr == 0) {
      if (print) cout << "skip to linear combination" << endl;
      // shift back
      for (int i=0; i<b_len; i++) {
        tmat[i] -= offset;
      }
      tmat -= offset;

      goto goto_lin_comb;
    }
    // if (offset == 3) {exit(0);}
    
    // if (0 && offset > 0 && sb_len > 0) {
    //   cout << endl; cout << "detect even more null coefficients from the DE" << endl;
    //   // subtract lam/eta from block
    //   cout << "subtract lam/eta from block" << endl;
    //   for (int lam, n=0; n<num_cum_eig; n++) {
    //     lam = cum_eig[n];

    //     if (mpc_lessthan_tol(eig_list[lam])) {
    //       for (int i=0; i<b_len; i++) {
    //         for (int j=0; j<b_len; j++) {
    //           poly_frac_set_pf(&block_m_lam[lam][i][j], &pfmat[offset+i][offset+j]);
    //         }
    //       }
    //       continue;
    //     }

    //     poly_frac_set_mpc(&lam_pf[lam], &eig_list[lam], nroots);
    //     poly_frac_mul_sym_pow(&lam_pf[lam], &lam_pf[lam], -1);
    //     poly_frac_neg(&lam_pf[lam]);

    //     for (int i=0; i<b_len; i++) {
    //       for (int j=0; j<b_len; j++) {
    //         if (i == j && block_log_len[lam] != 0) {
    //           poly_frac_add_pf(
    //             &block_m_lam[lam][i][j],
    //             &pfmat[offset+i][offset+j], &lam_pf[lam],
    //             roots, NULL
    //           );
    //         } else {
    //           poly_frac_set_pf(&block_m_lam[lam][i][j], &pfmat[offset+i][offset+j]);
    //         }
    //       }
    //     }
    //   }
    //   // if (offset == 1) {exit(0);}

    //   // CONTRIBUTION FROM PREV UN-TRANSFORMED SOL      
    //   // if (offset > 0) {
    //   if (1) {
    //     cout << "CONTRIBUTION FROM PREV UN-TRANSFORMED SOL" << endl;
    //     for (int lam, n=0; n<num_cum_eig; n++) {
    //       lam = cum_eig[n];
    //       if (sol_log_len[lam] == 0) {
    //         continue;
    //       }
    //       for (int i=0; i<offset; i++) {
    //         // cout << "lam, i = " << lam << ", " << i << endl;
    //         // cout << "log prof = " << log_prof[lam][i] << endl;
    //         poly_frac_build(&pf_sol[i]);
    //         if (log_prof[lam][i] == 0) {
    //           // cout << "set pf to zero" << endl;
    //           poly_frac_set_zero(&pf_sol[i]);
    //         } else {
    //           // cout << "set pf coeffs from poly" << endl;
    //           // print_poly(prev_sol[lam][0][i], tp_rank+1);
    //           poly_frac_set_coeffs(&pf_sol[i], sol[lam][0][i], 1);
    //           pf_sol[i].den_deg = 0;
    //           pf_sol[i].nroots = nroots;
    //           pf_sol[i].mults = new int[nroots];
    //           for (int k=0; k<nroots; k++) {
    //             pf_sol[i].mults[k] = 0;
    //           }
    //         }
    //         poly_frac_neg(&pf_sol[i]);
    //         // poly_frac_print(&pf_sol[i]);
    //       }
    //       // cout << "mul DE" << endl;
    //       // poly_frac_rk2_print(tmat, b_len, offset);
    //       poly_frac_rk2_mul_pf_rk1(DE_const_term[lam], pfmat+offset, pf_sol, b_len, offset, roots);

    //       // for (int i=0; i<b_len; i++) {
    //       //   cout << "lam, i = " << lam << ", " << i << endl;
    //       //   poly_frac_normal(&pf_const_term[lam][i]);
    //       //   poly_frac_print(&pf_const_term[lam][i]);
    //       //   const_term_pow_behav[lam][i] = poly_frac_pow_behav(&pf_const_term[lam][i]);
    //       //   cout << "behav = " << const_term_pow_behav[lam][i] << endl;
    //       // }
    //     }
    //   }

    //   pfmat += offset;
    //   for (int i=0; i<b_len; i++) {
    //     pfmat[i] += offset;
    //   }

    //   // CONTRIBUTION FROM PARTICULAR SOLUTION      
    //   // if (sb_len>0) {
    //   if (1) {
    //     cout << "CONTRIBUTION FROM PARTICULAR SOLUTION" << endl;
    //     cout << "diag mat times particular" << endl;
    //     for (int lamp, np=0; np<num_cum_eig; np++) {
    //       lamp = cum_eig[np];
    //       // cout << "lamp = " << lamp << endl;
    //       if (sol_log_len[lamp] == 0) {
    //         continue;
    //       }
          
    //       // cout << "mul mat" << endl;
    //       // poly_frac_rk2_print(tmat, b_len, offset);
    //       poly_frac_rk2_mul_pf_rk1(mat_times_tmat_times_psol[lamp], block_m_lam[lamp], pf_const_term[lamp], b_len, b_len, roots);

    //       // cumulate
    //       // if (sb_len > 0) {
    //       if (1) {
    //         for (int i=0; i<b_len; i++) {
    //           poly_frac_add_pf(
    //             &DE_const_term[lamp][i], &DE_const_term[lamp][i], &mat_times_tmat_times_psol[lamp][i],
    //             roots, NULL
    //           );
    //         }
    //       } else {
    //         for (int i=0; i<b_len; i++) {
    //           poly_frac_set_pf(&DE_const_term[lamp][i], &mat_times_tmat_times_psol[lamp][i]);
    //         }
    //       }
    //     }
    //   }

    //   // get power behavior
    //   cout << "get power behavior" << endl;
    //   for (int lamp, np=0; np<num_cum_eig; np++) {
    //     lamp = cum_eig[np];
    //     // cout << "lamp = " << lamp << endl;
    //     if (sol_log_len[lamp] == 0) {
    //       continue;
    //     }
        
    //     for (int i=0; i<b_len; i++) {
    //       DE_const_term_pow_behav[lamp][i] = poly_frac_pow_behav(&DE_const_term[lamp][i]);
    //     }
    //   }

    //   // CONTRIBUTION FROM HOMOGENEOUS SOLUTION
    //   cout << "CONTRIBUTION FROM HOMOGENEOUS SOLUTION" << endl;
    //   for (int lam, n=0; n<num_cum_eig; n++) {
    //     lam = cum_eig[n];
    //     if (block_log_len[lam] == 0) {
    //       continue;
    //     }

    //     for (int jc, j=0; j<match_constr_dim[lam]; j++) {
    //       jc = match_constr_idx[lam][j];

    //       // cout << "mul mat..." << endl;
    //       // poly_frac_rk2_mul_pf_rk1(tmat_times_hsol[lam], tmat, pf_sol+offset, b_len, b_len, roots);
    //       for (int i=0; i<b_len; i++) {
    //         // cout << "lam, i, match_constr_idx = " << lam << ", " << i << ", " << match_constr_idx[lam][j] << ", ";
    //         poly_frac_mul_pf(
    //           &mat_times_tmat_times_hsol[lam][i][j], &block_m_lam[lam][i][0], &tmat_times_hsol[lam][0][j]
    //         );

    //         for (int k=1; k<b_len; k++) {
    //           poly_frac_mul_pf(
    //             &pf_tmp, &block_m_lam[lam][i][k], &tmat_times_hsol[lam][k][j]
    //           );
    //           poly_frac_add_pf(
    //             &mat_times_tmat_times_hsol[lam][i][j], &mat_times_tmat_times_hsol[lam][i][j], &pf_tmp,
    //             roots,
    //             NULL
    //           );
    //         }

    //         tmp_pow_behav = poly_frac_pow_behav(&mat_times_tmat_times_hsol[lam][i][j]);
    //         pow_behav[lam][i][j] = poly_frac_pow_behav(&tmat_times_hsol[lam][i][j]);
    //         // cout << "pow behav = " << pow_behav[lam][i][j] << endl;
    //         if (tmp_pow_behav< DE_hom_pow_behav[lam][i]) {
    //           DE_hom_pow_behav[lam][i] = tmp_pow_behav;
    //         }
    //         // cout << "poly_frac[" << i << "][" << j << "]:" << endl;
    //         // poly_frac_print(&tmat_times_hsol[lam][i][j]);
    //       }
    //     }
    //   }

    //   // detect null coefficients
    //   int thresh_pow;
    //   for (int lam, n=0; n<num_cum_eig; n++) {
    //     lam = cum_eig[n];
    //     if (block_log_len[lam] == 0) {
    //       continue;
    //     }
        
    //     if (mpc_lessthan_tol(eig_list[lam])) {
    //       thresh_pow = 0;
    //     } else {
    //       thresh_pow = 1;
    //     }

    //     for (int i=0; i<b_len; i++) {
    //       cout << "lam, i = " << lam << ", " << i << endl;
    //       cout << "DE hom behav, DE const behav = " << DE_hom_pow_behav[lam][i] << ", " << DE_const_term_pow_behav[lam][i] << endl;
    //       if (DE_hom_pow_behav[lam][i] < DE_const_term_pow_behav[lam][i] && DE_hom_pow_behav[lam][i] < thresh_pow) {
    //         int num_non_zero = 0;
    //         for (int jc, j=0; j<match_constr_dim[lam]; j++) {
    //           jc = match_constr_idx[lam][j];
    //           if (coeff_is_null[jc] == 1) {continue;}
    //           // cout << "lam, i, match_constr_idx = " << lam << ", " << i << ", " << jc << ": ";
    //           if (poly_frac_pow_behav(&mat_times_tmat_times_hsol[lam][i][j]) == DE_hom_pow_behav[lam][i]) {
    //             num_non_zero++;
    //           }
    //         }
    //         if (num_non_zero == 1) {
    //           cout << "found new null coeff" << endl;
    //           num_null_coeff++;
    //           for (int jc, j=0; j<match_constr_dim[lam]; j++) {
    //             jc = match_constr_idx[lam][j];
    //             if (poly_frac_pow_behav(&tmat_times_hsol[lam][i][j]) == pow_behav_min[lam][i]) {
    //               coeff_is_null[jc] = 1;
    //               break;
    //             }
    //           }
    //         }
    //       }
    //     }
    //   }

    //   for (int i=0; i<b_len; i++) {
    //     pfmat[i] -= offset;
    //   }
    //   pfmat -= offset;

    //   num_constr = b_len - num_null_coeff;
    //   cout << "max num constraints = " << num_constr_max << endl;
    //   cout << "num null coeff = " << num_null_coeff << endl;
    //   for (int j=0; j<b_len; j++) {
    //     cout << coeff_is_null[j] << ", ";
    //   }
    //   cout << endl;
    //   cout << "num constr needed = " << num_constr << endl;

    //   if (num_constr == 0) {
    //     cout << "skip to linear combination" << endl;
    //     // shift back
    //     for (int i=0; i<b_len; i++) {
    //       tmat[i] -= offset;
    //     }
    //     tmat -= offset;
    //     goto goto_lin_comb;
    //   }

    // }

    if (print) cout << endl << "look for residual free constraints" << endl;
    // update min power behav
    // (after setting to zero some columns, the minimun might increase)    
    int first;
    for (int lam, n=0; n<num_cum_eig; n++) {
      lam = cum_eig[n];
      if (block_log_len[lam] == 0) {
        continue;
      }      
      for (int i=0; i<b_len; i++) {
        first = 1;
        for (int jc, j=0; j<match_constr_dim[lam]; j++) {
          jc = match_constr_idx[lam][j];
          if (coeff_is_null[jc] == 1) {
            continue;
          }
          if (first) {
            pow_behav_min[lam][i] = pow_behav[lam][i][j];
            first = 0;
          }
          if (pow_behav[lam][i][j] < pow_behav_min[lam][i]) {
            pow_behav_min[lam][i] = pow_behav[lam][i][j];
          }
        }
      }
    }


    // for (int i=0; i<num_constr_max; i++) {
    //   for (int j=0; j<b_len; j++) {
    //     free_constr_idx[i][j] = 0;
    //   }
    // }
    free_constr_const_term = new mpc_t[b_len];
    init_rk1_mpc(free_constr_const_term, b_len);
    for (int i=0; i<b_len; i++) {
      for (int lam, n=0; n<num_cum_eig; n++) {
        lam = cum_eig[n];
        if (print) {
        cout << "i, lam = " << i << ", " << lam << endl;
        cout << "block log len = " << block_log_len[lam] << endl;
        cout << "is in bound = " << is_in_bound[lam] << endl;
        }
        if (block_log_len[lam] == 0) {
          continue;
        }
        if (is_in_bound[lam] == 0) {
          continue;
        }
        if (print) {
        cout << "pow_behav_min = " << pow_behav_min[lam][i] << endl;
        cout << "bound_behav = " << bound_behav[offset+i][lam] << endl;
        }
        int LO_selected = 0;
        for (int p=pow_behav_min[lam][i]; p<bound_behav[offset+i][lam]; p++) {
          int found_non_zero = 0;
          for (int jc, j=0; j<match_constr_dim[lam]; j++) {
            jc = match_constr_idx[lam][j];
            if (print) {
            cout << "p, jc = " << p << ", " << jc << endl;
            cout << "part free constr = " << free_constr << endl;
            cout << "tmat times hsol:" << endl;
            poly_frac_print(&tmat_times_hsol[lam][i][j]);
            cout << "pow behav = " << poly_frac_pow_behav(&tmat_times_hsol[lam][i][j]) << endl;
            }
            if (coeff_is_null[jc] == 1) {
              if (print) cout << "coefficient is null" << endl;
              continue;
            }
            //////
            // EXTRACT HOMOGENEOUS COEFFICIENTS ASSOCIATED WITH SELECTED POWER BEHAVIOR
            //////
            if (poly_frac_pow_behav(&tmat_times_hsol[lam][i][j]) >= 0) {
              mpc_set_ui(constr[free_constr][jc], 0, MPFR_RNDN);
            } else if (p == poly_frac_pow_behav(&tmat_times_hsol[lam][i][j])) {
              // poly_frac_extract_LO(
              rel_err_poly_frac_extract_LO(
                &mpc_realref(constr[free_constr][jc]),
                &mpc_imagref(constr[free_constr][jc]),
                &tmat_times_hsol[lam][i][j], roots, -p
                , wp_bin
              );
              LO_selected = 1;
            } else if (p == poly_frac_pow_behav(&tmat_times_hsol[lam][i][j]) + 1) {
              // poly_frac_extract_NLO(
              rel_err_poly_frac_extract_NLO(
                &constr[free_constr][jc],
                &constr[free_constr-LO_selected][jc],
                &tmat_times_hsol[lam][i][j], roots, -p+1
                , wp_bin
              );
            } else if (p == poly_frac_pow_behav(&tmat_times_hsol[lam][i][j]) + 2) {
              // curretnly only LO and NLO extraction has been implemented;
              // extraction of NNLO and beyond of a poly_frac still not available
              printf("WARNING: trying to extract NNNLO behavior of a poly_frac while computing free constarints\n");
              perror("WARNING: trying to extract NNNLO behavior of a poly_frac while computing free constarints");
              exit(1);
            }
            // poly_frac_mul_sym_pow(&tmp_pf[free_constr], &tmat_times_hsol[lam][i][j], -p);
            // poly_frac_eval_sym_zero(&constr[free_constr][j], &tmp_pf[i], roots);
            // poly_frac_print(&tmp_pf[free_constr]);
            if (print) {
            print_mpc(&constr[free_constr][jc]); cout << endl;
            }
            // free_constr_idx[free_constr][jc] == 1;
            if (!mpc_lessthan_tol(constr[free_constr][jc])) {
              found_non_zero = 1;
            }
            // if (offset == 13) { continue;}
          }
          if (found_non_zero) { // perform csvd instead
            //////
            // EXTRACT CONSTANT TERM COEFFICIENT ASSOCIATED WITH SELECTED POWER BEHAVIOR
            //////
            if (print) cout << "EXTRACT CONSTANT TERM COEFFICIENT ASSOCIATED WITH SELECTED POWER BEHAVIOR" << endl;
            if (print) poly_frac_print(&pf_const_term[lam][i]);
            if (poly_frac_pow_behav(&pf_const_term[lam][i]) == 0) {
              if (print) cout << "pow behav of const term = 0" << endl;
              mpc_set_ui(free_constr_const_term[free_constr], 0, MPFR_RNDN);
            } else if (p == poly_frac_pow_behav(&pf_const_term[lam][i])) {
              if (print) cout << "extract LO of the constant term" << endl;
              // poly_frac_extract_LO(
              rel_err_poly_frac_extract_LO(
                &mpc_realref(free_constr_const_term[free_constr]),
                &mpc_imagref(free_constr_const_term[free_constr]),
                &pf_const_term[lam][i], roots, -p
                , wp_bin
              );
              LO_selected = 1;
            } else if (p == poly_frac_pow_behav(&pf_const_term[lam][i]) + 1) {
              if (print) cout << "extract NLO of the constant term" << endl;
              // poly_frac_extract_NLO(
              rel_err_poly_frac_extract_NLO(
                &free_constr_const_term[free_constr],
                &free_constr_const_term[free_constr-LO_selected],
                &pf_const_term[lam][i], roots, -p+1
                , wp_bin
              );
            } else if (p == poly_frac_pow_behav(&pf_const_term[lam][i]) + 2) {
              // curretnly only LO and NLO extraction has been implemented;
              // extraction of NNLO and beyond of a poly_frac still not available
              printf("WARNING: trying to extract NNLO behavior of a poly_frac while computing free constarints\n");
              perror("WARNING: trying to extract NNLO behavior of a poly_frac while computing free constarints");
              exit(1);
            }
            mpc_neg(free_constr_const_term[free_constr], free_constr_const_term[free_constr], MPFR_RNDN);
            if (print) {
            print_mpc(&free_constr_const_term[free_constr]); cout << endl;
            }

            if (free_constr >= 1) {
              //////
              // CHECK CONSISTENCY
              //////
              // prepare input for csvd
              int row_dim = free_constr+1;
              int col_dim = b_len-num_null_coeff;            
              mpc_t *csvd_mat = new mpc_t[row_dim*(col_dim+1)];
              init_rk1_mpc(csvd_mat, row_dim*(col_dim+1));
              mpc_t *csvd_mat0 = new mpc_t[row_dim*(col_dim+1)];
              init_rk1_mpc(csvd_mat0, row_dim*(col_dim+1));
              mpc_t *csvd_mat1 = new mpc_t[row_dim*(col_dim+1)];
              init_rk1_mpc(csvd_mat1, row_dim*(col_dim+1));
              build_csvd_mat(
                csvd_mat,
                row_dim, col_dim, constr, free_constr_const_term,
                b_len, coeff_is_null
              );

              // mpc_set_ui(csvd_mat[0+row_dim*0], 1, MPFR_RNDN);
              // mpc_set_ui(csvd_mat[1+row_dim*0], 1, MPFR_RNDN);
              // mpc_set_ui(csvd_mat[0+row_dim*3], 2, MPFR_RNDN);
              // mpc_set_ui(csvd_mat[1+row_dim*3], 2, MPFR_RNDN);

              if (print) {
              cout << "CSVD mat:" << endl;
              for (int i=0; i<row_dim; i++) {
                for (int j=0; j<col_dim+1; j++) {
                  cout << "i, j = " << i << ", " << j << ": ";
                  print_mpc(&csvd_mat[i+row_dim*j]); cout << endl;
                }
                cout << endl;
              }
              }

              for (int i=0; i<row_dim; i++) {
                for (int j=0; j<col_dim; j++) {
                  mpc_set(csvd_mat0[j+col_dim*i], csvd_mat[i+row_dim*j], MPFR_RNDN);
                }
              }

              if (print) {
              cout << "csvd_mat0:" << endl;
              for (int i=0; i<col_dim; i++) {
                for (int j=0; j<row_dim; j++) {
                  cout << "i, j = " << i << ", " << j << ": ";
                  print_mpc(&csvd_mat0[i+col_dim*j]); cout << endl;
                }
                cout << endl;
              }
              }

              for (int i=0; i<row_dim; i++) {
                for (int j=0; j<col_dim+1; j++) {
                  mpc_set(csvd_mat1[j+(col_dim+1)*i], csvd_mat[i+row_dim*j], MPFR_RNDN);
                }
              }
              
              if (print) {
              cout << "csvd_mat1:" << endl;
              for (int i=0; i<col_dim+1; i++) {
                for (int j=0; j<row_dim; j++) {
                  cout << "i, j = " << i << ", " << j << ": ";
                  print_mpc(&csvd_mat1[i+(col_dim+1)*j]); cout << endl;
                }
                cout << endl;
              }
              }

              // apply Rouché-Capelli theorem
              int rank0, rank1;
              // if (offset == 34) wp2 *= 0.90;
              // if (offset == 34) mpfr_tol_enlarge(0.80);
              rank0 = compute_rank(csvd_mat0, col_dim, row_dim);
              if (print) cout << "rank0 = " << rank0 << endl;
              rank1 = compute_rank(csvd_mat1, col_dim+1, row_dim);
              if (print) cout << "rank1 = " << rank1 << endl;
              if (rank0 != rank1) {
                printf("WARNING: found inconsistency while matching in zero: constraints are contradictory\n");
                perror("WARNING: found inconsistency while matching in zero: constraints are contradictory");
                LO_selected = 0;
              } else if (rank0 != free_constr+1) {
                // constraint linearly depedends on the previous ones
                LO_selected = 0;
              } else {
                // constraint does not linearly depedend on the previous ones
                free_constr++;
              }
            } else {
              free_constr++;
            }
          } else {
            LO_selected = 0;
          }
          if (free_constr == num_constr) {
            break;
          }
        }
        if (free_constr == num_constr) {
          break;
        }
      }
      if (free_constr == num_constr) {
        break;
      }
    }
    if (print) {
    cout << "free constr = " << free_constr << endl;
    cout << "free constr matrix:" << endl;
    for (int i=0; i<free_constr; i++) {
      for (int j=0; j<b_len; j++) {
        if (coeff_is_null[j] == 1) {continue;}
        cout << "i, j = " << i << ", " << j << ": ";
        print_mpc(&constr[i][j]); cout << endl;
      }
    }
    cout << "free constr const term:" << endl;
    for (int i=0; i<free_constr; i++) {
      cout << "i = " << i << ": ";
      print_mpc(&free_constr_const_term[i]); cout << endl;
    }
    }

    // if (0 && offset > 0 && sb_len > 0) {
    // // if (0) {
    //   cout << endl; cout << "look for even more residual free constraints from the DE" << endl;
    //   int thresh_pow;
    //   for (int i=0; i<b_len; i++) {
    //     for (int lam, n=0; n<num_cum_eig; n++) {
    //       lam = cum_eig[n];
    //       cout << "i, lam = " << i << ", " << lam << endl;
    //       if (block_log_len[lam] == 0) {
    //         continue;
    //       }
    //       cout << "DE hom behav, DE const behav = " << DE_hom_pow_behav[lam][i] << ", " << DE_const_term_pow_behav[lam][i] << endl;

    //       if (
    //         mpc_lessthan_tol(eig_list[lam]) &&\
    //         DE_hom_pow_behav[lam][i] == DE_const_term_pow_behav[lam][i] &&\
    //         DE_hom_pow_behav[lam][i] == 0
    //       ) {
    //         cout << "missing implementation to apply DE contraint on lam = 0" << endl;;
    //         continue;
    //       }
          
    //       // SET POWERS BELOW THRESHOLD TO ZERO
    //       cout << "SET POWERS BELOW THRESHOLD TO ZERO" << endl;
    //       int LO_selected = 0;
    //       for (int p=DE_hom_pow_behav[lam][i]; p<DE_const_term_pow_behav[lam][i]; p++) {
    //         int found_non_zero = 0;
    //         for (int jc, j=0; j<match_constr_dim[lam]; j++) {
    //           jc = match_constr_idx[lam][j];
    //           cout << "p, jc = " << p << ", " << jc << endl;
    //           // poly_frac_print(&mat_times_tmat_times_hsol[lam][i][j]);
    //           cout << "pow behav = " << poly_frac_pow_behav(&mat_times_tmat_times_hsol[lam][i][j]) << endl;
    //           if (coeff_is_null[jc] == 1) {
    //             cout << "coefficient is null" << endl;
    //             continue;
    //           }
    //           //////
    //           // EXTRACT HOMOGENEOUS COEFFICIENTS ASSOCIATED WITH SELECTED POWER BEHAVIOR
    //           //////
    //           if (poly_frac_pow_behav(&tmat_times_hsol[lam][i][j]) == 0) {
    //             mpc_set_ui(constr[free_constr][jc], 0, MPFR_RNDN);
    //           } else if (p == poly_frac_pow_behav(&mat_times_tmat_times_hsol[lam][i][j])) {
    //             poly_frac_extract_LO(
    //               &mpc_realref(constr[free_constr][jc]),
    //               &mpc_imagref(constr[free_constr][jc]),
    //               &mat_times_tmat_times_hsol[lam][i][j], roots, -p
    //             );
    //             LO_selected = 1;
    //           } else if (p == poly_frac_pow_behav(&mat_times_tmat_times_hsol[lam][i][j]) + 1) {
    //             poly_frac_extract_LO(
    //               &mpc_realref(constr[free_constr][jc]),
    //               &mpc_imagref(constr[free_constr][jc]),
    //               &mat_times_tmat_times_hsol[lam][i][j], roots, -p
    //             );
    //           } else if (p == poly_frac_pow_behav(&mat_times_tmat_times_hsol[lam][i][j]) + 2) {
    //             // curretnly only LO and NLO extraction has been implemented;
    //             // extraction of NNLO and beyond of a poly_frac still not available
    //             perror("WARNING: trying to extract NNNLO behavior of a poly_frac while computing free constarints");
    //             exit(1);
    //           }
    //           cout << "constraint:" << endl;
    //           print_mpc(&constr[free_constr][jc]); cout << endl;
              
    //           if (mpc_lessthan_tol(constr[free_constr][jc]) == 0) {
    //             found_non_zero = 1;
    //           }
    //         }

    //         if (found_non_zero) { // perform csvd instead
    //           //////
    //           // EXTRACT CONSTANT TERM COEFFICIENT ASSOCIATED WITH SELECTED POWER BEHAVIOR
    //           //////
    //           if (poly_frac_pow_behav(&DE_const_term[lam][i]) == 0) {
    //             mpc_set_ui(free_constr_const_term[free_constr], 0, MPFR_RNDN);
    //           } else if (p == DE_const_term_pow_behav[lam][i]) {
    //             poly_frac_extract_LO(
    //               &mpc_realref(free_constr_const_term[free_constr]),
    //               &mpc_imagref(free_constr_const_term[free_constr]),
    //               &DE_const_term[lam][i], roots, -p
    //             );
    //             LO_selected = 1;
    //           } else if (p == DE_const_term_pow_behav[lam][i] + 1) {
    //             poly_frac_extract_NLO(
    //               &free_constr_const_term[free_constr],
    //               &free_constr_const_term[free_constr-LO_selected],
    //               &DE_const_term[lam][i], roots, -p
    //             );
    //           } else if (p == DE_const_term_pow_behav[lam][i] + 2) {
    //             // curretnly only LO and NLO extraction has been implemented;
    //             // extraction of NNLO and beyond of a poly_frac still not available
    //             perror("WARNING: trying to extract NNNLO behavior of a poly_frac while computing free constarints");
    //             exit(1);
    //           }
    //           free_constr++;
    //         } else {
    //           LO_selected = 0;
    //         }

    //         if (free_constr == num_constr) {
    //           break;
    //         }
    //       }
          
    //       if (free_constr == num_constr) {
    //         break;
    //       }

    //       // SET THRESHOLD POWER
    //       cout << "SET THRESHOLD POWER" << endl;
    //       int found_non_zero = 0;
    //       for (int jc, j=0; j<match_constr_dim[lam]; j++) {
    //         jc = match_constr_idx[lam][j];
    //         cout << "p, jc = " << -1 << ", " << jc << endl;
    //         cout << "pow behav = " << poly_frac_pow_behav(&mat_times_tmat_times_hsol[lam][i][j]) << endl;
    //         if (coeff_is_null[jc] == 1) {
    //           cout << "coefficient is null" << endl;
    //           continue;
    //         }
            
    //         poly_frac_extract_LO(
    //           &mpc_realref(constr[free_constr][jc]),
    //           &mpc_imagref(constr[free_constr][jc]),
    //           &mat_times_tmat_times_hsol[lam][i][j], roots, 1
    //         );
            
    //         cout << "constraint:" << endl;
    //         print_mpc(&constr[free_constr][jc]); cout << endl;
            
    //         if (!mpc_lessthan_tol(constr[free_constr][jc])) {
    //           found_non_zero = 1;
    //         }
    //       }
          
    //       if (found_non_zero) { // perform csvd instead
    //         poly_frac_extract_LO(
    //           &mpc_realref(free_constr_const_term[free_constr]),
    //           &mpc_imagref(free_constr_const_term[free_constr]),
    //           &DE_const_term[lam][i], roots, 1
    //         );
    //         free_constr++;
    //       }

    //       if (free_constr == num_constr) {
    //         break;
    //       }

    //     }
    //     if (free_constr == num_constr) {
    //       break;
    //     }
    //   }
    //   cout << "num null coeff = " << num_null_coeff << endl;
    //   cout << "free constr = " << free_constr << endl;
    //   cout << "free constr matrix:" << endl;
    //   for (int i=0; i<free_constr; i++) {
    //     for (int j=0; j<b_len; j++) {
    //       if (coeff_is_null[j] == 1) {continue;}
    //       cout << "i, j = " << i << ", " << j << ": ";
    //       print_mpc(&constr[i][j]); cout << endl;
    //     }
    //   }
    //   cout << "free constr const term:" << endl;
    //   for (int i=0; i<free_constr; i++) {
    //     cout << "i = " << i << ": ";
    //     print_mpc(&free_constr_const_term[i]); cout << endl;
    //   }
    // }

    if (free_constr == num_constr) {
      // if (constr_const_term) {
      //   delete[] constr_const_term;
      // }
      // constr_const_term = free_constr_const_term;
      // shift back
      for (int i=0; i<b_len; i++) {
        tmat[i] -= offset;
      }
      tmat -= offset;
      goto goto_solve;
    }
  
    // if (offset == 5) {exit(0);}
        // int num_constr = b_len, free_constr = 0;
        // if (const_term_pow_behav[lam][i] < bound_behav[offset+i][lam]) {
        //   if (pow_behav[lam][i] >= bound_behav[offset+i][lam]) {
        //     // stop and signal error
        //   }
        //   // boundary condition not needed
        //   free_constr++;
        // } else if (const_term_pow_behav[lam][i] == bound_behav[offset+i][lam]) {
        //   if (pow_behav[lam][i] >= bound_behav[offset+i][lam]) {
        //     // stop and signal that a boundary condition is needed
        //   }
        //   // coefficient must be zero
        //   num_constr--;
        // } else if (const_term_pow_behav[lam][i] > bound_behav[offset+i][lam]) {
        //   if (pow_behav[lam][i] != bound_behav[offset+i][lam]) {
        //     // stop and signal error
        //   }
        //   // boundary condition needed
        // }
    
    if (print) cout << "evaluate HOMOGENEOUS solution" << endl;
    count_mi_idx = 0;
    mi_constr_idx = new int[b_len];
    if (free_constr == 0 && num_null_coeff == 0) {
      if (print) cout << "no null coeffs, no free constr" << endl;
    // if (free_constr + num_null_coeff < num_constr_max) {
      for (int i=0; i<b_len; i++) {
        // cout << "i = " << i << endl;
        // cout << "class lam = " << eq_class[offset+i] << endl;
        // select only MI indices corrensponding to non-null external constraints
        // i.e. those actually inserted as input in boundary generation
        for (int lam, n=0; n<num_cum_eig; n++) {
          lam = cum_eig[n];
          if (block_log_len[lam] == 0) {
            continue;
          }
          for (int jc, j=0; j<match_constr_dim[lam]; j++) {
            jc = match_constr_idx[lam][j];
            // cout << "lam, i, jc = " << lam << ", " << i << ", " << jc << ", " << endl;
            // cout << "store in " << count_i << ", " << j << endl;
            // poly_frac_print(&tmat_times_hsol[lam][i][jc]);
            poly_frac_mul_sym_pow(&tmp_pf[i], &tmat_times_hsol[lam][i][j], -bound_behav[offset+i][lam]);
            // cout << "pf to eval in zero:" << endl;
            // poly_frac_print(&tmp_pf[i]);
            // poly_frac_eval_sym_zero(&constr[i][match_constr_idx[lam][jc]], &tmp_pf[i], roots);
            poly_frac_eval_sym_zero(&constr[i][jc], &tmp_pf[i], roots);
          }
        }
      }
    } else if (0 < free_constr+num_null_coeff && free_constr+num_null_coeff < b_len) {
      // choose which MI index has to be used as constraint
      if (print) cout << "!!! choose which MI index has to be used as constraint" << endl;
      count_mi_idx = 0;
      // mi_constr_idx = new int[b_len];
      for (int i=0; i<b_len; i++) {
        if (print) {
        cout << "i = " << i << endl;
        cout << "external boundary:" << endl;
        print_mpc(&(*solutions)[offset+i][0]); cout << endl;
        }
        if (mpc_cmp_si((*solutions)[offset+i][0], 0) == 0) {
          continue;
        }
        if (print) cout << "mi_constr_idx = " << i << endl;
        mi_constr_idx[count_mi_idx] = i;
        for (int lam, n=0; n<num_cum_eig; n++) {
          lam = cum_eig[n];
          if (block_log_len[lam] == 0) {
            continue;
          }
          if (lam != lam0) {  // #suspicous: might not hold when more than one lambda is present in the boundary!
            continue;
          }
          for (int jc, j=0; j<match_constr_dim[lam]; j++) {
            jc = match_constr_idx[lam][j];
            if (coeff_is_null[jc] == 1) {continue;}
            if (print) {
            cout << "jc = " << jc << endl;
            cout << "bound_behav = " << bound_behav[offset+i][lam] << endl;
            cout << "store in " << free_constr+count_mi_idx << ", " << jc << endl;
            }
            // poly_frac_mul_sym_pow(&tmp_pf[i], &tmat_times_hsol[lam][i][j], -bound_behav[offset+i][lam]);
            // poly_frac_eval_sym_zero(&constr[free_constr+count_mi_idx][jc], &tmp_pf[i], roots);
            // poly_frac_extract_oder(
            rel_err_poly_frac_extract_oder(
              &constr[free_constr+count_mi_idx][jc],
              &tmat_times_hsol[lam][i][j], roots, nroots, -bound_behav[offset+i][lam]
              , wp_bin
            );
            if (print) {
            cout << "tmat times hsol:" << endl;
            poly_frac_print(&tmat_times_hsol[lam][i][j]);
            cout << "constr:" << endl;
            print_mpc(&constr[free_constr+count_mi_idx][jc]); cout << endl;
            }
          }
        }

        // check whether new constraint is independent
        // cout << "check whether new constraint is independent" << endl;
        int constr_found_non_zero = 0;
        for (int lam, n=0; n<num_cum_eig; n++) {
          lam = cum_eig[n];
          if (block_log_len[lam] == 0) {
            continue;
          }
          if (lam != lam0) {  // #suspicous: might not hold when more than one lambda is present in the boundary!
            continue;
          }
          for (int jc, j=0; j<match_constr_dim[lam]; j++) {
            jc = match_constr_idx[lam][j];
            if (coeff_is_null[jc] == 1) {continue;}
            if (!mpc_lessthan_tol(constr[free_constr+count_mi_idx][jc])) {
              constr_found_non_zero = 1;
              break;
            }
          }
        }
        
        if (constr_found_non_zero) {
          count_mi_idx++;
        } else {
          continue;
        }
        if (free_constr + count_mi_idx == num_constr) {
          break;
        }
      }
    }

    // // if (0 && block_log_len[lam0] > 0) {
    // if(1) {
    // // evaluate homogeneous solution
    // cout << "EVALUATE HOMOGENOUS SOLUTION IN ZERO" << endl;
    // // mpc_t *vec = new mpc_t[b_len];
    // // init_rk1_mpc(vec, b_len);
    // // #pair pf_limit_in_zero
    //   for (int lam, n=0; n<num_cum_eig; n++) {
    //     lam = cum_eig[n];
    //     cout << "lam = " << lam << endl;
    //     // if (lam != lam0) {
    //     //   continue;
    //     // }
    //     if (block_log_len[lam] == 0) {
    //       continue;
    //     }
    //     for (int j=0; j<match_constr_dim[lam]; j++) {
    //       cout << "j = " << j << endl;
    //       cout << "match_constr_idx = " << match_constr_idx[lam][j] << endl;
    //       cout << "pow behav = " << pow_behav[lam][j] << endl;
    //       if (pow_behav[lam][j] >= 0 || sol_log_len[lam0] > 0) {
    //         continue;
    //       }
    //       for (int i=0; i<b_len; i++) {
    //         cout << "i = " << i << endl;
    //         poly_frac_build(&pf_sol[offset+i]);
    //         print_poly(hsol[n][0][i][match_constr_idx[lam][j]], tp_rank+1);
    //         cout << "set coeffs..." << endl;
    //         poly_frac_set_coeffs(&pf_sol[offset+i], hsol[n][0][i][match_constr_idx[lam][j]], tp_rank+1);
    //         pf_sol[offset+i].den_deg = 0;
    //         pf_sol[offset+i].nroots = nroots;
    //         pf_sol[offset+i].mults = new int[nroots];
    //         cout << "set mults..." << endl;
    //         for (int k=0; k<nroots; k++) {
    //           pf_sol[offset+i].mults[k] = 0;
    //         }
    //       }
    //       cout << "mul tmat..." << endl;
    //       poly_frac_rk2_mul_pf_rk1(tmp_pf, tmat, pf_sol+offset, b_len, b_len, roots);
    //       cout << "eval..." << endl;
    //       for (int i=0; i<b_len; i++) {
    //         // poly_frac_add_pf(
    //         //   &tmp_pf[i],
    //         //   &tmp_pf[i],
    //         //   &tmat_times_sol[i],
    //         //   roots, NULL
    //         // );
    //         cout << "poly_frac to be evaluated:" << endl;
    //         poly_frac_print(&tmp_pf[i]);
    //         poly_frac_mul_sym_pow(&tmp_pf[i], &tmp_pf[i], -pow_behav[lam][j]);
    //         poly_frac_eval_sym_zero(&match_constr_mat[i][match_constr_idx[lam][j]], &tmp_pf[i], roots);
    //       }
    //     }
    //   }
    // }

    // for (int lam, n=0; n<num_cum_eig; n++) {
    //   lam = cum_eig[n];
    //   for (int j=0; j<match_constr_dim[lam]; j++) {
    //     if (pow_behav[lam][j] >= 0 || sol_log_len[lam0] > 0) {
    //       // set column to identity
    //       for (int i=0; i<b_len; i++) {
    //         if (i == match_constr_idx[lam][j]) {
    //           mpc_set_ui(match_constr_mat[i][match_constr_idx[lam][j]], 1, MPFR_RNDN);
    //         } else {
    //           mpc_set_ui(match_constr_mat[i][match_constr_idx[lam][j]], 0, MPFR_RNDN);
    //         }
    //       }
    //     }
    //   }
    // }

    // evaluate limit of particular solution
    // if (sb_len>0 && 0) {
    // // if (sb_len>0 && sol_log_len[lam0] > 0) {
    //   // cout << "EVALUATE PARTICULAR SOLUTION IN ZERO" << endl;
    //   for (int lamp, np=0; np<num_cum_eig; np++) {
    //     lamp = cum_eig[np];
    //     if (sol_log_len[lamp] == 0) {
    //       continue;
    //     }
    //     // if (lamp != lam0) {
    //     //   continue;
    //     // }
    //     // cout << "lamp = " << lamp << endl;
    //     for (int i=0; i<offset; i++) {
    //       // cout << "i = " << i << endl;
    //       poly_frac_build(&pf_sol[offset+i]);
    //       poly_frac_set_coeffs(&pf_sol[offset+i], psol[np][0][i][0], tp_rank+1);
    //       // poly_frac_print(&pf_sol[offset+i]);
    //       pf_sol[offset+i].den_deg = 0;
    //       pf_sol[offset+i].nroots = nroots;
    //       pf_sol[offset+i].mults = new int[nroots];
    //       for (int k=0; k<nroots; k++) {
    //         pf_sol[offset+i].mults[k] = 0;
    //       }
    //     }
    //     poly_frac_rk2_mul_pf_rk1(tmp_pf, tmat, pf_sol, b_len, b_len, roots);
    //     // for (int i=0; i<b_len; i++) {
    //     //   cout << "i = " << i << endl;
    //     //   poly_frac_print(&tmp_pf[i]);
    //     // //   // poly_frac_print(&tmat_times_sol[i]);
    //     // //   poly_frac_add_pf(
    //     // //     &tmp_pf[i],
    //     // //     &tmp_pf[i],
    //     // //     &tmat_times_sol[i],
    //     // //     roots, NULL
    //     // //   );
    //     // //   // poly_frac_print(&tmp_pf[i]);
    //     // }
    //   }
    // }

    // removed on 2024.06.03
    // // CONSTANT TERM FOR THE FREE CONSTRAINTS
    // cout << "build constant term for the free constraints" << endl;
    // count = 0;
    // for (int i=0; i<b_len; i++) {
    //   for (int lam, n=0; n<num_cum_eig; n++) {
    //     lam = cum_eig[n];
    //     if (block_log_len[lam] == 0 || is_in_bound[lam] == 0) {
    //       continue;
    //     }
    //     for (int p=pow_behav_min[lam][i]; p<bound_behav[offset+i][lam]; p++) {
    //       // cout << "lam, p, count = " << lam << ", " << p << ", " << count << endl;
    //       // poly_frac_print(&pf_const_term[lam][i]);
    //       poly_frac_mul_sym_pow(&pf_const_term[lam][i], &pf_const_term[lam][i], -p);
    //       // poly_frac_print(&tmat_times_psol[i]);
    //       poly_frac_eval_sym_zero(&tmpc, &pf_const_term[lam][i], roots);
    //       mpc_neg(constr_const_term[count], tmpc, MPFR_RNDN);
    //       count++;
    //     }
    //     if (count == free_constr) {
    //       break;
    //     }
    //   }
    //   if (count == free_constr) {
    //     break;
    //   }
    // }

    // CONSTANT TERM
    // initialize to external boundary
    for (int i=0; i<b_len; i++) {
      mpc_set(ext_constr_const_term[i], (*solutions)[offset+i][0], MPFR_RNDN);
    }

    // CONSTANT TERM FOR EXTERNAL CONSTRAINTS
    if (print) cout << "build constant term for the external constraints" << endl;
    // if (sb_len>0 && sol_log_len[lam0] > 0) {
    if ((sb_len > 0 || offset > 0) && num_constr > free_constr) {
      for (int lamp, np=0; np<num_cum_eig; np++) {
        lamp = cum_eig[np];
        if (sol_log_len[lamp] == 0 || is_in_bound[lamp] == 0) {
          continue;
        }
        for (int i=0; i<b_len; i++) { // sum only on selected external constraints
          if (print) cout << "i = " << i << endl;
          // #occhio
          if (const_term_pow_behav[lamp][i] == bound_behav[offset+i][lamp]) {
            poly_frac_mul_sym_pow(&tmp_pf[i], &pf_const_term[lamp][i], -const_term_pow_behav[lamp][i]);
            poly_frac_eval_sym_zero(&tmpc, &tmp_pf[i], roots);
            mpc_sub(
              ext_constr_const_term[i], ext_constr_const_term[i], tmpc,
              MPFR_RNDN
            );
          } else if (const_term_pow_behav[lamp][i] < bound_behav[offset+i][lamp]) {
            poly_frac_extract_oder(&tmpc, &pf_const_term[lamp][i], roots, nroots, -bound_behav[offset+i][lamp]);
            mpc_sub(
              ext_constr_const_term[i], ext_constr_const_term[i], tmpc,
              MPFR_RNDN
            );
            // perror("WARNING: while matching in zero: constant term (prev. sol. and part. sol.) \
            //   have power behaviour lower than external boundary. \
            //   External constraint cannot be satisfied, unless this is a free constraint."
            // );
          }
        }
      }
    }
    // if (offset == 7) {exit(0);}

    // shift back
    for (int i=0; i<b_len; i++) {
      tmat[i] -= offset;
    }
    tmat -= offset;

    //////
    // SOLVE MATCHING CONSTRAINTS
    //////
    goto_solve:

    // ELIMINATE COLUMNS ASSOCIATED WITH NULL COEFFS
    if (print) cout << "ELIMINATE COLUMNS ASSOCIATED WITH NULL COEFFS" << endl;
    int count_col_idx;
    mpc_t **constr_final;
    malloc_rk2_tens(constr_final, b_len, b_len);
    init_rk2_mpc(constr_final, b_len, b_len);
    for (int i=0; i<num_constr; i++) {
      count_col_idx = 0;
      for (int lam, n=0; n<num_cum_eig; n++) {
        lam = cum_eig[n];
        if (block_log_len[lam] == 0) {
          continue;
        }
        for (int jc, j=0; j<match_constr_dim[lam]; j++) {
          jc = match_constr_idx[lam][j];
          if (coeff_is_null[jc] == 1) {continue;}
          if (print) cout << "lam, i, j, count_col_idx = " << lam << ", " << i << ", " << j << ", " << count_col_idx << endl;
          mpc_set(constr_final[i][count_col_idx++], constr[i][jc], MPFR_RNDN);
        }
      }
    }

    count_col_idx = 0;
    for (int jc=0; jc<b_len; jc++) {
      if (coeff_is_null[jc] == 1) {continue;}
      col_idx[jc] = count_col_idx;
      count_col_idx++;
    }
    
    if (print) {
    cout << "final free constr mat:" << endl;
    print_rk2_mpc(constr_final, num_constr, num_constr);
    }

    if (print) cout << "linear system matrix:" << endl;
    if (print) print_rk2_mpc(constr_final, num_constr, num_constr);
    malloc_rk2_tens(constr_inv, num_constr, num_constr);
    init_rk2_mpc(constr_inv, num_constr, num_constr);
    mp_inverse(constr_inv, constr_final, num_constr);
    if (print) cout << "linear system inverse matrix:" << endl;
    if (print) print_rk2_mpc(constr_inv, num_constr, num_constr);
    mpc_t *final_const_term;
    // if (free_constr < num_constr) {
    // if (0 < free_constr + num_null_coeff && free_constr + num_null_coeff < b_len) {
    if (0 < free_constr+num_null_coeff && free_constr+num_null_coeff < b_len) {
      if (print) cout << "selecting both free and external constraints for constant term" << endl;
      if (print) cout << "select indices for the constant term" << endl;
      final_const_term = new mpc_t[b_len];
      init_rk1_mpc(final_const_term, b_len);

      for (int k=0; k<free_constr; k++) {
        mpc_set(final_const_term[k], free_constr_const_term[k], MPFR_RNDN);
      }

      // select MI indices for the constant term
      for (int k=0; k<num_constr-free_constr; k++) {
        if (print) {
        cout << "idx = " << free_constr+k << endl;
        cout << "mi_constr_idx = " << mi_constr_idx[k] << endl;
        }
        mpc_set(final_const_term[free_constr+k], ext_constr_const_term[mi_constr_idx[k]], MPFR_RNDN);
      }
    // } else {
    //   final_const_term = ext_constr_const_term;
    } else if (free_constr == num_constr) {
      if (print) cout << "selecting only free constraints for constant term" << endl;
      final_const_term = free_constr_const_term;
    } else if (free_constr == 0) {
      if (print) cout << "selecting only external constraints for constant term" << endl;
      final_const_term = ext_constr_const_term;
    }
    if (print) {
    cout << "final constant term:" << endl;
    for (int i=0; i<num_constr; i++) {
      print_mpc(&final_const_term[i]); cout << endl;
    }
    }

    // multiply the inverse times the constant term of the system
    // cout << "multiply the inverse times the constant term of the system" << endl;
    for (int i=0; i<num_constr; i++) {
      mpc_set_d(match_constr_coeffs[i], 0, MPFR_RNDN);
      // for (int lam, n=0; n<num_cum_eig; n++) {
      //   lam = cum_eig[n];
      //   if (block_log_len[lam] == 0 || is_in_bound[lam] == 0) {
      //     continue;
      //   }
      //   for (int j=0; j<match_constr_dim[lam]; j++) {
      for (int j=0; j<num_constr; j++) {
        mpc_fma(
          match_constr_coeffs[i],
          constr_inv[i][j],
          final_const_term[j],
          match_constr_coeffs[i],
          MPFR_RNDN
        );
        // }
      }
    }

    if (print) {
    cout << endl << "linear system solution:" << endl;
    for (int i=0; i<num_constr; i++) {
      print_mpc(&match_constr_coeffs[i]); cout << endl;
    }
    }

    goto_hardcoded_2L_box:
    goto_lin_comb:
    // build linear combination of the homogeneous
    if (print) cout << endl << "build linear combination of the homogeneous" << endl;
    for (int lam, n=0; n<num_cum_eig; n++) {
      lam = cum_eig[n];
      if (block_log_len[lam] == 0 || is_in_bound[lam] == 0) {
        continue;
      }
      if (print) cout << "adding lam = " << lam << endl;
      for (int jc, j=0; j<match_constr_dim[lam]; j++) {
        jc = match_constr_idx[lam][j];
        if (coeff_is_null[jc] == 1) {continue;}
        for (int l=0; l<sol_log_len[lam]-1; l++) {
        if (print) cout << "jc, l = " << jc << ", " << l << endl;
        if (print) cout << "col_idx = " << col_idx[jc] << endl;
          for (int k=0; k<=eta_ord; k++) {
            for (int i=0; i<b_len; i++) {
              mpc_fma(
                sol[lam][l][offset+i][k],
                match_constr_coeffs[col_idx[jc]],
                // match_constr_coeffs[jc],
                hsol[n][l][i][jc][k],
                sol[lam][l][offset+i][k],
                MPFR_RNDN
              );
            }
          }
        }
      }
      delete[] constr_idx[lam];
    }
    // del_rk2_tens(constr, num_constr);
    // del_rk2_tens(constr_inv, num_constr);

    // FREE
    poly_frac_rk1_free(pf_sol, offset+b_len);
    delete[] pf_sol;
    poly_frac_rk1_free(tmp_pf, b_len);
    delete[] tmp_pf;
    for (int lam, n=0; n<num_cum_eig; n++) {
      lam = cum_eig[n];
      if (block_log_len[lam] == 0) {
        continue;
      }
      poly_frac_rk2_free(tmat_times_hsol[lam], b_len, match_constr_dim[lam]);
      del_rk2_tens(tmat_times_hsol[lam], b_len);
    }
    delete[] tmat_times_hsol;
    for (int lam, n=0; n<num_cum_eig; n++) {
      lam = cum_eig[n];
      if (block_log_len[lam] == 0) {
        continue;
      }
      poly_frac_rk2_free(mat_times_tmat_times_hsol[lam], b_len, match_constr_dim[lam]);
      del_rk2_tens(mat_times_tmat_times_hsol[lam], b_len);
    }
    delete[] mat_times_tmat_times_hsol;
    for (int lam, n=0; n<num_cum_eig; n++) {
      lam = cum_eig[n];
      for (int i=0; i<b_len; i++) {
        poly_frac_free(&pf_const_term[lam][i]);
      }
      delete[] pf_const_term[lam];
    }
    delete[] pf_const_term;
    for (int lam, n=0; n<num_cum_eig; n++) {
      lam = cum_eig[n];
      for (int i=0; i<b_len; i++) {
        poly_frac_free(&tmat_times_psol[lam][i]);
      }
      delete[] tmat_times_psol[lam];
    }
    delete[] tmat_times_psol;
    for (int lam, n=0; n<num_cum_eig; n++) {
      lam = cum_eig[n];      
      for (int i=0; i<b_len; i++) {
        poly_frac_free(&mat_times_tmat_times_psol[lam][i]);
      }
      delete[] mat_times_tmat_times_psol[lam];
    }
    delete[] mat_times_tmat_times_psol;
    for (int lam, n=0; n<num_cum_eig; n++) {
      lam = cum_eig[n];
      for (int i=0; i<b_len; i++) {
        poly_frac_free(&mat_times_prevsol[lam][i]);
      }
      delete[] mat_times_prevsol[lam];
    }
    delete[] mat_times_prevsol;
    for (int lam, n=0; n<num_cum_eig; n++) {
      lam = cum_eig[n];
      for (int i=0; i<b_len; i++) {
        poly_frac_free(&DE_const_term[lam][i]);
      }
      delete[] DE_const_term[lam];
    }
    delete[] DE_const_term;

    poly_frac_free(&pf_tmp);

  }

  // solve the system
  if (match_in_zero == 0) {
    // cout << "wp_bin orig = " << wp_bin << endl;
    int wp_bin_red = wp_bin * 0.5;  // #hard-coded
    if (wp_bin_red < wp2/4) {wp_bin_red = wp_bin;}  // #hard-coded
    // cout << "wp_bin red = " << wp_bin_red << endl;
    //// invert the matrix
    if (print) {
    cout << endl; cout << "linear system matrix:" << endl;
    print_rk2_mpc(match_constr_mat, b_len, b_len);
    }
    mp_inverse(match_constr_inv_mat, match_constr_mat, b_len);
    if (print) {
    cout << endl; cout << "linear system inverse matrix:" << endl;
    print_rk2_mpc(match_constr_inv_mat, b_len, b_len);
    cout << endl; cout << "constant term:" << endl;
    for (int i=0; i<b_len; i++) {
      print_mpc(&(*solutions)[offset+i][0]); cout << endl;
    }
    }

    //// multiply the inverse times the constant term of the system
    // cout << "multiply the inverse times the constant term of the system" << endl;
    mpc_t tmpc; mpc_init3(tmpc, wp2, wp2);
    for (int i=0; i<b_len; i++) {
      mpc_set_d(match_constr_coeffs[i], 0, MPFR_RNDN);
      int j;
      for (j=0; j<b_len; j++) {
        // mpc_fma(match_constr_coeffs[i],
        //   match_constr_inv_mat[i][j],
        //   (*solutions)[offset+j][0],
        //   match_constr_coeffs[i],
        //   MPFR_RNDN);

        exp_rel_err_mpc_mul(tmpc, match_constr_inv_mat[i][j], (*solutions)[offset+j][0], MPFR_RNDN, wp_bin_red);
        exp_rel_err_mpc_add(match_constr_coeffs[i], match_constr_coeffs[i], tmpc, MPFR_RNDN, wp_bin_red);
      }

      // // rel err prune on last contribution
      // mpc_mul(tmpc, match_constr_inv_mat[i][j], (*solutions)[offset+j][0], MPFR_RNDN);
      // rel_err_mpc_add(match_constr_coeffs[i], match_constr_coeffs[i], tmpc, MPFR_RNDN);
    }
    mpc_clear(tmpc);
    if (print) {
      cout << "solution coefficients:" << endl;
      print_rk1_mpc(match_constr_coeffs, b_len);
    }
  }

  // cout << endl; cout << "linear system solution:" << endl;
  // for (int i=0; i<b_len; i++) {
  //   print_mpc(&match_constr_coeffs[i]); cout << endl;
  // }


  if (match_in_zero == 0) {
    // build linear combination of the homogeneous
    // cout << endl; cout << "build linear combination of the homogeneous" << endl;
    for (int lam, n=0; n<num_cum_eig; n++) {
      lam = cum_eig[n];
      if (block_log_len[lam] == 0) {
        continue;
      }
      for (int j=0; j<match_constr_dim[lam]; j++) {
        // if (match_in_zero) {
        //   cout << "pow behav = " << pow_behav[lam][j] << endl;
        //   if (pow_behav[lam][j] >= 0 || sol_log_len[lam0] > 0) {
        //     continue;
        //   }
        // }
        // cout << "adding lam = " << lam << endl;
        for (int l=0; l<sol_log_len[lam]-1; l++) {
          for (int k=0; k<=eta_ord; k++) {
            for (int i=0; i<b_len; i++) {
              // print_mpc(&sol[lam][l][offset+i][k]); cout << endl;
              // print_mpc(&match_constr_coeffs[match_constr_idx[lam][j]]); cout << endl;
              // print_mpc(&hsol[n][l][i][match_constr_idx[lam][j]][k]); cout << endl;
              mpc_fma(
                sol[lam][l][offset+i][k],
                match_constr_coeffs[match_constr_idx[lam][j]],
                hsol[n][l][i][match_constr_idx[lam][j]][k],
                sol[lam][l][offset+i][k],
                MPFR_RNDN
              );
            }
          }
        }
      }
    }
  }

  // add contribution from particular solution
  if (sb_len > 0) {
    // cout << endl; cout << "add contribution from particular solution" << endl;
    for (int lam, n=0; n<num_cum_eig; n++) {
      lam = cum_eig[n];
      if (sol_log_len[lam] == 0) {
        continue;
      }
      // cout << "lam = " << lam << endl;
      int higher_log = 1;
      for (int l=sol_log_len[lam]-2; l>=0; l--) {
        int k_start = 0;
        // ZERO ORDER
        if (higher_log && l>=rhs_log_len[lam]) {
          k_start = 1;
          int col_is_zero = 1;
          for (int i=0; i<b_len; i++) {
            // mpc_t tmpc; mpc_init3(tmpc, wp2, wp2);
            // mpc_add(
            //   tmpc,
            //   sol[lam][l][offset+i][0],
            //   psol[n][l][i][0][0],
            //   MPFR_RNDN
            // );

            exp_rel_err_mpc_add(
              sol[lam][l][offset+i][0],
              sol[lam][l][offset+i][0],
              psol[n][l][i][0][0],
              MPFR_RNDN,
              wp_bin
            );

            // cout << "rel err sum =  "; print_mpc(&sol[lam][l][offset+i][0]); cout << endl;
            // cout << "standard sum = "; print_mpc(&tmpc); cout << endl;

            // if (!mpc_equal_within_tol(tmpc, sol[lam][l][offset+i][0])) {
            //   cout << "rel err sum not corresponding with standard sum" << endl;
            //   perror("rel err sum not corresponding with standard sum");
            // }

            if (!mpc_zero_p(sol[lam][l][offset+i][0])) {
              col_is_zero = 0;
            } else {
              // cout << "lam, l, i = " << lam << ", " << l << ", " << i << endl;
              // cout << "rel err sum =  "; print_mpc(&sol[lam][l][offset+i][0]); cout << endl;
              // cout << "standard sum = "; print_mpc(&tmpc); cout << endl;

              // mpc_add(
              //   tmpc,
              //   sol[lam][l][offset+i][1],
              //   psol[n][l][i][0][1],
              //   MPFR_RNDN
              // );

              // cout << "standard order one = "; print_mpc(&tmpc); cout << endl;
            }
            // mpc_clear(tmpc);
          }
          
          if (col_is_zero) {
            // cout << "col is zero" << endl;
            for (int i=0; i<b_len; i++) {
              for (int k=1; k<=eta_ord; k++) {
                mpc_set_ui(sol[lam][l][offset+i][k], 0, MPFR_RNDN);
              }
            }
            continue;
          } else {
            higher_log = 0;
          }
        }



          // if (l >= rhs_log_len[lam]) {
          //   if (!mpfr_zero_p(mpc_realref(psol[n][l][i][0][0])) && !mpfr_zero_p(mpc_imagref(psol[n][l][i][0][0]))) {
          //     cout << "sol = "; print_mpc(&sol[lam][l][offset+i][0]); cout << endl;
          //     cout << "psol = "; print_mpc(&psol[n][l][i][0][0]); cout << endl;
          //     mpc_div(sol[lam][l][offset+i][0], sol[lam][l][offset+i][0], psol[n][l][i][0][0], MPFR_RNDN);
          //     cout << "ratio = "; print_mpc(&sol[lam][l][offset+i][0]); cout << endl;
          //     cout << "tol = "; mpfr_out_str(stdout, 10, 0, mpfr_tol, MPFR_RNDN); cout << endl;
          //     // mpfr_t tmpfr; mpfr_init2(tmpfr, wp2);
          //     // mpc_abs(tmpfr, sol[lam][l][offset+i][0], MPFR_RNDN);
          //     // cout << "abs = "; mpfr_out_str(stdout, 10, 0, tmpfr, MPFR_RNDN); cout << endl;
          //     // mpfr_clear(tmpfr);
          //   }
          //   if (mpc_lessthan_tol(sol[lam][l][offset+i][0])) {
          //     // if first eta-order is zero and there is no rhs,
          //     // all the eta-orders are zero
          //     for (int k=0; k<=eta_ord; k++) {
          //       mpc_set_ui(sol[lam][l][offset+i][k], 0, MPFR_RNDN);
          //     }
          //     continue;
          //   }
          // }

        for (int k=k_start; k<=eta_ord; k++) {
          for (int i=0; i<b_len; i++) {
          // if (k == 0) {
          //   cout << "l, i, k = " << l << ", " << i << ", " << k << endl;
          //   print_mpc(&sol[lam][l][offset+i][k]); cout << endl;
          //   print_mpc(&psol[n][l][i][0][k]); cout << endl;
          // }
            mpc_add(sol[lam][l][offset+i][k],
              sol[lam][l][offset+i][k],
              psol[n][l][i][0][k],
              MPFR_RNDN
            );
          }
        }
      }
    }
  }

  // cout << "delete" << endl;
  for (int lam, n=0; n<num_cum_eig; n++) {
    lam = cum_eig[n];
    if (block_log_len[lam] == 0) {
      continue;
    }
    delete[] match_constr_idx[lam];
  }
  delete[] match_constr_idx;
  delete[] match_constr_coeffs;
  delete[] match_constr_dim;
  mpc_rk2_clear(match_constr_mat, b_len, b_len);
  del_rk2_tens(match_constr_mat, b_len);
  mpc_rk2_clear(match_constr_inv_mat, b_len, b_len);
  del_rk2_tens(match_constr_inv_mat, b_len);
}


// #pair match_sol_around_zero
void pf_limit_in_zero(
  // OUTPUT
  mpc_t *final_sol,
  // INTPUT
  mpc_t ****sol, struct poly_frac **tmat,
  int dim, int lam0, int tp_rank,
  mpc_t *roots, int nroots,
  int rad_exp
) {
  // convert solution to poly_frac up to Poinc. rank of transformation
  // and keeping only the log power corresponing to the null eigenvalue,
  // since this is the only one that contributes to the limit
  // cout << "convert solution to poly_frac up to Poinc. rank of transformation" << endl;
  struct poly_frac *pf_sol;
  pf_sol = new struct poly_frac[dim];
  // cout << "solution:" << endl;
  for (int i=0; i<dim; i++) {
    // cout << "i = " << i << endl;
    // cout << "input poly:" << endl;
    // print_poly(sol[lam0][0][i], tp_rank+1);
    poly_frac_build(&pf_sol[i]);
    poly_frac_set_coeffs(&pf_sol[i], sol[lam0][0][i], tp_rank+1);
    pf_sol[i].den_deg = 0;
    pf_sol[i].nroots = nroots;
    pf_sol[i].mults = new int[nroots];
    for (int k=0; k<nroots; k++) {
      pf_sol[i].mults[k] = 0;
    }
    // poly_frac_print(&pf_sol[i]);
  }

  int wp_bin = - mpfr_log2_int(mpfr_tol);
  wp_bin *= 0.5;  // #hard-coded

  // apply matrix
  // cout << "apply matrix" << endl;
  // poly_frac_rk2_mul_pf_rk1(pf_sol, tmat, pf_sol, dim, dim, roots);
  rel_err_poly_frac_rk2_mul_pf_rk1(pf_sol, tmat, pf_sol, dim, dim, roots, wp_bin);

  // take limit in zero
  // cout << "take limit in zero" << endl;
  // poly_frac_rk1_prune_rel_tol(pf_sol, wp_bin, dim);
  poly_frac_rk1_prune_radius(pf_sol, wp_bin, rad_exp, dim);
  poly_frac_rk1_eval_sym_zero(final_sol, pf_sol, roots, dim);

  // FREE
  poly_frac_rk1_free(pf_sol, dim);
  delete[] pf_sol;
}


void pf_limit_in_zero_block(
  // OUTPUT
  mpc_t *final_sol,
  // IN-OUT
  struct poly_frac *pf_sol,
  // INTPUT
  mpc_t ****sol, struct poly_frac **tmat,
  int offset, int b_len, int lam0, int tp_rank,
  mpc_t *roots, int nroots,
  int rad_exp
) {
  struct poly_frac *pf_tmat_sol;
  pf_tmat_sol = new struct poly_frac[b_len];
  poly_frac_rk1_build(pf_tmat_sol, b_len);

  // fill poly_frac with the lam0 contribution to the block sol
  for (int i=offset; i<offset+b_len; i++) {    
    poly_frac_build(&pf_sol[i]);
    poly_frac_set_coeffs(&pf_sol[i], sol[lam0][0][i], tp_rank+1);
    pf_sol[i].den_deg = 0;
    pf_sol[i].nroots = nroots;
    pf_sol[i].mults = new int[nroots];
    for (int k=0; k<nroots; k++) {
      pf_sol[i].mults[k] = 0;
    }
  }

  int wp_bin = - mpfr_log2_int(mpfr_tol);
  wp_bin *= 0.5;  // #hard-coded
  
  rel_err_poly_frac_rk2_mul_pf_rk1(
    pf_tmat_sol,
    tmat+offset, pf_sol,
    b_len, offset+b_len,
    roots, wp_bin
  );

  poly_frac_rk1_prune_radius(pf_tmat_sol, wp_bin, rad_exp, b_len);
  poly_frac_rk1_eval_sym_zero(final_sol, pf_tmat_sol, roots, b_len);

  // FREE
  poly_frac_rk1_free(pf_tmat_sol, b_len);
  delete[] pf_tmat_sol;

}


void mpc_Log(
  // OUTPUT
  mpc_t *out,
  // INPUT
  mpc_t *in
) {
  mpc_t *in_copy;
  mpc_t in_copy_tmp;
  if (out == in) {
    mpc_init3(in_copy_tmp, wp2, wp2);
    mpc_set(in_copy_tmp, *in, MPFR_RNDN);
    in_copy = &in_copy_tmp;
  } else {
    in_copy = in;
  }
  cout << "in = "; mpc_out_str(stdout, 10, 0, *in_copy, MPFR_RNDN); cout << endl;

  mpfr_t tmpfr;
  mpfr_init2(tmpfr, wp2);
  mpfr_t mpfr_pi, mpfr_2pi;
  mpfr_init2(mpfr_pi, wp2);
  mpfr_init2(mpfr_2pi, wp2);
  mpfr_const_pi(mpfr_pi, MPFR_RNDN);
  mpfr_mul_ui(mpfr_2pi, mpfr_pi, 2, MPFR_RNDN);

  mpfr_t abs;
  mpfr_init2(abs, wp2);

  // real part
  mpc_abs(abs, *in_copy, MPFR_RNDN);
  // mpfr_out_str(stdout, 10, 0, abs, MPFR_RNDN); cout << endl;
  mpfr_log(mpc_realref(*out), abs, MPFR_RNDN);
  // print_mpc(out); cout << endl;

  // imaginary part
  mpfr_t *arg;
  arg = &mpc_imagref(*out);
  mpc_arg(*arg, *in_copy, MPFR_RNDN);
  cout << "arg = "; mpfr_out_str(stdout, 10, 0, *arg, MPFR_RNDN); cout << endl;
  // #TBD: update tol propagating error through the Arg function
  if (mpfr_sign_within_tol(*arg) < 0) {
    mpfr_add(tmpfr, *arg, mpfr_pi, MPFR_RNDN);
    if (!mpfr_lessthan_tol(tmpfr)) {
      mpfr_add(*arg, *arg, mpfr_2pi, MPFR_RNDN);
    }
  }
}


// #TBD: compute gamma_eps globally and use it to get gamma[eps-1]
void get_tadpole(
  // OUTPUT
  mpc_t *out,
  // INPUT
  mpc_t *mpc_eta, mpc_t *PS_ini, mpc_t *PS_fin, mpc_t *eps,
  struct poly_frac *tmat,
  mpc_t *eig_list, int lam, int analytic_continuation
) {
  // COMPUTE TADPOLE
  mpc_t omeps;
  mpc_init3(omeps, wp2, wp2);
  mpc_set_ui(omeps, 1, MPFR_RNDN);
  mpc_sub(omeps, omeps, *eps, MPFR_RNDN);
  cout << "one minus epsilon = "; print_mpc(&omeps); cout << endl;
  
  mpc_t mass;
  mpc_init3(mass, wp2, wp2);
  mpc_sub(mass, *PS_fin, *PS_ini, MPFR_RNDN);
  mpc_fma(mass, *mpc_eta, mass, *PS_ini, MPFR_RNDN);
  // mpc_sqr(mass, mass, MPFR_RNDN);
  
  // mpc_mul(mass, mass, *mpc_eta, MPFR_RNDN);
  cout << "mass = "; print_mpc(&mass); cout << endl;

  mpc_t pow;
  mpc_init3(pow, wp2, wp2);
  mpc_mul_ui(pow, omeps, 2, MPFR_RNDN);
  mpc_pow(*out, mass, pow, MPFR_RNDN);

  // mpfr_t gamma_eps;
  // mpfr_init2(gamma_eps, wp2);
  // mpfr_gamma(gamma_eps, mpc_realref(*eps), MPFR_RNDN);
  // mpc_mul_fr(*out, *out, gamma_eps, MPFR_RNDN);

  mpfr_t gamma;
  mpfr_init2(gamma, wp2);
  mpc_neg(omeps, omeps, MPFR_RNDN);
  mpfr_gamma(gamma, mpc_realref(omeps), MPFR_RNDN);
  cout << "gamma = "; mpfr_out_str(stdout, 10, 0, gamma, MPFR_RNDN); cout << endl;
  mpc_mul_fr(*out, *out, gamma, MPFR_RNDN);
  mpc_neg(*out, *out, MPFR_RNDN);
  cout << "original tadpole:"; cout << endl; print_mpc(out); cout << endl;

  // APPLY TRANSFORMATION MATRIX
  cout << "tmat:" << endl; poly_frac_print(tmat);
  struct poly_frac pf;
  poly_frac_build(&pf);
  poly_frac_mul_mpc(&pf, tmat, out);
  mpc_t *roots;
  roots = new mpc_t[1];
  mpc_init3(roots[0], wp2, wp2);
  mpc_set_ui(roots[0], 0, MPFR_RNDN);
  poly_frac_eval_value(out, &pf, roots, mpc_eta);
  cout << "transformed tadpole;"; cout << endl; print_mpc(out); cout << endl;

  // ANALYTIC CONTINUATION
  if (mpfr_sign_within_tol(mpc_realref(*mpc_eta)) > 0 ||\
    mpfr_sign_within_tol(mpc_imagref(*mpc_eta)) > 0) {
    analytic_continuation = 0;
  }
  mpc_t log_eta;
  mpc_init3(log_eta, wp2, wp2);
  mpc_t tmpc;
  mpc_init3(tmpc, wp2, wp2);
  if (analytic_continuation) {
    // cout << "ANALYTIC CONTINUATION: eta = ";
    // print_mpc(mpc_eta);
    // cout << endl;
    mpc_neg(*mpc_eta, *mpc_eta, MPFR_RNDN);
    mpc_log(log_eta, *mpc_eta, MPFR_RNDN);
    mpc_set_ui(tmpc, 0, MPFR_RNDN);
    mpfr_const_pi(mpc_imagref(tmpc), MPFR_RNDN);
    if (analytic_continuation == -1) {
      mpc_neg(tmpc, tmpc, MPFR_RNDN);
    }
    mpc_add(log_eta, log_eta, tmpc, MPFR_RNDN);
    // mpc_sub(log_eta, log_eta, tmpc, MPFR_RNDN);
    mpc_neg(*mpc_eta, *mpc_eta, MPFR_RNDN);
    // cout << "log = ";
    // print_mpc(&log_eta);
    // cout << endl;
  } else {
    mpc_log(log_eta, *mpc_eta, MPFR_RNDN);
  }
  mpc_t eta_lam;
  mpc_init3(eta_lam, wp2, wp2);
  if (analytic_continuation) {
    // analytic continuation
    //// eta^lam = exp{lam*log(eta)}
    // cout << "ANALYTIC CONTINUATION: eta = ";
    // print_mpc(mpc_eta);
    // cout << endl;
    mpc_sub(eta_lam, pow, eig_list[lam], MPFR_RNDN);
    mpc_mul(eta_lam, eta_lam, log_eta, MPFR_RNDN);
    mpc_exp(eta_lam, eta_lam, MPFR_RNDN);
  } else {
    mpc_pow(eta_lam, *mpc_eta, eig_list[lam], MPFR_RNDN);
  }

  // MULTIPLY TIMES ETA^(POW-LAM)
  mpc_mul(*out, *out, eta_lam, MPFR_RNDN);
}


void get_massless_bubble(
  // OUTPUT
  mpc_t *out,
  // INPUT
  mpc_t *mpc_eta, mpc_t *PS_ini, mpc_t *PS_fin, mpc_t *eps,
  struct poly_frac *tmat,
  mpc_t *eig_list, int lam, int analytic_continuation
) {
  // COMPUTE MASSLESS BUBBLE
  mpc_t omeps;
  mpc_init3(omeps, wp2, wp2);
  mpc_set_ui(omeps, 1, MPFR_RNDN);
  mpc_sub(omeps, omeps, *eps, MPFR_RNDN);
  cout << "one minus epsilon = "; print_mpc(&omeps); cout << endl;
  
  mpc_t p2;
  mpc_init3(p2, wp2, wp2);
  mpc_sub(p2, *PS_fin, *PS_ini, MPFR_RNDN);
  // mpc_mul(p2, p2, *mpc_eta, MPFR_RNDN);
  mpc_fma(p2, *mpc_eta, p2, *PS_ini, MPFR_RNDN);
  cout << "p2 = "; print_mpc(&p2); cout << endl;

  mpc_t pow;
  mpc_init3(pow, wp2, wp2);
  mpc_neg(pow, *eps, MPFR_RNDN);
  mpc_neg(*out, p2, MPFR_RNDN);
  mpc_pow_si(*out, *out, -1, MPFR_RNDN);
  // mpc_pow(*out, *out, pow, MPFR_RNDN);
  mpc_Log(out, out);
  cout << "log = "; print_mpc(out); cout << endl;
  mpc_mul(*out, *out, *eps, MPFR_RNDN);
  mpc_exp(*out, *out, MPFR_RNDN);

  mpfr_t gamma_eps;
  mpfr_init2(gamma_eps, wp2);
  mpfr_gamma(gamma_eps, mpc_realref(*eps), MPFR_RNDN);
  mpc_mul_fr(*out, *out, gamma_eps, MPFR_RNDN);
  cout << "gamma_eps = "; mpfr_out_str(stdout, 10, 0, gamma_eps, MPFR_RNDN); cout << endl;

  mpfr_t gamma1sq, gamma2;
  mpfr_init2(gamma1sq, wp2);
  mpfr_init2(gamma2, wp2);
  mpfr_gamma(gamma1sq, mpc_realref(omeps), MPFR_RNDN);
  mpfr_sqr(gamma1sq, gamma1sq, MPFR_RNDN);
  mpc_mul_ui(omeps, omeps, 2, MPFR_RNDN);
  mpfr_gamma(gamma2, mpc_realref(omeps), MPFR_RNDN);

  mpc_mul_fr(*out, *out, gamma1sq, MPFR_RNDN);
  mpc_div_fr(*out, *out, gamma2, MPFR_RNDN);
  // mpc_mul(*out, *out, pow, MPFR_RNDN);
  cout << "original bubble:"; cout << endl; print_mpc(out); cout << endl;

  // APPLY TRANSFORMATION MATRIX
  cout << "tmat:" << endl; poly_frac_print(tmat);
  struct poly_frac pf;
  poly_frac_build(&pf);
  poly_frac_mul_mpc(&pf, tmat, out);
  mpc_t *roots;
  roots = new mpc_t[1];
  mpc_init3(roots[0], wp2, wp2);
  mpc_set_ui(roots[0], 0, MPFR_RNDN);
  poly_frac_eval_value(out, &pf, roots, mpc_eta);
  cout << "transformed bubble:"; cout << endl; print_mpc(out); cout << endl;

  // ANALYTIC CONTINUATION
  if (mpfr_sign_within_tol(mpc_realref(*mpc_eta)) > 0 ||\
    mpfr_sign_within_tol(mpc_imagref(*mpc_eta)) > 0) {
    analytic_continuation = 0;
  }
  cout << "bubble analytic continuation: " << analytic_continuation << endl;
  mpc_t log_eta;
  mpc_init3(log_eta, wp2, wp2);
  mpc_t tmpc;
  mpc_init3(tmpc, wp2, wp2);
  if (analytic_continuation) {
    // cout << "ANALYTIC CONTINUATION: eta = ";
    // print_mpc(mpc_eta);
    // cout << endl;
    mpc_neg(*mpc_eta, *mpc_eta, MPFR_RNDN);
    mpc_log(log_eta, *mpc_eta, MPFR_RNDN);
    mpc_set_ui(tmpc, 0, MPFR_RNDN);
    mpfr_const_pi(mpc_imagref(tmpc), MPFR_RNDN);
    if (analytic_continuation == -1) {
      mpc_neg(tmpc, tmpc, MPFR_RNDN);
    }
    mpc_add(log_eta, log_eta, tmpc, MPFR_RNDN);
    // mpc_sub(log_eta, log_eta, tmpc, MPFR_RNDN);
    mpc_neg(*mpc_eta, *mpc_eta, MPFR_RNDN);
    // cout << "log = ";
    // print_mpc(&log_eta);
    // cout << endl;
  } else {
    mpc_log(log_eta, *mpc_eta, MPFR_RNDN);
  }
  mpc_t eta_lam;
  mpc_init3(eta_lam, wp2, wp2);
  if (analytic_continuation) {
    // analytic continuation
    //// eta^lam = exp{lam*log(eta)}
    // cout << "ANALYTIC CONTINUATION: eta = ";
    // print_mpc(mpc_eta);
    // cout << endl;
    // mpc_sub(eta_lam, pow, eig_list[lam], MPFR_RNDN);
    // mpc_mul(eta_lam, eta_lam, log_eta, MPFR_RNDN);
    mpc_mul(eta_lam, eig_list[lam], log_eta, MPFR_RNDN);
    mpc_exp(eta_lam, eta_lam, MPFR_RNDN);
  } else {
    mpc_pow(eta_lam, *mpc_eta, eig_list[lam], MPFR_RNDN);
  }

  // // MULTIPLY TIMES ETA^(POW-LAM)
  // mpc_mul(*out, *out, eta_lam, MPFR_RNDN);
  // MULTIPLY TIMES ETA^(-LAM)
  mpc_div(*out, *out, eta_lam, MPFR_RNDN);
}


void solve_zero(
  // IN-OUTPUT
  mpc_t ***solutions,
  // INPUT
  int dim, struct poly_frac **mat_ep,
  struct poly_frac **tmat, struct poly_frac **inv_tmat,
  int npoles, mpc_t *poles,
  mpc_t *matching_point,
  int nblocks, int **prof, int **sb_grid, int eta_ord,
  int num_classes, mpc_t *eig_list, int *eq_class, int *eig_grid,
  mpc_t *roots, int nroots,
  int cross, mpc_t *eta_target, int analytic_cont,
  int *is_mass, int *skip_inv, int ninvs, mpc_t *PS_ini, mpc_t *PS_fin, char *eps_str,
  int try_analytic,
  int **bound_behav, int **mi_eig, int *mi_eig_num,
  FILE *terminal
) {
  /*
  INPUT:
    - cross: boolean that activate singularity crossing:
      TRUE: cross singularity in zero and evaluate solution at the point provided by eta_target
        instead of computing the limit in zero;
      FALSE: no crossing, just the limit of the solution in zero;

      BOUND TO:
      - TRUE -> eta_target: must be non-NULL;
      - FALSE -> eta_target: suggested to be NULL;
  */
  int print = 0;
  int debug = 0;
  int prt_ord = 10;

  // check whether matching point is zero
  int match_in_zero;
  if (mpfr_zero_p(mpc_realref(*matching_point)) && mpfr_zero_p(mpc_imagref(*matching_point))) {
    match_in_zero = 1;
  } else {
    match_in_zero = 0;
  }
  if (print) cout << "match in zero: " << match_in_zero << endl;

  // find label of null eigenvalue
  mpfr_t tmp;
  mpfr_init2(tmp, wp2);
  int lam0 = 0;
  for (int lam=0; lam<num_classes; lam++) {
    // check if eigenvalue is zero
    mpfr_abs(tmp, mpc_realref(eig_list[lam]), MPFR_RNDN);
    if (mpfr_less_p(tmp, mpfr_tol)) {
      mpfr_abs(tmp, mpc_imagref(eig_list[lam]), MPFR_RNDN);
      if (mpfr_less_p(tmp, mpfr_tol)) {
        lam0 = lam;
        break; 
      }
    }
  }
  if (print) cout << "lam0 = " << lam0 << endl;

  int tp_rank;
  tp_rank = poly_frac_Poinc_rank(tmat, dim, dim);

  mpc_t **tmat_eval;
  malloc_rk2_tens(tmat_eval, dim, dim);
  init_rk2_mpc(tmat_eval, dim, dim);
  mpc_t *tmp_vec = new mpc_t[dim];
  init_rk1_mpc(tmp_vec, dim);

  if (print) {
  cout << "UNTRANSFORMED BOUNDARY:" << endl;
  for (int i=0; i<dim; i++) {
    print_mpc(&(*solutions)[i][0]); cout << endl;
  }
  }
  if (match_in_zero == 0) {
    // transform boundary
    // apply_tmat(*solutions, inv_tmat, dim, etav);
    poly_frac_rk2_eval_value(
      tmat_eval, inv_tmat, roots, matching_point,
      dim, dim
    );
    // cout << "tmat eval:" << endl;
    // print_rk2_mpc(tmat_eval, dim, dim);
    // mpc_t **prova;
    // malloc_rk2_tens(prova, dim, dim);
    // init_rk2_mpc(prova, dim, dim);
    // mpc_rk2_from_file((char*)"eval_tmat.txt", prova, dim, dim);
    // mpc_rk2_compare(prova, inv_tmat_eval, dim, dim);
    mpc_rk2_mul_mpc_rk2_slice(
      tmp_vec, tmat_eval, *solutions,
      dim, dim, 0
    );
    if (print) cout << "TRANSFORMED BOUNDARY:" << endl;
    for (int i=0; i<dim; i++) {
      mpc_set((*solutions)[i][0], tmp_vec[i], MPFR_RNDN);
      if (print) {print_mpc(&(*solutions)[i][0]); cout << endl;}
    }
  }

  // compute max block len
  int b_len_max = 0, b_len_tmp;
  for (int b=0; b<nblocks; b++) {
    b_len_tmp = prof[b][1] - prof[b][0] + 1;
    if (b_len_tmp > b_len_max) {
      b_len_max = b_len_tmp;
    }
  }

  mpc_t *****hsol, *****psol, ****sol;
  malloc_rk4_tens(sol, num_classes, b_len_max+2, dim, eta_ord+1);
  init_rk4_mpc(sol, num_classes, b_len_max+2, dim, eta_ord+1);
  // initialize to zero
  for (int i=0; i<dim; i++) {
    for (int k=0; k<eta_ord+1; k++) {
      mpc_set_d(sol[lam0][0][i][k], 0, MPFR_RNDN);
    }
  }

  mpfr_t mpfr_mem; mpfr_init2(mpfr_mem, wp2);
  // cout << endl; cout << "memory for the solution: ";
  // cout << num_classes << " x " << b_len_max+2 << " x " << dim << " x " << eta_ord+1 << " x ";
  // cout << mpfr_get_memory_usage(mpfr_mem)*2 << " bytes = ";
  // cout << num_classes*(b_len_max+2)*dim*(eta_ord+1)*mpfr_get_memory_usage(mpfr_mem)*2 << " bytes" << endl;
  mpfr_clear(mpfr_mem);
  //////
  // CHECK INPUTS
  //////
  // cout << endl; cout << "EIGENVALUE ANALYSIS:" << endl;
  // cout << "num classes = " << num_classes << endl;
  // cout << "eigenvalues list:" << endl;
  // for (int i=0; i<num_classes; i++) {
  //   cout << i << ": "; print_mpc(&eig_list[i]); cout << endl;
  // }
  // cout << endl;
  // // check normalization
  // cout << "eigenvalues grid:" << endl;
  // cout << "line, grid, class label" << endl;
  // for (int i=0; i<dim; i++) {
  //   cout << i << "\t" << eig_grid[i] << "\t" << eq_class[i] << endl;
  // }
  // cout << endl;
  // getchar();

  //////
  // BLOCK LOOP
  //////
  mpc_t eta_shift;
  mpc_init3(eta_shift, wp2, wp2);
  mpc_set_ui(eta_shift, 0, MPFR_RNDN);
  int b_len, sb_len, *b_idx, *sb_idx;
  int *lcm_deg, *block_max_deg, *num_max_deg, *den_max_deg;
  lcm_deg = new int[nblocks];
  block_max_deg = new int[nblocks];
  num_max_deg = new int[nblocks];
  den_max_deg = new int[nblocks];
  bool flag_lcm = true;
  // ex *gnc_lcm;
  // gnc_lcm = new ex[nblocks];
  mpc_t **lcm, ***block, ****const_term;
  lcm = new mpc_t*[nblocks];
  mpc_t **lcm_orig;
  lcm_orig = new mpc_t*[nblocks];
  int **lcm_roots;
  lcm_roots = new int*[nblocks];
  mpc_t ****mpc_lcm_mat_ep;
  malloc_rk3_tens(mpc_lcm_mat_ep, dim, dim, 2);
  int offset = 0;
  int num_cum_eig = 0; //, num_block_eig = 0;
  int *cum_eig, *block_eig;
  cum_eig = new int[num_classes];
  int *block_log_len = new int[num_classes];
  int *sblock_log_len = new int[num_classes];
  int *rhs_log_len = new int[num_classes];
  int *sol_log_len = new int[num_classes];
  int **log_prof = new int*[num_classes];
  int **block_log_prof = new int*[num_classes];  // #deprecated: not used anymore
  for (int i=0; i<num_classes; i++) {
    block_log_len[i] = 0;
    sblock_log_len[i] = 0;
    sol_log_len[i] = 0;
    rhs_log_len[i] = 0;
    log_prof[i] = new int[dim];
    for (int j=0; j<dim; j++) {
      log_prof[i][j] = 0;
    }
  }
  int chain_len, curr_eig;

  // copy mat_ep
  struct poly_frac **lcm_mat_ep;
  malloc_rk2_tens(lcm_mat_ep, dim, dim);
  poly_frac_rk2_build(lcm_mat_ep, dim, dim);
  poly_frac_rk2_set_pf_rk2(lcm_mat_ep, mat_ep, dim, dim);
  poly_frac_rk2_mul_sym_pow(lcm_mat_ep, 1, dim, dim);

  mpc_t **sol_at_target;
  mpc_t ***cross_zero;
  if (cross) {
    malloc_rk2_tens(sol_at_target, dim, 1);
    init_rk2_mpc(sol_at_target, dim, 1);
  }

  // needed for limit in zero
  mpfr_t tmpfr; mpfr_init2(tmpfr, wp2);
  mpc_abs(tmpfr, *matching_point, MPFR_RNDN);
  int rad_exp = mpfr_get_exp(tmpfr);
  struct poly_frac *pf_sol = new struct poly_frac[dim];

  for (int b=0; b<nblocks; b++) {
    // for (int b=0; b<1; b++) {
    if (b > 0) {fprintf(terminal, "\033[18D\033[K");}// fflush(terminal); usleep(sleep_time);}
    fprintf(terminal, "block %3d /%3d... ", b, nblocks); fflush(terminal); usleep(sleep_time);
  
    // build block indices
    build_block_indices(b_idx, sb_idx, &b_len, &sb_len, b, prof, sb_grid);
    if (print) {
    cout << endl; cout << "-----------------" << endl;
    cout << "b = " << b << endl;
    cout << "b_len = " << b_len << endl;
    cout << "sb_len = " << sb_len << endl;
    cout << "offset = " << offset << endl;
    // cout << "sb_idx:" << endl;
    // for (int a=0; a<sb_len; a++) {
    //   cout << a << ": " << sb_idx[a] << endl;
    // }
    }

    if (print) {
    for (int i=0; i<b_len; i++) {
      cout << "i, lam, grid = " << i << ", " << eq_class[offset+i] << ", " << eig_grid[offset+i] << endl;
    }
    }

    get_block_log_len(
      block_log_prof, block_log_len, rhs_log_len,
      b_len, num_classes, offset, eig_grid, eq_class
    );

    if (print) {
      cout << endl << "block_pf:" << endl;
      lcm_mat_ep += offset;
      for (int i=0; i<b_len; i++) {
        lcm_mat_ep[i] += offset;
      }
      poly_frac_rk2_print_to_math(lcm_mat_ep, b_len, b_len, roots);
      for (int i=0; i<b_len; i++) {
        lcm_mat_ep[i] -= offset;
      }
      lcm_mat_ep -= offset;
    }

    // coefficients of LCM, block and constant term
    sys_block_info(
      lcm[b], block, const_term,
      &flag_lcm, &lcm_deg[b], lcm_orig[b], lcm_roots[b],
      &block_max_deg[b], &num_max_deg[b], &den_max_deg[b],
      b_len, sb_len,
      b_idx, sb_idx,
      lcm_mat_ep, mpc_lcm_mat_ep,
      npoles, poles,
      eta_shift,
      eta_ord,
      sol, num_cum_eig, cum_eig, sblock_log_len, rhs_log_len, log_prof
    );

    if (print) {
    cout << endl; cout << "block_max_deg = " << block_max_deg[b] << endl;
    cout << "block:" << endl; print_rk3_mpc(block, b_len, b_len, block_max_deg[b]+1);
    cout << "lcm: " << endl; print_poly(lcm[b], lcm_deg[b]);
    if (sb_len > 0) {
      cout << endl; cout << "constant term: " << endl;
      for (int lam, n=0; n<num_cum_eig; n++) {
        lam = cum_eig[n];
        cout << "lam = " << lam << endl;
        for (int l=0; l<rhs_log_len[lam]; l++) {
          cout << "l = " << l << endl;
          print_rk2_mpc(const_term[n][l], b_len, prt_ord-1);
        }
      }
    }
    // if (b == 40) exit(0);

    cout << endl; cout << "previous solutions log-profile" << endl;
    for (int n=0; n<num_classes; n++) {
      cout << "eq class = " << n << endl;
      for (int i=0; i<offset; i++) {
        // if (b == 23) log_prof[0][i] = 1;
        cout << "log-prof[" << i << "] = " << log_prof[n][i] << endl;
      }
    }

    cout << endl; cout << "RHS log-lenghts" << endl;
    for (int i=0; i<num_classes; i++) {
      cout << "eq class = " << i << ", length = " << rhs_log_len[i] << endl;
    }
    }
    
    // get log-lengths of the solution
    // num_block_eig = 0;
    for (int lam=0; lam<num_classes; lam++) {
      sol_log_len[lam] = rhs_log_len[lam] + block_log_len[lam];
      if (block_log_len[lam] > 0 && rhs_log_len[lam] > 0) {
        // num_block_eig++;
        // x - 1 + y - 2 + 1 + 1 + 1
        sol_log_len[lam] =  rhs_log_len[lam] + block_log_len[lam];
      } else if (block_log_len[lam] > 0 && rhs_log_len[lam] == 0) {
        sol_log_len[lam] = block_log_len[lam];
      } else if (block_log_len[lam] == 0 && rhs_log_len[lam] > 0) {
        sol_log_len[lam] = rhs_log_len[lam] + 1;
      }
      for (int i=0; i<b_len; i++) {
        log_prof[lam][offset+i] = sol_log_len[lam];
        block_log_prof[lam][i] = block_log_len[lam];
      }

      if (sol_log_len[lam] > b_len_max+2) {
        fprintf(
          stderr,
          "b = %d, lam = %d: sol_log_len = %d greater than b_len_max + 2 = %d\n",
          b, lam, sol_log_len[lam], b_len_max+2
        );
        exit(1);
      }
    }
    if (print) {
    cout << endl; cout << "block log-lengths" << endl;
    for (int i=0; i<num_classes; i++) {
      cout << "eq class = " << i << ", length = " << block_log_len[i] << endl; 
    }
    cout << endl; cout << "solution log-lengths" << endl;
    for (int i=0; i<num_classes; i++) {
      cout << "eq class = " << i << ", length = " << sol_log_len[i] << endl; 
    }
    }

    // cumulate list of eigenvalues
    for (int i=0; i<b_len; i++) {
      if (eig_grid[offset+i] == 1) {
        // new eigenvalue appears
        cum_eig[num_cum_eig] = eq_class[offset+i];
        num_cum_eig++;
      }
    }
    if (print) {
    cout << endl; cout << "cumulated eigenvalues" << endl;
    for (int i=0; i<num_cum_eig; i++) {
      cout << cum_eig[i] << endl;
    }
    }

    // needed for TADPOLES AND BUBBLES
    // #TBD: make try_analytic a flag decided from chart, but always
    // activate it when it is necessary (i.e. when exiting from singularity)
    // int try_analytic = 1;
    int tad_mass;
    int *bub_mass = new int[2];
    int bub_p2;
    mpc_t hcsol;
    mpc_init3(hcsol, wp2, wp2);
    mpc_t eps;
    mpc_init3(eps, wp2, wp2);
    // if (try_analytic == 1 && b_len == 1 && (sb_len == 0 || sb_len == 2)){
    // if (try_analytic == 1 && b_len == 1 && (sb_len == 0)){
    if (try_analytic == 1 && b == 0){
      //////
      // DEPENDENCY MATRIX
      //////
      int **depmat;
      malloc_rk2_tens(depmat, dim, ninvs);
      // manually set depmat
      for (int i=0; i<dim; i++) {
        for (int s=0; s<ninvs; s++) {
          depmat[i][s] = 0;
        }
      }

      // gnc_to_mpc(&eps, (ex) *epv *1.);
      mpc_set_str_rat(&eps, eps_str);
      // cout << "epv = " << (ex) *epv << endl;
      // print_mpc(&eps); cout << endl;
      if (sb_len == 0) {
        // TADPOLES
        // build solution at the matching point
        // identify mass
        for (int s=0; s<ninvs; s++) {
          if (depmat[offset][s] == 1) {
            tad_mass = s;
            // cout << "tadpoles mass = " << tad_mass << endl;
            break;
          }
        }
        for (int lam, n=0; n<num_cum_eig; n++) {
          lam = cum_eig[n];
          if (block_log_len[lam] > 0) {
            // cout << "lam = " << lam << endl;
            if (skip_inv[tad_mass] == 0) {
              // compute tadpole
              get_tadpole(
                &sol[lam][0][offset][0],
                eta_target, &PS_ini[tad_mass], &PS_fin[tad_mass], &eps,
                &inv_tmat[offset][offset], eig_list, lam, analytic_cont
              );
            } else if (skip_inv[tad_mass] == 1) {
              mpc_set(sol[lam][0][offset][0], (*solutions)[offset][0], MPFR_RNDN);
            }
            // cout << "final tapdole:" << endl; print_mpc(&sol[lam][0][offset][0]); cout << endl;
            for (int k=1; k<=eta_ord; k++) {
              mpc_set_ui(sol[lam][0][offset][k], 0, MPFR_RNDN);
            }
          } else {
            for (int k=0; k<=eta_ord; k++) {
              mpc_set_ui(sol[lam][0][offset][k], 0, MPFR_RNDN);
            }
          }
          // fill other log-powers with zero
          for (int l=1; l<sol_log_len[lam]-1; l++) {
            for (int k=0; k<=eta_ord; k++) {
              mpc_set_ui(sol[lam][l][offset][k], 0, MPFR_RNDN);
            }
          }
        }
      } else if (sb_len == 2) {
        // BUBBLES
        // build solution at the matching point
        // identify masses and squared momentum
        int count_masses = 0;
        for (int s=0; s<ninvs; s++) {
          if (depmat[offset][s] == 1) {
            if (is_mass[s] == 1) {
              bub_mass[count_masses] = s;
              count_masses++;
              // cout << "found mass: " << s << endl;
            } else if (is_mass[s] == 0) {
              // cout << "found squared momentum: " << s << endl;
              bub_p2 = s;
            }
          }
        }
        if (count_masses == 1) {
          bub_mass[1] = bub_mass[0];
        }

        for (int lam, n=0; n<num_cum_eig; n++) {
          lam = cum_eig[n];
          if (block_log_len[lam] > 0) {
            // cout << "lam = " << lam << endl;
            if (skip_inv[bub_p2] == 0) {
              // compute bubble
              get_massless_bubble(
                &sol[lam][0][offset][0],
                eta_target, &PS_ini[bub_p2], &PS_fin[bub_p2], &eps,
                &inv_tmat[offset][offset], eig_list, lam, analytic_cont
              );
            } else if (skip_inv[bub_p2] == 1) {
              mpc_set(sol[lam][0][offset][0], (*solutions)[offset][0], MPFR_RNDN);
            }
            // cout << "final bubble:" << endl; print_mpc(&sol[lam][0][offset][0]); cout << endl;
            for (int k=1; k<=eta_ord; k++) {
              mpc_set_ui(sol[lam][0][offset][k], 0, MPFR_RNDN);
            }
          } else {
            for (int k=0; k<=eta_ord; k++) {
              mpc_set_ui(sol[lam][0][offset][k], 0, MPFR_RNDN);
            }
          }
          // fill other log-powers with zero
          for (int l=1; l<sol_log_len[lam]-1; l++) {
            for (int k=0; k<=eta_ord; k++) {
              mpc_set_ui(sol[lam][l][offset][k], 0, MPFR_RNDN);
            }
          }
        }
        // exit(0);
      }
    } else {
      // alloc memory for solutions: (# eigs) x (# log powers) x (# eta powers)
      hsol = new mpc_t****[num_cum_eig];
      psol = new mpc_t****[num_cum_eig];
      for (int lam, n=0; n<num_cum_eig; n++) {
        lam = cum_eig[n];
        if (sol_log_len[lam] == 0) {
          continue;
        }
        // #2BD: for the homogeneous, allocate only for block eigs
        malloc_rk4_tens(hsol[n], sol_log_len[lam], b_len, b_len, eta_ord+1);
        init_rk4_mpc(hsol[n], sol_log_len[lam], b_len, b_len, eta_ord+1);
        malloc_rk4_tens(psol[n], sol_log_len[lam], b_len, 1, eta_ord+1);
        init_rk4_mpc(psol[n], sol_log_len[lam], b_len, 1, eta_ord+1);
      }

      // HOMOGENEOUS SOLUTION
      // cout << endl; cout << "ENTER to solve homogeneous" << endl;
      // getchar();
      zero_solve_block(hsol,
        lcm[b], block,
        b_idx, b_len, sb_len,
        eta_ord, lcm_deg[b], block_max_deg[b],
        eig_list, cum_eig, num_cum_eig,
        block_log_len, sol_log_len, num_classes);
      
      if (print) {
      cout << endl; cout << "SOLUTION OF THE HOMOGENEOUS" << endl;
      if (1) {
        for (int lam, n=0; n<num_cum_eig; n++) {
          lam = cum_eig[n];
          if (block_log_len[lam] == 0) {
            continue;
          }
          cout << "lam = " << lam << endl;
          for (int l=0; l<sol_log_len[lam]; l++) {
            cout << "l = " << l << endl;
            print_rk3_mpc(hsol[n][l], b_len, b_len, prt_ord+1);
          }
        }
      }
      }

      // PARTICULAR SOLUTION
      if (sb_len > 0) {
        // cout << endl; cout << "ENTER to solve particular" << endl;
        // getchar();
        zero_solve_block(psol,
          lcm[b], block,
          b_idx, b_len, sb_len,
          eta_ord, lcm_deg[b], block_max_deg[b],
          eig_list, cum_eig, num_cum_eig,
          block_log_len, sol_log_len, num_classes,
          true, const_term, rhs_log_len
          // ,match_in_zero, lam0, *solutions + offset
        );

        if (print) {
        cout << endl; cout << "SOLUTION OF THE PARTICULAR" << endl;
        if (1) {
          for (int lam, n=0; n<num_cum_eig; n++) {
            lam = cum_eig[n];
            if (sol_log_len[lam] == 0) {
              continue;
            }
            cout << "lam = " << lam << endl;
            for (int l=0; l<sol_log_len[lam]; l++) {
              cout << "l = " << l << endl;
              print_rk3_mpc(psol[n][l], b_len, 1, prt_ord+1);
            }
          }
        }
        }
      }

      //////
      // CONSTRAINTS
      //////
      // LOG CONSTRAINTS
      if (sb_len > 0) {
        if (print) cout << endl << "enter to SOLVING LOG CONSTRAINTS" << endl;
        // getchar();
        solve_log_constraints(
          psol,
          eta_ord, b_len, offset, num_cum_eig, cum_eig, eq_class,
          block_log_len, sol_log_len, rhs_log_len, hsol
        );

        // if (b == 20) {
        if (print) {
          cout << "PARTICULAR after log matching" << endl;
          for (int lam, n=0; n<num_cum_eig; n++) {
            lam = cum_eig[n];
            if (sol_log_len[lam] == 0) {
              continue;
            }
            cout << "lam = " << lam << endl;
            for (int l=0; l<sol_log_len[lam]; l++) {
              cout << "l = " << l << endl;
              print_rk3_mpc(psol[n][l], b_len, 1, prt_ord+1);
            }
          }
        }

      }

      // MATCHING
      if (print) cout << endl << "enter to SOLVING MATCH CONSTRAINTS" << endl;
      // getchar();
      match_sol_around_zero(
        sol,
        eta_ord, b_len, sb_len, offset, num_cum_eig, cum_eig,
        num_classes, eq_class, eig_list,
        block_log_len, sol_log_len, rhs_log_len,
        hsol, psol,
        solutions, matching_point, lam0,
        analytic_cont,
        sol, tmat, nroots, roots,
        bound_behav, mi_eig, mi_eig_num, log_prof,
        mat_ep
      );

      for (int lam, n=0; n<num_cum_eig; n++) {
        lam = cum_eig[n];
        if (sol_log_len[lam] == 0) {
          continue;
        }
        mpc_rk4_clear(hsol[n], sol_log_len[lam], b_len, b_len, eta_ord+1);
        del_rk4_tens(hsol[n], sol_log_len[lam], b_len, b_len);
        mpc_rk4_clear(psol[n], sol_log_len[lam], b_len, 1, eta_ord+1);
        del_rk4_tens(psol[n], sol_log_len[lam], b_len, 1);
      }
      delete[] hsol;
      delete[] psol;
      
      for (int lam=0; lam<num_classes; lam++) {
        delete[] block_log_prof[lam];
      }

      if (sb_len > 0) {
        // delete constant term
        for (int lam, n=0; n<num_cum_eig; n++) {
          lam = cum_eig[n];
          if (rhs_log_len[lam] > 0) {
            mpc_rk3_clear(const_term[n],rhs_log_len[lam], b_len, eta_ord+1);
            del_rk3_tens(const_term[n], rhs_log_len[lam], b_len);
          }
        }
        delete[] const_term;
      }
    }

    if (cross) {
      // EVALUATE SOLUTIONS BEYOND ZERO
      if (print) cout << endl << "evaluate beyond zero" << endl;
      malloc_rk3_tens(cross_zero, num_cum_eig, b_len, 1);
      init_rk3_mpc(cross_zero, num_cum_eig, b_len, 1);
      for (int lam, n=0; n<num_cum_eig; n++) {
        lam = cum_eig[n];
        for (int l=0; l<sol_log_len[lam]; l++) {
          sol[lam][l] += offset;
        }
      }
      evaluate_series_around_zero(cross_zero, eta_target, NULL,
        eig_list, num_cum_eig, cum_eig, sol_log_len, b_len, 1, eta_ord, num_classes, lam0,
        true, sol,
        NULL, NULL, true,
        analytic_cont);
      if (print) cout << "evaluated series around zero" << endl;
      for (int lam, n=0; n<num_cum_eig; n++) {
        lam = cum_eig[n];
        for (int l=0; l<sol_log_len[lam]; l++) {
          sol[lam][l] -= offset;
        }
      }
      if (print) cout << "solution at target:" << endl;
      for (int i=0; i<b_len; i++) {
        if (print) cout << "i = " << i << endl;
        mpc_set_ui(sol_at_target[offset+i][0], 0, MPFR_RNDN);
        for (int lam, n=0; n<num_cum_eig; n++) {
          lam = cum_eig[n];
          if (print) cout << "lam = " << lam << endl;
          if (sol_log_len[lam] == 0) {
            continue;
          }
          if (print) {
          cout << "cross zero: "; print_mpc(&cross_zero[n][i][0]); cout << endl;
          }
          mpc_add(sol_at_target[offset+i][0], sol_at_target[offset+i][0], cross_zero[n][i][0], MPFR_RNDN);
        }
        if (print) {
        print_mpc(&sol_at_target[offset+i][0]); cout << endl;
        }
      }
      mpc_rk3_clear(cross_zero, num_cum_eig, b_len, 1);
      del_rk3_tens(cross_zero, num_cum_eig, b_len);

      // evaluate_series_around_zero(cross_zero, &eta_target, psol,
      //   eig_list, num_cum_eig, cum_eig, sol_log_len, b_len, 1, eta_ord, num_classes);
      
      // for (int i=0; i<b_len; i++) {
      //   mpc_set_ui(sol_at_target[offset+i], 0, MPFR_RNDN);
      //   for (int lam, n=0; n<num_cum_eig; n++) {
      //     lam = cum_eig[n];
      //     if (sol_log_len[lam] == 0) {
      //       continue;
      //     }
      //     cout << "lam = " << endl;
      //     mpc_add(sol_at_target[offset+i], sol_at_target[offset+i], cross_zero[n][i][0], MPFR_RNDN);
      //   }
      // }

      
      // tmat += offset;
      // for (int i=0; i<b_len; i++) {
      //   tmat[i] += offset;
      // }
      
      // if (b != 1 && b != 5){
      if (0) {
        // cout << endl; cout << "limit in zero:" << endl;
        // for (int lam, n=0; n<num_cum_eig; n++) {
        //   lam = cum_eig[n];
        //   for (int l=0; l<sol_log_len[lam]-1; l++) {
        //     sol[lam][l] += offset;
        //   }
        // }
        pf_limit_in_zero(
          tmp_vec,
          sol, tmat,
          b_len+offset, lam0, tp_rank,
          roots, npoles,
          1
        );
        // for (int i=0; i<b_len; i++) {
        //   print_mpc(&tmp_vec[offset+i]); cout << endl;
        // }
        // for (int lam, n=0; n<num_cum_eig; n++) {
        //   lam = cum_eig[n];
        //   for (int l=0; l<sol_log_len[lam]-1; l++) {
        //     sol[lam][l] -= offset;
        //   }
        // }
      }

      // cout << "transformed solution at target:" << endl;
      tmat += offset;
      poly_frac_rk2_eval_value(tmat_eval, tmat, roots, eta_target, b_len, offset+b_len);
      // if (b==6) {exit(0);}
      mpc_rk2_mul_mpc_rk2_slice(
        tmp_vec, tmat_eval, sol_at_target,
        b_len, b_len+offset, 0
      );
      // for (int i=0; i<b_len; i++) {
      //   // mpc_set(solutions[0][i][0], tmp_vec[i], MPFR_RNDN);
      //   print_mpc(&tmp_vec[i]); cout << endl;
      // }

      // for (int i=0; i<b_len; i++) {
      //   tmat[i] -= offset;
      // }
      tmat -= offset;
    }
    else {
      // cout << "COMPUTE LIMIT IN ZERO (BY BLOCK)" << endl;
      // pf_limit_in_zero_block(
      //   tmp_vec,
      //   pf_sol,
      //   sol, tmat,
      //   offset, b_len, lam0, tp_rank,
      //   roots, npoles,
      //   rad_exp
      // );
      // for (int i=0; i<b_len; i++) {
      //   mpc_set(solutions[0][offset+i][0], tmp_vec[i], MPFR_RNDN);
      // }
    }
    if (print) {
    cout << endl; cout << "enter to SOLUTION" << endl;
    // getchar();
    for (int lam, n=0; n<num_cum_eig; n++) {
      lam = cum_eig[n];
      cout << "lam = " << lam << endl;
      for (int l=0; l<sol_log_len[lam]-1; l++) {
        cout << "l = " << l << endl;
        print_rk2_mpc(sol[lam][l] + offset, b_len, prt_ord+1);
      }
    }
    }

    // cumulate log-lengths
    for (int lam, i=0; i<num_cum_eig; i++) {
      lam = cum_eig[i];
      if (sblock_log_len[lam] < sol_log_len[lam]) {
        sblock_log_len[lam] = sol_log_len[lam];
      }
    }
    // cout << endl; cout << "cumulated log-lengths" << endl;
    // for (int i=0; i<num_classes; i++) {
    //   cout << "eq class = " << i << ", length = " << sblock_log_len[i] << endl; 
    // }

    // FREE
    for (int i=0; i<b_len; i++) {
      for (int j=0; j<b_len; j++) {
        mpc_rk1_clear(mpc_lcm_mat_ep[b_idx[i]][b_idx[j]][0], block_max_deg[b]+1);
        delete[] mpc_lcm_mat_ep[b_idx[i]][b_idx[j]][0];
      }
      for (int j=0; j<sb_len; j++) {
        mpc_rk1_clear(mpc_lcm_mat_ep[b_idx[i]][sb_idx[j]][0], num_max_deg[b]+1);
        delete[] mpc_lcm_mat_ep[b_idx[i]][sb_idx[j]][0];
        mpc_rk1_clear(mpc_lcm_mat_ep[b_idx[i]][sb_idx[j]][1], den_max_deg[b]+1);
        delete[] mpc_lcm_mat_ep[b_idx[i]][sb_idx[j]][1];
      }
    }

    for (int i=0; i<b_len; i++) {
      for (int j=0; j<b_len; j++) {
        mpc_rk1_clear(block[i][j], block_max_deg[b]+1);
        delete[] block[i][j];
      }
      delete[] block[i];
    }
    delete[] block;

    mpc_rk1_clear(lcm[b], lcm_deg[b]+1);
    delete[] lcm[b];

    offset += b_len;
    // if (b == 18) {fflush(stdout); exit(0);}
    // if (b == 18) {break;}
  }
  fprintf(terminal, "\033[18D\033[K"); fflush(terminal); usleep(sleep_time);

  if (debug) {
  //////
  // EXPORT SOLUTIONS
  //////
  string solmath;
  FILE *fptr;
  char *tmp_filepath = NULL;
  char *sol_debug_name = (char*) malloc(MAX_PATH_LEN*sizeof(char));
  // char *filepath_debug_sol = strdup("debug/sol_lam");
  for (int lam, n=0; n<num_cum_eig; n++) {
    lam = cum_eig[n];
    for (int l=0; l<sol_log_len[lam]-1; l++) {
      mpc_rk2_to_math(&solmath, sol[lam][l], dim, eta_ord);
      snprintf(
        sol_debug_name, MAX_PATH_LEN, "%s%d%s%d",
        (char*)"sollam", lam, (char*)"l", l
      );
      make_dir((char*)"debug/");
      join_path(&tmp_filepath, (char*)"debug/", sol_debug_name);
      join_path(&tmp_filepath, tmp_filepath, (char*)".txt");
      cout << "export solution to " << tmp_filepath << endl;
      fptr = fopen(tmp_filepath, "w");
      fprintf(fptr, "%s", sol_debug_name);
      fprintf(fptr, "%s", (char*)"=");
      fprintf(fptr, "%s", solmath.c_str());
      fclose(fptr);
    }
  }
  }
  
  // #debug
  int nmi = dim, kplus = 1;
  struct poly_frac **pf_sol_dbg;
  malloc_rk2_tens(pf_sol_dbg, num_classes, dim);
  if (print) {
  cout << endl << "BEHAVIOUR OF THE MIs" << endl;
  // cout << "tp_rank = " << tp_rank << endl;
  for (int lam, n=0; n<num_cum_eig; n++) {
    lam = cum_eig[n];
    // cout << "lam = " << lam << endl;
    // print_mpc(&eig_list[n]); cout << endl;
    if (sol_log_len[lam] == 0) {
      // cout << "skip" << endl;
      continue;
    }
    for (int i=0; i<nmi; i++) {
      // cout << "i = " << i << endl;
      poly_frac_build(&pf_sol_dbg[lam][i]);
      if (log_prof[lam][i] == 0) {
        poly_frac_set_zero(&pf_sol_dbg[lam][i]);
        continue;
      }
      // cout << "input poly:" << endl;
      // print_poly(sol[lam][0][i], tp_rank+1+kplus);
      poly_frac_set_coeffs(&pf_sol_dbg[lam][i], sol[lam][0][i], tp_rank+1+kplus);
      pf_sol_dbg[lam][i].den_deg = 0;
      pf_sol_dbg[lam][i].nroots = nroots;
      pf_sol_dbg[lam][i].mults = new int[nroots];
      for (int k=0; k<nroots; k++) {
        pf_sol_dbg[lam][i].mults[k] = 0;
      }
      // poly_frac_print(&pf_sol_dbg[lam][i]);
    }

    // apply matrix
    // cout << "apply matrix" << endl;
    poly_frac_rk2_mul_pf_rk1(pf_sol_dbg[lam], tmat, pf_sol_dbg[lam], nmi, nmi, roots);
  }

  for (int m=0; m<dim; m++) {
    for (int lam, n=0; n<num_cum_eig; n++) {
      lam = cum_eig[n];
      if (sol_log_len[lam] == 0) {
        continue;
      }
      cout << "m, lam = " << m << ", " << lam << endl;
      poly_frac_print(&pf_sol_dbg[lam][m]);
    }
  }
  }

  if (cross) {
    // TRANSFORM BACK MI AT TARGET IN THE ORIGINAL GAUGE
    if (print) {
    cout << "TRANSFORM BACK MI AT TARGET IN THE ORIGINAL GAUGE" << endl;
    cout << "transformed solution at target:" << endl;
    for (int i=0; i<dim; i++) {
      cout << "MI n. " << i << endl;
      print_mpc(&sol_at_target[i][0]);
      cout << endl;
    }
    cout << endl;
    }
    poly_frac_rk2_eval_value(tmat_eval, tmat, roots, eta_target, dim, dim);
    mpc_rk2_mul_mpc_rk2_slice(
      tmp_vec, tmat_eval, sol_at_target,
      dim, dim, 0
    );
    // cout << "untransformed solution at target:" << endl;
    // the copy tmp_vec is unnecessary, we could avoid it if inside
    // pf_limit_in_zero we perform the rk2_mul_rk1 with a slice on the output
    for (int i=0; i<dim; i++) {
      mpc_set(solutions[0][i][0], tmp_vec[i], MPFR_RNDN);
      // print_mpc(&solutions[0][i][0]); cout << endl;
    }
    mpc_rk2_clear(sol_at_target, dim, 1);
    del_rk2_tens(sol_at_target, dim);
  } else {
    // COMPUTE LIMIT IN ZERO AND TRANSFORM BACK TO ORIGINAL GAUGE
    if (print) cout << "COMPUTE LIMIT IN ZERO" << endl;
    pf_limit_in_zero(
      tmp_vec,
      sol, tmat,
      dim, lam0, tp_rank,
      roots, npoles,
      rad_exp
    );
    // the copy tmp_vec is unnecessary, we could avoid it if inside
    // pf_limit_in_zero we perform the rk2_mul_rk1 with a slice on the output
    for (int i=0; i<dim; i++) {
      mpc_set(solutions[0][i][0], tmp_vec[i], MPFR_RNDN);
    }
  }

  // FREE
  mpfr_clear(tmpfr);
  delete[] cum_eig;
  delete[] block_log_len;
  for (int b=0; b<nblocks; b++) {
    mpc_rk1_clear(lcm_orig[b], lcm_deg[b]);
    delete[] lcm_orig[b];
  }
  delete[] lcm;
  delete[] lcm_orig;
  delete[] tmp_vec;
  mpc_rk2_clear(tmat_eval, dim, dim);
  del_rk2_tens(tmat_eval, dim);
  del_rk3_tens(mpc_lcm_mat_ep, dim, dim);
  poly_frac_rk2_free(lcm_mat_ep, dim, dim);
  del_rk2_tens(lcm_mat_ep, dim);
  mpc_rk4_clear(sol, num_classes, b_len_max+2, dim, eta_ord+1);
  del_rk4_tens(sol, num_classes, b_len_max+2, dim);

  return;
}


void center_around_pole(
  // OUTPUT
  struct poly_frac **sh_pfmat, mpc_t *sh_roots,
  // INPUT
  struct poly_frac **pfmat, mpc_t *roots, int label,
  int nroots, int dim1, int dim2
) {
  // check whether no shift is necessary (shift by zero)
  if (label == 0) {
    // copy matrix
    if (sh_pfmat != pfmat) {
      poly_frac_rk2_set_pf_rk2(sh_pfmat, pfmat, dim1, dim2);
    }

    // copy roots
    if (sh_roots != roots) {
      for (int k=0; k<nroots; k++) {
        mpc_set(sh_roots[k], roots[k], MPFR_RNDN);
      }
    }
    
    return;
  }

  //////
  // SHIFT MATRIX
  //////
  // cout << endl; cout << "shift matrix..." << endl;
  poly_frac_rk2_shift(
    sh_pfmat, pfmat, &roots[label], label,
    dim1, dim2
  );
  // poly_frac_rk2_print(sh_pfmat, dim, dim);

  //////
  // SHIFT POLES
  //////
  // cout << "shift poles..." << endl;
  // SWAP NULL AND SINGULAR ROOTS
  // the singular root becomes an exact zero and has to be moved in
  // poisition zero, therefore there is no need to do anything;
  // the null root shifts and becomes the opposite of the singular root,
  // therefore we just need to change a sign
  mpc_neg(sh_roots[label], roots[label], MPFR_RNDN);
  // SHIFT ALL THE OTHERS
  for (int k=1; k<nroots; k++) {
    if (k == label) {continue;}
    mpc_sub(sh_roots[k], roots[k], roots[label], MPFR_RNDN);
  }
  // cout << "shifted roots:" << endl;
  // print_poly(sh_roots, nroots-1);
}


void set_analytic_cont_sign(
  // OUTPUT
  int *analytic_cont,
  // INPUT
  mpc_t *m1, mpc_t *m2,
  int branch_deg, mpc_t *branch_poly
) {
  mpc_t zm, zp;
  mpc_init3(zm, wp2, wp2);
  mpc_init3(zp, wp2, wp2);
  mpfr_t argzm, argzp, argdiff;
  mpfr_init2(argzm, wp2);
  mpfr_init2(argzp, wp2);
  mpfr_init2(argdiff, wp2);
  mpfr_t mpfr_pi, mpfr_2pi;
  mpfr_init2(mpfr_pi, wp2);
  mpfr_init2(mpfr_2pi, wp2);
  mpfr_const_pi(mpfr_pi, MPFR_RNDN);
  mpfr_mul_ui(mpfr_2pi, mpfr_pi, 2, MPFR_RNDN);

  // evaluate branch polynomial in the matching points
  poly_eval(&zm, branch_poly, branch_deg, m1);
  poly_eval(&zp, branch_poly, branch_deg, m2);

  // cout << "zm = "; print_mpc(&zm); cout << endl;
  // cout << "zp = "; print_mpc(&zp); cout << endl;

  // mpc_set_si_si(zm, 1, 0, MPFR_RNDN);
  // mpc_arg(argzm, zm, MPFR_RNDN);
  // cout << "test argzm: "; mpfr_out_str(stdout, 10, 0, argzm, MPFR_RNDN); cout << endl;
  // mpc_set_si_si(zm, 0, 1, MPFR_RNDN);
  // mpc_arg(argzm, zm, MPFR_RNDN);
  // cout << "test argzm: "; mpfr_out_str(stdout, 10, 0, argzm, MPFR_RNDN); cout << endl;
  // mpc_set_si_si(zm, -1, 0, MPFR_RNDN);
  // mpc_arg(argzm, zm, MPFR_RNDN);
  // cout << "test argzm: "; mpfr_out_str(stdout, 10, 0, argzm, MPFR_RNDN); cout << endl;
  // mpc_set_si_si(zm, 0, -1, MPFR_RNDN);
  // mpc_arg(argzm, zm, MPFR_RNDN);
  // cout << "test argzm: "; mpfr_out_str(stdout, 10, 0, argzm, MPFR_RNDN); cout << endl;

  // #TBD:
  // 1. change mpfr_sign_within_tol in order to accept as input an arbitrary tolerance;
  // 2. pass as input the tolerance divided by the modulus of the complex number
  //    (error propagation for the argument of a complex number);

  // compute argument in such a way that arg is in [0, 2*Pi[, while mpc's arg is in ]-Pi, Pi]
  mpc_arg(argzm, zm, MPFR_RNDN);
  // cout << "initial argzm: "; mpfr_out_str(stdout, 10, 0, argzm, MPFR_RNDN); cout << endl;
  if (mpfr_sign_within_tol(argzm) < 0) {
    // mpfr_add(tmpfr, argzm, mpfr_pi, MPFR_RNDN);
    // if (mpfr_sign_within_tol(tmpfr) > 0) {
      mpfr_add(argzm, argzm, mpfr_2pi, MPFR_RNDN);
    // }
  // } else {
  //   mpfr_sub(tmpfr, argzm, mpfr_pi, MPFR_RNDN);
  //   } if (mpfr_sign_within_tol(tmpfr) == 0) {
  //   mpfr_neg(argzm, argzm, MPFR_RNDN);
  //   }
  }
  mpc_arg(argzp, zp, MPFR_RNDN);
  // cout << "initial argzp: "; mpfr_out_str(stdout, 10, 0, argzp, MPFR_RNDN); cout << endl;
  if (mpfr_sign_within_tol(argzp) < 0) {
    // mpfr_add(tmpfr, argzp, mpfr_pi, MPFR_RNDN);
    // if (mpfr_sign_within_tol(tmpfr) > 0) {
      mpfr_add(argzp, argzp, mpfr_2pi, MPFR_RNDN);
    // }
  // } else {
  //   mpfr_sub(tmpfr, argzp, mpfr_pi, MPFR_RNDN);
  //   } if (mpfr_sign_within_tol(tmpfr) == 0) {
  //     mpfr_neg(argzp, argzp, MPFR_RNDN);
  //   }
  }
  // cout << "arg(zm) = "; mpfr_out_str(stdout, 10, 0, argzm, MPFR_RNDN); cout << endl;
  // cout << "arg(zp) = "; mpfr_out_str(stdout, 10, 0, argzp, MPFR_RNDN); cout << endl;

  mpfr_sub(argzm, argzm, argzp, MPFR_RNDN);
  if (mpfr_sign_within_tol(argzm) < 0) {
    *analytic_cont = -1;
  } else if (mpfr_sign_within_tol(argzm) > 0) {
    *analytic_cont = 1;
  }
}


void propagate_infty(
  // IN-OUT
  mpc_t ***solutions,
  // INPUT
  char *eps_str,
  mpc_t *bound, mpc_t *path,
  int dim, poly_frac **pfmat,
  int nroots, mpc_t *roots,
  int **bound_behav, int **mi_eig, int *mi_eig_num,
  int nloops, int eta_ord,
  FILE *logfptr, FILE *terminal
) {
  int print = 0;

  fprintf(terminal, "infinity: "); fflush(terminal); usleep(sleep_time);

  //////
  // ENLARGE TOLERANCE
  //////
  double wp2_rel_decr = 0.80;
  double wp2_rel_decr_orig = wp2_rel_decr;
  double wp2_rel_decr_prune = 0.80;
  double wp2_rel_decr_prune_orig = wp2_rel_decr_prune;

  wp2_rel_decr *= wp2_rel_decr_orig;
  mpfr_tol_enlarge(wp2_rel_decr);
  // mpfr_to_gnc(&gnc_tol, &mpfr_tol);  // #uncomment-for-ginac
  // cout << "tol: "; mpfr_out_str(stdout, 10, 0, mpfr_tol, MPFR_RNDN); cout << endl;

  int wp_bin = - mpfr_log2_int(mpfr_tol);

  // check matrix in one point
  // mpc_t etapoint;
  // mpc_init3(etapoint, wp2, wp2);
  // mpc_set_ui(etapoint, 17, MPFR_RNDN);
  // mpc_t val;
  // mpc_init3(val, wp2, wp2);
  // cout << "{" << endl;
  // for (int i=0; i<dim; i++) {
  //   cout << "{" << endl;
  //   for (int j=0; j<dim; j++) {
  //     poly_frac_eval_value(&val, &pfmat[i][j], roots, &etapoint);
  //     // cout << "i, j, val = " << i << ", " << j << ", "; print_mpc(&val); cout << endl;
  //     mpfr_out_str(stdout, 10, 0, mpc_realref(val), MPFR_RNDN);
  //     if (j < dim-1) {
  //       cout << ",";
  //     }
  //     cout << endl;
  //   }
  //   cout << "}";
  //   if (i < dim-1) {
  //     cout << ",";
  //   }
  //   cout << endl;
  // }
  // cout << "}" << endl;

  //////
  // TRANSFORM eta -> 1/eta
  //////
  // dbg = 1;
  // fprintf(logfptr, "transform eta -> 1/eta\n"); fflush(logfptr);
  poly_frac **pfmat_infty;
  malloc_rk2_tens(pfmat_infty, dim, dim);
  poly_frac_rk2_build(pfmat_infty, dim, dim);
  int zero_label_infty;
  mpc_t *roots_infty = new mpc_t[nroots];
  init_rk1_mpc(roots_infty, nroots);
  // cout << "wp_bin = " << wp_bin << endl;
  rel_err_poly_frac_rk2_infty(
  // poly_frac_rk2_infty(
    pfmat_infty, roots_infty, &zero_label_infty,
    pfmat, roots,
    dim, dim, nroots
    , wp_bin
  );
  if (print) {
  cout << "matrix after eta -> 1/eta:" << endl;
  poly_frac_rk2_print(pfmat_infty, dim, dim);
  // poly_frac_rk2_print_to_math(pfmat_infty, dim, dim, roots_infty);
  }

  //////
  // FACTOR OUT BEHAVIOR (DEs FOR REGULAR SERIES)
  //////
  // fprintf(logfptr, "factor out behavior (DE for regular series)\n"); fflush(logfptr);
  mpc_t *pow_infty = new mpc_t[dim];
  init_rk1_mpc(pow_infty, dim);
  poly_frac pftmp;
  poly_frac_build(&pftmp);
  mpc_t nloop_eps_neg;
  mpc_init3(nloop_eps_neg, wp2, wp2);
  mpc_set_str_rat(&nloop_eps_neg, eps_str);  // eps
  mpc_mul_si(nloop_eps_neg, nloop_eps_neg, -nloops, MPFR_RNDN);  // -nloop*eps
  for (int i=0; i<dim; i++) {
    // cout << "i = " << i << endl;
    // cout << "input:" << endl;
    // poly_frac_print(&pfmat[i][i]);
    for (int j=0; j<dim; j++) {
      if (j == i) {continue;}
      // cout << "j = " << j << endl;
      // inv_tmat.DE.tmat = eta^-a.DE.eta^a
      // poly_frac_print(&pfmat_infty[0][i][j]);
      poly_frac_mul_sym_pow(
        &pfmat_infty[i][j],
        &pfmat_infty[i][j],
        - bound_behav[i][0] + bound_behav[j][0]
      );
      // poly_frac_print(&pfmat_infty[0][i][j]);
    }
    
    // cout << "diagonal:" << endl;
    // cout << "bound_behav = " << bound_behav[i][0] << endl;
    // -inv_tmat.D[tmat,eta] = -a/eta
    mpc_add_si(pow_infty[i], nloop_eps_neg, -bound_behav[i][0], MPFR_RNDN);  // -a(eps) = -nloop*eps - bound_behav
    // cout << "-nloop*eps - bound_behav = "; print_mpc(&pow_infty[i]); cout << endl;
    poly_frac_set_mpc(&pftmp, &pow_infty[i], nroots);
    poly_frac_mul_sym_pow(&pftmp, &pftmp, -1);  // -a(eps)/eta
    // cout << "1st addend:" << endl;
    // poly_frac_print(&pfmat_infty[i][i]);
    // cout << "2nd addend:" << endl;
    // poly_frac_print(&pftmp);
    rel_err_poly_frac_add_pf(
    // poly_frac_add_pf(
      &pfmat_infty[i][i],
      &pfmat_infty[i][i], &pftmp,
      roots_infty,
      NULL
      , wp_bin
    );
    // cout << "result:" << endl;
    // poly_frac_print(&pfmat_infty[i][i]);
  }
  // dbg = 0;

  int **bound_behav_reg;
  malloc_rk2_tens(bound_behav_reg, dim, 1);
  for (int i=0; i<dim; i++) {
    bound_behav_reg[i][0] = 0;
  }

  if (print) {
  cout << "matrix at infinity:" << endl;
  poly_frac_rk2_print(pfmat_infty, dim, dim);
  // poly_frac_rk2_print_to_math(pfmat_infty, dim, dim, roots_infty);
  }
  // cout << "roots at infity:" << endl;
  // print_poly(roots_infty, nroots-1);

  // // check matrix in one point
  // cout << "eval matrix at infinity:" << endl;
  // mpc_t etapoint;
  // mpc_init3(etapoint, wp2, wp2);
  // mpc_set_ui(etapoint, 17, MPFR_RNDN);
  // mpc_t val;
  // mpc_init3(val, wp2, wp2);
  // cout << "{" << endl;
  // for (int i=0; i<dim; i++) {
  //   cout << "{" << endl;
  //   for (int j=0; j<dim; j++) {
  //     poly_frac_eval_value(&val, &pfmat_infty[i][j], roots_infty, &etapoint);
  //     // cout << "i, j, val = " << i << ", " << j << ", "; print_mpc(&val); cout << endl;
  //     mpfr_out_str(stdout, 10, 0, mpc_realref(val), MPFR_RNDN);
  //     if (j < dim-1) {
  //       cout << ",";
  //     }
  //     cout << endl;
  //   }
  //   cout << "}";
  //   if (i < dim-1) {
  //     cout << ",";
  //   }
  //   cout << endl;
  // }
  // cout << "}" << endl;

  //////
  // SOLVE AT INFINITY
  //////
  fprintf(terminal, "singular: "); fflush(terminal); usleep(sleep_time);

  // prune
  wp2_rel_decr_prune *= wp2_rel_decr_prune_orig;
  // poly_frac_rk2_prune(pfmat_infty, dim, dim, wp2_rel_decr_prune);
  poly_frac_rk2_normal(pfmat_infty, dim, dim);

  // wp2 *= 0.90;

  // BLOCK PROFILE
  int nblocks, **prof, **sb_grid;
  poly_frac_rk2_find_profile(
    &prof, &sb_grid, &nblocks,
    pfmat_infty, dim
  );
  // cout << endl; cout << "BLOCK PROFILE:" << endl;
  // int tmp_offset = 0;
  // for (int b=0; b<nblocks; b++) {
  //   // select block
  //   int b_len = prof[b][1] - prof[b][0] + 1;
  //   // cout << "-------------------------" << endl;
  //   cout << "b = " << b << ", ";
  //   cout << "b_len = " << b_len << ", ";
  //   cout << "offset = " << tmp_offset << endl;
  //   tmp_offset += b_len;
  // }
  if (print) {
  cout << endl; cout << "BLOCK GRID:" << endl;
  for (int b=0; b<nblocks; b++) {
    for (int sb=0; sb<b; sb++) {
      cout << "b, sb = " << b << ", " << sb << ": ";
      cout << sb_grid[b][sb] << endl;
    }
    cout << endl;
  }
  }

  // NORMALIZE
  struct poly_frac **tmat, **inv_tmat;
  malloc_rk2_tens(tmat, dim, dim);
  poly_frac_rk2_build(tmat, dim, dim);
  malloc_rk2_tens(inv_tmat, dim, dim);
  poly_frac_rk2_build(inv_tmat, dim, dim);
  int num_classes;
  mpc_t *eig_list;
  int *eq_class = new int[dim];
  int *eig_grid = new int[dim];
  fprintf(logfptr, "\nnormalize at infinity...\n");
  pf_NormalizeMat(
    tmat, inv_tmat,
    &num_classes, eq_class, &eig_list, eig_grid,
    pfmat_infty,
    dim, nblocks, prof, sb_grid,
    roots_infty, nroots,
    terminal
  );
  if (print) {
  cout << "tmat=";
  poly_frac_rk2_print_to_math(tmat, dim, dim, roots_infty);
  // cout << "tmat:" << endl;
  // poly_frac_rk2_print(tmat, dim, dim);
  cout << "invtmat=";
  poly_frac_rk2_print_to_math(inv_tmat, dim, dim, roots_infty);
  // cout << "inv tmat:" << endl;
  // poly_frac_rk2_print(inv_tmat, dim, dim);
  cout << "mat=";
  poly_frac_rk2_print_to_math(pfmat_infty, dim, dim, roots_infty);

  cout << endl << "EIGENVALUES:" << endl;
  print_eigenvalues(
    dim, num_classes, eig_grid, eq_class,
    eig_list
  );
  }

  // cout << "mat:" << endl;
  // poly_frac_rk2_print(pfmat_infty, dim, dim);
  // exit(0);
  
  // // compute block profile for normalized matrix
  // int nblocks_norm, **prof_norm, **sb_grid_norm;
  // // cout << endl; cout << "computing block profile..." << endl;
  // poly_frac_rk2_find_profile(
  //   &prof_norm, &sb_grid_norm, &nblocks_norm, 
  //   pfmat_infty, dim
  // );
  // cout << endl; cout << "BLOCK PROFILE after normalization:" << endl;
  // tmp_offset = 0;
  // for (int b=0; b<nblocks_norm; b++) {
  //   // select block
  //   int b_len = prof_norm[b][1] - prof_norm[b][0] + 1;
  //   // cout << "-------------------------" << endl;
  //   cout << "b = " << b << ", ";
  //   cout << "b_len = " << b_len << ", ";
  //   cout << "offset = " << tmp_offset << endl;
  //   tmp_offset += b_len;
  // }

  // update sub-diagonal grid for normalized matrix
  int **sb_grid_norm;
  generate_sub_diag_grid(
    &sb_grid_norm,
    prof, nblocks, pfmat_infty
  );

  if (print) {
  cout << endl; cout << "BLOCK GRID after normalization:" << endl;
  for (int b=0; b<nblocks; b++) {
    for (int sb=0; sb<b; sb++) {
      cout << "b, sb = " << b << ", " << sb << ": ";
      cout << sb_grid_norm[b][sb] << endl;
    }
    cout << endl;
  }
  }

  // SOLVE
  mpc_t bound_pt, target_pt, eta_pow;
  mpc_init3(bound_pt, wp2, wp2);
  mpc_init3(target_pt, wp2, wp2);
  mpc_init3(eta_pow, wp2, wp2);

  // numeric epv = (numeric) eps_str;
  mpc_set_ui(bound_pt, 0, MPFR_RNDN);
  mpc_pow_si(target_pt, path[0], -1,  MPFR_RNDN);
  for (int i=0; i<dim; i++) {
    // mpc_init3(solutions[0][i][0], wp2, wp2);
    mpc_set(solutions[0][i][0], bound[i], MPFR_RNDN);
  }

  int ninvs = 1;
  int *is_mass = new int[1]; is_mass[0] = 0;
  int *skip_inv = new int[1]; skip_inv[0] = 0;
  mpc_t *PS_ini = new mpc_t[1];
  mpc_init3(PS_ini[0], wp2, wp2);
  mpc_set_ui(PS_ini[0], 0, MPFR_RNDN);
  mpc_t *PS_fin = new mpc_t[1];
  mpc_init3(PS_fin[0], wp2, wp2);
  mpc_set_ui(PS_fin[0], 1, MPFR_RNDN);
  fprintf(logfptr, "singular propagation at infinity...\n");
  solve_zero(
    solutions,
    dim, pfmat_infty,
    tmat, inv_tmat,
    nroots, roots_infty,
    &bound_pt,
    nblocks, prof, sb_grid_norm, eta_ord,
    num_classes, eig_list, eq_class, eig_grid,
    roots_infty, nroots,
    1, &target_pt, 1,
    is_mass, skip_inv, ninvs, PS_ini, PS_fin, eps_str,
    0,
    bound_behav_reg, mi_eig, mi_eig_num,
    terminal
  );
  // solve_zero(
  //   solutions,
  //   dim, pfmat_infty,
  //   tmat, inv_tmat,
  //   nroots, roots_infty,
  //   &bound_pt,
  //   nblocks, prof, sb_grid, eta_ord,
  //   num_classes, eig_list, eq_class, eig_grid,
  //   roots_infty, nroots,
  //   1, &target_pt, 1,
  //   NULL, is_mass, skip_inv, ninvs, PS_ini, PS_fin, eps_str,
  //   0,
  //   bound_behav_reg, mi_eig, mi_eig_num
  // );

  //////
  // MULTIPLY BACK BEHAVIOR
  //////
  fprintf(logfptr, "\nresult of singular propagation:\n");
  for (int i=0; i<dim; i++) {
    fprintf(logfptr, "i = %d: ", i);

    // multiply power behaviour
    mpc_pow(eta_pow, path[0], pow_infty[i], MPFR_RNDN);
    // cout << "pow_infty = "; print_mpc(&pow_infty[i]); cout << endl;
    // cout << "eta_pow = "; print_mpc(&eta_pow); cout << endl;
    mpc_mul(solutions[0][i][0], solutions[0][i][0], eta_pow, MPFR_RNDN);

    // print
    mpc_out_str(logfptr, 10, 0, solutions[0][i][0], MPFR_RNDN); fprintf(logfptr, "\n");
  }

  fprintf(terminal, "\033[10D\033[K"); fflush(terminal); usleep(sleep_time);
  fprintf(terminal, "\033[10D\033[K"); fflush(terminal); usleep(sleep_time);

  // FREE
  poly_frac_rk2_free(pfmat_infty, dim, dim);
  del_rk2_tens(pfmat_infty, dim);
  poly_frac_rk2_free(tmat, dim, dim);
  del_rk2_tens(tmat, dim);
  poly_frac_rk2_free(inv_tmat, dim, dim);
  del_rk2_tens(inv_tmat, dim);

}


void propagate_along_path(
  // OUTPUT
  mpc_t **sol_at_eps,
  // IN-OUT
  mpc_t ***solutions,
  // INPUT
  int dim, int eta_ord,
  int ep, char *eps_str,
  int neta_values, mpc_t *path, int *path_tags, int nsings,
  struct poly_frac **pfmat, int nroots, mpc_t *roots, int *sing_lab,
  int nblocks, int **prof, int **sb_grid,
  int nbranches, int *branch_deg, mpc_t **branch_poly, int **branch_sing_lab,
  int ninvs, mpc_t *PS_ini, mpc_t *PS_fin, char **symbols,
  int *is_mass, int *skip_inv,
  int **bound_behav, int **mi_eig, int *mi_eig_num,
  FILE *logfptr, FILE *terminal
) {
  int print = 0;
  mpc_t tmpc;
  mpc_init3(tmpc, wp2, wp2);

  // double wp2_rel_decr = 0.98;
  // double wp2_rel_decr = 0.95;
  // double wp2_rel_decr = 0.93;
  // double wp2_rel_decr = 0.90;
  // double wp2_rel_decr = 0.85;
  // double wp2_rel_decr = 0.81;
  // double wp2_rel_decr = 0.80;
  double wp2_rel_decr = 0.77;
  // double wp2_rel_decr = 0.70;
  // double wp2_rel_decr = 0.64;
  double wp2_rel_decr_orig = wp2_rel_decr;
  // double wp2_rel_decr_prune = 0.90;
  double wp2_rel_decr_prune = 0.80;
  // double wp2_rel_decr_prune = 0.77;
  // double wp2_rel_decr_prune = 0.70;
  double wp2_rel_decr_prune_orig = wp2_rel_decr_prune;

  // needed for analytic continuation
  int *analytic_cont;
  analytic_cont = new int[nsings];
  for (int k=0; k<nsings; k++) {
    analytic_cont[k] = 0;
  }

  // needed for shift
  struct poly_frac **sh_pfmat;
  malloc_rk2_tens(sh_pfmat, dim, dim);
  poly_frac_rk2_build(sh_pfmat, dim, dim);
  mpc_t *sh_roots;
  sh_roots = new mpc_t[nroots];
  init_rk1_mpc(sh_roots, nroots);
  mpc_set_ui(sh_roots[0], 0, MPFR_RNDN);
  // needed to check shift

  // needed for matrix normalization
  struct poly_frac **tmat, **inv_tmat;
  malloc_rk2_tens(tmat, dim, dim);
  poly_frac_rk2_build(tmat, dim, dim);
  malloc_rk2_tens(inv_tmat, dim, dim);
  poly_frac_rk2_build(inv_tmat, dim, dim);
  int num_classes;
  mpc_t *eig_list;
  int *eq_class = new int[dim];
  int *eig_grid = new int[dim];

  // needed for singular propagation
  mpc_t bound_pt, target_pt;
  mpc_init3(bound_pt, wp2, wp2);
  mpc_init3(target_pt, wp2, wp2);
  int crossl, crossr;
  mpc_t m1, m2;
  mpc_init3(m1, wp2, wp2);
  mpc_init3(m2, wp2, wp2);

  int solve_tag, count_sings = 0, nreg_steps = 1;

  // fprintf(logfptr, "{\n");
  for (int et=0; et<neta_values; et++) {
    if (et  > 0) {fprintf(terminal, "\033[23D\033[K");}// fflush(terminal); usleep(sleep_time);}
    fprintf(terminal, "path point %3d /%3d... ", et, neta_values); fflush(terminal); usleep(sleep_time);

    // if (et != neta_values - 1) continue;
    // FILE *bc_ptr = fopen("sol_at_eta-last_box1.txt", "r");
    // for (int i=0; i<dim; i++) {
    //   mpc_inp_str(solutions[0][i][0], bc_ptr, 0, 10, MPFR_RNDN);
    // }
    
    // if (et == 24) exit(0);
    fprintf(logfptr, "\n");
    fprintf(logfptr, "-----------------------------------------\n");
    fprintf(logfptr, "et = %d\n", et);
    fprintf(logfptr, "point = ");
    mpc_out_str(logfptr, 10, 0, path[et], MPFR_RNDN);
    fprintf(logfptr, "\n");
    fprintf(logfptr, "tag = %d\n", path_tags[et]);

    // print invariants at eta point
    // if (path_tags[et] == 0) {
    //   continue;
    // }
    // fprintf(logfptr, "{\n");
    // for (int s=0; s<ninvs; s++) {
    //   mpc_sub(tmpc, PS_fin[s], PS_ini[s], MPFR_RNDN);
    //   mpc_fma(tmpc, path[et], tmpc, PS_ini[s], MPFR_RNDN);
    //   fprintf(logfptr, "  %s", symbols[s+1]);
    //   if (is_mass[s]) {
    //     // fprintf(logfptr, "2");
    //     mpc_sqr(tmpc, tmpc, MPFR_RNDN);
    //   }
    //   fprintf(logfptr, " -> ");
    //   mpfr_out_str(logfptr, 10, 0, mpc_realref(tmpc), MPFR_RNDN);
    //   cout << " + I*(";
    //   mpfr_out_str(logfptr, 10, 0, mpc_imagref(tmpc), MPFR_RNDN);
    //   cout << ")";
    //   if (s < ninvs-1) {
    //     fprintf(logfptr, ",");
    //   }
    //   fprintf(logfptr, "\n");
    // }
    // fprintf(logfptr, "}");
    // if (et < neta_values-1) {
    //   fprintf(logfptr, ",");
    // }
    // fprintf(logfptr, "\n");
    fflush(logfptr);

    // process tag
    if (path_tags[et] == 0) {
      solve_tag = 0;
    } else if (path_tags[et] == 1) {
      if (et > 0 && (path_tags[et-1] == 2 || path_tags[et-1] == 1)) {
        solve_tag = 1;
      } else {
        continue;
      }
    } else if (path_tags[et] == 2) {
      nreg_steps++;
      continue;
    } else {
      fprintf(stderr, "tag not recognized");
      exit(EXIT_FAILURE);
    }
    // cout << "solve tag = " << solve_tag << endl;

    if (solve_tag == 1) {
      //////
      // REGULAR PROPAGATION
      //////
      fprintf(logfptr, "\nregular propagation...\n");
      fprintf(logfptr, "n. steps: %d\n", nreg_steps);
      fflush(logfptr);
      propagate_regular(
        solutions,
        nreg_steps+1, path+et-nreg_steps,
        dim, pfmat,
        nroots, roots,
        nblocks, prof, sb_grid, eta_ord
      );
      nreg_steps = 1;
      // print
      fprintf(logfptr, "\nresult of regular propagation:\n");
      for (int i=0; i<dim; i++) {
        fprintf(logfptr, "i = %d: ", i);
        mpc_out_str(logfptr, 10, 0, solutions[0][i][0], MPFR_RNDN); fprintf(logfptr, "\n");
        // print_mpc(&solutions[0][i][0]); cout << endl;
      }
      fflush(logfptr);
      // exit(0);
    } else if (solve_tag == 0) {
      fprintf(terminal, "singular: "); fflush(terminal); usleep(sleep_time);
      // update tolerance
      // cout << endl; cout << "update tol" << endl;
      // mpfr_mul_d(mpfr_tol, mpfr_tol, 1e30, MPFR_RNDN);
      // gnc_tol *= 1e30;
      wp2_rel_decr *= wp2_rel_decr_orig;
      // cout << "orig tol: "; mpfr_out_str(stdout, 10, 0, mpfr_tol, MPFR_RNDN); cout << endl;
      // cout << "wp2_rel_decr = " << wp2_rel_decr << endl;
      mpfr_tol_enlarge(wp2_rel_decr);
      // mpfr_to_gnc(&gnc_tol, &mpfr_tol);  // #uncomment-for-ginac
      // cout << "tol: "; mpfr_out_str(stdout, 10, 0, mpfr_tol, MPFR_RNDN); cout << endl;
      
      // int wp_bin = - mpfr_log2_int(mpfr_tol);

      // cout << endl; cout << "sing. point = "; print_mpc(&roots[sing_lab[count_sings]]); cout << endl;
      //////
      // SINGULAR PROPAGATION
      //////
      // # DECIDE TYPE OF CROSSING
      if (et == 0) {
        crossl = 0;
        crossr = 1;
      } else if (et == neta_values-1) {
        crossl = 1;
        crossr = 0;
      } else {
        crossl = 1;
        crossr = 1;
      }

      // # SET ANALYTIC CONTINUATION
      // cout << endl; cout << "set analytic continuation sign..." << endl;
      // cout << "sing lab: " << sing_lab[count_sings] << endl;
      if (crossl == 1) {
        mpc_set(m1, path[et-1], MPFR_RNDN);
      } else {
        // #suspicious
        mpc_sub(m1, path[et+1], path[et], MPFR_RNDN);
        mpc_sub(m1, path[et], m1, MPFR_RNDN);          
      }
      if (crossr == 1) {
        mpc_set(m2, path[et+1], MPFR_RNDN);
      } else {
        // #suspicious
        mpc_sub(m2, path[et], path[et-1], MPFR_RNDN);
        mpc_add(m2, path[et], m2, MPFR_RNDN);          
      }
      // cout << "m1 = "; print_mpc(&m1); cout << endl;
      // cout << "m2 = "; print_mpc(&m2); cout << endl;

      // if (cross == 1) {
      //   mpc_set(m1, path[et-1], MPFR_RNDN);
      //   mpc_set(m2, path[et+1], MPFR_RNDN);
      // } else {
      //   if (et == 0) {
      //     mpc_sub(m1, path[et+1], path[et], MPFR_RNDN);
      //     mpc_sub(m1, path[et], m1, MPFR_RNDN);
      //     mpc_set(m2, path[et+1], MPFR_RNDN);
      //   } else if (et == neta_values - 1) {
      //     mpc_set(m1, path[et-1], MPFR_RNDN);
      //     mpc_sub(m2, path[et], path[et-1], MPFR_RNDN);
      //     mpc_add(m2, path[et], m2, MPFR_RNDN);
      //   }
      // }
      // find associated branch point
      int found_branch = 0;
      for (int b=0; b<nbranches; b++) {
        for (int k=1; k<=branch_deg[b]; k++) {
          // cout << "branch sing lab: " << branch_sing_lab[b][k] << endl;
          if (branch_sing_lab[b][k] == sing_lab[count_sings]) {
            // this is a branch point
            found_branch = 1;
            // cout << "BRANCH POINT" << endl;
            set_analytic_cont_sign(
              &analytic_cont[count_sings],
              &m1, &m2, branch_deg[b], branch_poly[b]
            );
            break;
          } 
        }
      }
      if (found_branch == 0) {
        analytic_cont[count_sings] = 1;
      }
      // analytic_cont[count_sings] *= -1;
      // cout << "analytic continuation sign: " << analytic_cont[count_sings] << endl;

      // # SHIFT MATRIX AND POLES
      // cout << "unshifted matrix:" << endl;
      // poly_frac_rk2_print(pfmat, dim, dim);
      // cout << "shift by " << endl; print_mpc(&roots[sing_lab[count_sings]]); cout << endl;
      rel_err_center_around_pole(
        sh_pfmat, sh_roots,
        pfmat, roots, sing_lab[count_sings],
        nroots, dim, dim,
        wp2*(0.95)
      );
      if (print) {
      cout << "shifted matrix:" << endl;
      // poly_frac_rk2_print(sh_pfmat, dim, dim);
      poly_frac_rk2_print_to_math(sh_pfmat, dim, dim, sh_roots);
      }
      
      // shift boundary and target points
      // cout << "shift boundary and target points..." << endl;
      if (crossl == 1) {
        mpc_sub(bound_pt, path[et-1], roots[sing_lab[count_sings]], MPFR_RNDN);
      } else {
        mpc_set_ui(bound_pt, 0, MPFR_RNDN);
      }
      if (crossr == 1) {
        mpc_sub(target_pt, path[et+1], roots[sing_lab[count_sings]], MPFR_RNDN);
      }

      // prune
      wp2_rel_decr_prune *= wp2_rel_decr_prune_orig;
      // poly_frac_rk2_prune(sh_pfmat, dim, dim, wp2_rel_decr_prune);
      poly_frac_rk2_normal(sh_pfmat, dim, dim);

      // # NORMALIZE MATRIX
      fprintf(logfptr, "\nnormalize matrix...\n"); fflush(logfptr);
      // cout << "input matrix:" << endl;
      // poly_frac_rk2_print(sh_pfmat, dim, dim);
      // mpfr_mul_d(mpfr_tol, mpfr_tol, 1e10, MPFR_RNDN);
      // cout << "tol: "; mpfr_out_str(stdout, 10, 0, mpfr_tol, MPFR_RNDN); cout << endl;
      // if (et == 28) {
      // pf_check_NormalizeMat(
      //   tmat, inv_tmat,
      //   num_classes, eq_class, eig_list, eig_grid,
      //   sh_pfmat,
      //   dim, nblocks, prof, sb_grid,
      //   sh_roots, nroots
      // );
      // exit(0);
      // }
      // cout << "mat:" << endl;
      // poly_frac_rk2_print(sh_pfmat, dim, dim);
      // poly_frac_rk2_print_to_math(sh_pfmat, dim, dim, sh_roots);
      // cout << "tol: "; mpfr_out_str(stdout, 10, 0, mpfr_tol, MPFR_RNDN); cout << endl;
      pf_NormalizeMat(
        tmat, inv_tmat,
        &num_classes, eq_class, &eig_list, eig_grid,
        sh_pfmat,
        dim, nblocks, prof, sb_grid,
        sh_roots, nroots,
        terminal
      );
      if (print) {
      cout << "tmat=";
      poly_frac_rk2_print_to_math(tmat, dim, dim, sh_roots);
      // cout << "tmat:" << endl;
      // poly_frac_rk2_print(tmat, dim, dim);
      cout << "invtmat=";
      poly_frac_rk2_print_to_math(inv_tmat, dim, dim, sh_roots);
      // cout << "inv tmat:" << endl;
      // poly_frac_rk2_print(inv_tmat, dim, dim);
      cout << "mat=";
      poly_frac_rk2_print_to_math(sh_pfmat, dim, dim, sh_roots);

      cout << endl << "EIGENVALUES:" << endl;
      print_eigenvalues(
        dim, num_classes, eig_grid, eq_class,
        eig_list
      );
      }

      // cout << "mat:" << endl;
      // poly_frac_rk2_print(sh_pfmat, dim, dim);

      // // compute block profile for normalized matrix
      // int nblocks_norm, **prof_norm, **sb_grid_norm;
      // // cout << endl; cout << "computing block profile..." << endl;
      // poly_frac_rk2_find_profile(
      //   &prof_norm, &sb_grid_norm, &nblocks_norm, 
      //   sh_pfmat, dim
      // );
      // cout << endl; cout << "BLOCK PROFILE after normalization:" << endl;
      // tmp_offset = 0;
      // for (int b=0; b<nblocks_norm; b++) {
      //   // select block
      //   int b_len = prof_norm[b][1] - prof_norm[b][0] + 1;
      //   // cout << "-------------------------" << endl;
      //   cout << "b = " << b << ", ";
      //   cout << "b_len = " << b_len << ", ";
      //   cout << "offset = " << tmp_offset << endl;
      //   tmp_offset += b_len;
      // }

      // update sub-diagonal grid for normalized matrix
      int **sb_grid_norm;
      generate_sub_diag_grid(
        &sb_grid_norm,
        prof, nblocks, sh_pfmat
      );

      if (print) {
      cout << endl; cout << "BLOCK GRID after normalization:" << endl;
      for (int b=0; b<nblocks; b++) {
        for (int sb=0; sb<b; sb++) {
          cout << "b, sb = " << b << ", " << sb << ": ";
          cout << sb_grid_norm[b][sb] << endl;
        }
        cout << endl;
      }
      }
      // # SOLVE
      fprintf(logfptr, "singular propagation...\n"); fflush(logfptr);
      int try_analytic = 0;
      // if (et == 0) {try_analytic = 1;}
      solve_zero(
        solutions,
        dim, sh_pfmat,
        tmat, inv_tmat,
        nroots, sh_roots,
        &bound_pt,
        nblocks, prof, sb_grid_norm, eta_ord,
        num_classes, eig_list, eq_class, eig_grid,
        sh_roots, nroots,
        crossr, &target_pt, analytic_cont[count_sings],
        is_mass, skip_inv, ninvs, PS_ini, PS_fin, eps_str,
        try_analytic,
        bound_behav, mi_eig, mi_eig_num,
        terminal
      );
      // solve_zero(
      //   solutions,
      //   dim, sh_pfmat,
      //   tmat, inv_tmat,
      //   nroots, sh_roots,
      //   &bound_pt,
      //   nblocks, prof, sb_grid, eta_ord,
      //   num_classes, eig_list, eq_class, eig_grid,
      //   sh_roots, nroots,
      //   crossr, &target_pt, analytic_cont[count_sings],
      //   depmat, is_mass, skip_inv, ninvs, PS_ini, PS_fin, eps_str,
      //   try_analytic,
      //   bound_behav, mi_eig, mi_eig_num
      // );
      // print
      fprintf(logfptr, "\nresult of singular propagation:\n");
      for (int i=0; i<dim; i++) {
        fprintf(logfptr, "i = %d: ", i);
        mpc_out_str(logfptr, 10, 0, solutions[0][i][0], MPFR_RNDN); fprintf(logfptr, "\n");
        // print_mpc(&solutions[0][i][0]); cout << endl;
      }
      fflush(logfptr);
      // exit(0);

      count_sings++;
      // // update tolerance // DEBUG
      // mpfr_mul_d(mpfr_tol, mpfr_tol, 1e-30, MPFR_RNDN);
      // gnc_tol *= 1e-30;
      fprintf(terminal, "\033[10D\033[K"); fflush(terminal); usleep(sleep_time);
    }
  }
  fprintf(terminal, "\033[23D\033[K"); fflush(terminal); usleep(sleep_time);
  // copy to output
  for (int i=0; i<dim; i++) {
    mpc_set(sol_at_eps[i][ep], solutions[0][i][0], MPFR_RNDN);
  }

  fprintf(logfptr, "}\n"); fflush(logfptr);

  // FREE
  delete[] eq_class;
  delete[] eig_grid;
  poly_frac_rk2_free(tmat, dim, dim);
  del_rk2_tens(tmat, dim);
  poly_frac_rk2_free(inv_tmat, dim, dim);
  del_rk2_tens(inv_tmat, dim);
  poly_frac_rk2_free(sh_pfmat, dim, dim);
  del_rk2_tens(sh_pfmat, dim);
  mpc_rk1_clear(sh_roots, nroots);
  delete[] sh_roots;
}


void propagate_eps(
  // OUTPUT
  mpc_t **sol_at_eps,
  // INPUT
  int ep, int exit_sing, int nloops,
  poly_frac ***pfmat,
  int *zero_label, int *nroots, mpc_t **roots,
  int *neta_values, mpc_t **path, int **path_tags, int *nsings, int **sing_lab,
  mpc_t ***solutions, int dim, int eta_ord,
  int eps_num, char **eps_str,
  int ninvs, mpc_t *PS_ini, mpc_t *PS_fin, char **symbols,
  int *is_mass, int *skip_inv,
  int nbranches, int *branch_deg, mpc_t **branch_poly, int **branch_sing_lab,
  // int **bound_behav, int **mi_eig, int *mi_eig_num,
  int gen_bound, char *filepath_bound, mpc_t **bound,
  char *filepath_bound_build, char *filepath_bound_behav,
  // char *filepath_matrix, char *filepath_roots, // char *filepath_branch_sing_lab,
  char *filepath_path, char *filepath_path_tags, char *filepath_sol, char* dir_partial,
  char *file_ext, FILE *logfptr, int opt_write, int opt_checkpoint,
  FILE *terminal
) {
  int print = 0;
  char tmp_filepath[200];

  // struct poly_frac **pfmat;
  // malloc_rk2_tens(pfmat, dim, dim);
  // poly_frac_rk2_build(pfmat, dim, dim);
  // int nroots, zero_label;
  // mpc_t *roots;
  FILE *boundfptr;
  // ex gnc_tol_orig = gnc_tol;  // #uncomment-for-ginac
  mpfr_t mpfr_tol_orig;
  mpfr_init2(mpfr_tol_orig, wp2);
  mpfr_set(mpfr_tol_orig, mpfr_tol, MPFR_RNDN);

  int num_region_classes, **bound_behav, **mi_eig, *mi_eig_num;
  if (exit_sing) {
    //////
    // LOAD BOUNDARY BEHAVIOUR
    //////
    // cout << endl; cout << "read bound behav from file..." << endl;
    bound_behav_from_file(
      filepath_bound_behav,
      &bound_behav, &mi_eig, &mi_eig_num, &num_region_classes,
      dim
    );
  }

  //////
  // LOAD POLY_FRAC DE
  //////
  char tmp_filepath_pfmat[MAX_PATH_LEN];
  char tmp_filepath_roots[MAX_PATH_LEN];
  if (opt_checkpoint == 1) {
    snprintf(tmp_filepath_pfmat, sizeof(tmp_filepath_pfmat), "%s%s%d%s", dir_partial, "pfmat", ep, ".txt");
    snprintf(tmp_filepath_roots, sizeof(tmp_filepath_roots), "%s%s%d%s", dir_partial, "roots", ep, ".txt");
  } else {
    snprintf(tmp_filepath_pfmat, sizeof(tmp_filepath_pfmat), "%s%s%d%s", filepath_cache, "pfmat", ep, ".txt");
    snprintf(tmp_filepath_roots, sizeof(tmp_filepath_roots), "%s%s%d%s", filepath_cache, "roots", ep, ".txt");
  }
  // // MATRIX
  // cout << endl; cout << "reading poly_frac DE from file " << tmp_filepath_pfmat << endl;
  // malloc_rk2_tens(pfmat[ep], dim, dim);
  // poly_frac_rk2_build(pfmat[ep], dim, dim);
  // poly_frac_rk2_from_file(tmp_filepath_pfmat, pfmat[ep], dim, dim);
  // // ROOTS
  // cout << "reading mpc_t roots from file " << tmp_filepath_roots << endl;
  // nroots[ep] = count_lines(tmp_filepath_roots) - 1;
  // roots[ep] = new mpc_t[nroots[ep]];
  // init_rk1_mpc(roots[ep], nroots[ep]);
  // int_rk0_mpc_rk1_from_file(tmp_filepath_roots, roots[ep], nroots[ep], &zero_label[ep]);

  //////
  // LOAD BOUNDARY
  //////
  if (gen_bound == 0) {
    snprintf(tmp_filepath, sizeof(tmp_filepath), "%s%d%s", filepath_bound, ep, file_ext);
    fprintf(logfptr, "\nloading boundary from %s\n", tmp_filepath); fflush(logfptr);
    boundfptr = fopen(tmp_filepath, "r");
    if (boundfptr == NULL) {
      printf("Could not open file %s\n", tmp_filepath);
      exit(1);
    }
    for (int i=0; i<dim; i++) {
      // mpc_init3(solutions[0][i][0], wp2, wp2);
      mpc_inp_str(solutions[0][i][0], boundfptr, 0, 10, MPFR_RNDN);
    }
    fclose(boundfptr);
  } else {
    for (int i=0; i<dim; i++) {
      // mpc_init3(solutions[0][i][0], wp2, wp2);
      mpc_set(solutions[0][i][0], bound[ep][i], MPFR_RNDN);
    }
  }
  if (print) {
  cout << "BOUNDARY:" << endl;
  for (int i=0; i<dim; i++) {
    cout << "BC n." << i << ": "; print_mpc(&solutions[0][i][0]); cout << endl;
  }
  }

  if (exit_sing == -1) {
    //////
    // EXIT SINGULARITY AT INFINITY
    //////
    fprintf(logfptr, "\nEXIT SINGULARITY AT INFINITY...\n"); fflush(logfptr);
    propagate_infty(
      solutions,
      eps_str[ep],
      bound[ep], path[ep],
      dim, pfmat[ep],
      nroots[ep], roots[ep],
      bound_behav, mi_eig, mi_eig_num,
      nloops, eta_ord,
      logfptr, terminal
    );
  }

  //////
  // FIND BLOCK-PROFILE
  //////
  int nblocks, **prof, **sb_grid;   
  // cout << endl; cout << "computing block profile..." << endl;
  poly_frac_rk2_find_profile(
    &prof, &sb_grid, &nblocks,
    pfmat[ep], dim
    );
  // cout << endl; cout << "BLOCK PROFILE:" << endl;
  // int tmp_offset = 0;
  // for (int b=0; b<nblocks; b++) {
  //   // select block
  //   int b_len = prof[b][1] - prof[b][0] + 1;
  //   // cout << "-------------------------" << endl;
  //   cout << "b = " << b << ", ";
  //   cout << "b_len = " << b_len << ", ";
  //   cout << "offset = " << tmp_offset << endl;
  //   tmp_offset += b_len;
  // }
  if (print) {
  cout << endl; cout << "BLOCK GRID:" << endl;
  for (int b=0; b<nblocks; b++) {
    for (int sb=0; sb<b; sb++) {
      cout << "b, sb = " << b << ", " << sb << ": ";
      cout << sb_grid[b][sb] << endl;
    }
    cout << endl;
  }
  }

  //////
  // PROPAGATION
  //////
  fprintf(logfptr, "\nPROPAGATION...\n"); fflush(logfptr);
  // exit(0);
  propagate_along_path(
    sol_at_eps,
    solutions,
    dim, eta_ord,
    ep, eps_str[ep],
    neta_values[ep], path[ep], path_tags[ep], nsings[ep],
    pfmat[ep], nroots[ep], roots[ep], sing_lab[ep],
    nblocks, prof, sb_grid,
    nbranches, branch_deg, branch_poly, branch_sing_lab,
    ninvs, PS_ini, PS_fin, symbols,
    is_mass, skip_inv,
    bound_behav, mi_eig, mi_eig_num,
    logfptr, terminal
  );
  fprintf(logfptr, "\nresult of ep = %d\n", ep);
  for (int i=0; i<dim; i++) {
    // cout << "i = " << i << ": "; print_mpc(&solutions[0][i][0]); cout << endl;
    fprintf(logfptr, "MI n. %d: ", i); mpc_out_str(logfptr, 10, 0, sol_at_eps[i][ep], MPFR_RNDN); fprintf(logfptr, "\n");
  }
  fflush(logfptr);

  // SAVE RESULT TO FILE 
  if (opt_write > 0 || opt_write == -2) {
    snprintf(tmp_filepath, sizeof(tmp_filepath), "%s%d%s", filepath_sol, ep, file_ext);
    FILE *fptr = fopen(tmp_filepath, "w");
    for (int i=0; i<dim; i++) {
      mpc_out_str(fptr, 10, 0, sol_at_eps[i][ep], MPFR_RNDN); fprintf(fptr, "\n");
    }
    fclose(fptr);
  }

  mpfr_set(mpfr_tol, mpfr_tol_orig, MPFR_RNDN);
  // gnc_tol = gnc_tol_orig;  // #uncomment-for-ginac
  // wp2 = wp2_orig;
  // if (ep == 2) exit(0); // #dbg

  // // FREE
  // poly_frac_rk2_free(pfmat[ep], dim, dim);
  // del_rk2_tens(pfmat[ep], dim);
  // mpc_rk1_clear(roots[ep], nroots[ep]);
  // delete[] roots[ep];

}

