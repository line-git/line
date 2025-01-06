#include <iostream>
#include "mpc.h"

using namespace std;

#include "utils.h"
#include "tensor_utils.h"
extern "C" {
  #include "jordan.h"
}

#include "global_vars.h"


void jordan_triv_1(
  // OUTPUT
  mpc_t **jor, mpc_t **t_mat, int **ein_pos, int *num_ein,
  // INPUT
  mpc_t *ac, int N, mpc_t *prop, mpc_t *eig, int first_non_zero,
  int print
) {
  mpc_t tmpc; mpc_init3(tmpc, wp2, wp2);
  mpfr_t tmpfr; mpfr_init2(tmpfr, wp2);
  mpfr_t norm; mpfr_init2(norm, wp2);

  *num_ein = N;
  *ein_pos = (int*) malloc(*num_ein*sizeof(int));
  for (int i=0; i<*num_ein; i++) {
    (*ein_pos)[i] = i;
  }

  //////
  // JORDAN MATRIX
  //////
  ac += first_non_zero*N;

  if (!eig) {
    // get the only non-zero eigenvalue
    mpc_init3(*eig, wp2, wp2);
    mpc_set(*eig, ac[first_non_zero], MPFR_RNDN);
    for (int i=0; i<N-1-first_non_zero; i++) {
      mpc_fma(*eig, prop[first_non_zero+i], ac[first_non_zero+i+1], *eig, MPFR_RNDN);
    }
  }

  // set all elements to zero except for the (0,0) element
  for (int j=1; j<N; j++) {
    mpc_set_ui(jor[0][j], 0, MPFR_RNDN);
  }
  for (int i=1; i<N; i++) {
    for (int j=0; j<N; j++) {
      mpc_set_ui(jor[i][j], 0, MPFR_RNDN);
    }
  }
  // set (0,0) to the only non-null eigenvalue
  mpc_set(jor[0][0], *eig, MPFR_RNDN);

  //////
  // TRANSFORMATION MATRIX
  //////
  // 1st column
  mpc_norm(norm, ac[0], MPFR_RNDN);
  for (int i=1; i<N; i++) {
    mpc_norm(tmpfr, ac[i], MPFR_RNDN);
    mpfr_add(norm, norm, tmpfr, MPFR_RNDN);
  }
  mpfr_sqrt(norm, norm, MPFR_RNDN);
  // printf("norm = "); mpfr_out_str(stdout, 10, 0, norm, MPFR_RNDN); printf("\n");
  for (int i=0; i<N; i++) {
    mpc_div_fr(t_mat[i][0], ac[i], norm, MPFR_RNDN);
  }

  // columns associated with null columns
  if (print) printf("columns associated with null columns\n");
  for (int j=1; j<=first_non_zero; j++) {
    for (int i=0; i<N; i++) {
      if (i == j-1) {
        mpc_set_ui(t_mat[i][j], 1, MPFR_RNDN);
      } else {
        mpc_set_ui(t_mat[i][j], 0, MPFR_RNDN);
      }
    }
  }

  // remaining colums
  if (print) printf("remaining columns\n");
  for (int j=first_non_zero+1; j<N; j++) {
    for (int i=0; i<first_non_zero; i++) {
      mpc_set_ui(t_mat[i][j], 0, MPFR_RNDN);
    }

    mpc_norm(norm, prop[j-1], MPFR_RNDN);
    mpfr_add_ui(norm, norm, 1, MPFR_RNDN);
    mpfr_sqrt(norm, norm, MPFR_RNDN);
    mpc_div_fr(t_mat[first_non_zero][j], prop[j-1], norm, MPFR_RNDN);
    mpc_neg(t_mat[first_non_zero][j], t_mat[first_non_zero][j], MPFR_RNDN);

    for (int i=first_non_zero+1; i<N; i++) {
      if (i == j) {
        mpc_set_ui(t_mat[i][j], 1, MPFR_RNDN);
        mpc_div_fr(t_mat[i][j], t_mat[i][j], norm, MPFR_RNDN);
      } else {
        mpc_set_ui(t_mat[i][j], 0, MPFR_RNDN);
      }
    }
  }
  
  return;
}


void jordan_triv_2(
  // OUTPUT
  mpc_t **jor, mpc_t **t_mat, int **ein_pos, int *num_ein,
  // INPUT
  mpc_t *ac, int N, mpc_t *prop, int first_non_zero, int last_non_zero_row,
  int *is_zero,
  int print
) {
  mpc_t tmpc; mpc_init3(tmpc, wp2, wp2);
  mpfr_t tmpfr; mpfr_init2(tmpfr, wp2);
  mpfr_t norm; mpfr_init2(norm, wp2);

  *num_ein = N-1;
  *ein_pos = (int*) malloc(*num_ein*sizeof(int));
  for (int i=0; i<*num_ein; i++) {
    (*ein_pos)[i] = i;
  }

  //////
  // JORDAN MATRIX
  //////
  ac += first_non_zero*N;

  if (print) {
  printf("column:\n");
  for (int i=0; i<N; i++) {
    print_mpc(&ac[i]); printf("\n");
  }
  }

  // set all elements to zero except for the (N-2,N-1) element
  for (int i=0; i<N-2; i++) {
    for (int j=0; j<N; j++) {
      mpc_set_ui(jor[i][j], 0, MPFR_RNDN);
    }
  }
  for (int j=0; j<N-1; j++) {
    mpc_set_ui(jor[N-2][j], 0, MPFR_RNDN);
  }
  for (int j=0; j<N; j++) {
    mpc_set_ui(jor[N-1][j], 0, MPFR_RNDN);
  }

  // set (N-2,N-1) to one
  mpc_set_ui(jor[N-2][N-1], 1, MPFR_RNDN);

  //////
  // TRANSFORMATION MATRIX
  //////
  // set everything to zero
  for (int i=0; i<N; i++) {
    for (int j=0; j<N; j++) {
      mpc_set_ui(t_mat[i][j], 0, MPFR_RNDN);
    }
  }

  // find null proportionality coefficients  
  int *prop_is_zero = new int[N-1];
  for (int j=first_non_zero; j<N-1; j++) {
    if (mpc_lessthan_tol(prop[j])) {
      prop_is_zero[j] = 1;
    } else {
      prop_is_zero[j] = 0;
    }
    if (print) printf("prop_is_zero[%d] = %d\n", j, prop_is_zero[j]);
  }

  // find first non-zero prop coeff from the right
  int zero_right;
  for (zero_right=N-2; zero_right>=first_non_zero; zero_right--) {
    if (!prop_is_zero[zero_right]) {
      break;
    }
  }
  if (print) printf("zero_right = %d\n", zero_right);

  // NULL COLUMNS ON THE LEFT
  for (int j=0; j<first_non_zero-1; j++) {
    mpc_set_ui(t_mat[j+1][j], 1, MPFR_RNDN);
  }
  if (first_non_zero > 0 && first_non_zero <= zero_right) {
    mpc_norm(norm, prop[first_non_zero], MPFR_RNDN);
    mpfr_add_ui(norm, norm, 1, MPFR_RNDN);
    mpfr_sqrt(norm, norm, MPFR_RNDN);
    mpc_div_fr(
      t_mat[first_non_zero][first_non_zero-1],
      prop[first_non_zero], norm,
      MPFR_RNDN
    );
    mpc_set_si(t_mat[first_non_zero+1][first_non_zero-1], -1, MPFR_RNDN);
    mpc_div_fr(
      t_mat[first_non_zero+1][first_non_zero-1],
      t_mat[first_non_zero+1][first_non_zero-1], norm,
      MPFR_RNDN
    );
  }

  // NON-NULL COLUMN
  for (int p, j=first_non_zero; j<zero_right; j++) {
    if (print) printf("j = %d\n", j);
    // find first non-zero prop coeff
    for (p=j+1; p<N-1; p++) {
      if (prop_is_zero[p] == 0) {
        break;
      }
    }
    if (print) printf("p = %d\n", p);

    if (prop_is_zero[j]) {
      // place only 1
      if (print) printf("place only 1\n");
      mpc_set_ui(t_mat[j+1][j], 1, MPFR_RNDN);
    } else {
      // place both
      if (print) printf("place both\n");
      mpc_norm(norm, prop[p], MPFR_RNDN);
      mpc_norm(tmpfr, prop[j], MPFR_RNDN);
      mpfr_add(norm, norm, tmpfr, MPFR_RNDN);
      mpfr_sqrt(norm, norm, MPFR_RNDN);
      mpc_div_fr(t_mat[j+1][j], prop[p], norm, MPFR_RNDN);
      mpc_div_fr(t_mat[p+1][j], prop[j], norm, MPFR_RNDN);
      mpc_neg(t_mat[p+1][j], t_mat[p+1][j], MPFR_RNDN);
    }
  }

  int zero_right_plus = last_non_zero_row-1;
  if (zero_right_plus < zero_right) {
    zero_right_plus = zero_right;
  }

  // NULL COLUMNS ON THE RIGHT
  mpc_norm(norm, ac[0], MPFR_RNDN);
  for (int i=1; i<=zero_right; i++) {
    if (print) printf("i = %d\n", i);
    mpc_norm(tmpfr, ac[i], MPFR_RNDN);
    mpfr_add(norm, norm, tmpfr, MPFR_RNDN);
    if (print) {printf("cum norm2 = "); print_mpc(&ac[i]); printf("\n");}
  }
  for (int j=zero_right; j<zero_right_plus; j++) {
    if (print) printf("j = %d\n", j);
    mpc_norm(tmpfr, ac[j+1], MPFR_RNDN);
    mpfr_add(norm, norm, tmpfr, MPFR_RNDN);
    if (print) {
    printf("cum norm2 = ");
    mpfr_out_str(stdout, 10, 0, norm, MPFR_RNDN); printf("\n");
    }
    if (is_zero[j+2] == 1 && j+1 != zero_right) {
      mpc_set_ui(t_mat[j+2][j], 1, MPFR_RNDN);
      continue;
    }
    mpfr_sqrt(tmpfr, norm, MPFR_RNDN);
    for (int i=0; i<=j+1; i++) {
      if (print) printf("i = %d\n", i);
      mpc_div_fr(t_mat[i][j], ac[i], tmpfr, MPFR_RNDN);
      mpc_neg(t_mat[i][j], t_mat[i][j], MPFR_RNDN);
      if (print) {print_mpc(&t_mat[i][j]); printf("\n");}
    }
  }

  // deal with the case where last elements of the column are zero
  for (int j=zero_right_plus; j<N-2; j++) {
    mpc_set_ui(t_mat[j+2][j], 1, MPFR_RNDN);
  }

  // SECOND TO LAST COLUMN
  ac -= first_non_zero*N;
  for (int i=0; i<=last_non_zero_row; i++) {
    mpc_set(t_mat[i][N-2], ac[i+(zero_right+1)*N], MPFR_RNDN);
  }
  
  // LAST COLUMN
  mpc_set_ui(t_mat[zero_right+1][N-1], 1, MPFR_RNDN);

  return;
}


void jordan_wrp(
  // OUTPUT
  mpc_t **jor, mpc_t **t_mat, int **ein_pos, int *num_ein,
  // INPUT
  mpfr_t *ar, mpfr_t *ai, int N, int nm, int debug
) {
  int print;
  double wp2_rel_decr = 0.7;
  if (debug < 0 || debug >= 80) {
    print = 1;
  } else {
    print = 0;
  }

  //////
  // CHECK FOR TRIVIAL JORDAN
  //////
  if (print) printf("CHECK FOR TRIVIAL JORDAN\n");
  // fill complex matrix
  mpc_t *ac = (mpc_t*) malloc(N*nm*sizeof(mpc_t));
  for (int i=0; i<N*nm; i++) {
    mpc_init3(ac[i], wp2, wp2);
    mpc_set_fr_fr(ac[i], ar[i], ai[i], MPFR_RNDN);
  }

  // check if columns are all proporional to each other
  mpfr_t mpfr_tol_orig; mpfr_init2(mpfr_tol_orig, wp2);
  mpfr_set(mpfr_tol_orig, mpfr_tol, MPFR_RNDN);
  mpfr_tol_enlarge(wp2_rel_decr);
  // cout << "JORDAN WRP: mpfr_tol_orig = "; mpfr_out_str(stdout, 10, 0, mpfr_tol_orig, MPFR_RNDN); cout << endl;
  // cout << "JORDAN WRP: mpfr_tol = "; mpfr_out_str(stdout, 10, 0, mpfr_tol, MPFR_RNDN); cout << endl;
  int triv_jor = 1;
  int first_non_zero, first_non_zero_row, last_non_zero_row;
  int *is_zero;
  mpc_t *prop;
  mpc_t eig; mpc_init3(eig, wp2, wp2);
  if (nm != N || N == 1) {
    triv_jor = 0;
  } else {
    prop = (mpc_t*) malloc((N-1)*sizeof(mpc_t));

    // find null elements
    is_zero = new int[N*N];
    for (int j=0; j<N; j++) {
      for (int i=0; i<N; i++) {
        // print_mpc(&ac[i+N*j]); printf("\n");
        if (mpc_lessthan_tol(ac[i+N*j])) {
          is_zero[i+N*j] = 1;
        } else {
          is_zero[i+N*j] = 0;
        }
        // cout << is_zero[i+N*j] << endl;
      }
    }

    // identify first non-zero column
    int zero_col;
    first_non_zero = 0;
    for (int j=0; j<N; j++) {
      zero_col = 1;
      for (int i=0; i<N; i++) {
        if (is_zero[i+N*j] == 0) {
          zero_col = 0;
          break;
        }
      }
      if (zero_col) {
        first_non_zero++;
        continue;
      } else {
        break;
      }
    }
    if (print) printf("first non-zero col = %d\n", first_non_zero);

    if (first_non_zero == N) {
      // INPUT IS ZERO
      if (print) printf("INPUT IS ZERO\n");
      for (int i=0; i<N; i++) {
        for (int j=0; j<N; j++) {
          mpc_set_ui(jor[i][j], 0, MPFR_RNDN);
          mpc_set_ui(t_mat[i][j], 0, MPFR_RNDN);
        }
        mpc_set_ui(t_mat[i][i], 1, MPFR_RNDN);
      }

      *num_ein = N;
      *ein_pos = (int*) malloc(*num_ein*sizeof(int));
      for (int i=0; i<*num_ein; i++) {
        (*ein_pos)[i] = i;
      }
      mpfr_set(mpfr_tol, mpfr_tol_orig, MPFR_RNDN);
      // cout << "JORDAN WRP: return: mpfr_tol = "; mpfr_out_str(stdout, 10, 0, mpfr_tol, MPFR_RNDN); cout << endl;
      return;
    }

    // check for proportionality starting from first non-zero column
    if (print) printf("check for proportionality starting from first non-zero column\n");
    for (int i=first_non_zero; i<N-1; i++) {
      if (print) printf("i = %d\n", i);
      mpc_init3(prop[i], wp2, wp2);
      if (print) {
      printf("1st vec:\n");
      print_rk1_mpc(ac+first_non_zero*N, N);
      printf("2nd vec:\n");
      print_rk1_mpc(ac+(i+1)*N, N);
      }
      if (
        !mpc_rk1_prop(
            &prop[i],
            ac+first_non_zero*N, ac+(i+1)*N, N,
            is_zero+first_non_zero*N, is_zero+(i+1)*N
        )
      ) {
        triv_jor = 0;
        break;
      }
      if (print) {printf("prop[%d] = ", i); print_mpc(&prop[i]); printf("\n");}
    }
    
    // find the first non-zero row
    first_non_zero_row = 0;
    int zero_row;
    for (int i=0; i<N; i++) {
      zero_row = 1;
      for (int j=0; j<N; j++) {
        if (!is_zero[i+j*N]) {
          zero_row = 0;
          break;
        }
      }
      if (zero_row) {
        first_non_zero_row++;
        continue;
      } else {
        break;
      }
    }
    if (print) printf("first non-zero row = %d\n", first_non_zero_row);

    if (first_non_zero_row > 0) {
      //////
      // TURN-ON PERMUTATION
      //////
      if (print) printf("TURN-ON PERMUTATION\n");
      // permute first row with first non-zero row
      mpfr_rk2f_switch_rows(ar, 0, first_non_zero_row, N);
      mpfr_rk2f_switch_rows(ai, 0, first_non_zero_row, N);
      
      // permute first col with first non-zero col
      mpfr_rk2f_switch_cols(ar, 0, first_non_zero_row, N, N);
      mpfr_rk2f_switch_cols(ai, 0, first_non_zero_row, N, N);

      // call itself
      mpfr_set(mpfr_tol, mpfr_tol_orig, MPFR_RNDN);
      // cout << "JORDAN WRP: call itself: mpfr_tol = "; mpfr_out_str(stdout, 10, 0, mpfr_tol, MPFR_RNDN); cout << endl;
      jordan_wrp(
        jor, t_mat, ein_pos, num_ein,
        ar, ai, N, nm, debug
      );

      if (print)
      {
      printf("permuted tmat:\n");
      print_rk2_mpc(t_mat, N, N);
      }

      // permute rows in the transformation matrix
      mpc_rk2_switch_rows(t_mat, 0, first_non_zero_row, N);

      // unpermute
      mpfr_rk2f_switch_cols(ar, 0, first_non_zero_row, N, N);
      mpfr_rk2f_switch_cols(ai, 0, first_non_zero_row, N, N);
      mpfr_rk2f_switch_rows(ar, 0, first_non_zero_row, N);
      mpfr_rk2f_switch_rows(ai, 0, first_non_zero_row, N);

      return;
    }

    if (triv_jor) {
      // find the last non-zero row
      last_non_zero_row = N-1;
      for (int i=N-1; i>=0; i--) {
        zero_row = 1;
        for (int j=0; j<N; j++) {
          if (!is_zero[i+j*N]) {
            zero_row = 0;
            break;
          }
        }
        if (zero_row) {
          last_non_zero_row--;
          continue;
        } else {
          break;
        }
      }
      if (print) printf("last non-zero row = %d\n", last_non_zero_row);

      // get the only non-zero eigenvalue
      if (print) printf("get the only non-zero eigenvalue\n");
      ac += first_non_zero*N;
      mpc_set(eig, ac[first_non_zero], MPFR_RNDN);
      for (int i=0; i<N-1-first_non_zero; i++) {
        mpc_fma(eig, prop[first_non_zero+i], ac[first_non_zero+i+1], eig, MPFR_RNDN);
      }
      if (print) {printf("eig = "); print_mpc(&eig); printf("\n");}
      ac -= first_non_zero*N;

      // check whether the only non-zero eigenvalues accidentaly
      // happens to be zero as well
      if (mpc_lessthan_tol(eig)) {
        triv_jor = 2;
      }
    }
  }

  mpfr_set(mpfr_tol, mpfr_tol_orig, MPFR_RNDN);

  // triv_jor = 0;
  if (triv_jor == 1) {
    if (print) printf("TRIVIAL JORDAN 1\n");
    jordan_triv_1(
      jor, t_mat, ein_pos, num_ein,
      ac, N, prop, &eig, first_non_zero,
      print
    );
    free(prop);
  } else if (triv_jor == 2) {
    if (print) printf("TRIVIAL JORDAN 2\n");
    jordan_triv_2(
      jor, t_mat, ein_pos, num_ein,
      ac, N, prop, first_non_zero, last_non_zero_row,
      is_zero+first_non_zero*N,
      print
    );
  } else {
    if (print) printf("STANDARD JORDAN\n");
    jordan(
      jor, t_mat, ein_pos, num_ein,
      ar, ai, N, nm, debug
    );
  }

  free(ac);
}

