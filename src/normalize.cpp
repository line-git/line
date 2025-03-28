#include <iostream>
#include "mpc.h"

using namespace std;

#include "global_vars.h"
#include "malloc_defs.h"
#include "utils.h"
// #include "tensor_utils.h"
// #include "conversions.h"
// #include "setup.h"
#include "tensor_utils.h"
// #include "system_analyzer.h"
#include "jordan_wrp.h"
extern "C" {
  #include "algebra.h"
  #include "jordan.h"
  #include "rel_err_mpc.h"
}


// needed for Jordan decomposition
static int c__0 = 0;


void RLeeAlg1(
  // OUTPUT
  int *k0, int *dimS, int *&S, mpc_t **Delta,
  // INPUT
  mpc_t **L0, int num_ein, int r, int *chain_lengths
) {
  int print = 0;
	// TOLERANCE PARAMETERS
	// mpfr_t tol, told, tmp;
  mpfr_t told, tmp;
  mpfr_init2(told, wp2);
  // mpfr_init2(tol, wp2);
  // mpfr_set_d(tol, 10, MPFR_RNDN);
  // mpfr_pow_si(tol, tol, 8, MPFR_RNDN);
  mpfr_init2(tmp, wp2);
  // mpfr_set_d(tmp, 2, MPFR_RNDN);
	// mpfr_pow_si(tmp, tmp, - ((int) wp2), MPFR_RNDN);
  // mpfr_mul(tol, tol, tmp, MPFR_RNDN);
  // cout << "tolerance = ";
  // mpfr_out_str(stdout, 10, 0, mpfr_tol, MPFR_RNDN);
  // cout << endl;
  // getchar();

  int wp_bin = - mpfr_log2_int(mpfr_tol);
  mpc_t tmp_fma; mpc_init3(tmp_fma, wp2, wp2);	

  // {
  // int n = 3, m = 2;
  // int actual_n = 2;
  // static int c__0 = 0;
  // mpc_t *tu = new mpc_t[n*n];
  // for (int i=0; i<n*n; i++) {
  //   mpc_init3(tu[i], wp2, wp2);
  // }
  // mpc_t *tv= new mpc_t[m*m];
  // for (int i=0; i<m*m; i++) {
  //   mpc_init3(tv[i], wp2, wp2);
  // }
  // mpfr_t *ts = new mpfr_t[min(n, m)];
  // for (int i=0; i<min(n, m); i++) {
  //   mpfr_init2(ts[i], wp2);
  // }
  // mpc_t *tmat = new mpc_t[n*m];
  // for (int i=0; i<n*m; i++) {
  //   mpc_init3(tmat[i], wp2, wp2);
  // }
  // mpc_set_d_d(tmat[0], 1, 0, MPFR_RNDN);
  // mpc_set_d_d(tmat[1], 2, 0, MPFR_RNDN);
  // mpc_set_d_d(tmat[2], 5, 0, MPFR_RNDN);
  // mpc_set_d_d(tmat[3], 3, 0, MPFR_RNDN);
  // mpc_set_d_d(tmat[4], 4, 0, MPFR_RNDN);
  // mpc_set_d_d(tmat[5], 5, 0, MPFR_RNDN);
  
  // cout << "tmat before: " << endl;
  // for (int i=0; i<n*m; i++) {
  //   print_mpc(&tmat[i]);
  //   cout << ", ";
  // }
  // getchar();
  // csvd1_(tmat, &n, &m, &actual_n, &m, &c__0, &actual_n, &m, ts, tu, tv, min(actual_n, m));
  // cout << "tmat after: " << endl;
  // for (int i=0; i<n*m; i++) {
  //   print_mpc(&tmat[i]);
  //   cout << ", ";
  // }
  // cout << endl << endl;
  // cout << "s: " << endl;
  // for (int i=0; i<min(n, m); i++) {
  //   mpfr_out_str(stdout, 10, 0, ts[i], MPFR_RNDN);
  //   cout << ", ";
  // }
  // cout << endl;
  // cout << "u: " << endl;
  // for (int i=0; i<n*n; i++) {
  //   mpc_out_str(stdout, 10, 0, tu[i], MPFR_RNDN);
  //   cout << ", ";
  // }
  // cout << endl;
  // cout << "v: " << endl;
  // for (int i=0; i<m*m; i++) {
  //   mpc_out_str(stdout, 10, 0, tv[i], MPFR_RNDN);
  //   cout << ", ";
  // }
  // cout << endl << endl;
  // getchar();
  // }

  //////
  // PREPARATION
  //////
  // cout << "I am in RLeeAlg1" << endl;

  //// initialize Delta to zero
  for (int i=0; i<num_ein; i++) {
    for (int j=0; j<num_ein; j++) {
      mpc_set_d_d(Delta[i][j], 0, 0, MPFR_RNDN);
    }
  }
  // static int c__0 = 0;
  int row_dim;
  *dimS = 0;
  S = new int[num_ein];
  // int *row_idx, *col_idx;
  // row_idx = new int[num_ein];
  // col_idx = new int[num_ein];
  int count_idx = 0;

  mpc_t *L0t = new mpc_t[num_ein*num_ein];
  init_rk1_mpc(L0t, num_ein*num_ein);
  // cout << "allocated L0t" << endl;

  // // LINE 5
  // for (int k=0; k++; k<num_ein) {
  //   if (k == S[count_S]) {
  //     count_S++;
  //     continue;
  //   }
  //   row_idx[count_idx] = k;
  //   count_idx++;
  // }

  //////
  // OUTER LOOP
  //////
  // LINE 4
  bool isInS = false, next_col;
  int col_rank, new_col_rank, num_ein2 = num_ein * num_ein, col_idx;
  int row_dim2, col_dim, col_dim2, min_dim, null_idx;
  int *rows_to_change = new int[num_ein];
  mpc_t *u = new mpc_t[num_ein2];
  for (int i=0; i<num_ein2; i++) {
    mpc_init3(u[i], wp2, wp2);
  }
  mpc_t *v= new mpc_t[num_ein2];
  for (int i=0; i<num_ein2; i++) {
    mpc_init3(v[i], wp2, wp2);
  }
  mpfr_t *s = new mpfr_t[num_ein];
  for (int i=0; i<num_ein; i++) {
    mpfr_init2(s[i], wp2);
  }
  while(1) {
    // if only one index is left out of S, add it to S and break
    if (*dimS == num_ein - 1) {
      // find complement to S
      for (int i=0; i<num_ein; i++) {
        for (int ip=0; ip<*dimS; ip++) {
          if (S[ip] == i) {
            isInS = true;
            break;
          }
        }
        if (!isInS) {
          S[num_ein-1] = i;
          break;
        }
      }
    }
    
    row_dim = num_ein - *dimS;
    if (print || dbg) cout << "row_dim = " << row_dim << endl;

    //////
    // INNER LOOP
    //////
    col_rank = 0;
    col_idx = 0;
    next_col = false;
    col_dim = 0;
    for (int j=0; j<num_ein; j++) {
      if (print || dbg) {cout << endl; cout << "j = " << j << endl;}
      // getchar();

      isInS = false;
      for (int l=0; l<*dimS; l++) {
        if (j == S[l]) {
          isInS = true;
          break;
        }
      }
      if (isInS) {
        if (print || dbg) cout << "is in S" << endl;
        continue;
      }
      rows_to_change[col_dim] = j;
      col_dim++;
      if (print || dbg) cout << "col_dim = " << col_dim << endl;

      // LINE 5
      //// remove lines from L0 while preparing the input for the SVD
      for (int jt, jp=0; jp<col_dim; jp++) {
        jt = rows_to_change[jp];
        count_idx = 0;
        for (int k=0; k<num_ein; k++) {
          isInS = false;
          for (int l=0; l<*dimS; l++) {
            // cout << "l = " << l << endl;
            if (k == S[l]) {
              isInS = true;
              break;
            }
          }
          if (isInS) {
            continue;
          }
          // cout << "count_idx, jp = " << count_idx << ", " << jp << endl;
          // mpc_init3(L0t[count_idx + jp*row_dim], wp2, wp2);
          mpc_set(L0t[count_idx + jp*row_dim], L0[k][jt], MPFR_RNDN);
          count_idx++;
        }
      }

      //// CSVD
      row_dim2 = row_dim * row_dim;
      // col_dim = j+1;
      col_dim2 = col_dim * col_dim;
      min_dim = min(row_dim, col_dim);

      if (print || dbg) {
      cout << "entering csvd with:" << endl;
      cout << "num_ein = " << num_ein << endl;
      cout << "row_dim = " << row_dim << endl;
      cout << "col_dim = " << col_dim << endl;
      cout << "L0t before: " << endl;
      for (int i=0; i<row_dim; i++) {
        for (int jp=0; jp<col_dim; jp++) {
          print_mpc(&L0t[i + row_dim*jp]);
          cout << ", ";
        }
        cout << endl; 
      }
      cout << endl << endl;
      }
      // // getchar();
      csvd1_(L0t, &row_dim, &col_dim, &row_dim, &col_dim, &c__0, &row_dim, &col_dim, s, u, v, min_dim, 1);
      // cout << "CSVD done" << endl;
      
      if (print || dbg) {
      cout << "L0t after: " << endl;
      for (int i=0; i<row_dim; i++) {
        for (int jp=0; jp<col_dim; jp++) {
          print_mpc(&L0t[i + row_dim*jp]);
          cout << ", ";
        }
        cout << endl;
      }
      cout << endl << endl;
      
      cout << "s: " << endl;
      for (int i=0; i<min_dim; i++) {
        mpfr_out_str(stdout, 10, 0, s[i], MPFR_RNDN);
        cout << ", ";
      }
      cout << endl << endl;
      
      cout << "u: " << endl;
      for (int jp=0; jp<row_dim; jp++) {
        for (int i=0; i<row_dim; i++) {
          mpc_out_str(stdout, 10, 0, u[i + jp*row_dim], MPFR_RNDN);
          cout << ", ";
        }
        cout << endl;
      }
      cout << endl << endl;
    
      cout << "v: " << endl;
      for (int jp=0; jp<col_dim; jp++) {
        for (int i=0; i<col_dim; i++) {
          mpc_out_str(stdout, 10, 0, v[i + jp*col_dim], MPFR_RNDN);
          cout << ", ";
        }
        cout << endl;
      }
      cout << endl;
      }
      // // getchar();

      //// determine rank
      new_col_rank = 0;
      mpfr_mul_ui(told, mpfr_tol, row_dim*col_dim, MPFR_RNDN);
      for (int i=0; i<min_dim; i++) {
        mpfr_abs(tmp, s[i], MPFR_RNDN);
        if (mpfr_greater_p(tmp, told)) {
          new_col_rank++;
        }
      }
      // cout << "new_col_rank = " << new_col_rank << endl;
      // getchar();

      // LINE 6
      col_idx = j;
      if (new_col_rank == col_rank) {
        // check whether column index is already in S
        isInS = false;
        for (int l=0; l<*dimS; l++) {
          if (S[l] == col_idx) {
            isInS = true;
            break;
          }
        }
        if (isInS) {
          continue;
        }
        
        // if not, then the search is over
        break;
      } else {
        col_rank = new_col_rank;
        continue;
      }
      // getchar();
    }
    // cout << "col_idx = " << col_idx << endl;
    // getchar();
    
    // find coefficients from null space
    //// select vector of null space that has a non-zero last component
    for (int j=0; j<col_dim - col_rank; j++) {
      if (v[col_dim-1 + j*col_dim] != 0) {
        null_idx = j;
        // normalize coefficients so that last component equals one
        for (int i=0; i<col_dim-1; i++) {
          // mpc_div(
          exp_rel_err_mpc_div(
            v[i + j*col_dim], v[i + j*col_dim], v[col_dim-1 + j*col_dim], MPFR_RNDN
            , wp_bin
          );
        }
      }
    }
    // cout << "null_idx = " << null_idx << endl;

    // LINE 9
    //// set to zero the selected column of L0
    for (int i=0; i<num_ein; i++) {
      mpc_set_d_d(L0[i][col_idx], 0, 0, MPFR_RNDN);
    }
    //// change rows of L0 whose Jordan chains have the same
    //// length of the one associated with the selected column
    for (int it, i=0; i<col_dim-1; i++ ) {
      it = rows_to_change[i];
      if (chain_lengths[it] == chain_lengths[col_idx]) {
        // subtract to the row c_i times the row of L0 with row index equal
        // to the index of the selected column
        for (int j=0; j<num_ein; j++ ) {
          mpc_neg(v[i + null_idx*col_dim], v[i + null_idx*col_dim], MPFR_RNDN);
          // mpc_fma(L0[it][j], v[i + null_idx*col_dim], L0[col_idx][j], L0[it][j], MPFR_RNDN);
          exp_rel_err_mpc_mul(tmp_fma, v[i + null_idx*col_dim], L0[col_idx][j], MPFR_RNDN, wp_bin);
          exp_rel_err_mpc_add(L0[it][j], L0[it][j], tmp_fma, MPFR_RNDN, wp_bin);
          mpc_neg(v[i + null_idx*col_dim], v[i + null_idx*col_dim], MPFR_RNDN);
        }
      }
    }
    
    if (print || dbg) {
    cout << "new L0: " << endl;
    for (int i=0; i<num_ein; i++) {
      for (int jp=0; jp<num_ein; jp++) {
        print_mpc(&L0[i][jp]);
        cout << ", ";
      }
      cout << endl;
    }
    cout << endl << endl;
    }

    // LINE 10
    // update Delta matrix
    //// set selected column to itself minus a linear combination
    //// of previuous columns with coefficients given by the null space
    for (int i=0; i<num_ein; i++) {
      for (int jt, j=0; j<col_dim-1; j++) {
        jt = rows_to_change[j];
        // mpc_fma(Delta[i][col_idx], v[j + null_idx*col_dim], Delta[i][jt], Delta[i][col_idx], MPFR_RNDN);
        exp_rel_err_mpc_mul(tmp_fma, v[j + null_idx*col_dim], Delta[i][jt], MPFR_RNDN, wp_bin);
        exp_rel_err_mpc_add(Delta[i][col_idx], Delta[i][col_idx], tmp_fma, MPFR_RNDN, wp_bin);
      }
    }
    //// finally, add coefficients of the null space to selected column
    // cout << "check on v" << endl;
    for (int it, i=0; i<col_dim-1; i++) {
      it = rows_to_change[i];
      // print_mpc(&v[i + null_idx*col_dim]);
      // cout << endl;
      exp_rel_err_mpc_add(Delta[it][col_idx], Delta[it][col_idx], v[i + null_idx*col_dim], MPFR_RNDN, wp_bin);
    }

    // LINE 11
    // update S
    // cout << "adding " << col_idx << " to S" << endl;
    S[*dimS] = col_idx;
    (*dimS)++;
    // cout << "dimS = " << *dimS << endl;
    if (*dimS > num_ein) {
      fprintf(stdout, "error in RLee algorithm: dimS > # eigenvalues\n");
      fprintf(stderr, "error in RLee algorithm: dimS > # eigenvalues\n");
      exit(1);
    }

    // cout << "new Delta: " << endl;
    // for (int i=0; i<num_ein; i++) {
    //   for (int jp=0; jp<num_ein; jp++) {
    //     print_mpc(&Delta[i][jp]);
    //     cout << ", ";
    //   }
    //   cout << endl;
    // }
    // cout << endl << endl;
    // getchar();
    
    // cout << "internal S: " << endl;
    // for (int a=0; a<*dimS; a++) {
    //   cout << S[a] << ", ";
    // }
    // cout << endl << endl;
    // getchar();

    if (col_idx <= r-1) {
      break;
    }
  }

  // LINE 13
  // remove the last selected column index from S
  *k0 = S[*dimS-1];
  (*dimS)--;

  // cout << "final S: " << endl;
  // for (int a=0; a<*dimS; a++) {
  //   cout << S[a] << ", ";
  // }
  // cout << endl << endl;
  
  // cout << "k0 = " << *k0 << endl << endl;

  // cout << "final Delta: " << endl;
  // for (int i=0; i<num_ein; i++) {
  //   for (int jp=0; jp<num_ein; jp++) {
  //     print_mpc(&Delta[i][jp]);
  //     cout << ", ";
  //   }
  //   cout << endl;
  // }
  // cout << endl << endl;
  // getchar();

  // delete[] u;
  // delete[] v;
  // delete[] s;

  // FREE
  mpc_clear(tmp_fma);
}


void find_projector(
  // OUTPUT
  mpc_t **proj,
  // INPUT
  mpfr_t *A0r, mpfr_t *A0i, mpc_t **A1, int dim
) {
  int print = 0;
  int jordan_debug = 0;
  if (print == 1) {
    jordan_debug = -1;
  }
  #define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
  int temp;

  int wp_bin = - mpfr_log2_int(mpfr_tol);

  // JORDAN DECOMPOSITION OF A0
  mpc_t **jor, **t_mat;
  malloc_rk2_tens(jor, dim, dim);
  init_rk2_mpc(jor, dim, dim);
  malloc_rk2_tens(t_mat, dim, dim);
  init_rk2_mpc(t_mat, dim, dim);
  int num_ein, *ein_pos;
  jordan_wrp(jor, t_mat, &ein_pos, &num_ein, A0r, A0i, dim, dim, jordan_debug);
  mpc_rk2_prune_abs_tol(jor, dim, dim);
  mpc_rk2_prune_abs_tol(t_mat, dim, dim);

  if (print || dbg) {
  cout << "1st jordan tmat:" << endl;
  print_rk2_mpc(t_mat, dim, dim);
  cout << "1st jordan mat:" << endl;
  print_rk2_mpc(jor, dim, dim);
  }

  // cout << "num_ein = " << num_ein << endl;
  // cout << "ein_pos: " << endl;
  // for (int i=0; i<num_ein; i++) {
  //   cout << ein_pos[i] << ", ";
  // }
  // cout << endl << endl;

  // order Jordan chains from longest to shortest
  int *chain_lengths = new int[num_ein];
  for (int i=0; i<num_ein-1; i++) {
    chain_lengths[i] = ein_pos[i+1] - ein_pos[i];
  }
  chain_lengths[num_ein-1] = dim - ein_pos[num_ein-1];

  // cout << "chain_lengths: " << endl;
  // for (int i=0; i<num_ein; i++) {
  //   cout << chain_lengths[i] << ", ";
  // }
  // cout << endl << endl;

  int *perm_ein = new int[num_ein];
  for (int i=0; i<num_ein; i++) {
    perm_ein[i] = i;
  }
  for (int i=0; i<num_ein; i++) {
    for (int j=i; j<num_ein; j++) {
      if (chain_lengths[j] > chain_lengths[i]) {
        SWAP(chain_lengths[i], chain_lengths[j]);
        SWAP(perm_ein[i], perm_ein[j]);
      }
    }
  }

  // cout << "perm_ein: " << endl;
  // for (int i=0; i<num_ein; i++) {
  //   cout << perm_ein[i] << ", ";
  // }
  // cout << endl << endl;
  // getchar();

  //// invert transformation matrix
  mpc_t **t_mat_inv;
  mpc_t **copy_t_mat;
  malloc_rk2_tens(t_mat_inv, dim, dim);
  init_rk2_mpc(t_mat_inv, dim, dim);
  malloc_rk2_tens(copy_t_mat, dim, dim);
  init_rk2_mpc(copy_t_mat, dim, dim);
  copy_rk2_mpc(copy_t_mat, t_mat, dim, dim);
  mp_inverse(t_mat_inv, copy_t_mat, dim);
  mpc_rk2_prune_abs_tol(t_mat_inv, dim, dim);
  if (print || dbg) {
    cout << "1st jordan inv tmat:" << endl;
    print_rk2_mpc(t_mat_inv, dim, dim);
  }

  // BUILD L0 AND L1
  mpc_t **L0, **L1;
  malloc_rk2_tens(L0, num_ein, num_ein);
  init_rk2_mpc(L0, num_ein, num_ein);
  malloc_rk2_tens(L1, num_ein, num_ein);
  init_rk2_mpc(L1, num_ein, num_ein);
  mpc_t tmp;
  mpc_init3(tmp, wp2, wp2);
  mpc_t tmpflop;
  mpc_init3(tmpflop, wp2, wp2);

  for (int row, col, i=0; i<num_ein; i++) {
    // the row index picks the position of the head of the Jordan chain;
    if (perm_ein[i] == num_ein -1) {
      row = dim - 1;
    } else {
      row = ein_pos[perm_ein[i]] + chain_lengths[i] - 1;
    }
    
    for (int j=0; j<num_ein; j++) {
      col = ein_pos[perm_ein[j]];
      mpc_set_d_d(L0[i][j], 0, 0, MPFR_RNDN);
      mpc_set_d_d(L1[i][j], 0, 0, MPFR_RNDN);
      for (int h=0; h<dim; h++) {
        // mpc_fma(L1[i][j], t_mat_inv[row][h], t_mat[h][col], L1[i][j], MPFR_RNDN);
        // cout << "t_mat_inv[row][h] = "; print_mpc(&t_mat_inv[row][h]); cout << endl;
        // cout << "t_mat[h][col] = "; print_mpc(&t_mat[h][col]); cout << endl;
        exp_rel_err_mpc_mul(tmpflop, t_mat_inv[row][h], t_mat[h][col], MPFR_RNDN, wp_bin);
        // cout << "tmpflop = "; print_mpc(&tmpflop); cout << endl;
        exp_rel_err_mpc_add(L1[i][j], L1[i][j], tmpflop, MPFR_RNDN, wp_bin);
        for (int k=0; k<dim; k++) {
          exp_rel_err_mpc_mul(tmp, A1[h][k], t_mat[k][col], MPFR_RNDN, wp_bin);
          // cout << "tmp = "; print_mpc(&tmp); cout << endl;
          // mpc_fma(L0[i][j], t_mat_inv[row][h], tmp, L0[i][j], MPFR_RNDN);
          exp_rel_err_mpc_mul(tmpflop, t_mat_inv[row][h], tmp, MPFR_RNDN, wp_bin);
          // cout << "tmpflop = "; print_mpc(&tmpflop); cout << endl;
          exp_rel_err_mpc_add(L0[i][j], L0[i][j], tmpflop, MPFR_RNDN, wp_bin);
        }
      }
    }
  }

  if (print || dbg) {
  cout << "L0: " << endl; print_rk2_mpc(L0, num_ein, num_ein);
  cout << "L1: " << endl; print_rk2_mpc(L1, num_ein, num_ein);
  }

  // count non-zero elements on the L1 diagonal
  int r = 0;
  for (int i=0; i<num_ein; i++) {
    // if (mpfr_zero_p(mpc_realref(L1[i][i]))) {
    if (mpc_lessthan_tol(L1[i][i])) {
      r++;
    }
  }
  if (print || dbg) cout << "r = " << r << endl;

  // Apply Roman Lee Algorithm 1 to build matrix Delta
  int k0, *S, dimS;
  mpc_t **Delta;
  malloc_rk2_tens(Delta, num_ein, num_ein);
  init_rk2_mpc(Delta, num_ein, num_ein);
  RLeeAlg1(&k0, &dimS, S, Delta, L0, num_ein, r, chain_lengths);
  if (print || dbg) cout << "dimS = " << dimS << endl;
  if (print || dbg) cout << "k0 = " << k0 << endl;
  if (print || dbg) {cout << "Delta: " << endl; print_rk2_mpc(Delta, num_ein, num_ein);}

  // Transform matrix U according to U -> U*(1 + E) and U^-1 accordingly
  //// (for U^-1 we have the same operations but with the rows instead of
  //// the columns, and with a negative flop instead of a positive one)
  //
  //// first transform every column by superposing previous
  //// columns with proper indices;
  ////
  //// start with the column indices in S, from the greatest to the 
  //// smallest: in this way we can dinamically change the columns of U
  //// since every columns is superposed only with the previous ones;
  
  //// order elements of S from greatest to smallest
  for (int i=0; i<dimS; i++) {
    for (int j=i; j<dimS; j++) {
      if (S[j] > S[i]) {
        SWAP(S[i], S[j]);
      }
    }
  }

  //// linearly combine elements of previous columns
  //// with the coefficients stored in Delta
  for (int k, kp=0; kp<dimS; kp++) {
    k = ein_pos[perm_ein[S[kp]]];
    // build the linear combination component by component
    for (int i=0; i<dim; i++) {
      for (int l, lp=0; lp<S[kp]; lp++) {
        l = ein_pos[perm_ein[lp]];
        // cout << "k, i, l = " << k << ", " <<  i << ", " << l << endl;
        // mpc_fma(t_mat[i][k], Delta[lp][S[kp]], t_mat[i][l], t_mat[i][k], MPFR_RNDN);
        exp_rel_err_mpc_mul(tmpflop, Delta[lp][S[kp]], t_mat[i][l], MPFR_RNDN, wp_bin);
        exp_rel_err_mpc_add(t_mat[i][k], t_mat[i][k], tmpflop, MPFR_RNDN, wp_bin);
        mpc_neg(t_mat_inv[l][i], t_mat_inv[l][i], MPFR_RNDN);
        // mpc_fma(t_mat_inv[l][i], Delta[lp][S[kp]], t_mat_inv[k][i], t_mat_inv[l][i], MPFR_RNDN);
        exp_rel_err_mpc_mul(tmpflop, Delta[lp][S[kp]], t_mat_inv[k][i], MPFR_RNDN, wp_bin);
        exp_rel_err_mpc_add(t_mat_inv[l][i], t_mat_inv[l][i], tmpflop, MPFR_RNDN, wp_bin);
        mpc_neg(t_mat_inv[l][i], t_mat_inv[l][i], MPFR_RNDN);
      }
    }
  }
  //// now we just need to process k0 (if column index is not
  //// in S U {k0} then the corresponding column in Delta is zero)
  ////
  //// for k0 we need to combine every indices corresponding to
  //// every element of the Jordan chain, sice the chain of k0
  //// has more than one element
  for (int k, kp=0; kp<chain_lengths[k0]; kp++) {
    k = ein_pos[perm_ein[k0]] + kp;
    // build the linear combination component by component
    for (int i=0; i<dim; i++) {
      for (int l, lp=0; lp<k0; lp++) {
        l = ein_pos[perm_ein[lp]] + kp;
        // mpc_fma(t_mat[i][k], Delta[lp][k0], t_mat[i][l], t_mat[i][k], MPFR_RNDN);
        exp_rel_err_mpc_mul(tmpflop, Delta[lp][k0], t_mat[i][l], MPFR_RNDN, wp_bin);
        exp_rel_err_mpc_add(t_mat[i][k], t_mat[i][k], tmpflop, MPFR_RNDN, wp_bin);
        mpc_neg(t_mat_inv[l][i], t_mat_inv[l][i], MPFR_RNDN);
        // mpc_fma(t_mat_inv[l][i], Delta[k0][lp], t_mat_inv[k][i], t_mat_inv[l][i], MPFR_RNDN);
        exp_rel_err_mpc_mul(tmpflop, Delta[k0][lp], t_mat_inv[k][i], MPFR_RNDN, wp_bin);
        exp_rel_err_mpc_add(t_mat_inv[l][i], t_mat_inv[l][i], tmpflop, MPFR_RNDN, wp_bin);
        mpc_neg(t_mat_inv[l][i], t_mat_inv[l][i], MPFR_RNDN);
      }
    }
  }


  // // test if inverse is really the inverse
  // mpc_t **test;
  // malloc_rk2_tens(test, dim, dim);
  // init_rk2_mpc(test, dim, dim);
  // cout << "check the inverse:" << endl;
  // getchar();
  // for (int i=0; i<dim; i++) {
  //     for (int j=0; j<dim; j++) {
  //       mpc_set_d_d(test[i][j], 0, 0, MPFR_RNDN);
  //         for (int k=0; k<dim; k++) {
  //           mpc_fma(test[i][j], t_mat[i][k], t_mat_inv[k][j], test[i][j], MPFR_RNDN);
  //         }
  //         print_mpc(&test[i][j]);
  //         cout << endl;
  //     }
  //     cout << endl;
  // }
  // getchar();


  // Construct the projector
  int k;
  for (int i=0; i<dim; i++) {
    for (int j=0; j<dim; j++) {
      // project over k0
      k = ein_pos[perm_ein[k0]];
      // if (i==0 && j==0)
      // cout << "k0, perm_ein, k = " << k0 << ", " << perm_ein[k0] << ", " << k << endl;
      exp_rel_err_mpc_mul(
        proj[i][j], t_mat[i][k], t_mat_inv[k][j], MPFR_RNDN
        , wp_bin
      );
      // project over elements of S
      for (int kp=0; kp<dimS; kp++) {
        k = ein_pos[perm_ein[S[kp]]];
        // if (i==0 && j==0)
        // cout << "S[kp], perm_ein, k = " << S[kp] << ", " << perm_ein[S[kp]] << ", " << k << endl;
        // mpc_fma(proj[i][j], t_mat[i][k], t_mat_inv[k][j], proj[i][j], MPFR_RNDN);
        exp_rel_err_mpc_mul(
          tmp, t_mat[i][k], t_mat_inv[k][j], MPFR_RNDN
          , wp_bin
        );
        exp_rel_err_mpc_add(proj[i][j], proj[i][j], tmp, MPFR_RNDN, wp_bin);
      }
    }
  }

  // FREE
  mpc_clear(tmp);
  mpc_clear(tmpflop);

  #undef SWAP
}


void find_eigen_duplicates(
  // IN-OUT
  int* eig_grid,
  //OUTPUT
  mpc_t *eigenvalues, int dim
) {
  mpfr_t abs_diff;
  mpfr_init2(abs_diff, wp2);
  for (int i=0; i<dim; i++) {
    if (eig_grid[i] == 0) {
      continue;
    }
    for (int j=i-1; j>=0; j--) {
      if (eig_grid[j] == 0) {
        continue;
      }
      // mpfr_sub(abs_diff, mpc_realref(eigenvalues[i]), mpc_realref(eigenvalues[j]), MPFR_RNDN);
      // mpfr_abs(abs_diff, abs_diff, MPFR_RNDN);
      // if (mpfr_less_p(abs_diff, mpfr_tol)) {
      if (mpc_equal_within_tol(eigenvalues[i], eigenvalues[j])) {
        // +1 is added by convention because the negative of zero stays zero
        eig_grid[i] = -(j+1);
        break;
      }
    }
  }
}


void decompose_eigenvalues(
  // OUTPUT
  int *eig_int_diff,
  // INPUT
  int dim, mpc_t *eigenvalues, int *eig_grid,
  int *eq_class, int *class_min
) {
  mpfr_t diff;
  mpfr_init2(diff, wp2);
  for (int i=0; i<dim; i++) {
    // cout << "i = " << i << endl;
    // if this is the minimum, no need to manually decompose
    if (i == class_min[eq_class[i]]) {
      eig_int_diff[i] = 0;
      continue;
    }
    if (eig_grid[i] == 0) {
      continue;
      }
    if (eig_grid[i] < 0) {
      // assign integer difference to duplicates
      // cout << "assign integer difference to duplicates" << endl;
      // cout << "eig_grid[" << i << "] = " << eig_grid[i] << endl;
      // cout << "eig_int_diff = " << eig_int_diff[-eig_grid[i]-1] << endl;
      eig_int_diff[i] = eig_int_diff[-eig_grid[i]-1];
      continue;
    }
    
    // find the difference w.r.t. the class minimum
    mpfr_sub(diff, mpc_realref(eigenvalues[i]), mpc_realref(eigenvalues[class_min[eq_class[i]]]), MPFR_RNDN); 
    eig_int_diff[i] = (int) mpfr_get_si(diff, MPFR_RNDN);
  }

}


void group_eigenvalues(
  // OUTPUT
  int *&out_class_min, int *num_classes,
  // INPUT
  int dim, mpc_t *eigenvalues, int *eig_grid,
  int *eq_class, int *eig_int_diff
) {
  int print = 0;
  // find eigenvalues with the same real part
  // cout << "find duplicates" << endl;
  find_eigen_duplicates(eig_grid, eigenvalues, dim);
  // for (int i=0; i<dim; i++) {
  //   cout << "eig_grid[" << i << "] = " << eig_grid[i] << endl;
  // }

  // cout << "assign classes" << endl;
  int *class_min = new int[dim];
  bool found = false;
  *num_classes = 0;
  for (int i=0; i<dim; i++) {
    // cout << "i = " << i << endl;
    if (eig_grid[i] == 0) {
      // assign class to Jordan vectors that are not pure eigenstates
      //// find the chain bottom (first appearence of a positive or of a negative value
      //// in eig_grid going backward) and copy its class label
      //// (for sure the bottom of the chain has already been processed)
      for (int j=i-1; j>=0; j--) {
        if (eig_grid[j] >= 1 || eig_grid[j] < 0) {
          eq_class[i] = eq_class[j];
          break;
        }
      }
      continue;
    }
    if (eig_grid[i] < 0) {
      // assign class to duplicates
      //// rember that +1 was added by convention, so it hase to be subtracted back
      eq_class[i] = eq_class[-eig_grid[i]-1];
      // update minimum of the class
      if (mpfr_less_p(mpc_realref(eigenvalues[i]), mpc_realref(eigenvalues[class_min[eq_class[i]]]))) {
        class_min[eq_class[i]] = i;
      }
      continue;
    }
    // check if eigenvalue belongs to a class already appeared before
    found = false;
    for (int j=i-1; j>=0; j--) {
      if (eig_grid[j] == 0) {
        continue;
      }
      if (mpfr_is_diff_int(mpc_realref(eigenvalues[i]), mpc_realref(eigenvalues[j])) &&\
        mpfr_equal_within_tol(mpc_imagref(eigenvalues[i]), mpc_imagref(eigenvalues[j]))
      ) {
        found = true;
        eq_class[i] = eq_class[j];
        // update minimum of the class
        //// the mininum of the class labeled by eq_class[j] is stored 
        //// into the position given by class_min[...]
        if (mpfr_less_p(mpc_realref(eigenvalues[i]), mpc_realref(eigenvalues[class_min[eq_class[j]]]))) {
          class_min[eq_class[j]] = i;
        }
        break;
      }
    }

    // if this is the first appearance of the eigenvalue,
    // then a new class is found
    if (!found) {
      // assign the new label
      eq_class[i] = *num_classes;
      // assign position of the minimum
      class_min[*num_classes] = i;
      // increas number of classes
      (*num_classes)++;
    }
  }

  // alloc and prepare the output containing the position of the minima
  out_class_min = new int[*num_classes];
  for (int i=0; i<*num_classes; i++) {
    out_class_min[i] = class_min[i];
  }

  // DECOMPOSE EIGENVALUES
  // cout << "decompose" << endl;
  decompose_eigenvalues(
    eig_int_diff,
    dim, eigenvalues, eig_grid,
    eq_class, class_min
  );

  // print for debug
  if (print || dbg) {
  cout << "INPUT OF EIGENVALUE NORMALIZATION:" << endl;
  cout << "line, grid, class label, min pos, int diff, eigenvalues" << endl;
  for (int i=0; i<dim; i++) {
    cout << i << "\t" << eig_grid[i] << "\t" << eq_class[i] << "\t" << class_min[eq_class[i]] << "\t" << eig_int_diff[i] << "\t";
    print_mpc(&eigenvalues[i]);
    cout << endl;
  }
  cout << endl;
  }

  // // always use lam=0 as eq class = 0
  
  // // find label of lam=0
  // mpfr_t mpfr_zero;
  // mpfr_init2(mpfr_zero, wp2);
  // mpfr_set_ui(mpfr_zero, 0, MPFR_RNDN);
  // int lam0;
  // for (int i=0; i<dim; i++) {
  //   if (mpfr_lessthan_tol(mpc_imagref(eigenvalues[i])) && mpfr_is_diff_int(mpc_realref(eigenvalues[i]), mpfr_zero)){
  //     cout << "label of lam=0 is " << eq_class[i] << endl;
  //     lam0 = eq_class[i];
  //     break;
  //   }
  // }

  // // change eq class
  // for (int i=0; i<dim; i++) {
  //   if (eq_class[i] == lam0) {
  //     eq_class[i] = 0;
  //   } else if (eq_class[i] == 0) {
  //     eq_class[i] = lam0;
  //   }
  // }

  // // change class min
  // int tmp;
  // tmp = class_min[0];
  // class_min[0] = class_min[lam0];
  // class_min[lam0] = tmp;

  // // print for debug
  // cout << "after using eq class = 0 for lam = 0:" << endl;
  // cout << "line, grid, class label, min pos, int diff, eigenvalues" << endl;
  // for (int i=0; i<dim; i++) {
  //   cout << i << "\t" << eig_grid[i] << "\t" << eq_class[i] << "\t" << class_min[eq_class[i]] << "\t" << eig_int_diff[i] << "\t";
  //   print_mpc(&eigenvalues[i]);
  //   cout << endl;
  // }
  // cout << endl;
}

