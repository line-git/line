#include <iostream>
#include "mpc.h"

using namespace std;

#include "global_vars.h"
#include "malloc_defs.h"
#include "utils.h"
#include "tensor_utils.h"
#include "system_analyzer.h"
#include "poly_frac.h"
#include "normalize.h"
#include "jordan_wrp.h"
extern "C" {
  #include "algebra.h"
  // #include "jordan.h"
  #include "cpoly.h"
  #include "rel_err_mpc.h"
  // #include "in_out.h"
}

// needed for Jordan decomposition
static int c__0 = 0;


//////
// POLY_FRAC NORMALIZE
/////
void pf_to_Fuchsian(
  // OUTPUT
  struct poly_frac **bal,
  struct poly_frac **inv_bal,
  // IN-OUT
  struct poly_frac **matb,
  // INPUT
  int b_len,
  mpc_t *roots, int nroots
) {
  /*
  Tranform the input square matrix to Fuchsian form, providing in output
  the direct and the inverse transformation matrices. The transformation
  law is:
  
    M' = invT*M*T - invT*Derivate[T]
    
  OUTPUT:
    - bal, inv_bal: rk2 tensors of polynomial fractions corresponding to
      the balance transformation and its inverse, respectively;
  IN-OUT:
    - matb: rk2 of polynomial fractions to be transfromed, overwritten
      by the transformed matrix;
  INPUT:
    - b_len: dimension of the square matrices;
    - roots: array containing the poles of the global matrix (where the
      1st element is zero);
    - nroots: number of poles of the global matrix;
  */
  int print = 0;
  int b_len2 = b_len*b_len;
  int p_rank;
  int index;
  // leading and next-to-leading orders
  mpfr_t *A0r, *A0i;
  mpc_t **A1;
  A0r = new mpfr_t[b_len2];
  A0i = new mpfr_t[b_len2];
  init_rk1_mpfr(A0r, b_len2);
  init_rk1_mpfr(A0i, b_len2);
  malloc_rk2_tens(A1, b_len, b_len);
  init_rk2_mpc(A1, b_len, b_len);
  // projector
  mpc_t **proj;
  malloc_rk2_tens(proj, b_len, b_len);
  init_rk2_mpc(proj, b_len, b_len);
  // utilities
  mpc_t mpc_one;
  mpc_init3(mpc_one, wp2, wp2);
  mpc_set_ui(mpc_one, 1, MPFR_RNDN);
  struct poly_frac pf_one;
  poly_frac_build(&pf_one);
  poly_frac_set_ui(&pf_one, 1, nroots);
  // contributions to the transformation matrices
  struct poly_frac **proj_mat, **mat_proj, **proj_mat_proj;
  malloc_rk2_tens(proj_mat, b_len, b_len);
  poly_frac_rk2_build(proj_mat, b_len, b_len);
  malloc_rk2_tens(mat_proj, b_len, b_len);
  poly_frac_rk2_build(mat_proj, b_len, b_len);
  malloc_rk2_tens(proj_mat_proj, b_len, b_len);
  poly_frac_rk2_build(proj_mat_proj, b_len, b_len);
  struct poly_frac **der_contr;
  malloc_rk2_tens(der_contr, b_len, b_len);
  poly_frac_rk2_build(der_contr, b_len, b_len);

  int wp_bin = - mpfr_log2_int(mpfr_tol);

  int num_it = 0;
  while (1) {
    // compute Poinc. rank
    // p_rank = mp_Poinc_rank(mults, b_len, b_len);
    poly_frac_rk2_prune_rel_tol(matb, wp_bin, b_len, b_len);
    p_rank = poly_frac_Poinc_rank(matb, b_len, b_len);
    if (print) {
    cout << "p_rank = " << p_rank << endl;
    cout << "pf block:" << endl;
    poly_frac_rk2_print(matb, b_len, b_len);
    }
    // getchar();
    if (p_rank <= 0) {
      if (num_it == 0) {
        // set transformation matrices to identity
        poly_frac_rk2_set_id(bal, b_len, nroots);
        poly_frac_rk2_set_id(inv_bal, b_len, nroots);
      }
      break;
    }

    // EXTRACT LO AND NLO
    poly_frac_rk2_prune_rel_tol(matb, wp_bin, b_len, b_len);
    rel_err_poly_frac_mat_extract_LO(
      A0r, A0i,
      matb, roots, p_rank+1,
      b_len, b_len
      , wp_bin
    );
    if (print || dbg) {
    cout << "Leading Order:" << endl;
    for (int i=0; i<b_len; i++) {
      cout << "i = " << i << endl;
      for (int j=0; j<b_len; j++) {
        cout << "j = " << j << endl;
        index = i + b_len*j;
        mpfr_out_str(stdout, 10, 0, A0r[index], MPFR_RNDN);
        cout << endl;
        mpfr_out_str(stdout, 10, 0, A0i[index], MPFR_RNDN);
        cout << endl;
      }
    }
    cout << endl;
    }
    // poly_frac_mat_extract_NLOc(
    rel_err_poly_frac_mat_extract_NLOc(
      A1,
      A0r, A0i, matb, roots, p_rank+1,
      b_len, b_len
      , wp_bin
    );
    if (print || dbg) {
    cout << "A1:" << endl; print_rk2_mpc(A1, b_len, b_len); cout << endl;
    }

    // cout << "finding projector..." << endl;
    find_projector(proj, A0r, A0i, A1, b_len);
    if (print || dbg) {
    cout << "projector:" << endl;
    print_rk2_mpc(proj, b_len, b_len);
    }

    // apply balance transformation
    //// build (eta-1)*proj*mat
    // cout << endl;
    // cout << "build proj*mat" << endl;
    rel_err_mpc_rk2_mul_poly_frac_rk2(
      proj_mat, proj, matb, b_len,
      roots,
      wp_bin
    );
    //// build mat*proj
    // cout << endl;
    // cout << "build mat*proj" << endl;
    rel_err_poly_frac_rk2_mul_mpc_rk2(
      mat_proj, matb, proj, b_len,
      roots,
      wp_bin
    );
    //// build proj*mat*proj
    // cout << endl;
    // cout << "build proj*mat*proj" << endl;
    rel_err_poly_frac_rk2_mul_mpc_rk2(
      proj_mat_proj, proj_mat, proj, b_len,
      roots,
      wp_bin
    );
    //// multiply f(eta) factors
    // cout << endl;
    // cout << "multiply eta factors" << endl;
    for (int i=0; i<b_len; i++) {
      for (int j=0; j<b_len; j++) {
        // cout << "i, j = " << i << ", " << j << endl;
        // (eta-1) * proj*mat
        poly_frac_mul_root(&proj_mat[i][j], &mpc_one, 1);
        // cout << "done proj*mat" << endl;
        // -(eta-1)/eta * mat*proj
        poly_frac_mul_root(&mat_proj[i][j], &mpc_one, 1);
        poly_frac_mul_sym_pow(&mat_proj[i][j], &mat_proj[i][j], -1);
        poly_frac_neg(&mat_proj[i][j]);
        // cout << "done mat*proj" << endl;
        // -(eta-1)^2/eta * proj*mat*proj
        poly_frac_mul_root(&proj_mat_proj[i][j], &mpc_one, 2);
        poly_frac_mul_sym_pow(&proj_mat_proj[i][j], &proj_mat_proj[i][j], -1);
        poly_frac_neg(&proj_mat_proj[i][j]);
        // cout << "done proj*mat*proj" << endl;
      }
    }
    //// construct derivative contribution
    // cout << endl;
    // cout << "construct derivative contribution" << endl;
    for (int i=0; i<b_len; i++) {
      for (int j=0; j<b_len; j++) {
        // convert projector to poly_frac
        poly_frac_set_mpc(&der_contr[i][j], &proj[i][j], nroots);
        // multiply it times 1/eta
        poly_frac_mul_sym_pow(&der_contr[i][j], &der_contr[i][j], -1);
        // poly_frac_neg(&der_contr[i][j]);
        // poly_frac_print(&der_contr[i][j], nroots);
      }
    }
    //// sum all the contributions
    // cout << endl << endl;
    // cout << "sum all the contributions" << endl;
    for (int i=0; i<b_len; i++) {
      for (int j=0; j<b_len; j++) {
        // cout << "i, j = " << i << ", " << j << endl;
        // poly_frac_print(&matb[i][j], nroots);
        // cout << "adding:" << endl;
        // poly_frac_print(&proj_mat[i][j], nroots);
        rel_err_poly_frac_add_pf(
          &matb[i][j], &matb[i][j], &proj_mat[i][j],
          roots, NULL,
          wp_bin
        );
        // cout << "summed proj*mat" << endl;
        // poly_frac_print(&matb[i][j], nroots);
        // cout << "adding:" << endl;
        // poly_frac_print(&mat_proj[i][j], nroots);
        rel_err_poly_frac_add_pf(
          &matb[i][j], &matb[i][j], &mat_proj[i][j],
          roots, NULL,
          wp_bin
        );
        // cout << "summed mat*proj" << endl;
        // poly_frac_print(&matb[i][j], nroots);
        // cout << "adding:" << endl;
        // poly_frac_print(&proj_mat_proj[i][j], nroots);
        rel_err_poly_frac_add_pf(
          &matb[i][j], &matb[i][j], &proj_mat_proj[i][j],
          roots, NULL,
          wp_bin
        );
        // cout << "summed proj*mat*proj" << endl;
        // poly_frac_print(&matb[i][j], nroots);
        // cout << "adding:" << endl;
        // poly_frac_print(&der_contr[i][j], nroots);
        rel_err_poly_frac_add_pf(
          &matb[i][j], &matb[i][j], &der_contr[i][j],
          roots, NULL,
          wp_bin
        );
        // cout << "summed derivative contribution" << endl;
        // poly_frac_print(&matb[i][j]);
      }
    }
    // cout << "applied balance" << endl;

    // cumulate balance transformation
    // cout << "cumulate balance transformation" << endl;
    if (num_it == 0) {
      // direct transformation
      poly_frac_rk2_set_mpc_rk2(bal, proj, b_len, b_len, nroots);
      // cout << "converted projector" << endl;
      // poly_frac_rk2_print(bal, b_len, b_len, nroots);
      poly_frac_rk2_mul_root(
        bal, &mpc_one, 1,
        b_len, b_len
      );
      // cout << "multiplied root" << endl;
      poly_frac_rk2_mul_sym_pow(bal, -1, b_len, b_len);
      // cout << "multiplied sym pow" << endl;
      poly_frac_rk2_neg(bal, b_len, b_len);
      // cout << "done negative" << endl;
      for (int i=0; i<b_len; i++) {
        // cout << "i = " << i << endl;
        // poly_frac_print(&bal[i][i]);
        // poly_frac_print(&pf_one, nroots);
        rel_err_poly_frac_add_pf(
          &bal[i][i], &bal[i][i], &pf_one,
          roots, NULL,
          wp_bin
        );
        // poly_frac_print(&bal[i][i]);
      }
      // cout << "added identity" << endl;
      // poly_frac_print(&pf_one, nroots);

      // inverse transformation
      poly_frac_rk2_set_mpc_rk2(inv_bal, proj, b_len, b_len, nroots);
      poly_frac_rk2_mul_root(
        inv_bal, &mpc_one, 1,
        b_len, b_len
      );
      for (int i=0; i<b_len; i++) {
        // cout << "i = " << i << endl;
        // poly_frac_print(&inv_bal[i][i]);
        // poly_frac_print(&pf_one);
        rel_err_poly_frac_add_pf(
          &inv_bal[i][i], &inv_bal[i][i], &pf_one,
          roots, NULL,
          wp_bin
        );
        // poly_frac_print(&inv_bal[i][i]);
      }
    } else {
      // cumulate direct transformation
      // cout << "cumulate direct transformation" << endl;
      rel_err_poly_frac_rk2_mul_mpc_rk2(
        mat_proj, bal, proj, b_len,
        roots,
        wp_bin
      );
      poly_frac_rk2_mul_root(
        mat_proj, &mpc_one, 1,
        b_len, b_len
      );
      poly_frac_rk2_mul_sym_pow(mat_proj, 1, b_len, b_len);
      poly_frac_rk2_neg(mat_proj, b_len, b_len);
      rel_err_poly_frac_rk2_add_pf_rk2(
        bal, bal, mat_proj,
        roots, NULL,
        b_len, b_len,
        wp_bin
      );

      // cumulate inverse transformation
      // cout << "cumulate inverse transformation" << endl;
      rel_err_poly_frac_rk2_mul_mpc_rk2(
        mat_proj, inv_bal, proj, b_len,
        roots,
        wp_bin
      );
      poly_frac_rk2_mul_root(
        mat_proj, &mpc_one, 1,
        b_len, b_len
      );
      rel_err_poly_frac_rk2_add_pf_rk2(
        inv_bal, inv_bal, mat_proj,
        roots, NULL,
        b_len, b_len,
        wp_bin
      );
    }
    
    poly_frac_rk2_prune_rel_tol(matb, wp_bin, b_len, b_len);
    num_it++;
    if (num_it > 1000) {
      perror("too many iterations in diagonal fuchsianization");
      exit(1);
    }
  }
  delete[] A0r;
  delete[] A0i;
  del_rk2_tens(A1, b_len);
  del_rk2_tens(proj, b_len);
  poly_frac_rk2_free(proj_mat, b_len, b_len);
  poly_frac_rk2_free(mat_proj, b_len, b_len);
  poly_frac_rk2_free(proj_mat_proj, b_len, b_len);
  poly_frac_rk2_free(der_contr, b_len, b_len);

  // cout << "TRANSFORMED MATRIX" << endl;
  // poly_frac_rk2_print(matb, b_len, b_len);
  // cout << "BALANCE TRANSFORMATION" << endl;
  // poly_frac_rk2_print(bal, b_len, b_len);
  // cout << "INVERSE BALANCE TRANSFORMATION" << endl;
  // poly_frac_rk2_print(inv_bal, b_len, b_len);
  // getchar();

}


void print_eigenvalues(
  int dim, int num_classes, int *eig_grid, int *eq_class,
  mpc_t *eig_list
) {
  cout << "num classes = " << num_classes << endl;
  cout << "line, grid, class label, int diff, eigenvalues" << endl;
  for (int i=0; i<dim; i++) {
    cout << i << "\t" << eig_grid[i] << "\t" << eq_class[i]  << "\t";
    print_mpc(&eig_list[eq_class[i]]); cout << endl;
  }
}


void pf_NormalizeDiagonal(
  // OUTPUT
  struct poly_frac **tmat, struct poly_frac **inv_tmat,
  mpc_t **eig_list, int *num_classes, int *eq_class, int *eig_grid,
  // IN-OUT
  struct poly_frac **mat_ep,
  // INPUT
  int dim, int nblocks, int **prof, int **sb_grid,
  mpc_t *roots, int nroots,  
  FILE *terminal
) {
  int print = 0;
  int jordan_debug = 0;
  if (print == 1) {
   jordan_debug = -1;
  }
  mpfr_t tmpfr;
  mpfr_init2(tmpfr, wp2);
  int b_len, sb_len, b_len2;
  int p_rank;

  // needed as input for the Jordan decomposition
  int num_ein, *eig_pos;
  mpfr_t *LOre, *LOim;
  // array of matrices to hold the output of the Jordan decomposition of
  // the diagonal blocks
  mpc_t ***jor, ***jor_tmat, ***jor_inv_tmat;
  jor_tmat = new mpc_t**[nblocks];
  jor_inv_tmat = new mpc_t**[nblocks];
  jor = new mpc_t**[nblocks];

  // needed for eigenvalues
  mpc_t *eigenvalues;
  eigenvalues = new mpc_t[dim];
  init_rk1_mpc(eigenvalues, dim);
  int *floor_eig;
  floor_eig = new int[dim];
  int *class_min, *eig_int_diff, offset = 0;
  eig_int_diff = new int[dim];
  int *chain_len = new int[dim];
  // int *eig_grid;
  // eig_grid = new int[dim];
  // int *eq_class;
  // eq_class = new int[dim];
  // int num_classes;

  int wp_bin = - mpfr_log2_int(mpfr_tol);

  if (print) cout << "FUCHSIANIZE DIAGONAL BLOCKS" << endl;
  fprintf(terminal, "diagonal: "); fflush(terminal); usleep(sleep_time);
  for (int b=0; b<nblocks; b++) {
    if (b > 0) {fprintf(terminal, "\033[18D\033[K");}// fflush(terminal); usleep(sleep_time);}
    fprintf(terminal, "block %3d /%3d... ", b, nblocks-1); fflush(terminal); usleep(sleep_time);
    b_len = prof[b][1] - prof[b][0] + 1;
    b_len2 = b_len*b_len;
    if (0 && b == 40) dbg = 1;
    if (print || dbg) {
    cout << "-------------------------" << endl;
    cout << "b = " << b << endl;
    cout << "b_len = " << b_len << endl;
    cout << endl;
    }

    // select block
    poly_frac_rk2_slice(&mat_ep, &mat_ep, prof[b][0], prof[b][1], prof[b][0]);
    poly_frac_rk2_slice(&tmat, &tmat, prof[b][0], prof[b][1], prof[b][0]);
    poly_frac_rk2_slice(&inv_tmat, &inv_tmat, prof[b][0], prof[b][1], prof[b][0]);

    // if (print) cout << "compute avarage exp2" << endl;
    // int block_p_rank = poly_frac_Poinc_rank(mat_ep, b_len, b_len);
    // int num_deg_max = poly_frac_k2_num_deg_max(mat_ep, b_len, b_len);
    // cout << "p_rank = " << block_p_rank << endl;
    // cout << "num_deg_max = " << num_deg_max << endl;
    // if (num_deg_max != -1) {
    //   double *num_exp2_avg = new double[num_deg_max+p_rank+2];
    //   poly_frac_k2_num_exp2_avg(num_exp2_avg, mat_ep, b_len, b_len, block_p_rank, num_deg_max);
    //   if (print) {     
    //     for (int k=-block_p_rank-1; k<=num_deg_max; k++) {
    //       cout << "pow " << k << ": " << num_exp2_avg[k] << endl;
    //     }
    //     cout << endl;
    //   }
    // }

    // // deselect block
    // poly_frac_rk2_unslice(&mat_ep, prof[b][0], prof[b][1], prof[b][0]);
    // poly_frac_rk2_unslice(&tmat, prof[b][0], prof[b][1], prof[b][0]);
    // poly_frac_rk2_unslice(&inv_tmat, prof[b][0], prof[b][1], prof[b][0]);

    // continue;


    // print block
    if (dbg) {
    cout << "block before to_Fuchsian:" << endl;
    poly_frac_rk2_print(mat_ep, b_len, b_len);
    }
    // getchar();

    if (print) cout << "enter to_Fuchsian" << endl;
    pf_to_Fuchsian(tmat, inv_tmat, mat_ep, b_len, roots, nroots);
    if (print) cout << "exit to_Fuchsian:" << endl;
    if (dbg) {
    cout << "block after to_Fuchsian:" << endl;
    poly_frac_rk2_print(mat_ep, b_len, b_len);
    }
    // if (b==19) {
    //   cout << "inv tmat after to_Fuchsian" << endl;
    //   poly_frac_rk2_print(inv_tmat, b_len, b_len);
    // }
    // cout << "MAT:" << endl;
    // poly_frac_rk2_print(mat_ep, b_len, b_len);

    // CONVERT LO TO JORDAN NORMAL FORM
    LOre = new mpfr_t[b_len2];
    init_rk1_mpfr(LOre, b_len2);
    LOim = new mpfr_t[b_len2];
    init_rk1_mpfr(LOim, b_len2);
    // cout << "extract new LO" << endl;
    // poly_frac_rk2_print(mat_ep, b_len, b_len);
    poly_frac_rk2_prune_rel_tol(mat_ep, wp_bin, b_len, b_len);
    // poly_frac_rk2_print(mat_ep, b_len, b_len);
    rel_err_poly_frac_mat_extract_LO(
      LOre, LOim,
      mat_ep, roots, 1,
      b_len, b_len
      , wp_bin
    );
    // cout << "new LO extracted" << endl;
    malloc_rk2_tens(jor[b], b_len, b_len);
    init_rk2_mpc(jor[b], b_len, b_len);
    malloc_rk2_tens(jor_tmat[b], b_len, b_len);
    init_rk2_mpc(jor_tmat[b], b_len, b_len);
    malloc_rk2_tens(jor_inv_tmat[b], b_len, b_len);
    init_rk2_mpc(jor_inv_tmat[b], b_len, b_len);
    if (print || dbg) {
    cout << "enter Jordan decomposition" << endl;
    cout << "input:" << endl;
    for (int i=0; i<b_len; i++) {
      for (int j=0; j<b_len; j++) {
        cout << "(";
        mpfr_out_str(stdout, 10, 0, LOre[i+b_len*j], MPFR_RNDN);
        cout << " ";
        mpfr_out_str(stdout, 10, 0, LOim[i+b_len*j], MPFR_RNDN);
        cout << ")" << endl;
      }
      cout << endl;
    }
    }
    jordan_wrp(jor[b], jor_tmat[b], &eig_pos, &num_ein, LOre, LOim, b_len, b_len, jordan_debug);
    // if (b == 40) {
    //   // mpc_rk2_to_file((char*)"jor_tmat.txt", jor_tmat[b], b_len, b_len);
    //   mpc_rk2_from_file((char*)"jor_tmat.txt", jor_tmat[b], b_len, b_len);
    // }
    mpc_rk2_prune_abs_tol(jor[b], b_len, b_len);
    mpc_rk2_prune_abs_tol(jor_tmat[b], b_len, b_len);
    if (print || dbg) {
    cout << "exit Jordan decomposition" << endl;
    cout << "output:" << endl;
    print_rk2_mpc(jor[b], b_len, b_len);
    cout << "jor tmat:" << endl;
    print_rk2_mpc(jor_tmat[b], b_len, b_len);
    }
    
    // ANALYZE EIGENVALUES
    //// collect eigenvalues from the diagonal of the Jordan block
    // cout << "eigenvalues:" << endl;
    for (int i=0; i<b_len; i++) {
      mpc_set(eigenvalues[prof[b][0]+i], jor[b][i][i], MPFR_RNDN);
      // print_mpc(&jor[b][i][i]);
      // cout << endl;
      // also initialize to zero the eigenvalues grid and the chain lenghts
      eig_grid[prof[b][0]+i] = 0;
      chain_len[prof[b][0]+i] = 0;
    }
    // also initialize to zero the eigenvalues grid and the chain lenghts
    for (int npos, i=0; i<num_ein; i++) {
      eig_grid[prof[b][0]+eig_pos[i]] = 1;
      if (i < num_ein-1) {
        npos = eig_pos[i + 1];
      } else {
        npos = b_len;
      }
      chain_len[prof[b][0]+eig_pos[i]] = npos - eig_pos[i];
    }

    free(eig_pos);
    delete[] LOre;
    delete[] LOim;

    // deselect block
    poly_frac_rk2_unslice(&mat_ep, prof[b][0], prof[b][1], prof[b][0]);
    poly_frac_rk2_unslice(&tmat, prof[b][0], prof[b][1], prof[b][0]);
    poly_frac_rk2_unslice(&inv_tmat, prof[b][0], prof[b][1], prof[b][0]);
  }
  fprintf(terminal, "\033[18D\033[K"); fflush(terminal); usleep(sleep_time);
  fprintf(terminal, "\033[10D\033[K"); fflush(terminal); usleep(sleep_time);
  dbg = 0;

  //////
  // NORMALIZE EIGENVALUES
  //////
  if (print) cout << "grouping eigenvalues..." << endl;
  // exit(0);
  group_eigenvalues(
    class_min, num_classes,
    dim, eigenvalues, eig_grid,
    eq_class, eig_int_diff
  );
  // cout << "finshed grouping" << endl;
  // cout << "INPUT OF NORMALIZATION" << endl;
  // cout << "line, grid, class label, min pos, int diff, eigenvalues" << endl;
  // for (int i=0; i<dim; i++) {
  //   cout << i << "\t" << eig_grid[i] << "\t" << eq_class[i] << "\t" << class_min[eq_class[i]] << "\t" << eig_int_diff[i] << "\t";
  //   print_mpc(&eigenvalues[i]);
  //   cout << endl;
  // }
  // cout << endl;
  for (int i=0; i<dim; i++) {
    mpfr_add(tmpfr, mpc_realref(eigenvalues[i]), mpfr_tol, MPFR_RNDN);
    floor_eig[i] = mpfr_get_si(tmpfr, MPFR_RNDD); 
  }

  // needed for permutations
  int *eig_perm, *eig_grid_copy, *eq_class_copy, *eig_int_diff_copy, *chain_len_copy, *floor_eig_copy;
  mpc_t mpc_diff;
  mpc_init3(mpc_diff, wp2, wp2);
  mpfr_t mpfr_abs_diff;
  mpfr_init2(mpfr_abs_diff, wp2);

  // struct poly_frac **shearing, **inv_shearing;
  struct poly_frac **pf_tmp_mat;
  struct poly_frac pf_one_over_eta;
  poly_frac_build(&pf_one_over_eta);
  poly_frac_set_ui(&pf_one_over_eta, 1, nroots);
  poly_frac_mul_sym_pow(&pf_one_over_eta, &pf_one_over_eta, -1);
  struct poly_frac pf_minus_one_over_eta;
  poly_frac_build(&pf_minus_one_over_eta);
  poly_frac_set_pf(&pf_minus_one_over_eta, &pf_one_over_eta);
  poly_frac_neg(&pf_minus_one_over_eta);
  int eig_normalized;
  int it;
  int *new_chain_len, *already_selected;
  if (print) cout << "NORMALIZE EIGENVALUES" << endl;
  // exit(0);
  // mpfr_mul_d(mpfr_tol, mpfr_tol, 1e70, MPFR_RNDN);
  for (int b=0; b<nblocks; b++) {
    b_len = prof[b][1] - prof[b][0] + 1;
    b_len2 = b_len*b_len;
    if (0 && b == 40) dbg = 1;
    if (print || dbg) {
    cout << "-------------------------" << endl;
    cout << "b = " << b << endl;
    cout << "b_len = " << b_len << endl;
    cout << endl; cout << "eigenvalues:" << endl;
    cout << "line, grid, class label, min pos, int diff, eigenvalues" << endl;
    for (int i, ip=0; ip<b_len; ip++) {
      i = offset+ip;
      cout << i << "\t" << eig_grid[i] << "\t" << eq_class[i] << "\t" << class_min[eq_class[i]] << "\t" << eig_int_diff[i] << "\t";
      print_mpc(&eigenvalues[i]); cout << endl;
    }
    }

    poly_frac_rk2_slice(&mat_ep, &mat_ep, prof[b][0], prof[b][1], prof[b][0]);
    poly_frac_rk2_slice(&tmat, &tmat, prof[b][0], prof[b][1], prof[b][0]);
    poly_frac_rk2_slice(&inv_tmat, &inv_tmat, prof[b][0], prof[b][1], prof[b][0]);
    
    //////
    // NORMALIZE BLOCK
    //////
    malloc_rk2_tens(pf_tmp_mat, b_len, b_len);
    poly_frac_rk2_build(pf_tmp_mat, b_len, b_len);
    eig_normalized = 1;
    it = 0;
    while (1) {
      // cout << "eigen floors:" << endl;
      // for (int i=0; i<b_len; i++) {
      //   cout << "i: " << i << ", floor = " << floor_eig[offset+i] << endl;
      // }
      // first check if block is already normalized
      eig_normalized = 1;
      for (int i=0; i<b_len; i++) {
        // if (eig_int_diff[offset+i] != 0) {
        if (floor_eig[offset+i] != 0) {
          eig_normalized = 0;
          break;
        }
      }
      if (eig_normalized) {
        // turn the Leading Order to Jordan Normal Form
        //// apply Jordan transformation
        // cout << "mul jor on the right" << endl;
        rel_err_poly_frac_rk2_mul_mpc_rk2(
          pf_tmp_mat, mat_ep, jor_tmat[b], b_len,
          roots,
          wp_bin
        );
        //// cumulate Jordan transformation
        // cout << "cumulate Jordan transformation" << endl;
        rel_err_poly_frac_rk2_mul_mpc_rk2(
          tmat, tmat, jor_tmat[b], b_len,
          roots,
          wp_bin
        );
        //// apply inverse Jordan transformation
        mp_inverse(jor_inv_tmat[b], jor_tmat[b], b_len);
        mpc_rk2_prune_abs_tol(jor_inv_tmat[b], b_len, b_len);
        //// cumulate inverse Jordan transformation
        rel_err_mpc_rk2_mul_poly_frac_rk2(
          mat_ep, jor_inv_tmat[b], pf_tmp_mat, b_len,
          roots,
          wp_bin
        );
        // if (b == 19) {
        //   cout << "tmat inv before jor:" << endl;
        //   poly_frac_rk2_print(inv_tmat, b_len, b_len);
        //   cout << "inv jor tmat:" << endl;
        //   print_rk2_mpc(jor_inv_tmat[b], b_len, b_len);
        // }
        rel_err_mpc_rk2_mul_poly_frac_rk2(
          inv_tmat, jor_inv_tmat[b], inv_tmat, b_len,
          roots,
          wp_bin
        );
        // if (b == 19) {
        //   cout << "tmat inv after jor:" << endl;
        //   poly_frac_rk2_print(inv_tmat, b_len, b_len);
        // }

        // store new eigenvalues
        for (int i=0; i<b_len; i++) {
          mpc_set(eigenvalues[offset+i], jor[b][i][i], MPFR_RNDN);
        }

        break;
      }

      it++;
      if (it>1000) {
        perror("too many iterations in eigenvalue normalization");
        exit(1);
      }
      if (print) cout << "it = " << it << endl;

      //////
      // SHEARING TRANSFORMATION
      //////
      // // build shearing transformation
      // malloc_rk2_tens(shearing, b_len, b_len);
      // poly_frac_rk2_build(shearing, b_len, b_len);
      // malloc_rk2_tens(inv_shearing, b_len, b_len);
      // poly_frac_rk2_build(inv_shearing, b_len, b_len);
      // poly_frac_rk2_set_mpc_rk2(
      //   shearing, jor_tmat[b], b_len, b_len, nroots
      // );
      // poly_frac_rk2_set_mpc_rk2(
      //   inv_shearing, jor_inv_tmat[b], b_len, b_len, nroots
      // );

      // apply Jordan transformation
      if (print || dbg) {
      cout << "apply Jordan transformation" << endl;
      cout << "block before: " << endl;
      poly_frac_rk2_print(mat_ep, b_len, b_len);
      cout << "jordan before: " << endl;
      print_rk2_mpc(jor_tmat[b], b_len, b_len);
      cout << "mat*jor" << endl;
      }
      rel_err_poly_frac_rk2_mul_mpc_rk2(
        pf_tmp_mat, mat_ep, jor_tmat[b], b_len,
        roots,
        wp_bin
      );
      if (dbg) {
      cout << "intermediate block 1: " << endl;
      poly_frac_rk2_print(pf_tmp_mat, b_len, b_len);
      // cumulate Jordan transformation
      cout << "tmat*jor_tmat" << endl;
      }
      rel_err_poly_frac_rk2_mul_mpc_rk2(
        tmat, tmat, jor_tmat[b], b_len,
        roots,
        wp_bin
      );
      if (print || dbg) {
      cout << "cumulated tmat: " << endl;
      poly_frac_rk2_print(tmat, b_len, b_len);
      // apply inverse Jordan transformation
      cout << "apply inverse Jordan transformation" << endl;
      }
      mp_inverse(jor_inv_tmat[b], jor_tmat[b], b_len);
      mpc_rk2_prune_abs_tol(jor_inv_tmat[b], b_len, b_len);
      if (print || dbg) {
      cout << "inverse jordan before: " << endl;
      print_rk2_mpc(jor_inv_tmat[b], b_len, b_len);
      }
      rel_err_mpc_rk2_mul_poly_frac_rk2(
        mat_ep, jor_inv_tmat[b], pf_tmp_mat, b_len,
        roots,
        wp_bin
      );
      if (print || dbg) {
      cout << "intermediate block 2: " << endl;
      poly_frac_rk2_print(mat_ep, b_len, b_len);
      // cumulate inverse Jordan transformation
      // cout << "jor_inv_tmat*inv_tmat" << endl;
      }
      rel_err_mpc_rk2_mul_poly_frac_rk2(
        inv_tmat, jor_inv_tmat[b], inv_tmat, b_len,
        roots,
        wp_bin
      );
      if (print || dbg) {
      cout << "cumulated inv_tmat: " << endl;
      poly_frac_rk2_print(inv_tmat, b_len, b_len);
      cout << "block before: " << endl;
      poly_frac_rk2_print(mat_ep, b_len, b_len);
      }

      // cout << endl; cout << "current eigenvalues:" << endl;
      // for (int i=0; i<b_len; i++) {
      //   print_mpc(&jor[b][i][i]); cout << endl;
      // }
      // cout << endl;

      // multiply times eta (eta^-1) the columns (rows) of the (inverse)
      // Jordan transformation indicated by the eigenvalue decomposition,
      // i.e. those associated with positive integer differences
      int shearing_direction;
      if (max_ints(floor_eig+offset, b_len) > 0) {
        shearing_direction = -1;
        if (print || dbg) cout << "shear dir = " << shearing_direction << endl;
        for (int j=0; j<b_len; j++) {
          if (print || dbg) cout << "j = " << j << endl;
          if (print || dbg) cout << eig_int_diff[offset+j] << endl;
          if (print || dbg) cout << "floor: " << floor_eig[offset+j] << endl;
          // if (eig_int_diff[offset+j] != 0) {
            if (floor_eig[offset+j] > 0) {
            if (print || dbg) cout << "selected" << endl;
            // lower by 1 the value of the decomposition
            // eig_int_diff[offset+j]--;
            floor_eig[offset+j]--;
            // and of the eigenvalue as well
            mpc_sub_ui(eigenvalues[offset+j], eigenvalues[offset+j], 1, MPFR_RNDN);

            for (int i=0; i<b_len; i++) {
              // multiply column times eta
              poly_frac_mul_sym_pow(&mat_ep[i][j], &mat_ep[i][j], 1);
              poly_frac_mul_sym_pow(&tmat[i][j], &tmat[i][j], 1);

              // multiply row of the inverse times 1/eta
              poly_frac_mul_sym_pow(&mat_ep[j][i], &mat_ep[j][i], -1);
              poly_frac_mul_sym_pow(&inv_tmat[j][i], &inv_tmat[j][i], -1);
            }

            // derivative contribution
            rel_err_poly_frac_add_pf(
              &mat_ep[j][j], &mat_ep[j][j], &pf_minus_one_over_eta,
              roots, NULL,
              wp_bin
            );
            if (print || dbg) {
              poly_frac_print(&mat_ep[j][j]);
            }
          }
        }
      } else {
        shearing_direction = 1;
        if (print || dbg) cout << "shear dir = " << shearing_direction << endl;
        for (int j=0; j<b_len; j++) {
          if (print || dbg) cout << "j = " << j << endl;
          if (print || dbg) cout << eig_int_diff[offset+j] << endl;
          if (print || dbg) cout << "floor: " << floor_eig[offset+j] << endl;
          // if (eig_int_diff[offset+j] != 0) {
            if (floor_eig[offset+j] < 0) {
            if (print || dbg) cout << "selected" << endl;
            // lower by 1 the value of the decomposition
            // eig_int_diff[offset+j]--;
            floor_eig[offset+j]++;
            // and of the eigenvalue as well
            mpc_add_ui(eigenvalues[offset+j], eigenvalues[offset+j], 1, MPFR_RNDN);

            for (int i=0; i<b_len; i++) {
              // multiply column times eta
              poly_frac_mul_sym_pow(&mat_ep[i][j], &mat_ep[i][j], -1);
              poly_frac_mul_sym_pow(&tmat[i][j], &tmat[i][j], -1);

              // multiply row of the inverse times 1/eta
              poly_frac_mul_sym_pow(&mat_ep[j][i], &mat_ep[j][i], 1);
              poly_frac_mul_sym_pow(&inv_tmat[j][i], &inv_tmat[j][i], 1);
            }

            // derivative contribution
            // cout << "addend1: " << endl;
            // poly_frac_print(&mat_ep[j][j]);
            // cout << "addend2: " << endl;
            // poly_frac_print(&pf_one_over_eta);
            rel_err_poly_frac_add_pf(
              &mat_ep[j][j], &mat_ep[j][j], &pf_one_over_eta,
              roots, NULL,
              wp_bin
            );
            if (print || dbg) {
              poly_frac_print(&mat_ep[j][j]);
            }
          }
        }
      }
      if (print || dbg) {
      cout << "block after: " << endl;
      poly_frac_rk2_print(mat_ep, b_len, b_len);
      }
      // getchar();


      // CONVERT LO TO JORDAN NORMAL FORM
      LOre = new mpfr_t[b_len2];
      init_rk1_mpfr(LOre, b_len2);
      LOim = new mpfr_t[b_len2];
      init_rk1_mpfr(LOim, b_len2);
      // cout << "compute new leading order" << endl;
      poly_frac_rk2_prune_rel_tol(mat_ep, wp_bin, b_len, b_len);
      rel_err_poly_frac_mat_extract_LO(
        LOre, LOim,
        mat_ep, roots, 1,
        b_len, b_len
        , wp_bin
      );
      // cout << "done" << endl;
      if (print) {
      cout << "enter Jordan decomposition" << endl;
      cout << "input:" << endl;
      for (int i=0; i<b_len; i++) {
        for (int j=0; j<b_len; j++) {
          cout << "(";
          mpfr_out_str(stdout, 10, 0, LOre[i+b_len*j], MPFR_RNDN);
          cout << ", ";
          mpfr_out_str(stdout, 10, 0, LOim[i+b_len*j], MPFR_RNDN);
          cout << ")" << endl;
        }
        cout << endl;
      }
      }
      jordan_wrp(jor[b], jor_tmat[b], &eig_pos, &num_ein, LOre, LOim, b_len, b_len, jordan_debug);
      // if (b == 40) {
      //   // mpc_rk2_to_file((char*)"jor_tmat2.txt", jor_tmat[b], b_len, b_len);
      //   mpc_rk2_from_file((char*)"jor_tmat2.txt", jor_tmat[b], b_len, b_len);
      // }
      mpc_rk2_prune_abs_tol(jor[b], b_len, b_len);
      mpc_rk2_prune_abs_tol(jor_tmat[b], b_len, b_len);
      if (print) {
      cout << "exit Jordan decomposition" << endl;
      cout << "output:" << endl;
      print_rk2_mpc(jor[b], b_len, b_len);
      cout << "jor tmat:" << endl;
      print_rk2_mpc(jor_tmat[b], b_len, b_len);
      cout << "new eigenvalues:" << endl;
      for (int j=0; j<num_ein; j++) {
        cout << "j = " << j << ": ";
        print_mpc(&jor[b][eig_pos[j]][eig_pos[j]]);
        cout << endl;
      }
      }

      // //// find the permutation of the eigenvectors
      // if (b_len > 1) {
      //   // cout << "find the permutation of the eigenvectors" << endl;
      //   eig_perm = new int[b_len];
      //   eq_class_copy = new int[b_len];
      //   eig_grid_copy = new int[b_len];
      //   eig_int_diff_copy = new int[b_len];
      //   chain_len_copy = new int[b_len];
      //   already_selected = new int[num_ein];
      //   floor_eig_copy = new int[b_len];
      //   //// get new chain lenghts
      //   new_chain_len = new int[num_ein];
      //   for (int j=0; j<num_ein-1; j++) {
      //     new_chain_len[j] = eig_pos[j+1] - eig_pos[j];
      //   }
      //   new_chain_len[num_ein-1] = b_len - eig_pos[num_ein-1];
      //   for (int j=0; j<num_ein; j++) {
      //     already_selected[j] = 0;
      //   }

      //   for (int i=0; i<b_len; i++) {
      //     // cout << "i = " << i << endl;
      //     if (eig_grid[offset+i] == 0) {
      //       continue;
      //     }
      //     // store old_values
      //     eig_grid_copy[i] = eig_grid[offset+i];
      //     eq_class_copy[i] = eq_class[offset+i];
      //     eig_int_diff_copy[i] = eig_int_diff[offset+i];
      //     chain_len_copy[i] = chain_len[offset+i];
      //     floor_eig_copy[i] = floor_eig[offset+i];
      //     //// (and set to zero the new grid and chain lengths)
      //     eig_grid[offset+i] = 0;
      //     chain_len[offset+i] = 0;
          
      //     // consider the i-th eigenvalue and find its new position
      //     // cout << "consider the i-th eigenvalue and find its new position" << endl;
      //     for (int j=0; j<num_ein; j++) {
      //       // cout << "j = " << j << endl;
      //       // if chain length is different, skip
      //       if (chain_len_copy[i] != new_chain_len[j]) {
      //         continue;
      //       }
      //       // if already selected, skip
      //       if (already_selected[j]) {
      //         continue;
      //       }
      //       mpc_sub(mpc_diff, eigenvalues[offset+i], jor[b][eig_pos[j]][eig_pos[j]], MPFR_RNDN);
      //       mpc_abs(mpfr_abs_diff, mpc_diff, MPFR_RNDN);
      //       if (mpfr_less_p(mpfr_abs_diff, mpfr_tol)) {
      //         eig_perm[j] = i;
      //         already_selected[j] = 1;
      //         // if it is a minimum, update its position
      //         if (class_min[eq_class[offset+i]] == offset+i) {
      //           class_min[eq_class[offset+i]] = offset+eig_pos[j];
      //         }
      //         break;
      //       }
      //     }
      //   }
      //   // cout << "apply the permutation" << endl;
      //   //// apply the permutation
      //   for (int j=0; j<num_ein; j++) {
      //     // cout << "j = " << j << endl;
      //     // move the entire Jordan chain in the proper position
      //     for (int jp=0; jp<new_chain_len[j]; jp++) {
      //       // cout << "jp = " << jp << endl;
      //       mpc_set(eigenvalues[offset+eig_perm[j]+jp], jor[b][eig_pos[j]+jp][eig_pos[j]+jp], MPFR_RNDN);
      //       // mpfr_add(tmpfr, mpc_realref(jor[b][eig_pos[j]+jp][eig_pos[j]+jp]), mpfr_tol, MPFR_RNDN);
      //       // floor_eig[offset+eig_perm[j]+jp] = mpfr_get_si(tmpfr, MPFR_RNDD);
      //     }
      //     // permute all the classification properties
      //     // cout << "permute all the classification properties" << endl;
      //     eig_grid[offset+eig_perm[j]] = eig_grid_copy[eig_pos[j]];
      //     eq_class[offset+eig_perm[j]] = eq_class_copy[eig_pos[j]];
      //     eig_int_diff[offset+eig_perm[j]] = eig_int_diff_copy[eig_pos[j]];
      //     chain_len[offset+eig_perm[j]] = chain_len_copy[eig_pos[j]];
      //     floor_eig[offset+eig_perm[j]] = floor_eig_copy[eig_pos[j]];

      //     // assign properties to non-pure eigenstates
      //     // cout << "assign properties to non-pure eigenstates" << endl;
      //     for (int jp=1; jp<new_chain_len[jp]; jp++) {
      //         eq_class[offset+eig_perm[j]+jp] = eq_class[offset+eig_perm[j]];
      //         class_min[offset+eig_perm[j]+jp] = class_min[offset+eig_perm[j]];
      //         eig_int_diff[offset+eig_perm[j]+jp] = eig_int_diff[offset+eig_perm[j]];
      //         floor_eig[offset+eig_perm[j]+jp] = floor_eig_copy[offset+eig_perm[j]];
      //     }
      //   }
      //   // ////// deal with non-pure eigenstates
      //   // for (int i=0; i<b_len; i++) {
      //   //   if (eig_grid[offset+i] == 0) {
      //   //     for (int j=i-1; j>=0; j--) {
      //   //       if (eig_grid[j] == 1 || eig_grid[j] < 0) {
      //   //         eq_class[offset+i] = eq_class[offset+j];
      //   //         class_min[offset+i] = class_min[offset+j];
      //   //         eig_int_diff[offset+i] = eig_int_diff[offset+j];
      //   //         break;
      //   //       }
      //   //     }
      //   //   }
      //   // }
      // }

      if (b_len > 1) {
        // EIGENVALUES GRID AND CHAIN LENS
        for (int i=0; i<b_len; i++) {
          eig_grid[prof[b][0]+i] = 0;
          chain_len[prof[b][0]+i] = 0;
        }
        for (int npos, i=0; i<num_ein; i++) {
          eig_grid[prof[b][0]+eig_pos[i]] = 1;
          if (i < num_ein-1) {
            npos = eig_pos[i + 1];
          } else {
            npos = b_len;
          }
          chain_len[prof[b][0]+eig_pos[i]] = npos - eig_pos[i];
        }

        for (int i=0; i<num_ein-1; i++) {
          mpfr_add(tmpfr, mpc_realref(jor[b][eig_pos[i]][eig_pos[i]]), mpfr_tol, MPFR_RNDN);
          floor_eig[offset+eig_pos[i]] = mpfr_get_si(tmpfr, MPFR_RNDD);
          for (int j=1; j<eig_pos[i+1]-eig_pos[i]; j++) {
            floor_eig[offset+eig_pos[i]+j] = floor_eig[offset+eig_pos[i]];
          }
        }
        mpfr_add(tmpfr, mpc_realref(jor[b][eig_pos[num_ein-1]][eig_pos[num_ein-1]]), mpfr_tol, MPFR_RNDN);
        floor_eig[offset+eig_pos[num_ein-1]] = mpfr_get_si(tmpfr, MPFR_RNDD);
        for (int j=1; j<b_len-eig_pos[num_ein-1]; j++) {
          floor_eig[offset+eig_pos[num_ein-1]+j] = floor_eig[offset+eig_pos[num_ein-1]];
        }
      }
    }
    // cout << "FINAL BLOCK" << endl;
    // poly_frac_rk2_print(mat_ep, b_len, b_len);
    // cout << "CUMULATED TMAT" << endl;
    // poly_frac_rk2_print(tmat, b_len, b_len);
    // cout << "CUMULATED INV TMAT" << endl;
    // poly_frac_rk2_print(inv_tmat, b_len, b_len);
    // getchar();

    // p_rank = poly_frac_Poinc_rank(mat_ep, b_len, b_len);
    // if (p_rank > 0) {
    //   cout << "ERROR: Poinc. rank = " << p_rank;
    //   cout << " after normalization of block b = " << b << endl;
    //   exit(0);
    // }

    // deselct blocks
    poly_frac_rk2_unslice(&mat_ep, prof[b][0], prof[b][1], prof[b][0]);
    poly_frac_rk2_unslice(&tmat, prof[b][0], prof[b][1], prof[b][0]);
    poly_frac_rk2_unslice(&inv_tmat, prof[b][0], prof[b][1], prof[b][0]);

    // FREE
    poly_frac_rk2_free(pf_tmp_mat, b_len, b_len);
    del_rk2_tens(pf_tmp_mat, b_len);

    if (!eig_normalized) {
      delete[] LOre;
      delete[] LOim;
      delete[] eig_pos;
      // if (b_len > 1) {
      // delete[] eig_perm;
      // delete[] eq_class_copy;
      // delete[] eig_grid_copy;
      // delete[] eig_int_diff_copy;
      // delete[] new_chain_len;
      // delete[] already_selected;
      // delete[] chain_len_copy;
      // }
    }

    offset += b_len;
  }
  dbg = 0;
  
  // check normalization
  // cout << endl; cout << "OUTPUT OF NORMALIZATION" << endl;
  // cout << "line, grid, class label, min pos, int diff, eigenvalues" << endl;
  // for (int i=0; i<dim; i++) {
  //   cout << i << "\t" << eig_grid[i] << "\t" << eq_class[i] << "\t" << class_min[eq_class[i]] << "\t" << eig_int_diff[i] << "\t";
  //   print_mpc(&eigenvalues[i]);
  //   cout << endl;
  // }
  // cout << "grouping on the output" << endl;
  // #2BD: possiamo rimuoverlo? In teoria si', controlliamo
  // dopo aver aggiunto shearing shift up bisogna prima modificare il codice
  // e soltando dopo rimuovere questo grouping
  group_eigenvalues(
    class_min, num_classes,
    dim, eigenvalues, eig_grid,
    eq_class, eig_int_diff
  );
  if (print) {
  cout << endl; cout << "OUTPUT OF FINAL GROUPING" << endl;
  cout << "num classes = " << *num_classes << endl;
  cout << "line, grid, class label, min pos, int diff, eigenvalues" << endl;
  for (int i=0; i<dim; i++) {
    cout << i << "\t" << eig_grid[i] << "\t" << eq_class[i] << "\t" << class_min[eq_class[i]] << "\t" << eig_int_diff[i] << "\t";
    print_mpc(&eigenvalues[i]);
    cout << endl;
  }
  }

  // PREPARE OUTPUT EIGENVALUES
  (*eig_list) = new mpc_t[*num_classes];
  init_rk1_mpc(*eig_list, *num_classes);
  for (int i=0; i<*num_classes; i++) {
    mpc_set((*eig_list)[i], eigenvalues[class_min[i]], MPFR_RNDN);
  }

  // FREE
  for (int b=0; b<nblocks; b++) {
    b_len = prof[b][1] - prof[b][0] + 1;
    mpc_rk2_clear(jor[b], b_len, b_len);
    del_rk2_tens(jor[b], b_len);
    mpc_rk2_clear(jor_tmat[b], b_len, b_len);
    del_rk2_tens(jor_tmat[b], b_len);
    mpc_rk2_clear(jor_inv_tmat[b], b_len, b_len);
    del_rk2_tens(jor_inv_tmat[b], b_len);
  }
  delete[] jor;
  delete[] jor_tmat;
  delete[] jor_inv_tmat;

  if (print) cout << endl << "TRANSFORM SUB-DIAGONAL BLOCKS" << endl;
  // TRANSFORM SUB-DIAGONAL BLOCK
  // poly_frac_rk2_print(mat_ep, dim, dim);
  // getchar();
  // struct poly_frac **matsb, **prev_bal, **inv_bal;
  for (int b=0; b<nblocks; b++) {
    b_len = prof[b][1] - prof[b][0] + 1;
    if (0 && b == 40) dbg = 1;
    if (print || dbg) {
    cout << "-------------------------" << endl;
    cout << "b = " << b << endl;
    cout << "b_len = " << b_len << endl;
    cout << endl;
    }

    // select current balance
    poly_frac_rk2_slice(&inv_tmat, &inv_tmat, prof[b][0], prof[b][1], prof[b][0]);
    // cout << "input inv bal" << endl;
    // poly_frac_rk2_print(inv_tmat, b_len, b_len);
    for (int sb=0; sb<b; sb++) {
      if (sb_grid[b][sb] == 0) {
        continue;
      }
      sb_len = prof[sb][1] - prof[sb][0] + 1;
      // if (b == 40 && sb == 27) dbg = 1;
      if (print || dbg) {
      cout << "-------------------------" << endl;
      cout << "sb = " << sb << endl;
      cout << "sb_len = " << sb_len << endl;
      cout << endl;
      }

      // select sub-block, previous balance
      poly_frac_rk2_slice(&mat_ep, &mat_ep, prof[b][0], prof[b][1], prof[sb][0]);
      // cout << "selected sub-diagonal block" << endl;
      // poly_frac_rk2_print(mat_ep, b_len, sb_len);
      poly_frac_rk2_slice(&tmat, &tmat, prof[sb][0], prof[sb][1], prof[sb][0]);
      
      // apply transformation
      // matsb -> inv_bal * matsb * prev_bal
      if (dbg) {
      cout << "matsb -> inv_bal * matsb * prev_bal" << endl;
      cout << "input block:" << endl;
      poly_frac_rk2_print_to_math(mat_ep, b_len, sb_len, roots);
      cout << "multiply on the right" << endl;
      cout << "input prev bal:" << endl;
      poly_frac_rk2_print_to_math(tmat, sb_len, sb_len, roots);
      }
      rel_err_poly_frac_rk2_mul_pf_rk2(mat_ep, mat_ep, tmat, b_len, sb_len, sb_len, roots, wp_bin);
      if (dbg) {      
      cout << "out:" << endl;
      poly_frac_rk2_print_to_math(mat_ep, b_len, sb_len, roots);
      cout << "multiply on the left" << endl;
      cout << "input inv tmat:" << endl;
      poly_frac_rk2_print_to_math(inv_tmat, b_len, b_len, roots);
      }
      rel_err_poly_frac_rk2_mul_pf_rk2(mat_ep, inv_tmat, mat_ep, b_len, b_len, sb_len, roots, wp_bin);
      if (dbg) {
      cout << "out:" << endl;
      poly_frac_rk2_print_to_math(mat_ep, b_len, sb_len, roots);
      }

      // deselect sub-block, previous balance
      poly_frac_rk2_unslice(&mat_ep, prof[b][0], prof[b][1], prof[sb][0]);
      poly_frac_rk2_unslice(&tmat, prof[sb][0], prof[sb][1], prof[sb][0]);
      // dbg = 0;
    }
    dbg = 0;
    // deselct current balance
    poly_frac_rk2_unslice(&inv_tmat, prof[b][0], prof[b][1], prof[b][0]);
  }

  // set to zero top off-diagonal elements
  // cout << "set to zero top off-diagonal elements" << endl;
  for (int b=0; b<nblocks; b++) {
    // cout << "b = " << b << endl;
    for (int i=prof[b][0]; i<=prof[b][1]; i++) {
      for (int j=prof[b][1]+1; j<dim; j++) {
        // cout << "i, j = " << i << ", " << j << endl;
        poly_frac_set_zero(&tmat[i][j]);
        poly_frac_set_zero(&inv_tmat[i][j]);
      }
    }
  }

}


void pf_apply_sub_diagonal_p_reduction(
  // IN-OUT
  struct poly_frac **mat,
  // INPUT
  int p_rank, int b, int sb, int nblocks, int **prof,
  mpc_t *sb_sys_sol, mpc_t *roots, int nroots
) {
  int b_len, sb_len, h_len, k_len;
  b_len = prof[b][1] - prof[b][0] + 1;
  sb_len = prof[sb][1] - prof[sb][0] + 1;

  struct poly_frac pf_tmp;
  poly_frac_build(&pf_tmp);

  int wp_bin = - mpfr_log2_int(mpfr_tol);

  // apply inverse transformation on the left
  // (modify only the elements of the b-th row up to the sb-th element)
  for (int k=0; k<=sb; k++) {
    k_len = prof[k][1] - prof[k][0] + 1;
    // cout << "k = " << k << endl;
    // update sub-diagonal block (b, k)
    for (int i=0; i<b_len; i++) {
      for (int j=0; j<k_len; j++) {
        // cout << "i, j = " << i << ", " << j << endl;
        for (int a=0; a<sb_len; a++) {
          // cout << "a = " << a << endl;
          poly_frac_mul_sym_pow(&pf_tmp, &mat[prof[sb][0] + a][prof[k][0] + j], -p_rank);
          poly_frac_neg(&pf_tmp);
          rel_err_poly_frac_add_pf(
            &mat[prof[b][0] + i][prof[k][0] + j],
            &mat[prof[b][0] + i][prof[k][0] + j],
            &pf_tmp,
            roots,
            &sb_sys_sol[i*sb_len + a],
            wp_bin
          );
          // poly_frac_print(&mat[prof[b][0] + i][prof[k][0] + j]);
        }
      }
    }
  }
  

  // apply direct transformation on the right
  // (modify only the elements of the sb-th column starting from b-th element)
  for (int h=b; h<nblocks; h++) {
    h_len = prof[h][1] - prof[h][0] + 1;
    // cout << "h = " << h << endl;
    // update sub-diagonal block (h, sb)
    for (int i=0; i<h_len; i++) {
      for (int j=0; j<sb_len; j++) {
        // cout << "i, j = " << i << ", " << j << endl;
        for (int a=0; a<b_len; a++) {
          // cout << "a = " << a << endl;
          poly_frac_mul_sym_pow(&pf_tmp, &mat[prof[h][0] + i][prof[b][0] + a], -p_rank);
          rel_err_poly_frac_add_pf(
            &mat[prof[h][0] + i][prof[sb][0] + j],
            &mat[prof[h][0] + i][prof[sb][0] + j],
            &pf_tmp,
            roots,
            &sb_sys_sol[a*sb_len + j],
            wp_bin
          );
          // poly_frac_print(&mat[prof[h][0] + i][prof[sb][0] + j]);
        }
      }
    }
  }

  // finally, update the (b, sb) sub-diagonal block
  for (int i=0; i<b_len; i++) {
    for (int j=0; j<sb_len; j++) {
      // cout << "i, j = " << i << ", " << j << endl;
      // cout << "mpc:" << endl;
      // print_mpc(&sb_sys_sol[i*sb_len + j]);
      // cout << endl;
      poly_frac_set_mpc(&pf_tmp, &sb_sys_sol[i*sb_len + j], nroots);
      poly_frac_mul_sym_pow(&pf_tmp, &pf_tmp, -p_rank-1);
      // cout << "converted to poly and diveded by pow:" << endl;
      // poly_frac_print(&pf_tmp);
      poly_frac_mul_ui(&pf_tmp, &pf_tmp, p_rank);
      // cout << "multiplied by p_rank" << endl;
      // poly_frac_print(&pf_tmp);
      // cout << "1st addend:" << endl;
      // poly_frac_print(&mat[prof[b][0] + i][prof[sb][0] + j]);
      rel_err_poly_frac_add_pf(
        &mat[prof[b][0] + i][prof[sb][0] + j],
        &mat[prof[b][0] + i][prof[sb][0] + j],
        &pf_tmp,
        roots,
        NULL,
        wp_bin
      );
      // poly_frac_print(&mat[prof[b][0] + i][prof[sb][0] + j]);
    }
  }
}


void pf_to_Fuchsian_global(
  // IN-OUT
  struct poly_frac **mat_ep,
  struct poly_frac **tmat, struct poly_frac **inv_tmat,
  // INPUT
  int dim, int nblocks, int **prof, int **sb_grid,
  mpc_t *roots, int nroots,
  FILE *terminal
) {
  int print = 0;
  int p_rank, big_dim, b_len, sb_len, k_len;
  int *b_idx, *sb_idx, sb_len_tot;
  // int **prof, **sb_grid, nblocks, sb_len_tot, *b_idx, *sb_idx;
  // cout << "finding new profile..." << endl;
  // pf_find_profile(prof, sb_grid, nblocks, mat_ep, dim);
  // cout << "done" << endl;

  mpc_t **sb_sys_mat, **sb_sys_inv_mat, *sb_sys_RHS, *sb_sys_sol;
  mpc_t mpc_tmp;
  mpc_init3(mpc_tmp, wp2, wp2);
  struct poly_frac pf_tmp;
  poly_frac_build(&pf_tmp);
  struct poly_frac **matsb;
  int num_it;
  if (print) {
    cout << "to_Fuchsian_global tol: ";
    mpfr_out_str(stdout, 10, 0, mpfr_tol, MPFR_RNDN); cout << endl;
  }
  int wp_bin = - mpfr_log2_int(mpfr_tol);
  fprintf(terminal, "off-diagonal: "); fflush(terminal); usleep(sleep_time);
  for (int b=0; b<nblocks; b++) {
    if (b > 0) {fprintf(terminal, "\033[18D\033[K");}// fflush(terminal); usleep(sleep_time);}
    fprintf(terminal, "block %3d /%3d... ", b, nblocks-1); fflush(terminal); usleep(sleep_time);
    b_len = prof[b][1] - prof[b][0] + 1;
    if (print) {
    cout << "################################################################# b = " << b << endl;
    cout << "b_len = " << b_len << endl;
    }
    // if (b == 40) dbg = 1;
    // matbR = new ex*[b_len];
    // rk2_subtens(matbR, mat_ep, prof[b][0], prof[b][1], prof[b][0]);
    // p_rank_debugR = Poinc_rank(matbR, b_len, b_len, g_et);
    // cout << "p_rank of matbR (b = " << b << ") = " << p_rank_debugR << endl;
    // getchar();

    // initialize to zero sub-diagoal row of the transfomation matrices
    // cout << "initialize to zero sub-diagoal row of the transfomation matrices" << endl;
    build_block_indices(b_idx, sb_idx, &b_len, &sb_len_tot, b, prof, sb_grid);
    for (int i=0; i<b_len; i++) {
      // for (int j=0; j<sb_len_tot; j++) {
        for (int j=0; j<prof[b][1]; j++) {
        // for (int j=0; j<dim; j++) {
        // cout << "i, j = " << i << ", " << j << endl;
        // cout << "b_idx, sb_idx = " << b_idx[i] << ", " << sb_idx[j] << endl;
        // poly_frac_set_zero(&tmat[b_idx[i]][sb_idx[j]]);
        // poly_frac_set_zero(&inv_tmat[b_idx[i]][sb_idx[j]]);
        poly_frac_set_zero(&tmat[b_idx[i]][j]);
        poly_frac_set_zero(&inv_tmat[b_idx[i]][j]);
      }
    }

    for (int sb=b-1; sb>=0; sb--) {
      sb_len = prof[sb][1] - prof[sb][0] + 1;
      big_dim = b_len*sb_len;
      if (print) {
      cout << "--------------- sb = " << sb << endl;
      cout << "sb_len = " << sb_len << endl;
      }
      // cout << "b_len*sb_len = " << big_dim << endl;
      // matbC = new ex*[sb_len];
      // rk2_subtens(matbC, mat_ep, prof[sb][0], prof[sb][1], prof[sb][0]);
      // p_rank_debugC = Poinc_rank(matbC, sb_len, sb_len, g_et);
      // cout << "p_rank of matbC (sb = " << sb << ") = " << p_rank_debugC << endl;
      // if (!sb_grid[b][sb]) {
      //   // cout << "skip" << endl;
      //   continue;
      // }
      
      matsb = new struct poly_frac*[b_len];
      poly_frac_rk2_subtens(&matsb, &mat_ep, prof[b][0], prof[b][1], prof[sb][0]);
      // if (b == 19 && (sb == 6 || sb == 0)) {
      // cout << "sub-diagonal block:" << endl;
      // poly_frac_rk2_print(matsb, b_len, sb_len);
      // }

      num_it = 0;
      while (1) {
        if (print || dbg) cout << "it = " << num_it << endl;
        if (print || dbg) cout << "matsb:" << endl;
        if (print || dbg) poly_frac_rk2_print(matsb, b_len, sb_len);
        // compute Poinc. Rank of the sub-block
        poly_frac_rk2_prune_rel_tol(matsb, wp_bin, b_len, sb_len);
        poly_frac_rk2_normal(matsb, b_len, sb_len);
        p_rank = poly_frac_Poinc_rank(matsb, b_len, sb_len);
        if (print || dbg) cout << "matsb pruned:" << endl;
        if (print || dbg) poly_frac_rk2_print(matsb, b_len, sb_len);
        if (print || dbg) cout << "p_rank = " << p_rank << endl;
        // getchar();
        if (p_rank <= 0) {
          break;
        }
        // poly_frac_rk2_print(matsb, b_len, sb_len);
        // fill sub-block system matrix and RHS
        // cout << "fill sub-block system matrix and RHS" << endl;
        malloc_rk2_tens(sb_sys_mat, big_dim, big_dim);
        // initialize to zero system matrix
        for (int i=0; i<big_dim; i++) {
          for (int j=0; j<big_dim; j++) {
            mpc_init3(sb_sys_mat[i][j], wp2, wp2);
            mpc_set_ui(sb_sys_mat[i][j], 0, MPFR_RNDN);
          }
        }
        malloc_rk2_tens(sb_sys_inv_mat, big_dim, big_dim);
        init_rk2_mpc(sb_sys_inv_mat, big_dim, big_dim);
        sb_sys_RHS = new mpc_t[big_dim];
        init_rk1_mpc(sb_sys_RHS, big_dim);
        sb_sys_sol = new mpc_t[big_dim];
        init_rk1_mpc(sb_sys_sol, big_dim);
        for (int iR=0; iR<b_len; iR++) {
          // fill RHS
          if (print || dbg) cout << "fill RHS" << endl;
          for (int jC=0; jC<sb_len; jC++) {
            // cout << "index = " << sb_len*iR+jC << endl;
            // poly_frac_print(&matsb[iR][jC]);
            rel_err_poly_frac_extract_LO(
              &mpc_realref(sb_sys_RHS[sb_len*iR+jC]),
              &mpc_imagref(sb_sys_RHS[sb_len*iR+jC]),
              &matsb[iR][jC],
              roots, p_rank+1
              , wp_bin
            );
            mpc_neg(sb_sys_RHS[sb_len*iR+jC], sb_sys_RHS[sb_len*iR+jC], MPFR_RNDN);
            if (print || dbg) {print_mpc(&sb_sys_RHS[sb_len*iR+jC]); cout << endl;}
          }
          for (int jR=0; jR<b_len; jR++) {
            for (int iC=0; iC<sb_len; iC++) {
              // contribution from diagonal_block[b] * sub-diagnoal block[b, sb]
              if (print || dbg) cout << "contribution from diagonal_block[b] * sub-diagnoal block[b, sb]" << endl;
              // poly_frac_print(&mat_ep[prof[b][0]+iR][prof[b][0]+jR]);
              rel_err_poly_frac_extract_LO(
                &mpc_realref(sb_sys_mat[sb_len*iR+iC][sb_len*jR+iC]),
                &mpc_imagref(sb_sys_mat[sb_len*iR+iC][sb_len*jR+iC]),
                &mat_ep[prof[b][0]+iR][prof[b][0]+jR],
                roots, 1
                , wp_bin
              );
              if(print) {print_mpc(&sb_sys_mat[sb_len*iR+iC][sb_len*jR+iC]); cout << endl;}
            }
            if (iR == jR) {
              for (int iC=0; iC<sb_len; iC++) {
                // contribution from p_rank * sub-diagonal block[b, sb]
                if (print || dbg) cout << "contribution from p_rank * sub-diagonal block[b, sb]" << endl;
                mpc_add_ui(
                  sb_sys_mat[sb_len*iR+iC][sb_len*iR+iC],
                  sb_sys_mat[sb_len*iR+iC][sb_len*iR+iC],
                  p_rank,
                  MPFR_RNDN
                );
                // if (
                //   mpfr_get_exp(mpc_imagref(sb_sys_mat[sb_len*iR+iC][sb_len*iR+iC])) < -140 && !mpfr_zero_p(mpc_imagref(sb_sys_mat[sb_len*iR+iC][sb_len*iR+iC])) ||\
                //   mpfr_get_exp(mpc_realref(sb_sys_mat[sb_len*iR+iC][sb_len*iR+iC])) < -140 && !mpfr_zero_p(mpc_realref(sb_sys_mat[sb_len*iR+iC][sb_len*iR+iC]))
                // ) {
                //   cout << "iR, jR, iC = " << iR << ", " << jR << ", " << iC << endl;
                //   cout << "mat:" << endl;
                //   print_mpc(&sb_sys_mat[sb_len*iR+iC][sb_len*iR+iC]); cout << endl;
                //   fflush(stdout);
                //   perror("constr mat construction: cancellation missed");
                //   exit(1);
                // }
                // print_mpc(&sb_sys_mat[sb_len*iR+iC][sb_len*iR+iC]);
                // cout << endl;
                for (int jC=0; jC<sb_len; jC++) {
                  // contribution from sub-diagnoal block[b, sb] * diagonal_block[sb]
                  // cout << "contribution from sub-diagnoal block[b, sb] * diagonal_block[sb]" << endl;
                  // cout << "pf:" << endl;
                  // poly_frac_print(&mat_ep[prof[sb][0]+jC][prof[sb][0]+iC]);
                  rel_err_poly_frac_extract_LO(
                    &mpc_realref(mpc_tmp),
                    &mpc_imagref(mpc_tmp),
                    &mat_ep[prof[sb][0]+jC][prof[sb][0]+iC],
                    roots, 1
                    , wp_bin
                  );
                  if (print || dbg) {cout << "L0 = "; print_mpc(&mpc_tmp); cout << endl;}
                  // cout << "original matrix element:" << endl;
                  // print_mpc(&sb_sys_mat[sb_len*iR+iC][sb_len*iR+jC]);
                  // cout << endl;
                  mpc_neg(mpc_tmp, mpc_tmp, MPFR_RNDN);
                  exp_rel_err_mpc_add(
                    sb_sys_mat[sb_len*iR+iC][sb_len*iR+jC],
                    sb_sys_mat[sb_len*iR+iC][sb_len*iR+jC],
                    mpc_tmp,
                    MPFR_RNDN,
                    wp_bin
                  );
                  if (print || dbg) {cout << "output:" << endl; print_mpc(&sb_sys_mat[sb_len*iR+iC][sb_len*iR+jC]); cout << endl;}
                }
              }
            }
          }
        }

        // solve linear system
        if (print || dbg) {
        cout << "solve linear system" << endl;
        cout << "check RHS:" << endl;
        print_poly(sb_sys_RHS, big_dim-1);
        cout << "check matrix:" << endl;
        print_rk2_mpc(sb_sys_mat, big_dim, big_dim);
        }
        mp_inverse(sb_sys_inv_mat, sb_sys_mat, big_dim);
        if (print || dbg) {
        cout << "check inverse matrix:" << endl;
        print_rk2_mpc(sb_sys_inv_mat, big_dim, big_dim);
        }
        // cout << "inverse computed" << endl;
        for (int i=0; i<big_dim; i++) {
          mpc_set_ui(sb_sys_sol[i], 0, MPFR_RNDN);
          int j;
          for (j=0; j<big_dim; j++) {
            // mpc_fma(
            //   sb_sys_sol[i],
            //   sb_sys_inv_mat[i][j], sb_sys_RHS[j], sb_sys_sol[i],
            //   MPFR_RNDN
            // );
            exp_rel_err_mpc_mul(mpc_tmp, sb_sys_inv_mat[i][j], sb_sys_RHS[j], MPFR_RNDN, wp_bin);
            exp_rel_err_mpc_add(sb_sys_sol[i], sb_sys_sol[i], mpc_tmp, MPFR_RNDN, wp_bin);
          }
          // // rel err prune on last contribution
          // mpc_mul(mpc_tmp, sb_sys_inv_mat[i][j], sb_sys_RHS[j], MPFR_RNDN);
          // rel_err_mpc_add(sb_sys_sol[i], sb_sys_sol[i], mpc_tmp, MPFR_RNDN);
        }
        if (print || dbg) {
        cout << "solution:" << endl;
        print_poly(sb_sys_sol, big_dim-1);
        }
        // getchar();

        // apply transformation
        // cout << "apply transformation" << endl;
        pf_apply_sub_diagonal_p_reduction(
          mat_ep,
          p_rank, b, sb, nblocks, prof,
          sb_sys_sol, roots, nroots
        );

        // cumulate transformation
        // cout << "cumulate transformation" << endl;
        for (int i=0; i<b_len; i++) {
          for (int j=0; j<sb_len; j++) {
            poly_frac_set_mpc(&pf_tmp, &sb_sys_sol[i*sb_len + j], nroots);
            poly_frac_mul_sym_pow(&pf_tmp, &pf_tmp, -p_rank);
            rel_err_poly_frac_add_pf(
              &tmat[prof[b][0]+i][prof[sb][0]+j],
              &tmat[prof[b][0]+i][prof[sb][0]+j],
              &pf_tmp,
              roots,
              NULL,
              wp_bin
            );
          }
        }

        // cumulate inverse transformation
        // cout << "cumulate inverse transformation" << endl;
        // (modify only the elements of the b-th row up to the sb-th element)
        for (int k=0; k<sb; k++) {
          // cout << "k = " << k << endl;
          k_len = prof[k][1] - prof[k][0] + 1;
          if (!sb_grid[sb][k]) {
            // cout << "skip" << endl;
            continue;
          }
          // update sub-diagonal block (b, k)
          for (int i=0; i<b_len; i++) {
            for (int j=0; j<k_len; j++) {
              // cout << "i, j = " << i << ", " << j << endl;
              for (int a=0; a<sb_len; a++) {
                // cout << "a = " << a << endl;
                // cout << "it, jt = " << prof[sb][0] + a << ", " << prof[k][0] + j << endl;
                // poly_frac_print(&inv_tmat[prof[sb][0] + a][prof[k][0] + j]);
                poly_frac_mul_sym_pow(&pf_tmp, &inv_tmat[prof[sb][0] + a][prof[k][0] + j], -p_rank);
                poly_frac_neg(&pf_tmp);
                rel_err_poly_frac_add_pf(
                  &inv_tmat[prof[b][0] + i][prof[k][0] + j],
                  &inv_tmat[prof[b][0] + i][prof[k][0] + j],
                  &pf_tmp,
                  roots,
                  &sb_sys_sol[i*sb_len + a],
                  wp_bin
                );
              }
            }
          }
        }
        for (int i=0; i<b_len; i++) {
          for (int j=0; j<sb_len; j++) {
            poly_frac_set_mpc(&pf_tmp, &sb_sys_sol[i*sb_len + j], nroots);
            poly_frac_mul_sym_pow(&pf_tmp, &pf_tmp, -p_rank);
            poly_frac_neg(&pf_tmp);
            rel_err_poly_frac_add_pf(
              &inv_tmat[prof[b][0] + i][prof[sb][0] + j],
              &inv_tmat[prof[b][0] + i][prof[sb][0] + j],
              &pf_tmp,
              roots,
              NULL,
              wp_bin
            );
          }
        }

        // FREE
        mpc_rk2_clear(sb_sys_mat, big_dim, big_dim);
        del_rk2_tens(sb_sys_mat, big_dim);
        mpc_rk2_clear(sb_sys_inv_mat, big_dim, big_dim);
        del_rk2_tens(sb_sys_inv_mat, big_dim);
        mpc_rk1_clear(sb_sys_RHS, big_dim);
        mpc_rk1_clear(sb_sys_sol, big_dim);
        delete[] sb_sys_sol;

        num_it++;
        // if (num_it > 10000) {
        if (num_it > 10) {
          perror("too many iterations in global fuchsianization");
          exit(1);
        }
      }

      if (num_it > 0) {
        del_rk2_tens(sb_sys_mat, big_dim);
        del_rk2_tens(sb_sys_inv_mat, big_dim);
      }
      delete[] matsb;
    }
  }
  fprintf(terminal, "\033[18D\033[K"); fflush(terminal); usleep(sleep_time);
  fprintf(terminal, "\033[14D\033[K"); fflush(terminal); usleep(sleep_time);
  dbg = 0;

  // FREE
  mpc_clear(mpc_tmp);
}


void pf_NormalizeMat(
  // OUTPUT
  struct poly_frac **tmat, struct poly_frac **inv_tmat,
  int *num_classes, int *eq_class, mpc_t **eig_list, int *eig_grid, 
  // IN-OUT
  struct poly_frac **mat_ep,
  // INPUT
  int dim, int nblocks, int **prof, int **sb_grid,
  mpc_t *roots, int nroots,
  FILE *terminal
) {
  int print = 0;
  fprintf(terminal, "normalize: "); usleep(sleep_time);

  // if (print) cout << "compute avarage exp2" << endl;
  // int p_rank = poly_frac_Poinc_rank(mat_ep, dim, dim);
  // int num_deg_max = poly_frac_k2_num_deg_max(mat_ep, dim, dim);
  // double *num_exp2_avg = new double[num_deg_max+p_rank+2];
  // poly_frac_k2_num_exp2_avg(num_exp2_avg, mat_ep, dim, dim, p_rank, num_deg_max);
  // if (print) {
  //   for (int k=-p_rank-1; k<=num_deg_max; k++) {
  //     cout << "pow " << k << ": " << num_exp2_avg[k] << endl;
  //   }
  //   cout << endl;
  // }

  if (print) cout << "enter NormalizeDiagonal" << endl;
  pf_NormalizeDiagonal(
    tmat, inv_tmat,
    eig_list, num_classes, eq_class, eig_grid,
    mat_ep,
    dim, nblocks, prof, sb_grid,
    roots, nroots,
    terminal
  );
  if (print) cout << "exit NormalizeDiagonal" << endl;
  // int p_rank = poly_frac_Poinc_rank(mat_ep, dim, dim);
  // cout << "Poinc. rank = " << p_rank << endl;
  // for (int b_len, b=0; b<nblocks; b++) {
  //   b_len = prof[b][1] - prof[b][0] + 1;
  //   poly_frac_rk2_slice(&mat_ep, &mat_ep, prof[b][0], prof[b][1], prof[b][0]);
  //   p_rank = poly_frac_Poinc_rank(mat_ep, b_len, b_len);
  //   if (p_rank > 0) {
  //     cout << "ERROR: Poinc. rank = " << p_rank;
  //     cout << " for block b = " << b << endl;
  //     exit(0);
  //   }
  //   poly_frac_rk2_unslice(&mat_ep, prof[b][0], prof[b][1], prof[b][0]);
  // }
  // getchar();
  // exit(0);
  // cout << "inv tmat after NormalizeDiagonal" << endl;
  // poly_frac_rk2_print(inv_tmat, dim, dim);

  // find new profile
  // pf_find_profile(prof, sb_grid, nblocks, mat_ep, dim);
  struct poly_frac **tmatsb, **inv_tmatsb;
  malloc_rk2_tens(tmatsb, dim, dim);
  poly_frac_rk2_build(tmatsb, dim, dim);
  malloc_rk2_tens(inv_tmatsb, dim, dim);
  poly_frac_rk2_build(inv_tmatsb, dim, dim);
  if (print) cout << "enter to_Fuchsian_global" << endl;
  pf_to_Fuchsian_global(
    mat_ep,
    tmatsb, inv_tmatsb,
    dim, nblocks, prof, sb_grid,
    roots, nroots,
    terminal
  );
  if (print) cout << "exit to_Fuchian_global" << endl;
  // p_rank = poly_frac_Poinc_rank(mat_ep, dim, dim);
  // cout << "Poinc. rank = " << p_rank << endl;
  // for (int b_len, b=0; b<nblocks; b++) {
  //   b_len = prof[b][1] - prof[b][0] + 1;
  //   for (int sb_len, sb=b; sb>=0; sb--) {
  //     sb_len = prof[sb][1] - prof[sb][0] + 1;
  //     poly_frac_rk2_slice(&mat_ep, &mat_ep, prof[b][0], prof[b][1], prof[sb][0]);
  //     p_rank = poly_frac_Poinc_rank(mat_ep, b_len, sb_len);
  //     if (p_rank > 0) {
  //       cout << "ERROR: Poinc. rank = " << p_rank;
  //       cout << " for block (b, sb) = (" << b << ", " << sb << ")" << endl;
  //       cout << "i, j = " << prof[b][0] << ", " << prof[sb][0] << endl;
  //       cout << "address = " << &mat_ep[0][0] << endl;
  //       cout << "block:" << endl;
  //       poly_frac_rk2_print(mat_ep, b_len, sb_len);
  //       exit(0);
  //     }
  //     poly_frac_rk2_unslice(&mat_ep, prof[b][0], prof[b][1], prof[sb][0]);
  //   }
  // }
  // getchar();
  // exit(0);
  // cout << "MATRIX:" << endl;
  // poly_frac_rk2_print(pfmat, dim, dim);
  // cout << "TMAT SB:" << endl;
  // poly_frac_rk2_print(tmatsb, dim, dim);
  // cout << "INV TMAT SB:" << endl;
  // poly_frac_rk2_print(inv_tmatsb, dim, dim);

  // cumulate transformations
  struct poly_frac pf_tmp;
  poly_frac_build(&pf_tmp);
  int b_len, sb_len;
  for (int b=0; b<nblocks; b++) {
    b_len = prof[b][1] - prof[b][0] + 1;
    // cout << "-------------------------" << endl;
    // cout << "b = " << b << endl;
    // cout << "b_len = " << b_len << endl;
    // cout << endl;
    for (int sb=b-1; sb>=0; sb--) {
      sb_len = prof[sb][1] - prof[sb][0] + 1;
      // cout << "-------------------------" << endl;
      // cout << "sb = " << sb << endl;
      // cout << "sb_len = " << sb_len << endl;
      // cout << endl;
      for (int i=0; i<b_len; i++) {
        for (int j=0; j<sb_len; j++) {
          // cout << "i, j = " << i << ", " << j << endl;
          // direct transformation
          // cout << "direct transformation" << endl;
          // cout << "1st factor:" << endl;
          // poly_frac_print(&tmat[prof[b][0]+i][prof[b][0]]);
          // cout << "2nd factor:" << endl;
          // poly_frac_print(&tmatsb[prof[b][0]][prof[sb][0]+j]);
          poly_frac_mul_pf(
            &tmat[prof[b][0]+i][prof[sb][0]+j],
            &tmat[prof[b][0]+i][prof[b][0]],
            &tmatsb[prof[b][0]][prof[sb][0]+j]
          );
          // cout << "result:" << endl;
          // poly_frac_print(&tmat[prof[b][0]+i][prof[sb][0]+j]);
          for (int k=1; k<b_len; k++) {
            // cout << "k = " << k << endl;
            // cout << "1st factor:" << endl;
            // poly_frac_print(&tmat[prof[b][0]+i][prof[b][0]+k]);
            // cout << "2nd factor:" << endl;
            // poly_frac_print(&tmatsb[prof[b][0]+k][prof[sb][0]+j]);
            poly_frac_mul_pf(
              &pf_tmp,
              &tmat[prof[b][0]+i][prof[b][0]+k],
              &tmatsb[prof[b][0]+k][prof[sb][0]+j]
            );
            // cout << "2nd addend:" << endl;
            // poly_frac_print(&pf_tmp);
            // cout << "1st addend:" << endl;
            // poly_frac_print(&tmat[prof[b][0]+i][prof[sb][0]+j]);
            poly_frac_add_pf(
              &tmat[prof[b][0]+i][prof[sb][0]+j],
              &tmat[prof[b][0]+i][prof[sb][0]+j],
              &pf_tmp,
              roots,
              NULL
            );
            // cout << "result:" << endl;
            // poly_frac_print(&tmat[prof[b][0]+i][prof[sb][0]+j]);
          }
          // inverse transformation
          // cout << "inverse transformation" << endl;
          // cout << "1st factor:" << endl;
          // poly_frac_print(&inv_tmatsb[prof[b][0]+i][prof[sb][0]]);
          // cout << "2nd factor:" << endl;
          // poly_frac_print(&inv_tmat[prof[sb][0]][prof[sb][0]+j]);
          poly_frac_mul_pf(
            &inv_tmat[prof[b][0]+i][prof[sb][0]+j],
            &inv_tmatsb[prof[b][0]+i][prof[sb][0]],
            &inv_tmat[prof[sb][0]][prof[sb][0]+j]
          );
          // cout << "result:" << endl;
          // poly_frac_print(&inv_tmat[prof[b][0]+i][prof[sb][0]+j]);
          for (int k=1; k<sb_len; k++) {
            // cout << "k = " << k << endl;
            // cout << "1st factor:" << endl;
            // poly_frac_print(&inv_tmatsb[prof[b][0]+i][prof[sb][0]+k]);
            // cout << "2nd factor:" << endl;
            // poly_frac_print(&inv_tmat[prof[sb][0]+k][prof[sb][0]+j]);
            poly_frac_mul_pf(
              &pf_tmp,
              &inv_tmatsb[prof[b][0]+i][prof[sb][0]+k],
              &inv_tmat[prof[sb][0]+k][prof[sb][0]+j]
            );
            // cout << "2nd addend:" << endl;
            // poly_frac_print(&pf_tmp);
            // cout << "1st addend:" << endl;
            // poly_frac_print(&inv_tmat[prof[b][0]+i][prof[sb][0]+j]);
            poly_frac_add_pf(
              &inv_tmat[prof[b][0]+i][prof[sb][0]+j],
              &inv_tmat[prof[b][0]+i][prof[sb][0]+j],
              &pf_tmp,
              roots,
              NULL
            );
            // cout << "result:" << endl;
            // poly_frac_print(&inv_tmat[prof[b][0]+i][prof[sb][0]+j]);
          }
        }
      }
    }
  }
  // cout << "FINAL TMAT:" << endl;
  // poly_frac_rk2_print(tmat, dim, dim);
  // cout << "FINAL INV TMAT:" << endl;
  // poly_frac_rk2_print(inv_tmat, dim, dim);

  // FREE
  poly_frac_rk2_free(tmatsb, dim, dim);
  del_rk2_tens(tmatsb, dim);
  poly_frac_rk2_free(inv_tmatsb, dim, dim);
  del_rk2_tens(inv_tmatsb, dim);
  poly_frac_free(&pf_tmp);

  fprintf(terminal, "\033[11D\033[K"); fflush(terminal);
  usleep(sleep_time);
}


void pf_check_NormalizeMat(
  struct poly_frac **tmat, struct poly_frac **inv_tmat,
  int num_classes, int *eq_class, mpc_t *eig_list, int *eig_grid,
  struct poly_frac **mat_ep,
  int dim, int nblocks, int **prof, int **sb_grid,
  mpc_t *roots, int nroots,
  FILE *terminal
) {
  int p_rank = poly_frac_Poinc_rank(mat_ep, dim, dim);
  cout << "input Poinc. rank = " << p_rank << endl;

  struct poly_frac **pfmat_tmp;
  malloc_rk2_tens(pfmat_tmp, dim, dim);
  poly_frac_rk2_build(pfmat_tmp, dim, dim);
  struct poly_frac **pfmat_tmp1;
  malloc_rk2_tens(pfmat_tmp1, dim, dim);
  poly_frac_rk2_build(pfmat_tmp1, dim, dim);
  struct poly_frac **mat_ep_backup;
  malloc_rk2_tens(mat_ep_backup, dim, dim);
  poly_frac_rk2_build(mat_ep_backup, dim, dim);
  poly_frac_rk2_set_pf_rk2(mat_ep_backup, mat_ep, dim, dim);
  // cout << "mat = " << endl;
  // poly_frac_rk2_print_to_math(mat_ep, dim, dim, roots);
  // cout << "start NormalizeMat" << endl;
  pf_NormalizeMat(
    tmat, inv_tmat,
    &num_classes, eq_class, &eig_list, eig_grid,
    mat_ep,
    dim, nblocks, prof, sb_grid,
    roots, nroots,
    terminal
  );
  cout << "end NormalizeMat" << endl;
  p_rank = poly_frac_Poinc_rank(mat_ep, dim, dim);
  cout << "Poinc. rank = " << p_rank << endl;
  // cout << "nmat = " << endl;
  // poly_frac_rk2_print_to_math(mat_ep, dim, dim, roots);
  // cout << "T = " << endl;
  // poly_frac_rk2_print_to_math(tmat, dim, dim, roots);
  // cout << "invT = " << endl;
  // poly_frac_rk2_print_to_math(inv_tmat, dim, dim, roots);
  // exit(0);

  // CHECK RESULTS
  // print_rk2_tens_to_mat(mat_ep, dim, dim);
  // print_rk2_tens_to_mat(inv_tmat, dim, dim);
  // getchar();
  cout << "apply inverse transformation to whole matrix" << endl;
  // cout << "multiply on the right" << endl;
  poly_frac_rk2_mul_pf_rk2(
    mat_ep, mat_ep, inv_tmat,
    dim, dim, dim,
    roots
  );
  // poly_frac_rk2_print(mat_ep, dim, dim);
  // cout << "multiply on the left" << endl;
  poly_frac_rk2_mul_pf_rk2(
    mat_ep, tmat, mat_ep,
    dim, dim, dim,
    roots
  );
  // poly_frac_rk2_print(mat_ep, dim, dim);
  // cout << "derivative contribution:" << endl;
  for (int i=0; i<dim; i++) {
    for (int j=0; j<dim; j++) {
      // cout << "i, j = " << i << ", " << j << endl;
      // cout << "input:" << endl;
      // poly_frac_print(&inv_tmat[i][j]);
      poly_frac_derivative(&pfmat_tmp[i][j], &inv_tmat[i][j], roots);
      // cout << "derivative:" << endl;
      // poly_frac_print(&pfmat_tmp[i][j]);
      poly_frac_neg(&pfmat_tmp[i][j]);
    }
  }
  // cout << "derInvT = " << endl;
  // poly_frac_rk2_print_to_math(pfmat_tmp, dim, dim, roots);
  // exit(0);
  // poly_frac_rk2_print(pfmat_tmp, dim, dim);
  poly_frac_rk2_mul_pf_rk2(
    pfmat_tmp1, tmat, pfmat_tmp,
    dim, dim, dim,
    roots
  );
  for (int i=0; i<dim; i++) {
    for (int j=0; j<dim; j++) {
      poly_frac_add_pf(
        &mat_ep[i][j], &mat_ep[i][j], &pfmat_tmp1[i][j],
        roots, NULL
      );
    }
  }
  cout << "applied" << endl;
  cout << "check difference" << endl;
  for (int i=0; i<dim; i++) {
    for (int j=0; j<dim; j++) {
      poly_frac_set_pf(&pfmat_tmp[i][j], &mat_ep[i][j]);
      poly_frac_neg(&pfmat_tmp[i][j]);
      poly_frac_add_pf(&pfmat_tmp[i][j], &pfmat_tmp[i][j], &mat_ep_backup[i][j], roots, NULL);
      poly_frac_prune(&pfmat_tmp[i][j]);
      if (pfmat_tmp[i][j].num_vdeg != -1) {
        cout << "#### (i, j) = " << i << ", " << j << endl;
        cout << "original:" << endl;
        poly_frac_print(&mat_ep_backup[i][j]);
        cout << "reconstructed:" << endl;
        poly_frac_print(&mat_ep[i][j]);
        cout << "difference:" << endl;
        poly_frac_print(&pfmat_tmp[i][j]);
      }
    }
  }
  cout << "done" << endl;
}

