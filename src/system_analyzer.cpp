#include <iostream>
#include <cstring>
// #include <complex>
#include "mpc.h"

using namespace std;

#include "global_vars.h"
#include "malloc_defs.h"
#include "utils.h"
// #include "conversions.h"
#include "eq_solver.h"
#include "tensor_utils.h"
#include "poly_frac.h"
#include "setup.h"
#include "codec.h"
extern "C" {
  #include "cpoly.h"
  #include "mp_roots_poly.h"
  #include "in_out.h"
  #include "ex_tree.h"
}


void generate_sub_diag_grid(
  // OUTPUT
  int ***subblock,
  // INPUT
  int **profile, int nblocks, struct poly_frac **pfmat
) {
  int i, j;

  // grid to store the positions of non-vanishing sub-blocks
  (*subblock) = new int * [nblocks];
  for (i = 0; i < nblocks; i++) {
    (*subblock)[i] = new int [nblocks];
  }
  // initialize to zero every element
  for (i = 0; i < nblocks; i++) {
    for (j = 0; j < nblocks; j++) { 
      (*subblock)[i][j] = 0;
    }
  }

  // find non-vanishing sub-blocks
  int found;
  for (i = 0; i < nblocks; i++) {
    for (j = 0; j < i; j++) {
      found = 0;
      // select sub-block and look for non-vanishing elements
      for (int kk = profile[i][0]; kk <= profile[i][1]; kk++) {
        for (int ll = profile[j][0]; ll <= profile[j][1]; ll++) {
          if(pfmat[kk][ll].num_vdeg != -1) {
            found=1;
            break;
          }
        }
        if(found == 1) {
          (*subblock)[i][j] = 1;
          break;
        }
      }
    }
  }
}


// #PAIR
void poly_frac_rk2_find_profile(
  // OUTPUT
  int ***profile, int ***subblock, int *nblocks,
  // INPUT
  struct poly_frac **pfmat, int g_dim
) {
  int i, j;

  // temporary profile that stores the position
  // of the last non-zero element of each row
  int *tmpprof;
  tmpprof = new int [g_dim];
  for (i = 0; i < g_dim; i++) {
    tmpprof[i] = i;
    for (j = g_dim-1; j >= i; j--) {
      if(pfmat[i][j].num_vdeg != -1) {
        tmpprof[i] = j;
        break;
      }
    }
  }

  // extract number of blocks and compute block ranges
  *nblocks = 0;
  int **inner_profile = new int*[g_dim];
  for (i=0; i<g_dim; i++) {
    inner_profile[i] = new int[2];
  }
  int col_idx_max = 0;
  inner_profile[0][0] = 0;
  for (i = 0; i < g_dim; i++) {
    if (tmpprof[i] > col_idx_max) {
      // update max column index for the block
      col_idx_max = tmpprof[i];
    }

    if (i == col_idx_max) {
      // end of block
      inner_profile[(*nblocks)++][1] = col_idx_max;
      if (i == g_dim-1) {
        break;
      }
      // start profile for next block
      inner_profile[*nblocks][0] = i+1;
    }
  }

  // *nblocks X 2 matrix to store block ranges in the output
  (*profile) = new int * [*nblocks];
  for (i = 0; i < *nblocks; i++) {
    (*profile)[i] = new int [2];
    (*profile)[i][0] = inner_profile[i][0];
    (*profile)[i][1] = inner_profile[i][1];
  }

  // // compute block ranges from the temporary profile
  // profile[0][0] = 0;
  // profile[0][1] = tmpprof[0];
  // int k = 1;
  // for (i = 1; i < g_dim; i++) {
  //   if(tmpprof[i] == i) {
  //     profile[k][0] = profile[k-1][1] + 1;
  //     profile[k][1] = i;
  //     k++;
  //   }
  // }
  delete[] tmpprof;

  // generte grid to store the positions of non-vanishing sub-blocks
  generate_sub_diag_grid(
    subblock,
    *profile, *nblocks, pfmat
  );
  
}


void update_new_roots(
  // OUTPUT
  int *label, int *num_out_roots, mpc_t *out_roots, mpfr_t *out_tols,
  // INPUT
  int *num_in_roots, mpc_t *in_roots, mpfr_t *in_tols,
  mpc_t *root, mpfr_t *tol
) {
  mpc_t tmpc;
  mpc_init3(tmpc, wp2, wp2);
  mpfr_t tmpfr, tmpfr1;
  mpfr_init2(tmpfr, wp2);
  mpfr_init2(tmpfr1, wp2);
  mpc_t *old_root;
  mpfr_t *old_tol;
  int new_root = 1;
  for (int kp=0; kp<*num_in_roots+*num_out_roots; kp++) {
    // cout << "compare with kp = " << kp << endl;
    if (kp <*num_in_roots) {
      old_root = in_roots + kp;
      old_tol = in_tols + kp;
    } else {
      old_root = out_roots + kp - *num_in_roots;
      old_tol = out_tols + kp - *num_in_roots;
    }
    mpc_sub(tmpc, *root, *old_root, MPFR_RNDN);
    if (!mpc_lessthan_tol(*old_root)) {
      // if (dbg) {
      //   cout << "dividing" << endl;
      //   cout << "root:" << endl;
      //   print_mpc(root); cout << endl;
      //   cout << "old_root:" << endl;
      //   print_mpc(old_root); cout << endl;
      //   cout << "tmpc:" << endl;
      //   print_mpc(&tmpc); cout << endl;
      // }
      mpc_div(tmpc, tmpc, *old_root, MPFR_RNDN);
    }
    mpc_abs(tmpfr, tmpc, MPFR_RNDN);
    mpfr_max(tmpfr1, *old_tol, *tol, MPFR_RNDN);
    // cout << "old_root:" << endl;
    // print_mpc(old_root); cout << endl;
    // cout << "old_tol:" << endl;
    // mpfr_out_str(stdout, 10, 0, *old_tol, MPFR_RNDN); cout << endl;
    // if (kp == 2) {
    // cout << "rel difference:" << endl;
    // mpfr_out_str(stdout, 10, 0, tmpfr, MPFR_RNDN); cout << endl;
    // cout << "new tol:" << endl;
    // mpfr_out_str(stdout, 10, 0, *tol, MPFR_RNDN); cout << endl;
    // cout << "old tol:" << endl;
    // mpfr_out_str(stdout, 10, 0, *old_tol, MPFR_RNDN); cout << endl;
    // }
    if (mpfr_lessequal_p(tmpfr, tmpfr1)) {
      // cout << "it corresponds to root n. " << kp << endl;
      *label = kp;
      if (mpfr_less_p(*tol, *old_tol)) {
        mpfr_set(*old_tol, *tol, MPFR_RNDN);
        mpc_set(*old_root, *root, MPFR_RNDN);
      }
      new_root = 0;
      break;
    }
  }
  if (new_root) {
    // cout << "found new root:" << endl;
    // print_mpc(root); cout << endl;
    mpc_init3(out_roots[*num_out_roots], wp2, wp2);
    // cout << "init new root" << endl;
    mpc_set(out_roots[*num_out_roots], *root, MPFR_RNDN);
    // cout << "set new root" << endl;
    mpfr_init2(out_tols[*num_out_roots], wp2);
    mpfr_set(out_tols[*num_out_roots], *tol, MPFR_RNDN);
    // cout << "set new tol" << endl;
    *label = *num_in_roots + *num_out_roots;
    // cout << "set new label" << endl;
    (*num_out_roots)++;
  }
  // getchar();

}


// <<linearize_masses_nd>>
void linearize_masses_nd(
	// IN-OUT
	struct lnode *nd,
	// INPUT
	int *is_mass
) {
	// change m2 -> m^2
	if (nd->op == 's') {
    if (nd->number == 0) {
      return;
    }
    if (is_mass[nd->number-1]) {
      // cout << "found mass symbol" << endl;
      // change m2 -> m^2 at node level
      nd->op = '^';
      nd->son = (struct lnode*) malloc(sizeof(struct lnode));
      create_lnode(nd->son, 0, 's');
      nd->son->bro = NULL;
      nd->son->number = nd->number;
      nd->number = 2;
    }
	} else {
		struct lnode *son = nd->son;
		// printf("\nn operand = %d:\n", root->n);
		while (son) {
			linearize_masses_nd(son, is_mass);
			son = son->bro;
		}
	}
}


// <<process_branch_points>>
void process_branch_points(
  // OUTPUT
  int *num_branch_roots, mpc_t *branch_roots, mpfr_t *branch_tols,
  int **branch_sing_lab,
  int *branch_deg, mpc_t **branch_poly, 
  // IN-OUT
  int *nroots, mpc_t **roots, mpfr_t **tols,
  // INPUT
  int nbranches, char **branch_cut_str,
  int ninvs, char **symbols, int *is_mass,
  int *skip_inv, char ***ep_kin, struct poly_frac *pspf
) {

  // LINEARIZE MASS SYMBOLS
  // branch cuts input expressions depends on masses
  // instead of squared masses
  char **lin_symbols;
  linearize_masses_str(&lin_symbols, ninvs, symbols, is_mass);
  // cout << "linearized symbols:" << endl;
  // for (int s=0; s<=ninvs; s++) {
  //   cout << lin_symbols[s] << endl;
  // }

  // ANALYSE BRANCH CUTS
  struct lnode nd; lnode_build(&nd);
  char *nd_str;
  struct poly_frac pf;
  poly_frac_build(&pf);

  // int dummy_nroots = 1;
	// mpc_t *roots = new mpc_t[1];
	// mpc_init3(roots[0], wp2, wp2);
	// mpc_set_ui(roots[0], 0, MPFR_RNDN);
	// mpfr_t *tols = new mpfr_t[1];
	// mpfr_init2(tols[0], wp2);
	// mpfr_set(tols[0], mpfr_tol, MPFR_RNDN);

  int ldeg;
  mpc_t *single_branch_roots;
  mpfr_t *single_branch_tols;
  *num_branch_roots = 0;
  for (int b=0; b<nbranches; b++) {
    // cout << "b = " << b << endl;
    // cout << "branch cut:" << branch_cut_str[b] << endl;

    // parse expression
    lnode_free(&nd);
    lnode_parse_expression(&nd, branch_cut_str[b], (const char**) lin_symbols, ninvs+1, 0);
    lnode_expand(&nd);
    // lnode_print(&nd, (const char**) lin_symbols); printf("\n");
    nd_str = lnode_to_str(&nd, (char*)", ");
    // printf("%s\n", nd_str);

    // decode expression
    decode_tree_pf(
      &pf,
      nroots, roots, tols,
      pspf,
      ep_kin, skip_inv-1, is_mass-1,
      &nd_str, (char*) ", ",
      wp2, mpfr_tol, 0
    );
    // cout << "branch cut poly frac:" << endl;
    // poly_frac_print(&pf);

    // branch cut polynomial
    branch_deg[b] = pf.num_deg;
    branch_poly[b] = new mpc_t[pf.num_deg+1];
    ldeg = pf.num_deg-pf.num_vdeg;
    for (int k=0; k<ldeg; k++) {
      mpc_init3(branch_poly[b][k], wp2, wp2);
      mpc_set_ui(branch_poly[b][k], 0, MPFR_RNDN);
    }
    for (int k=ldeg; k<=branch_deg[b]; k++) {
      mpc_init3(branch_poly[b][k], wp2, wp2);
      mpc_set(branch_poly[b][k], pf.coeffs[k-ldeg], MPFR_RNDN);
    }

    // cout << "coeffs:" << endl;
    // print_poly(branch_poly[b], branch_deg[b]);

    // find branch cut roots
    single_branch_roots = new mpc_t[branch_deg[b]+1];
    init_rk1_mpc(single_branch_roots, branch_deg[b]+1);
    single_branch_tols = new mpfr_t[branch_deg[b]+1];
    init_rk1_mpfr(single_branch_tols, branch_deg[b]+1);
    branch_sing_lab[b] = new int[branch_deg[b]+1];
    // cout << "finding roots..." << endl;
    // mp_zroots_impr(branch_poly[b], branch_deg[b], single_branch_roots, 1, wp2, mpfr_tol, single_branch_tols, 1);
    mp_zroots_impr(
      single_branch_roots, single_branch_tols,
      branch_poly[b], branch_deg[b], NULL,
      1, wp2, mpfr_tol, 1
    );
    for (int k=1; k<=branch_deg[b]; k++) {
      // cout << "root:"; print_mpc(&single_branch_roots[k]); cout << endl;
      update_new_roots(
        &branch_sing_lab[b][k], num_branch_roots, branch_roots, branch_tols,
        nroots, *roots, *tols,
        &single_branch_roots[k], &single_branch_tols[k]
      );
      // cout << "branch sing lab = " << branch_sing_lab[b][k] << endl;
      // cout << "num branch_roots = " << *num_branch_roots << endl;
    }
    delete[] single_branch_roots;
    delete[] single_branch_tols;
  }

  mpc_rk1_append(
    nroots, roots,
    *num_branch_roots, branch_roots
  );
  mpfr_rk1_append(
    nroots, tols,
    *num_branch_roots, branch_tols
  );
  *nroots += *num_branch_roots;

  // FREE
  lnode_free(&nd);

}


// <<linearize_mass_DE>>
void linearize_mass_DE(
  // IN-OUT
  struct lnode *nd,
  // INPUT
  int sym_idx
) {
  // multiply mass matrices times 2*m: A_m = 2*m*A_m2

  if (nd->op == '-' && nd->son->op == '*') {
    // if we have the negative of a product,
    // then just linearize the inner product
    linearize_mass_DE(nd->son, sym_idx);
    return;
  }

  // first check if m^2 and/or 2 are present in the denominator,
  // in this case remove them;
  // only works for expanded products
  int mul_m = 1, mul_2 = 1; // 0: exit, 1: multiply m, 2: multiply 2, 3: multiply both
  
  if (nd->op == '*') {
    struct lnode *son = nd->son, *tmp_nd;
    while (son && (mul_m  || mul_2)) {
      if (son->op == '/') {
        if (son->son->op == 's' && son->son->number == sym_idx) {
          mul_m = 0;
          // remove node
          lnode_remove_son(nd, son);
        } else if (son->son->op == 'n' && mpz_cmp_ui(son->son->mpz, 2) == 0) {
          mul_2 = 0;
          // remove node
          tmp_nd = son;
          son = son->bro;
          lnode_remove_son(nd, tmp_nd);
          continue;
        }
      } else if (son->op == '^') {
        if (son->son->op == 's' && son->son->number == sym_idx) {
          mul_m = 0;
          if (son->number == -1) {
            // remove node
            tmp_nd = son;
            son = son->bro;
            lnode_remove_son(nd, tmp_nd);
            continue;
          } else {
            son->number += 1;
          }
        } else if (son->son->op == 'n' && mpz_cmp_ui(son->son->mpz, 2) == 0) {
          mul_2 = 0;
          if (son->number == -1) {
            // remove node
            tmp_nd = son;
            son = son->bro;
            lnode_remove_son(nd, tmp_nd);
          } else {
            son->number += 1;
          }
        }
      }
      son = son->bro;
    }
  }

  if (mul_m == 0 && mul_2 == 0) {
    return;
  } else if (mul_m == 1 && mul_2 == 0) {
    // multiply mass

    struct lnode nd_copy = *nd; 

    // create product node
    // struct lnode out;
    nd->op = '*';
    nd->n = 2;
    nd->son = (struct lnode*) malloc(sizeof(struct lnode));
    *(nd->son)= nd_copy;
    nd->bro = NULL;

    // crate symbol son node: lnode = m_{sym_idx}
    nd->son->bro = (struct lnode*) malloc(sizeof(struct lnode));
    create_lnode(nd->son->bro, 0, 's');
    nd->son->bro->number = sym_idx;
    nd->son->bro->bro = NULL;

  } else if (mul_m == 0 && mul_2 == 1) {
    // multiply 2

    struct lnode nd_copy = *nd; 

    // create product node
    // struct lnode out;
    nd->op = '*';
    nd->n = 2;
    nd->son = (struct lnode*) malloc(sizeof(struct lnode));
    *(nd->son)= nd_copy;
    nd->bro = NULL;

    // crate numeric son node: lnode = 2
    nd->son->bro = (struct lnode*) malloc(sizeof(struct lnode));
    create_lnode(nd->son->bro, 0, 'n');
    mpz_init(nd->son->bro->mpz);
    mpz_set_ui(nd->son->bro->mpz, 2);
    nd->son->bro->bro = NULL;
    
  } else if (mul_m == 1 && mul_2 == 1) {
    // multiply both

    struct lnode nd_copy = *nd; 

    // create product node
    // struct lnode out;
    nd->op = '*';
    nd->n = 3;
    nd->son = (struct lnode*) malloc(sizeof(struct lnode));
    *(nd->son)= nd_copy;
    nd->bro = NULL;

    // crate numeric son node: lnode = 2
    nd->son->bro = (struct lnode*) malloc(sizeof(struct lnode));
    create_lnode(nd->son->bro, 0, 'n');
    mpz_init(nd->son->bro->mpz);
    mpz_set_ui(nd->son->bro->mpz, 2);

    // crate symbol son node: lnode = m_{sym_idx}
    nd->son->bro->bro = (struct lnode*) malloc(sizeof(struct lnode));
    create_lnode(nd->son->bro->bro, 0, 's');
    nd->son->bro->bro->number = sym_idx;
    nd->son->bro->bro->bro = NULL;
  }

  return;
}


void wrt_cmp_DE(
  // OUTPUT
  int **perm,
  // INPUT
  int *zero_label, int *nroots, mpc_t **roots,
  poly_frac ***pfmat, int eps_num, int dim,
  double wp2_rel_decr,
  char *file_ext, char *filepath_matrix, char *filepath_roots,
  int opt_write
) {
  char tmp_filepath[200];
  for (int ep=0; ep<eps_num; ep++) {
    if (opt_write) {
      //////
      // WRITE TO FILE
      //////
      
      // matrix
      snprintf(tmp_filepath, sizeof(tmp_filepath), "%s%d%s", filepath_matrix, ep, file_ext);
      cout << endl; cout << "writing to " << tmp_filepath << endl;
      poly_frac_rk2_to_file(
        tmp_filepath,
        pfmat[ep], dim, dim
      );
      
      // roots
      snprintf(tmp_filepath, sizeof(tmp_filepath), "%s%d%s", filepath_roots, ep, file_ext);
      cout << "writing to " << tmp_filepath << endl;
      // mpc_rk1_to_file(tmp_filepath, roots, nroots);
      int_rk0_mpc_rk1_to_file(tmp_filepath, roots[ep], nroots[ep], zero_label[ep]);
      // branch sing labs
      // snprintf(tmp_filepath, sizeof(tmp_filepath), "%s%d%s", filepath_branch_sing_lab, ep, file_ext);
      // int_rk2_to_file(tmp_filepath, branch_sing_lab_ep, 3, 2);
    } else {
      //////
      // COMPARE
      //////
      
      // enlarge tolerance
      mpfr_t mpfr_tol_orig;
      mpfr_init2(mpfr_tol_orig, wp2);
      mpfr_set(mpfr_tol_orig, mpfr_tol, MPFR_RNDN);
      // cout << endl; cout << "pre-enlarged mpfr tol:" << endl;
      // mpfr_out_str(stdout, 10, 0, mpfr_tol, MPFR_RNDN); cout << endl;
      mpfr_tol_enlarge(wp2_rel_decr);
      // cout << endl; cout << "check with enlarged mpfr tol:" << endl;
      // mpfr_out_str(stdout, 10, 0, mpfr_tol, MPFR_RNDN); cout << endl;

      // ROOTS
      mpc_t *bench_roots;
      int bench_nroots, bench_zero_label;
      snprintf(tmp_filepath, sizeof(tmp_filepath), "%s%d%s", filepath_roots, ep, file_ext);
      cout << endl; cout << "reading mpc_t roots from file " << tmp_filepath << endl;
      bench_nroots = count_lines(tmp_filepath) - 1;
      bench_roots = new mpc_t[bench_nroots];
      init_rk1_mpc(bench_roots, bench_nroots);
      int_rk0_mpc_rk1_from_file(tmp_filepath, bench_roots, bench_nroots, &bench_zero_label);
      cout << "CHECK roots..." << endl;
      if (bench_nroots != nroots[ep]) {
        cout << "CHECK: FAIL" << endl;
        cout << "nroots = " << nroots[ep] << endl;
        cout << "benchmark nroots = " << bench_nroots << endl;
        cout << "benchmark roots:" << endl;
        print_poly(bench_roots, bench_nroots-1);
        exit(1);
      }
      if (bench_zero_label != zero_label[ep]) {
        cout << "CHECK: FAIL" << endl;
        cout << "zero_label = " << zero_label[ep] << endl;
        cout << "benchmark zero_label = " << bench_zero_label << endl;
        exit(1);
      }
      perm[ep] = new int[nroots[ep]];
      mpc_rk1_compare_perm(perm[ep], bench_nroots, bench_roots, roots[ep]);
      // cout << "permutation:" << endl;
      // for (int k=0; k<bench_nroots; k++) {
      //   cout << "k[" << k << "] = " << perm[k] << endl;
      // }

      // POLY_FRAC MATRICES
      struct poly_frac **bench_pfmat;
      malloc_rk2_tens(bench_pfmat, dim, dim);
      poly_frac_rk2_build(bench_pfmat, dim, dim);
      snprintf(tmp_filepath, sizeof(tmp_filepath), "%s%d%s", filepath_matrix, ep, file_ext);
      cout << endl; cout << "reading poly_frac DE from file " << tmp_filepath << endl;
      poly_frac_rk2_from_file(tmp_filepath, bench_pfmat, dim, dim);
      cout << "CHECK pfmat..." << endl;
      poly_frac_rk2_perm_roots(bench_pfmat, perm[ep], dim, dim);
      poly_frac_rk2_prune(bench_pfmat, dim, dim, wp2_rel_decr);
      poly_frac_rk2_normal(bench_pfmat, dim, dim);
      poly_frac_rk2_compare(bench_pfmat, pfmat[ep], dim, dim, roots[ep]);

      // restore original tolerance
      mpfr_set(mpfr_tol, mpfr_tol_orig, MPFR_RNDN);
    }
  }
}


// <<generate_poly_frac_DE>>
void generate_poly_frac_DE(
  // OUTPUT
  poly_frac ***pfmat,
  int *zero_label, int *nroots, mpc_t **roots, mpfr_t **tols,
  // INPUT
  int nroots_branch, mpc_t *roots_branch, mpfr_t *tols_branch,
  int ninvs, char **symbols, int *is_mass,
  poly_frac *pspf,
  int *skip_inv, char ***ep_kin,
  int dim, char ****mats_str,
  int nbranches, int *branch_deg,
  int eps_num, char **eps_str,
  FILE *terminal
  // ,
  // char *file_ext, char *filepath_matrix, char *filepath_roots,
  // int opt_write
) {
	// int nroots = 1;
	// mpc_t *roots = new mpc_t[1];
	// mpc_init3(roots[0], wp2, wp2);
	// mpc_set_ui(roots[0], 0, MPFR_RNDN);
	// mpfr_t *tols = new mpfr_t[1];
	// mpfr_init2(tols[0], wp2);
	// mpfr_set(tols[0], mpfr_tol, MPFR_RNDN);

  poly_frac ***pfmats;
  malloc_rk3_tens(pfmats, ninvs, dim, dim);

  // int **branch_sing_lab_ep;
  // malloc_rk2_tens(branch_sing_lab_ep, nbranches, 3);

	// double wp2_rel_decr = 1.0;
  double wp2_rel_decr = 0.93;
  // double wp2_rel_decr = 0.60;

  // // enlarge tolerance for decoding
  // mpfr_t mpfr_tol_orig;
  // mpfr_init2(mpfr_tol_orig, wp2);
  // mpfr_set(mpfr_tol_orig, mpfr_tol, MPFR_RNDN);
  // mpfr_tol_enlarge(wp2_rel_decr);

  //////
  // PARSE
  //////
  cout << endl; cout << "PARSING..." << endl;
  struct lnode ***mats_nd;
  malloc_rk3_tens(mats_nd, ninvs, dim, dim);
  for (int s=0; s<ninvs; s++) {
    // cout << "----------------------------------------" << endl;
    // cout << "s = " << s << endl;
    if (skip_inv[s] == 1) {
      continue;
    }
    for (int i=0; i<dim; i++) {
      for (int j=0; j<dim; j++) {
        //////
        // ENCODE
        //////
        // cout << "i, j = " << i << ", " << j << endl;
        // encode with parser
        // if (s == 0 && i == 26 && j == 22) {
        // cout << "i, j = " << i << ", " << j << endl;
        // cout << endl; cout << "input expresion (for parser):" << endl;
        // cout << mats_str[s][i][j] << endl;
        // }
        struct lnode nd;
        // cout << "parsing expression..." << endl;
        // fflush(stdout);
        lnode_parse_expression(&mats_nd[s][i][j], mats_str[s][i][j], (const char**) symbols, ninvs+1, 0);
        // if (s == 0 && i == 26 && j == 22) {
        // cout << "parsed tree:" << endl;
        // lnode_print(&mats_nd[s][i][j], (const char**) symbols); printf("\n");
        // }
        linearize_masses_nd(&mats_nd[s][i][j], is_mass);
        // if (s == 0 && i == 26 && j == 22) {
        // cout << "linearized masses:" << endl;
        // lnode_print(&mats_nd[s][i][j], (const char**) symbols); printf("\n");
        // }
        // cout << "expanding..." << endl;
        lnode_expand(&mats_nd[s][i][j]);
        // if (s == 0 && i == 26 && j == 22) {
        // printf("expanded:\n");
        // lnode_print(&mats_nd[s][i][j], (const char**) symbols); printf("\n");
        // }
        if (is_mass[s]) {
          // A_m = 2*m*A_m2
          linearize_mass_DE(&mats_nd[s][i][j], s+1);
          // if (s == 0 && i == 26 && j == 22) {
          // printf("DE mass linearized:\n");
          // lnode_print(&mats_nd[s][i][j], (const char**) symbols); printf("\n");
          // }
        }
        // printf("encoded:\n");
        // printf("%s\n", lnode_to_str(&nd, (char*)", "));
      }
    }
  }

  //////
  // DECODE
  //////
  cout << "DECODING..." << endl;
  // int nroots_ini = *nroots;

  fprintf(terminal, "process mat: "); fflush(terminal); usleep(sleep_time);
  for (int ep=0; ep<eps_num; ep++) {
    if (ep > 0) {fprintf(terminal, "\033[22D\033[K");}// fflush(terminal); usleep(sleep_time);}
    fprintf(terminal, "eps value %3d /%3d... ", ep, eps_num-1); fflush(terminal); usleep(sleep_time);
    fprintf(stdout, "\n############################################## ep = %d\n", ep);
		
		ep_kin[0][0] = strdup(eps_str[ep]);
		// cout << "eps = " << ep_kin[0][0] << endl;
		poly_frac_build(&pspf[0]);
		generate_PS_pf(&pspf[0], ep_kin, 0, 0, 0, 1, wp2);
		// cout << "eps pf:" << endl;
		// poly_frac_print(&pspf[0]);

    // reset roots to branch point roots
    // *nroots = nroots_ini;
    nroots[ep] = nroots_branch;
    roots[ep] = new mpc_t[nroots_branch];
    init_rk1_mpc(roots[ep], nroots_branch);
    copy_rk1_mpc(roots[ep], roots_branch, nroots_branch);
    tols[ep] = new mpfr_t[nroots_branch];
    init_rk1_mpfr(tols[ep], nroots_branch);
    copy_rk1_mpfr(tols[ep], tols_branch, nroots_branch);

    char *nd_str;
		for (int s=0; s<ninvs; s++) {
			// cout << "----------------------------------------" << endl;
			// cout << "s = " << s << endl;
      if (skip_inv[s] == 1) {
        continue;
      }
			poly_frac_rk2_build(pfmats[s], dim, dim);
			for (int i=0; i<dim; i++) {
				for (int j=0; j<dim; j++) {
          // cout << "i, j = " << i << ", " << j << endl;
					nd_str = lnode_to_str(&mats_nd[s][i][j], (char*)", ");
          // if (s == 0 && i == 49 && j == 0) {
          // cout << "i, j = " << i << ", " << j << endl;
					// printf("%s\n", nd_str);
          // dbg = 1;
          // }
          // fflush(stdout);
					// continue;

					//////
					// DECODE
					//////
					decode_tree_pf(
						&pfmats[s][i][j],
						&nroots[ep], &roots[ep], &tols[ep],
						pspf,
						ep_kin, skip_inv-1, is_mass-1,
						&nd_str, (char*) ", ",
						wp2, mpfr_tol, 0
					);
          // dbg = 0;
          // if (s == 0 && i == 26 && j == 22) {
					// cout << "poly frac:" << endl;
					// poly_frac_print(&pfmats[s][i][j]);
          // }

					free(nd_str);
				}
			}
      // cout << "root results:" << endl;
      // print_poly(roots, nroots-1);
		}
    if (mpc_lessthan_tol(roots[ep][0])) {
      mpc_set_ui(roots[ep][0], 0, MPFR_RNDN);
    }
		// cout << endl; cout << "ROOTS:" << endl;
		// cout << "nroots = " << nroots[ep] << endl;
		// print_poly(roots[ep], nroots[ep]-1);

    // cout << "POLY_FRAC MATRICES" << endl;
    // for (int s=0; s<ninvs; s++) {
    //   cout << "s = " << s << endl;
    //   poly_frac_rk2_print(pfmats[s], dim, dim);
    // }

		//////
		// COMBINE POLY_FRAC MATRICES TOGETHER
		//////
    int wp_bin = - mpfr_log2_int(mpfr_tol);
    for (int i=0; i<dim; i++) {
      for (int j=0; j<dim; j++) {
        poly_frac_set_zero(&pfmat[ep][i][j]);
        for (int s=0; s<ninvs; s++) {
          if (skip_inv[s]) {
            continue;
          }
					poly_frac_roots_append_zero(&pfmats[s][i][j], nroots[ep]);
          rel_err_poly_frac_add_pf(
            &pfmat[ep][i][j], &pfmat[ep][i][j], &pfmats[s][i][j],
            roots[ep], &pspf[1+s].coeffs[pspf[1+s].num_vdeg],
            wp_bin
          );
        }
      }
    }

    // // restore tolerance
    // mpfr_set(mpfr_tol, mpfr_tol_orig, MPFR_RNDN);

		// PRUNE
		// poly_frac_rk2_prune(pfmat[ep], dim, dim, wp2_rel_decr);
    // poly_frac_rk2_normal(pfmat[ep], dim, dim);
    poly_frac_rk2_prune_rel_tol_real_max(pfmat[ep], wp_bin, dim, dim);
    mpc_rk1_prune_re_im(roots[ep], wp_bin, nroots[ep]);  // #2BD prune each single root as soon as it is discovered

    // zero_label
    zero_label[ep] = -1;
    for (int i=0; i<dim; i++) {
      for (int j=0; j<dim; j++) {
        if (pfmat[ep][i][j].mults) {
          if (pfmat[ep][i][j].mults[0] != 0) {
            zero_label[ep] = 0;
            break;
          }
        }
      }
      if (zero_label[ep] == 0) {
        break;
      }
    }
  }
  fprintf(terminal, "\033[22D\033[K"); fflush(terminal); usleep(sleep_time);
  fprintf(terminal, "\033[13D\033[K"); fflush(terminal); usleep(sleep_time);

}


void build_block_indices(int *&b_idx, int *&sb_idx, int *b_len, int *sb_len, int b, int** prof, int **sb_grid) {
  int row_start = prof[b][0], row_end = prof[b][1];
  *b_len = row_end - row_start + 1;
  *sb_len = 0;
  for (int sb=0; sb<b; sb++) {
    if (sb_grid[b][sb]) {
      *sb_len = *sb_len + prof[sb][1] - prof[sb][0] + 1;
    }
  }

  b_idx = new int[*b_len];
  sb_idx = new int[*sb_len];
  // build block indices
  for (int i=0; i<row_end-row_start+1; i++) {
    b_idx[i] = prof[b][0] + i;
  }
  // build sub-block indices
  int l=0;
  for (int sb=0; sb<b; sb++) {
    if (sb_grid[b][sb]) {
      for (int j=prof[sb][0]; j<=prof[sb][1]; j++) {
        sb_idx[l] = j;
        l++;
      }
    }
  }
}


template <typename type>
void select_block(type **&block, int b, int sb, type **mat, int** prof) {
  int row_start = prof[b][0], row_end = prof[b][1];
  int col_start = prof[sb][0], col_end = prof[sb][1];

  // copy by reference
  block = new type*[row_end - row_start + 1];
  for (int i=0; i<row_end-row_start+1; i++) {
    block[i] = new type[col_end - col_start +1];
    block[i] = &(mat[row_start + i][col_start]);
  }

  // // copy by values
  // block = new ex*[end - start + 1];
  // for (int i=0; i<end-start+1; i++) {
  //   block[i] = new ex[end - start +1];
  //   for (int j=0; j<end-start+1; j++) {
  //     block[i][j] = mat[row_start + i][col_start + j];
  //   }
  // }
}


void generate_constant_term(
  // OUTPUT
  mpc_t **const_term,
  // INPUT
  mpc_t ***sblock_num, mpc_t ***sblock_den, mpc_t **solutions,
  int b_len, int sb_len, int *sb_idx, int num_max_deg, int den_max_deg, int eta_ord,
  int log_sb_idx_dim = 0, int *log_sb_idx = NULL
) {
  bool solve_in_zero = true;
  if (log_sb_idx == NULL) {
    solve_in_zero = false;
    log_sb_idx_dim = sb_len;
    log_sb_idx = (int *) new int[sb_len];
    for (int j=0; j<sb_len; j++) {
      log_sb_idx[j] = j;
    }
  }
  // cout << "solve in zero = " << solve_in_zero << endl;
  mpc_t *tmp_pol = new mpc_t[eta_ord+1];
  init_poly(tmp_pol, eta_ord);
  for (int i=0; i<b_len; i++) {
    // cout << "i = " << i << endl;
    set_null_poly(const_term[i], eta_ord);
    for (int j=0; j<log_sb_idx_dim; j++) {
      // cout << "j = " << j << endl;
      // cout << "log_sb_idx = " << log_sb_idx[j] << endl;
      // cout << "sb_idx = " << sb_idx[log_sb_idx[j]] << endl;
      
      // div_poly(tmp_pol, sblock_num[i][j], sblock_den[i][j], eta_ord, num_max_deg, den_max_deg);
      // mul_poly(tmp_pol, tmp_pol, solutions[sb_idx[j]], eta_ord, eta_ord);
      // if (solve_in_zero) {
      //   cout << "solutions = " << endl;
      //   print_poly(solutions[sb_idx[log_sb_idx[j]]], 10);
      //   cout << "sblock num = " << endl;
      //   print_poly(sblock_num[i][log_sb_idx[j]], num_max_deg);
      //   cout << "sblock den = " << endl;
      //   print_poly(sblock_den[i][log_sb_idx[j]], den_max_deg);
      // }
      // cout << "1st factor (solutions):" << endl;
      // print_poly(solutions[sb_idx[log_sb_idx[j]]], 10);
      // cout << "2nd factor (numerators):" << endl;
      // print_poly(sblock_num[i][log_sb_idx[j]], num_max_deg);
      mul_poly(tmp_pol, solutions[sb_idx[log_sb_idx[j]]], sblock_num[i][log_sb_idx[j]], eta_ord, num_max_deg);
      // cout << "done mul" << endl;
      // print_poly(tmp_pol, 1);
      div_poly(tmp_pol, tmp_pol, sblock_den[i][log_sb_idx[j]], eta_ord, den_max_deg);
      // cout << "done div" << endl;
      // print_poly(tmp_pol, 0);
      add_poly(const_term[i], const_term[i], tmp_pol, eta_ord, eta_ord);
      // cout << "done add" << endl;
      // print_poly(const_term[i], 0);
    }
  }
  delete[] tmp_pol;
}


void sys_block_info(
  // OUTPUT
  mpc_t *&lcm, mpc_t ***&block, mpc_t ****&const_term,
  // IN-OUT
  bool *flag_lcm, int *lcm_deg, mpc_t *&lcm_orig, int *&lcm_roots,
  int *block_max_deg, int *num_max_deg, int *den_max_deg,
  // INPUT
  int b_len, int sb_len,
  int *b_idx, int *sb_idx,
  struct poly_frac **lcm_mat_ep, mpc_t ****mpc_lcm_mat_ep,
  int npoles, mpc_t *poles,
  mpc_t eta_shift,
  int eta_ord,
  mpc_t ****solutions, int num_prev_eig=-1, int *prev_eig=NULL,
  int *prev_log_len=NULL, int *sol_log_len=NULL, int **log_prof=NULL
) {
  //////
  // INPUT
  // - to load info for solving in regular points, pass the address of a rank-2 tensor as
  //    "solutions", the 1st entry of a rank-3 tensor as "const_term" and
  //    leave "num_prev_eig" and "prev_eig" as default;

  //////
  // LCM
  //////
  int **nrem_roots;
  // ex num_den, gnc_expr;
  mpc_t tmpc;
  mpc_init3(tmpc, wp2, wp2);
  int ***rem_mults;
  int count;
  int *lcm_mults;
  int ldeg;
  if (*flag_lcm) {
    // get multiplicites of the lcm and compute its degree
    lcm_mults = new int[npoles];
    for (int k=0; k<npoles; k++) {
      lcm_mults[k] = 0;
    }
    for (int i=0; i<b_len; i++) {
      for (int j=0; j<b_len; j++) {
        if (lcm_mat_ep[b_idx[i]][b_idx[j]].num_vdeg == -1) {
          continue;
        }
        for (int k=0; k<npoles; k++) {
          lcm_mults[k] = max(lcm_mults[k], lcm_mat_ep[b_idx[i]][b_idx[j]].mults[k]);
        }
      }
    }
    *lcm_deg = 0;
    for (int k=0; k<npoles; k++) {
      *lcm_deg += lcm_mults[k];
    }
    // cout << "lcm deg = " << *lcm_deg << endl;
    // compute lcm coefficients
    lcm_roots = new int[*lcm_deg];
    lcm_orig = new mpc_t[*lcm_deg + 1];
    init_rk1_mpc(lcm_orig, *lcm_deg + 1);
    mpc_set_d(lcm_orig[0], 1, MPFR_RNDN);
    int tmp_lcm_deg = 0;
    for (int k=0; k<npoles; k++) {
      for (int kp=0; kp<lcm_mults[k]; kp++) {
        lcm_roots[tmp_lcm_deg] = k;
        poly_mul_root(lcm_orig, tmp_lcm_deg, poles[k]);
        tmp_lcm_deg++;
      }
    }
    // cout << "lcm roots:" << endl;
    // for (int k=0; k<*lcm_deg; k++) {
    //   print_mpc(&poles[out_lcm_roots[k]]);
    //   cout << endl;
    // }
    // cout << "mpc lcm:" << endl;
    // print_poly(lcm_orig, *lcm_deg);

    // compute numerator max degree
    malloc_rk2_tens(nrem_roots, b_len, b_len);
    malloc_rk3_tens(rem_mults, b_len, b_len, npoles);
    *block_max_deg = 0;
    for (int i=0; i<b_len; i++) {
      for (int j=0; j<b_len; j++) {
        nrem_roots[i][j] = 0;
        if (lcm_mat_ep[b_idx[i]][b_idx[j]].num_vdeg == -1) {
          for (int k=0; k<npoles; k++) {
            rem_mults[i][j][k] = 0;
            nrem_roots[i][j] += 0;
          }
        } else {
          for (int k=0; k<npoles; k++) {
            rem_mults[i][j][k] = lcm_mults[k] - lcm_mat_ep[b_idx[i]][b_idx[j]].mults[k];
            nrem_roots[i][j] += rem_mults[i][j][k];
          }
          *block_max_deg = max(*block_max_deg, lcm_mat_ep[b_idx[i]][b_idx[j]].num_deg + nrem_roots[i][j]);
        }
      }
    }
    // cout << "max degree = " << *block_max_deg << endl;
    
    // turn into coefficients and multiply times the lcm
    for (int i=0; i<b_len; i++) {
      for (int j=0; j<b_len; j++) {
        // cout << "i, j = " << i << ", " << j << endl;
        mpc_lcm_mat_ep[b_idx[i]][b_idx[j]][0] = new mpc_t[*block_max_deg+1];
        init_rk1_mpc(mpc_lcm_mat_ep[b_idx[i]][b_idx[j]][0], *block_max_deg+1);
        // turn into coefficients
        // cout << "turn into coefficients" << endl;
        if (lcm_mat_ep[b_idx[i]][b_idx[j]].num_vdeg == -1) {
          for (int k=0; k<=*block_max_deg; k++) {
            mpc_set_ui(
              mpc_lcm_mat_ep[b_idx[i]][b_idx[j]][0][k], 0, MPFR_RNDN
            );
          }
        } else {
          ldeg = lcm_mat_ep[b_idx[i]][b_idx[j]].num_deg - lcm_mat_ep[b_idx[i]][b_idx[j]].num_vdeg;
          for (int k=0; k<ldeg; k++) {
            mpc_set_d(mpc_lcm_mat_ep[b_idx[i]][b_idx[j]][0][k], 0, MPFR_RNDN);
          }
          for (int k=0; k<=lcm_mat_ep[b_idx[i]][b_idx[j]].num_vdeg; k++) {  
            mpc_set(
              mpc_lcm_mat_ep[b_idx[i]][b_idx[j]][0][ldeg+k],
              lcm_mat_ep[b_idx[i]][b_idx[j]].coeffs[k],
              MPFR_RNDN
            );
          }
          // multiply times roots not present in the denominator
          // cout << "multiply times roots not present in the denominator" << endl;
          count = 0;
          for (int k=0; k<npoles; k++) {
            for (int m=0; m<rem_mults[i][j][k]; m++) {
              // cout << "k, m = " << k << ", " << m << endl;
              poly_mul_root(
                mpc_lcm_mat_ep[b_idx[i]][b_idx[j]][0],
                lcm_mat_ep[b_idx[i]][b_idx[j]].num_deg + count,
                poles[k]
              );
              count++;
            }
          }
          // set to zero remaining powers
          // cout << "set to zero remaining powers" << endl;
          for (int k=lcm_mat_ep[b_idx[i]][b_idx[j]].num_deg+nrem_roots[i][j]+1; k<=*block_max_deg; k++) {
            mpc_set_d(mpc_lcm_mat_ep[b_idx[i]][b_idx[j]][0][k], 0, MPFR_RNDN);
          }
        }
      }
    }
    del_rk2_tens(nrem_roots, b_len);
    del_rk3_tens(rem_mults, b_len, b_len);
    // getchar();
    // cout << "done flag bock" << endl;
  }

  //////
  // BLOCK
  //////
  // select block numerators and turn them into coefficients
  malloc_rk3_tens(block, b_len, b_len, *block_max_deg + 1);
  init_rk3_mpc(block, b_len, b_len, *block_max_deg + 1);
  // mpc_t eta_rescale;
  // mpc_init3(eta_rescale, wp2, wp2);
  // gnc_to_mpc(&eta_rescale, pow(denom(etav), *lcm_deg));
  // cout << "mpc coefficients after lcm:" << endl;
  for (int i=0; i<b_len; i++) {
    for (int j=0; j<b_len; j++) {
      // cout << "i, j = " << i << ", " << j << endl;
      for (int k=0; k<=*block_max_deg; k++) {
      // for (int k=0; k<=mat_deg[b_idx[i]][b_idx[j]][0]+nrem_roots[i][j]; k++) {
        // print_mpc(&mpc_lcm_mat_ep[b_idx[i]][b_idx[j]][0][k]);
        // cout << endl;
        mpc_set(block[i][j][k], mpc_lcm_mat_ep[b_idx[i]][b_idx[j]][0][k], MPFR_RNDN);
        // print_mpc(&block[i][j][k]);
        // cout << endl;
      }
      // shift
      poly_shift(eta_shift, block[i][j], *block_max_deg);
      // // rescale
      // for (int k=0; k<=*block_max_deg; k++) {
      //   mpc_mul(block[i][j][k], block[i][j][k], eta_rescale, MPFR_RNDN);
      // }
    }
  }

  // // print block
  // cout << "block coefficients:" << endl;
  // print_rk3_mpc(block, b_len, b_len, (*block_max_deg)+1);
  // getchar();
  
  // rk2gnc_to_rk3mpc(block, lcm_mat_ep, b_idx, b_idx, b_len, b_len, *block_max_deg + 1, g, etav, *lcm_deg);
  // // print block
  // cout << "block coefficients:" << endl;
  // print_rk3_mpc(block, b_len, b_len, (*block_max_deg)+1);
  // getchar();

  // SHIFT LCM
  // ex shifted_g;
  // shifted_g = (*g).subs(lst{g_et == g_et + etav});
  // cout << "ginac shifted LCM: " << shifted_g << endl;
  // shifted_g = (shifted_g) * pow(denom(etav), *lcm_deg);
  // cout << "LCM of traslation: " << shifted_g << endl;
  // cout << "LCM degree = " << *lcm_deg << endl;
  lcm = new mpc_t[*lcm_deg + 1];
  for (int k=0; k<=(*lcm_deg); k++) {
    // cout << "k = " << k << endl;
    mpc_init3(lcm[k], wp2, wp2);
    mpc_set(lcm[k], lcm_orig[k], MPFR_RNDN);
  }
  poly_shift(eta_shift, lcm, *lcm_deg);
  // for (int k=0; k<=*lcm_deg; k++) {
  //   mpc_mul(lcm[k], lcm[k], eta_rescale, MPFR_RNDN);
  //   // cout << coeff(shifted_g.expand(), g_et, k) << endl;
  // }
  // cout << "mpc shifted lcm" << endl;
  // print_poly(lcm, *lcm_deg);
  // getchar();

  // cout << "denom eta = " << denom(etav) << endl;
  // cout << "pow denom eta = " << pow(denom(etav), *lcm_deg) << endl;
  // ex etav_num_den = numer_denom(etav);
  // mpc_init3(mpc_etav_den, wp2, wp2);
  // gnc_to_mpc(&mpc_etav, etav_num_den[0]);
  // gnc_to_mpc(&mpc_etav_den, etav_num_den[1]);
  // mpc_div(mpc_etav, mpc_etav, mpc_etav_den, MPFR_RNDN);
  // // cout << etav << endl;
  // // print_mpc(&mpc_etav);
  // // cout << endl;
  // mpc_t eta_rescale;
  // for (int k=0; k<=*lcm_deg; k++) {
  //   mpc_div(lcm[k], lcm[k], eta_rescale, MPFR_RNDN);
  // }
  // poly_shift(eta_shift, lcm, *lcm_deg);
  // gnc_to_mpc(&mpc_etav, pow(denom(etav), *lcm_deg));
  // for (int k=0; k<=*lcm_deg; k++) {
  //   mpc_mul(lcm[k], lcm[k], eta_rescale, MPFR_RNDN);
  //   cout << coeff(shifted_g.expand(), g_et, k) << endl;
  // }
  // print_poly(lcm, *lcm_deg);
  // getchar();

  // if sb_len is zero, no constant term is needed
  if (sb_len == 0) {
    return;
  }

  //////
  // SUB-BLOCKS
  //////
  int nden_rem_roots;
  if (*flag_lcm) {
    malloc_rk2_tens(nrem_roots, b_len, sb_len);
    malloc_rk3_tens(rem_mults, b_len, sb_len, npoles);
    *num_max_deg = 0;
    *den_max_deg = 0;
    for (int i=0; i<b_len; i++) {
      for (int j=0; j<sb_len; j++) {
        // cout << "i, j = " << i << ", " << j << endl;
        // cout << "lcm mults" << endl;
        // for (int k=0; k<npoles; k++) {
        //   cout << "k = " << k << ": " << lcm_mults[k] << endl;
        // }
        // cout << "den mults" << endl;
        // for (int k=0; k<npoles; k++) {
        //   cout << "k = " << k << ": " << lcm_mat_ep[b_idx[i]][sb_idx[j]].mults[k] << endl;
        // }
        // cout << "den deg = " << mat_deg[b_idx[i]][sb_idx[j]][1] << endl;
        nrem_roots[i][j] = 0;
        if (lcm_mat_ep[b_idx[i]][sb_idx[j]].num_vdeg == -1) {
          continue;
          // for (int k=0; k<npoles; k++) {
          //   rem_mults[i][j][k] = 0;
          //   nrem_roots[i][j] += 0;
          // }
        } else {
          for (int k=0; k<npoles; k++) {
            // cout << "k = " << k << endl;
            rem_mults[i][j][k] = lcm_mults[k] - lcm_mat_ep[b_idx[i]][sb_idx[j]].mults[k];
            if (rem_mults[i][j][k] < 0) {
              rem_mults[i][j][k] = 0;
            } else if (rem_mults[i][j][k] > 0) {
              nrem_roots[i][j] += rem_mults[i][j][k];
            }
          }
          *num_max_deg = max(*num_max_deg, lcm_mat_ep[b_idx[i]][sb_idx[j]].num_deg + nrem_roots[i][j]);
          *den_max_deg = max(*den_max_deg, lcm_mat_ep[b_idx[i]][sb_idx[j]].den_deg + nrem_roots[i][j] - *lcm_deg);
        }
      }
    }
    // cout << "num max degree = " << *num_max_deg << endl;
    // cout << "den max degree = " << *den_max_deg << endl;

    // turn into coefficients and multiply times the LCM
    for (int i=0; i<b_len; i++) {
      for (int j=0; j<sb_len; j++) {
        // cout << "i, j = " << i << ", " << j << endl;
        // NUMERATOR
        mpc_lcm_mat_ep[b_idx[i]][sb_idx[j]][0] = new mpc_t[*num_max_deg+1];
        init_rk1_mpc(mpc_lcm_mat_ep[b_idx[i]][sb_idx[j]][0], *num_max_deg+1);
        // turn numerator into coefficients
        // cout << "turn numerator into coefficients" << endl;
        if (lcm_mat_ep[b_idx[i]][sb_idx[j]].num_vdeg == -1) {
          for (int k=0; k<=*num_max_deg; k++) {
            mpc_set_ui(
              mpc_lcm_mat_ep[b_idx[i]][sb_idx[j]][0][k], 0, MPFR_RNDN
            );
          }
        } else {
          ldeg = lcm_mat_ep[b_idx[i]][sb_idx[j]].num_deg - lcm_mat_ep[b_idx[i]][sb_idx[j]].num_vdeg;
          for (int k=0; k<ldeg; k++) {
            mpc_set_d(mpc_lcm_mat_ep[b_idx[i]][sb_idx[j]][0][k], 0, MPFR_RNDN);
          }
          for (int k=0; k<=lcm_mat_ep[b_idx[i]][sb_idx[j]].num_vdeg; k++) {
            mpc_set(
              mpc_lcm_mat_ep[b_idx[i]][sb_idx[j]][0][ldeg+k],
              lcm_mat_ep[b_idx[i]][sb_idx[j]].coeffs[k],
              MPFR_RNDN
            );
          }
          // multiply numerator times the lcm
          // cout << "multiply numerator times the lcm" << endl;
          count = 0;
          for (int k=0; k<npoles; k++) {
            for (int m=0; m<rem_mults[i][j][k]; m++) {
              // cout << "k, m = " << k << ", " << m << endl;
              // cout << "poly in:" << endl;
              // print_poly(mpc_lcm_mat_ep[b_idx[i]][sb_idx[j]][0], lcm_mat_ep[b_idx[i]][sb_idx[j]].num_deg + count);
              poly_mul_root(
                mpc_lcm_mat_ep[b_idx[i]][sb_idx[j]][0],
                lcm_mat_ep[b_idx[i]][sb_idx[j]].num_deg + count,
                poles[k]
              );
              count++;
            }
          }
          // set to zero remaining powers of the numenator
          // cout << "set to zero remaining powers of the numenator" << endl;
          for (int k=lcm_mat_ep[b_idx[i]][sb_idx[j]].num_deg+nrem_roots[i][j]+1; k<=*num_max_deg; k++) {
            mpc_set_d(mpc_lcm_mat_ep[b_idx[i]][sb_idx[j]][0][k], 0, MPFR_RNDN);
          }
        }


        // DENOMINATOR
        // cout << "DENOMINATOR" << endl;
        mpc_lcm_mat_ep[b_idx[i]][sb_idx[j]][1] = new mpc_t[*den_max_deg+1];
        init_rk1_mpc(mpc_lcm_mat_ep[b_idx[i]][sb_idx[j]][1], *den_max_deg+1);
        // turn denominator into coefficients
        mpc_set_ui(mpc_lcm_mat_ep[b_idx[i]][sb_idx[j]][1][0], 1, MPFR_RNDN);
        if (lcm_mat_ep[b_idx[i]][sb_idx[j]].num_vdeg == -1) {
          for (int k=1; k<=*den_max_deg; k++) {
            mpc_set_ui(mpc_lcm_mat_ep[b_idx[i]][sb_idx[j]][1][k], 0, MPFR_RNDN);
          }
        }
        else {
          count = 0;
          for (int k=0; k<npoles; k++) {
            for (int m=0; m<lcm_mat_ep[b_idx[i]][sb_idx[j]].mults[k]-lcm_mults[k]; m++) {
              // cout << "k, m = " << k << ", " << m << endl;
              poly_mul_root(
                mpc_lcm_mat_ep[b_idx[i]][sb_idx[j]][1],
                count,
                poles[k]
              );
              count++;
            }
          }
          // set to zero remaining powers of the denominator
          for (int k=lcm_mat_ep[b_idx[i]][sb_idx[j]].den_deg+nrem_roots[i][j]-*lcm_deg+1; k<=*den_max_deg; k++) {
            mpc_set_d(mpc_lcm_mat_ep[b_idx[i]][sb_idx[j]][1][k], 0, MPFR_RNDN);
          } 
        }
      }
    }
    del_rk2_tens(nrem_roots, b_len);
    del_rk3_tens(rem_mults, b_len, sb_len);
  }

  // find max degree of sub-block numerators and denominators
  // int num_max_deg, den_max_deg;
  // find_num_max_degree(&num_max_deg, lcm_mat_ep, &shifted_g, etav, b_len, sb_len, b_idx, sb_idx);
  // find_den_max_degree(&den_max_deg, lcm_mat_ep, &shifted_g, etav, b_len, sb_len, b_idx, sb_idx);

  // select sub-block numerators and denominators and turn them into coefficient
  mpc_t ***sblock_num, ***sblock_den;
  malloc_rk3_tens(sblock_num, b_len, sb_len, *num_max_deg+1);
  init_rk3_mpc(sblock_num, b_len, sb_len, *num_max_deg+1);
  malloc_rk3_tens(sblock_den, b_len, sb_len, *den_max_deg+1);
  init_rk3_mpc(sblock_den, b_len, sb_len, *den_max_deg+1);
  for (int i=0; i<b_len; i++) {
    for (int j=0; j<sb_len; j++) {
      // cout << "i,j=" << i << "," << j << endl;
      // cout << "NUM" << endl;
      // print_poly(mpc_lcm_mat_ep[b_idx[i]][sb_idx[j]][0], *num_max_deg);
      for (int k=0; k<=*num_max_deg; k++) {
      // for (int k=0; k<=mat_deg[b_idx[i]][sb_idx[j]][0]+nrem_roots[i][j]; k++) {
        // cout << "k=" << k <<  endl;
        mpc_set(sblock_num[i][j][k], mpc_lcm_mat_ep[b_idx[i]][sb_idx[j]][0][k], MPFR_RNDN);
      }
      // shift
      poly_shift(eta_shift, sblock_num[i][j], *num_max_deg);
      // // rescale
      // for (int k=0; k<=*num_max_deg; k++) {
      //   mpc_mul(sblock_num[i][j][k], sblock_num[i][j][k], eta_rescale, MPFR_RNDN);
      // }
      // for (int k=0; k<=*den_max_deg; k++) {
      // cout << "DEN" << endl;
      for (int k=0; k<=*den_max_deg; k++) {
      // for (int k=0; k<=mat_deg[b_idx[i]][sb_idx[j]][1]+nrem_roots[i][j]-*lcm_deg; k++) {
        // cout << "k=" << k <<  endl;
        mpc_set(sblock_den[i][j][k], mpc_lcm_mat_ep[b_idx[i]][sb_idx[j]][1][k], MPFR_RNDN);
      }
      // shift
      poly_shift(eta_shift, sblock_den[i][j], *den_max_deg);
    }
  }

  // // print sub-block
  // cout << "sub-block numerators:" << endl;
  // print_rk3_mpc(sblock_num, b_len, sb_len, *num_max_deg+1);
  // cout << "sub-block denominators:" << endl;
  // print_rk3_mpc(sblock_den, b_len, sb_len, *den_max_deg+1);

  // rk2gnc_to_rk3mpc_num_den(sblock_num, sblock_den, lcm_mat_ep, b_idx, sb_idx, b_len, sb_len, *num_max_deg+1, *den_max_deg+1, g, etav, *lcm_deg);

  //////
  // CONSTANT TERM
  /////
  // cout << "GENERATE CONSTANT TERM" << endl;
  if (num_prev_eig == -1) {
    malloc_rk2_tens(const_term[0][0], b_len, eta_ord+1);
    generate_constant_term(
      const_term[0][0],
      sblock_num, sblock_den, *(*solutions),
      b_len, sb_len, sb_idx, *num_max_deg, *den_max_deg, eta_ord);
  } else {
    int *log_sb_idx = new int[sb_len];
    int log_sb_idx_dim;
    const_term = new mpc_t***[num_prev_eig];
    for (int lam, n=0; n<num_prev_eig; n++) {
      lam = prev_eig[n];
      // cout << "lam = " << lam << endl;
      // #2BD: se salviamo le massime potenze associate ad ogni lam e ad ogni posizione di MI, si dovrebbe poter evitare il count
      // e lo si rende anche più preciso, perché si può prevedere se una certa potenza non compare per via dell'annullarsi
      // di un elemento del sotto-blocco
      sol_log_len[lam] = count_log_len(mpfr_tol, solutions[lam], prev_log_len[lam], sb_len, eta_ord+1, sb_idx);
      // for (int i=0; i<sb_len; i++) {
      //   if(sol_log_len[lam] < log_prof[lam][sb_idx[i]]) {
      //     sol_log_len[lam] = log_prof[lam][sb_idx[i]];
      //   }
      // }
      if (sol_log_len[lam] == 0) {
        continue;
      }
      malloc_rk3_tens(const_term[n], sol_log_len[lam], b_len, eta_ord+1);
      for (int l=0; l<sol_log_len[lam]; l++) {
        // cout << "l = " << l << endl;
        // generate an alternative sb_idx using the log profile, in order to avoid
        // using the NaN's in the previous solutions
        log_sb_idx_dim = 0;
        for (int i=0; i<sb_len; i++) {
          // cout << "i = " << i << endl;
          // cout << "sb_idx = " << sb_idx[i] << endl;
          // cout << "log_profile = " << log_prof[lam][sb_idx[i]] << endl;
          if (log_prof[lam][sb_idx[i]] - 1 > l) {
            log_sb_idx[log_sb_idx_dim] = i;
            log_sb_idx_dim++;
          }
        }
        // for (int i=0; i<log_sb_idx_dim; i++) {
        //   cout << "log_sb_idx[" << i << "] = " << log_sb_idx[i] << endl;
        // }
        generate_constant_term(
          const_term[n][l],
          sblock_num, sblock_den, solutions[lam][l],
          b_len, sb_len, sb_idx, *num_max_deg, *den_max_deg, eta_ord,
          log_sb_idx_dim, log_sb_idx);
        // cout << "generated constant term:" << endl;
        // print_rk2_mpc(const_term[lam][l], b_len, eta_ord+1);
      }
    }
    delete[] log_sb_idx;
  }
  // delete numerators and denominators
  del_rk3_tens(sblock_num, b_len, sb_len);
  del_rk3_tens(sblock_den, b_len, sb_len);

}

