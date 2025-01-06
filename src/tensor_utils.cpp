#include <iostream>
#include "mpc.h"

using namespace std;

#include "setup.h"
#include "utils.h"
#include "global_vars.h"
// #include "malloc_defs.h"
#include "poly_frac.h"


//////
// PRINT
//////
void print_rk1_mpc(mpc_t *tens, int dim) {
  for (int i=0; i<dim; i++) {
    cout << "i = " << i << ": "; print_mpc(&tens[i]); cout << endl;
  }
}


void print_rk2_mpc(mpc_t **tens, int dim1, int dim2) {
  for (int i1=0; i1<dim1; i1++) {
    for (int i2=0; i2<dim2; i2++) {
      cout << "(i1, i2) = " << i1 << ", " << i2 << endl;
      print_mpc(&tens[i1][i2]);
      cout << endl;
    }
    cout << endl;
  }
}


void print_rk3_mpc(mpc_t ***tens, int dim1, int dim2, int dim3) {
  for (int i1=0; i1<dim1; i1++) {
    cout << "i1 = " << i1 << endl;
    for (int i2=0; i2<dim2; i2++) {
      for (int i3=0; i3<dim3; i3++) {
        cout << "(i2, i3) = " << i2 << ", " << i3 << endl;
        print_mpc(&tens[i1][i2][i3]);
        cout << endl;
      }
      cout << endl;
    }
    cout << endl;
  }
}


//////
// INITIALIZE WITH MULTIPLE PRECISION
//////
void init_rk1_mpfr(mpfr_t *tens, int dim1) {
  for (int i1=0; i1<dim1; i1++) {
      mpfr_init2(tens[i1], wp2);
  }
}


void init_rk1_mpc(mpc_t *tens, int dim1) {
  for (int i1=0; i1<dim1; i1++) {
      mpc_init3(tens[i1], wp2, wp2);
  }
}


void init_rk1_mpc_wp2(int wp2, mpc_t *tens, int dim1) {
  for (int i1=0; i1<dim1; i1++) {
      mpc_init3(tens[i1], wp2, wp2);
  }
}


void init_rk2_mpfr(mpfr_t **tens, int dim1, int dim2) {
  for (int i1=0; i1<dim1; i1++) {
    for (int i2=0; i2<dim2; i2++) {
      mpfr_init2(tens[i1][i2], wp2);
    }
  }
}


void init_rk2_mpc(mpc_t **tens, int dim1, int dim2) {
  for (int i1=0; i1<dim1; i1++) {
    for (int i2=0; i2<dim2; i2++) {
      mpc_init3(tens[i1][i2], wp2, wp2);
    }
  }
}


void init_rk3_mpc(mpc_t ***tens, int dim1, int dim2, int dim3) {
  for (int i1=0; i1<dim1; i1++) {
    for (int i2=0; i2<dim2; i2++) {
      for (int i3=0; i3<dim3; i3++) {
        mpc_init3(tens[i1][i2][i3], wp2, wp2);
      }
    }
  }
}


void init_rk4_mpc(mpc_t ****tens, int dim1, int dim2, int dim3, int dim4) {
  for (int i1=0; i1<dim1; i1++) {
    for (int i2=0; i2<dim2; i2++) {
      for (int i3=0; i3<dim3; i3++) {
        for (int i4=0; i4<dim4; i4++) {
          mpc_init3(tens[i1][i2][i3][i4], wp2, wp2);
        }
      }
    }
  }
}


//////
// CLEAR
//////
void mpq_rk1_clear(
  mpq_t *tens, int dim
) {
  for (int i=0; i<dim; i++) {
    mpq_clear(tens[i]);
  }
}

void mpc_rk1_clear(
  mpc_t *tens, int dim
) {
  for (int i=0; i<dim; i++) {
    mpc_clear(tens[i]);
  }
}


//////
// COPY
//////
void copy_rk1_int(int *copy, int *tens, int dim1) {
  for (int i1=0; i1<dim1; i1++) {
    copy[i1] = tens[i1];
  }
}

void copy_rk1_mpc(mpc_t *copy, mpc_t *tens, int dim1) {
  for (int i1=0; i1<dim1; i1++) {
    mpc_set(copy[i1], tens[i1], MPFR_RNDN);
  }
}

void copy_rk1_mpfr(mpfr_t *copy, mpfr_t *tens, int dim1) {
  for (int i1=0; i1<dim1; i1++) {
    mpfr_set(copy[i1], tens[i1], MPFR_RNDN);
  }
}


void copy_rk2_mpc(mpc_t **copy, mpc_t **tens, int dim1, int dim2) {
  for (int i1=0; i1<dim1; i1++) {
    for (int i2=0; i2<dim2; i2++) {
      mpc_set(copy[i1][i2], tens[i1][i2], MPFR_RNDN);
    }
  }
}


//////
// COMPARE
//////
void int_rk1_compare(int *tens1, int *tens2, int dim) {
  int found_diff = 0;
  for (int i=0; i<dim; i++) {
    // mpc_sub(tmp, tens1[i], tens2[i], MPFR_RNDN);
    // if (!mpc_lessthan_tol(tmp)) {
    if (tens1[i] != tens2[i]) {
      found_diff = 1;
      cout << "i = " << i << endl;
      cout << "1st: " << endl;
      cout << tens1[i] << endl;
      cout << "2nd: " << endl;
      cout << tens2[i] << endl;
    }
  }

  if (!found_diff) {
    cout << "CHECK: PASS" << endl;
    cout << endl;
  } else {
    cout << "CHECK: FAIL" << endl;
    cout << endl;
    exit(1);
  }
}

void int_rk1_compare_perm(
  int *perm, int *tens1, int *tens2, int dim
) {
  int found_diff = 0;
  for (int i=0; i<dim; i++) {
    // mpc_sub(tmp, tens1[i], tens2[i], MPFR_RNDN);
    // if (!mpc_lessthan_tol(tmp)) {
    if (perm[tens1[i]] != tens2[i]) {
      found_diff = 1;
      cout << "i = " << i << endl;
      cout << "1st: " << endl;
      cout << tens1[i] << endl;
      cout << "perm[1st] = " << endl;
      cout << perm[tens1[i]] << endl;
      cout << "2nd: " << endl;
      cout << tens2[i] << endl;
    }
  }

  if (!found_diff) {
    cout << "CHECK: PASS" << endl;
  } else {
    cout << "CHECK: FAIL" << endl;
    exit(1);
  }
}


void mpc_rk1_compare(mpc_t *tens1, mpc_t *tens2, int dim) {
  mpc_t tmp;
  mpc_init3(tmp, wp2, wp2);
  int found_diff = 0;
  for (int i=0; i<dim; i++) {
    // mpc_sub(tmp, tens1[i], tens2[i], MPFR_RNDN);
    // if (!mpc_lessthan_tol(tmp)) {
    if (!mpc_equal_within_tol(tens1[i], tens2[i])) {
      found_diff = 1;
      cout << "i = " << i << endl;
      cout << "1st: " << endl;
      print_mpc(&tens1[i]);
      cout << endl;
      cout << "2nd: " << endl;
      print_mpc(&tens2[i]);
      cout << endl;
      mpc_sub(tmp, tens1[i], tens2[i], MPFR_RNDN);
      cout << "difference:" << endl;
      print_mpc(&tmp);
      cout << endl;
    }
  }

  if (!found_diff) {
    cout << "CHECK: PASS" << endl;
    cout << endl;
  } else {
    cout << "CHECK: FAIL" << endl;
    cout << endl;
    exit(1);
  }
}


void mpc_rk1_compare_double(mpc_t *tens1, mpc_t *tens2, int dim) {
  mpc_t tmp;
  mpc_init3(tmp, wp2, wp2);
  int found_diff = 0;
  double d1r, d1i, d2r, d2i;
  double tol = 1e-10, cmpr, cmpi;
  for (int i=0; i<dim; i++) {
    d1r = mpfr_get_d(mpc_realref(tens1[i]), MPFR_RNDN);
    d1i = mpfr_get_d(mpc_imagref(tens1[i]), MPFR_RNDN);
    d2r = mpfr_get_d(mpc_realref(tens2[i]), MPFR_RNDN);
    d2i = mpfr_get_d(mpc_imagref(tens2[i]), MPFR_RNDN);

    cmpr = fabs(d1r - d2r);
    if (d1r != 0.) {cmpr /= fabs(d1r);}
    cmpi = fabs(d1i - d2i);
    if (d1i != 0.) {cmpi /= fabs(d1i);}

    if (cmpr > tol || cmpi > tol) {
      found_diff = 1;
      cout << "i = " << i << endl;
      cout << "1st: " << endl;
      print_mpc(&tens1[i]); cout << endl;
      printf(" double: (%.27g %.27g)\n", d1r, d1i);
      // cout << "  double: (" << d1r << " " << d1i <<  ")" << endl;
      cout << "2nd: " << endl;
      print_mpc(&tens2[i]); cout << endl;
      mpc_sub(tmp, tens1[i], tens2[i], MPFR_RNDN);
      printf(" double: (%.27g %.27g)\n", d2r, d2i);
      // cout << "  double: (" << d2r << " " << d2i <<  ")" << endl;
      cout << "difference:" << endl;
      printf(" double: (%.27g %.27g)\n", d1r-d2r, d1i-d2i);
      // cout << "  double: (" << d1r-d2r << " " << d2i-d2i <<  ")" << endl;
      print_mpc(&tmp); cout << endl;
    }
  }

  if (!found_diff) {
    cout << "CHECK: PASS" << endl;
    cout << endl;
  } else {
    cout << "CHECK: FAIL" << endl;
    cout << endl;
    exit(1);
  }
}


void mpc_rk1_compare_perm(
	// OUTPUT
	int *perm,
	// INPUT
	int dim, mpc_t *tens1, mpc_t *tens2
) {
	int found_perm, found_diff = 0;
	for (int i=0; i<dim; i++) {
    // cout << "i = " << i << endl;
    // cout << "looking for perm idx of:" << endl;
    // print_mpc(&tens1[i]); cout << endl;
		perm[i] = -1;
		for (int j=0; j<dim; j++) {
      // cout << "candidate j = " << j << ":" << endl;
      // print_mpc(&tens2[j]); cout << endl;
			if (mpc_equal_within_tol(tens1[i], tens2[j])) {
				found_perm = 1;
				perm[i] = j; 
				break;
			}
		}
		if (perm[i] == -1) {
			found_diff = 1;
			cout << "permutation index not found for i = " << i << endl;
			cout << "element:" << endl;
      print_mpc(&tens1[i]); cout << endl;
		}
	}

  if (!found_diff) {
    cout << "CHECK: PASS" << endl;
  } else {
    cout << "CHECK: FAIL" << endl;
    exit(1);
  }
}


int mpc_rk1_prop(
  // OUTPUT
  mpc_t *cf,
  // INPUT
  mpc_t *v1, mpc_t *v2, int dim,
  int *is_zero_v1, int *is_zero_v2
) {
  mpc_t tmpc; mpc_init3(tmpc, wp2, wp2);
  mpc_set_ui(*cf, 0, MPFR_RNDN);
  mpc_t *cfs = (mpc_t*) malloc(dim*sizeof(mpc_t));

  int is_zero_null;
  if (!is_zero_v1) {
    is_zero_null = 1;
    is_zero_v1 = (int*) malloc(dim*sizeof(int));
    is_zero_v2 = (int*) malloc(dim*sizeof(int));
    for (int i=0; i<dim; i++) {
      if (mpc_lessthan_tol(v1[i])) {
        is_zero_v1[i] = 1;
      } else {
        is_zero_v1[i] = 0;
      }
      if (mpc_lessthan_tol(v2[i])) {
        is_zero_v2[i] = 1;
      } else {
        is_zero_v2[i] = 0;
      }
    }
  } else {
    is_zero_null = 0;
    // cout << "is_zero_v1:" << endl;
    // for (int i=0; i<dim; i++) {
    //   cout << is_zero_v1[i] << endl;
    // }
    // cout << "is_zero_v2:" << endl;
    // for (int i=0; i<dim; i++) {
    //   cout << is_zero_v2[i] << endl;
    // }
  }

  int count_legit = 0;
  for (int i=0; i<dim; i++) {
    // cout << "i = " << i << endl;
    mpc_init3(cfs[i], wp2, wp2);

    // get ratio
    if (is_zero_v2[i]) {
      if (is_zero_v1[i]) {
        continue;
      } else {
        mpc_set_ui(cfs[count_legit], 0, MPFR_RNDN);
      }
    } else {
      if (is_zero_v1[i]) {
        // cout << "v2 non zero, v1 zero" << endl;
        return 0;
      } else {
        mpc_div(cfs[count_legit], v2[i], v1[i], MPFR_RNDN);
      }
    }
    // cout << "cfs = "; print_mpc(&cfs[count_legit]); cout << endl;

    // cumulate for avarage
    mpc_add(*cf, *cf, cfs[count_legit], MPFR_RNDN);

    if (count_legit == 0) {
      count_legit++;
      continue;
    }

    // check whether coeffs is compatible with previous ones
    if (!mpc_equal_within_tol(cfs[count_legit], cfs[count_legit-1])) {
			free(cfs);
      if (is_zero_null) {
        free(is_zero_v1);
        free(is_zero_v2);
      }
      return 0;
    }

    count_legit++;
  }

  // vector are proportional, return avarage coeff
  // printf("count_legit = %d\n", count_legit);
  mpc_div_ui(*cf, *cf, count_legit, MPFR_RNDN);
	free(cfs);
  if (is_zero_null) {
    free(is_zero_v1);
    free(is_zero_v2);
  }
  return 1;
}


void mpc_rk2_compare(mpc_t **tens1, mpc_t **tens2, int dim1, int dim2) {
  mpc_t tmp;
  mpc_init3(tmp, wp2, wp2);
  int found_diff = 0;
  for (int i1=0; i1<dim1; i1++) {
    for (int i2=0; i2<dim2; i2++) {
      if (!mpc_equal_within_tol(tens1[i1][i2], tens2[i1][i2])) {
        found_diff = 1;
        cout << "i1, i2 = " << i1 << ", " << i2 << endl;
        cout << "1st: "; print_mpc(&tens1[i1][i2]); cout << endl;
        cout << "2nd: "; print_mpc(&tens2[i1][i2]); cout << endl;
        mpc_sub(tmp, tens1[i1][i2], tens2[i1][i2], MPFR_RNDN);
        cout << "difference:"; print_mpc(&tmp); cout << endl;
      }
    }
  }

  if (!found_diff) {
    cout << "CHECK: PASS" << endl;
  } else {
    cout << "CHECK: FAIL" << endl;
    cout << endl;
    exit(1);
  }
}


void poly_frac_rk2_compare(
  struct poly_frac **tens1,
  struct poly_frac **tens2,
  int dim1, int dim2,
  mpc_t *roots
) {
  struct poly_frac tmp;
  poly_frac_build(&tmp);
  int found_diff = 0;
  for (int i1=0; i1<dim1; i1++) {
    for (int i2=0; i2<dim2; i2++) {
      if (!poly_frac_cmp_pf(&tens1[i1][i2], &tens2[i1][i2])) {
        found_diff = 1;
        cout << "i1, i2 = " << i1 << ", " << i2 << endl;
        cout << "1st: " << endl; poly_frac_print(&tens1[i1][i2]); cout << endl;
        cout << "2nd: " << endl; poly_frac_print(&tens2[i1][i2]); cout << endl;
        cout << "difference: " << endl;
        poly_frac_set_pf(&tmp, &tens2[i1][i2]);
        poly_frac_neg(&tmp);
        poly_frac_add_pf(&tmp, &tens1[i1][i2], &tmp, roots, NULL);
        poly_frac_print(&tmp);
      };
    }
  }

  if (!found_diff) {
    cout << "CHECK: PASS" << endl;
  } else {
    cout << "CHECK: FAIL" << endl;
    exit(1);
  }
}


//////
// APPEND
//////
void int_rk1_append_zero(
  // IN-OUT
  int *dim, int **vec,
  // INPUT
  int ap_dim
) {
  if (ap_dim == 0) {
    return;
  }

  int new_dim = *dim + ap_dim;

  // hold old values
  int *holder;
  if (*dim > 0) {
    holder = *vec;
  }
  
  // allocate
  // cout << "allocate" << endl;
  (*vec) = new int[new_dim];
  // copy old values
  // cout << "copy old values" << endl;
  for (int i=0; i<*dim; i++) {
    (*vec)[i] = holder[i];
  }
  // delete holder
  // cout << "delete holder" << endl;
  if (*dim > 0) {
    delete[] holder;
  }
  
  // copy new values
  // cout << "copy new values" << endl;
  for (int i=*dim; i<new_dim; i++) {
    (*vec)[i] = 0;
  }

  *dim = new_dim;

}


void int_rk1_append(
  // IN-OUT
  int *dim, int **vec,
  // INPUT
  int ap_dim, int *ap_vec
) {
  if (ap_dim == 0) {
    return;
  }

  int new_dim = *dim + ap_dim;

  // hold old values
  int *holder;
  if (*dim > 0) {
    holder = *vec;
  }
  
  // allocate
  // cout << "allocate" << endl;
  (*vec) = new int[new_dim];
  // copy old values
  // cout << "copy old values" << endl;
  for (int i=0; i<*dim; i++) {
    (*vec)[i] = holder[i];
  }
  // delete holder
  // cout << "delete holder" << endl;
  if (*dim > 0) {
    delete[] holder;
  }

  // copy new values
  // cout << "copy new values" << endl;
  for (int i=*dim; i<new_dim; i++) {
    (*vec)[i] = ap_vec[i-*dim];
  }
  *dim = new_dim;
}


void mpc_rk1_append(
  // IN-OUT
  int *dim, mpc_t **vec,
  // INPUT
  int ap_dim, mpc_t *ap_vec
) {
  if (ap_dim == 0) {
    return;
  }
  // hold old values
  mpc_t *holder;
  if (*dim > 0) {
    holder = *vec;
  }
  
  // allocate
  // cout << "allocate" << endl;
  (*vec) = new mpc_t[*dim+ap_dim];
  // copy old values
  // cout << "copy old values" << endl;
  for (int i=0; i<*dim; i++) {
    mpc_init3((*vec)[i], wp2, wp2);
    mpc_set((*vec)[i], holder[i], MPFR_RNDN);
  }
  // delete holder
  // cout << "delete holder" << endl;
  if (*dim > 0) {
    delete[] holder;
  }
  // copy new values
  // cout << "copy new values" << endl;
  for (int i=*dim; i<*dim+ap_dim; i++) {
    mpc_init3((*vec)[i], wp2, wp2);
    mpc_set((*vec)[i], ap_vec[i-*dim], MPFR_RNDN);
  }
}


void mpfr_rk1_append(
  // IN-OUT
  int *dim, mpfr_t **vec,
  // INPUT
  int ap_dim, mpfr_t *ap_vec
) {
  if (ap_dim == 0) {
    return;
  }
  // hold old values
  mpfr_t *holder;
  if (*dim > 0) {
    holder = *vec;
  }
  // allocate
  (*vec) = new mpfr_t[*dim+ap_dim];
  // copy old values
  for (int i=0; i<*dim; i++) {
    mpfr_init2((*vec)[i], wp2);
    mpfr_set((*vec)[i], holder[i], MPFR_RNDN);
  }
  // delete holder
  if (*dim > 0) {
    delete[] holder;
  }
  // copy new values
  for (int i=*dim; i<*dim+ap_dim; i++) {
    mpfr_init2((*vec)[i], wp2);
    mpfr_set((*vec)[i], ap_vec[i-*dim], MPFR_RNDN);
  }
}


//////
// TRANSFORM
//////
void mpfr_rk2f_switch_rows(
  // IN-OUT
  mpfr_t *tens,
  // INPUT
  int idx1, int idx2, int dim
) {
  mpfr_t tmpfr; mpfr_init2(tmpfr, wp2);
  for (int j=0; j<dim; j++) {
    mpfr_set(tmpfr, tens[idx1+j*dim], MPFR_RNDN);
    mpfr_set(tens[idx1+j*dim], tens[idx2+j*dim], MPFR_RNDN);
    mpfr_set(tens[idx2+j*dim], tmpfr, MPFR_RNDN);
  }
}


void mpfr_rk2f_switch_cols(
  // IN-OUT
  mpfr_t *tens,
  // INPUT
  int idx1, int idx2, int dim1, int dim2
) {
  mpfr_t tmpfr; mpfr_init2(tmpfr, wp2);
  for (int i=0; i<dim1; i++) {
    mpfr_set(tmpfr, tens[i+idx1*dim2], MPFR_RNDN);
    mpfr_set(tens[i+idx1*dim2], tens[i+idx2*dim2], MPFR_RNDN);
    mpfr_set(tens[i+idx2*dim2], tmpfr, MPFR_RNDN);
  }
}


void mpc_rk2_switch_rows(
  // IN-OUT
  mpc_t **tens,
  int idx1, int idx2, int dim
) {
  mpc_t tmpc; mpc_init3(tmpc, wp2, wp2);
  for (int j=0; j<dim; j++) {
    mpc_set(tmpc, tens[idx1][j], MPFR_RNDN);
    mpc_set(tens[idx1][j], tens[idx2][j], MPFR_RNDN);
    mpc_set(tens[idx2][j], tmpc, MPFR_RNDN);
  }
}


void mpc_rk2_switch_cols(
  // IN-OUT
  mpc_t **tens,
  int idx1, int idx2, int dim
) {
  mpc_t tmpc; mpc_init3(tmpc, wp2, wp2);
  for (int i=0; i<dim; i++) {
    mpc_set(tmpc, tens[i][idx1], MPFR_RNDN);
    mpc_set(tens[i][idx1], tens[i][idx2], MPFR_RNDN);
    mpc_set(tens[i][idx2], tmpc, MPFR_RNDN);
  }
}

//////
// TENSOR ALGEBRA
//////
void mpc_rk2_mul_mpc_rk1(
  // OUTPUT
  mpc_t *out_tens,
  // INPUT
  mpc_t **tens1, mpc_t *tens2,
  int dim1, int dim2
) {
  for (int i1=0; i1<dim1; i1++) {
    mpc_set_ui(out_tens[i1], 0, MPFR_RNDN);
    for (int i2=0; i2<dim2; i2++) {
      mpc_fma(
        out_tens[i1],
        tens1[i1][i2], tens2[i2], out_tens[i1],
        MPFR_RNDN
      );
    }
  }
}


void mpc_rk2_mul_mpc_rk2(
  // OUTPUT
  mpc_t **out_tens,
  // INPUT
  mpc_t **tens1, mpc_t **tens2,
  int dim1, int dim, int dim2
) {
  for (int i1=0; i1<dim1; i1++) {
    for (int i2=0; i2<dim2; i2++) {
      mpc_set_ui(out_tens[i1][i2], 0, MPFR_RNDN);
      for (int i=0; i<dim; i++) {
        mpc_fma(
          out_tens[i1][i2],
          tens1[i1][i], tens2[i][i2], out_tens[i1][i2],
          MPFR_RNDN
        );
      }
    }
  }
}


void mpc_rk2_mul_mpc_rk2_slice(
  // OUTPUT
  mpc_t *out_tens,
  // INPUT
  mpc_t **tens1, mpc_t **tens2,
  int dim1, int dim2, int slice
) {
  for (int i1=0; i1<dim1; i1++) {
    mpc_set_ui(out_tens[i1], 0, MPFR_RNDN);
    for (int i2=0; i2<dim2; i2++) {
      mpc_fma(
        out_tens[i1],
        tens1[i1][i2], tens2[i2][slice], out_tens[i1],
        MPFR_RNDN
      );
    }
  }
}


//////
// OPERATIONS
//////
void int_rk1_to_str(
  // OUTPUT
  char **out,
  // INPUT
  int dim, int *vec
) {
  if (*out) {
    free(*out);
  }
  (*out) = (char*) malloc(((MAX_POW_DIGITS+1)*dim+2+1)*sizeof(char));

  int c = 0, cc;
  (*out)[c++] = '[';

  char *value_str = (char*) malloc((MAX_POW_DIGITS + 1)*sizeof(char));
  for (int i=0; i<dim; i++) {
    // attach value
    snprintf(value_str, sizeof(value_str), "%d", vec[i]);
    cc = 0;
    while (value_str[cc] != '\0') {
      (*out)[c++] = value_str[cc++];
    }
    (*out)[c++] = ',';
  }
  (*out)[c-1] = ']';
  (*out)[c] = '\0';

}


int int_rk1_count_postivie(
  // INPUT
  int *vec, int dim
) {
  int count = 0;
  for (int i=0; i<dim; i++) {
    if (vec[i] > 0) {
      count++;
    }
  }
  return count;
}


int int_rk1_sum_postivie(
  // INPUT
  int *vec, int dim
) {
  int sum = 0;
  for (int i=0; i<dim; i++) {
    if (vec[i] > 0) {
      sum += vec[i];
    }
  }
  return sum;
}


int int_rk1_sum_negative_abs(
  // INPUT
  int *vec, int dim
) {
  int sum = 0;
  for (int i=0; i<dim; i++) {
    if (vec[i] < 0) {
      sum += -vec[i];
    }
  }
  return sum;
}


int int_rk1_is_positive_to_decimal(
  int *vec, int dim
) {
  int sec = 0, two = 1;
  for (int p=0; p<dim; p++) {
    if (vec[p] > 0) {
      sec += two;
    }
    two *= 2;
  }
  return sec;
}


int int_rk1_is_positive_to_decimal_mirror(
  int *vec, int dim
) {
  int sec = 0, two = 1;
  for (int p=dim-1; p>=0; p--) {
    if (vec[p] > 0) {
      sec += two;
    }
    two *= 2;
  }
  return sec;
}


int int_rk1_is_zero_to_decimal(
  int *vec, int dim
) {
  int sec = 0, two = 1;
  for (int p=0; p<dim; p++) {
    if (vec[p] == 0) {
      sec += two;
    }
    two *= 2;
  }
  return sec;
}

