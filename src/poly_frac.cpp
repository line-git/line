#include <iostream>
#include "mpc.h"

using namespace std;

#include "global_vars.h"
#include "utils.h"
#include "tensor_utils.h"
// #include "conversions.h"
#include "malloc_defs.h"
#include "poly_frac.h"
extern "C" {
  #include "cpoly.h"
  // #include "mp_roots_poly.h"
  #include "rel_err_mpc.h"
}


//////
// FUNCTIONS FOR POLYNOMIALS SAVED BY ROOTS
//////
void poly_eval_sym_zero_by_roots_wo_0(
  // OUTPUT
  mpc_t *value,
  // INPUT,
  mpc_t *roots, int *mults, int nroots
) {
  // #2BD: decide if it is worth using multiplicites and mpc_pow_ui
  //    pros: better than mul precision-wise;
  //    cons: need to store into temporary variable even when mults[k]
  //      equals one;
  // a similar argument can be made for the sign variable (swapping a bit
  // at every step vs considering if mults[k] is odd or even)
  mpc_set_ui(*value, 1, MPFR_RNDN);
  int sign = 1;
  for (int k=1; k<nroots; k++) {
    for (int m=0; m<mults[k]; m++) {
      // cout << "k, m = " << k << ", " << m << endl;
      mpc_mul(*value, *value, roots[k], MPFR_RNDN);
      sign = -sign;
    }
  }
  if (sign < 0) {
    mpc_neg(*value, *value, MPFR_RNDN);
  }
}


void poly_coeff1_by_roots_wo_0(
  // OUPUT
  mpc_t *value,
  // INPUT,
  mpc_t *roots, int *mults, int nroots
) {
  mpc_set_ui(*value, 0, MPFR_RNDN);
  mpc_t tmp0;
  mpc_init3(tmp0, wp2, wp2);
  mpc_t tmp1;
  mpc_init3(tmp1, wp2, wp2);
  for (int k=1; k<nroots; k++) {
    if (mults[k] == 0) {
      continue;
    }
    mpc_pow_ui(tmp0, roots[k], mults[k]-1, MPFR_RNDN);
    if (mults[k]%2 == 0) {
      mpc_neg(tmp0, tmp0, MPFR_RNDN);
    }
    for (int kp=1; kp<nroots; kp++) {
      if (kp == k || mults[kp] == 0) {
        continue;
      }
      mpc_pow_ui(tmp1, roots[kp], mults[kp], MPFR_RNDN);
      mpc_mul(tmp0, tmp0, tmp1, MPFR_RNDN);
      if (mults[k]%2 == 1) {
        mpc_neg(tmp0, tmp0, MPFR_RNDN);
      }
    }
    mpc_mul_ui(tmp0, tmp0, mults[k], MPFR_RNDN);
    mpc_add(*value, *value, tmp0, MPFR_RNDN);
  }
}


void poly_eval_by_roots(
  // OUTPUT
  mpc_t *out,
  // INPUT
  mpc_t *roots, int *mults, int nroots, mpc_t *val
) {
  mpc_t tmp;
  mpc_init3(tmp, wp2, wp2);
  mpc_set_ui(*out, 1, MPFR_RNDN);
  for (int k=0; k<nroots; k++) {
    if (mults[k] == 0) {
      continue;
    }
    mpc_sub(tmp, *val, roots[k], MPFR_RNDN);
    mpc_pow_ui(tmp, tmp, mults[k], MPFR_RNDN);
    mpc_mul(*out, *out, tmp, MPFR_RNDN);
  }
}


//////
// FUCNTIONS FOR POLYNOMIALS SAVED BY COEFFICIENTS
//////
void poly_eval(
  // OUTPUT
  mpc_t *out,
  // INPUT
  mpc_t *pol, int deg, mpc_t *val
) {
  mpc_set(*out, pol[deg], MPFR_RNDN);
  if (deg == 0) {
    return;
  }
  mpc_fma(*out, *val, pol[deg], pol[deg-1], MPFR_RNDN);
  for (int k=deg-1; k>0; k--) {
      mpc_fma(*out, *val, *out, pol[k-1], MPFR_RNDN);
  }
}


void mul_vpoly(
  mpc_t *out_pol, mpc_t *pol1, mpc_t *pol2,
  int deg1, int deg2,
  int vdeg1, int vdeg2,
  int ldeg1, int ldeg2
){
  /*
  Perform multiplication of two polynomials:

    pol1 = \sum_{i=ldeg1}^{deg1} a_i * x^i
    pol2 = \sum_{j=ldeg2}^{deg2} b_j * x^j

    pol1*pol2 = \sum_{k=ldeg1+ldeg2}^{deg1+deg2} \sum_{i=\max{0, k-deg2}}^{k} a_i * b_{k-i}
  */
  int i, k, i_min = 0, i_max, kp;
  // mpc_t tmp;
  // mpc_init3(tmp, wp2, wp2);
  
  int ldeg = ldeg1+ldeg2;
  for (k=0; k<=vdeg1+vdeg2; k++) {
    kp = k + ldeg;
    // cout << "k, kp = " << k << ", " << kp << endl;
    mpc_set_d(out_pol[k], 0, MPFR_RNDN);
    if (ldeg1 > kp - deg2) {
      i_min = ldeg1;
    } else {
      i_min = kp - deg2;
    }
    if (deg1 < kp - ldeg2) {
      i_max = deg1;
    } else {
      i_max = kp - ldeg2;
    }
    // cout << "i_min, i_max = " << i_min << ", " << i_max << endl;
    for (i=i_min; i<=i_max; i++) {
      // cout << "i = " << i << endl;
      mpc_fma(out_pol[k], pol1[i-ldeg1], pol2[kp-i-ldeg2], out_pol[k], MPFR_RNDN);
    }
  }
}


void poly_mul_root_new_out(
  mpc_t *out, mpc_t *coeffs, int deg, mpc_t root
) {
  mpc_set(out[deg+1], coeffs[deg], MPFR_RNDN);
  for (int k=deg; k>0; k--) {
    mpc_neg(out[k], coeffs[k], MPFR_RNDN);
    mpc_fma(out[k], out[k], root, coeffs[k-1], MPFR_RNDN);
  }
  mpc_mul(out[0], coeffs[0], root, MPFR_RNDN);
  mpc_neg(out[0], out[0], MPFR_RNDN);
}


//////
// POLY_FRAC STRUCTURE FUNDAMENTALS
//////

void poly_frac_build(
  struct poly_frac *pf
) {
  (*pf).coeffs = NULL;
  (*pf).mults = NULL;
}


void poly_frac_rk1_build(
  struct poly_frac *pf_tens,
  int dim1
) {
  for (int i1=0; i1<dim1; i1++) {
    pf_tens[i1].coeffs = NULL;
    pf_tens[i1].mults = NULL;
  }
}


void poly_frac_rk2_build(
  struct poly_frac **pf_tens,
  int dim1, int dim2
) {
  for (int i1=0; i1<dim1; i1++) {
    for (int i2=0; i2<dim2; i2++) {
      pf_tens[i1][i2].coeffs = NULL;
      pf_tens[i1][i2].mults = NULL;
    }
  }
}


void poly_frac_rk3_build(
  struct poly_frac ***pf_tens,
  int dim1, int dim2, int dim3
) {
  for (int i1=0; i1<dim1; i1++) {
    for (int i2=0; i2<dim2; i2++) {
      for (int i3=0; i3<dim3; i3++) {
        pf_tens[i1][i2][i3].coeffs = NULL;
        pf_tens[i1][i2][i3].mults = NULL;
      }
    }
  }
}


void poly_frac_free(
  struct poly_frac *pf
) {
  if (pf->mults) {
    delete[] pf->mults;
  }
  if (pf->coeffs) {
    for (int k=0; k<=pf->num_vdeg; k++) {
      mpc_clear(pf->coeffs[k]);
    }
    delete[] pf->coeffs;
  }
}


void poly_frac_rk1_free(
  struct poly_frac *pf_tens, int dim1
) {
  for (int i1=0; i1<dim1; i1++) {
    poly_frac_free(&pf_tens[i1]);
  }
}


void poly_frac_rk2_free(
  struct poly_frac **pf_tens, int dim1, int dim2
) {
  for (int i1=0; i1<dim1; i1++) {
    for (int i2=0; i2<dim2; i2++) {
      poly_frac_free(&pf_tens[i1][i2]);
    }
  }
}


void poly_frac_rk3_free(
  struct poly_frac ***pf_tens, int dim1, int dim2, int dim3
) {
  for (int i1=0; i1<dim1; i1++) {
    for (int i2=0; i2<dim2; i2++) {
      for (int i3=0; i3<dim3; i3++) {
        poly_frac_free(&pf_tens[i1][i2][i3]);
      }
    }
  }
}


void poly_frac_alloc_init_coeffs(
  struct poly_frac *pf
) {
  (*pf).coeffs = new mpc_t[(*pf).num_vdeg+1];
  init_rk1_mpc((*pf).coeffs, (*pf).num_vdeg+1);
}


void del_rk1_poly_frac(
  struct poly_frac *pftens,
  int dim  
) {
  for (int i=0; i<dim; i++) {
    if (pftens[i].coeffs) {
      delete[] pftens[i].coeffs;
    }
    if (pftens[i].mults) {
      delete[] pftens[i].mults;
    }
  }
  delete[] pftens;
}


void del_rk2_poly_frac(
  struct poly_frac **pftens,
  int dim1, int dim2
) {
  for (int i1=0; i1<dim1; i1++) {
    for (int i2=0; i2<dim2; i2++) {
      if (pftens[i1][i2].coeffs) {
        delete[] pftens[i1][i2].coeffs;
      }
      if (pftens[i1][i2].mults) {
        delete[] pftens[i1][i2].mults;
      }
    }
    delete[] pftens[i1];
  }
  delete[] pftens;
}


//////
// PRINT FUNCTIONS
//////
void poly_frac_print(
  struct poly_frac *pf
) {
  if ((*pf).num_vdeg == -1) {
    cout << "0" << endl;
    return;
  }

  cout << "num deg: " << (*pf).num_deg << " ("<< (*pf).num_vdeg << ")" << endl;
  if ((*pf).coeffs) {
    for (int k=0; k<=(*pf).num_vdeg; k++) {
      cout << "pow " << (*pf).num_deg - (*pf).num_vdeg + k << ": ";
      print_mpc(&(*pf).coeffs[k]);
      cout << endl;
    }
  }
  cout << "den deg: " << (*pf).den_deg << endl;
  if ((*pf).mults) {
    cout << "{";
    for (int k=0; k<(*pf).nroots; k++) {
      if ((*pf).mults[k] == 0) {
        continue;
      }
      cout << "r" << k << ": " << (*pf).mults[k] << ", ";
    }
    cout << "}" << endl;
  }
}


void poly_frac_rk2_print(
  struct poly_frac **pf_tens,
  int dim1, int dim2
) {
  for (int i1=0; i1<dim1; i1++) {
    for (int i2=0; i2<dim2; i2++) {
      cout << "(i1, i2) = " << i1 << ", " << i2 << endl;
      poly_frac_print(&pf_tens[i1][i2]);
    }
    cout << endl;
  }
}


void poly_frac_print_to_math(
  struct poly_frac *pf,
  mpc_t *roots
) {
  if ((*pf).num_vdeg == -1) {
    cout << "0";
    return;
  }
  string output = "(";

  char *mantissa_real = NULL, *mantissa_imag = NULL;
  string sign_real = "", sign_imag = "", id_mat = "", mpc_str;
  long exponent_real, exponent_imag;
  int bor, boi;  
  for (int k=0; k<=(*pf).num_vdeg; k++) {
    // cout << "k = " << k <<  endl;
    mantissa_real = mpfr_get_str(NULL, &exponent_real, 10, 0, mpc_realref((*pf).coeffs[k]), MPFR_RNDN);
    // cout << mantissa_real << endl;
    mantissa_imag = mpfr_get_str(NULL, &exponent_imag, 10, 0, mpc_imagref((*pf).coeffs[k]), MPFR_RNDN);
    // cout << mantissa_imag << endl;
    sign_real = "";
    bor = 0;
    if (mpfr_sgn(mpc_realref((*pf).coeffs[k])) <= 0) {
      sign_real = "-";
      bor = 1;
    }
    sign_imag = "";
    boi = 0;
    if (mpfr_sgn(mpc_imagref((*pf).coeffs[k])) <= 0) {
      sign_imag = "-";
      boi = 1;
    }
    // mpc_str = "("+sign_real+"0."+((string) (mantissa_real+bor))+"*^"+to_string(exponent_real)+\
    //   " + I*("+sign_imag+"0."+((string) (mantissa_imag+boi))+"*^"+to_string(exponent_imag)+"))";
    // if (mpc_lessthan_tol((*pf).coeffs[k])) {
    if (0 && mpc_lessthan_tol((*pf).coeffs[k])) {
      mpc_str="0";
    } else {
      mpc_str = "(";
      // if (!mpfr_lessthan_tol(mpc_realref(((*pf).coeffs[k])))) {
      if (1 || !mpfr_lessthan_tol(mpc_realref(((*pf).coeffs[k])))) {
        mpc_str += sign_real+"0."+((string) (mantissa_real+bor))+"*^"+to_string(exponent_real);
      }
      // if (!mpfr_lessthan_tol(mpc_imagref(((*pf).coeffs[k])))) {
      if (1 || !mpfr_lessthan_tol(mpc_imagref(((*pf).coeffs[k])))) {
        mpc_str += " + I*("+sign_imag+"0."+((string) (mantissa_imag+boi))+"*^"+to_string(exponent_imag)+")";
      }
      mpc_str+=")";
    }
    // cout << mpc_str << endl;
    output += "+"+mpc_str+"*eta^"+to_string(k + (*pf).num_deg - (*pf).num_vdeg);
    // cout << output << endl;
  }
  output += ")";
  if ((*pf).mults[0] > 0) {
    output += "/eta^"+to_string((*pf).mults[0]);
  }
  // cout << output << endl;
  for (int k=1; k<(*pf).nroots; k++) {
    if ((*pf).mults[k] == 0) {
      continue;
    }
    // cout << "k = " << k << endl;
    mantissa_real = mpfr_get_str(NULL, &exponent_real, 10, 0, mpc_realref(roots[k]), MPFR_RNDN);
    mantissa_imag = mpfr_get_str(NULL, &exponent_imag, 10, 0, mpc_imagref(roots[k]), MPFR_RNDN);
    sign_real = "";
    bor = 0;
    if (mpfr_sgn(mpc_realref(roots[k])) < 0) {
      sign_real = "-";
      bor = 1;
    }
    sign_imag = "";
    boi = 0;
    if (mpfr_sgn(mpc_imagref(roots[k])) < 0) {
      sign_imag = "-";
      boi = 1;
    }
    // mpc_str = "("+sign_real+"0."+((string) (mantissa_real+bor))+"*^"+to_string(exponent_real)+\
    //   " + I*("+sign_imag+"0."+((string) (mantissa_imag+boi))+"*^"+to_string(exponent_imag)+"))";
    // if (mpc_lessthan_tol(roots[k])) {
    if (0 && mpc_lessthan_tol(roots[k])) {
      mpc_str="0";
    } else {
      mpc_str = "(";
      // if (!mpfr_lessthan_tol(mpc_realref((roots[k])))) {
      if (1 || !mpfr_lessthan_tol(mpc_realref((roots[k])))) {
        mpc_str += sign_real+"0."+((string) (mantissa_real+bor))+"*^"+to_string(exponent_real);
      }
      // if (!mpfr_lessthan_tol(mpc_imagref((roots[k])))) {
      if (1 || !mpfr_lessthan_tol(mpc_imagref((roots[k])))) {
        mpc_str += " + I*("+sign_imag+"0."+((string) (mantissa_imag+boi))+"*^"+to_string(exponent_imag)+")";
      }
      mpc_str+=")";
    }
    output += "/(eta-"+mpc_str+")^"+to_string((*pf).mults[k]);
    // cout << output << endl;
  }
  cout << output;
}


void poly_frac_rk2_print_to_math(
  struct poly_frac **pf_tens,
  int dim1, int dim2,
  mpc_t *roots
) {
  cout << "{" << endl;
  for (int i1=0; i1<dim1; i1++) {
    cout << "{";
    for (int i2=0; i2<dim2; i2++) {
      poly_frac_print_to_math(&pf_tens[i1][i2], roots);
      if (i2 < dim2-1) {
        cout << ", ";
      }
    }
    cout << "}";
    if (i1 < dim1-1) {
      cout << ", ";
    }
    cout << endl;
  }
  cout << "}" << endl;
}


//////
// ASSIGNMENT FUNCTIONS
//////
void poly_frac_set_info(
  struct poly_frac *pf,
  int num_deg, int num_vdeg, int den_deg
) {
  (*pf).num_deg = num_deg;
  (*pf).num_vdeg = num_vdeg;
  (*pf).den_deg = den_deg;
}


void poly_frac_set_info_pf(
  struct poly_frac *pfout,
  struct poly_frac *pfin
) {
  (*pfout).num_deg = (*pfin).num_deg;
  (*pfout).num_vdeg = (*pfin).num_vdeg;
  (*pfout).den_deg = (*pfin).den_deg;
  (*pfout).nroots = (*pfin).nroots;
}


void poly_frac_set_mults_pf(
  struct poly_frac *pfout,
  struct poly_frac *pfin
) {
  if ((*pfin).num_vdeg == -1) {
    return;
  }
  
  if (!(*pfout).mults) {
    (*pfout).mults = new int[(*pfin).nroots];
  }
  for (int k=0; k<(*pfin).nroots; k++) {
    (*pfout).mults[k] = (*pfin).mults[k];
  }
}


void poly_frac_set_zero(
    struct poly_frac *pf
) {
  (*pf).num_deg = -1;
  (*pf).num_vdeg = -1;
  (*pf).den_deg = 0;
  (*pf).nroots = 0;
  if((*pf).coeffs) {
    // cout << "delete coeffs" << endl;
    delete[] (*pf).coeffs;
    // cout << "set coeffs to NULL" << endl;
    (*pf).coeffs = NULL;
  }
  if((*pf).mults) {
    // cout << "delete mults" << endl;
    delete[] (*pf).mults;
    // cout << "set mults to NULL" << endl;
    (*pf).mults = NULL;
  }
  // cout << "done" << endl;
}


void poly_frac_set_ui(
  struct poly_frac *pf,
  int ui_num,
  int nroots
) {
  if (ui_num == 0) {
    poly_frac_set_zero(pf);
    return;
  }
  (*pf).num_deg = 0;
  (*pf).num_vdeg = 0;
  (*pf).den_deg = 0;
  (*pf).nroots = nroots;
  if ((*pf).nroots > 0) {
    if (!(*pf).mults) {
      (*pf).mults = new int[nroots];
    }
    for (int k=0; k<nroots; k++) {
      (*pf).mults[k] = 0;
    }
  }
  if ((*pf).coeffs) {
    delete[] (*pf).coeffs;
  }
  (*pf).coeffs = new mpc_t[1];
  mpc_init3((*pf).coeffs[0], wp2, wp2);
  mpc_set_ui((*pf).coeffs[0], ui_num, MPFR_RNDN);
}


void poly_frac_set_si(
  struct poly_frac *pf,
  int si_num,
  int nroots
) {
  if (si_num >=0) {
    poly_frac_set_ui(pf, si_num, nroots);
  } else {
    poly_frac_set_ui(pf, -si_num, nroots);
    poly_frac_neg(pf);
  }
}


void poly_frac_set_ui_wp2(
  int wp2,
  struct poly_frac *pf,
  int ui_num,
  int nroots
) {
  if (ui_num == 0) {
    poly_frac_set_zero(pf);
    return;
  }
  (*pf).num_deg = 0;
  (*pf).num_vdeg = 0;
  (*pf).den_deg = 0;
  (*pf).nroots = nroots;
  if ((*pf).nroots > 0) {
    if (!(*pf).mults) {
      (*pf).mults = new int[nroots];
    }
    for (int k=0; k<nroots; k++) {
      (*pf).mults[k] = 0;
    }
  }
  if ((*pf).coeffs) {
    delete[] (*pf).coeffs;
  }
  (*pf).coeffs = new mpc_t[1];
  mpc_init3((*pf).coeffs[0], wp2, wp2);
  mpc_set_ui((*pf).coeffs[0], ui_num, MPFR_RNDN);
}


void poly_frac_set_mpc(
  // OUTPUT
  struct poly_frac *pf,
  // INPUT
  mpc_t *mpcin,
  int nroots
) {
  // if mpc input is zero, assign poly_frac
  if (mpc_lessthan_tol(*mpcin)) {
    poly_frac_set_zero(pf);
    // (*pf).num_vdeg = -1;
    // (*pf).coeffs = NULL;
    // (*pf).den_deg = 0;
    return;
  }

  if (nroots > 0) {
    if (!(*pf).mults) {
      (*pf).mults = new int[nroots];
    }
  }
  if ((*pf).coeffs) {
    mpc_rk1_clear((*pf).coeffs, (*pf).num_vdeg+1);
    delete[] (*pf).coeffs;
  }
  (*pf).num_deg = 0;
  (*pf).num_vdeg = 0;
  (*pf).den_deg = 0;
  (*pf).nroots = nroots;

  (*pf).coeffs = new mpc_t[1];
  mpc_init3((*pf).coeffs[0], wp2, wp2);
  mpc_set((*pf).coeffs[0], *mpcin, MPFR_RNDN);
  for (int k=0; k<nroots; k++) {
    (*pf).mults[k] = 0;
  }
}


void poly_frac_rk2_set_mpc_rk2(
  struct poly_frac **pf_tens,
  mpc_t **mpc_tens,
  int dim1, int dim2,
  int nroots
) {
  for (int i1=0; i1<dim1; i1++) {
    for (int i2=0; i2<dim2; i2++) {
      poly_frac_set_mpc(&pf_tens[i1][i2], &mpc_tens[i1][i2], nroots);
    }
  }  
}


void poly_frac_set_pf(
  struct poly_frac *pfout,
  struct poly_frac *pfin
) {
  // check if input is equal to output
  if (pfout == pfin) {
    return;
  }

  // check if input is zero
  if ((*pfin).num_vdeg == -1) {
    poly_frac_set_zero(pfout);
    return;
  }

  if ((*pfin).nroots > 0) {
    // cout << "alloc mults" << endl;
    if (!(*pfout).mults) {
      (*pfout).mults = new int[(*pfin).nroots];
    }
    // cout << "set mults" << endl;
    for (int k=0; k<(*pfin).nroots; k++) {
      (*pfout).mults[k] = (*pfin).mults[k];
    }
  }
  // cout << "alloc coeffs" << endl;
  if ((*pfout).coeffs) {
    mpc_rk1_clear((*pfout).coeffs, (*pfout).num_vdeg+1);
    delete[] (*pfout).coeffs;
  }
  (*pfout).coeffs = new mpc_t[(*pfin).num_vdeg+1];
  // cout << "set coeffs" << endl;
  for (int k=0; k<=(*pfin).num_vdeg; k++) {
    // cout << "k = " << k << endl;
    mpc_init3((*pfout).coeffs[k], wp2, wp2);
    // print_mpc(&(*pfin).coeffs[k]);
    // cout << endl;
    mpc_set((*pfout).coeffs[k], (*pfin).coeffs[k], MPFR_RNDN);
  }

  // cout << "set degrees" << endl;
  (*pfout).num_deg = (*pfin).num_deg;
  (*pfout).num_vdeg = (*pfin).num_vdeg;
  (*pfout).den_deg = (*pfin).den_deg;
  (*pfout).nroots = (*pfin).nroots;

  // cout << "done" << endl;
}


void poly_frac_rk2_set_pf_rk2(
  struct poly_frac **pf_tens_out,
  struct poly_frac **pf_tens_in,
  int dim1, int dim2
) {
  for (int i1=0; i1<dim1; i1++) {
    for (int i2=0; i2<dim2; i2++) {
      poly_frac_set_pf(&pf_tens_out[i1][i2], &pf_tens_in[i1][i2]);
    }
  }
}


void poly_frac_set_coeffs(
  struct poly_frac *pf,
  mpc_t *coeffs, int deg
) {
  // find lower degree
  int ldeg = 0;
  for (int k=0; k<=deg; k++) {
    if (mpc_lessthan_tol(coeffs[k])) {
      ldeg++;
    } else {
      break;
    }
  }
  if (ldeg == deg+1) {
    // input is zero
    poly_frac_set_zero(pf);
    return;
  }

  if ((*pf).coeffs) {
    mpc_rk1_clear((*pf).coeffs, (*pf).num_vdeg+1);
    delete[] (*pf).coeffs;
  }
  (*pf).num_deg = deg;
  (*pf).num_vdeg = deg - ldeg;
  (*pf).den_deg = 0;
  (*pf).coeffs = new mpc_t[(*pf).num_vdeg+1];
  for (int k=0; k<=(*pf).num_vdeg; k++) {
    mpc_init3((*pf).coeffs[k], wp2, wp2);
    mpc_set((*pf).coeffs[k], coeffs[ldeg+k], MPFR_RNDN);
  }
}


void poly_frac_set_coeffs_wp2(
  struct poly_frac *pf,
  mpc_t *coeffs, int deg,
  int wp2
) {
  // find lower degree
  int ldeg = 0;
  for (int k=0; k<=deg; k++) {
    if (mpc_lessthan_tol(coeffs[k])) {
      ldeg++;
    } else {
      break;
    }
  }
  if (ldeg == deg+1) {
    // input is zero
    poly_frac_set_zero(pf);
    return;
  }

  if ((*pf).coeffs) {
    delete[] (*pf).coeffs;
    mpc_rk1_clear((*pf).coeffs, (*pf).num_vdeg+1);
  }
  (*pf).num_deg = deg;
  (*pf).num_vdeg = deg - ldeg;
  (*pf).den_deg = 0;
  (*pf).coeffs = new mpc_t[(*pf).num_vdeg+1];
  for (int k=0; k<=(*pf).num_vdeg; k++) {
    mpc_init3((*pf).coeffs[k], wp2, wp2);
    mpc_set((*pf).coeffs[k], coeffs[ldeg+k], MPFR_RNDN);
  }
}


void poly_frac_rk2_set_id(
  struct poly_frac **pf_tens,
  int dim,
  int nroots
) {
  for (int i1=0; i1<dim; i1++) {
    poly_frac_set_ui(&pf_tens[i1][i1], 1, nroots);
    for (int i2=0; i2<i1; i2++) {
      poly_frac_set_zero(&pf_tens[i1][i2]);
    }
    for (int i2=i1+1; i2<dim; i2++) {
      poly_frac_set_zero(&pf_tens[i1][i2]);
    }
  }
}


//////
// FUNCTIONS FOR SYMBOLIC MANIPULATION
//////
int poly_frac_pow_behav(
  struct poly_frac *pf
) {
  // check if input is null
  if ((*pf).num_vdeg == -1) {
    return 0;
  }
  return (*pf).num_deg - (*pf).num_vdeg - (*pf).mults[0];
}


//////
// FUNCTIONS FOR SYMBOLIC MANIPULATION
//////
void poly_frac_normal(
  struct poly_frac *pf
) {
  int num_ldeg = (*pf).num_deg - (*pf).num_vdeg;
  if (num_ldeg == 0 || (*pf).mults[0] == 0) {
    return;
  }
  int pow_behav = num_ldeg - (*pf).mults[0];
  if (pow_behav < 0) {
    (*pf).num_deg = (*pf).num_vdeg;
    (*pf).den_deg -= num_ldeg;
    (*pf).mults[0] = -pow_behav;
  } else if (pow_behav > 0) {
    (*pf).num_deg -= (*pf).mults[0];
    (*pf).den_deg -= (*pf).mults[0];
    (*pf).mults[0] = 0;
  } else if (pow_behav == 0) {
    (*pf).num_deg = (*pf).num_vdeg;
    (*pf).den_deg -= (*pf).mults[0];
    (*pf).mults[0] = 0;
  }
}


void poly_frac_prune(
  struct poly_frac *pf
) {
  // #TBD: propagate error when comparing higher powers and tolerance
  // check if input is zero
  if ((*pf).num_vdeg == -1) {
    return;
  }

  // find real lower and upper degree
  int ldeg = 0, udeg = 0;
  for (int k=0; k<=(*pf).num_vdeg; k++) {
    if(mpc_lessthan_tol((*pf).coeffs[k])) {
      ldeg++;
    } else {
      break;
    }
  }

  for (int k=(*pf).num_vdeg; k>=0; k--) {
    if(mpc_lessthan_tol((*pf).coeffs[k])) {
      udeg++;
    } else {
      break;
    }
  }

  // prune
  if (ldeg == 0 && udeg == 0) {
    // nothing to be pruned
    return;
  } else if (ldeg > pf->num_vdeg) {
    // prune all coefficients
    poly_frac_set_zero(pf);
    return;
  } else {
    for (int k=ldeg; k<=pf->num_vdeg-udeg; k++) {
      mpc_set((*pf).coeffs[k-ldeg], (*pf).coeffs[k], MPFR_RNDD);
    }
    (*pf).num_vdeg -= ldeg + udeg;
    (*pf).num_deg -= udeg;
  }

}


void poly_frac_prune_idx(
  // IN-OUT
  struct poly_frac *pf,
  // INPUT
  int idx
) {
  if (pf->num_vdeg == -1) {
    return;
  }

  mpc_t *holder = pf->coeffs;
  pf->coeffs = new mpc_t[pf->num_vdeg+1-idx];
  int k;
  for (k=idx; k<=pf->num_vdeg; k++) {
    mpc_init3(pf->coeffs[k-idx], wp2, wp2);;
    mpc_set(pf->coeffs[k-idx], holder[k], MPFR_RNDN);
    mpc_clear(holder[k]);
  }
  for (k=0; k<idx; k++) {
    mpc_clear(holder[k]);
  }
  delete[] holder;
  pf->num_vdeg -= idx;

}


void poly_frac_prune_tol(
  struct poly_frac *pf, double wp2_rel_decr
) {
  // enlarge tolerance
  mpfr_t mpfr_tol_orig;
  mpfr_init2(mpfr_tol_orig, wp2);
  mpfr_set(mpfr_tol_orig, mpfr_tol, MPFR_RNDN);
  mpfr_tol_enlarge(wp2_rel_decr);
  // cout << "tol: "; mpfr_out_str(stdout, 10, 0, mpfr_tol, MPFR_RNDN); cout << endl;

  // prune
  poly_frac_prune(pf);

  // restore tolerance
  mpfr_set(mpfr_tol, mpfr_tol_orig, MPFR_RNDN);
}


int mpc_get_exp2(
  mpc_t in
) {
  if (mpfr_zero_p(mpc_realref(in))) {
    return mpfr_get_exp(mpc_imagref(in));
  } else if (mpfr_zero_p(mpc_imagref(in))) {
    return mpfr_get_exp(mpc_realref(in));
  } else {
    int exp2re =  mpfr_get_exp(mpc_realref(in));
    int exp2im =  mpfr_get_exp(mpc_imagref(in));
    if (exp2re > exp2im) {
      return exp2re;
    } else {
      return exp2im;
    }
  }
}


void poly_frac_prune_rel_tol_old(
  struct poly_frac *pf
) {
  if (pf->num_vdeg < 1) {
    poly_frac_prune(pf);
    return;
  }

  // multiply tolerance times 2^exp2, exp2 being the exponent in base 2
  // of the abs of the NLO coeff
  int exp2re, exp2im, exp2;
  exp2re = mpfr_get_exp(mpc_realref(pf->coeffs[1]));
  exp2im = mpfr_get_exp(mpc_imagref(pf->coeffs[1]));
  if (exp2re > exp2im) {
    exp2 = exp2re;
  } else {
    exp2 = exp2im;
  }
  for (int k=2; k<=pf->num_vdeg; k++) {
    exp2re = mpfr_get_exp(mpc_realref(pf->coeffs[k]));
    exp2im = mpfr_get_exp(mpc_imagref(pf->coeffs[k]));
    if (exp2re > exp2im) {
      if (exp2re > exp2) {
        exp2 = exp2re;
      }
    } else {
      if (exp2im > exp2) {
        exp2 = exp2im;
      }
    }

    // if (exp2 < mpfr_get_exp(mpc_realref(pf->coeffs[0]))) {
    //   if (exp2 < mpfr_get_exp(mpc_imagref(pf->coeffs[0]))) {
    //     // previous order was compatible with zero, use next one
    //     continue;
    //   } else {
    //     break;
    //   }
    // } else {
    //   break;
    // }

  }
  mpfr_mul_2si(mpfr_tol, mpfr_tol, exp2, MPFR_RNDN);
  // cout << "coeff = "; print_mpc(&pf->coeffs[k]); cout << endl;
  // cout << "exp2 = " << exp2 << endl;  
  // cout << "mpfr_tol = "; mpfr_out_str(stdout, 10, 0, mpfr_tol, MPFR_RNDN); cout << endl;

  // prune
  poly_frac_prune(pf);

  // restore tolerance
  mpfr_mul_2si(mpfr_tol, mpfr_tol, -exp2, MPFR_RNDN);
}


void poly_frac_prune_rel_tol_no_zero_p(
  struct poly_frac *pf, int wp_bin
) {
  if (pf->num_vdeg < 1) {
    poly_frac_prune(pf);
    return;
  }

  int tol_exp2 = mpfr_get_exp(mpfr_tol);
  int exp2 = mpc_get_exp2(pf->coeffs[0]);
  int exp2_next, exp2_diff, k, drop_idx = -1, drop_diff;
  for (k=1; k<=pf->num_vdeg; k++) {
    exp2_next = mpc_get_exp2(pf->coeffs[k]);
    // if (exp2_next > exp2) {
      exp2_diff = exp2_next - exp2;
    // } else {
    //   exp2_diff = exp2 - exp2_next;
    // }
    
    if (exp2_diff > wp_bin) {
      if (drop_idx == -1) {
        break;
      } else {
        drop_idx = -1;
      }
    }

    if (exp2_diff < 0) {
      if (-exp2_diff > wp_bin) {
        drop_idx = k;
      }
    }

    exp2 = exp2_next;
    continue;
    
  }

  if (k == pf->num_vdeg+1) {
    // no need to prune
    return;
  }
  
  printf("poly_frac_prune_rel_tol: pruning with idx = %d\n", k);
  poly_frac_print(pf);
  poly_frac_prune_idx(pf, k);
  poly_frac_normal(pf);

}


void poly_frac_prune_rel_tol(
  struct poly_frac *pf, int wp_bin
) {
  if (pf->num_vdeg < 1) {
    poly_frac_prune(pf);
    return;
  }

  int zero_p = mpc_zero_p(pf->coeffs[0]);
  int zero_p_next;
  int exp2 = mpc_get_exp2(pf->coeffs[0]);
  int exp2_next, exp2_diff, k, drop_idx = -1, drop_diff;
  for (k=1; k<=pf->num_vdeg; k++) {
    zero_p_next = mpc_zero_p(pf->coeffs[k]);    
    exp2_next = mpc_get_exp2(pf->coeffs[k]);
    
    if (zero_p && zero_p_next) {
      continue;
    }

    exp2_diff = exp2_next - exp2;    
    
    if (exp2_diff >= 0) {
      if (exp2_diff > wp_bin) {
        if (!zero_p_next) {
          if (drop_idx == -1) {
            break;
          } else {
            drop_idx = -1;
          }
        }
      } else {
        if (zero_p) {
          if (drop_idx == -1) {
            break;
          } else {
            drop_idx = -1;
          }
        } else if (zero_p_next) {
          drop_idx = k;
        }
      }
    } else {
      exp2_diff *= -1;
      if (exp2_diff > wp_bin) {
        if (!zero_p) {
          drop_idx = k;
        }
      } else {
        if (zero_p) {
          if (drop_idx == -1) {
            break;
          } else {
            drop_idx = -1;
          }          
        } else if (zero_p_next) {
          drop_idx = k;
        }
      }
    }

    zero_p = zero_p_next;
    exp2 = exp2_next;
    continue;
    
  }

  if (k == pf->num_vdeg+1) {
    // no need to prune
    return;
  }

  poly_frac_prune_idx(pf, k);
  poly_frac_normal(pf);

}


void poly_frac_trim_zero_p(
  struct poly_frac *pf
) {
  if (pf->num_vdeg == -1) {
    return;
  }

  int k;

  // find lower degree
  int ldeg = 0;
  for (k=0; k<=pf->num_vdeg; k++) {
    if(mpc_zero_p(pf->coeffs[k])) {
      ldeg++;
    } else {
      break;
    }
  }

  // find upper degree
  int udeg = 0;
  for (k=pf->num_vdeg; k>=0; k--) {
    if(mpc_zero_p(pf->coeffs[k])) {
      udeg++;
    } else {
      break;
    }
  }

  // prune
  if (ldeg == 0 && udeg == 0) {
    // nothing to be pruned
    return;
  } else if (ldeg > pf->num_vdeg) {
    // prune all coefficients
    poly_frac_set_zero(pf);
    return;
  }

  mpc_t *holder = pf->coeffs;
  pf->coeffs = new mpc_t[pf->num_vdeg-ldeg-udeg+1];
  for (k=ldeg; k<=pf->num_vdeg-udeg; k++) {
    mpc_init3(pf->coeffs[k-ldeg], wp2, wp2);;
    mpc_set(pf->coeffs[k-ldeg], holder[k], MPFR_RNDN);
    mpc_clear(holder[k]);
  }
  for (k=0; k<ldeg; k++) {
    mpc_clear(holder[k]);
  }
  for (k=pf->num_vdeg-udeg+1; k<=pf->num_vdeg; k++) {
    mpc_clear(holder[k]);
  }
  delete[] holder;
  pf->num_vdeg -= ldeg + udeg;
  pf->num_deg -= udeg;

  poly_frac_normal(pf);

}


void poly_frac_prune_rel_tol_real_max(
  struct poly_frac *pf, int wp_bin
) {
  if (pf->num_vdeg < 1) {
    poly_frac_prune(pf);
    return;
  }

  // cout << "PRUNE:" << endl;
  // poly_frac_print(pf);

  int k;

  // extract exponents and find maximum
  mp_exp_t *exps_re = new mp_exp_t[pf->num_vdeg+1];
  mp_exp_t *exps_im = new mp_exp_t[pf->num_vdeg+1];
  mp_exp_t exp_max;
  exps_re[0] = mpfr_get_exp(mpc_realref(pf->coeffs[0]));
  exps_im[0] = mpfr_get_exp(mpc_imagref(pf->coeffs[0]));
  exp_max = exps_re[0];
  if (exps_im[0] > exp_max) {
    exp_max = exps_im[0];
  }
  for (k=1; k<=pf->num_vdeg; k++) {
    exps_re[k] = mpfr_get_exp(mpc_realref(pf->coeffs[k]));
    exps_im[k] = mpfr_get_exp(mpc_imagref(pf->coeffs[k]));

    if (exps_re[k] > exp_max) {
      exp_max = exps_re[k];
    }
    if (exps_im[k] > exp_max) {
      exp_max = exps_im[k];
    }
  }

  // cout << "exp max = " << exp_max << endl;

  // PRUNE IN RELATIVE TO MAX
  for (k=0; k<=pf->num_vdeg; k++) {
    if (exp_max - exps_re[k] > wp_bin ) {
      // cout << "prune real of k = " << k << endl;
      mpfr_set_ui(mpc_realref(pf->coeffs[k]), 0, MPFR_RNDN);
    }
    if (exp_max - exps_im[k] > wp_bin ) {
      // cout << "prune imag of k = " << k << endl;
      mpfr_set_ui(mpc_imagref(pf->coeffs[k]), 0, MPFR_RNDN);
    }
  }

  // TRIM
  poly_frac_trim_zero_p(pf);

}


void poly_frac_rk2_trim_zero_p(
  struct poly_frac **pf,
  int dim1, int dim2
) {
  for (int i1=0; i1<dim1; i1++) {
    for (int i2=0; i2<dim2; i2++) {
      poly_frac_trim_zero_p(&pf[i1][i2]);
    }
  }
}


void poly_frac_rk2_prune_rel_tol_real_max(
  struct poly_frac **pf, int wp_bin,
  int dim1, int dim2
) {
  for (int i1=0; i1<dim1; i1++) {
    for (int i2=0; i2<dim2; i2++) {
      poly_frac_prune_rel_tol_real_max(&pf[i1][i2], wp_bin);
    }
  }
}


void poly_frac_rk1_prune_rel_tol(
  struct poly_frac *pf, int wp_bin,
  int dim1
) {
  // cout << "before prune:" << endl;
  // poly_frac_rk1_print(pf, dim1);
  for (int i1=0; i1<dim1; i1++) {
    poly_frac_prune_rel_tol(&pf[i1], wp_bin);
  }
  // cout << "after prune:" << endl;
  // poly_frac_rk1_print(pf, dim1);
}


void poly_frac_prune_radius_real(
  struct poly_frac *pf, int wp_bin, int rad_exp
) {
  if (pf->num_vdeg <= 0) {
    return;
  }
  
  int zero_p = mpfr_zero_p(mpc_realref(pf->coeffs[0]));
  int zero_p_next;
  int exp2 = mpfr_get_exp(mpc_realref(pf->coeffs[0]));
  int exp2_next, exp2_diff, k, drop_idx = -1, drop_diff;
  for (k=1; k<=pf->num_vdeg; k++) {
    zero_p_next = mpfr_zero_p(mpc_realref(pf->coeffs[k]));
    exp2_next = mpfr_get_exp(mpc_realref(pf->coeffs[k]));
    
    if (zero_p && zero_p_next) {
      continue;
    }
    
    exp2_diff = exp2_next + rad_exp - exp2;
    if (exp2_diff >= 0) {
      if (exp2_diff > wp_bin) {
        if (!zero_p_next) {
          if (drop_idx == -1) {
            break;
          } else {
            drop_idx = -1;
          }
        }
      } else {
        if (zero_p) {
          if (drop_idx == -1) {
            break;
          } else {
            drop_idx = -1;
          }
        } else if (zero_p_next) {
          drop_idx = k;
        }
      }
    } else {
      exp2_diff *= -1;
      if (exp2_diff > wp_bin) {
        if (!zero_p) {
          drop_idx = k;
        }
      } else {
        if (zero_p) {
          if (drop_idx == -1) {
            break;
          } else {
            drop_idx = -1;
          }          
        } else if (zero_p_next) {
          drop_idx = k;
        }
      }
    }

    zero_p = zero_p_next;
    exp2 = exp2_next;
    continue;
    
  }

  if (k == pf->num_vdeg+1) {
    // no need to prune
    return;
  }

  for (int kp=0; kp<k; kp++) {
    mpfr_set_ui(mpc_realref(pf->coeffs[kp]), 0,  MPFR_RNDN);
  }

}


void poly_frac_prune_radius_imag(
  struct poly_frac *pf, int wp_bin, int rad_exp
) {
  if (pf->num_vdeg <= 0) {
    return;
  }

  int zero_p = mpfr_zero_p(mpc_imagref(pf->coeffs[0]));
  int zero_p_next;
  int exp2 = mpfr_get_exp(mpc_imagref(pf->coeffs[0]));
  int exp2_next, exp2_diff, k, drop_idx = -1, drop_diff;
  for (k=1; k<=pf->num_vdeg; k++) {
    zero_p_next = mpfr_zero_p(mpc_imagref(pf->coeffs[k]));
    exp2_next = mpfr_get_exp(mpc_imagref(pf->coeffs[k]));
    
    if (zero_p && zero_p_next) {
      continue;
    }
    
    exp2_diff = exp2_next + rad_exp - exp2;
    if (exp2_diff >= 0) {
      if (exp2_diff > wp_bin) {
        if (!zero_p_next) {
          if (drop_idx == -1) {
            break;
          } else {
            drop_idx = -1;
          }
        }
      } else {
        if (zero_p) {
          if (drop_idx == -1) {
            break;
          } else {
            drop_idx = -1;
          }
        } else if (zero_p_next) {
          drop_idx = k;
        }
      }
    } else {
      exp2_diff *= -1;
      if (exp2_diff > wp_bin) {
        if (!zero_p) {
          drop_idx = k;
        }
      } else {
        if (zero_p) {
          if (drop_idx == -1) {
            break;
          } else {
            drop_idx = -1;
          }          
        } else if (zero_p_next) {
          drop_idx = k;
        }
      }
    }

    zero_p = zero_p_next;
    exp2 = exp2_next;
    continue;
    
  }

  if (k == pf->num_vdeg+1) {
    // no need to prune
    return;
  }

  for (int kp=0; kp<k; kp++) {
    mpfr_set_ui(mpc_imagref(pf->coeffs[kp]), 0,  MPFR_RNDN);
  }

}


void poly_frac_prune_radius_real_imag(
  struct poly_frac *pf, int wp_bin, int rad_exp, int start_real
) {
  if (pf->num_vdeg <= 0) {
    return;
  }
  
  int zero_p, exp2, real = start_real;
  if (real) {
    zero_p = mpfr_zero_p(mpc_realref(pf->coeffs[0]));
    exp2 = mpfr_get_exp(mpc_realref(pf->coeffs[0]));
    real = 0;
  } else {
    zero_p = mpfr_zero_p(mpc_imagref(pf->coeffs[0]));
    exp2 = mpfr_get_exp(mpc_imagref(pf->coeffs[0]));
    real = 1;
  }
  int zero_p_next;
  int exp2_next, exp2_diff, k, drop_idx = -1, drop_diff;
  for (k=1; k<=pf->num_vdeg; k++) {
    if (real) {
      zero_p_next = mpfr_zero_p(mpc_realref(pf->coeffs[k]));
      exp2_next = mpfr_get_exp(mpc_realref(pf->coeffs[k]));
      real = 0;
    } else {
      zero_p_next = mpfr_zero_p(mpc_imagref(pf->coeffs[k]));
      exp2_next = mpfr_get_exp(mpc_imagref(pf->coeffs[k]));
      real = 1;
    }
    
    
    if (zero_p && zero_p_next) {
      continue;
    }
    
    exp2_diff = exp2_next + rad_exp - exp2;
    if (exp2_diff >= 0) {
      if (exp2_diff > wp_bin) {
        if (!zero_p_next) {
          if (drop_idx == -1) {
            break;
          } else {
            drop_idx = -1;
          }
        }
      } else {
        if (zero_p) {
          if (drop_idx == -1) {
            break;
          } else {
            drop_idx = -1;
          }
        } else if (zero_p_next) {
          drop_idx = k;
        }
      }
    } else {
      exp2_diff *= -1;
      if (exp2_diff > wp_bin) {
        if (!zero_p) {
          drop_idx = k;
        }
      } else {
        if (zero_p) {
          if (drop_idx == -1) {
            break;
          } else {
            drop_idx = -1;
          }          
        } else if (zero_p_next) {
          drop_idx = k;
        }
      }
    }

    zero_p = zero_p_next;
    exp2 = exp2_next;
    continue;
    
  }

  if (k == pf->num_vdeg+1) {
    // no need to prune
    return;
  }

  real = start_real;
  for (int kp=0; kp<k; kp++) {
    if (real) {
      mpfr_set_ui(mpc_realref(pf->coeffs[kp]), 0,  MPFR_RNDN);
      real = 0;
    } else {
      mpfr_set_ui(mpc_imagref(pf->coeffs[kp]), 0,  MPFR_RNDN);
      real = 1;
    }
  }

}


void poly_frac_prune_radius(
  struct poly_frac *pf, int wp_bin, int rad_exp
) {
  // cout << "before prune radius:" << endl;
  // poly_frac_print(pf);

  // prune real part
  poly_frac_prune_radius_real(pf, wp_bin, rad_exp);
  // cout << "after prune radius real:" << endl;
  // poly_frac_print(pf);
  // prune imag part
  poly_frac_prune_radius_imag(pf, wp_bin, rad_exp);
  // cout << "after prune radius imag:" << endl;
  // poly_frac_print(pf);

  // prune alternate
  poly_frac_prune_radius_real_imag(pf, wp_bin, rad_exp, 1);
  poly_frac_prune_radius_real_imag(pf, wp_bin, rad_exp, 0);

  // trim
  poly_frac_trim_zero_p(pf);

  // cout << "after trim:" << endl;
  // poly_frac_print(pf);
  
}


void poly_frac_rk1_prune_radius(
  struct poly_frac *pf, int wp_bin, int rad_exp,
  int dim1
) {
  for (int i1=0; i1<dim1; i1++) {
    // cout << "i1 = " << i1 << endl;
    poly_frac_prune_radius(&pf[i1], wp_bin, rad_exp);
  }
}


void poly_frac_rk2_prune_radius(
  struct poly_frac **pf, int wp_bin, int rad_exp,
  int dim1, int dim2
) {
  for (int i1=0; i1<dim1; i1++) {
    for (int i2=0; i2<dim2; i2++) {
      poly_frac_prune_radius(&pf[i1][i2], wp_bin, rad_exp);
    }
  }
}


void poly_frac_rk2_normal(
  struct poly_frac **pf,
  int dim1, int dim2
) {
  for (int i1=0; i1<dim1; i1++) {
    for (int i2=0; i2<dim2; i2++) {
      poly_frac_normal(&pf[i1][i2]);
    }
  }
}


void poly_frac_rk2_prune(
  struct poly_frac **pf,
  int dim1, int dim2, double wp2_rel_decr
) {
  // enlarge tolerance
  mpfr_t mpfr_tol_orig;
  mpfr_init2(mpfr_tol_orig, wp2);
  mpfr_set(mpfr_tol_orig, mpfr_tol, MPFR_RNDN);
  // int wp2_orig = wp2;
  mpfr_tol_enlarge(wp2_rel_decr);

  // cout << endl; cout << "pruning rk2 with:" << endl;
  // cout << "wp2_rel_decr_prune = " << wp2_rel_decr << endl;
  // cout << "wp2 = " << (int) (wp2_rel_decr * wp2) << endl;
  // cout << "tol: "; mpfr_out_str(stdout, 10, 0, mpfr_tol, MPFR_RNDN); cout << endl;

  // prune elements
  for (int i1=0; i1<dim1; i1++) {
    for (int i2=0; i2<dim2; i2++) {
      poly_frac_prune(&pf[i1][i2]);
    }
  }

  // restore tolerance
  mpfr_set(mpfr_tol, mpfr_tol_orig, MPFR_RNDN);
  // wp2 = wp2_orig;

}


void poly_frac_rk2_prune_rel_tol(
  struct poly_frac **pf, int wp_bin,
  int dim1, int dim2
) {
  // cout << "before prune:" << endl;
  // poly_frac_rk2_print(pf, dim1, dim2);
  for (int i1=0; i1<dim1; i1++) {
    for (int i2=0; i2<dim2; i2++) {
      poly_frac_prune_rel_tol(&pf[i1][i2], wp_bin);
    }
  }
  // cout << "after prune:" << endl;
  // poly_frac_rk2_print(pf, dim1, dim2);
}


void poly_frac_shift(
  struct poly_frac *pfout,
  struct poly_frac *pfin,
  mpc_t *shift,
  int swap_root_lab
) {
  if ((*pfin).num_vdeg == -1) {
    poly_frac_set_zero(pfout);
    return;
  }

  if (mpc_lessthan_tol(*shift)) {
    if (pfout != pfin) {
      poly_frac_set_pf(pfout, pfin);
    } else {
      return;
    }
  }

  // SHIFT COEFFICIENTS
  if ((*pfout).coeffs) {
    mpc_rk1_clear((*pfout).coeffs, (*pfout).num_vdeg+1);
    delete[] (*pfout).coeffs;
  }
  (*pfout).coeffs = new mpc_t[(*pfin).num_vdeg+1];
  init_rk1_mpc((*pfout).coeffs, (*pfin).num_vdeg+1);
  for (int k=0; k<=(*pfin).num_vdeg; k++) {
    mpc_set((*pfout).coeffs[k], (*pfin).coeffs[k], MPFR_RNDN);
  }
  poly_shift(*shift, (*pfout).coeffs, (*pfin).num_vdeg);

  // poly_frac_set_pf(pfout, pfin);
  poly_frac_set_info_pf(pfout, pfin);
  poly_frac_set_mults_pf(pfout, pfin);

  // SWAP MULTIPLICITIES
  (*pfout).mults[0] = (*pfin).mults[swap_root_lab];
  (*pfout).mults[swap_root_lab] = (*pfin).mults[0];

  // SHIFT eta^ldeg
  mpc_t root;
  mpc_init3(root, wp2, wp2);
  mpc_neg(root, *shift, MPFR_RNDN);
  int ldeg = (*pfin).num_deg-(*pfin).num_vdeg;
  for (int k=0; k<ldeg; k++) {
    poly_mul_root((*pfout).coeffs, (*pfout).num_vdeg+k, root);
  }
  (*pfout).num_vdeg = (*pfout).num_deg;

  // cout << "before pruning:" << endl;
  // poly_frac_print(pfout);
  // poly_frac_prune(pfout);
  // cout << "after pruning:" << endl;
  // poly_frac_print(pfout);
  // cout << "after normal:" << endl;
  // poly_frac_normal(pfout);
  // poly_frac_print(pfout);
}


void poly_frac_rk2_shift(
  struct poly_frac **pfout,
  struct poly_frac **pfin,
  mpc_t *shift,
  int swap_root_lab,
  int dim1, int dim2
) {
  for (int i1=0; i1<dim1; i1++) {
    for (int i2=0; i2<dim2; i2++) {
      poly_frac_shift(&pfout[i1][i2], &pfin[i1][i2], shift, swap_root_lab);
    }
  }
}


void poly_frac_perm_roots(
  // IN-OUT
  struct poly_frac *pf,
  // INPUT
  int *perm
) {
  // check if input is zero
  if ((*pf).num_vdeg == -1) {
    return;
  }

  if ((*pf).mults) {
    int *mults_unperm = new int[(*pf).nroots];
    for (int k=0; k<(*pf).nroots; k++) {
      mults_unperm[k] = (*pf).mults[k];
    }
    for (int k=0; k<(*pf).nroots; k++) {
      (*pf).mults[perm[k]] = mults_unperm[k];
    }
    delete[] mults_unperm;
  }
}


void poly_frac_rk2_perm_roots(
  // IN-OUT
  struct poly_frac **pf,
  // INPUT
  int *perm, int dim1, int dim2
) {
  for (int i1=0; i1<dim1; i1++) {
    for (int i2=0; i2<dim2; i2++) {
      poly_frac_perm_roots(&pf[i1][i2], perm);
    }
  }
}


void poly_frac_arg_inv(
  // OUTPUT
  struct poly_frac *pfout,
  // INPUT
  struct poly_frac *pfin, mpc_t *roots
) {
  // cout << "pf in:" << endl;
  // poly_frac_print(pfin);
  // cout << "nroots = " << (*pfin).nroots << endl;

  if ((*pfin).num_vdeg == -1) {
    poly_frac_set_zero(pfout);
    return;
  }
  
  // MULTIPLICITIES (BESIDES NULL ROOT)    
  if ((*pfout).mults) {
    if ((*pfout).nroots != (*pfin).nroots) {
      delete[] (*pfout).mults;
      (*pfout).mults = new int[(*pfin).nroots];
    }
  } else {
    (*pfout).mults = new int[(*pfin).nroots];
  }
  for (int k=1; k<(*pfin).nroots; k++) {
    // cout << "k = " << k << endl;
    // cout << "pfin.mult = " << (*pfin).mults[k] << endl;
    (*pfout).mults[k] = (*pfin).mults[k];
  }
  (*pfout).nroots = (*pfin).nroots;

  // DECIDE DEGREES AND NULL ROOT MULTIPLICITY
  if ((*pfout).coeffs) {
    if ((*pfout).num_vdeg != (*pfin).num_vdeg) {
      mpc_rk1_clear((*pfout).coeffs, (*pfout).num_vdeg+1);
      delete[] (*pfout).coeffs;
      (*pfout).coeffs = new mpc_t[(*pfin).num_vdeg+1];
      init_rk1_mpc((*pfout).coeffs, (*pfin).num_vdeg+1);
    }
  } else {
    (*pfout).coeffs = new mpc_t[(*pfin).num_vdeg+1];
    init_rk1_mpc((*pfout).coeffs, (*pfin).num_vdeg+1);
  }
  (*pfout).num_vdeg = (*pfin).num_vdeg;

  int pow_behav = (*pfin).den_deg - (*pfin).num_deg;  
  // cout << "pow behav = " << pow_behav << endl;
  if (pow_behav >= 0) {
    (*pfout).num_deg = (*pfin).num_vdeg + pow_behav;
    (*pfout).mults[0] = 0;
    (*pfout).den_deg = (*pfin).den_deg - (*pfin).mults[0];
  } else if (pow_behav < 0) {
    (*pfout).num_deg = (*pfin).num_vdeg;
    (*pfout).mults[0] = -pow_behav;
    (*pfout).den_deg = (*pfin).den_deg - (*pfin).mults[0] - pow_behav;
  }

  // NORMALIZATION FACTOR
  mpc_t norm;
  mpc_init3(norm, wp2, wp2);
  mpc_set_ui(norm, 1, MPFR_RNDN);
  for (int k=1; k<(*pfin).nroots; k++) {
    for (int m=0; m<(*pfin).mults[k]; m++) {      
      mpc_mul(
        norm, norm, roots[k], MPFR_RNDN  // , wp_bin
      );
      // cout << "root = "; print_mpc(&roots[k]); cout << endl;
      // cout << "norm = "; print_mpc(&norm); cout << endl;
    }
  }
  if (((*pfin).den_deg - (*pfin).mults[0])%2 != 0) {
    mpc_neg(norm, norm, MPFR_RNDN);
  }

  // NUMERATOR (revert order and normalize)

  for (int k=0; k<=(*pfout).num_vdeg; k++) {
    // cout << "k = " << k << endl;
    // cout << "coeff = "; print_mpc(&(*pfin).coeffs[(*pfout).num_vdeg-k]); cout << endl;
    // cout << "norm = "; print_mpc(&norm); cout << endl;
    mpc_div((*pfout).coeffs[k], (*pfin).coeffs[(*pfout).num_vdeg-k], norm, MPFR_RNDN);
    // cout << "division = "; print_mpc(&(*pfout).coeffs[k]); cout << endl;
  }
  
}


void poly_frac_rk2_infty(
  // OUTPUT
  struct poly_frac **pfout, mpc_t *out_roots, int *out_zero_label,
  // INPUT
  struct poly_frac **pfin, mpc_t *in_roots,
  int dim1, int dim2, int nroots
) {
  // invert roots
  mpc_set_ui(out_roots[0], 0, MPFR_RNDN);
  for (int k=1; k<nroots; k++) {
    mpc_pow_si(out_roots[k], in_roots[k], -1, MPFR_RNDN);
  }

  // invert matrix elements and set new zero_label
  // *out_zero_label = -1;
  for (int i1=0; i1<dim1; i1++) {
    for (int i2=0; i2<dim2; i2++) {
      // cout << "i1, i2 = " << i1 << ", " << i2 << endl;
      poly_frac_arg_inv(&pfout[i1][i2], &pfin[i1][i2], in_roots);
      // cout << "pfout:" << endl;
      // poly_frac_print(&pfout[i1][i2]);
      // if (pfout[i1][i2].mults) {
      //   if (pfout[i1][i2].mults[0] != 0) {
      //     *out_zero_label = 0;
      //   }
      // }
      // mpc_t val, eta_val, inv_eta_val;
      // mpc_init3(val, wp2, wp2);
      // mpc_init3(eta_val, wp2, wp2);
      // mpc_init3(inv_eta_val, wp2, wp2);
      // mpc_set_ui(eta_val, 17, MPFR_RNDN);
      // mpc_pow_si(inv_eta_val, eta_val, -1, MPFR_RNDN);
      // cout << "eval pfin in eta = 17:" << endl;
      // poly_frac_eval_value(&val, &pfin[i1][i2], in_roots, &eta_val);
      // print_mpc(&val); cout << endl;
      // cout << "eval pfout in 1/eta = 1/17:" << endl;
      // poly_frac_eval_value(&val, &pfout[i1][i2], out_roots, &inv_eta_val);
      // print_mpc(&val); cout << endl;

      // chain-rule factor
      if (pfout[i1][i2].num_vdeg != -1) {
        poly_frac_neg(&pfout[i1][i2]);
        poly_frac_mul_sym_pow(&pfout[i1][i2], &pfout[i1][i2], -2);
      }
    }
  }

  *out_zero_label = 0;

}


//////
// BASIC ARITHMETIC FUNCTIONS
//////
void poly_frac_neg(
  struct poly_frac *pf
) {
  for (int k=0; k<=(*pf).num_vdeg; k++) {
    mpc_neg((*pf).coeffs[k], (*pf).coeffs[k], MPFR_RNDN);
  }
}


void poly_frac_rk2_neg(
  struct poly_frac **pfmat, int dim1, int dim2
) {
  for (int i1=0; i1<dim1; i1++) {
    for (int i2=0; i2<dim2; i2++) {
      poly_frac_neg(&pfmat[i1][i2]);
    }
  }
}


void poly_frac_mul_ui(
  struct poly_frac *pfout,
  struct poly_frac *pfin,
  int ui
) {
  // if input polynomial fraction is zero, set output to zero
  if ((*pfin).num_vdeg == -1) {
    if (pfout != pfin) {
      poly_frac_set_zero(pfout);
    }
    return;
  }

  // if input mpc is zero, set output to zero
  if (ui == 0) {
    poly_frac_set_zero(pfout);
    return;
  }

  if (pfout != pfin) {
    poly_frac_set_info_pf(pfout, pfin);
    if (!(*pfout).mults) {
      (*pfout).mults = new int[(*pfin).nroots];
    }
    if ((*pfout).coeffs) {
      mpc_rk1_clear((*pfout).coeffs, (*pfout).num_vdeg+1);
      delete[] (*pfout).coeffs;
    }
    poly_frac_alloc_init_coeffs(pfout);
    for (int k=0; k<(*pfin).nroots; k++) {
      (*pfout).mults[k] = (*pfin).mults[k];
    }
  }

  for (int k=0; k<=(*pfin).num_vdeg; k++) {
    mpc_mul_ui((*pfout).coeffs[k], (*pfin).coeffs[k], ui, MPFR_RNDN);
  }
}


void poly_frac_mul_si(
  struct poly_frac *pfout,
  struct poly_frac *pfin,
  int si
) {
  if (si>=0) {
    poly_frac_mul_ui(pfout, pfin, si);
  } else {
    poly_frac_mul_ui(pfout, pfin, -si);
    poly_frac_neg(pfout);
  }
}


void poly_frac_mul_mpc(
  struct poly_frac *pfout,
  struct poly_frac *pfin,
  mpc_t *mpcin
) {
  // if input polynomial fraction is zero, set output to zero
  if ((*pfin).num_vdeg == -1) {
    if (pfout != pfin) {
      poly_frac_set_zero(pfout);
    }
    return;
  }

  // if input mpc is zero, set output to zero
  if (mpc_lessthan_tol(*mpcin)) {
    // cout << "about to set zero" << endl;
    poly_frac_set_zero(pfout);
    // if ((*pfout).coeffs) {
    //   delete[] (*pfout).coeffs;
    //   (*pfout).coeffs = NULL;
    // }
    // (*pfout).num_vdeg = -1;
    // for (int k=0; k<nroots; k++) {
    //   (*pfout).mults[k] = 0;
    // }
    return;
  }

  if (pfout != pfin) {
    if (!(*pfout).mults) {
      (*pfout).mults = new int[(*pfin).nroots];
    }
    if ((*pfout).coeffs) {
      mpc_rk1_clear((*pfout).coeffs, (*pfout).num_vdeg+1);
      delete[] (*pfout).coeffs;
    }
    poly_frac_set_info_pf(pfout, pfin);
    poly_frac_alloc_init_coeffs(pfout);
    for (int k=0; k<(*pfin).nroots; k++) {
      (*pfout).mults[k] = (*pfin).mults[k];
    }
  }

  for (int k=0; k<=(*pfin).num_vdeg; k++) {
    mpc_mul((*pfout).coeffs[k], (*pfin).coeffs[k], *mpcin, MPFR_RNDN);
  }
}


void poly_frac_div_mpc(
  struct poly_frac *pfout,
  struct poly_frac *pfin,
  mpc_t *mpcin
) {
  // if input polynomial fraction is zero, set output to zero
  if ((*pfin).num_vdeg == -1) {
    if (pfout != pfin) {
      poly_frac_set_zero(pfout);
    }
    return;
  }

  // if input mpc is zero, give warning
  if (mpc_lessthan_tol(*mpcin)) {
    cout << endl; cout << "WARNING: dividing by quantity lower than tolerance" << endl; cout << endl;
  }

  if (pfout != pfin) {
    poly_frac_set_info_pf(pfout, pfin);
    if (!(*pfout).mults) {
      (*pfout).mults = new int[(*pfin).nroots];
    }
    if ((*pfout).coeffs) {
      mpc_rk1_clear((*pfout).coeffs, (*pfout).num_vdeg+1);
      delete[] (*pfout).coeffs;
    }
    poly_frac_alloc_init_coeffs(pfout);
    for (int k=0; k<(*pfin).nroots; k++) {
      (*pfout).mults[k] = (*pfin).mults[k];
    }
  }

  for (int k=0; k<=(*pfin).num_vdeg; k++) {
    mpc_div((*pfout).coeffs[k], (*pfin).coeffs[k], *mpcin, MPFR_RNDN);
  }
}


void poly_frac_mul_pf(
  struct poly_frac *pfout,
  struct poly_frac *pfin1,
  struct poly_frac *pfin2
) {
  // if input is zero, set output to zero
  if ((*pfin1).num_vdeg == -1) {
    if (pfout != pfin1) {
      poly_frac_set_zero(pfout);
    }
    return;
  }
  if ((*pfin2).num_vdeg == -1) {
    if (pfout != pfin2) {
      poly_frac_set_zero(pfout);
    }
    return;
  }

  // deal with overwriting
  // cout << "deal with overwriting" << endl;
  struct poly_frac *pf1, *pf2;
  int ovrwrt;
  if (pfout == pfin1 && pfout != pfin2) {
    ovrwrt = 1;
    pf1 = new struct poly_frac[1];
    poly_frac_build(pf1);
    poly_frac_set_pf(pf1, pfin1);
    pf2 = pfin2;
  } else if (pfout != pfin1 && pfout == pfin2) {
    ovrwrt = 2;
    pf2 = new struct poly_frac[1];
    poly_frac_build(pf2);
    poly_frac_set_pf(pf2, pfin2);
    pf1 = pfin1;
  } else if (pfout == pfin1 && pfout == pfin2) {
    ovrwrt = 3;
    pf1 = new struct poly_frac[1];
    pf2 = new struct poly_frac[1];
    poly_frac_build(pf1);
    poly_frac_set_pf(pf1, pfin1);
    poly_frac_build(pf2);
    poly_frac_set_pf(pf2, pfin2);
  } else {
    ovrwrt = 0;
    pf1 = pfin1;
    pf2 = pfin2;
  }

  // multiply numerators
  // cout << "multiply numerators" << endl;
  int ldeg1 = (*pf1).num_deg - (*pf1).num_vdeg;
  int ldeg2 = (*pf2).num_deg - (*pf2).num_vdeg;
  if ((*pfout).coeffs) {
    mpc_rk1_clear((*pfout).coeffs, (*pfout).num_vdeg+1);
    delete[] (*pfout).coeffs;
  }
  (*pfout).num_deg = (*pfin1).num_deg + (*pfin2).num_deg;
  (*pfout).num_vdeg =  (*pfout).num_deg - ldeg1 - ldeg2;
  (*pfout).coeffs = new mpc_t[(*pfout).num_vdeg+1];
  init_rk1_mpc((*pfout).coeffs, (*pfout).num_vdeg+1);
  // cout << "num_deg = " << (*pfout).num_deg << endl;
  // cout << "num_vdeg = " << (*pfout).num_vdeg << endl;
  mul_vpoly(
    (*pfout).coeffs, (*pf1).coeffs, (*pf2).coeffs,
    (*pf1).num_deg, (*pf2).num_deg,
    (*pf1).num_vdeg, (*pf2).num_vdeg,
    ldeg1, ldeg2
  );
  // cout << "numerator:" << endl;
  // print_poly((*pfout).coeffs, (*pfout).num_vdeg);

  // multiply denominators
  // cout << "multiply denominators" << endl;
  // cout << "den degs: " << (*pf1).den_deg << ", " << (*pf2).den_deg << endl;
  if (!(*pfout).mults) {
    (*pfout).mults = new int[(*pf1).nroots];
  }
  (*pfout).nroots = (*pf1).nroots;
  (*pfout).den_deg = 0;
  for (int k=0; k<(*pf1).nroots; k++) {
    (*pfout).mults[k] = (*pf1).mults[k] + (*pf2).mults[k];
    (*pfout).den_deg += (*pfout).mults[k];
  }
  // if ((*pf1).den_deg != 0 && (*pf2).den_deg != 0) {
  //   if (!(*pfout).mults) {
  //     (*pfout).mults = new int[(*pf1).nroots];
  //   }
  //   (*pfout).nroots = (*pf1).nroots;
  //   (*pfout).den_deg = 0;
  //   for (int k=0; k<(*pf1).nroots; k++) {
  //     (*pfout).mults[k] = (*pf1).mults[k] + (*pf2).mults[k];
  //     (*pfout).den_deg += (*pfout).mults[k];
  //   }
  // } else if ((*pf1).den_deg != 0 && (*pf2).den_deg == 0) {
  //   if (!(*pfout).mults) {
  //     (*pfout).mults = new int[(*pf1).nroots];
  //   }
  //   (*pfout).nroots = (*pf1).nroots;
  //   (*pfout).den_deg = 0;
  //   for (int k=0; k<(*pf1).nroots; k++) {
  //     (*pfout).mults[k] = (*pf1).mults[k];
  //     (*pfout).den_deg += (*pfout).mults[k];
  //   }
  // } else if ((*pf1).den_deg == 0 && (*pf2).den_deg != 0) {
  //   if (!(*pfout).mults) {
  //     (*pfout).mults = new int[(*pf2).nroots];
  //   }
  //   (*pfout).nroots = (*pf2).nroots;
  //   (*pfout).den_deg = 0;
  //   for (int k=0; k<(*pf2).nroots; k++) {
  //     (*pfout).mults[k] = (*pf2).mults[k];
  //     (*pfout).den_deg += (*pfout).mults[k];
  //   }
  // } else if ((*pf1).den_deg == 0 && (*pf2).den_deg == 0) {
  //   cout << "both den degs are zero" << endl;
  //   (*pfout).nroots = 0;
  //   (*pfout).den_deg = 0;
  //   cout << "about to return:" << endl;
  //   poly_frac_print(pfout);
  //   cout << "nroot = " << (*pfout).nroots << endl;
  //   cout << pfout << endl;
  //   cout << "delete internal pf" << endl;
  //   if (pfout == pfin1) {
  //     cout << "delete 1st one" << endl;
  //     delete[] (*pf1).coeffs;
  //     if ((*pf1).mults) {
  //       delete[] (*pf1).mults;
  //     }
  //   }
  //   if (pfout == pfin2) {
  //     cout << "delete 2nd one" << endl;
  //     delete[] (*pf2).coeffs;
  //     if ((*pf2).mults) {
  //       delete[] (*pf2).mults;
  //     }
  //   }

  //   cout << "return" << endl;
  //   return;
  // }

  // normal
  poly_frac_normal(pfout);

  if (ovrwrt == 1) {
    poly_frac_free(pf1);
    delete[] pf1;
  } else if (ovrwrt == 2) {
    poly_frac_free(pf2);
    delete[] pf2;
  } else if (ovrwrt == 3) {
    poly_frac_free(pf1);
    delete[] pf1;
    poly_frac_free(pf2);
    delete[] pf2;
  }
}


void poly_frac_mul_sym_pow(
  struct poly_frac *pfout,
  struct poly_frac *pfin,
  int pow
) {
  // if input is zero, set output to zero
  if ((*pfin).num_vdeg == -1) {
    if (pfout == pfin) {
      return;
    }
    poly_frac_set_zero(pfout);
    return;
  }

  if (pfout != pfin) {
    poly_frac_set_pf(pfout, pfin);
  }

  // if power is zero, do nothing
  if (pow == 0) {
    return;
  }

  int pow_behav = (*pfin).num_deg - (*pfin).num_vdeg - (*pfin).mults[0] + pow;
  if (pow_behav >= 0) {
    (*pfout).num_deg = (*pfin).num_vdeg + pow_behav;
    (*pfout).den_deg = (*pfin).den_deg - (*pfin).mults[0];
    (*pfout).mults[0] = 0;
  } else if (pow_behav < 0) {
    (*pfout).num_deg = (*pfin).num_vdeg;
    (*pfout).den_deg = (*pfin).den_deg - pow_behav - (*pfin).mults[0];
    (*pfout).mults[0] = -pow_behav;
  }
}


void poly_frac_rk2_mul_sym_pow(
  struct poly_frac **pfmat, int pow,
  int dim1, int dim2
) {
  for (int i1=0; i1<dim1; i1++) {
    for (int i2=0; i2<dim2; i2++) {
      poly_frac_mul_sym_pow(&pfmat[i1][i2], &pfmat[i1][i2], pow);
    }
  }
}


void poly_frac_mul_root(
  struct poly_frac *pf,
  mpc_t *root, int mult
) {
  // if the input is zero, do nothing
  if ((*pf).num_vdeg == -1) {
    return;
  }

  mpc_t *original_coeffs;
  original_coeffs = (*pf).coeffs;
  (*pf).coeffs = new mpc_t[(*pf).num_vdeg+mult+1];
  init_rk1_mpc((*pf).coeffs, (*pf).num_vdeg+mult+1);
  // cout << "allocated new coeffs" << endl;
  // cout << "old memory:" << endl;
  // print_poly(original_coeffs, (*pf).num_vdeg);
  // cout << ";" << endl;
  poly_mul_root_new_out((*pf).coeffs, original_coeffs, (*pf).num_vdeg, *root);
  // cout << "multiplied 1st root" << endl;
  for (int k=1; k<mult; k++) {
    // cout << "k = " << k << endl;
    poly_mul_root((*pf).coeffs, (*pf).num_vdeg+k, *root);
    // cout << "multiplied other root" << endl;
  }
  mpc_rk1_clear(original_coeffs, (*pf).num_vdeg+1);
  delete[] original_coeffs;
  (*pf).num_deg += mult;
  (*pf).num_vdeg += mult;  
}


void poly_frac_rk2_mul_root(
  struct poly_frac **pfmat,
  mpc_t *root, int mult,
  int dim1, int dim2
) {
  for (int i1=0; i1<dim1; i1++) {
    for (int i2=0; i2<dim2; i2++) {
      poly_frac_mul_root(&pfmat[i1][i2], root, mult);
    }
  }
}


void poly_frac_den_to_poly(
  // OUTPUT
  mpc_t *coeffs,
  // INPUT
  struct poly_frac *pfin, mpc_t *roots
) {
  mpc_set_ui(coeffs[0], 1, MPFR_RNDN);
  // for cycle to mul mults
  int count=0;
  for (int k=0; k<(*pfin).nroots; k++) {
    for (int m=0; m<(*pfin).mults[k]; m++) {
      poly_mul_root(coeffs, count, roots[k]);
      count++;
    }
  }
}


void poly_frac_den_to_num(
  // OUPUT
  struct poly_frac *pfout,
  // INPUT
  struct poly_frac *pfin,
  mpc_t *roots,
  int wp2
) {
  //////
  // BUILD NUMERATOR
  //////
  if ((*pfout).coeffs) {
    mpc_rk1_clear((*pfout).coeffs, (*pfout).num_vdeg+1);
    delete[] (*pfout).coeffs;
    // cout << "deleting" << endl;
  }
  (*pfout).num_deg = (*pfin).den_deg;
  if ((*pfin).den_deg > 0) {
    (*pfout).num_vdeg = (*pfin).den_deg - (*pfin).mults[0]; 
  } else {
    (*pfout).num_vdeg = 0;
  }
  (*pfout).coeffs = new mpc_t[(*pfout).num_vdeg+1];
  // cout << "den to num wp2 = " << wp2 << endl;
  for (int k=0; k<=(*pfout).num_vdeg; k++) {
    mpc_init3((*pfout).coeffs[k], wp2, wp2);
  }
  // init_rk1_mpc((*pfout).coeffs, (*pfout).num_vdeg+1);
  if ((*pfout).num_vdeg > 0) {
    poly_frac_den_to_poly((*pfout).coeffs, pfin, roots);
  } else {
    mpc_set_ui((*pfout).coeffs[0], 1, MPFR_RNDN);
  }
}


void poly_frac_pow_ui(
  struct poly_frac *pfout,
  struct poly_frac *pfin,
  int pow
) {
  if (pow == 0) {
    poly_frac_set_ui(pfout, 1, (*pfin).nroots);
    return;
  }

  if (pow == 1) {
    if (pfout != pfin) {
      poly_frac_set_pf(pfout, pfin);
      return;
    }
  }

  // deal with overwriting
  struct poly_frac *pf;
  int ovrwrt;
  if (pfout == pfin) {
    ovrwrt = 1;
    pf = new struct poly_frac[1];
    poly_frac_build(pf);
    poly_frac_set_pf(pf, pfin);
  } else {
    ovrwrt = 0;
    pf = pfin;
  }

  // poly_frac_set_ui(pfout, 1, (*pfin).nroots);
  poly_frac_set_pf(pfout, pf);
  for (int i=1; i<pow; i++) {
    poly_frac_mul_pf(pfout, pfout, pf);
  }

  if (ovrwrt) {
    poly_frac_free(pf);
    delete[] pf;
  }


}


void poly_frac_pow_si(
  // OUTPUT
  struct poly_frac *pfout,
  mpc_t **roots, mpfr_t **tols,
  // INPUT
  struct poly_frac *pfin,
  int pow,
  mpc_t *in_roots
) {
  if (pow >= 0) {
    poly_frac_pow_ui(pfout, pfin, pow);
    return;
  }

  // poly_frac_inv(pfout, roots, tols, pfin, in_roots);

  // for (int i=1; i<-pow; i++) {
  //   poly_frac_mul_pf(pfout, pfout, pfout);
  // }
}


void poly_frac_add_pf(
  struct poly_frac *pfout,
  struct poly_frac *pfin1,
  struct poly_frac *pfin2,
  mpc_t *roots,
  mpc_t *flop_mpc
) {
  // check whether 1st is zero
  if ((*pfin1).num_vdeg == -1) {
    // cout << "1st is zero" << endl;
    // if the mpc exists, check if it is null
    if (flop_mpc != NULL) {
      if (mpc_lessthan_tol(*flop_mpc)) {
        // return the null polynomial fraction
        poly_frac_set_zero(pfout);
        // if ((*pfout).coeffs) {
        //   delete[] (*pfout).coeffs;
        //   (*pfout).coeffs = NULL;
        // }
        // (*pfout).num_vdeg = -1;
        return;
      } else {
        poly_frac_mul_mpc(pfout, pfin2, flop_mpc);
        return;
      }
    }
    // return the 2nd
    if (pfout != pfin2) {
      // cout << "set output to 2nd addend" << endl;
      // poly_frac_print(pfin2, nroots);
      poly_frac_set_pf(pfout, pfin2);
      // cout << "done" << endl;
    }
    return;
  }
  // check whether 2nd is zero
  if ((*pfin2).num_vdeg == -1) {
    if (pfout != pfin1) {
      poly_frac_set_pf(pfout, pfin1);
    }
    return;
  }
  if (flop_mpc != NULL) {
    if (mpc_lessthan_tol(*flop_mpc)) {
      if (pfout != pfin1) {
        poly_frac_set_pf(pfout, pfin1);
      }
      return;
    }
  }

  int k, m, count;
  // compute LCM between denominators, and degrees
  (*pfout).nroots = (*pfin1).nroots;
  if (!(*pfout).mults) {
    (*pfout).mults = new int[(*pfin1).nroots];
  }

  // deal with overwriting
  int *mults, ovwrt_mults;
  if (pfout == pfin1 || pfout == pfin2) {
    mults = new int[(*pfin1).nroots];
    ovwrt_mults = 1;
  } else {
    mults = (*pfout).mults;
    ovwrt_mults = 0;
  }

  int num_deg1 = (*pfin1).num_deg, num_deg2 = (*pfin2).num_deg;
  int num_ldeg1 = (*pfin1).num_deg - (*pfin1).num_vdeg;
  int num_ldeg2 = (*pfin2).num_deg - (*pfin2).num_vdeg;
  (*pfout).den_deg = 0;
  int tmp_int;
  //// null root
  if (tmp_int = (*pfin1).mults[0] - (*pfin2).mults[0], tmp_int > 0) {
    mults[0] = (*pfin1).mults[0];
    num_deg2 += tmp_int;
    num_ldeg2 += tmp_int;
  } else if (tmp_int = (*pfin2).mults[0] - (*pfin1).mults[0], tmp_int > 0) {
    mults[0] = (*pfin2).mults[0];
    num_deg1 += tmp_int;
    num_ldeg1 += tmp_int;
  } else if ((*pfin1).mults[0] == (*pfin2).mults[0]) {
    mults[0] = (*pfin1).mults[0];
  }
  (*pfout).den_deg += mults[0];
  //// other roots
  for (k=1; k<(*pfin1).nroots; k++) {
    if ((*pfin1).mults[k] > (*pfin2).mults[k]) {
      mults[k] = (*pfin1).mults[k];
      num_deg2 += (*pfin1).mults[k] - (*pfin2).mults[k];
    } else if ((*pfin1).mults[k] < (*pfin2).mults[k]) {
      mults[k] = (*pfin2).mults[k];
      num_deg1 += (*pfin2).mults[k] - (*pfin1).mults[k];
    } else if ((*pfin1).mults[k] == (*pfin2).mults[k]) {
      mults[k] = (*pfin1).mults[k];
    }
    (*pfout).den_deg += mults[k];
  }
  int num_vdeg1 = num_deg1 - num_ldeg1;
  int num_vdeg2 = num_deg2 - num_ldeg2;
  (*pfout).num_deg = max(num_deg1, num_deg2);
  int num_ldeg = min(num_ldeg1, num_ldeg2);
  int num_vdeg = (*pfout).num_deg - num_ldeg;
  // cout << "num_ldeg1 = " << num_ldeg1 << endl;
  // cout << "num_deg1 = " << num_deg1 << endl;
  // cout << "num_ldeg2 = " << num_ldeg2 << endl;
  // cout << "num_deg2 = " << num_deg2 << endl;

  // compute coefficients of the numerator
  //// contribution from 1st addend
  //// (the coefficients are already stored into the output)
  mpc_t *num1 = new mpc_t[num_vdeg1+1];
  init_rk1_mpc(num1, num_vdeg1+1);
  copy_poly(num1, (*pfin1).coeffs, (*pfin1).num_vdeg);
  count = 0;
  for (k=1; k<(*pfout).nroots; k++) {
    // cout << "k = " << k << endl;
    for (m=0; m<mults[k]-(*pfin1).mults[k]; m++) {
      // cout << "m = " << m << endl;
      poly_mul_root(num1, (*pfin1).num_vdeg+count, roots[k]);
      count++;
    }
  }
  // cout << "preparato 1st" << endl;
  // print_poly(num1, num_vdeg1);
  //// contribution from 2nd addend
  mpc_t *num2 = new mpc_t[num_vdeg2+1];
  init_rk1_mpc(num2, num_vdeg2+1);
  if (flop_mpc == NULL) {
    copy_poly(num2, (*pfin2).coeffs, (*pfin2).num_vdeg);
  } else {
    for (k=0; k<=(*pfin2).num_vdeg; k++) {
      mpc_mul(num2[k], (*pfin2).coeffs[k], *flop_mpc, MPFR_RNDN);
    }
  }
  count = 0;
  for (k=1; k<(*pfout).nroots; k++) {
    for (m=0; m<mults[k]-(*pfin2).mults[k]; m++) {
      poly_mul_root(num2, (*pfin2).num_vdeg+count, roots[k]);
      count++;
    }
  }
  // cout << "preparato 2nd" << endl;
  // print_poly(num2, num_vdeg2);

  if (ovwrt_mults) {
    for (k=0; k<(*pfout).nroots; k++) {
      (*pfout).mults[k] = mults[k];
    }
    delete[] mults;
  }
  //// sum of the two contributions
  if ((*pfout).coeffs) {
    mpc_rk1_clear((*pfout).coeffs, (*pfout).num_vdeg+1);
    delete[] (*pfout).coeffs;
  }
  (*pfout).coeffs = new mpc_t[num_vdeg+1];
  init_rk1_mpc((*pfout).coeffs, num_vdeg+1);
  int found1, found2, kp;
  for (k=num_ldeg; k<=(*pfout).num_deg; k++) {
    kp = k - num_ldeg;
    // cout << "k, kp = " << k << ", " << kp << endl;
    if (k >= num_ldeg1 && k <= num_deg1) {
      found1 = 1;
    } else {
      found1 = 0;
    }
    if (k >= num_ldeg2 && k <= num_deg2) {
      found2 = 1;
    } else {
      found2 = 0;
    }

    if (found1 == 0 && found2 == 0) {
      mpc_set_ui((*pfout).coeffs[kp], 0, MPFR_RNDN);
    } else if (found1 && found2 == 0) {
      mpc_set((*pfout).coeffs[kp], num1[k-num_ldeg1], MPFR_RNDN);
    } else if (found1 == 0 && found2) {
      mpc_set((*pfout).coeffs[kp], num2[k-num_ldeg2], MPFR_RNDN);
    } else if (found1 && found2) {
      mpc_add((*pfout).coeffs[kp], num1[k-num_ldeg1], num2[k-num_ldeg2], MPFR_RNDN);
    }

    // cout << "found1, found2 = " << found1 << ", " << found2 << endl;
    // print_mpc(&(*pfout).coeffs[kp]);
    // cout << endl;
  }

  // lower padding
  int lpad = 0, upad = num_vdeg;
  for (k=0; k<=num_vdeg; k++) {
    // if(mpc_lessthan_tol((*pfout).coeffs[k])) {
    if(mpc_zero_p((*pfout).coeffs[k])) {
      lpad++;
      continue;
    } else {
      break;
    }
  }
  // cout << "lpad = " << lpad << endl;
  if (lpad == num_vdeg + 1) {
    mpc_rk1_clear((*pfout).coeffs, num_vdeg+1);
    delete[] (*pfout).coeffs;
    (*pfout).coeffs = NULL;
    (*pfout).num_vdeg = -1;
    (*pfout).den_deg = 0;
    return;
  }

  // upper padding
  for (k=num_vdeg; k>=0; k--) {
    // if(mpc_lessthan_tol((*pfout).coeffs[k])) {
    if(mpc_zero_p((*pfout).coeffs[k])) {
      upad--;
      continue;
    } else {
      break;
    }
  }
  // cout << "upad = " << upad << endl;

  // prune null coefficients
  if (lpad>0 || upad<num_vdeg1) {
    int old_vdeg = num_vdeg;
    mpc_t *unpruned_coeffs;
    unpruned_coeffs = (*pfout).coeffs;
    (*pfout).num_deg -= num_vdeg - upad;
    num_vdeg = upad-lpad;
    (*pfout).coeffs = new mpc_t[num_vdeg+1];
    for (k=0; k<=num_vdeg; k++) {
      mpc_init3((*pfout).coeffs[k], wp2, wp2);
      mpc_set((*pfout).coeffs[k], unpruned_coeffs[lpad+k], MPFR_RNDN);
    }
    mpc_rk1_clear(unpruned_coeffs, old_vdeg+1);
    delete[] unpruned_coeffs;
  }

  // (in case it is really needed)
  // simplify roots from numerator and LCM

  // (*pfout).coeffs = new mpc_t[num_vdeg+1];
  // for (k=0; k<=num_vdeg; k++) {
  //   mpc_init3((*pfout).coeffs[k], wp2, wp2);
  //   mpc_set((*pfout).coeffs[k], num[k], MPFR_RNDN);
  // }
  (*pfout).num_vdeg = num_vdeg;
  // (*pfout).num_deg = num_deg;

  // normal
  poly_frac_normal(pfout);

  mpc_rk1_clear(num1, num_vdeg1+1); delete[] num1;
  mpc_rk1_clear(num2, num_vdeg2+1); delete[] num2;
}


void poly_frac_rk2_add_pf_rk2(
  struct poly_frac **pf_tens_out,
  struct poly_frac **pf_tens_in1,
  struct poly_frac **pf_tens_in2,
  mpc_t *roots, mpc_t *flop_mpc,
  int dim1, int dim2
) {
  for (int i1=0; i1<dim1; i1++) {
    for (int i2=0; i2<dim2; i2++) {
      poly_frac_add_pf(
        &pf_tens_out[i1][i2], &pf_tens_in1[i1][i2], &pf_tens_in2[i1][i2],
        roots, flop_mpc
      );
    }
  }
}


//////
// EVALUATION AND DIFFERENTIATION FUNCTIONS
//////
int poly_frac_eval_sym_zero(
  // OUTPUT
  mpc_t *value,
  // INPUT
  struct poly_frac *pf, mpc_t *roots
) {
  // if input is zero, return zero
  if ((*pf).num_vdeg == -1) {
    mpc_set_ui(*value, 0, MPFR_RNDN);
    return 0;
  }
  int pow_behav = (*pf).num_deg - (*pf).num_vdeg - (*pf).mults[0];
  // cout << "num_deg = " << (*pf).num_deg << endl;
  // cout << "num_vdeg = " << (*pf).num_vdeg << endl;
  // cout << "mults[0] = " << (*pf).mults[0] << endl;
  // cout << "pow_behav" << pow_behav << endl;
  if (pow_behav < 0) {
    fprintf(stderr, "Trying to evaluate in zero singular fraction of polynomials! Exiting...\n");
    printf("poly_frac:\n");
    poly_frac_print(pf);
    // return -1;
    exit(EXIT_FAILURE);
    // return -1;
  } else if (pow_behav > 0) {
    mpc_set_ui(*value, 0, MPFR_RNDN);
    return 0;
  } else if (pow_behav == 0) {
    poly_eval_sym_zero_by_roots_wo_0(
      value,
      roots, (*pf).mults, (*pf).nroots
    );
    mpc_div(*value, (*pf).coeffs[0], *value, MPFR_RNDN);
    return 1;
  }
  return -1;
}


void poly_frac_rk1_eval_sym_zero(
  // OUTPUT
  mpc_t *mpc_tens,
  // INPUT
  struct poly_frac *pf_tens, mpc_t *roots, int dim
) {
  for (int i=0; i<dim; i++) {
    // cout << "eval in zero i = " << i << endl;
    // poly_frac_print(&pf_tens[i]);
    if (poly_frac_eval_sym_zero(&mpc_tens[i], &pf_tens[i], roots) == -1) {
      cout << "singular pf at i == " << i << endl;
    };
  }
}


void poly_frac_eval_value(
  // OUTPUT
  mpc_t *out,
  // INPUT
  struct poly_frac *pf, mpc_t *roots, mpc_t *val
) {
  // if input is zero, return zero
  // poly_frac_print(pf);
  // cout << "vdeg = " << pf->num_vdeg << endl;

  if ((*pf).num_vdeg == -1) {
    // cout << "set zero" << endl;
    mpc_set_ui(*out, 0, MPFR_RNDN);
    return;
  }

  // if point is zero, call proper routine
  if (mpc_lessthan_tol(*val)) {
    // cout << "going to dedicated routine" << endl;
    poly_frac_eval_sym_zero(out, pf, roots);
    return;
  }
  
  poly_frac_normal(pf);
  mpc_t tmp;
  mpc_init3(tmp, wp2, wp2);

  // EVALUATE NUMERATOR
  poly_eval(out, (*pf).coeffs, (*pf).num_vdeg, val);
  // cout << "eval num = "; print_mpc(out); cout << endl;
  // multiply lower degree
  if ((*pf).num_vdeg < (*pf).num_deg) {
    int ldeg = (*pf).num_deg - (*pf).num_vdeg;
    mpc_pow_ui(tmp, *val, ldeg, MPFR_RNDN);
    mpc_mul(*out, *out, tmp, MPFR_RNDN);
  }
  // cout << "num" << endl; print_mpc(out); cout << endl;

  // EVALUATE DENOMINATOR
  poly_eval_by_roots(&tmp, roots, (*pf).mults, (*pf).nroots, val);
  // cout << "den" << endl; print_mpc(&tmp); cout << endl;

  // DIVIDE NUMERATOR BY DENOMINATOR
  mpc_div(*out, *out, tmp, MPFR_RNDN);
  // cout << "out" << endl;
  // print_mpc(out); cout << endl;
  
}


void poly_frac_rk2_eval_value(
  // OUTPUT
  mpc_t **mpc_tens,
  // INPUT
  struct poly_frac **pf_tens, mpc_t *roots, mpc_t *val,
  int dim1, int dim2
) {
  // cout << "val = "; print_mpc(val); cout << endl;
  for (int i1=0; i1<dim1; i1++) {
    for (int i2=0; i2<dim2; i2++) {
      // cout << "i1, i2 = " << i1 << ", " << i2 << endl;
      poly_frac_eval_value(
        &mpc_tens[i1][i2],
        &pf_tens[i1][i2], roots, val
      );
      // if (i1 == 31) {
      //   cout << "matrix element eval:" << endl;
      //   print_mpc(&mpc_tens[i1][i2]); cout << endl;
      // }
    }
  }
}


void poly_frac_extract_LO(
  // OUTPUT
  mpfr_t *out_re, mpfr_t *out_im,
  // INPUT
  struct poly_frac *pf, mpc_t *roots, int pow
) {
  // if input is zero, return zero
  if ((*pf).num_vdeg == -1) {
    mpfr_set_ui(*out_re, 0, MPFR_RNDN);
    mpfr_set_ui(*out_im, 0, MPFR_RNDN);
    return;
  }
  int pow_behav = (*pf).num_deg - (*pf).num_vdeg - (*pf).mults[0] + pow;
  // cout << "pow_behav = " << pow_behav << endl;
  if (pow_behav < 0) {
    fprintf(stderr, "Trying to evaluate in zero singular fraction of polynomials! Exiting...\n");
    printf("poly_frac:\n");
    poly_frac_print(pf);
    exit(EXIT_FAILURE);
  } else if (pow_behav > 0) {
    mpfr_set_ui(*out_re, 0, MPFR_RNDN);
    mpfr_set_ui(*out_im, 0, MPFR_RNDN);
  } else if (pow_behav == 0) {
    mpc_t tmp;
    mpc_init3(tmp, wp2, wp2);
    // cout << "evaluate denominator in zero" << endl;
    poly_eval_sym_zero_by_roots_wo_0(
      &tmp,
      roots, (*pf).mults, (*pf).nroots
    );
    mpc_div(tmp, (*pf).coeffs[0], tmp, MPFR_RNDN);
    mpc_real(*out_re, tmp, MPFR_RNDN);
    mpc_imag(*out_im, tmp, MPFR_RNDN);
  }
}


void poly_frac_mat_extract_LO(
  // OUTPUT
  mpfr_t *LOre, mpfr_t *LOim,
  // INPUT
  struct poly_frac  **pfmat, mpc_t *roots, int pow,
  int dim1, int dim2
) {
  for (int i1=0; i1<dim1; i1++) {
    for (int i2=0; i2<dim2; i2++) {
      // cout << "i1, i2 = " << i1 << ", " << i2 << endl;
      // poly_frac_print(&pfmat[i1][i2]);
      poly_frac_extract_LO(
        &LOre[i1 + dim1*i2], &LOim[i1 + dim1*i2],
        &pfmat[i1][i2], roots, pow
      );
    }
  }
}


void poly_frac_mat_extract_LOc(
  // OUTPUT
  mpc_t **LO,
  // INPUT
  struct poly_frac  **pfmat, mpc_t *roots, int pow,
  int dim1, int dim2
) {
  for (int i1=0; i1<dim1; i1++) {
    for (int i2=0; i2<dim2; i2++) {
      poly_frac_extract_LO(
        &mpc_realref(LO[i1][i2]), &mpc_imagref(LO[i1][i2]),
        &pfmat[i1][i2], roots, pow
      );
    }
  }
}


void poly_frac_extract_NLO(
  // OUTPUT
  mpc_t *NLO,
  // INPUT
  mpc_t *LO, struct poly_frac *pf, mpc_t *roots, int pow
) {
  if ((*pf).num_vdeg == -1) {
    mpc_set_ui(*NLO, 0, MPFR_RNDN);
    return;
  }
  
  int pow_behav = (*pf).num_deg - (*pf).num_vdeg - (*pf).mults[0] + pow;
  if (pow_behav < 0) {
    fprintf(stderr, "Trying to evaluate in zero singular fraction of polynomial! Exiting...\n");
    printf("poly_frac:\n");
    poly_frac_print(pf);
    exit(EXIT_FAILURE);
  } else if (pow_behav > 1) {
    mpc_set_ui(*NLO, 0, MPFR_RNDN);
  } else if (pow_behav == 1) {
    poly_eval_sym_zero_by_roots_wo_0(
      NLO,
      roots, (*pf).mults, (*pf).nroots
    );
    mpc_div(*NLO, (*pf).coeffs[0], *NLO, MPFR_RNDN);
  } else if (pow_behav == 0) {
    poly_coeff1_by_roots_wo_0(
      NLO,
      roots, (*pf).mults, (*pf).nroots
    );
    mpc_mul(*NLO, *NLO, *LO, MPFR_RNDN);
    if ((*pf).num_deg > 0) {
      mpc_sub(*NLO, (*pf).coeffs[1], *NLO, MPFR_RNDN);
    } else {
      mpc_neg(*NLO, *NLO, MPFR_RNDN);
    }
    mpc_t tmp;
    mpc_init3(tmp, wp2, wp2);
    poly_eval_sym_zero_by_roots_wo_0(
      &tmp,
      roots, (*pf).mults, (*pf).nroots
    );
    mpc_div(*NLO, *NLO, tmp, MPFR_RNDN);

    // FREE
    mpc_clear(tmp);
  }
}


void poly_frac_mat_extract_NLOc(
  // OUTPUT
  mpc_t **NLO,
  // INPUT
  mpfr_t *LOre, mpfr_t *LOim,
  struct poly_frac  **pfmat, mpc_t *roots, int pow,
  int dim1, int dim2
) {
  int index;
  mpc_t tmp;
  mpc_init3(tmp, wp2, wp2);
  for (int i1=0; i1<dim1; i1++) {
    for (int i2=0; i2<dim2; i2++) {
      index =  i1 + dim1*i2;
      mpc_set_fr_fr(tmp, LOre[index], LOim[index], MPFR_RNDN);
      poly_frac_extract_NLO(
        &NLO[i1][i2],
        &tmp, &pfmat[i1][i2], roots, pow
      );
    }
  }
}


void poly_frac_extract_oder(
  // OUTPUT
  mpc_t *out,
  // INPUT
  struct poly_frac *pf, mpc_t *roots, int nroots, int order
) {
  int pow_behav = poly_frac_pow_behav(pf);
  struct poly_frac pfin;
  poly_frac_build(&pfin);
  poly_frac_set_pf(&pfin, pf);
  struct poly_frac pflead;
  poly_frac_build(&pflead);
  poly_frac_mul_sym_pow(&pfin, &pfin, -pow_behav);
  // poly_frac_print(&pfin);

  for (int p=pow_behav; p<=order; p++) {
    poly_frac_eval_sym_zero(out, &pfin, roots);
    if (p<order) {
      mpc_neg(*out, *out, MPFR_RNDN);    
      poly_frac_set_mpc(&pflead, out, nroots);
      poly_frac_add_pf(&pfin, &pfin, &pflead, roots, NULL);
      poly_frac_mul_sym_pow(&pfin, &pfin, -1);
    }
  }

}


void poly_frac_derivative(
  struct poly_frac *pfout,
  struct poly_frac *pfin,
  mpc_t *roots
) {
  // if input is zero, return input
  if ((*pfin).num_vdeg == -1) {
    if (pfout != pfin) {
      poly_frac_set_zero(pfout);
    }
    return;
  }

  // BUILD 1ST CONTRIBUTION: Der[num]/den
  // cout <<  "BUILD 1ST CONTRIBUTION: Der[num]/den" << endl;
  struct poly_frac pf1;
  poly_frac_build(&pf1);
  if ((*pfin).num_deg == 0) {
    poly_frac_set_zero(&pf1);
  } else {
    int deg_offset;
    int ldeg1 = (*pfin).num_deg - (*pfin).num_vdeg;;
    pf1.num_deg = (*pfin).num_deg - 1;
    if (ldeg1 == 0) {
      pf1.num_vdeg = (*pfin).num_vdeg - 1;
      deg_offset = 1;
    } else {
      pf1.num_vdeg = (*pfin).num_vdeg;
      deg_offset = 0;
    }
    pf1.den_deg = (*pfin).den_deg;
    pf1.nroots = (*pfin).nroots;
    pf1.mults = new int[(*pfin).nroots];
    for (int k=0; k<(*pfin).nroots; k++) {
      pf1.mults[k] = (*pfin).mults[k];
    }
    pf1.coeffs = new mpc_t[pf1.num_vdeg+1];
    for (int k=0; k<=pf1.num_vdeg; k++) {
      mpc_init3(pf1.coeffs[k], wp2, wp2);
      mpc_mul_ui(pf1.coeffs[k], (*pfin).coeffs[k+deg_offset], ldeg1+k+deg_offset, MPFR_RNDN);
    }    
  }
  // cout << "1st addend:" << endl;
  // poly_frac_print(&pf1);


  // BUILD 2ND CONTRIBUTION: -num*Der[den]/den^2
  // cout <<  "BUILD 2ND CONTRIBUTION: -num*Der[den]/den^2" << endl;
  struct poly_frac pf2;
  poly_frac_build(&pf2);
  if ((*pfin).den_deg == 0) {
    poly_frac_set_pf(pfout, &pf1);
    return;
  }
  pf2.num_deg = (*pfin).num_deg + (*pfin).den_deg - 1;
  pf2.num_vdeg = pf2.num_deg - (*pfin).mults[0] + 1;
  pf2.den_deg = 2 * (*pfin).den_deg;
  pf2.nroots = (*pfin).nroots;
  pf2.mults = new int[(*pfin).nroots];
  // cout << "num_deg = " << pf2.num_deg << endl;
  // cout << "num_vdeg = " << pf2.num_vdeg << endl;
  // cout << "den_deg = " << pf2.den_deg << endl;
  for (int k=0; k<(*pfin).nroots; k++) {
    pf2.mults[k] = 2*(*pfin).mults[k];
  }
  pf2.coeffs = new mpc_t[pf2.num_vdeg+1];
  // init_rk1_mpc(pf2.coeffs, pf2.num_vdeg+1);
  for (int k=0; k<=pf2.num_vdeg; k++) {
    mpc_init3(pf2.coeffs[k], wp2, wp2);
    mpc_set_ui(pf2.coeffs[k], 0, MPFR_RNDN);
  }

  mpc_t *der_den = new mpc_t[(*pfin).den_deg-(*pfin).mults[0]+1];
  init_rk1_mpc(der_den, (*pfin).den_deg-(*pfin).mults[0]+1);
  mpc_t *num_der_den = new mpc_t[pf2.num_vdeg+1];
  init_rk1_mpc(num_der_den, pf2.num_vdeg+1);
  int count = 0;
  for (int ldeg, k=0; k<(*pfin).nroots; k++) {
    if ((*pfin).mults[k] == 0) {
      continue;
    }
    // cout << "k = " << k << endl;
    // cout << "root:" << endl;
    // print_mpc(&roots[k]);
    // cout << endl;
    // compute contribution to the derivative due to the k-th root
    mpc_set_ui(der_den[0], (*pfin).mults[k], MPFR_RNDN);
    count = 0;
    for (int mt, kp=1; kp<(*pfin).nroots; kp++) {
      if (kp == k) {
        mt = (*pfin).mults[kp] - 1;
      } else {
        mt = (*pfin).mults[kp];
      }
      for (int mp=0; mp<mt; mp++) {
        poly_mul_root(der_den, count, roots[kp]);
        count++;
      }
    }
    // cout << "D[den]:" << endl;
    // print_poly(der_den, (*pfin).den_deg-(*pfin).mults[0]);

    // multiply the contribution times the numerator
    if (k == 0) {
      ldeg = (*pfin).mults[0] - 1;
    } else {
      ldeg = (*pfin).mults[0];
    }
    mul_vpoly(
      num_der_den, (*pfin).coeffs, der_den,
      (*pfin).num_deg, (*pfin).den_deg - 1,
      (*pfin).num_vdeg, (*pfin).den_deg - 1 - ldeg,
      (*pfin).num_deg - (*pfin).num_vdeg, ldeg
    );
    // cout << "num*D[den]:" << endl;
    // print_poly(num_der_den, pf2.num_vdeg);

    // cumulate
    add_poly(pf2.coeffs, pf2.coeffs, num_der_den, pf2.num_vdeg, pf2.num_vdeg);
  }
  poly_frac_normal(&pf2);
  poly_frac_neg(&pf2);
  // cout << "2nd addend" << endl;
  // poly_frac_print(&pf2);

  // SUM CONTRIBUTIONS
  // cout << "SUM CONTRIBUTIONS" << endl;
  poly_frac_add_pf(pfout, &pf1, &pf2, roots, NULL);
}


//////
// MATRIX MULTIPLICATION FUNCTIONS
//////
void mpc_rk2_mul_poly_frac_rk2(
  struct poly_frac **pf_tens_out,
  mpc_t **mpc_tens,
  struct poly_frac **pf_tens_in,
  int dim,
  mpc_t *roots
) {
  struct poly_frac **pf_tens;
  int to_be_del;
  if (pf_tens_out == pf_tens_in) {
    to_be_del = 1;
    // copy input
    malloc_rk2_tens(pf_tens, dim, dim);
    for (int i1=0; i1<dim; i1++) {
      for (int i2=0; i2<dim; i2++) {
        poly_frac_build(&pf_tens[i1][i2]);
        poly_frac_set_pf(&pf_tens[i1][i2], &pf_tens_in[i1][i2]);
      }
    }
  } else {
    to_be_del = 0;
    pf_tens = pf_tens_in;
  }

  // row by column product
  for (int i1=0; i1<dim; i1++) {
    for (int i2=0; i2<dim; i2++) {
      // pf_tens_out[i1][i2].mults = new int[nroots];
      poly_frac_mul_mpc(
        &pf_tens_out[i1][i2], &pf_tens[0][i2], &mpc_tens[i1][0]
      );
      for (int i3=1; i3<dim; i3++) {
        poly_frac_add_pf(
          &pf_tens_out[i1][i2], &pf_tens_out[i1][i2], &pf_tens[i3][i2],
          roots,
          &mpc_tens[i1][i3]
        );
      }
    }
  }

  if (to_be_del) {
    poly_frac_rk2_free(pf_tens, dim, dim);
  }
}


void poly_frac_rk2_mul_mpc_rk2(
  struct poly_frac **pf_tens_out,
  struct poly_frac **pf_tens_in,
  mpc_t **mpc_tens,
  int dim,
  mpc_t *roots
) {
  struct poly_frac **pf_tens;
  int to_be_del;
  if (pf_tens_out == pf_tens_in) {
    to_be_del = 1;
    // copy input
    malloc_rk2_tens(pf_tens, dim, dim);
    for (int i1=0; i1<dim; i1++) {
      for (int i2=0; i2<dim; i2++) {
        poly_frac_build(&pf_tens[i1][i2]);
        poly_frac_set_pf(&pf_tens[i1][i2], &pf_tens_in[i1][i2]);
      }
    }
  } else {
    to_be_del = 0;
    pf_tens = pf_tens_in;
  }
  
  // row by column product
  for (int i1=0; i1<dim; i1++) {
    for (int i2=0; i2<dim; i2++) {
      // cout << "i, j = " << i1 << ", " << i2 << endl;
      // cout << "1st term:" << endl;
      // poly_frac_print(&pf_tens_in[i1][0]);
      // cout << "2nd term" << endl;
      // print_mpc(&mpc_tens[0][i2]);
      // cout << endl;
      poly_frac_mul_mpc(
        &pf_tens_out[i1][i2], &pf_tens[i1][0], &mpc_tens[0][i2]
      );
      // cout << "result" << endl;
      // poly_frac_print(&pf_tens_out[i1][i2]);
      // cout << "mul" << endl;
      for (int i3=1; i3<dim; i3++) {
        // cout << "k = " << i3 << endl;
        // cout << "1st addend" << endl;
        // poly_frac_print(&pf_tens_out[i1][i2]);
        // cout << "2nd addend" << endl;
        // poly_frac_print(&pf_tens[i1][i3]);
        // cout << "times" << endl;
        // print_mpc(&mpc_tens[i3][i2]);
        // cout << endl;
        poly_frac_add_pf(
          &pf_tens_out[i1][i2], &pf_tens_out[i1][i2], &pf_tens[i1][i3],
          roots,
          &mpc_tens[i3][i2]
        );
        // cout << "result" << endl;
        // poly_frac_print(&pf_tens_out[i1][i2]);
      }
      // cout << "flop" << endl;
    }
  }

  if (to_be_del) {
    // cout << "trying to delete" << endl;
    poly_frac_rk2_free(pf_tens, dim, dim);
    // cout << "done" << endl;
  }
}


void poly_frac_rk2_mul_pf_rk1(
  struct poly_frac *pf_tens_out,
  struct poly_frac **pf_tens_in1,
  struct poly_frac *pf_tens_in2,
  int dim1, int dim2,
  mpc_t *roots
) {
  struct poly_frac **pf_tens1, *pf_tens2;
  int to_be_del;
  if (pf_tens_out == pf_tens_in2) {
    to_be_del = 2;
    // copy input
    pf_tens2 = new struct poly_frac[dim2];
    for (int i2=0; i2<dim2; i2++) {
      poly_frac_build(&pf_tens2[i2]);
      poly_frac_set_pf(&pf_tens2[i2], &pf_tens_in2[i2]);
    }
    pf_tens1 = pf_tens_in1;
  } else {
    to_be_del = 0;
    pf_tens1 = pf_tens_in1;
    pf_tens2 = pf_tens_in2;
  }

  // row by column product
  struct poly_frac pf_tmp;
  poly_frac_build(&pf_tmp);
  for (int i1=0; i1<dim1; i1++) {
    poly_frac_mul_pf(
      &pf_tens_out[i1], &pf_tens1[i1][0], &pf_tens2[0]
    );
    for (int i2=1; i2<dim2; i2++) {
      // cout << "i1, i2 = " << i1 << ", " << i2 << endl;
      poly_frac_mul_pf(
        &pf_tmp, &pf_tens1[i1][i2], &pf_tens2[i2]
      );
      poly_frac_add_pf(
        &pf_tens_out[i1], &pf_tens_out[i1], &pf_tmp,
        roots,
        NULL
      );
    }
  }

  if (to_be_del == 2) {
    // cout << "trying to delete" << endl;
    del_rk1_poly_frac(pf_tens2, dim2);
    // cout << "done" << endl;
  }
}


void poly_frac_rk2_mul_pf_rk2(
  struct poly_frac **pf_tens_out,
  struct poly_frac **pf_tens_in1,
  struct poly_frac **pf_tens_in2,
  int dim1, int dim, int dim2,
  mpc_t *roots
) {
  struct poly_frac **pf_tens1, **pf_tens2;
  int to_be_del;
  if (pf_tens_out == pf_tens_in1) {
    to_be_del = 1;
    // copy input
    malloc_rk2_tens(pf_tens1, dim1, dim);
    for (int i1=0; i1<dim1; i1++) {
      for (int i2=0; i2<dim; i2++) {
        poly_frac_build(&pf_tens1[i1][i2]);
        poly_frac_set_pf(&pf_tens1[i1][i2], &pf_tens_in1[i1][i2]);
      }
    }
    pf_tens2 = pf_tens_in2;
  } else if (pf_tens_out == pf_tens_in2) {
    to_be_del = 2;
    // copy input
    malloc_rk2_tens(pf_tens2, dim, dim2);
    for (int i1=0; i1<dim; i1++) {
      for (int i2=0; i2<dim2; i2++) {
        poly_frac_build(&pf_tens2[i1][i2]);
        poly_frac_set_pf(&pf_tens2[i1][i2], &pf_tens_in2[i1][i2]);
      }
    }
    pf_tens1 = pf_tens_in1;
  } else {
    to_be_del = 0;
    pf_tens1 = pf_tens_in1;
    pf_tens2 = pf_tens_in2;
  }
  
  // row by column product
  struct poly_frac pf_tmp;
  poly_frac_build(&pf_tmp);
  for (int i1=0; i1<dim1; i1++) {
    for (int i2=0; i2<dim2; i2++) {
      // cout << "i, j = " << i1 << ", " << i2 << endl;
      // cout << "1st factor" << endl;
      // poly_frac_print(&pf_tens1[i1][0]);
      // cout << "2nd factor" << endl;
      // poly_frac_print(&pf_tens2[0][i2]);
      poly_frac_mul_pf(
        &pf_tens_out[i1][i2], &pf_tens1[i1][0], &pf_tens2[0][i2]
      );
      // cout << "result" << endl;
      // poly_frac_print(&pf_tens_out[i1][i2]);
      for (int i3=1; i3<dim; i3++) {
        // cout << "k = " << i3 << endl;
        // cout << "1st factor" << endl;
        // poly_frac_print(&pf_tens1[i1][i3]);
        // cout << "2nd factor" << endl;
        // poly_frac_print(&pf_tens2[i3][i2]);
        poly_frac_mul_pf(
          &pf_tmp, &pf_tens1[i1][i3], &pf_tens2[i3][i2]
        );
        // cout << "result" << endl;
        // poly_frac_print(&pf_tmp);
        poly_frac_add_pf(
          &pf_tens_out[i1][i2], &pf_tens_out[i1][i2], &pf_tmp,
          roots,
          NULL
        );
        // cout << "result" << endl;
        // poly_frac_print(&pf_tens_out[i1][i2]);
      }
    }
  }

  if (to_be_del == 1) {
    // cout << "trying to delete" << endl;
    poly_frac_rk2_free(pf_tens1, dim1, dim);
    // cout << "done" << endl;
  } else if (to_be_del == 2) {
    // cout << "trying to delete" << endl;
    poly_frac_rk2_free(pf_tens2, dim, dim2);
    // cout << "done" << endl;
  }
}


//////
// OTHER MATRIX FUNCTIONS
//////
void pf_find_profile(
  // OUTPUT
  int **&profile, int **&subblock, int &nblocks,
  // INPUT
  struct poly_frac **mat, int g_dim
) {
  int i, j;

  // temporary profile that stores the position
  // of the last non-zero element of each row
  int *tmpprof;
  tmpprof = new int [g_dim];
  for (i = 0; i < g_dim; i++) {
    tmpprof[i] = i;
    for (j = g_dim-1; j >= i; j--) {
      if(mat[i][j].num_vdeg != -1) {
        tmpprof[i] = j;
        break;
      }
    }
  }

  // extract number of blocks
  nblocks=0;
  for (i = 0; i < g_dim; i++)
    if(tmpprof[i] == i) {
      nblocks++;
    }

  // nblocks X 2 matrix to store block ranges
  profile = new int * [nblocks];
  for (i = 0; i < nblocks; i++) {
    profile[i] = new int [2];
  }

  // compute block ranges from the temporary profile
  profile[0][0] = 0;
  profile[0][1] = tmpprof[0];
  int k = 1;
  for (i = 1; i < g_dim; i++) {
    if(tmpprof[i] == i) {
      profile[k][0] = profile[k-1][1] + 1;
      profile[k][1] = i;
      k++;
    }
  }
  delete[] tmpprof;
  
  // grid to store the positions of non-vanishing sub-blocks
  subblock = new int * [nblocks];
  for (i = 0; i < nblocks; i++) {
    subblock[i] = new int [nblocks];
  }
  // initialize to zero every element
  for (i = 0; i < nblocks; i++) {
    for (j = 0; j < nblocks; j++) { 
      subblock[i][j] = 0;
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
          if(mat[kk][ll].num_vdeg != -1) {
            found=1;
            break;
          }
        }
        if(found == 1) {
          subblock[i][j] = 1;
          break;
        }
      }
    }
  }
}


int poly_frac_Poinc_rank(
  struct poly_frac **pfmat, int dim1, int dim2
) {
  int p_rank = 0;
  for (int i1=0; i1<dim1; i1++) {
    for (int i2=0; i2<dim2; i2++) {
      if (pfmat[i1][i2].num_vdeg == -1) {
        continue;
      }
      p_rank = max(p_rank, pfmat[i1][i2].mults[0]);
    }
  }

  return p_rank - 1;
}


int poly_frac_Poinc_rank_sh(
  struct poly_frac **pfmat, int dim1, int dim2, int sh
) {
  int p_rank = 0;
  for (int i1=0; i1<dim1; i1++) {
    for (int i2=0; i2<dim2; i2++) {
      if (pfmat[i1][i2].num_vdeg == -1) {
        continue;
      }
      p_rank = max(p_rank, pfmat[i1][i2].mults[sh]);
    }
  }

  return p_rank - 1;
}


int poly_frac_k2_num_deg_max(
  struct poly_frac **pfmat, int dim1, int dim2
) {
  int num_deg_max = -1;
  for (int i1=0; i1<dim1; i1++) {
    for (int i2=0; i2<dim2; i2++) {
      if (pfmat[i1][i2].num_vdeg == -1) {
        continue;
      }
      if (pfmat[i1][i2].num_deg > num_deg_max) {
        cout << "new num_deg_max at i1, i2 = " << i1 << ", " << i2 << endl;
        poly_frac_print(&pfmat[i1][i2]);
        num_deg_max = pfmat[i1][i2].num_deg;
      }
    }
  }

  return num_deg_max;
}


void poly_frac_k2_num_exp2_avg(
  // OUTPUT
  double *num_exp2_avg,
  // INPUT
  struct poly_frac **pfmat, int dim1, int dim2,
  int Poinc_rank, int num_deg_max
) {
  int ldeg, k;

  int *ncontr = new int[num_deg_max+Poinc_rank+2];

  // initialize to zero
  for (k=0; k<num_deg_max+Poinc_rank+2; k++) {
    num_exp2_avg[k] = 0;
    ncontr[k] = 0;
  }

  // cumulate
  for (int i1=0; i1<dim1; i1++) {
    for (int i2=0; i2<dim2; i2++) {
      if (pfmat[i1][i2].num_vdeg == -1) {
        continue;
      }
      ldeg = pfmat[i1][i2].num_deg - pfmat[i1][i2].num_vdeg; 
      for (k=0; k<=pfmat[i1][i2].num_vdeg; k++) {
        cout << "k = " << k << endl;
        if (mpc_zero_p(pfmat[i1][i2].coeffs[k])) {
          continue;
        }
        num_exp2_avg[Poinc_rank+1-pfmat[i1][i2].mults[0]+ldeg+k] += mpc_get_exp2(pfmat[i1][i2].coeffs[k]);
        ncontr[Poinc_rank+1-pfmat[i1][i2].mults[0]+ldeg+k]++;
      }
    }
  }

  // normalize
  cout << "normalize" << endl;
  for (k=0; k<num_deg_max+Poinc_rank+2; k++) {
    if (ncontr[k] > 0) {
      num_exp2_avg[k] /= ncontr[k];
    } else {
      num_exp2_avg[k] = -1;
    }
  }
}


//////
// UTILS
//////
void poly_frac_roots_append_zero(
  // IN-OUT
  struct poly_frac *pf,
  // INPUT
  int nroots
) {
  if ((*pf).nroots < nroots) {
    // cout << "update nroots of inner pf n." << n << endl;
    // cout << "old nroots = " << inner_pf[i].nroots << endl;
    // cout << "new nroots = " << *nroots << endl;
    if ((*pf).mults) {
      // cout << "append" << endl;
      // cout << "old vec:" << endl;
      // for (int k=0; k<inner_pf[i].nroots; k++) {
      // 	cout << inner_pf[i].mults[k] << endl;
      // }
      int_rk1_append_zero(&(*pf).nroots, &(*pf).mults, nroots - (*pf).nroots);
      // cout << "new vec:" << endl;
      // for (int k=0; k<inner_pf[i].nroots; k++) {
      // 	cout << inner_pf[i].mults[k] << endl;
      // }
    } else {
      // cout << "alloc" << endl;
      (*pf).mults = new int[nroots];
      for (int k=0; k<nroots; k++) {
        (*pf).mults[k] = 0;
      }
      (*pf).nroots = nroots;
    }
  }
}


void poly_frac_rk2_slice(
  struct poly_frac ***sub_tens,
  struct poly_frac ***tens,
  int dim1_start, int dim1_end, int dim2_start
) {
  // cout << "before starting" << endl;
  // cout << *tens << endl;
  // cout << *sub_tens << endl;
  *sub_tens = *tens + dim1_start;
  // cout << "temporaneo" << endl;
  // cout << *sub_tens << endl;
  // cout << "dopo primo shift" << endl;
  // *sub_tens += dim1_start;
  // cout << (*tens) + dim1_start << endl;
  // cout << *sub_tens << endl;
  // cout << "singoli shift" << endl;
  for (int i1=0; i1<=dim1_end-dim1_start; i1++) {
    // cout << "i = " << i1 << endl;
    // cout << "in1: " << (*sub_tens)[i1] << endl;
    // cout << "in2: " << (*tens)[i1+dim1_start] << endl;
    // cout << "goal: " << (*tens)[dim1_start+i1] + dim2_start << endl;
    // cout << "try : " << (*sub_tens)[i1] + dim2_start << endl;
    (*sub_tens)[i1] += dim2_start;
    // cout << "original: " << (*tens)[dim1_start+i1] + dim2_start << endl;
    // cout << "selected: " << (*sub_tens)[i1] << endl;
  }
}


void poly_frac_rk2_unslice(
  struct poly_frac ***tens,
  int dim1_start, int dim1_end, int dim2_start
) {
  for (int i1=0; i1<=dim1_end-dim1_start; i1++) {
    (*tens)[i1] -= dim2_start;
  }
  *tens -= dim1_start;
}


void poly_frac_rk2_subtens(
  struct poly_frac ***sub_tens,
  struct poly_frac ***tens,
  int dim1_start, int dim1_end, int dim2_start
) {
  for (int i1=0; i1<=dim1_end-dim1_start; i1++) {
    (*sub_tens)[i1] = (*tens)[dim1_start+i1] + dim2_start;
  }
}


int poly_frac_cmp_pf(
  struct poly_frac *pf1,
  struct poly_frac *pf2
) {
  // check if they are both zeros
  if ((*pf1).num_vdeg == -1 && (*pf1).num_vdeg == -1) {
    return 1;
  }
  // compare properties
  if ((*pf1).num_deg != (*pf2).num_deg) {
    cout << "different num deg" << endl;
    return 0;
  }
  if ((*pf1).num_vdeg != (*pf2).num_vdeg) {
    cout << "different num virtual deg" << endl;
    return 0;
  }
  if ((*pf1).den_deg != (*pf2).den_deg) {
    cout << "different den deg" << endl;
    return 0;
  }
  if ((*pf1).nroots != (*pf2).nroots) {
    cout << "different nroots: " << (*pf1).nroots << ", " << (*pf2).nroots << endl;
    return 0;
  }

  // compare multiplicities
  for (int k=0; k<(*pf1).nroots; k++) {
    if ((*pf1).mults[k] != (*pf2).mults[k]) {
      cout << "different multiplicities k = " << k << ": " << (*pf1).mults[k] << ", " << (*pf2).mults[k] << endl;
      return 0;
    }
  }

  // compare coefficients
  for (int k=0; k<(*pf1).num_vdeg; k++) {
    if (!mpc_equal_within_tol((*pf1).coeffs[k], (*pf2).coeffs[k])) {
      cout << "different coefficients for k = " << k << endl;
      return 0;
    }
  }

  return 1;
}


//////
// CONVERSIONS
//////

void root_prof_to_poly_frac_den(
  struct poly_frac *pf,
  int *root_prof, int deg,
  int nroots //, int zero_label
) {
  if ((*pf).mults) {
    delete[] (*pf).mults;
  }
  (*pf).mults = new int[nroots];
  (*pf).nroots = nroots;

  for (int k=0; k<nroots; k++) {
    (*pf).mults[k] = 0;
  }

  // #DEPRECATED
  // zero_label was used when the list of roots was not
  // initialized with zero at the beginning, thus root_prof
  // had to be incremented by one when adding zero at the
  // beginning of the list
  //
  // if (zero_label == -1) {
  //   for (int k=0; k<deg; k++) {
  //     (*pf).mults[root_prof[k]+1]++;
  //   }
  // }
  // else {
  //   for (int k=0; k<deg; k++) {
  //     if (root_prof[k] < zero_label) {
  //       (*pf).mults[root_prof[k]+1]++;
  //     } else if (root_prof[k] == zero_label) {
  //       (*pf).mults[0]++;
  //     } else if (root_prof[k] > zero_label) {
  //       (*pf).mults[root_prof[k]]++;
  //     }
  //   }
  // }
  for (int k=0; k<deg; k++) {    
    (*pf).mults[root_prof[k]]++;
  }

}


void root_prof_rk2_to_poly_frac_den_rk2(
  struct poly_frac **pf_tens,
  int ***root_prof,
  int nroots,// int zero_label,
  int dim1, int dim2
) {
  for (int i1=0; i1<dim1; i1++) {
    for (int i2=0; i2<dim2; i2++) {
      // if (i1 == 0 && i2 == 0) {
      //   cout << "den deg = " << pf_tens[i1][i2].den_deg << endl;
      //   for (int k=0; k<pf_tens[i1][i2].den_deg; k++) {
      //     cout << "root_prof = " << root_prof[i1][i2][k] << endl;
      //   }
      // }
      root_prof_to_poly_frac_den(
        &pf_tens[i1][i2], root_prof[i1][i2],
        pf_tens[i1][i2].den_deg, nroots //, zero_label
      );
      // if (i1 == 9 && i2 == 0) {
      //   poly_frac_print(&pf_tens[i1][i2]);
      // }
    }
  }
}



//////
// READ/WRITE
//////
void poly_frac_rk2_to_file(
  char *file_name,
  struct poly_frac **pfmat, int dim1, int dim2
) {
  FILE *fptr = fopen(file_name, "w");
  for (int i=0; i<dim1; i++) {
    for (int j=0; j<dim2; j++) {
      // virtual degree
      fprintf(fptr, "%d\n", pfmat[i][j].num_vdeg);
      if (pfmat[i][j].num_vdeg == -1) {continue;}
      // numerator degree
      fprintf(fptr, "%d\n", pfmat[i][j].num_deg);
      // denominator degree
      fprintf(fptr, "%d\n", pfmat[i][j].den_deg);
      // number of roots
      fprintf(fptr, "%d\n", pfmat[i][j].nroots);
      // numerator coefficients
      for (int k=0; k<=pfmat[i][j].num_vdeg; k++) {
        mpc_out_str(fptr, 10, 0, pfmat[i][j].coeffs[k], MPFR_RNDN);
        fprintf(fptr, "\n");
      }
      // denominator multiplicities
      for (int k=0; k<pfmat[i][j].nroots; k++) {
        fprintf(fptr, "%d\n", pfmat[i][j].mults[k]);
      }
    }
    fprintf(fptr, "\n");
  }
  fclose(fptr);
}


void poly_frac_rk2_from_file(
  char *file_name,
  struct poly_frac **pfmat, int dim1, int dim2
) {
  FILE *fptr = fopen(file_name, "r");
  
  // check if file exists
  if (fptr == NULL) {
    printf("Could not open file %s\n", file_name);
    exit(1);
  }

  // load matrix
  for (int i=0; i<dim1; i++) {
    for (int j=0; j<dim2; j++) {
      if (pfmat[i][j].coeffs) {
        mpc_rk1_clear(pfmat[i][j].coeffs, pfmat[i][j].num_vdeg+1);
        delete[] pfmat[i][j].coeffs;
      }
      // virtual degree
      fscanf(fptr, "%d", &pfmat[i][j].num_vdeg);
      if (pfmat[i][j].num_vdeg == -1) {
        poly_frac_set_zero(&pfmat[i][j]);
        continue;
      }
      // numerator degree
      fscanf(fptr, "%d\n", &pfmat[i][j].num_deg);
      // denominator degree
      fscanf(fptr, "%d\n", &pfmat[i][j].den_deg);
      // number of roots
      fscanf(fptr, "%d\n", &pfmat[i][j].nroots);
      // numerator coefficients
      pfmat[i][j].coeffs = new mpc_t[pfmat[i][j].num_vdeg+1];
      for (int k=0; k<=pfmat[i][j].num_vdeg; k++) {
        mpc_init3(pfmat[i][j].coeffs[k], wp2, wp2);
        mpc_inp_str(pfmat[i][j].coeffs[k], fptr, 0, 10, MPFR_RNDN);
      }
      // denominator multiplicities
      if (pfmat[i][j].mults) {
        delete[] pfmat[i][j].mults;
      }
      pfmat[i][j].mults = new int[pfmat[i][j].nroots];
      for (int k=0; k<pfmat[i][j].nroots; k++) {
        fscanf(fptr, "%d\n", &pfmat[i][j].mults[k]);
      }
    }
  }
  fclose(fptr);
}


//////
// LOCAL PRECISION
//////
void poly_frac_set_pf_wp2(
	int wp2,
  struct poly_frac *pfout,
  struct poly_frac *pfin
) {
  // check if input is equal to outpu
  if (pfout == pfin) {
    return;
  }

  // check if input is zero
  if ((*pfin).num_vdeg == -1) {
    poly_frac_set_zero(pfout);
    return;
  }

  // cout << "set degrees" << endl;
  (*pfout).num_deg = (*pfin).num_deg;
  (*pfout).den_deg = (*pfin).den_deg;
  (*pfout).nroots = (*pfin).nroots;
  if ((*pfin).nroots > 0) {
    // cout << "alloc mults" << endl;
    if (!(*pfout).mults) {
      (*pfout).mults = new int[(*pfin).nroots];
    }
    // cout << "set mults" << endl;
    for (int k=0; k<(*pfin).nroots; k++) {
      (*pfout).mults[k] = (*pfin).mults[k];
    }
  }
  // cout << "alloc coeffs" << endl;
  if ((*pfout).coeffs) {
    mpc_rk1_clear((*pfout).coeffs, (*pfout).num_vdeg+1);
    delete[] (*pfout).coeffs;
  }
  (*pfout).num_vdeg = (*pfin).num_vdeg;
  (*pfout).coeffs = new mpc_t[(*pfin).num_vdeg+1];
  // cout << "set coeffs" << endl;
  for (int k=0; k<=(*pfin).num_vdeg; k++) {
    // cout << "k = " << k << endl;
    mpc_init3((*pfout).coeffs[k], wp2, wp2);
    // print_mpc(&(*pfin).coeffs[k]);
    // cout << endl;
    mpc_set((*pfout).coeffs[k], (*pfin).coeffs[k], MPFR_RNDN);
  }
  // cout << "done" << endl;
}


void poly_frac_add_pf_wp2(
	int wp2,
  struct poly_frac *pfout,
  struct poly_frac *pfin1,
  struct poly_frac *pfin2,
  mpc_t *roots,
  mpc_t *flop_mpc
) {
  // check whether 1st is zero
  if ((*pfin1).num_vdeg == -1) {
    // cout << "1st is zero" << endl;
    // if the mpc exists, check if it is null
    if (flop_mpc != NULL) {
      if (mpc_lessthan_tol(*flop_mpc)) {
        // return the null polynomial fraction
        poly_frac_set_zero(pfout);
        // if ((*pfout).coeffs) {
        //   delete[] (*pfout).coeffs;
        //   (*pfout).coeffs = NULL;
        // }
        // (*pfout).num_vdeg = -1;
        return;
      } else {
        poly_frac_mul_mpc(pfout, pfin2, flop_mpc);
        return;
      }
    }
    // return the 2nd
    if (pfout != pfin2) {
      // cout << "set output to 2nd addend" << endl;
      // poly_frac_print(pfin2, nroots);
      poly_frac_set_pf_wp2(wp2, pfout, pfin2);
      // cout << "done" << endl;
    }
    return;
  }
  // check whether 2nd is zero
  if ((*pfin2).num_vdeg == -1) {
    if (pfout != pfin1) {
      poly_frac_set_pf_wp2(wp2, pfout, pfin1);
    }
    return;
  }
  if (flop_mpc != NULL) {
    if (mpc_lessthan_tol(*flop_mpc)) {
      if (pfout != pfin1) {
        poly_frac_set_pf_wp2(wp2, pfout, pfin1);
      }
      return;
    }
  }

  int k, m, count;
  // compute LCM between denominators, and degrees
  (*pfout).nroots = (*pfin1).nroots;
  if (!(*pfout).mults) {
      (*pfout).mults = new int[(*pfin1).nroots];
  }

  // deal with overwriting
  int *mults, ovwrt_mults;
  if (pfout == pfin1 || pfout == pfin2) {
    mults = new int[(*pfin1).nroots];
    ovwrt_mults = 1;
  } else {
    mults = (*pfout).mults;
    ovwrt_mults = 0;
  }

  int num_deg1 = (*pfin1).num_deg, num_deg2 = (*pfin2).num_deg;
  int num_ldeg1 = (*pfin1).num_deg - (*pfin1).num_vdeg;
  int num_ldeg2 = (*pfin2).num_deg - (*pfin2).num_vdeg;
  (*pfout).den_deg = 0;
  int tmp_int;
  //// null root
  if (tmp_int = (*pfin1).mults[0] - (*pfin2).mults[0], tmp_int > 0) {
    mults[0] = (*pfin1).mults[0];
    num_deg2 += tmp_int;
    num_ldeg2 += tmp_int;
  } else if (tmp_int = (*pfin2).mults[0] - (*pfin1).mults[0], tmp_int > 0) {
    mults[0] = (*pfin2).mults[0];
    num_deg1 += tmp_int;
    num_ldeg1 += tmp_int;
  } else if ((*pfin1).mults[0] == (*pfin2).mults[0]) {
    mults[0] = (*pfin1).mults[0];
  }
  (*pfout).den_deg += mults[0];
  //// other roots
  for (k=1; k<(*pfin1).nroots; k++) {
    if ((*pfin1).mults[k] > (*pfin2).mults[k]) {
      mults[k] = (*pfin1).mults[k];
      num_deg2 += (*pfin1).mults[k] - (*pfin2).mults[k];
    } else if ((*pfin1).mults[k] < (*pfin2).mults[k]) {
      mults[k] = (*pfin2).mults[k];
      num_deg1 += (*pfin2).mults[k] - (*pfin1).mults[k];
    } else if ((*pfin1).mults[k] == (*pfin2).mults[k]) {
      mults[k] = (*pfin1).mults[k];
    }
    (*pfout).den_deg += mults[k];
  }
  int num_vdeg1 = num_deg1 - num_ldeg1;
  int num_vdeg2 = num_deg2 - num_ldeg2;
  (*pfout).num_deg = max(num_deg1, num_deg2);
  int num_ldeg = min(num_ldeg1, num_ldeg2);
  int num_vdeg = (*pfout).num_deg - num_ldeg;
  // cout << "num_ldeg1 = " << num_ldeg1 << endl;
  // cout << "num_deg1 = " << num_deg1 << endl;
  // cout << "num_ldeg2 = " << num_ldeg2 << endl;
  // cout << "num_deg2 = " << num_deg2 << endl;

  // compute coefficients of the numerator
  //// contribution from 1st addend
  //// (the coefficients are already stored into the output)
  mpc_t *num1 = new mpc_t[num_vdeg1+1];
  init_rk1_mpc_wp2(wp2, num1, num_vdeg1+1);
  copy_poly(num1, (*pfin1).coeffs, (*pfin1).num_vdeg);
  count = 0;
  for (k=1; k<(*pfout).nroots; k++) {
    // cout << "k = " << k << endl;
    for (m=0; m<mults[k]-(*pfin1).mults[k]; m++) {
      // cout << "m = " << m << endl;
      poly_mul_root(num1, (*pfin1).num_vdeg+count, roots[k]);
      count++;
    }
  }
  // cout << "preparato 1st" << endl;
  // print_poly(num1, num_vdeg1);
  //// contribution from 2nd addend
  mpc_t *num2 = new mpc_t[num_vdeg2+1];
  init_rk1_mpc_wp2(wp2, num2, num_vdeg2+1);
  if (flop_mpc == NULL) {
    copy_poly(num2, (*pfin2).coeffs, (*pfin2).num_vdeg);
  } else {
    for (k=0; k<=(*pfin2).num_vdeg; k++) {
      mpc_mul(num2[k], (*pfin2).coeffs[k], *flop_mpc, MPFR_RNDN);
    }
  }
  count = 0;
  for (k=1; k<(*pfout).nroots; k++) {
    for (m=0; m<mults[k]-(*pfin2).mults[k]; m++) {
      poly_mul_root(num2, (*pfin2).num_vdeg+count, roots[k]);
      count++;
    }
  }
  // cout << "preparato 2nd" << endl;
  // print_poly(num2, num_vdeg2);

  if (ovwrt_mults) {
    for (k=0; k<(*pfout).nroots; k++) {
      (*pfout).mults[k] = mults[k];
    }
    delete[] mults;
  }
  //// sum of the two contributions
  if ((*pfout).coeffs) {
    mpc_rk1_clear((*pfout).coeffs, (*pfout).num_vdeg+1);
    delete[] (*pfout).coeffs;
  }
  (*pfout).coeffs = new mpc_t[num_vdeg+1];
  init_rk1_mpc_wp2(wp2, (*pfout).coeffs, num_vdeg+1);
  int found1, found2, kp;
  for (k=num_ldeg; k<=(*pfout).num_deg; k++) {
    kp = k - num_ldeg;
    // cout << "k, kp = " << k << ", " << kp << endl;
    if (k >= num_ldeg1 && k <= num_deg1) {
      found1 = 1;
    } else {
      found1 = 0;
    }
    if (k >= num_ldeg2 && k <= num_deg2) {
      found2 = 1;
    } else {
      found2 = 0;
    }

    if (found1 == 0 && found2 == 0) {
      mpc_set_ui((*pfout).coeffs[kp], 0, MPFR_RNDN);
    } else if (found1 && found2 == 0) {
      mpc_set((*pfout).coeffs[kp], num1[k-num_ldeg1], MPFR_RNDN);
    } else if (found1 == 0 && found2) {
      mpc_set((*pfout).coeffs[kp], num2[k-num_ldeg2], MPFR_RNDN);
    } else if (found1 && found2) {
      mpc_add((*pfout).coeffs[kp], num1[k-num_ldeg1], num2[k-num_ldeg2], MPFR_RNDN);
    }

    // cout << "found1, found2 = " << found1 << ", " << found2 << endl;
    // print_mpc(&(*pfout).coeffs[kp]);
    // cout << endl;
  }

  // lower padding
  int lpad = 0, upad = num_vdeg;
  for (k=0; k<=num_vdeg; k++) {
    if(mpc_lessthan_tol((*pfout).coeffs[k])) {
      lpad++;
      continue;
    } else {
      break;
    }
  }
  // cout << "lpad = " << lpad << endl;
  if (lpad == num_vdeg + 1) {
    mpc_rk1_clear((*pfout).coeffs, num_vdeg+1);
    delete[] (*pfout).coeffs;
    (*pfout).coeffs = NULL;
    (*pfout).num_vdeg = -1;
    (*pfout).den_deg = 0;
    return;
  }

  // upper padding
  for (k=num_vdeg; k>=0; k--) {
    if(mpc_lessthan_tol((*pfout).coeffs[k])) {
      upad--;
      continue;
    } else {
      break;
    }
  }
  // cout << "upad = " << upad << endl;

  // prune null coefficients
  if (lpad>0 || upad<num_vdeg1) {
    int old_vdeg = num_vdeg;
    mpc_t *unpruned_coeffs;
    unpruned_coeffs = (*pfout).coeffs;
    (*pfout).num_deg -= num_vdeg - upad;
    num_vdeg = upad-lpad;
    (*pfout).coeffs = new mpc_t[num_vdeg+1];
    for (k=0; k<=num_vdeg; k++) {
      mpc_init3((*pfout).coeffs[k], wp2, wp2);
      mpc_set((*pfout).coeffs[k], unpruned_coeffs[lpad+k], MPFR_RNDN);
    }
    mpc_rk1_clear(unpruned_coeffs, old_vdeg+1);
    delete[] unpruned_coeffs;
  }

  // (in case it is really needed)
  // simplify roots from numerator and LCM

  // (*pfout).coeffs = new mpc_t[num_vdeg+1];
  // for (k=0; k<=num_vdeg; k++) {
  //   mpc_init3((*pfout).coeffs[k], wp2, wp2);
  //   mpc_set((*pfout).coeffs[k], num[k], MPFR_RNDN);
  // }
  (*pfout).num_vdeg = num_vdeg;
  // (*pfout).num_deg = num_deg;

  // normal
  poly_frac_normal(pfout);

  mpc_rk1_clear(num1, num_vdeg1+1); delete[] num1;
  mpc_rk1_clear(num2, num_vdeg2+1); delete[] num2;
}


void poly_frac_mul_pf_wp2(
	int wp2,
  struct poly_frac *pfout,
  struct poly_frac *pfin1,
  struct poly_frac *pfin2
) {
  // if input is zero, set output to zero
  if ((*pfin1).num_vdeg == -1) {
    if (pfout != pfin1) {
      poly_frac_set_zero(pfout);
    }
    return;
  }
  if ((*pfin2).num_vdeg == -1) {
    if (pfout != pfin2) {
      poly_frac_set_zero(pfout);
    }
    return;
  }

  // deal with overwriting
  // cout << "deal with overwriting" << endl;
  struct poly_frac *pf1, *pf2;
  int ovrwrt;
  if (pfout == pfin1 && pfout != pfin2) {
    ovrwrt = 1;
    pf1 = new struct poly_frac[1];
    poly_frac_build(pf1);
    poly_frac_set_pf_wp2(wp2, pf1, pfin1);
    pf2 = pfin2;
  } else if (pfout != pfin1 && pfout == pfin2) {
    ovrwrt = 2;
    pf2 = new struct poly_frac[1];
    poly_frac_build(pf2);
    poly_frac_set_pf_wp2(wp2, pf2, pfin2);
    pf1 = pfin1;
  } else if (pfout == pfin1 && pfout == pfin2) {
    ovrwrt = 3;
    pf1 = new struct poly_frac[1];
    pf2 = new struct poly_frac[1];
    poly_frac_build(pf1);
    poly_frac_set_pf_wp2(wp2, pf1, pfin1);
    poly_frac_build(pf2);
    poly_frac_set_pf_wp2(wp2, pf2, pfin2);
  } else {
    ovrwrt = 0;
    pf1 = pfin1;
    pf2 = pfin2;
  }

  // multiply numerators
  // cout << "multiply numerators" << endl;
  int ldeg1 = (*pf1).num_deg - (*pf1).num_vdeg;
  int ldeg2 = (*pf2).num_deg - (*pf2).num_vdeg;
  if ((*pfout).coeffs) {
    mpc_rk1_clear((*pfout).coeffs, (*pfout).num_vdeg+1);
    delete[] (*pfout).coeffs;
  }
  (*pfout).num_deg = (*pfin1).num_deg + (*pfin2).num_deg;
  (*pfout).num_vdeg =  (*pfout).num_deg - ldeg1 - ldeg2;
  (*pfout).coeffs = new mpc_t[(*pfout).num_vdeg+1];
  init_rk1_mpc_wp2(wp2, (*pfout).coeffs, (*pfout).num_vdeg+1);
  // cout << "num_deg = " << (*pfout).num_deg << endl;
  // cout << "num_vdeg = " << (*pfout).num_vdeg << endl;
  mul_vpoly(
    (*pfout).coeffs, (*pf1).coeffs, (*pf2).coeffs,
    (*pf1).num_deg, (*pf2).num_deg,
    (*pf1).num_vdeg, (*pf2).num_vdeg,
    ldeg1, ldeg2
  );
  // cout << "numerator:" << endl;
  // print_poly((*pfout).coeffs, (*pfout).num_vdeg);

  // multiply denominators
  // cout << "multiply denominators" << endl;
  // cout << "den degs: " << (*pf1).den_deg << ", " << (*pf2).den_deg << endl;
  if (!(*pfout).mults) {
    (*pfout).mults = new int[(*pf1).nroots];
  }
  (*pfout).nroots = (*pf1).nroots;
  (*pfout).den_deg = 0;
  for (int k=0; k<(*pf1).nroots; k++) {
    (*pfout).mults[k] = (*pf1).mults[k] + (*pf2).mults[k];
    (*pfout).den_deg += (*pfout).mults[k];
  }
  // if ((*pf1).den_deg != 0 && (*pf2).den_deg != 0) {
  //   if (!(*pfout).mults) {
  //     (*pfout).mults = new int[(*pf1).nroots];
  //   }
  //   (*pfout).nroots = (*pf1).nroots;
  //   (*pfout).den_deg = 0;
  //   for (int k=0; k<(*pf1).nroots; k++) {
  //     (*pfout).mults[k] = (*pf1).mults[k] + (*pf2).mults[k];
  //     (*pfout).den_deg += (*pfout).mults[k];
  //   }
  // } else if ((*pf1).den_deg != 0 && (*pf2).den_deg == 0) {
  //   if (!(*pfout).mults) {
  //     (*pfout).mults = new int[(*pf1).nroots];
  //   }
  //   (*pfout).nroots = (*pf1).nroots;
  //   (*pfout).den_deg = 0;
  //   for (int k=0; k<(*pf1).nroots; k++) {
  //     (*pfout).mults[k] = (*pf1).mults[k];
  //     (*pfout).den_deg += (*pfout).mults[k];
  //   }
  // } else if ((*pf1).den_deg == 0 && (*pf2).den_deg != 0) {
  //   if (!(*pfout).mults) {
  //     (*pfout).mults = new int[(*pf2).nroots];
  //   }
  //   (*pfout).nroots = (*pf2).nroots;
  //   (*pfout).den_deg = 0;
  //   for (int k=0; k<(*pf2).nroots; k++) {
  //     (*pfout).mults[k] = (*pf2).mults[k];
  //     (*pfout).den_deg += (*pfout).mults[k];
  //   }
  // } else if ((*pf1).den_deg == 0 && (*pf2).den_deg == 0) {
  //   cout << "both den degs are zero" << endl;
  //   (*pfout).nroots = 0;
  //   (*pfout).den_deg = 0;
  //   cout << "about to return:" << endl;
  //   poly_frac_print(pfout);
  //   cout << "nroot = " << (*pfout).nroots << endl;
  //   cout << pfout << endl;
  //   cout << "delete internal pf" << endl;
  //   if (pfout == pfin1) {
  //     cout << "delete 1st one" << endl;
  //     delete[] (*pf1).coeffs;
  //     if ((*pf1).mults) {
  //       delete[] (*pf1).mults;
  //     }
  //   }
  //   if (pfout == pfin2) {
  //     cout << "delete 2nd one" << endl;
  //     delete[] (*pf2).coeffs;
  //     if ((*pf2).mults) {
  //       delete[] (*pf2).mults;
  //     }
  //   }

  //   cout << "return" << endl;
  //   return;
  // }

  // normal
  poly_frac_normal(pfout);

  // FREE
  if (ovrwrt == 1) {
    poly_frac_free(pf1);
    delete[] pf1;
  } else if (ovrwrt == 2) {
    poly_frac_free(pf2);
    delete[] pf2;
  } else if (ovrwrt == 3) {
    poly_frac_free(pf1);
    delete[] pf1;
    poly_frac_free(pf2);
    delete[] pf2;
  }
}


void poly_frac_pow_ui_wp2(
	int wp2,
  struct poly_frac *pfout,
  struct poly_frac *pfin,
  int pow
) {
  if (pow == 0) {
    poly_frac_set_ui_wp2(wp2, pfout, 1, (*pfin).nroots);
    return;
  }

  if (pow == 1) {
    if (pfout != pfin) {
      poly_frac_set_pf_wp2(wp2, pfout, pfin);
      return;
    }
  }

  // deal with overwriting
  struct poly_frac *pf;
  int ovrwrt;
  if (pfout == pfin) {
    ovrwrt = 1;
    pf = new struct poly_frac[1];
    poly_frac_build(pf);
    poly_frac_set_pf_wp2(wp2, pf, pfin);
  } else {
    ovrwrt = 0;
    pf = pfin;
  }

  // poly_frac_set_ui(pfout, 1, (*pfin).nroots);
  poly_frac_set_pf_wp2(wp2, pfout, pf);
  for (int i=1; i<pow; i++) {
    poly_frac_mul_pf_wp2(wp2, pfout, pfout, pf);
  }

  if (ovrwrt) {
    poly_frac_free(pf);
    delete[] pf;
  }


}


void poly_frac_set_mpc_wp2(
	int wp2,
  // OUTPUT
  struct poly_frac *pf,
  // INPUT
  mpc_t *mpcin,
  int nroots
) {
  // if mpc input is zero, assign poly_frac
  if (mpc_lessthan_tol(*mpcin)) {
    poly_frac_set_zero(pf);
    // (*pf).num_vdeg = -1;
    // (*pf).coeffs = NULL;
    // (*pf).den_deg = 0;
    return;
  }

  if (nroots > 0) {
    if (!(*pf).mults) {
      (*pf).mults = new int[nroots];
    }
  }
  if ((*pf).coeffs) {
    mpc_rk1_clear((*pf).coeffs, (*pf).num_vdeg+1);
    delete[] (*pf).coeffs;
  }
  (*pf).num_deg = 0;
  (*pf).num_vdeg = 0;
  (*pf).den_deg = 0;
  (*pf).nroots = nroots;

  (*pf).coeffs = new mpc_t[1];
  mpc_init3((*pf).coeffs[0], wp2, wp2);
  mpc_set((*pf).coeffs[0], *mpcin, MPFR_RNDN);
  for (int k=0; k<nroots; k++) {
    (*pf).mults[k] = 0;
  }
}


//////
// RELATIVE ERROR FUNTIONS
//////
void rel_err_poly_mul_root(mpc_t *coeffs, int deg, mpc_t root, int wp_bin) {
    /*
    Perform the product between a polynomial

        \sum_{i=0}^{deg} a_i * x^i

    (whose coefficients are given as input) and the polynomial

        x - root

    */
    mpc_set(coeffs[deg+1], coeffs[deg], MPFR_RNDN);
    for (int k=deg; k>0; k--) {
        // printf("k = %d\n", k);
        mpc_neg(coeffs[k], coeffs[k], MPFR_RNDN);
        // mpc_fma(coeffs[k], coeffs[k], root, coeffs[k-1], MPFR_RNDN);
        // printf("coeff = "); print_mpc(&coeffs[k]); printf("\n");
        exp_rel_err_mpc_mul(coeffs[k], coeffs[k], root, MPFR_RNDN, wp_bin);
        exp_rel_err_mpc_add(coeffs[k], coeffs[k], coeffs[k-1], MPFR_RNDN, wp_bin);
    }
    // printf("k = 0\n");
    // printf("coeff = "); print_mpc(&coeffs[0]); printf("\n");
    exp_rel_err_mpc_mul(coeffs[0], coeffs[0], root, MPFR_RNDN, wp_bin);
    mpc_neg(coeffs[0], coeffs[0], MPFR_RNDN);
}


void rel_err_poly_mul_root_cp(mpc_t *out, mpc_t *coeffs, int deg, mpc_t root, int wp_bin) {
    /*
    Perform the product between a polynomial

        \sum_{i=0}^{deg} a_i * x^i

    (whose coefficients are given as input) and the polynomial

        x - root

    */

    mpc_set(out[deg+1], coeffs[deg], MPFR_RNDN);
    for (int k=deg; k>0; k--) {
        // ERROR: coeff should not change sign
        // mpc_neg(out[k], coeffs[k], MPFR_RNDN);
        // mpc_fma(out[k], out[k], root, coeffs[k-1], MPFR_RNDN);
        mpc_neg(coeffs[k], coeffs[k], MPFR_RNDN);
        // mpc_fma(out[k], coeffs[k], root, coeffs[k-1], MPFR_RNDN);
        exp_rel_err_mpc_mul(out[k], coeffs[k], root, MPFR_RNDN, wp_bin);
        exp_rel_err_mpc_add(out[k], out[k], coeffs[k-1], MPFR_RNDN, wp_bin);
    }
    exp_rel_err_mpc_mul(out[0], coeffs[0], root, MPFR_RNDN, wp_bin);
    mpc_neg(out[0], out[0], MPFR_RNDN);
}


void rel_err_mul_vpoly(
  mpc_t *out_pol, mpc_t *pol1, mpc_t *pol2,
  int deg1, int deg2,
  int vdeg1, int vdeg2,
  int ldeg1, int ldeg2,
  int wp_bin
) {
  /*
  Perform multiplication of two polynomials:

    pol1 = \sum_{i=ldeg1}^{deg1} a_i * x^i
    pol2 = \sum_{j=ldeg2}^{deg2} b_j * x^j

    pol1*pol2 = \sum_{k=ldeg1+ldeg2}^{deg1+deg2} \sum_{i=\max{0, k-deg2}}^{k} a_i * b_{k-i}
  */
  int i, k, i_min = 0, i_max, kp;
  // mpc_t tmp;
  // mpc_init3(tmp, wp2, wp2);
  
  int ldeg = ldeg1+ldeg2;
  mpc_t tmpc; mpc_init3(tmpc, wp2, wp2);
  for (k=0; k<=vdeg1+vdeg2; k++) {
    kp = k + ldeg;
    // cout << "k, kp = " << k << ", " << kp << endl;
    mpc_set_ui(out_pol[k], 0, MPFR_RNDN);
    if (ldeg1 > kp - deg2) {
      i_min = ldeg1;
    } else {
      i_min = kp - deg2;
    }
    if (deg1 < kp - ldeg2) {
      i_max = deg1;
    } else {
      i_max = kp - ldeg2;
    }
    // cout << "i_min, i_max = " << i_min << ", " << i_max << endl;
    for (i=i_min; i<=i_max; i++) {
      // cout << "i = " << i << endl;
      // mpc_fma(out_pol[k], pol1[i-ldeg1], pol2[kp-i-ldeg2], out_pol[k], MPFR_RNDN);
      // cout << "pol1[" << i-ldeg1 << "] = "; print_mpc(&pol1[i-ldeg1]); cout << endl;
      // cout << "pol2[" << kp-i-ldeg2 << "] = "; print_mpc(&pol2[kp-i-ldeg2]); cout << endl;
      // exp_rel_err_mpc_mul(
      exp_rel_err_mpc_mul(
        tmpc, pol1[i-ldeg1], pol2[kp-i-ldeg2], MPFR_RNDN
        , wp_bin
      );
      exp_rel_err_mpc_add(
        out_pol[k], out_pol[k], tmpc, MPFR_RNDN
        , wp_bin
      );
    }
  }
}


void rel_err_poly_frac_mul_mpc(
  struct poly_frac *pfout,
  struct poly_frac *pfin,
  mpc_t *mpcin,
  int wp_bin
) {
  // if input polynomial fraction is zero, set output to zero
  if ((*pfin).num_vdeg == -1) {
    if (pfout != pfin) {
      poly_frac_set_zero(pfout);
    }
    return;
  }

  // if input mpc is zero, set output to zero
  if (mpc_lessthan_tol(*mpcin)) {
    // cout << "about to set zero" << endl;
    poly_frac_set_zero(pfout);
    // if ((*pfout).coeffs) {
    //   delete[] (*pfout).coeffs;
    //   (*pfout).coeffs = NULL;
    // }
    // (*pfout).num_vdeg = -1;
    // for (int k=0; k<nroots; k++) {
    //   (*pfout).mults[k] = 0;
    // }
    return;
  }

  if (pfout != pfin) {
    if (!(*pfout).mults) {
      (*pfout).mults = new int[(*pfin).nroots];
    }
    if ((*pfout).coeffs) {
      mpc_rk1_clear((*pfout).coeffs, (*pfout).num_vdeg+1);
      delete[] (*pfout).coeffs;
    }
    poly_frac_set_info_pf(pfout, pfin);
    poly_frac_alloc_init_coeffs(pfout);
    for (int k=0; k<(*pfin).nroots; k++) {
      (*pfout).mults[k] = (*pfin).mults[k];
    }
  }

  for (int k=0; k<=(*pfin).num_vdeg; k++) {
    exp_rel_err_mpc_mul((*pfout).coeffs[k], (*pfin).coeffs[k], *mpcin, MPFR_RNDN, wp_bin);
  }
}


void rel_err_poly_frac_add_pf(
  struct poly_frac *pfout,
  struct poly_frac *pfin1,
  struct poly_frac *pfin2,
  mpc_t *roots,
  mpc_t *flop_mpc,
  int wp_bin
) {
  // check whether 1st is zero
  if ((*pfin1).num_vdeg == -1) {
    // cout << "1st is zero" << endl;
    // if the mpc exists, check if it is null
    if (flop_mpc != NULL) {
      if (mpc_lessthan_tol(*flop_mpc)) {
        // return the null polynomial fraction
        poly_frac_set_zero(pfout);
        // if ((*pfout).coeffs) {
        //   delete[] (*pfout).coeffs;
        //   (*pfout).coeffs = NULL;
        // }
        // (*pfout).num_vdeg = -1;
        return;
      } else {
        // poly_frac_mul_mpc(
        rel_err_poly_frac_mul_mpc(
          pfout, pfin2, flop_mpc
          , wp_bin
        );
        return;
      }
    }
    // return the 2nd
    if (pfout != pfin2) {
      // cout << "set output to 2nd addend" << endl;
      // poly_frac_print(pfin2, nroots);
      poly_frac_set_pf(pfout, pfin2);
      // cout << "done" << endl;
    }
    return;
  }
  // check whether 2nd is zero
  if ((*pfin2).num_vdeg == -1) {
    if (pfout != pfin1) {
      poly_frac_set_pf(pfout, pfin1);
    }
    return;
  }
  if (flop_mpc != NULL) {
    if (mpc_lessthan_tol(*flop_mpc)) {
      if (pfout != pfin1) {
        poly_frac_set_pf(pfout, pfin1);
      }
      return;
    }
  }

  int k, m, count;
  // compute LCM between denominators, and degrees
  (*pfout).nroots = (*pfin1).nroots;
  if (!(*pfout).mults) {
    (*pfout).mults = new int[(*pfin1).nroots];
  }

  // deal with overwriting
  int *mults, ovwrt_mults;
  if (pfout == pfin1 || pfout == pfin2) {
    mults = new int[(*pfin1).nroots];
    ovwrt_mults = 1;
  } else {
    mults = (*pfout).mults;
    ovwrt_mults = 0;
  }

  int num_deg1 = (*pfin1).num_deg, num_deg2 = (*pfin2).num_deg;
  int num_ldeg1 = (*pfin1).num_deg - (*pfin1).num_vdeg;
  int num_ldeg2 = (*pfin2).num_deg - (*pfin2).num_vdeg;
  (*pfout).den_deg = 0;
  int tmp_int;
  //// null root
  if (tmp_int = (*pfin1).mults[0] - (*pfin2).mults[0], tmp_int > 0) {
    mults[0] = (*pfin1).mults[0];
    num_deg2 += tmp_int;
    num_ldeg2 += tmp_int;
  } else if (tmp_int = (*pfin2).mults[0] - (*pfin1).mults[0], tmp_int > 0) {
    mults[0] = (*pfin2).mults[0];
    num_deg1 += tmp_int;
    num_ldeg1 += tmp_int;
  } else if ((*pfin1).mults[0] == (*pfin2).mults[0]) {
    mults[0] = (*pfin1).mults[0];
  }
  (*pfout).den_deg += mults[0];
  //// other roots
  for (k=1; k<(*pfin1).nroots; k++) {
    if ((*pfin1).mults[k] > (*pfin2).mults[k]) {
      mults[k] = (*pfin1).mults[k];
      num_deg2 += (*pfin1).mults[k] - (*pfin2).mults[k];
    } else if ((*pfin1).mults[k] < (*pfin2).mults[k]) {
      mults[k] = (*pfin2).mults[k];
      num_deg1 += (*pfin2).mults[k] - (*pfin1).mults[k];
    } else if ((*pfin1).mults[k] == (*pfin2).mults[k]) {
      mults[k] = (*pfin1).mults[k];
    }
    (*pfout).den_deg += mults[k];
  }
  int num_vdeg1 = num_deg1 - num_ldeg1;
  int num_vdeg2 = num_deg2 - num_ldeg2;
  (*pfout).num_deg = max(num_deg1, num_deg2);
  int num_ldeg = min(num_ldeg1, num_ldeg2);
  int num_vdeg = (*pfout).num_deg - num_ldeg;
  // cout << "num_ldeg1 = " << num_ldeg1 << endl;
  // cout << "num_deg1 = " << num_deg1 << endl;
  // cout << "num_ldeg2 = " << num_ldeg2 << endl;
  // cout << "num_deg2 = " << num_deg2 << endl;

  // compute coefficients of the numerator
  //// contribution from 1st addend
  //// (the coefficients are already stored into the output)
  mpc_t *num1 = new mpc_t[num_vdeg1+1];
  init_rk1_mpc(num1, num_vdeg1+1);
  copy_poly(num1, (*pfin1).coeffs, (*pfin1).num_vdeg);
  count = 0;
  for (k=1; k<(*pfout).nroots; k++) {
    // cout << "k = " << k << endl;
    for (m=0; m<mults[k]-(*pfin1).mults[k]; m++) {
      // cout << "m = " << m << endl;
      rel_err_poly_mul_root(
      // poly_mul_root(
        num1, (*pfin1).num_vdeg+count, roots[k]
        , wp_bin
      );
      count++;
    }
  }
  // cout << "preparato 1st" << endl;
  // print_poly(num1, num_vdeg1);
  //// contribution from 2nd addend
  mpc_t *num2 = new mpc_t[num_vdeg2+1];
  init_rk1_mpc(num2, num_vdeg2+1);
  if (flop_mpc == NULL) {
    copy_poly(num2, (*pfin2).coeffs, (*pfin2).num_vdeg);
  } else {
    for (k=0; k<=(*pfin2).num_vdeg; k++) {
      exp_rel_err_mpc_mul(
      // mpc_mul(
        num2[k], (*pfin2).coeffs[k], *flop_mpc, MPFR_RNDN
        , wp_bin
      );
    }
  }
  count = 0;
  for (k=1; k<(*pfout).nroots; k++) {
    for (m=0; m<mults[k]-(*pfin2).mults[k]; m++) {
      // cout << "k, m = " << k << ", " << m << endl;
      // cout << "root = : "; print_mpc(&roots[k]); cout << endl;
      rel_err_poly_mul_root(
      // poly_mul_root(
        num2, (*pfin2).num_vdeg+count, roots[k]
        , wp_bin
      );
      count++;
    }
  }
  // cout << "preparato 2nd" << endl;
  // print_poly(num2, num_vdeg2);

  if (ovwrt_mults) {
    for (k=0; k<(*pfout).nroots; k++) {
      (*pfout).mults[k] = mults[k];
    }
    delete[] mults;
  }
  //// sum of the two contributions
  if ((*pfout).coeffs) {
    mpc_rk1_clear((*pfout).coeffs, (*pfout).num_vdeg+1);
    delete[] (*pfout).coeffs;
  }
  (*pfout).coeffs = new mpc_t[num_vdeg+1];
  init_rk1_mpc((*pfout).coeffs, num_vdeg+1);
  int found1, found2, kp;
  for (k=num_ldeg; k<=(*pfout).num_deg; k++) {
    kp = k - num_ldeg;
    // cout << "k, kp = " << k << ", " << kp << endl;
    if (k >= num_ldeg1 && k <= num_deg1) {
      found1 = 1;
    } else {
      found1 = 0;
    }
    if (k >= num_ldeg2 && k <= num_deg2) {
      found2 = 1;
    } else {
      found2 = 0;
    }

    if (found1 == 0 && found2 == 0) {
      mpc_set_ui((*pfout).coeffs[kp], 0, MPFR_RNDN);
    } else if (found1 && found2 == 0) {
      mpc_set((*pfout).coeffs[kp], num1[k-num_ldeg1], MPFR_RNDN);
    } else if (found1 == 0 && found2) {
      mpc_set((*pfout).coeffs[kp], num2[k-num_ldeg2], MPFR_RNDN);
    } else if (found1 && found2) {
      // mpc_add(
      exp_rel_err_mpc_add(
        (*pfout).coeffs[kp], num1[k-num_ldeg1], num2[k-num_ldeg2], MPFR_RNDN
        , wp_bin
      );
    }

    // cout << "found1, found2 = " << found1 << ", " << found2 << endl;
    // print_mpc(&(*pfout).coeffs[kp]);
    // cout << endl;
  }

  // lower padding
  int lpad = 0, upad = num_vdeg;
  for (k=0; k<=num_vdeg; k++) {
    // if(mpc_lessthan_tol((*pfout).coeffs[k])) {
    if(mpc_zero_p((*pfout).coeffs[k])) {
      lpad++;
      continue;
    } else {
      break;
    }
  }
  // cout << "lpad = " << lpad << endl;
  if (lpad == num_vdeg + 1) {
    mpc_rk1_clear((*pfout).coeffs, num_vdeg+1);
    delete[] (*pfout).coeffs;
    (*pfout).coeffs = NULL;
    (*pfout).num_vdeg = -1;
    (*pfout).den_deg = 0;
    return;
  }

  // upper padding
  for (k=num_vdeg; k>=0; k--) {
    // if(mpc_lessthan_tol((*pfout).coeffs[k])) {
    if(mpc_zero_p((*pfout).coeffs[k])) {
      upad--;
      continue;
    } else {
      break;
    }
  }
  // cout << "upad = " << upad << endl;

  // prune null coefficients
  if (lpad>0 || upad<num_vdeg1) {
    int old_vdeg = num_vdeg;
    mpc_t *unpruned_coeffs;
    unpruned_coeffs = (*pfout).coeffs;
    (*pfout).num_deg -= num_vdeg - upad;
    num_vdeg = upad-lpad;
    (*pfout).coeffs = new mpc_t[num_vdeg+1];
    for (k=0; k<=num_vdeg; k++) {
      mpc_init3((*pfout).coeffs[k], wp2, wp2);
      mpc_set((*pfout).coeffs[k], unpruned_coeffs[lpad+k], MPFR_RNDN);
    }
    mpc_rk1_clear(unpruned_coeffs, old_vdeg+1);
    delete[] unpruned_coeffs;
  }

  // (in case it is really needed)
  // simplify roots from numerator and LCM

  // (*pfout).coeffs = new mpc_t[num_vdeg+1];
  // for (k=0; k<=num_vdeg; k++) {
  //   mpc_init3((*pfout).coeffs[k], wp2, wp2);
  //   mpc_set((*pfout).coeffs[k], num[k], MPFR_RNDN);
  // }
  (*pfout).num_vdeg = num_vdeg;
  // (*pfout).num_deg = num_deg;

  // normal
  poly_frac_normal(pfout);

  mpc_rk1_clear(num1, num_vdeg1+1); delete[] num1;
  mpc_rk1_clear(num2, num_vdeg2+1); delete[] num2;
}


void rel_err_poly_frac_mul_pf(
  struct poly_frac *pfout,
  struct poly_frac *pfin1,
  struct poly_frac *pfin2,
  int wp_bin
) {
  // if input is zero, set output to zero
  if ((*pfin1).num_vdeg == -1) {
    if (pfout != pfin1) {
      poly_frac_set_zero(pfout);
    }
    return;
  }
  if ((*pfin2).num_vdeg == -1) {
    if (pfout != pfin2) {
      poly_frac_set_zero(pfout);
    }
    return;
  }

  // deal with overwriting
  // cout << "deal with overwriting" << endl;
  struct poly_frac *pf1, *pf2;
  int ovrwrt;
  if (pfout == pfin1 && pfout != pfin2) {
    ovrwrt = 1;
    pf1 = new struct poly_frac[1];
    poly_frac_build(pf1);
    poly_frac_set_pf(pf1, pfin1);
    pf2 = pfin2;
  } else if (pfout != pfin1 && pfout == pfin2) {
    ovrwrt = 2;
    pf2 = new struct poly_frac[1];
    poly_frac_build(pf2);
    poly_frac_set_pf(pf2, pfin2);
    pf1 = pfin1;
  } else if (pfout == pfin1 && pfout == pfin2) {
    ovrwrt = 3;
    pf1 = new struct poly_frac[1];
    pf2 = new struct poly_frac[1];
    poly_frac_build(pf1);
    poly_frac_set_pf(pf1, pfin1);
    poly_frac_build(pf2);
    poly_frac_set_pf(pf2, pfin2);
  } else {
    ovrwrt = 0;
    pf1 = pfin1;
    pf2 = pfin2;
  }

  // multiply numerators
  // cout << "multiply numerators" << endl;
  int ldeg1 = (*pf1).num_deg - (*pf1).num_vdeg;
  int ldeg2 = (*pf2).num_deg - (*pf2).num_vdeg;
  if ((*pfout).coeffs) {
    mpc_rk1_clear((*pfout).coeffs, (*pfout).num_vdeg+1);
    delete[] (*pfout).coeffs;
  }
  (*pfout).num_deg = (*pfin1).num_deg + (*pfin2).num_deg;
  (*pfout).num_vdeg =  (*pfout).num_deg - ldeg1 - ldeg2;
  (*pfout).coeffs = new mpc_t[(*pfout).num_vdeg+1];
  init_rk1_mpc((*pfout).coeffs, (*pfout).num_vdeg+1);
  // cout << "num_deg = " << (*pfout).num_deg << endl;
  // cout << "num_vdeg = " << (*pfout).num_vdeg << endl;
  rel_err_mul_vpoly(
    (*pfout).coeffs, (*pf1).coeffs, (*pf2).coeffs,
    (*pf1).num_deg, (*pf2).num_deg,
    (*pf1).num_vdeg, (*pf2).num_vdeg,
    ldeg1, ldeg2,
    wp_bin
  );
  // cout << "numerator:" << endl;
  // print_poly((*pfout).coeffs, (*pfout).num_vdeg);

  // multiply denominators
  // cout << "multiply denominators" << endl;
  // cout << "den degs: " << (*pf1).den_deg << ", " << (*pf2).den_deg << endl;
  if (!(*pfout).mults) {
    (*pfout).mults = new int[(*pf1).nroots];
  }
  (*pfout).nroots = (*pf1).nroots;
  (*pfout).den_deg = 0;
  for (int k=0; k<(*pf1).nroots; k++) {
    (*pfout).mults[k] = (*pf1).mults[k] + (*pf2).mults[k];
    (*pfout).den_deg += (*pfout).mults[k];
  }
  // if ((*pf1).den_deg != 0 && (*pf2).den_deg != 0) {
  //   if (!(*pfout).mults) {
  //     (*pfout).mults = new int[(*pf1).nroots];
  //   }
  //   (*pfout).nroots = (*pf1).nroots;
  //   (*pfout).den_deg = 0;
  //   for (int k=0; k<(*pf1).nroots; k++) {
  //     (*pfout).mults[k] = (*pf1).mults[k] + (*pf2).mults[k];
  //     (*pfout).den_deg += (*pfout).mults[k];
  //   }
  // } else if ((*pf1).den_deg != 0 && (*pf2).den_deg == 0) {
  //   if (!(*pfout).mults) {
  //     (*pfout).mults = new int[(*pf1).nroots];
  //   }
  //   (*pfout).nroots = (*pf1).nroots;
  //   (*pfout).den_deg = 0;
  //   for (int k=0; k<(*pf1).nroots; k++) {
  //     (*pfout).mults[k] = (*pf1).mults[k];
  //     (*pfout).den_deg += (*pfout).mults[k];
  //   }
  // } else if ((*pf1).den_deg == 0 && (*pf2).den_deg != 0) {
  //   if (!(*pfout).mults) {
  //     (*pfout).mults = new int[(*pf2).nroots];
  //   }
  //   (*pfout).nroots = (*pf2).nroots;
  //   (*pfout).den_deg = 0;
  //   for (int k=0; k<(*pf2).nroots; k++) {
  //     (*pfout).mults[k] = (*pf2).mults[k];
  //     (*pfout).den_deg += (*pfout).mults[k];
  //   }
  // } else if ((*pf1).den_deg == 0 && (*pf2).den_deg == 0) {
  //   cout << "both den degs are zero" << endl;
  //   (*pfout).nroots = 0;
  //   (*pfout).den_deg = 0;
  //   cout << "about to return:" << endl;
  //   poly_frac_print(pfout);
  //   cout << "nroot = " << (*pfout).nroots << endl;
  //   cout << pfout << endl;
  //   cout << "delete internal pf" << endl;
  //   if (pfout == pfin1) {
  //     cout << "delete 1st one" << endl;
  //     delete[] (*pf1).coeffs;
  //     if ((*pf1).mults) {
  //       delete[] (*pf1).mults;
  //     }
  //   }
  //   if (pfout == pfin2) {
  //     cout << "delete 2nd one" << endl;
  //     delete[] (*pf2).coeffs;
  //     if ((*pf2).mults) {
  //       delete[] (*pf2).mults;
  //     }
  //   }

  //   cout << "return" << endl;
  //   return;
  // }

  // normal
  poly_frac_normal(pfout);

  if (ovrwrt == 1) {
    poly_frac_free(pf1);
    delete[] pf1;
  } else if (ovrwrt == 2) {
    poly_frac_free(pf2);
    delete[] pf2;
  } else if (ovrwrt == 3) {
    poly_frac_free(pf1);
    delete[] pf1;
    poly_frac_free(pf2);
    delete[] pf2;
  }
}


void rel_err_poly_frac_rk2_mul_pf_rk1(
  struct poly_frac *pf_tens_out,
  struct poly_frac **pf_tens_in1,
  struct poly_frac *pf_tens_in2,
  int dim1, int dim2,
  mpc_t *roots,
  int wp_bin
) {
  struct poly_frac **pf_tens1, *pf_tens2;
  int to_be_del;
  if (pf_tens_out == pf_tens_in2) {
    to_be_del = 2;
    // copy input
    pf_tens2 = new struct poly_frac[dim2];
    for (int i2=0; i2<dim2; i2++) {
      poly_frac_build(&pf_tens2[i2]);
      poly_frac_set_pf(&pf_tens2[i2], &pf_tens_in2[i2]);
    }
    pf_tens1 = pf_tens_in1;
  } else {
    to_be_del = 0;
    pf_tens1 = pf_tens_in1;
    pf_tens2 = pf_tens_in2;
  }

  // row by column product
  struct poly_frac pf_tmp;
  poly_frac_build(&pf_tmp);
  for (int i1=0; i1<dim1; i1++) {
    rel_err_poly_frac_mul_pf(
      &pf_tens_out[i1], &pf_tens1[i1][0], &pf_tens2[0],
      wp_bin
    );
    for (int i2=1; i2<dim2; i2++) {
      // cout << "i1, i2 = " << i1 << ", " << i2 << endl;
      rel_err_poly_frac_mul_pf(
        &pf_tmp, &pf_tens1[i1][i2], &pf_tens2[i2],
        wp_bin
      );
      rel_err_poly_frac_add_pf(
        &pf_tens_out[i1], &pf_tens_out[i1], &pf_tmp,
        roots,
        NULL,
        wp_bin
      );
    }
  }

  if (to_be_del == 2) {
    // cout << "trying to delete" << endl;
    del_rk1_poly_frac(pf_tens2, dim2);
    // cout << "done" << endl;
  }
}


void rel_err_poly_frac_rk2_mul_pf_rk2(
  struct poly_frac **pf_tens_out,
  struct poly_frac **pf_tens_in1,
  struct poly_frac **pf_tens_in2,
  int dim1, int dim, int dim2,
  mpc_t *roots,
  int wp_bin
) {
  struct poly_frac **pf_tens1, **pf_tens2;
  int to_be_del;
  if (pf_tens_out == pf_tens_in1) {
    to_be_del = 1;
    // copy input
    malloc_rk2_tens(pf_tens1, dim1, dim);
    for (int i1=0; i1<dim1; i1++) {
      for (int i2=0; i2<dim; i2++) {
        poly_frac_build(&pf_tens1[i1][i2]);
        poly_frac_set_pf(&pf_tens1[i1][i2], &pf_tens_in1[i1][i2]);
      }
    }
    pf_tens2 = pf_tens_in2;
  } else if (pf_tens_out == pf_tens_in2) {
    to_be_del = 2;
    // copy input
    malloc_rk2_tens(pf_tens2, dim, dim2);
    for (int i1=0; i1<dim; i1++) {
      for (int i2=0; i2<dim2; i2++) {
        poly_frac_build(&pf_tens2[i1][i2]);
        poly_frac_set_pf(&pf_tens2[i1][i2], &pf_tens_in2[i1][i2]);
      }
    }
    pf_tens1 = pf_tens_in1;
  } else {
    to_be_del = 0;
    pf_tens1 = pf_tens_in1;
    pf_tens2 = pf_tens_in2;
  }
  
  // row by column product
  struct poly_frac pf_tmp;
  poly_frac_build(&pf_tmp);
  for (int i1=0; i1<dim1; i1++) {
    for (int i2=0; i2<dim2; i2++) {
      // cout << "i, j = " << i1 << ", " << i2 << endl;
      // cout << "1st factor" << endl;
      // poly_frac_print(&pf_tens1[i1][0]);
      // cout << "2nd factor" << endl;
      // poly_frac_print(&pf_tens2[0][i2]);
      rel_err_poly_frac_mul_pf(
        &pf_tens_out[i1][i2], &pf_tens1[i1][0], &pf_tens2[0][i2],
        wp_bin
      );
      // cout << "result" << endl;
      // poly_frac_print(&pf_tens_out[i1][i2]);
      for (int i3=1; i3<dim; i3++) {
        // cout << "k = " << i3 << endl;
        // cout << "1st factor" << endl;
        // poly_frac_print(&pf_tens1[i1][i3]);
        // cout << "2nd factor" << endl;
        // poly_frac_print(&pf_tens2[i3][i2]);
        rel_err_poly_frac_mul_pf(
          &pf_tmp, &pf_tens1[i1][i3], &pf_tens2[i3][i2],
          wp_bin
        );
        // cout << "result" << endl;
        // poly_frac_print(&pf_tmp);
        rel_err_poly_frac_add_pf(
          &pf_tens_out[i1][i2], &pf_tens_out[i1][i2], &pf_tmp,
          roots,
          NULL,
          wp_bin
        );
        // cout << "result" << endl;
        // poly_frac_print(&pf_tens_out[i1][i2]);
      }
    }
  }

  if (to_be_del == 1) {
    // cout << "trying to delete" << endl;
    poly_frac_rk2_free(pf_tens1, dim1, dim);
    // cout << "done" << endl;
  } else if (to_be_del == 2) {
    // cout << "trying to delete" << endl;
    poly_frac_rk2_free(pf_tens2, dim, dim2);
    // cout << "done" << endl;
  }
}


void rel_err_poly_frac_rk2_mul_mpc_rk2(
  struct poly_frac **pf_tens_out,
  struct poly_frac **pf_tens_in,
  mpc_t **mpc_tens,
  int dim,
  mpc_t *roots,
  int wp_bin
) {
  struct poly_frac **pf_tens;
  int to_be_del;
  if (pf_tens_out == pf_tens_in) {
    to_be_del = 1;
    // copy input
    malloc_rk2_tens(pf_tens, dim, dim);
    for (int i1=0; i1<dim; i1++) {
      for (int i2=0; i2<dim; i2++) {
        poly_frac_build(&pf_tens[i1][i2]);
        poly_frac_set_pf(&pf_tens[i1][i2], &pf_tens_in[i1][i2]);
      }
    }
  } else {
    to_be_del = 0;
    pf_tens = pf_tens_in;
  }
  
  // row by column product
  for (int i1=0; i1<dim; i1++) {
    for (int i2=0; i2<dim; i2++) {
      // cout << "i, j = " << i1 << ", " << i2 << endl;
      // cout << "1st term:" << endl;
      // poly_frac_print(&pf_tens_in[i1][0]);
      // cout << "2nd term" << endl;
      // print_mpc(&mpc_tens[0][i2]);
      // cout << endl;
      rel_err_poly_frac_mul_mpc(
      // poly_frac_mul_mpc(
        &pf_tens_out[i1][i2], &pf_tens[i1][0], &mpc_tens[0][i2]
        , wp_bin
      );
      // cout << "result" << endl;
      // poly_frac_print(&pf_tens_out[i1][i2]);
      // cout << "mul" << endl;
      for (int i3=1; i3<dim; i3++) {
        // cout << "k = " << i3 << endl;
        // cout << "1st addend" << endl;
        // poly_frac_print(&pf_tens_out[i1][i2]);
        // cout << "2nd addend" << endl;
        // poly_frac_print(&pf_tens[i1][i3]);
        // cout << "times" << endl;
        // print_mpc(&mpc_tens[i3][i2]);
        // cout << endl;
        rel_err_poly_frac_add_pf(
          &pf_tens_out[i1][i2], &pf_tens_out[i1][i2], &pf_tens[i1][i3],
          roots,
          &mpc_tens[i3][i2],
          wp_bin
        );
        // cout << "result" << endl;
        // poly_frac_print(&pf_tens_out[i1][i2]);
      }
      // cout << "flop" << endl;
    }
  }

  if (to_be_del) {
    // cout << "trying to delete" << endl;
    poly_frac_rk2_free(pf_tens, dim, dim);
    // cout << "done" << endl;
  }
}


void rel_err_mpc_rk2_mul_poly_frac_rk2(
  struct poly_frac **pf_tens_out,
  mpc_t **mpc_tens,
  struct poly_frac **pf_tens_in,
  int dim,
  mpc_t *roots,
  int wp_bin
) {
  struct poly_frac **pf_tens;
  int to_be_del;
  if (pf_tens_out == pf_tens_in) {
    to_be_del = 1;
    // copy input
    malloc_rk2_tens(pf_tens, dim, dim);
    for (int i1=0; i1<dim; i1++) {
      for (int i2=0; i2<dim; i2++) {
        poly_frac_build(&pf_tens[i1][i2]);
        poly_frac_set_pf(&pf_tens[i1][i2], &pf_tens_in[i1][i2]);
      }
    }
  } else {
    to_be_del = 0;
    pf_tens = pf_tens_in;
  }

  // row by column product
  for (int i1=0; i1<dim; i1++) {
    for (int i2=0; i2<dim; i2++) {
      // pf_tens_out[i1][i2].mults = new int[nroots];
      rel_err_poly_frac_mul_mpc(
      // poly_frac_mul_mpc(
        &pf_tens_out[i1][i2], &pf_tens[0][i2], &mpc_tens[i1][0]
        , wp_bin
      );
      for (int i3=1; i3<dim; i3++) {
        rel_err_poly_frac_add_pf(
          &pf_tens_out[i1][i2], &pf_tens_out[i1][i2], &pf_tens[i3][i2],
          roots,
          &mpc_tens[i1][i3],
          wp_bin
        );
      }
    }
  }

  if (to_be_del) {
    poly_frac_rk2_free(pf_tens, dim, dim);
  }
}


void rel_err_poly_frac_rk2_add_pf_rk2(
  struct poly_frac **pf_tens_out,
  struct poly_frac **pf_tens_in1,
  struct poly_frac **pf_tens_in2,
  mpc_t *roots, mpc_t *flop_mpc,
  int dim1, int dim2,
  int wp_bin
) {
  for (int i1=0; i1<dim1; i1++) {
    for (int i2=0; i2<dim2; i2++) {
      rel_err_poly_frac_add_pf(
        &pf_tens_out[i1][i2], &pf_tens_in1[i1][i2], &pf_tens_in2[i1][i2],
        roots, flop_mpc,
        wp_bin
      );
    }
  }
}


void rel_err_poly_eval_sym_zero_by_roots_wo_0(
  // OUTPUT
  mpc_t *value,
  // INPUT,
  mpc_t *roots, int *mults, int nroots,
  int wp_bin
) {
  // #2BD: decide if it is worth using multiplicites and mpc_pow_ui
  //    pros: better than mul precision-wise;
  //    cons: need to store into temporary variable even when mults[k]
  //      equals one;
  // a similar argument can be made for the sign variable (swapping a bit
  // at every step vs considering if mults[k] is odd or even)
  mpc_set_ui(*value, 1, MPFR_RNDN);
  int sign = 1;
  for (int k=1; k<nroots; k++) {
    for (int m=0; m<mults[k]; m++) {
      // cout << "k, m = " << k << ", " << m << endl;
      exp_rel_err_mpc_mul(*value, *value, roots[k], MPFR_RNDN, wp_bin);
      sign = -sign;
    }
  }
  if (sign < 0) {
    mpc_neg(*value, *value, MPFR_RNDN);
  }
}


void rel_err_poly_frac_extract_LO(
  // OUTPUT
  mpfr_t *out_re, mpfr_t *out_im,
  // INPUT
  struct poly_frac *pf, mpc_t *roots, int pow,
  int wp_bin
) {
  // if input is zero, return zero
  if ((*pf).num_vdeg == -1) {
    mpfr_set_ui(*out_re, 0, MPFR_RNDN);
    mpfr_set_ui(*out_im, 0, MPFR_RNDN);
    return;
  }
  int pow_behav = (*pf).num_deg - (*pf).num_vdeg - (*pf).mults[0] + pow;
  // cout << "pow_behav = " << pow_behav << endl;
  if (pow_behav < 0) {
    fprintf(stderr, "Trying to evaluate in zero singular fraction of polynomials! Exiting...\n");
    printf("poly_frac:\n");
    poly_frac_print(pf);
    exit(EXIT_FAILURE);
  } else if (pow_behav > 0) {
    mpfr_set_ui(*out_re, 0, MPFR_RNDN);
    mpfr_set_ui(*out_im, 0, MPFR_RNDN);
  } else if (pow_behav == 0) {
    mpc_t tmp;
    mpc_init3(tmp, wp2, wp2);
    // cout << "evaluate denominator in zero" << endl;
    rel_err_poly_eval_sym_zero_by_roots_wo_0(
      &tmp,
      roots, (*pf).mults, (*pf).nroots,
      wp_bin
    );
    // mpc_div(
    exp_rel_err_mpc_div(
      tmp, (*pf).coeffs[0], tmp, MPFR_RNDN
      , wp_bin
    );
    mpc_real(*out_re, tmp, MPFR_RNDN);
    mpc_imag(*out_im, tmp, MPFR_RNDN);
  }
}


void rel_err_poly_frac_mat_extract_LO(
  // OUTPUT
  mpfr_t *LOre, mpfr_t *LOim,
  // INPUT
  struct poly_frac  **pfmat, mpc_t *roots, int pow,
  int dim1, int dim2,
  int wp_bin
) {
  for (int i1=0; i1<dim1; i1++) {
    for (int i2=0; i2<dim2; i2++) {
      // cout << "i1, i2 = " << i1 << ", " << i2 << endl;
      // poly_frac_print(&pfmat[i1][i2]);
      rel_err_poly_frac_extract_LO(
        &LOre[i1 + dim1*i2], &LOim[i1 + dim1*i2],
        &pfmat[i1][i2], roots, pow,
        wp_bin
      );
      continue;  // comment to check for missed cancellation

      if (
        mpfr_get_exp(LOre[i1+dim1*i2]) < -140 && !mpfr_zero_p(LOre[i1+dim1*i2]) ||\
        mpfr_get_exp(LOim[i1+dim1*i2]) < -140 && !mpfr_zero_p(LOim[i1+dim1*i2])
      ) {
        cout << "i1, i2 = " << i1 << ", " << i2 << endl;
        cout << "pf:" << endl;
        poly_frac_print(&pfmat[i1][i2]);
        cout << "LOre:" << endl;
        mpfr_out_str(stdout, 10, 0, LOre[i1+dim1*i2], MPFR_RNDN); cout << endl;
        cout << "LOim:" << endl;
        mpfr_out_str(stdout, 10, 0, LOim[i1+dim1*i2], MPFR_RNDN); cout << endl;
        fflush(stdout);
        perror("LO extraction: cancellation missed");
        exit(1);
      }
    }
  }
}


void rel_err_poly_coeff1_by_roots_wo_0(
  // OUPUT
  mpc_t *value,
  // INPUT,
  mpc_t *roots, int *mults, int nroots,
  int wp_bin
) {
  mpc_set_ui(*value, 0, MPFR_RNDN);
  mpc_t tmp0;
  mpc_init3(tmp0, wp2, wp2);
  mpc_t tmp1;
  mpc_init3(tmp1, wp2, wp2);
  for (int k=1; k<nroots; k++) {
    if (mults[k] == 0) {
      continue;
    }
    mpc_pow_ui(tmp0, roots[k], mults[k]-1, MPFR_RNDN);
    if (mults[k]%2 == 0) {
      mpc_neg(tmp0, tmp0, MPFR_RNDN);
    }
    for (int kp=1; kp<nroots; kp++) {
      if (kp == k || mults[kp] == 0) {
        continue;
      }
      mpc_pow_ui(tmp1, roots[kp], mults[kp], MPFR_RNDN);
      exp_rel_err_mpc_mul(tmp0, tmp0, tmp1, MPFR_RNDN, wp_bin);
      if (mults[k]%2 == 1) {
        mpc_neg(tmp0, tmp0, MPFR_RNDN);
      }
    }
    mpc_mul_ui(tmp0, tmp0, mults[k], MPFR_RNDN);
    exp_rel_err_mpc_add(*value, *value, tmp0, MPFR_RNDN, wp_bin);
  }
}


void rel_err_poly_frac_extract_NLO(
  // OUTPUT
  mpc_t *NLO,
  // INPUT
  mpc_t *LO, struct poly_frac *pf, mpc_t *roots, int pow,
  int wp_bin
) {
  if ((*pf).num_vdeg == -1) {
    mpc_set_ui(*NLO, 0, MPFR_RNDN);
    return;
  }
  
  int pow_behav = (*pf).num_deg - (*pf).num_vdeg - (*pf).mults[0] + pow;
  if (pow_behav < 0) {
    fprintf(stderr, "Trying to evaluate in zero singular fraction of polynomial! Exiting...\n");
    printf("poly_frac:\n");
    poly_frac_print(pf);
    exit(EXIT_FAILURE);
  } else if (pow_behav > 1) {
    mpc_set_ui(*NLO, 0, MPFR_RNDN);
  } else if (pow_behav == 1) {
    rel_err_poly_eval_sym_zero_by_roots_wo_0(
      NLO,
      roots, (*pf).mults, (*pf).nroots,
      wp_bin
    );
    // mpc_div(
    exp_rel_err_mpc_div(
      *NLO, (*pf).coeffs[0], *NLO, MPFR_RNDN
      , wp_bin
    );
  } else if (pow_behav == 0) {
    rel_err_poly_coeff1_by_roots_wo_0(
      NLO,
      roots, (*pf).mults, (*pf).nroots
      , wp_bin
    );
    exp_rel_err_mpc_mul(*NLO, *NLO, *LO, MPFR_RNDN, wp_bin);
    if ((*pf).num_deg > 0) {
      mpc_neg(*NLO, *NLO, MPFR_RNDN);
      exp_rel_err_mpc_add(*NLO, (*pf).coeffs[1], *NLO, MPFR_RNDN, wp_bin);
    } else {
      mpc_neg(*NLO, *NLO, MPFR_RNDN);
    }
    mpc_t tmp;
    mpc_init3(tmp, wp2, wp2);
    rel_err_poly_eval_sym_zero_by_roots_wo_0(
      &tmp,
      roots, (*pf).mults, (*pf).nroots,
      wp_bin
    );
    // mpc_div(
    exp_rel_err_mpc_div(
      *NLO, *NLO, tmp, MPFR_RNDN
      , wp_bin
    );
    
    // FREE
    mpc_clear(tmp);
  }
}


void rel_err_poly_frac_mat_extract_NLOc(
  // OUTPUT
  mpc_t **NLO,
  // INPUT
  mpfr_t *LOre, mpfr_t *LOim,
  struct poly_frac  **pfmat, mpc_t *roots, int pow,
  int dim1, int dim2,
  int wp_bin
) {
  int index;
  mpc_t tmp;
  mpc_init3(tmp, wp2, wp2);
  for (int i1=0; i1<dim1; i1++) {
    for (int i2=0; i2<dim2; i2++) {
      index =  i1 + dim1*i2;
      mpc_set_fr_fr(tmp, LOre[index], LOim[index], MPFR_RNDN);
      rel_err_poly_frac_extract_NLO(
        &NLO[i1][i2],
        &tmp, &pfmat[i1][i2], roots, pow,
        wp_bin
      );
      continue;  // comment to check for missed cancellation
      
      if (
        mpfr_get_exp(mpc_imagref(NLO[i1][i2])) < -140 && !mpfr_zero_p(mpc_imagref(NLO[i1][i2])) ||\
        mpfr_get_exp(mpc_realref(NLO[i1][i2])) < -140 && !mpfr_zero_p(mpc_realref(NLO[i1][i2]))
      ) {
        cout << "i1, i2 = " << i1 << ", " << i2 << endl;
        cout << "pf:" << endl;
        poly_frac_print(&pfmat[i1][i2]);
        cout << "NLO:" << endl;
        print_mpc(&NLO[i1][i2]); cout << endl;
        fflush(stdout);
        perror("NLO extraction: cancellation missed");
        exit(1);
      }
    }
  }
}


int rel_err_poly_frac_eval_sym_zero(
  // OUTPUT
  mpc_t *value,
  // INPUT
  struct poly_frac *pf, mpc_t *roots,
  int wp_bin
) {
  // if input is zero, return zero
  if ((*pf).num_vdeg == -1) {
    mpc_set_ui(*value, 0, MPFR_RNDN);
    return 0;
  }
  int pow_behav = (*pf).num_deg - (*pf).num_vdeg - (*pf).mults[0];
  // cout << "num_deg = " << (*pf).num_deg << endl;
  // cout << "num_vdeg = " << (*pf).num_vdeg << endl;
  // cout << "mults[0] = " << (*pf).mults[0] << endl;
  // cout << "pow_behav" << pow_behav << endl;
  if (pow_behav < 0) {
    fprintf(stderr, "Trying to evaluate in zero singular fraction of polynomials! Exiting...\n");
    printf("poly_frac:\n");
    poly_frac_print(pf);
    // return -1;
    exit(EXIT_FAILURE);
    // return -1;
  } else if (pow_behav > 0) {
    mpc_set_ui(*value, 0, MPFR_RNDN);
    return 0;
  } else if (pow_behav == 0) {
    rel_err_poly_eval_sym_zero_by_roots_wo_0(
      value,
      roots, (*pf).mults, (*pf).nroots,
      wp_bin
    );
    exp_rel_err_mpc_div(*value, (*pf).coeffs[0], *value, MPFR_RNDN, wp_bin);
    return 1;
  }
  return -1;
}


void rel_err_poly_frac_extract_oder(
  // OUTPUT
  mpc_t *out,
  // INPUT
  struct poly_frac *pf, mpc_t *roots, int nroots, int order,
  int wp_bin
) {
  int pow_behav = poly_frac_pow_behav(pf);
  struct poly_frac pfin;
  poly_frac_build(&pfin);
  poly_frac_set_pf(&pfin, pf);
  struct poly_frac pflead;
  poly_frac_build(&pflead);
  poly_frac_mul_sym_pow(&pfin, &pfin, -pow_behav);
  // poly_frac_print(&pfin);

  for (int p=pow_behav; p<=order; p++) {
    rel_err_poly_frac_eval_sym_zero(out, &pfin, roots, wp_bin);
    if (p<order) {
      mpc_neg(*out, *out, MPFR_RNDN);    
      poly_frac_set_mpc(&pflead, out, nroots);
      rel_err_poly_frac_add_pf(&pfin, &pfin, &pflead, roots, NULL, wp_bin);
      poly_frac_mul_sym_pow(&pfin, &pfin, -1);
    }
  }

}


void rel_err_poly_frac_arg_inv(
  // OUTPUT
  struct poly_frac *pfout,
  // INPUT
  struct poly_frac *pfin, mpc_t *roots,
  int wp_bin
) {
  // cout << "pf in:" << endl;
  // poly_frac_print(pfin);
  // cout << "nroots = " << (*pfin).nroots << endl;

  if ((*pfin).num_vdeg == -1) {
    poly_frac_set_zero(pfout);
    return;
  }
  
  // MULTIPLICITIES (BESIDES NULL ROOT)    
  if ((*pfout).mults) {
    if ((*pfout).nroots != (*pfin).nroots) {
      delete[] (*pfout).mults;
      (*pfout).mults = new int[(*pfin).nroots];
    }
  } else {
    (*pfout).mults = new int[(*pfin).nroots];
  }
  for (int k=1; k<(*pfin).nroots; k++) {
    // cout << "k = " << k << endl;
    // cout << "pfin.mult = " << (*pfin).mults[k] << endl;
    (*pfout).mults[k] = (*pfin).mults[k];
  }
  (*pfout).nroots = (*pfin).nroots;

  // DECIDE DEGREES AND NULL ROOT MULTIPLICITY
  if ((*pfout).coeffs) {
    if ((*pfout).num_vdeg != (*pfin).num_vdeg) {
      mpc_rk1_clear((*pfout).coeffs, (*pfout).num_vdeg+1);
      delete[] (*pfout).coeffs;
      (*pfout).coeffs = new mpc_t[(*pfin).num_vdeg+1];
      init_rk1_mpc((*pfout).coeffs, (*pfin).num_vdeg+1);
    }
  } else {
    (*pfout).coeffs = new mpc_t[(*pfin).num_vdeg+1];
    init_rk1_mpc((*pfout).coeffs, (*pfin).num_vdeg+1);
  }
  (*pfout).num_vdeg = (*pfin).num_vdeg;

  int pow_behav = (*pfin).den_deg - (*pfin).num_deg;  
  // cout << "pow behav = " << pow_behav << endl;
  if (pow_behav >= 0) {
    (*pfout).num_deg = (*pfin).num_vdeg + pow_behav;
    (*pfout).mults[0] = 0;
    (*pfout).den_deg = (*pfin).den_deg - (*pfin).mults[0];
  } else if (pow_behav < 0) {
    (*pfout).num_deg = (*pfin).num_vdeg;
    (*pfout).mults[0] = -pow_behav;
    (*pfout).den_deg = (*pfin).den_deg - (*pfin).mults[0] - pow_behav;
  }

  // NORMALIZATION FACTOR
  mpc_t norm;
  mpc_init3(norm, wp2, wp2);
  mpc_set_ui(norm, 1, MPFR_RNDN);
  for (int k=1; k<(*pfin).nroots; k++) {
    for (int m=0; m<(*pfin).mults[k]; m++) {      
      // cout << "norm = "; print_mpc(&norm); cout << endl;
      // cout << "root = "; print_mpc(&roots[k]); cout << endl;
      exp_rel_err_mpc_mul(
      // mpc_mul(
        norm, norm, roots[k], MPFR_RNDN
        , wp_bin
      );
    }
  }
  if (((*pfin).den_deg - (*pfin).mults[0])%2 != 0) {
    mpc_neg(norm, norm, MPFR_RNDN);
  }
  // cout << "final norm = "; print_mpc(&norm); cout << endl;

  // NUMERATOR (revert order and normalize)

  for (int k=0; k<=(*pfout).num_vdeg; k++) {
    // cout << "k = " << k << endl;
    // cout << "coeff = "; print_mpc(&(*pfin).coeffs[(*pfout).num_vdeg-k]); cout << endl;
    // cout << "norm = "; print_mpc(&norm); cout << endl;
    // mpc_div(
    exp_rel_err_mpc_div(
      (*pfout).coeffs[k], (*pfin).coeffs[(*pfout).num_vdeg-k], norm, MPFR_RNDN
      , wp_bin
    );
    // cout << "division = "; print_mpc(&(*pfout).coeffs[k]); cout << endl;
  }
  
}


void rel_err_poly_frac_rk2_infty(
  // OUTPUT
  struct poly_frac **pfout, mpc_t *out_roots, int *out_zero_label,
  // INPUT
  struct poly_frac **pfin, mpc_t *in_roots,
  int dim1, int dim2, int nroots,
  int wp_bin
) {
  // invert roots
  mpc_set_ui(out_roots[0], 0, MPFR_RNDN);
  for (int k=1; k<nroots; k++) {
    mpc_pow_si(out_roots[k], in_roots[k], -1, MPFR_RNDN);
  }

  // invert matrix elements and set new zero_label
  // *out_zero_label = -1;
  for (int i1=0; i1<dim1; i1++) {
    for (int i2=0; i2<dim2; i2++) {
      // cout << "i1, i2 = " << i1 << ", " << i2 << endl;
      // if (i1 == 137 && i2 == 0) dbg = 1;
      rel_err_poly_frac_arg_inv(&pfout[i1][i2], &pfin[i1][i2], in_roots, wp_bin);
      // dbg = 0;
      // cout << "pfout:" << endl;
      // poly_frac_print(&pfout[i1][i2]);
      // if (pfout[i1][i2].mults) {
      //   if (pfout[i1][i2].mults[0] != 0) {
      //     *out_zero_label = 0;
      //   }
      // }
      // mpc_t val, eta_val, inv_eta_val;
      // mpc_init3(val, wp2, wp2);
      // mpc_init3(eta_val, wp2, wp2);
      // mpc_init3(inv_eta_val, wp2, wp2);
      // mpc_set_ui(eta_val, 17, MPFR_RNDN);
      // mpc_pow_si(inv_eta_val, eta_val, -1, MPFR_RNDN);
      // cout << "eval pfin in eta = 17:" << endl;
      // poly_frac_eval_value(&val, &pfin[i1][i2], in_roots, &eta_val);
      // print_mpc(&val); cout << endl;
      // cout << "eval pfout in 1/eta = 1/17:" << endl;
      // poly_frac_eval_value(&val, &pfout[i1][i2], out_roots, &inv_eta_val);
      // print_mpc(&val); cout << endl;

      // chain-rule factor
      if (pfout[i1][i2].num_vdeg != -1) {
        poly_frac_neg(&pfout[i1][i2]);
        poly_frac_mul_sym_pow(&pfout[i1][i2], &pfout[i1][i2], -2);
      }
    }
  }

  *out_zero_label = 0;

}


void rel_err_poly_shift(mpc_t cnst, mpc_t *d, int n, int wp_bin)
{
  mpc_t tmp_flop; mpc_init3(tmp_flop, wp2, wp2);
	int k,j;
	for (j=0; j<n; j++) {
	  for (k=n-1; k>=j; k--) {
      // mpc_fma(d[k], cnst, d[k+1], d[k], MPFR_RNDN);
      // mpc_mul(tmp_flop, cnst, d[k+1], MPFR_RNDN);
      exp_rel_err_mpc_mul(tmp_flop, cnst, d[k+1], MPFR_RNDN, wp_bin);
      // cout << "cnst = "; mpc_out_str(stdout, 10, 0, cnst, MPFR_RNDN); cout << endl;
      // cout << "d = "; print_mpc(&d[k+1]); cout << endl;
      // cout << "flop = "; print_mpc(&tmp_flop); cout << endl;

      // cout << "d = "; print_mpc(&d[k]); cout << endl;
      // cout << "flop = "; print_mpc(&tmp_flop); cout << endl;
      exp_rel_err_mpc_add(d[k], d[k], tmp_flop, MPFR_RNDN, wp_bin);
      // mpc_add(d[k], d[k], tmp_flop, MPFR_RNDN);
      // cout << "d = "; print_mpc(&d[k]); cout << endl;
    }
  }
  mpc_clear(tmp_flop);
}


void rel_err_poly_frac_shift(
  struct poly_frac *pfout,
  struct poly_frac *pfin,
  mpc_t *shift,
  int swap_root_lab,
  int wp_bin
) {
  if ((*pfin).num_vdeg == -1) {
    poly_frac_set_zero(pfout);
    return;
  }

  if (mpc_lessthan_tol(*shift)) {
    if (pfout != pfin) {
      poly_frac_set_pf(pfout, pfin);
    } else {
      return;
    }
  }

  // SHIFT COEFFICIENTS
  if ((*pfout).coeffs) {
    mpc_rk1_clear((*pfout).coeffs, (*pfout).num_vdeg+1);
    delete[] (*pfout).coeffs;
  }
  (*pfout).coeffs = new mpc_t[(*pfin).num_vdeg+1];
  init_rk1_mpc((*pfout).coeffs, (*pfin).num_vdeg+1);
  for (int k=0; k<=(*pfin).num_vdeg; k++) {
    mpc_set((*pfout).coeffs[k], (*pfin).coeffs[k], MPFR_RNDN);
  }
  // poly_shift(*shift, (*pfout).coeffs, (*pfin).num_vdeg);
  rel_err_poly_shift(*shift, (*pfout).coeffs, (*pfin).num_vdeg, wp_bin);

  // poly_frac_set_pf(pfout, pfin);
  poly_frac_set_info_pf(pfout, pfin);
  poly_frac_set_mults_pf(pfout, pfin);

  // SWAP MULTIPLICITIES
  (*pfout).mults[0] = (*pfin).mults[swap_root_lab];
  (*pfout).mults[swap_root_lab] = (*pfin).mults[0];

  // SHIFT eta^ldeg
  mpc_t root;
  mpc_init3(root, wp2, wp2);
  mpc_neg(root, *shift, MPFR_RNDN);
  int ldeg = (*pfin).num_deg-(*pfin).num_vdeg;
  if (ldeg > 0) {
    mpc_t *holder = (*pfout).coeffs;
    (*pfout).coeffs = new mpc_t[(*pfout).num_deg+1];
    init_rk1_mpc((*pfout).coeffs, (*pfout).num_deg+1);
    rel_err_poly_mul_root_cp((*pfout).coeffs, holder, (*pfout).num_vdeg, root, wp_bin);
    mpc_rk1_clear(holder, (*pfout).num_vdeg+1);
    delete[] holder;
  }
  for (int k=1; k<ldeg; k++) {
    rel_err_poly_mul_root((*pfout).coeffs, (*pfout).num_vdeg+k, root, wp_bin);
    // rel_err_poly_mul_root_cp((*pfout).coeffs, (*pfout).num_vdeg+k, root, wp_bin);
  }
  (*pfout).num_vdeg = (*pfout).num_deg;

  // cout << "before pruning:" << endl;
  // poly_frac_print(pfout);
  // poly_frac_prune(pfout);
  // cout << "after pruning:" << endl;
  // poly_frac_print(pfout);
  // cout << "after normal:" << endl;
  // poly_frac_normal(pfout);
  // poly_frac_print(pfout);
}


void rel_err_poly_frac_rk2_shift(
  struct poly_frac **pfout,
  struct poly_frac **pfin,
  mpc_t *shift,
  int swap_root_lab,
  int dim1, int dim2,
  int wp_bin
) {
  for (int i1=0; i1<dim1; i1++) {
    for (int i2=0; i2<dim2; i2++) {
      // cout << "i1, i2 = " << i1 << ", " << i2 << endl;
      rel_err_poly_frac_shift(&pfout[i1][i2], &pfin[i1][i2], shift, swap_root_lab, wp_bin);
    }
  }
}


//////
// RELATIVE ERROR FUNTIONS
//////
void rel_err_center_around_pole(
  // OUTPUT
  struct poly_frac **sh_pfmat, mpc_t *sh_roots,
  // INPUT
  struct poly_frac **pfmat, mpc_t *roots, int label,
  int nroots, int dim1, int dim2,
  int wp_bin
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
  rel_err_poly_frac_rk2_shift(
    sh_pfmat, pfmat, &roots[label], label,
    dim1, dim2,
    wp_bin
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
    // mpc_sub(sh_roots[k], roots[k], roots[label], MPFR_RNDN);    
    exp_rel_err_mpc_add(sh_roots[k], roots[k], sh_roots[label], MPFR_RNDN, wp_bin);
  }
  // cout << "shifted roots:" << endl;
  // print_poly(sh_roots, nroots-1);
}


void mpc_prune_re_im(
  // IN-OUT
  mpc_t *in,
  // INPUT
  int wp_bin
) {
  // if input is real or imaginary, do nothing
  if (mpfr_zero_p(mpc_imagref(*in))) {
    return;
  }
  if (mpfr_zero_p(mpc_realref(*in))) {
    return;
  }

  int exp_diff = mpfr_get_exp(mpc_realref(*in)) - mpfr_get_exp(mpc_imagref(*in));
  if (exp_diff > wp_bin) {
    mpfr_set_ui(mpc_imagref(*in), 0, MPFR_RNDN);
  } else if (exp_diff < -wp_bin) {
    mpfr_set_ui(mpc_realref(*in), 0, MPFR_RNDN);
  }
}


void mpc_rk1_prune_re_im(
  // IN-OUT
  mpc_t *vec,
  // INPUT
  int wp_bin, int dim
) {
  for (int k=0; k<dim; k++) {
    mpc_prune_re_im(&vec[k], wp_bin);
  }
}


void rel_err_mpc_rk2_mul_mpc_rk1(
  // OUTPUT
  mpc_t *out_tens,
  // INPUT
  mpc_t **tens1, mpc_t *tens2,
  int dim1, int dim2,
  int wp_bin
) {
  mpc_t tmpc; mpc_init3(tmpc, wp2, wp2);
  for (int i1=0; i1<dim1; i1++) {
    mpc_set_ui(out_tens[i1], 0, MPFR_RNDN);
    for (int i2=0; i2<dim2; i2++) {
      cout << "i1, i2 = " << i1 << ", " << i2 << endl;
      // mpc_fma(
      //   out_tens[i1],
      //   tens1[i1][i2], tens2[i2], out_tens[i1],
      //   MPFR_RNDN
      // );
      // mpc_mul(tmpc, tens1[i1][i2], tens2[i2], MPFR_RNDN);
      exp_rel_err_mpc_mul(tmpc, tens1[i1][i2], tens2[i2], MPFR_RNDN, wp_bin);
      cout << "cum = "; print_mpc(&out_tens[i1]); cout << endl;
      cout << "fac1 = "; print_mpc(&tens1[i1][i2]); cout << endl;
      cout << "fac2 = "; print_mpc(&tens2[i2]); cout << endl;
      cout << "tmpc = "; print_mpc(&tmpc); cout << endl;
      // mpc_add(out_tens[i1], out_tens[i1], tmpc, MPFR_RNDN);
      exp_rel_err_mpc_add(out_tens[i1], out_tens[i1], tmpc, MPFR_RNDN, wp_bin);
    }
    cout << "tot = "; print_mpc(&out_tens[i1]); cout << endl;
  }
  mpc_clear(tmpc);
}

