#include <iostream>
#include <cstring>
// #include <complex>
#include "mpc.h"

using namespace std;

#include "malloc_defs.h"
#include "utils.h"
#include "tensor_utils.h"
#include "setup.h"
#include "topology.h"
#include "codec.h"
#include "poly_frac.h"
#include "system_analyzer.h"
#include "kira.h"
extern "C" {
  #include "ex_tree.h"
  #include "cpoly.h"
  #include "mp_roots_poly.h"
  #include "in_out.h"
}


void mpz_rk1_sort_pows(
  // OUTPUT
  int *deg, mpz_t **coeffs,
  // INPUT
  int nterms, mpz_t *terms, int *pows
) {
  int k;

  // find degree as max power
  *deg = 0;
  for (k=0; k<nterms; k++) {
    if (pows[k] > *deg) {
      *deg = pows[k];
    }
  }

  // allocate
  *coeffs = new mpz_t[*deg+1];
  for (k=0; k<=*deg; k++) {
    mpz_init((*coeffs)[k]);
    mpz_set_ui((*coeffs)[k], 0);
  }

  // sort and sum
  for (k=0; k<nterms; k++) {
    mpz_add((*coeffs)[pows[k]], (*coeffs)[pows[k]], terms[k]);
  }

}


void mpq_rk1_sort_pows(
  // OUTPUT
  int *deg, mpq_t **coeffs,
  // INPUT
  int nterms, mpq_t *terms, int *pows
) {
  int k;

  // find degree as max power
  *deg = 0;
  for (k=0; k<nterms; k++) {
    if (pows[k] > *deg) {
      *deg = pows[k];
    }
  }

  // allocate
  *coeffs = new mpq_t[*deg+1];
  for (k=0; k<=*deg; k++) {
    mpq_init((*coeffs)[k]);
    mpq_set_ui((*coeffs)[k], 0, 1);
  }

  // sort and sum
  for (k=0; k<nterms; k++) {
    mpq_add((*coeffs)[pows[k]], (*coeffs)[pows[k]], terms[k]);
  }

}


void mpz_poly_eval_mpq(
  // OUTPUT
  mpq_t *out,
  // INPUT
  int deg, mpz_t *coeffs, mpz_t *npows, mpz_t *dpows
) {
  // cout << "deg = " << deg << endl;
  mpz_t tmp; mpz_init(tmp);
  // NUMERATOR
  mpz_set_ui(mpq_numref(*out), 0);
  for (int k=0; k<=deg; k++) {
    // cout << "k = " << k << endl;
    if (mpz_sgn(coeffs[k]) == 0) {
      continue;
    }
    // cout << "npows = "; mpz_out_str(stdout, 10, npows[k]); cout << endl;
    // cout << "dpows = "; mpz_out_str(stdout, 10, dpows[deg-k]); cout << endl;
    // cout << "coeffs = "; mpz_out_str(stdout, 10, coeffs[k]); cout << endl;
    mpz_mul(tmp, npows[k], dpows[deg-k]);
    mpz_mul(tmp, tmp, coeffs[k]);
    mpz_add(mpq_numref(*out), mpq_numref(*out), tmp);
  }

  // DENOMINATOR
  mpz_set(mpq_denref(*out), dpows[deg]);

  // SIMPLIFY
  mpq_canonicalize(*out);
  
  // FREE
  mpz_clear(tmp);
}


typedef struct BitField {
  unsigned int bit : 1;
} BitField;


char* read_file_into_string(const char* filename) {
	FILE* file = fopen(filename, "r");
	if (file == NULL) {
		perror("Errore apertura file");
		return NULL;
	}

	// go to the end of file to compute size
	fseek(file, 0, SEEK_END);
	long file_size = ftell(file);
	rewind(file);

	// alloc memory for the string
	char* buffer = (char*)malloc((file_size + 1) * sizeof(char));
	if (buffer == NULL) {
		perror("Errore allocazione memoria");
		fclose(file);
		return NULL;
	}

	// read file and close
	size_t read_size = fread(buffer, sizeof(char), file_size, file);
	buffer[read_size] = '\0';  // assicura il terminatore nullo

	fclose(file);
	return buffer;
}


typedef struct pfac {
  int deg;  // polynomial degree
  mpq_t *coeffs;  // exact integer coefficients
  int *labs;  // root labels
} pfac;


void pfac_build(pfac *fac) {
  // fac->nroots = 0;
  fac->labs = NULL;
  fac->deg = 0;
  fac->coeffs = NULL;
}


void pfac_free(
  pfac *fac
) {
  if (fac->coeffs) {
    for (int k=0; k<=fac->deg; k++) {
      mpq_clear(fac->coeffs[k]);
    }
    delete[] fac->coeffs;
  }
  if (fac->labs) {
    delete[] fac->labs;
  }
}


void pfac_print(pfac *fac) {
  cout << "deg = " << fac->deg << endl;
  cout << "coeffs:" << endl;
  for (int k=0; k<=fac->deg; k++) {
    cout << "pow " << k << ": ";
    mpq_out_str(stdout, 10, fac->coeffs[k]); cout << endl;
  }
  cout << "labels: {";
  if (fac->labs) {
    for (int k=0; k<fac->deg; k++) {
      cout << "r" << fac->labs[k] << ", ";
    }
  cout << "}" << endl;
  }
}


typedef struct pfaclist {
  int nfacs;  // number of items, -1 if not first item
  pfac *fac;  // pfac item
  pfaclist *next;  // next pfac item, NULL for last item
} pfaclist;


void pfaclist_build(pfaclist *pfacl) {
  pfacl->nfacs = 0;
  pfacl->fac = NULL;
  pfacl->next = NULL;
}


void pfaclist_free(pfaclist *pfacl) {
  if (pfacl->fac) {
    pfac_free(pfacl->fac);
    free(pfacl->fac);
  }
  if (pfacl->next) {
    pfaclist_free(pfacl->next);
    free(pfacl->next);
  }
}


typedef struct pden {
  int deg;  // degree of denominator
  int ldeg;  // lower degree of denominator
  int nroots;  // global number of roots
  int *mults;  // multiplicity of each root
  int nfacs;  // number of factors
  int *facs;  // array of factor labels
  int *fmults;  // multiplicity of each factor
  mpq_t norm;  // normalization constant
} pden;


void pden_build(
  pden *den
) {
  den->deg = 0;
  den->ldeg = 0;
  den->nroots = 0;
  den->mults = NULL;
  den->nfacs = 0;
  den->facs = NULL;
  den->fmults = NULL;
  mpq_init(den->norm);
  mpq_set_ui(den->norm, 1, 1);
}


void pden_free(
  pden *den
) {
  if (den->mults) {
    delete[] den->mults;
  }
  if (den->facs) {
    delete[] den->facs;
  }
  if (den->fmults) {
    delete[] den->fmults;
  }
  mpq_clear(den->norm);
}


void pden_print(pden *den) {
  cout << "deg, ldeg = " << den->deg << ", " << den->ldeg << endl;
  cout << "roots: {";
  if (den->mults) {
    for (int k=0; k<den->nroots; k++) {
      if (den->mults[k] > 0) {
        cout << "r" << k << ": " << den->mults[k] << ", ";
      }
    }
  }
  cout << "}" << endl;
  cout << "factors: {";
  if (den->facs) {
    for (int f=0; f<den->nfacs; f++) {
      cout << "f" << den->facs[f] << ": " << den->fmults[f] << ", ";
    }
  }
  cout << "}" << endl;
  cout << "norm: "; mpq_out_str(stdout, 10, den->norm); cout << endl;
}


void pden_set_ui(
  pden *den, int ui
) {
  den->deg = 0;
  den->ldeg = 0;
  den->nfacs = 0;
  mpq_set_ui(den->norm, ui, 1);
}


void pden_set_pden(
  pden *out, pden *in
) {
  // check if input is equal to output
  if (out == in) {
    return;
  }

  // check if input is uninitialized
  if (in->deg == -1) {
    out->deg = -1;
    return;
  }
  
  // copy degree
  out->deg = in->deg;
  out->ldeg = in->ldeg;

  // copy root mults
  if (out->nroots != in->nroots) {
    out->nroots = in->nroots;
    if (out->mults) {
      delete[] out->mults;
    }
    if (in->mults) {
      out->mults = new int[in->nroots];
    }
  }
  for (int k=0; k<in->nroots; k++) {
    out->mults[k] = in->mults[k];
  }

  // copy factor labels and mults
  if (out->nfacs != in->nfacs) {
    out->nfacs = in->nfacs;
    if (out->facs) {
      delete[] out->facs;
    }
    if (in->facs) {
      out->facs = new int[in->nfacs];
    }
    if (out->fmults) {
      delete[] out->fmults;
    }
    if (in->fmults) {
      out->fmults = new int[in->nfacs];
    }
  }
  for (int k=0; k<in->nfacs; k++) {
    out->facs[k] = in->facs[k];
    out->fmults[k] = in->fmults[k];
  }

  // copy normalization
  mpq_set(out->norm, in->norm);

}


void pfac_set_coeffs(
  // OUTPUT
  pfac *fac,
  // INPUT
  mpq_t *coeffs, int deg
) {
  fac->deg = deg;

  // if (fac->subs) {
  //   delete[] fac->subs;
  // }
  // fac->subs = new int*[deg+1];

  // if (fac->sdim) {
  //   delete[] fac->sdim;
  // }
  // fac->sdim = new int[deg+1];

  
  // if (fac->mults) {
  //   delete[] fac->mults;
  // }
  // fac->mults = new int[deg+1];

  if (fac->coeffs) {
    delete[] fac->coeffs;
  }
  fac->coeffs = new mpq_t[deg+1];

  for (int k=0; k<=deg; k++) {
    mpq_init(fac->coeffs[k]);
    mpq_set(fac->coeffs[k], coeffs[k]);
  }

}


void pfac_set_mults(
  // OUTPUT
  pfac *fac,
  // INPUT
  int *mults, int nroots
) {
  int count = 0;
  for (int k=1; k<nroots; k++) {
    for (int m=0; m<mults[k]; m++) {
      fac->labs[count++] = k;
    }
  }
}


int coeffs_div_pfac(
  // OUTPUT
  // int *mults_rem,
  // IN-OUT
  mpq_t *coeffs,
  // INPUT
  int deg, int *mults, int nroots,
  pfac *fac
) {
  int fully_in = 0;
  int *mults_cp = new int[nroots];
  for (int k=0; k<nroots; k++) {
    mults_cp[k] = mults[k];
  }

  int deg_in = 0;
  for (int k=0; k<fac->deg; k++) {
    // cout << "k, lab = " << k << ", " << fac->labs[k] << endl;
    if (mults_cp[fac->labs[k]] > 0) {
      mults_cp[fac->labs[k]]--;
      deg_in++;
    }
  }

  if (deg_in > 0) {
    if (deg_in == fac->deg) {
      fully_in = 1;
    } else { 
      fully_in = -1;
    }
  }

  if (dbg) cout << "fully_in = " << fully_in << endl;
  if (fully_in == 1) {
    // // lower multiplicites
    // for (int k=0; k<fac->deg; k++) {
    //   mults[fac->labs[k]]--;
    // }

    // perform exact division
    if (dbg) cout << "perform exact division" << endl;
    if (div_poly_mpq(coeffs, coeffs, fac->coeffs, deg, fac->deg)) {
      perror("exact polynomial division gave remainder");
      exit(1);
    }
    if (dbg) {
    cout << "result:" << endl;
    for (int k=0; k<=deg-fac->deg; k++) {
      cout << "pow " << k << ": ";
      mpq_out_str(stdout, 10, coeffs[k]); cout << endl;
    }
    }
  }

  // FREE
  delete[] mults_cp;

  return fully_in*deg_in;
}


void pfaclist_insert_empty(
  // IN-OUT
  pfaclist *facl
) {
  pfaclist *holder = facl->next;
	facl->next = (pfaclist*) malloc(sizeof(struct pfaclist));
  pfaclist_build(facl->next);
	facl->next->next = holder;
}


void pfaclist_remove(
  // IN-OUT
  pfaclist *facs, pfac *fac
) {
  pfaclist *pfacl = facs;
  while (pfacl) {
    if (!pfacl->next) {
      break;
    }
    if (pfacl->next->fac == fac) {
      pfac_free(pfacl->next->fac);
      pfacl->next = pfacl->next->next;
      facs->nfacs--;
    }
    pfacl = pfacl->next;
  }

}


void pfac_check(
  pfac *fac, mpc_t *roots
) {
  // reconstruct coeffs from roots
  mpc_t *rec_cfs = new mpc_t[fac->deg+1];
  init_rk1_mpc(rec_cfs, fac->deg+1);
  mpc_set_ui(rec_cfs[0], 1, MPFR_RNDN);
  for (int k=0; k<fac->deg; k++) {
    poly_mul_root(rec_cfs, k, roots[fac->labs[k]]);
  }

  // convert coeffs from rational
  mpc_t *conv_cfs = new mpc_t[fac->deg+1];
  init_rk1_mpc(conv_cfs, fac->deg+1);
  for (int k=0; k<=fac->deg; k++) {
    mpfr_set_q(mpc_realref(conv_cfs[k]), fac->coeffs[k], MPFR_RNDN);
    mpfr_set_ui(mpc_imagref(conv_cfs[k]), 0, MPFR_RNDN);
  }

  // check
  mpc_rk1_compare_double(conv_cfs, rec_cfs, fac->deg+1);

  delete[] rec_cfs;
  delete[] conv_cfs;
}


int pfac_cmp(
  // OUTPUT
  int **is_in,
  // INPUT
  pfac *fac1, pfac *fac2
) {
  if (*is_in) {
    delete[] *is_in;
  }
  *is_in = new int[fac2->deg];

  int found, k1, k2;
  for (k2=0; k2<fac2->deg; k2++) {
    (*is_in)[k2] = 0;
  }
  for (k1=0; k1<fac1->deg; k1++) {
    found = 0;
    for (k2=0; k2<fac2->deg; k2++) {
      if ((*is_in)[k2] == 1) {continue;}
      if (fac2->labs[k2] == fac1->labs[k1]) {
        found = 1;
        (*is_in)[k2] = 1;
        break;
      }
    }
    if (found) {
      continue;
    } else {
      return 0;
    }
  }

  return 1;
}


// <<pfac_split>>
int pfac_split(
  // IN-OUT
  pfaclist *pfacl, pfaclist *facs,
  // INPUT
  int sub_deg, mpq_t *sub_coeffs,
  int *mults, int nroots
) {
  // insert empty node for the new factor
  // cout << "insert empty node for the new factor" << endl;
  pfaclist_insert_empty(pfacl);
  pfacl->next->fac = (pfac*) malloc(sizeof(pfac));
  pfac_build(pfacl->next->fac);

  // set coeffs of polynomial division
  // cout << "set coeffs of polynomial division" << endl;
  div_poly_mpq(
    pfacl->fac->coeffs, pfacl->fac->coeffs, sub_coeffs,
    pfacl->fac->deg, sub_deg
  );

  // set coeffs of new factor
  // cout << "set coeffs of new factor" << endl;
  pfacl->next->fac->coeffs = new mpq_t[sub_deg+1];
  for (int k=0; k<=sub_deg; k++) {
    mpq_init(pfacl->next->fac->coeffs[k]);
    mpq_set(pfacl->next->fac->coeffs[k], sub_coeffs[k]);
  }

  // copy old root labels
  // cout << "copy old root labels" << endl;
  int deg_old = pfacl->fac->deg;
  int *labs_old = new int[pfacl->fac->deg];
  copy_rk1_int(labs_old, pfacl->fac->labs, pfacl->fac->deg);

  // set degrees
  pfacl->fac->deg -= sub_deg;
  pfacl->next->fac->deg = sub_deg;

  // sort root labels between the two sub-factors
  // cout << "sort root labels between the two sub-factors" << endl;
  delete[] pfacl->fac->labs;
  pfacl->fac->labs = new int[pfacl->fac->deg];
  pfacl->next->fac->labs = new int[sub_deg];

  int count1 = 0, count2 = 0;
  for (int idx, k=0; k<deg_old; k++) {
    idx = labs_old[k];
    // cout << "k, idx, mult = " << k << ", " << idx << ", " << mults[idx] << endl;
    if (mults[idx] > 0) {
      pfacl->next->fac->labs[count1++] = idx;
      mults[idx]--;
    } else {
      pfacl->fac->labs[count2++] = idx;
    }
  }

  // // check for merge
  // pfaclist *pfaclp = facs;
  // for (int fp=0; fp<facs->nfacs; fp++) {
  //   if (pfac_cmp(pfaclp->fac, pfacl->fac)) {
  //     return fp;
  //   }
  //   if (pfaclp->next) {
  //     pfaclp = pfaclp->next;
  //   }
  // }

  return -1;

}


int pfac_cut(
  // IN-OUT
  pfaclist *pfacl, pfaclist *facs,
  // INPUT
  int sub_deg, mpq_t *sub_coeffs,
  int *mults, int nroots
) {
  // set coeffs of polynomial division
  // cout << "set coeffs of polynomial division" << endl;
  div_poly_mpq(
    pfacl->fac->coeffs, pfacl->fac->coeffs, sub_coeffs,
    pfacl->fac->deg, sub_deg
  );

  // copy old root labels
  // cout << "copy old root labels" << endl;
  int deg_old = pfacl->fac->deg;
  int *labs_old = new int[pfacl->fac->deg];
  copy_rk1_int(labs_old, pfacl->fac->labs, pfacl->fac->deg);

  // set degree
  pfacl->fac->deg -= sub_deg;

  delete[] pfacl->fac->labs;
  pfacl->fac->labs = new int[pfacl->fac->deg];

  int count = 0;
  for (int idx, k=0; k<deg_old; k++) {
    idx = labs_old[k];
    // cout << "k, idx, mult = " << k << ", " << idx << ", " << mults[idx] << endl;
    if (mults[idx] > 0) {
      mults[idx]--;
    } else {
      pfacl->fac->labs[count++] = idx;
    }
  }

  // // check for merge
  // pfaclist *pfaclp = facs;
  // for (int fp=0; fp<facs->nfacs; fp++) {
  //   if (pfac_cmp(pfaclp->fac, pfacl->fac)) {
  //     return fp;
  //   }
  //   if (pfaclp->next) {
  //     pfaclp = pfaclp->next;
  //   }
  // }

  return -1;

}


void pfac_split_pfac(
  // IN-OUT
  pfaclist *pfacl,
  // INPUT
  int sub_deg, mpq_t *sub_coeffs,
  int *is_in
) {
  // // inserty empty node for the new factor
  // // cout << "insert empty node for the new factor" << endl;
  // pfaclist_insert_empty(pfacl);
  // pfacl->next->fac = (pfac*) malloc(sizeof(pfac));
  // pfac_build(pfacl->next->fac);

  // set coeffs of polynomial division
  // cout << "set coeffs of polynomial division" << endl;
  div_poly_mpq(
    pfacl->fac->coeffs, pfacl->fac->coeffs, sub_coeffs,
    pfacl->fac->deg, sub_deg
  );

  // // set coeffs of new factor
  // // cout << "set coeffs of new factor" << endl;
  // pfacl->next->fac->coeffs = new mpq_t[sub_deg+1];
  // for (int k=0; k<=sub_deg; k++) {
  //   mpq_init(pfacl->next->fac->coeffs[k]);
  //   mpq_set(pfacl->next->fac->coeffs[k], sub_coeffs[k]);
  // }

  // copy old root labels
  // cout << "copy old root labels" << endl;
  int deg_old = pfacl->fac->deg;
  int *labs_old = new int[pfacl->fac->deg];
  copy_rk1_int(labs_old, pfacl->fac->labs, pfacl->fac->deg);

  // set degrees
  pfacl->fac->deg -= sub_deg;
  // pfacl->next->fac->deg = sub_deg;
  // make factor monic
  for (int k=0; k<pfacl->fac->deg; k++) {
    mpq_div(pfacl->fac->coeffs[k], pfacl->fac->coeffs[k], pfacl->fac->coeffs[pfacl->fac->deg]);
  }
  mpq_set_ui(pfacl->fac->coeffs[pfacl->fac->deg], 1, 1);


  // sort root labels between the two sub-factors
  // cout << "sort root labels between the two sub-factors" << endl;
  delete[] pfacl->fac->labs;
  pfacl->fac->labs = new int[pfacl->fac->deg];
  // pfacl->next->fac->labs = new int[sub_deg];

  // int count1 = 0, count2 = 0;
  // for (int k=0; k<deg_old; k++) {
  //   // cout << "k, idx, mult = " << k << ", " << idx << ", " << mults[idx] << endl;
  //   if (is_in[k] > 0) {
  //     pfacl->next->fac->labs[count1++] = labs_old[k];
  //   } else {
  //     pfacl->fac->labs[count2++] = labs_old[k];
  //   }
  // }

  int count = 0;
  for (int k=0; k<deg_old; k++) {
    // cout << "k, idx, mult = " << k << ", " << idx << ", " << mults[idx] << endl;
    if (is_in[k] == 0) {
      pfacl->fac->labs[count++] = labs_old[k];
    }
  }
}


int pfaclist_split(
  // IN-OUT
  pfaclist *facs,
  // INPUT
  pfac *fac, mpc_t *roots
) {
  pfaclist *pfacl = facs;
  int *is_in = NULL, remove_fac = 0;
  while (pfacl) {
    if (pfacl->fac == fac) {
      pfacl = pfacl->next;
      continue;
    }
    
    // check whether fac is contained
    if (pfac_cmp(&is_in, fac, pfacl->fac)) {
      // fac is contained
      if (fac->deg == pfacl->fac->deg) {
        // fac is exactly the same
        // mark factor as to be removed
        remove_fac = 1;
      } else {
        // split factor
        // cout << "factor before split:" << endl;
        // pfac_print(pfacl->fac);
        pfac_split_pfac(pfacl, fac->deg, fac->coeffs, is_in);
        // cout << "factors after split:" << endl;
        // pfac_print(pfacl->fac);
        // pfac_check(pfacl->fac, roots);
        // pfac_print(fac);
        // pfac_check(fac, roots);
        
        // recursively split new factor list
        pfaclist_split(facs, pfacl->fac, roots);
      }
    }

    pfacl = pfacl->next;
  }

  if (remove_fac) {
    // cout << "remove factor:" << endl;
    // pfac_print(fac);
    pfaclist_remove(facs, fac);
  }

  return remove_fac;
}


void pden_roots_append_zero(
  // IN-OUT
  pden *den,
  // INPUT
  int nroots
) {
  if (den->nroots < nroots) {
    if (den->mults) {
      int_rk1_append_zero(&den->nroots, &den->mults, nroots - den->nroots);
    } else {
      // cout << "alloc" << endl;
      den->mults = new int[nroots];
      for (int k=0; k<nroots; k++) {
        den->mults[k] = 0;
      }

      den->nroots = nroots;
    }
  }
}


void root_prof_update_mults(
  // IN-OUT
  int *mults,
  // INPUT
  int *root_prof, int deg, int nroots
) {
  for (int k=0; k<deg; k++) {
    mults[root_prof[k]]++;
  }
}


void root_prof_update_fac(
  // IN-OUT
  pfac *fac, int idx,
  // INPUT
  int *root_prof, int deg, int nroots
) {  
  for (int k=0; k<deg; k++) {
    fac->labs[idx+k] = root_prof[k];
  }
}


// <<poly_mpq_test_roots>>
void poly_mpq_test_roots(
  // OUTPUT
  int *mults2, int *mults, int *npass,
  // INPUT
  mpq_t *coeffs, int deg,
  mpc_t *roots, int nroots,
  int *root_labs
) {
  int k;

  // needed to perform the actual root test
	mpfr_t EPSS;
	mpfr_init2(EPSS, wp2);
	mpfr_mul_ui(EPSS, mpfr_tol, 10, MPFR_RNDN);
	mpfr_t EPSS_rel_increase;
	mpfr_init2(EPSS_rel_increase, wp2);
	mpfr_set_ui(EPSS_rel_increase, 10, MPFR_RNDN);
	mpfr_pow_si(EPSS_rel_increase, EPSS_rel_increase, wp2*3/10/20, MPFR_RNDN);
  mpfr_mul(EPSS, EPSS, EPSS_rel_increase, MPFR_RNDN);
	mpfr_t err_tmp, tmpfr_tmp, abx_tmp;
	mpfr_init2(err_tmp, wp2);
	mpfr_init2(tmpfr_tmp, wp2);
	mpfr_init2(abx_tmp, wp2);
	mpc_t b_tmp, d_tmp, f_tmp;
	mpc_init3(b_tmp, wp2, wp2);
	mpc_init3(d_tmp, wp2, wp2);
	mpc_init3(f_tmp, wp2, wp2);

  // convert rational coefficients to wp
  mpc_t *coeffs_mpc = new mpc_t[deg+1];
  // *npass = 0;
  for (k=0; k<=deg; k++) {
    mpc_init3(coeffs_mpc[k], wp2, wp2);
    mpc_set_q(coeffs_mpc[k], coeffs[k], MPFR_RNDN);
  }

  // test roots
	mpc_t *defl = new mpc_t[deg+1];
	for (k=0; k<=deg; k++) {
		mpc_init3(defl[k], wp2, wp2);
	}

  int mul;
  for (int r, kp=0; kp<nroots; kp++) {
    r = root_labs[kp];
    if (dbg) {
    cout << "r = " << r << ", ";
    cout << endl << "eval in "; print_mpc(&roots[r]); cout << endl;
    }
    if (mults[r] > 0) {
      if (dbg) cout << "mul = " << mults[r] << " already" << endl;
      continue;
    }

    // restore deflation to original polynomial
    for (k=0; k<=deg; k++) {
      mpc_set(defl[k], coeffs_mpc[k], MPFR_RNDN);
    }

    // find multiplicity
    mul = 0;
    while (1) {
      // cout << "cum mul = " << mul << endl;
      if (mul == deg) {break;}

      if (mul > 0) {
        // perform deflation
        deflation(defl, deg-mul+1, &roots[r], wp2);
      }

      // test root
      if (!test_root(
        defl, deg-mul, &roots[r], wp2, EPSS, &err_tmp,
        &b_tmp, &d_tmp, &f_tmp, &abx_tmp, &tmpfr_tmp
      )) {
        break;
      }
      mul++;
    }
    if (dbg) cout << "mul = " << mul << endl;
    mults[r] += mul;
    if (mults2) {
      mults2[r] +=mul;
    }
    *npass += mul;
  }

  // FREE
  mpc_rk1_clear(coeffs_mpc, deg+1); delete[] coeffs_mpc;
  mpc_rk1_clear(defl, deg+1); delete[] defl;
  mpfr_clear(EPSS);
  mpfr_clear(EPSS_rel_increase);
  mpfr_clear(err_tmp); mpfr_clear(tmpfr_tmp); mpfr_clear(abx_tmp);
  mpc_clear(b_tmp); mpc_clear(d_tmp); mpc_clear(f_tmp);
  
}


void pden_facs_update(
  // IN-OUT
  pden *den,
  // INPUT
  int *split_idx, int num_split
  // int *split_idx, int *merge_idx, int num_split
) {
  // // find index of split factor
  // int f_split = -1;
  // for (int f=0; f<den->nfacs; f++) {
  //   if (den->facs[f] == split_idx) {
  //     f_split = f;
  //     break;
  //   }
  // }
  // if (f_split == -1) {
  //   // split factor is not present
  //   return;
  // }

  // // shift one place up every factor after the split
  // for (int f=den->nfacs; f>f_split; f--) {
  //   den->facs[f+1] = den->facs[f]+1;
  //   den->fmults[f+1] = den->fmults[f];
  // }
  // den->facs[f_split+1] = split_idx+1;
  // den->fmults[f_split+1] = den->fmults[f_split];
  // den->nfacs++;

  if (split_idx[0] < 0) {
    return;
  }

  int count_split = 0;
  for (int c=0; c<num_split; c++) {
    for (int f=0; f<den->nfacs; f++) {
      if (den->facs[f] == split_idx[c]) {
        count_split += den->fmults[f];
        break;
      }
    }
  }

  int f_split = -1;  // index of split factor
  for (int f=0; f<den->nfacs; f++) {
    if (den->facs[f] > split_idx[0]) {
      den->facs[f]++;
    } else if (den->facs[f] == split_idx[0]) {
      f_split = f;
    }
  }

  if (count_split > 0) {
    den->facs[den->nfacs] = split_idx[0];
    den->fmults[den->nfacs] = count_split;
    den->nfacs++;
  }

}


void pden_group_roots(
  // IN-OUT
  pden *den,
  // INPUT
  pfaclist *facs
) {
  int f, k, kp, fac_is_in, lab, high_mul;
  // clear denominator
  for (f=0; f<den->nfacs; f++) {
    den->fmults[f] = 0;
  }
  den->nfacs = 0;

  
  int *mults_cp = new int[den->nroots];
  copy_rk1_int(mults_cp, den->mults, den->nroots);
  int *mults_neg = new int[den->nroots];
  for (k=0; k<den->nroots; k++) {
    mults_neg[k] = 0;
  }

  pfaclist *pfacl = facs;
  for (f=0; f<facs->nfacs; f++) {
    fac_is_in = 1;
    high_mul = 0;
    while(fac_is_in) {
      // check multiplicities
      for (k=0; k<pfacl->fac->deg; k++) {
        lab = pfacl->fac->labs[k];
        if (mults_cp[lab] - mults_neg[lab] > 0) {
          mults_neg[lab]++;
        } else {
          fac_is_in = 0;
          break;
        }
      }
      
      if (fac_is_in) {
        // lower multiplicites
        for (k=0; k<pfacl->fac->deg; k++) {
          lab = pfacl->fac->labs[k];
          mults_cp[lab]--;
        }

        // increase denominator factor multiplicity
        if (high_mul == 0) {
          den->facs[den->nfacs] = f;
          high_mul = 1;
        }
        den->fmults[den->nfacs]++;
      } else {
        if (high_mul) {
          den->nfacs++;
        }
      }

      for (k=0; k<pfacl->fac->deg; k++) {
        lab = pfacl->fac->labs[k];
        mults_neg[lab] = 0;
      }

    }

    if (pfacl->next) {
      pfacl = pfacl->next;
    };
  }

  // FREE
  delete[] mults_cp;
  delete[] mults_neg;
  
}


int denom_coeffs_to_facs(
  // OUTPUT
  pden *den, int *num_split, int **split_idx,
  // IN-OUT
  pfaclist *facs,
  int *nroots, mpc_t **roots, mpfr_t **tols,
  // INPUT
  mpq_t *coeffs, int deg
) {
  // DEGREE AND LOWER DEGREE
  den->deg = deg;
  den->ldeg = 0;
	for (int k=0; k<=deg; k++) {
		// cout << "sign of "; mpq_out_str(stdout, 10, coeffs_mpq[k]); cout << endl;
		if (mpq_sgn(coeffs[k]) == 0) {
			den->ldeg++;
		}	else {
			break;
		}
	}
  // cout << "ldeg = " << den->ldeg << endl;

  den->nroots = *nroots;
  if (den->mults) {
    delete[] den->mults;
  }
  den->mults = new int[*nroots];
  den->mults[0] = den->ldeg;
  int *mults_rem = new int[*nroots];
  for (int k=1; k<*nroots; k++) {
    den->mults[k] = 0;
    mults_rem[k] = 0;
  }

  deg -= den->ldeg;
  coeffs += den->ldeg;

  //////
  // FIND PREVIOUS ROOTS
  //////
  // evaluate polynomial in previous roots
  // cout << "evaluate polynomial in previous roots" << endl;
  int nprev_roots = 0;
  // poly_mpq_test_roots(
  //   den->mults+1, &nprev_roots, coeffs, deg,
  //   *roots+1, *nroots-1,
  // );
  // if (dbg) cout << "num prev roots = " << nprev_roots << endl;

  pfaclist *pfacl = facs;
  // // update global nroots
  // pfaclist *pfacl = facs;
  // if (nprev_roots < deg) {
  //   // cout << "update global nroots" << endl;
  //   // update nroots of denominator
  //   pden_roots_append_zero(den, *nroots + deg - nprev_roots);
  // }

  //////
  // FIND FACTORS
  //////
  // divide previous factors
  // cout << "divide previous factors" << endl;
  den->nfacs = 0;
  if (den->facs) {
    delete[] den->facs;
  }
  den->facs = new int[deg+1];
  if (den->fmults) {
    delete[] den->fmults;
  }
  den->fmults = new int[deg+1];
  for (int k=0; k<=deg; k++) {
    den->fmults[k] = 0;
  }
  int nrem_roots = deg;
  pfacl = facs;
  int *fac_incl = new int[facs->nfacs];
  for (int is_fac, f=0; f<facs->nfacs; f++) {
    if (dbg) cout << "f = " << f << endl;
    if (dbg) pfac_print(pfacl->fac);

    poly_mpq_test_roots(
      mults_rem, den->mults, &nprev_roots, coeffs, deg,
      *roots, pfacl->fac->deg, pfacl->fac->labs
    );

    is_fac = 0;
    while ((fac_incl[f] = coeffs_div_pfac(
      coeffs, deg, mults_rem, *nroots,
      pfacl->fac)) > 0) {
      if (dbg) cout << "factor is divided" << endl;
      // lower multiplicites
      for (int k=0; k<pfacl->fac->deg; k++) {
        mults_rem[pfacl->fac->labs[k]]--;
      }

      // append factor labael
      den->facs[den->nfacs] = f;
      den->fmults[den->nfacs]++;
      nrem_roots -= pfacl->fac->deg;
      nprev_roots -= pfacl->fac->deg;
      deg -= pfacl->fac->deg;  
      is_fac = 1;
    }
    if (is_fac) {
      den->nfacs++;
    }    
    if (pfacl->next) {
      pfacl = pfacl->next;
    }
  }
  if (dbg) cout << "num rem roots = " << nrem_roots << endl;
  if (dbg) cout << "num prev roots = " << nprev_roots << endl;

  // update global nroots
  if (nprev_roots < deg) {
    // cout << "update global nroots" << endl;
    // update nroots of denominator
    pden_roots_append_zero(den, *nroots + deg - nprev_roots);
  }

  // // partially included factors
  // *num_split = 0;
  // pfacl = facs;
  // for (int f=0; f<facs->nfacs; f++) {
  //   if (fac_incl[f] < 0 && fac_incl[f] == -nrem_roots) {
  //     (*num_split)++;
  //   }
  //   if (pfacl->next) {
  //     pfacl = pfacl->next;
  //   }
  // }
  // if (*num_split > 0) {
  //   (*split_idx) = new int[*num_split];
  //   // (*merge_idx) = new int[*num_split];
  // } else {
  //   (*split_idx) = new int[1];
  //   // (*merge_idx) = new int[1];
  // }

  // int count_split = 0;
  // (*split_idx)[0] = -1;
  // pfacl = facs;
  // for (int f=0; f<facs->nfacs; f++) {
  //   if (dbg) cout << "f, incl = " << f << ", " << fac_incl[f] << endl;
  //   if (fac_incl[f] < 0 && fac_incl[f] == -nrem_roots) {
  //     // split factor      
  //     if (count_split == 0) {
  //       if (dbg) cout << "factor to be split" << endl;
  //       pfac_print(pfacl->fac);
  //       cout << "" << endl;
  //     } else {
  //       if (dbg) cout << "cut factor" << endl;
  //       // suppose denominator remainder is {r1} and there are two partially included
  //       // factors f=f1: {r1,r2}, f=f2: {r1,r3}, then the first one is split in
  //       // {r1},{r2}, while the second must be cut in {r3}
  //       cout << "factor before cut (f = " << f << "):" << endl;
  //       pfac_print(pfacl->fac);
  //       pfac_cut(pfacl, facs, nrem_roots, coeffs, mults_rem, *nroots);
  //       // (*merge_idx)[count_split] = pfac_cut(pfacl, facs, nrem_roots, coeffs, mults_rem, *nroots);
  //       cout << "factor after cut:" << endl;
  //       pfac_print(pfacl->fac);
  //       cout << "" << endl;
  //     }
  //     (*split_idx)[count_split++] = f;      
  //   }
  //   if (pfacl->next) {
  //     pfacl = pfacl->next;
  //   }
  // }

  // // split factor
  // if (*num_split > 0) {
  //   if (dbg) cout << "split factor" << endl;
  //   pfacl = facs;
  //   for (int f=0; f<facs->nfacs; f++) {
  //     if (f == (*split_idx)[0]) {
  //       cout << "factor before split (f = " << f << "):" << endl;
  //       pfac_print(pfacl->fac);
  //       pfac_split(pfacl, facs, nrem_roots, coeffs, mults_rem, *nroots);
  //       // (*merge_idx)[0] = pfac_split(pfacl, facs, nrem_roots, coeffs, mults_rem, *nroots);
  //       facs->nfacs++;
  //       // assign new factor to denominator
  //       den->facs[den->nfacs] = f+1;
  //       den->fmults[den->nfacs]++;
  //       cout << "factors after split:" << endl;
  //       pfac_print(pfacl->fac);
  //       pfac_print(pfacl->next->fac);
  //       // update factor labels of denominator
  //       for (int fp=0; fp<den->nfacs; fp++) {
  //         if (den->facs[fp] > f) {
  //           if (dbg) cout << "fp, facs = " << fp << ", " << den->facs[fp] << endl;
  //           den->facs[fp]++;
  //         }
  //       }
  //       den->nfacs++;
  //       return;
  //     }
  //     if (pfacl->next) {
  //       pfacl = pfacl->next;
  //     }
  //   }
  // }  
  // if (dbg) cout << "num rem roots = " << nrem_roots << endl;

  //////
  // SET NORMALIZATION
  //////
  mpq_set(den->norm, coeffs[nrem_roots]);
  if (dbg) {cout << "norm = "; mpq_out_str(stdout, 10, den->norm); cout << endl;}
  if (nrem_roots == 0) {
    // store normalization constant
    // mpq_set(den->norm, coeffs[0]);
    delete[] mults_rem;
    delete[] fac_incl;
    return 0;
  }
  for (int k=0; k<nrem_roots; k++) {
    mpq_div(coeffs[k], coeffs[k], coeffs[nrem_roots]);
  }
  mpq_set_ui(coeffs[nrem_roots], 1, 1);

  //////
  // UPDATE GLOBAL LIST OF FACTORS WITH RESIDUAL COEFFS
  //////
  // cout << "UPDATE GLOBAL LIST OF FACTORS WITH RESIDUAL COEFFS" << endl;
  if (pfacl->fac) {
    pfacl->next = (pfaclist*) malloc(sizeof(pfaclist));
    pfacl = pfacl->next;
    pfaclist_build(pfacl);
  }
  pfacl->fac = (pfac*) malloc(sizeof(pfac));
  pfac_build(pfacl->fac);
  pfacl->fac->deg = nrem_roots;
  pfacl->fac->labs = new int[nrem_roots];
  pfac_set_coeffs(pfacl->fac, coeffs, nrem_roots);
  // pfac_set_mults(pfacl->fac, mults_rem, den->nroots);  // #next: move these deflation mults into root finder
  if (dbg) {
  cout << "factor coeffs:" << endl;
  for (int k=0; k<=pfacl->fac->deg; k++) {
    cout << "k = " << k << ": ";
    mpq_out_str(stdout, 10, pfacl->fac->coeffs[k]); cout << endl;
  }
  }

  // append new factor to denominator
  den->facs[den->nfacs] = facs->nfacs;
  den->fmults[den->nfacs++] = 1;
  facs->nfacs++;
  if (dbg) cout << "nfacs = " << facs->nfacs << endl;

  if (nrem_roots == nprev_roots) {
    //////
    // SPLIT PARTIALLY INCLUDED FACTORS
    //////
    pfac_set_mults(pfacl->fac, mults_rem, den->nroots);
    if (pfaclist_split(facs, pfacl->fac, *roots)) {
      return 1;
    } else {
      // cout << "new factor: nfacs = " << facs->nfacs << endl;
      // pfac_print(pfacl->fac);
      // pfac_check(pfacl->fac, *roots);
      // cout << "cumulated roots: nroots = " << *nroots << endl;
      // print_poly(*roots, *nroots-1);
      return 0;
    };
  }

  // if (dbg) {
  // cout << "factor multiplicites:" << endl;
  // for (int k=0; k<pfacl->fac->deg; k++) {
  //   cout << pfacl->fac->labs[k] << endl;
  // }
  // }

  //////
  // FIND NEW ROOTS
  //////
  if (dbg) cout << "FIND NEW ROOTS" << endl;

  // HANDLE MULTIPLICITIES OF PREVIOUS ROOTS
  int *root_prof = new int[deg];
  int count = 0;
  for (int k=1; k<*nroots; k++) {
    for (int m=0; m<mults_rem[k]; m++) {
      root_prof[count] = k;
      count++;
    }
  }

  // FIND ROOTS FROM RATIONAL COEFFS
  if (dbg) cout << "FIND ROOTS FROM RATIONAL COEFFS" << endl;
  mpc_t *den_roots;
  mpfr_t *den_tols;
  int den_vdeg;
  int found_mul_roots = 0;
  poly_mpq_find_roots(
    &found_mul_roots,
    &den_vdeg, &den_roots, &den_tols,
    coeffs, deg,
    wp2, mpfr_tol, 0, 1,
    *nroots, mults_rem, *roots
  );
  // int num_new_roots_max = den_vdeg - nprev_roots;
  int num_new_roots_max = den_vdeg;
  if (dbg) {
  cout << "den_vdeg = " << den_vdeg << endl;
  cout << "max num new roots = " << num_new_roots_max << endl;
  cout << "roots from rational coeffs:" << endl;
  print_poly(den_roots+1, den_vdeg-1);
  }
  
  // if (found_mul_roots) {
  //   nprev_roots = 0;
  // }

  // update list of global roots and find root profile
  if (dbg) cout << "update list of global roots and find root profile" << endl;
  int num_new_roots = 0;
  mpc_t *new_roots = new mpc_t[num_new_roots_max+1];
  mpfr_t *new_tols = new mpfr_t[num_new_roots_max+1];
  // for (int k=1; k<=num_new_roots_max; k++) {
  for (int k=nprev_roots+1; k<=den_vdeg; k++) {
    update_new_roots(
      &root_prof[k-1], &num_new_roots, new_roots, new_tols,
      nroots, *roots, *tols,
      &den_roots[k] , &den_tols[k]
    );
  }
  if (dbg) {
  cout << "root profile:" << endl;
  for (int k=0; k<deg; k++) {
    cout << "k, rp = " << k << ", " << root_prof[k] << endl;
  }
  cout << "roots: " << endl;
  print_poly(*roots, *nroots-1);
  cout << "new roots: " << endl;
  print_poly(new_roots, num_new_roots-1);
  }

  mpc_rk1_append(
    nroots, roots,
    num_new_roots, new_roots
  );
  mpfr_rk1_append(
    nroots, tols,
    num_new_roots, new_tols
  );
  *nroots += num_new_roots;
  if (dbg) cout << "nroots = " << *nroots << endl;
  if (dbg) print_poly(*roots, *nroots-1);
  
  // update multiplicites
  if (dbg) cout << "update mults:" << endl;
  root_prof_update_mults(den->mults, root_prof+nprev_roots, deg-nprev_roots, *nroots);
  root_prof_update_fac(
    pfacl->fac, 0,
    root_prof, deg, *nroots
  );

  // cout << "new factor: nfacs = " << facs->nfacs << endl;
  // pfac_print(pfacl->fac);
  // pfac_check(pfacl->fac, *roots);
  // cout << "cumulated roots: nroots = " << *nroots << endl;
  // print_poly(*roots, *nroots-1);

  // cout << "den:" << endl;
  // pden_print(den);

  // FREE
  delete[] mults_rem;
  delete[] fac_incl;
  delete[] root_prof;
  mpc_rk1_clear(den_roots, den_vdeg+1);
  delete[] den_roots;
  mpfr_rk1_clear(den_tols, den_vdeg+1);
  delete[] den_tols;
  mpc_rk1_clear(new_roots, num_new_roots);
  delete[] new_roots;
  mpfr_rk1_clear(new_tols, num_new_roots);
  delete[] new_tols;

  return 0;

}


typedef struct rpfrac {
  int num_deg;
  mpq_t *coeffs;
  pden den;
} rpfrac;


void rpfrac_build(
  rpfrac *rpf
) {
  rpf->coeffs = NULL;
  pden_build(&rpf->den);
}


void rpfrac_free(
  rpfrac *rpf
) {
  if (rpf->coeffs) {
    for (int k=0; k<=rpf->num_deg; k++) {
      mpq_clear(rpf->coeffs[k]);
    }
    delete[] rpf->coeffs;
  }
  pden_free(&rpf->den);
}


void rpfrac_print(
  rpfrac *rpf
) {
  if (rpf->num_deg == -1) {
    cout << "0" << endl;
    return;
  }

  cout << "num deg: " << rpf->num_deg << endl;
  if (rpf->coeffs) {
    for (int k=0; k<=rpf->num_deg; k++) {
      cout << "pow " << k << ": ";
      mpq_out_str(stdout, 10, rpf->coeffs[k]); cout << endl;
    }
  }

  cout << "pden:" << endl;
  pden_print(&rpf->den);
  
}


void rpfrac_neg(
  rpfrac *rpf
) {
  mpq_neg(rpf->den.norm, rpf->den.norm);
}


void rpfrac_set_zero(
  rpfrac *rpf
) {
  rpf->num_deg = -1;
  pden_set_ui(&rpf->den, 1);
}


void rpfrac_set_ui(
  rpfrac *rpf,
  int ui_num,
  int nroots
) {
  if (ui_num == 0) {
    rpfrac_set_zero(rpf);
    return;
  }
  (*rpf).num_deg = 0;
  (*rpf).den.deg = 0;
  (*rpf).den.nroots = nroots;
  (*rpf).den.nfacs = 0;
  if ((*rpf).den.facs) {
    delete[] (*rpf).den.facs;
    (*rpf).den.facs = NULL;
  }
  if ((*rpf).den.fmults) {
    delete[] (*rpf).den.fmults;
    (*rpf).den.fmults = NULL;
  }
  if ((*rpf).den.nroots > 0) {
    if (!(*rpf).den.mults) {
      (*rpf).den.mults = new int[nroots];
    }
    for (int k=0; k<nroots; k++) {
      (*rpf).den.mults[k] = 0;
    }
  }

  // COEFFICIENTS
  if ((*rpf).coeffs) {
    delete[] (*rpf).coeffs;
  }
  (*rpf).coeffs = new mpq_t[1];
  mpq_init((*rpf).coeffs[0]);
  mpq_set_ui((*rpf).coeffs[0], ui_num, 1);

  // NORMALIZATION
  mpq_set_ui(rpf->den.norm, 1, 1);
}


void rpfrac_set_si(
  rpfrac *rpf,
  int si_num,
  int nroots
) {
  if (si_num >=0) {
    rpfrac_set_ui(rpf, si_num, nroots);
  } else {
    rpfrac_set_ui(rpf, -si_num, nroots);
    rpfrac_neg(rpf);
  }
}


void rpfrac_set_rpf(
  rpfrac *out, rpfrac *in
) {
  // check if input is equal to output
  if (out == in) {
    return;
  }

  // check if input is zero
  if (in->num_deg == -1) {
    rpfrac_set_zero(out);
    return;
  }

  // SET NUMERATOR
  out->num_deg = in->num_deg;
  if (out->coeffs) {
    delete[] out->coeffs;
  }
  out->coeffs = new mpq_t[in->num_deg+1];
  for (int k=0; k<=in->num_deg; k++) {
    mpq_init(out->coeffs[k]);
    mpq_set(out->coeffs[k], in->coeffs[k]);
  }

  // SET DENOMINATOR
  pden_set_pden(&out->den, &in->den);

}


void rpfrac_add_rpf(
  // OUTPUT
  rpfrac *out,
  // INPUT
  rpfrac *in1, rpfrac *in2,
  pfaclist *facl, pfac **fac_addr
) {
  // check whether 1st is zero
  if (in1->num_deg == -1) {
    // return the 2nd
    if (out != in2) {
      rpfrac_set_rpf(out, in2);
    }
    return;
  }
  // check whether 2nd is zero
  if (in2->num_deg == -1) {
    // return the 1st
    if (out != in1) {
      rpfrac_set_rpf(out, in1);
    }
    return;
  }

  // deal with overwriting
  int ovwrt, *facs = NULL, *fmults = NULL;
  if (out == in1 || out == in2) {
    ovwrt = 1;
  } else {
    ovwrt = 0;
    facs = out->den.facs;
    fmults = out->den.fmults;
  }
  if (facs) {
    delete[] facs;
  }
  facs = new int[in1->den.nfacs+in2->den.nfacs];
  if (fmults) {
    delete[] fmults;
  }
  fmults = new int[in1->den.nfacs+in2->den.nfacs];

  //////
  // COMPUTE LCM BETWEEN DENOMINATORS
  //////
  // cout << "COMPUTE LCM BETWEEN DENOMINATORS" << endl;
  int nfacs = in1->den.nfacs;
  int *common = new int[in2->den.nfacs];
  for (int f=0; f<in2->den.nfacs; f++) {
    common[f] = -1;
  }
  int *fmults_mul1 = new int[facl->nfacs];
  int *fmults_mul2 = new int[facl->nfacs];
  int num_deg1 = in1->num_deg, num_deg2 = in2->num_deg;
  for (int f=0; f<facl->nfacs; f++) {
    fmults_mul1[f] = 0;
    fmults_mul2[f] = 0;
  }
  for (int f=0; f<in1->den.nfacs; f++) {
    // cout << "f = " << f << endl;
    facs[f] = in1->den.facs[f];
    fmults[f] = in1->den.fmults[f];
    fmults_mul2[facs[f]] = fmults[f];
    for (int fp=0; fp<in2->den.nfacs; fp++) {
      if (in2->den.facs[fp] == facs[f]) {
        common[fp] = f;
        // cout << "common[" << fp << "]" << endl;
        if (in2->den.fmults[fp] >= fmults[f]) {
          fmults[f] = in2->den.fmults[fp];
          fmults_mul2[facs[f]] = 0;
        } else {
          fmults_mul2[facs[f]] -= in2->den.fmults[fp];
        }
        break;
      }
    }
    num_deg2 += fmults_mul2[facs[f]]*fac_addr[facs[f]]->deg;
    // cout << "facs = " << facs[f] << endl;
    // cout << "fmults = " << fmults[f] << endl;
    // cout << "fmults_mul2 = " << fmults_mul2[facs[f]] << endl;
    // cout << "num_deg2 = " << num_deg2 << endl;
  }

  for (int f=0; f<in2->den.nfacs; f++) {
    // cout << "f = " << f << endl;
    if (common[f] >= 0) {
      if (in2->den.fmults[f] > in1->den.fmults[common[f]]) {
        fmults_mul1[in2->den.facs[f]] = in2->den.fmults[f] - in1->den.fmults[common[f]];
        num_deg1 += fmults_mul1[in2->den.facs[f]]*fac_addr[in2->den.facs[f]]->deg;
      } else {
        fmults_mul1[in2->den.facs[f]] = 0;
      }
      continue;
    } else {
      fmults_mul1[in2->den.facs[f]] += in2->den.fmults[f];
      facs[nfacs] = in2->den.facs[f];
      fmults[nfacs] = in2->den.fmults[f];
      num_deg1 += fmults_mul1[in2->den.facs[f]]*fac_addr[facs[nfacs]]->deg;
      // cout << "facs = " << facs[nfacs] << endl;
      // cout << "fmults = " << fmults[nfacs] << endl;
      nfacs++;
    }
    // cout << "fmults_mul1 = " << fmults_mul1[in2->den.facs[f]] << endl;
    // cout << "num_deg1 = " << num_deg1 << endl;
  }

  // LOWER DEGREES
  int ldeg_mul1;
  int ldeg_mul2;
  if (in1->den.ldeg <= in2->den.ldeg) {
    ldeg_mul1 = in2->den.ldeg - in1->den.ldeg;
    ldeg_mul2 = 0;
    out->den.ldeg = in2->den.ldeg;
  } else {
    ldeg_mul1 = 0;
    ldeg_mul2 = in1->den.ldeg - in2->den.ldeg;
    out->den.ldeg = in1->den.ldeg;
  }
  num_deg1 += ldeg_mul1;
  num_deg2 += ldeg_mul2;
  // cout << "num_deg1, num_deg2 = " << num_deg1 << ", " << num_deg2 << endl;

  //////
  // NUMERATORS TIMES FACTORS
  //////
  // cout << "NUMERATORS TIMES FACTORS" << endl;
  
  // 1ST NUMERATOR
  mpq_t *pol11, *pol12;
  pol11 = new mpq_t[num_deg1+1];
  pol12 = new mpq_t[num_deg1+1];
  for (int k=0; k<=num_deg1; k++) {
    mpq_init(pol11[k]);
    mpq_init(pol12[k]);
  }
  mpq_t *pout = pol11, *num1 = pol12, *ptmp;
  
  for (int k=0; k<ldeg_mul1; k++) {
    mpq_set_ui(num1[k], 0, 1);
  }
  for (int k=ldeg_mul1; k<=ldeg_mul1+in1->num_deg; k++) {
    mpq_set(num1[k], in1->coeffs[k-ldeg_mul1]);
  }

  int count = ldeg_mul1+in1->num_deg;
  int nmuls = 0;
  pfaclist *cursor = facl;
  for (int f=0; f<facl->nfacs; f++) {
    for (int m=0; m<fmults_mul1[f]; m++) {
      // cout << "f, m = " << f << ", " << m << endl;
      // cout << "1st fac:" << endl;
      // for (int k=0; k<=count; k++) {
      //   mpq_out_str(stdout, 10, num1[k]); cout << endl;
      // }
      // cout << "2nd fac:" << endl;
      // pfac_print(cursor->fac);
      mul_poly_mpq_tot_wrp(pout, num1, cursor->fac->coeffs, count, cursor->fac->deg);
      count += cursor->fac->deg;
      nmuls++;
      // swap pointers
      ptmp = pout;
      pout = num1;
      num1 = ptmp;
    }
    cursor = cursor->next;
  }
  // cout << "1st num:" << endl;
  // for (int k=0; k<=num_deg1; k++) {
  //   mpq_out_str(stdout, 10, num1[k]); cout << endl;
  // }

  // 2nd NUMERATOR
  // cout << "2nd NUMERATOR" << endl;
  mpq_t *pol21 = new mpq_t[num_deg2+1];
  mpq_t *pol22 = new mpq_t[num_deg2+1];
  for (int k=0; k<=num_deg2; k++) {
    mpq_init(pol21[k]);
    mpq_init(pol22[k]);
  }
  pout = pol21;
  mpq_t *num2 = pol22;

  for (int k=0; k<ldeg_mul2; k++) {
    mpq_set_ui(num2[k], 0, 1);
  }
  for (int k=ldeg_mul2; k<=ldeg_mul2+in2->num_deg; k++) {
    mpq_set(num2[k], in2->coeffs[k-ldeg_mul2]);
  }

  count = ldeg_mul2+in2->num_deg;
  nmuls = 0;
  cursor = facl;
  for (int f=0; f<facl->nfacs; f++) {
    for (int m=0; m<fmults_mul2[f]; m++) {
      // cout << "f, m = " << f << ", " << m << endl;
      // cout << "1st fac:" << endl;
      // for (int k=0; k<=count; k++) {
      //   mpq_out_str(stdout, 10, num2[k]); cout << endl;
      // }
      // cout << "2nd fac:" << endl;
      // pfac_print(cursor->fac);
      mul_poly_mpq_tot_wrp(pout, num2, cursor->fac->coeffs, count, cursor->fac->deg);
      count += cursor->fac->deg;
      nmuls++;
      // swap pointers
      ptmp = pout;
      pout = num2;
      num2 = ptmp;
    }
    cursor = cursor->next;
  }
  // cout << "2nd num:" << endl;
  // for (int k=0; k<=num_deg2; k++) {
  //   mpq_out_str(stdout, 10, num2[k]); cout << endl;
  // }

  //////
  // ADD NUMERATORS
  //////
  // cout << "ADD NUMERATORS" << endl;
  int mindeg, maxdeg;
  mpq_t norm_tmp;
  mpq_init(norm_tmp);
  if (num_deg1 <= num_deg2) {
    mindeg = num_deg1;
    maxdeg = num_deg2;
    // normalization factor
    mpq_div(norm_tmp, in2->den.norm, in1->den.norm);
    for (int k=0; k<=mindeg; k++) {
      mpq_mul(num1[k], num1[k], norm_tmp);
    }
    mpq_set(out->den.norm, in2->den.norm);
  } else {
    mindeg = num_deg2;
    maxdeg = num_deg1;
    // normalization factor
    mpq_div(norm_tmp, in1->den.norm, in2->den.norm);
    for (int k=0; k<=mindeg; k++) {
      mpq_mul(num2[k], num2[k], norm_tmp);
    }
    mpq_set(out->den.norm, in1->den.norm);
  }
  // cout << "norm 1 = "; mpq_out_str(stdout, 10, in1->den.norm); cout << endl;
  // cout << "norm 2 = "; mpq_out_str(stdout, 10, in2->den.norm); cout << endl;
  // cout << "norm tmp = "; mpq_out_str(stdout, 10, norm_tmp); cout << endl;
  if (out->coeffs) {
    mpq_rk1_clear(out->coeffs, out->num_deg+1);
    delete[] out->coeffs;
  }
  out->coeffs = new mpq_t[maxdeg+1];

  // common powers
  // cout << "common powers" << endl;
  for (int k=0; k<=mindeg; k++) {
    mpq_init(out->coeffs[k]);
    // normalization factors
    mpq_add(out->coeffs[k], num1[k], num2[k]);
    // mpq_out_str(stdout, 10, num1[k]); cout << " + ";
    // mpq_out_str(stdout, 10, num2[k]); cout << " = ";
    // mpq_out_str(stdout, 10, out->coeffs[k]); cout << endl;
  }

  // higher powers
  if (maxdeg == num_deg1) {
    for (int k=mindeg+1; k<=maxdeg; k++) {
      mpq_init(out->coeffs[k]);
      mpq_set(out->coeffs[k], num1[k]);
    }
  } else if (maxdeg == num_deg2) {
    for (int k=mindeg+1; k<=maxdeg; k++) {
      mpq_init(out->coeffs[k]);
      mpq_set(out->coeffs[k], num2[k]);
    }
  }

  //////
  // ASSIGN DENOMINATOR DEG, NFACS, FACS, FMULTS, MULTS AND CLEAN
  //////
  out->num_deg = maxdeg;
  out->den.nfacs = nfacs;
  out->den.deg = out->den.ldeg;
  for (int f=0; f<nfacs; f++) {
    out->den.deg += fmults[f]*fac_addr[facs[f]]->deg;
  }

  // assign facs and fmults if overwriting
  if (ovwrt) {
    if (out->den.facs) {
      delete[] out->den.facs;
    }
    out->den.facs = new int[nfacs];
    if (out->den.fmults) {
      delete[] out->den.fmults;
    }
    out->den.fmults = new int[nfacs];

    for (int f=0; f<nfacs; f++) {
      out->den.facs[f] = facs[f];
      out->den.fmults[f] = fmults[f];
    }
    // delete
    delete[] facs;
    delete[] fmults;
  }

  // FREE
  delete[] common;
  delete[] fmults_mul1;
  delete[] fmults_mul2;
  mpq_rk1_clear(pol11, num_deg1+1);
  delete[] pol11;
  mpq_rk1_clear(pol12, num_deg1+1);
  delete[] pol12;
  mpq_rk1_clear(pol21, num_deg2+1);
  delete[] pol21;
  mpq_rk1_clear(pol22, num_deg2+1);
  delete[] pol22;
  mpq_clear(norm_tmp);
  
}


void rpfrac_mul_ui_ovrwrt(
  rpfrac *rpf,
  int ui
) {
  mpq_t mpq_ui;
  mpq_init(mpq_ui);
  mpq_set_ui(mpq_ui, 1, ui);
  mpq_mul(rpf->den.norm, rpf->den.norm, mpq_ui);
  mpq_clear(mpq_ui);
}


void rpfrac_mul_si_ovrwrt(
  rpfrac *rpf,
  int si
) {
  if (si > 0) {
    rpfrac_mul_ui_ovrwrt(rpf, si);
  } else {
    rpfrac_mul_ui_ovrwrt(rpf, -si);
    rpfrac_neg(rpf);
  }
}


void poly_frac_set_rpf(
  poly_frac *pf, rpfrac *rpf, pfaclist *facl, pfac **fac_addr, int nroots
) {
  //////
  // POWER BEHAVIOR AND DEGREES
  //////
  // detect lower degree of the numerator
  int rpf_num_ldeg = 0;
  for (int k=0; k<rpf->num_deg; k++) {
    if (mpq_sgn(rpf->coeffs[k]) == 0) {
      rpf_num_ldeg++;
    } else {
      break;
    }
  }
  int pow_behav = rpf_num_ldeg - rpf->den.ldeg;
  int num_ldeg, den_ldeg;
  if (pow_behav > 0) {
    num_ldeg = pow_behav;
    den_ldeg = 0;
    pf->num_deg = rpf->num_deg - rpf->den.ldeg;
    pf->den_deg = rpf->den.deg - rpf->den.ldeg;
  } else if (pow_behav <= 0) {
    num_ldeg = 0;
    den_ldeg = -pow_behav;
    pf->num_deg = rpf->num_deg - rpf_num_ldeg;
    pf->den_deg = rpf->den.deg - rpf_num_ldeg;
  }
  pf->num_vdeg = pf->num_deg - num_ldeg;

  if (pf->coeffs) {
    delete[] pf->coeffs;
  }
  pf->coeffs = new mpc_t[pf->num_vdeg+1];

  //////
  // NORMALIZATION CONSTANT
  //////
  mpq_t norm;
  mpq_init(norm);
  mpq_set(norm, rpf->den.norm);
  for (int f=0; f<rpf->den.nfacs; f++) {
    for (int m=0; m<rpf->den.fmults[f]; m++) {
      mpq_mul(norm, norm, fac_addr[rpf->den.facs[f]]->coeffs[fac_addr[rpf->den.facs[f]]->deg]);
    }
  }

  //////
  // NUMERATOR COEFFICIENTS
  //////
  mpq_t tmpq;
  mpq_init(tmpq);
  for (int k=0; k<=pf->num_vdeg; k++) {
    mpc_init3(pf->coeffs[k], wp2, wp2);
    mpq_div(tmpq, rpf->coeffs[rpf_num_ldeg+k], norm);
    mpc_set_q(pf->coeffs[k], tmpq, MPFR_RNDN);
    // mpfr_set_ui(mpc_imagref(pf->coeffs[k]), 0, MPFR_RNDN);
  }

  //////
  // MULTIPLICITIES
  //////
  pf->nroots = nroots;
  if (pf->mults) {
    delete[] pf->mults;
  }
  pf->mults = new int[nroots];
  for (int k=0; k<nroots; k++) {
    pf->mults[k] = 0;
  }

  int *fmults_glob = new int[facl->nfacs];
  for (int f=0; f<facl->nfacs; f++) {
    fmults_glob[f] = 0;
  }
  for (int f=0; f<rpf->den.nfacs; f++) {
    // cout << "f, fac, mult = " << f << ", " << rpf->den.facs[f] << ", " << rpf->den.fmults[f] << endl;
    fmults_glob[rpf->den.facs[f]] += rpf->den.fmults[f];
  }

  pfaclist *cursor = facl;
  for (int f=0; f<facl->nfacs; f++) {
    if (fmults_glob[f] == 0) {cursor = cursor->next; continue;}
    // cout << "f = " << f << endl;
    // pfac_print(cursor->fac);
    // cycle through root labels
    for (int r=0; r<cursor->fac->deg; r++) {
      // cout << "r, lab = " << r << ", " << cursor->fac->labs[r] << endl;
      pf->mults[cursor->fac->labs[r]] += fmults_glob[f];
    }
    cursor = cursor->next;
  }

  // null root
  pf->mults[0] = den_ldeg;

  // FREE
  mpq_clear(norm);
  mpq_clear(tmpq);
  delete[] fmults_glob;
}


// typedef struct lnode_pf {
//   int sign;
//   int info;  // 1: pf is one, -1: pf is zero, 2: other
//   int den_is_one;  // 1: den is one
//   int num_nterms;  // number of terms in the numerator
//   int *num_pows;  // power of terms in the numerator
//   lnode *num_terms;  // term trees in the numerator
//   int den_nterms;  // number of terms in the denominator
//   int *den_pows;  // power of terms in the denominator
//   lnode *den_terms;  // trees in the denominator
// } lnode_pf;


void lnode_pf_build(lnode_pf *ndpf) {
  ndpf->sign = 0;
  ndpf->info = 0;
  ndpf->den_is_one = 1;
  ndpf->num_nterms = 0;
  ndpf->num_pows = NULL;
  ndpf->num_terms = NULL;
  ndpf->den_nterms = 0;
  ndpf->den_pows = NULL;
  ndpf->den_terms = NULL;
}


void lnode_pf_free(lnode_pf *ndpf) {
  if(ndpf->num_pows) {
    free(ndpf->num_pows);
  }  
  if(ndpf->num_terms) {
    for (int k=0; k<ndpf->num_nterms; k++) {
      lnode_free(&ndpf->num_terms[k]);
    }
    free(ndpf->num_terms);
  }
  if (ndpf->den_pows) {
    free(ndpf->den_pows);
  }
  if (ndpf->den_terms) {
    for (int k=0; k<ndpf->den_nterms; k++) {
      lnode_free(&ndpf->den_terms[k]);
    }
    free(ndpf->den_terms);
  }
}


void lnode_poly_coeffs_eval(
  // IN-OUT
  mpq_t **mpq_coeffs, int *deg,
  // INPUT
  lnode *nd_coeffs, int *pows, int nterms, mpq_t *val
) {
  mpq_t *mpq_terms = new mpq_t[nterms];
  *deg = pows[0];
  char *nd_str = NULL;
  if (dbg) cout << "nterms = " << nterms << endl;
  if (dbg) cout << "lnode terms:" << endl;
  for (int k=0; k<nterms; k++) {
    if (dbg) cout << "k, pow = " << k << ", " << pows[k] << ": ";
    if (dbg) {cout << endl; lnode_print(&nd_coeffs[k], NULL); cout << endl;}

    // update max power
    if (*deg < pows[k]) {
      *deg = pows[k];
    }

    // decode coefficient
    mpq_init(mpq_terms[k]);
    if(nd_str) {free(nd_str);}
    nd_str = lnode_to_str(&nd_coeffs[k], (char*)", ");
    decode_tree_mpq(&mpq_terms[k], val, &nd_str, (char*)", ");
    if (dbg) {cout << "  mpq: "; mpq_out_str(stdout, 10, mpq_terms[k]); cout << endl;}
  }
  if(nd_str) {free(nd_str);}

  // cumulate coefficients
  *mpq_coeffs = new mpq_t[*deg+1];
  for (int k=0; k<=*deg; k++) {
    mpq_init((*mpq_coeffs)[k]);
    mpq_set_ui((*mpq_coeffs)[k], 0, 1);
  }
  for (int k=0; k<nterms; k++) {
    mpq_add((*mpq_coeffs)[pows[k]], (*mpq_coeffs)[pows[k]], mpq_terms[k]);
  }

  // FREE
  mpq_rk1_clear(mpq_terms, nterms);
  delete[] mpq_terms;
  

}


void process_kira_IBPs(
  // OUTPUT
  mpz_t ******coeffs_num_den, int ******pows_num_den,
  int *****nterms_num_den, lnode_pf ***der_ndpf,
  char ****kira_str,
  // INPUT
  int ninvs, char **symbols, int *is_mass,
  int dim_eta, int max_num_contr,
  LI* MI_eta, char *dir_amflow,
  FILE *terminal
) {
  int print = 0;
  
  char *topo_name = strdup("MIeta");
  char *kira_results_filepath = NULL;
  join_path(&kira_results_filepath, dir_amflow, (char*)"kira/results/");
  join_path(&kira_results_filepath, kira_results_filepath, topo_name);
  join_path(&kira_results_filepath, kira_results_filepath, (char*)"/kira_target.kira");

  char *tmp_str = NULL;
  size_t len = 0;
  char *line = NULL;
  cout << endl; cout << "reading Kira results from " << kira_results_filepath << endl;
  FILE *fptr = fopen(kira_results_filepath, "r");
  int mi_idx, mi_count, c = 0, line_count = 0;
  int topo_len = strlen(topo_name);

  lnode ******nd_num_den;
  malloc_rk4_tens(nd_num_den, dim_eta, max_num_contr, dim_eta, 2);
  lnode *nd_pointer;
  BitField sign, sign1;

  int m_stop = -2, d_stop = -2, n_stop = -2;
  // int m_pfc_start = 105;
  // int m_pfc_start = 109;
  int m_pfc_start = 0;

  //////
  // READ KIRA RESULTS FILE AND EXTRACT COEFFICIENT TREES
  //////
  cout << endl; cout << "READ KIRA RESULTS FILE AND EXTRACT COEFFICIENT TREES" << endl;
  struct lnode nd; lnode_build(&nd);
  fprintf(terminal, "process Kira results: "); fflush(terminal); usleep(sleep_time);
  for (int m=0; m<dim_eta; m++) {
    if (m > 0) {fprintf(terminal, "\033[15D\033[K");}// fflush(terminal); usleep(sleep_time);}
    fprintf(terminal, "MI %3d /%3d... ", m, dim_eta); fflush(terminal); usleep(sleep_time);
    // if (m == 105 || m == 115) {
    //   cout << "line count = " << line_count << endl;
    //   if (m == 115) exit(0);
    // }
    for (int d=0; d<MI_eta[m].nmass; d++) {
      if (print) {cout << endl; cout << "m, d  = " << m << ", " << d << endl;}

      // read coefficients
      if(tmp_str) {free(tmp_str);}
      tmp_str = (char*) malloc(((MAX_POW_DIGITS+1)*dim_eta+2)*sizeof(char));

      // check whether contribution is a master
      mi_idx = LI_rk1_get_idx(MI_eta[m].mder_contr[d].pows_str, MI_eta, dim_eta);
      if (mi_idx != -1) {
        // contribution is a master
        if (print) cout << "contribution is a master: " << MI_eta[m].mder_contr[d].pows_str << endl;
        if (print) cout << "mi_idx = " << mi_idx << endl;
        der_ndpf[m][d][mi_idx].info = 1;
        continue;  
      }

      // skip first line
      getline(&line, &len, fptr);
      line_count++;
      // cout << "line count = " << line_count << endl;
      // cout << "first line: " << line << endl;
      // check whether target is the expected one
      end_at_char(line+topo_len, '*');
      if (strcmp(line+topo_len, MI_eta[m].mder_contr[d].pows_str) != 0) {
        perror("target is not the expected one");
        cout << "line " << line_count << endl;
        cout << "read: " << line+topo_len << endl;
        cout << "expected: " << MI_eta[m].mder_contr[d].pows_str << endl;
        cout << "from master: " << MI_eta[m].pows_str << endl;
        exit(1);
      }

      mi_count = 0;
      while (getline(&line, &len, fptr)) {
        line_count++;
        // cout << "line count = " << line_count << endl;
        // cout << "mi_count = " << mi_count << endl;
        if (strcmp(line, "\n") == 0 || mi_count == dim_eta) {
          break;
        }
        // continue;
        // cout << "line: "<< line << endl;

        // skip topology name
        c = topo_len;
      
        // read power list
        while (line[c-1] != ']') {
          tmp_str[c-topo_len] = line[c];
          c++;
        }
        tmp_str[c-topo_len] = '\0';
        // cout << "power list: " << tmp_str << endl;

        // find corresponding MI
        mi_idx = LI_rk1_get_idx(tmp_str, MI_eta, dim_eta);
        if (mi_idx == -1) {
          cout << "error for m, d  = " << m << ", " << d << endl;
          perror("master of kira results do not correspond to MIs");
          exit(1);
        }
        // if (m == 51) cout << "found mi_idx = " << mi_idx << endl;
        if (print) {cout << endl; cout << "m, mi_idx, d = " << m << ", " << mi_idx << ", " << d << endl;}
        if (m == m_stop && d == d_stop && mi_idx == n_stop) {dbg = 1; cout << "dbg = 1" << endl;}

        // mi_idx++; continue;

        // go to math expression
        while (line[c] != '*') {
          c++;
        }
        // store math expression
        line[c] = '+';  //change '*' with '+' to cumulate
        end_at_char(line+c, '\n');
        // cout << "cumulate: " << line+c << endl;
        // struct lnode nd;
        // if (!(m == 55 && d == 0 && mi_idx == 1)) {continue;}
        if (print) cout << "math string:" << endl << line+c << endl;
        // if (!(m == 137 && d == 3 && mi_idx == 1)) {continue;}
        // if (m == 105 && d == 0 && mi_idx == 0) {
        //   cout << "math string:" << endl << line+c << endl;
        //   exit(0);
        // }
        // char* file_content = read_file_into_string((char*)"processed_expression.txt");
        // lnode_parse_expression(&nd, file_content+1, (const char**) symbols, ninvs+1, 0);
        lnode_free(&nd);
        lnode_parse_expression(&nd, line+c+1, (const char**) symbols, ninvs+1, 0);
        // exit(0);

        //////
        // EXTRACT DENOMINATOR COEFFICIENTS
        //////
        if (print) cout << "EXTRACT DENOMINATOR COEFFICIENTS" << endl;
        struct lnode *nd_ptr;

        // remove sign
        if (nd.op == '-') {
          der_ndpf[m][d][mi_idx].sign = -1;
          nd_ptr = nd.son;
        } else {
          nd_ptr = &nd;
          der_ndpf[m][d][mi_idx].sign = 1;
        }
        if (print) cout << "sign = " << der_ndpf[m][d][mi_idx].sign << endl;

        int num_den = 0;
        if (nd_ptr->op == '*') {
          if (nd_ptr->son->bro->op == '/') {
            num_den = 1;
          }
        }
        if (print) cout << "num_den = " << num_den << endl;

        if (num_den == 1) {
          // cout << "nd:" << endl; lnode_print(nd.son->bro->son, NULL); cout << endl;
          der_ndpf[m][d][mi_idx].den_is_one = 0,
          lnode_coefficient_list_helper(
            &der_ndpf[m][d][mi_idx].den_terms,
            &der_ndpf[m][d][mi_idx].den_pows, &der_ndpf[m][d][mi_idx].den_nterms,
            nd_ptr->son->bro->son, 1
          );

          if (print) {
          cout << "nterms = " << der_ndpf[m][d][mi_idx].den_nterms << endl;
          for (int k=0; k<der_ndpf[m][d][mi_idx].den_nterms; k++) {
            cout << "k, pow = " << k << ", " << der_ndpf[m][d][mi_idx].den_pows[k] << endl;
            lnode_print(&der_ndpf[m][d][mi_idx].den_terms[k], NULL); cout << endl;
          }
          }

          nd_num_den[m][d][mi_idx][1] = new lnode*[der_ndpf[m][d][mi_idx].den_nterms];
          pows_num_den[m][d][mi_idx][1] = new int*[der_ndpf[m][d][mi_idx].den_nterms];
          nterms_num_den[m][d][mi_idx][1] = new int[der_ndpf[m][d][mi_idx].den_nterms];
          coeffs_num_den[m][d][mi_idx][1] = new mpz_t*[der_ndpf[m][d][mi_idx].den_nterms];
          for (int nt=0; nt<der_ndpf[m][d][mi_idx].den_nterms; nt++) {
            if (print) cout << "nt = " << nt << endl;
            nd_num_den[m][d][mi_idx][1][nt] = NULL;
            if (der_ndpf[m][d][mi_idx].den_terms[nt].op == '-') {
              nd_pointer = der_ndpf[m][d][mi_idx].den_terms[nt].son;
              sign.bit = 1;
            } else {
              nd_pointer = &der_ndpf[m][d][mi_idx].den_terms[nt];
              sign.bit = 0;
            }
            if (nd_pointer->op == '*') {
              if (nd_pointer->n == 1) {
                nd_pointer = nd_pointer->son;
              }
            }

            lnode_coefficient_list_helper(
              &nd_num_den[m][d][mi_idx][1][nt],
              &pows_num_den[m][d][mi_idx][1][nt], &nterms_num_den[m][d][mi_idx][1][nt],
              nd_pointer, 0
            );
            
            if (print) cout << "eps nterms = " << nterms_num_den[m][d][mi_idx][1][nt] << endl;
            coeffs_num_den[m][d][mi_idx][1][nt] = new mpz_t[nterms_num_den[m][d][mi_idx][1][nt]];
            for (int a=0; a<nterms_num_den[m][d][mi_idx][1][nt]; a++) {
              if (print) cout << "a, pow = " << a << ", " << pows_num_den[m][d][mi_idx][1][nt][a] << endl;
              if (print) {lnode_print(&nd_num_den[m][d][mi_idx][1][nt][a], NULL); cout << endl;}
              
              sign1.bit = sign.bit;
              if (nd_num_den[m][d][mi_idx][1][nt][a].op == '-') {
                nd_pointer = nd_num_den[m][d][mi_idx][1][nt][a].son;
                sign1.bit ^= 1;
              } else {
                nd_pointer = &nd_num_den[m][d][mi_idx][1][nt][a];
              }
              if (nd_pointer->op == '*') {
                if (nd_pointer->n == 1) {
                  nd_pointer = nd_pointer->son;
                }
              }
              mpz_init(coeffs_num_den[m][d][mi_idx][1][nt][a]);
              mpz_set(coeffs_num_den[m][d][mi_idx][1][nt][a], nd_pointer->mpz);
              if (sign1.bit) {
                mpz_neg(coeffs_num_den[m][d][mi_idx][1][nt][a], coeffs_num_den[m][d][mi_idx][1][nt][a]);
              }
              if (print) {mpz_out_str(stdout, 10, coeffs_num_den[m][d][mi_idx][1][nt][a]); cout << endl;}

              // free
              lnode_free(&nd_num_den[m][d][mi_idx][1][nt][a]);
            }
            lnode_free(&der_ndpf[m][d][mi_idx].den_terms[nt]);
            free(nd_num_den[m][d][mi_idx][1][nt]);
          }
          free(der_ndpf[m][d][mi_idx].den_terms);
          // free(der_ndpf[m][d][mi_idx].den_pows);
          delete[] nd_num_den[m][d][mi_idx][1];
        } else {
          // cout << "set den to one" << endl;
          der_ndpf[m][d][mi_idx].den_is_one = 1;
        }

        //////
        // EXTRACT NUMERATOR COEFFICIENTS
        //////
        if (print) cout << "EXTRACT NUMERATOR COEFFICIENTS" << endl;
        der_ndpf[m][d][mi_idx].info = 2;
        if (m == m_stop+1 && d == 0) break;
        // mi_count++;
        // continue;
        if (num_den == 1) {
          if (print) {cout << "nd:" << endl; lnode_print(nd_ptr->son, NULL); cout << endl;}
          lnode_coefficient_list_helper(
            &der_ndpf[m][d][mi_idx].num_terms,
            &der_ndpf[m][d][mi_idx].num_pows, &der_ndpf[m][d][mi_idx].num_nterms,
            nd_ptr->son, 1
          );
        } else if (num_den == 0) {
          if (print) {cout << "nd:" << endl; lnode_print(nd_ptr, NULL); cout << endl;}
          lnode_coefficient_list_helper(
            &der_ndpf[m][d][mi_idx].num_terms,
            &der_ndpf[m][d][mi_idx].num_pows, &der_ndpf[m][d][mi_idx].num_nterms,
            nd_ptr, 1
          );
        }

        if (print) {
        cout << "nterms = " << der_ndpf[m][d][mi_idx].num_nterms << endl;
        for (int k=0; k<der_ndpf[m][d][mi_idx].num_nterms; k++) {
          cout << "k, pow = " << k << ", " << der_ndpf[m][d][mi_idx].num_pows[k] << endl;
          lnode_print(&der_ndpf[m][d][mi_idx].num_terms[k], NULL); cout << endl;
        }
        }

        if (der_ndpf[m][d][mi_idx].sign == -1) {
          // change sign
          for (int k=0; k<der_ndpf[m][d][mi_idx].num_nterms; k++) {
            lnode_neg(&der_ndpf[m][d][mi_idx].num_terms[k]);
          }
        }


        nd_num_den[m][d][mi_idx][0] = new lnode*[der_ndpf[m][d][mi_idx].num_nterms];
        pows_num_den[m][d][mi_idx][0] = new int*[der_ndpf[m][d][mi_idx].num_nterms];
        nterms_num_den[m][d][mi_idx][0] = new int[der_ndpf[m][d][mi_idx].num_nterms];
        coeffs_num_den[m][d][mi_idx][0] = new mpz_t*[der_ndpf[m][d][mi_idx].num_nterms];
        for (int nt=0; nt<der_ndpf[m][d][mi_idx].num_nterms; nt++) {
          if (print) cout << "nt = " << nt << endl;
          nd_num_den[m][d][mi_idx][0][nt] = NULL;
          if (der_ndpf[m][d][mi_idx].num_terms[nt].op == '-') {
            nd_pointer = der_ndpf[m][d][mi_idx].num_terms[nt].son;
            sign.bit = 1;
          } else {
            nd_pointer = &der_ndpf[m][d][mi_idx].num_terms[nt];
            sign.bit = 0;
          }
          if (nd_pointer->op == '*') {
            if (nd_pointer->n == 1) {
              nd_pointer = nd_pointer->son;
            }
          }

          lnode_coefficient_list_helper(
            &nd_num_den[m][d][mi_idx][0][nt],
            &pows_num_den[m][d][mi_idx][0][nt], &nterms_num_den[m][d][mi_idx][0][nt],
            nd_pointer, 0
          );
          
          if (print) cout << "eps nterms = " << nterms_num_den[m][d][mi_idx][0][nt] << endl;
          coeffs_num_den[m][d][mi_idx][0][nt] = new mpz_t[nterms_num_den[m][d][mi_idx][0][nt]];
          for (int a=0; a<nterms_num_den[m][d][mi_idx][0][nt]; a++) {
            if (print) cout << "a, pow = " << a << ", " << pows_num_den[m][d][mi_idx][0][nt][a] << endl;
            if (print) {lnode_print(&nd_num_den[m][d][mi_idx][0][nt][a], NULL); cout << endl;}

            sign1.bit = sign.bit;
            if (nd_num_den[m][d][mi_idx][0][nt][a].op == '-') {
              nd_pointer = nd_num_den[m][d][mi_idx][0][nt][a].son;
              sign1.bit ^= 1;
            } else {
              nd_pointer = &nd_num_den[m][d][mi_idx][0][nt][a];
            }
            if (nd_pointer->op == '*') {
              if (nd_pointer->n == 1) {
                nd_pointer = nd_pointer->son;
              }
            }
            mpz_init(coeffs_num_den[m][d][mi_idx][0][nt][a]);
            mpz_set(coeffs_num_den[m][d][mi_idx][0][nt][a], nd_pointer->mpz);
            if (sign1.bit) {
              mpz_neg(coeffs_num_den[m][d][mi_idx][0][nt][a], coeffs_num_den[m][d][mi_idx][0][nt][a]);
            }
            if (print) {mpz_out_str(stdout, 10, coeffs_num_den[m][d][mi_idx][0][nt][a]); cout << endl;}

            // free
            lnode_free(&nd_num_den[m][d][mi_idx][0][nt][a]);
          }
          lnode_free(&der_ndpf[m][d][mi_idx].num_terms[nt]);
          free(nd_num_den[m][d][mi_idx][0][nt]);
        }
        free(der_ndpf[m][d][mi_idx].num_terms);
        // free(der_ndpf[m][d][mi_idx].num_pows);
        delete[] nd_num_den[m][d][mi_idx][0];

        mi_count++;
        continue;  // comment to proceed with poly_frac check
        if (m < m_pfc_start) {continue;}
        kira_str[m][d][mi_idx] = lnode_to_str(&nd, (char*)", ");
      }
      if (m == m_stop+1 && d == 0) break;
    }
    if (m == m_stop+1) break;
  }
  fprintf(terminal, "\033[15D\033[K"); fflush(terminal); usleep(sleep_time);
  fclose(fptr);

  // FREE
  if (kira_results_filepath) free(kira_results_filepath);
  if (line) free(line);
  lnode_free(&nd);
  // for (int m=0; m<dim_eta; m++) {
  //   for (int d=0; d<max_num_contr; d++) {
  //     for (int n=0; n<dim_eta; n++) {
  //       lnode_pf_free(&der_ndpf[m][d][n]);
  //     }
  //   }
  // }
  // del_rk3_tens(der_ndpf, dim_eta, max_num_contr);
  del_rk4_tens(nd_num_den, dim_eta, max_num_contr, dim_eta);
  // return;

}


void kira_IBPs_to_DE_pf(
  // OUTPUT
  poly_frac ***pfmat,
  int *zero_label, int *nroots, mpc_t **roots, mpfr_t **tols,
  // INPUT
  int nroots_branch, mpc_t *roots_branch, mpfr_t *tols_branch,
  int ninvs, char **symbols, int *is_mass,
  poly_frac *pspf,
  int *skip_inv, char ***ep_kin,
  int dim_eta, char ****mats_str,
  int nbranches, int *branch_deg,
  int eps_num, char **eps_str,
  LI* MI_eta, char *dir_amflow,
  int max_num_contr,
  mpz_t ******coeffs_num_den, int ******pows_num_den,
  int *****nterms_num_den, lnode_pf ***der_ndpf,
  char ****kira_str,
  FILE *terminal
) {
  // dbg = 1;
  int print = 0;

  double wp2_rel_decr = 0.93;

  poly_frac tmpf;
  poly_frac_build(&tmpf);

  mpc_t tmpc;
  mpc_init3(tmpc, wp2, wp2);

  // // uncomment to proceed with poly_frac check
  // malloc_rk3_tens(kira_str, dim_eta, max_num_contr, dim_eta);
  // for (int m=0; m<dim_eta; m++) {
  //   for (int d=0; d<max_num_contr; d++) {
  //     for (int n=0; n<dim_eta; n++) {
  //       kira_str[m][d][n] = NULL;
  //     }
  //   }
  // }

  // int m_stop = 142, d_stop = 0, n_stop = 0;
  // int m_stop = 126, d_stop = 0, n_stop = 0;
  // int m_stop = 126, d_stop = 2, n_stop = 0;
  // int m_stop = 114, d_stop = 0, n_stop = 0;
  // int m_stop = 105, d_stop = 3, n_stop = 0;
  // int m_stop = 61, d_stop = 0, n_stop = 0;
  // int m_stop = 55, d_stop = 0, n_stop = 6;
  // int m_stop = 55, d_stop = 0, n_stop = 0;
  // int m_stop = 11, d_stop = 1, n_stop = 1;
  // int m_stop = 9, d_stop = 2, n_stop = 1;
  // int m_stop = 3, d_stop = 0, n_stop = 0;
  int m_stop = -2, d_stop = -2, n_stop = -2;
  // int m_pfc_start = 105;
  // int m_pfc_start = 109;
  int m_pfc_start = 0;

  char *kira_str_cp = NULL;
  dbg = 0;
  for (int ep=0; ep<eps_num; ep++) {
    if (ep > 0) {fprintf(terminal, "\033[22D\033[K");}// fflush(terminal); usleep(sleep_time);}
    fprintf(terminal, "eps value %3d /%3d... ", ep, eps_num-1); fflush(terminal); usleep(sleep_time);
    fprintf(stdout, "\n############################################## ep = %d\n", ep);

		ep_kin[0][0] = strdup(eps_str[ep]);
		// cout << "eps = " << ep_kin[0][0] << endl;
		poly_frac_build(&pspf[0]);
		generate_PS_pf(&pspf[0], ep_kin, 0, 0, 0, 1, wp2);

    // rational space-time dimension
    mpq_t *stdim = new mpq_t[1];
    mpq_init(stdim[0]);
    mpq_t mpq_si; 
    mpq_init(mpq_si);
    mpq_set_si(mpq_si, -2, 1);
    mpq_set_str(stdim[0], eps_str[ep], 10);
    mpq_mul(stdim[0], stdim[0], mpq_si);
    mpq_set_si(mpq_si, 4, 1);
    mpq_add(stdim[0], mpq_si, stdim[0]);

    // space-time dimension powers
    int stdim_npows = 200;  // #hard-coded
    mpz_t *stdim_num_pows = new mpz_t[stdim_npows];
    mpz_t *stdim_den_pows = new mpz_t[stdim_npows];
    mpz_init(stdim_num_pows[0]);
    mpz_set_ui(stdim_num_pows[0], 1);
    mpz_init(stdim_den_pows[0]);
    mpz_set_ui(stdim_den_pows[0], 1);
    for (int p=1; p<stdim_npows; p++) {
      mpz_init(stdim_num_pows[p]);
      mpz_mul(stdim_num_pows[p], stdim_num_pows[p-1], mpq_numref(stdim[0]));
      mpz_init(stdim_den_pows[p]);
      mpz_mul(stdim_den_pows[p], stdim_den_pows[p-1], mpq_denref(stdim[0]));
    }

    // reset roots to branch point roots
    nroots[ep] = nroots_branch;
    roots[ep] = new mpc_t[nroots_branch];
    init_rk1_mpc(roots[ep], nroots_branch);
    copy_rk1_mpc(roots[ep], roots_branch, nroots_branch);
    tols[ep] = new mpfr_t[nroots_branch];
    init_rk1_mpfr(tols[ep], nroots_branch);
    copy_rk1_mpfr(tols[ep], tols_branch, nroots_branch);
    // copy for benchmark
    int nroots_bench = nroots_branch;
    mpc_t *roots_bench = new mpc_t[nroots_branch];
    init_rk1_mpc(roots_bench, nroots_branch);
    copy_rk1_mpc(roots_bench, roots_branch, nroots_branch);
    mpfr_t *tols_bench = new mpfr_t[nroots_branch];
    init_rk1_mpfr(tols_bench, nroots_branch);
    copy_rk1_mpfr(tols_bench, tols_branch, nroots_branch);

    // initialize to zero
    for (int m=0; m<dim_eta; m++) {
      for (int n=0; n<dim_eta; n++) {
        poly_frac_set_zero(&pfmat[ep][m][n]);
      }
    }

    pfaclist facs;
    pfaclist_build(&facs);
    rpfrac ***rpf_der;
    malloc_rk3_tens(rpf_der, dim_eta, max_num_contr, dim_eta);

    // initialize every possible contribution
    for (int m=0; m<dim_eta; m++) {
      for (int d=0; d<MI_eta[m].nmass; d++) {
        for (int n=0; n<dim_eta; n++) {
          // cout << "m, d, n = " << m << ", " << d << ", " << n << endl;
          rpfrac_build(&rpf_der[m][d][n]);
          rpf_der[m][d][n].num_deg = -1;
          rpf_der[m][d][n].den.deg = -1;
        }
      }
    }

    // needed for check with poly_frac
    poly_frac ***pf_der = NULL;
    // // uncomment to proceed with poly_frac check
    // int wp2_orig = wp2;
    // mpfr_t mpfr_tol_orig; mpfr_init2(mpfr_tol_orig, wp2);
    // mpfr_set(mpfr_tol_orig, mpfr_tol, MPFR_RNDN);
    // wp2 *= 3;
    // mpfr_clear(mpfr_tol); mpfr_tol_set();
    // malloc_rk3_tens(pf_der, dim_eta, max_num_contr, dim_eta);
    // for (int m=0; m<dim_eta; m++) {
    //   for (int d=0; d<MI_eta[m].nmass; d++) {
    //     for (int n=0; n<dim_eta; n++) {
    //       poly_frac_build(&pf_der[m][d][n]);
    //       poly_frac_set_zero(&pf_der[m][d][n]);
    //     }
    //   }
    // }
    
    dbg = 0;
    int num_split, *split_idx, den_deg;
    mpq_t *den_mpq_coeffs = NULL;
    rpfrac **rpf_mat;
    malloc_rk2_tens(rpf_mat, dim_eta, dim_eta);
    poly_frac pf_check;
    poly_frac_build(&pf_check);
    int roots_bench_cp = 0;
    mpz_t *sorted_coeffs_num;
    mpz_t *sorted_coeffs_den;
    int sorted_deg_num;
    int sorted_deg_den;
    mpq_t *unsorted_coeffs_num;
    mpq_t *unsorted_coeffs_den;
    int unsorted_deg_num;
    int unsorted_deg_den;
    for (int m=0; m<dim_eta; m++) {
      if (m > 0) {fprintf(terminal, "\033[15D\033[K");}// fflush(terminal); usleep(sleep_time);}
      fprintf(terminal, "MI %3d /%3d... ", m, dim_eta); fflush(terminal); usleep(sleep_time);
      // CONVERT DERIVATIVE CONTRIBUTIONS
      // cout << "CONVERT DERIVATIVE CONTRIBUTIONS" << endl;      
      for (int d=0; d<MI_eta[m].nmass; d++) {
        for (int n=0; n<dim_eta; n++) {
          if (dbg) {break;}
          if (print) {cout << endl; cout << "m, d, n  = " << m << ", " << d << ", " << n << endl;}
          if (m == m_stop && d == d_stop && n == n_stop) {dbg = 1; cout << "dbg = 1" << endl;}
          // if (m == 55 && d == 0 && n == 0) {dbg = 1; cout << "dbg = 1" << endl;}
          // if (!(m == 55 && d == 0 && n == 1)) {continue;}
          if (der_ndpf[m][d][n].info == -1) {
            continue;
          } else if (der_ndpf[m][d][n].info == 1) {
            rpfrac_set_si(&rpf_der[m][d][n], 1, nroots[ep]);
            continue;  // comment to proceed with poly_frac check
            if (m < m_pfc_start) {continue;}
            poly_frac_set_si(&pf_der[m][d][n], 1, nroots[ep]);
            continue;
          }

          //////
          // DECODE NUMERATOR
          //////
          if (print || dbg) cout << "DECODE NUMERATOR" << endl;
          if (der_ndpf[m][d][n].info == 2) {
            // find degree and decode terms
            // lnode_poly_coeffs_eval(
            //   &rpf_der[m][d][n].coeffs, &rpf_der[m][d][n].num_deg,
            //   der_ndpf[m][d][n].num_terms,
            //   der_ndpf[m][d][n].num_pows, der_ndpf[m][d][n].num_nterms, stdim
            // );

            unsorted_coeffs_num = new mpq_t[der_ndpf[m][d][n].num_nterms];
            for (int nt=0; nt<der_ndpf[m][d][n].num_nterms; nt++) {
              // SORT STDIM POWERS
              mpz_rk1_sort_pows(
                &sorted_deg_num, &sorted_coeffs_num,
                nterms_num_den[m][d][n][0][nt], coeffs_num_den[m][d][n][0][nt], pows_num_den[m][d][n][0][nt]
              );
              
              // SUM
              mpq_init(unsorted_coeffs_num[nt]);
              mpz_poly_eval_mpq(
                &unsorted_coeffs_num[nt],
                sorted_deg_num, sorted_coeffs_num, stdim_num_pows, stdim_den_pows
              );

              for (int k=0; k<=sorted_deg_num; k++) {
                mpz_clear(sorted_coeffs_num[k]);
              }
              delete[] sorted_coeffs_num;

            }

            if (print) {
            cout << "numerator rational terms:" << endl;
            for (int k=0; k<der_ndpf[m][d][n].num_nterms; k++) {
              cout << "k = " << k << ": ";
              mpq_out_str(stdout, 10, unsorted_coeffs_num[k]); cout << endl;
            }
            }

            // SORT ETA POWERS
            mpq_rk1_sort_pows(
              &rpf_der[m][d][n].num_deg, &rpf_der[m][d][n].coeffs,
              der_ndpf[m][d][n].num_nterms, unsorted_coeffs_num, der_ndpf[m][d][n].num_pows
            );

            if (print) {
            cout << "numerator rational coeffs:" << endl;
            for (int k=0; k<=rpf_der[m][d][n].num_deg; k++) {
              cout << "k = " << k << ": ";
              mpq_out_str(stdout, 10, rpf_der[m][d][n].coeffs[k]); cout << endl;
            }
            }

            for (int nt=0; nt<der_ndpf[m][d][n].num_nterms; nt++) {
              mpq_clear(unsorted_coeffs_num[nt]);
            }
            delete[] unsorted_coeffs_num;
          }

          //////
          // DECODE DENOMINATOR
          //////
          if (print || dbg) cout << "DECODE DENOMINATOR" << endl;
          // if (m == 142 && d == 0 && n == 0) {
          //   dbg = 1;
          // } else {
          //   dbg = 0;
          // }
          if (der_ndpf[m][d][n].den_is_one == 1) {
            // set to one
            pden_set_ui(&rpf_der[m][d][n].den, 1);
          } else if (der_ndpf[m][d][n].den_is_one == 0) {
            // find degree and decode terms
            // lnode_poly_coeffs_eval(
            //   &den_mpq_coeffs, &den_deg,
            //   der_ndpf[m][d][n].den_terms,
            //   der_ndpf[m][d][n].den_pows, der_ndpf[m][d][n].den_nterms, stdim
            // );


            unsorted_coeffs_den = new mpq_t[der_ndpf[m][d][n].den_nterms];
            for (int nt=0; nt<der_ndpf[m][d][n].den_nterms; nt++) {
              // SORT STDIM POWERS
              mpz_rk1_sort_pows(
                &sorted_deg_den, &sorted_coeffs_den,
                nterms_num_den[m][d][n][1][nt], coeffs_num_den[m][d][n][1][nt], pows_num_den[m][d][n][1][nt]
              );
              
              // SUM
              mpq_init(unsorted_coeffs_den[nt]);
              mpz_poly_eval_mpq(
                &unsorted_coeffs_den[nt],
                sorted_deg_den, sorted_coeffs_den, stdim_num_pows, stdim_den_pows
              );

              for (int k=0; k<=sorted_deg_den; k++) {
                mpz_clear(sorted_coeffs_den[k]);
              }
              delete[] sorted_coeffs_den;

            }

            if (print) {
            cout << "denominator rational terms:" << endl;
            for (int k=0; k<der_ndpf[m][d][n].den_nterms; k++) {
              cout << "k = " << k << ": ";
              mpq_out_str(stdout, 10, unsorted_coeffs_den[k]); cout << endl;
            }
            }

            // SORT ETA POWERS
            mpq_rk1_sort_pows(
              &den_deg, &den_mpq_coeffs,
              der_ndpf[m][d][n].den_nterms, unsorted_coeffs_den, der_ndpf[m][d][n].den_pows
            );

            if (print) {
            cout << "denominator rational coeffs:" << endl;
            for (int k=0; k<=den_deg; k++) {
              cout << "k = " << k << ": ";
              mpq_out_str(stdout, 10, den_mpq_coeffs[k]); cout << endl;
            }
            }

            for (int nt=0; nt<der_ndpf[m][d][n].den_nterms; nt++) {
              mpq_clear(unsorted_coeffs_den[nt]);
            }
            delete[] unsorted_coeffs_den;

            if (dbg) {
            cout << "cumulated roots:" << endl;
            print_rk1_mpc(roots[ep], nroots[ep]);
            cout << "denominator rational coefficients:" << endl;
            for (int k=0; k<=den_deg; k++) {
              cout << "k = " << k << ": ";
              mpq_out_str(stdout, 10, den_mpq_coeffs[k]); cout << endl;
            }
            }

            denom_coeffs_to_facs(
              &rpf_der[m][d][n].den, &num_split, &split_idx,
              &facs, &nroots[ep], &roots[ep], &tols[ep],
              den_mpq_coeffs, den_deg
            );            
            
            // cout << "denom:" << endl;
            // pden_print(&rpf_der[m][d][n].den);
            // if (m == 55 && d == 0 && n == 0) {dbg = 0; cout << "dbg = 1" << endl;}

            // // update previous denominators
            // // cout << "update previous denominators" << endl;
            // for (int mp=0; mp<=m; mp++) {
            //   for (int dp=0; dp<MI_eta[mp].nmass; dp++) {
            //     for (int np=0; np<dim_eta; np++) {
            //       // cout << "mp, dp, np = " << mp << ", " << dp << ", " << np << endl;
            //       if (mp == m) {
            //         if (dp > d) {
            //           // denominator not processed yet
            //           continue;
            //         }
            //         if (dp == d && np >= n) {
            //           // current denominator is already up to date (np == n)
            //           // or denominator not processed yes (np > n)
            //           continue;
            //         }
            //       }
            //       if (rpf_der[mp][dp][np].den.deg == -1) {
            //         // empty contribution
            //         continue;
            //       }
            //       pden_facs_update(&rpf_der[mp][dp][np].den, split_idx, num_split);
            //     }
            //   }
            // }

            for (int k=0; k<=den_deg; k++) {
              mpq_clear(den_mpq_coeffs[k]);
            }
            delete[] den_mpq_coeffs;
          }

          continue;  // comment to proceed with poly_frac check
          if (m < m_pfc_start) {continue;}
          if (!roots_bench_cp) {
            nroots_bench = nroots[ep];
            roots_bench = new mpc_t[nroots_bench];
            init_rk1_mpc(roots_bench, nroots_bench);
            copy_rk1_mpc(roots_bench, roots[ep], nroots_bench);
            tols_bench = new mpfr_t[nroots_bench];
            init_rk1_mpfr(tols_bench, nroots_bench);
            copy_rk1_mpfr(tols_bench, tols[ep], nroots_bench);
            roots_bench_cp = 1;
          }

          //////
          // POLY_FRAC
          //////
          int old_nroots = nroots_bench;
          kira_str_cp = strdup(kira_str[m][d][n]);
          if (m == 9 && d == 2 && n == 1) {
            cout << "encoded string:" << endl;
            cout << kira_str_cp << endl;
            cout << "roots:" << endl;
            print_rk1_mpc(roots_bench, nroots_bench);
          }
          decode_tree_pf(
						&pf_der[m][d][n],
						&nroots_bench, &roots_bench, &tols_bench,
						pspf,
						ep_kin, skip_inv-1, is_mass-1,
						&kira_str_cp, (char*) ", ",
						// wp2, mpfr_tol, 0
						wp2, mpfr_tol, 1
					);
          if (m == 9 && d == 2 && n == 1) {
            cout << "output poly_frac:" << endl;
            poly_frac_print(&pf_der[m][d][n]);
            cout << "roots:" << endl;
            print_rk1_mpc(roots_bench, nroots_bench);
          }
          if (kira_str_cp) {free(kira_str_cp);}
          // cout << "m, n, d = " << m << ", " << n << ", " << d << endl;
          // cout << "nroots: old, new = " << old_nroots << ", " << nroots[ep] << endl;
          // if (old_nroots != nroots[ep]) {
          //   cout << nroots_bench - old_nroots << " new roots:" << endl;
          //   print_poly(roots_bench+old_nroots, nroots_bench-old_nroots-1);
          // 
          if (dbg) break;
        }
        if (dbg) break;
      }
      dbg = 0;



      // //////
      // // ANTICIPATED FREE MEMORY
      // //////

      // // rpfrac derivative contributions
      // for (int d=0; d<MI_eta[m].nmass; d++) {
      //   for (int n=0; n<dim_eta; n++) {
      //     rpfrac_free(&rpf_der[m][d][n]);
      //   }
      // }
      // del_rk2_tens(rpf_der[m], max_num_contr);
      
      // // rpfrac partial sums
      // for (int n=0; n<dim_eta; n++) {
      //   rpfrac_free(&rpf_mat[m][n]);
      // }
      // delete[] rpf_mat[m];

      // // poly_frac derivative contributions
      // if (pf_der) {
      //   poly_frac_rk2_free(pf_der[m], MI_eta[m].nmass, dim_eta);
      //   del_rk2_tens(pf_der[m], max_num_contr);
      // }
      // continue;

      

      // group denominator roots into factors
      // cout << "group denominator roots into factors" << endl;
      // if (m >= m_pfc_start) {print = 1;}
      for (int d=0; d<MI_eta[m].nmass; d++) {
        for (int n=0; n<dim_eta; n++) {
          if (print) {cout << endl; cout << "m, d, n  = " << m << ", " << d << ", " << n << endl;}
          if (m == m_stop && d == d_stop && n == n_stop) {dbg = 1; cout << "dbg = 1" << endl;}
          if (der_ndpf[m][d][n].info == -1) {
            // empty contribution
            continue;
          } else if (der_ndpf[m][d][n].info == 1) {
            // denominator is one
            pden_roots_append_zero(&rpf_der[m][d][n].den, nroots[ep]);
            continue;
          }
          
          pden_roots_append_zero(&rpf_der[m][d][n].den, nroots[ep]);
          if (dbg) {
          cout << "den before grouping:" << endl;
          pden_print(&rpf_der[m][d][n].den);
          }
          pden_group_roots(&rpf_der[m][d][n].den, &facs);
          if (dbg) {
          cout << "den after grouping:" << endl;
          pden_print(&rpf_der[m][d][n].den);
          }

          continue;  // comment to proceed with poly_frac check
          if (m < m_pfc_start) {continue;}
          // compare with rpfrac
          poly_frac_roots_append_zero(&pf_der[m][d][n], nroots[ep]);
          if (print || dbg) {
          cout << "compare rpfrac with poly_frac..." << endl;
          cout << "benchmark:" << endl;
          poly_frac_print(&pf_der[m][d][n]);
          cout << "rpfrac:" << endl;
          rpfrac_print(&rpf_der[m][d][n]);
          cout << "bench roots:" << endl;
          print_rk1_mpc(roots_bench, nroots_bench);
          }
          //// store addresses of factor list for fast access
          pfac **fac_addr = new pfac*[facs.nfacs];
          pfaclist *cursor = &facs;
          for (int f=0; f<facs.nfacs; f++) {
            fac_addr[f] = cursor->fac;
            cursor = cursor->next;
          }
          poly_frac tmpf;
          poly_frac_build(&tmpf);
          poly_frac_set_rpf(&tmpf, &rpf_der[m][d][n], &facs, fac_addr, nroots[ep]);
          if (print || dbg) cout << "converted poly_frac:" << endl;
          if (print || dbg) poly_frac_print(&tmpf);
          //// check
          poly_frac_neg(&tmpf);
          poly_frac_add_pf(&tmpf, &tmpf, &pf_der[m][d][n], roots_bench, NULL);
          if (print || dbg) {
          cout << "difference:" << endl;
          poly_frac_print(&tmpf);
          }
          poly_frac_prune_tol(&tmpf, 0.50);
          if (tmpf.num_vdeg != -1) {
            printf("rpfrac not corresponding to poly_frac\n");
            cout << "roots:" << endl;
            print_rk1_mpc(roots[ep], nroots[ep]);
            perror("rpfrac not corresponding to poly_frac");
            exit(1);
          }

          // free memory
          poly_frac_free(&tmpf);
          delete[] fac_addr;

          if (dbg) break;
        }
        if (dbg) break;
      }
      // print = 0;

      // store addresses of factor list for fast access
      pfac **fac_addr = new pfac*[facs.nfacs];
      pfaclist *cursor = &facs;
      if (dbg || print) cout << "FACTORS:" << endl;
      for (int f=0; f<facs.nfacs; f++) {
        if (dbg || print) cout << "f = " << f << endl;
        if (dbg || print) pfac_print(cursor->fac);
        fac_addr[f] = cursor->fac;
        cursor = cursor->next;
      }

      //////
      // SUM TOGETHER DERIVATIVE CONTRIBUTIONS
      //////      
      // cout << "SUM TOGETHER DERIVATIVE CONTRIBUTIONS" << endl;
      for (int n=0; n<dim_eta; n++) {
        rpfrac_build(&rpf_mat[m][n]);
        rpfrac_set_zero(&rpf_mat[m][n]);
        poly_frac pf_psum;
        poly_frac_build(&pf_psum);
        poly_frac_set_zero(&pf_psum);
        for (int d=0; d<MI_eta[m].nmass; d++) {
          if (print) cout << "m, n, d = " << m << ", " << n << ", " << d << endl;
          if (m == m_stop && d == d_stop && n == n_stop) {dbg = 1; cout << "dbg = 1" << endl;}
          // multiply times chain rule factor
          rpfrac_mul_si_ovrwrt(&rpf_der[m][d][n], MI_eta[m].pows[MI_eta[m].mder_idx[d]]);
          if (print || dbg) {
          cout << "1st addend:" << endl;
          rpfrac_print(&rpf_mat[m][n]);
          cout << "2nd addend:" << endl;
          rpfrac_print(&rpf_der[m][d][n]);
          }
          rpfrac_add_rpf(
            &rpf_mat[m][n], &rpf_mat[m][n], &rpf_der[m][d][n],
            &facs, fac_addr
          );
          if (print || dbg) {
          cout << "result:" << endl;
          rpfrac_print(&rpf_mat[m][n]);
          }
          continue;  // comment to proceed with poly_frac check
          if (m < m_pfc_start) {continue;}

          poly_frac_set_rpf(&pfmat[ep][m][n], &rpf_mat[m][n], &facs, fac_addr, nroots[ep]);
          if (print || dbg) {
          cout << "converted:" << endl;
          poly_frac_print(&pfmat[ep][m][n]);
          }

          //////
          // POLY_FRAC
          //////
          poly_frac_roots_append_zero(&pf_der[m][d][n], nroots[ep]);
          poly_frac_mul_si(&pf_der[m][d][n], &pf_der[m][d][n], MI_eta[m].pows[MI_eta[m].mder_idx[d]]);
          if (print || dbg) {
          cout << "1st poly_frac addend:" << endl;
          poly_frac_print(&pf_psum);
          cout << "2nd poly_frac addend:" << endl;
          poly_frac_print(&pf_der[m][d][n]);
          }
          poly_frac_add_pf(
            &pf_psum, &pf_psum, &pf_der[m][d][n],
            roots[ep], NULL
          );
          if (print || dbg) {
          cout << "poly_frac result:" << endl;
          poly_frac_print(&pf_psum);
          }

          // check
          poly_frac_neg(&pf_psum);
          poly_frac_add_pf(&pf_check, &pf_psum, &pfmat[ep][m][n], roots[ep], NULL);
          if (print || dbg) {
          cout << "difference:" << endl;
          poly_frac_print(&pf_check);
          }
          poly_frac_prune_tol(&pf_check, 0.50);
          if (pf_check.num_vdeg != -1) {
            printf("rpfrac not corresponding to poly_frac\n");
            cout << "rpfrac:" << endl;
            rpfrac_print(&rpf_mat[m][n]);
            cout << "converted:" << endl;
            poly_frac_print(&pfmat[ep][m][n]);
            cout << "pf:" << endl;
            poly_frac_neg(&pf_psum);
            poly_frac_print(&pf_psum);
            perror("rpfrac not corresponding to poly_frac");
            exit(1);
          }
          poly_frac_neg(&pf_psum);
          
          if (dbg) {break;}
        }
        // cout << "converting..." << endl;
        // cout << "rpfrac:" << endl;
        // rpfrac_print(&rpf_mat[m][n]);
        poly_frac_set_rpf(&pfmat[ep][m][n], &rpf_mat[m][n], &facs, fac_addr, nroots[ep]);
        // cout << "poly_frac:" << endl;
        // poly_frac_print(&pfmat[ep][m][n]);
        if (dbg) {break;}
      }
      
      //////
      // FREE MEMORY
      //////
      delete[] fac_addr;

      // rpfrac derivative contributions
      for (int d=0; d<MI_eta[m].nmass; d++) {
        for (int n=0; n<dim_eta; n++) {
          rpfrac_free(&rpf_der[m][d][n]);
        }
      }
      del_rk2_tens(rpf_der[m], max_num_contr);
      
      // rpfrac partial sums
      for (int n=0; n<dim_eta; n++) {
        rpfrac_free(&rpf_mat[m][n]);
      }
      delete[] rpf_mat[m];

      // poly_frac derivative contributions
      if (pf_der) {
        poly_frac_rk2_free(pf_der[m], MI_eta[m].nmass, dim_eta);
        del_rk2_tens(pf_der[m], max_num_contr);
      }

      if (dbg) break;
    }
    fprintf(terminal, "\033[15D\033[K"); fflush(terminal); usleep(sleep_time);

    // for (int m=0; m<dim_eta; m++) {
    //   for (int d=0; d<max_num_contr; d++) {
    //     for (int n=0; n<dim_eta; n++) {
    //       lnode_pf_free(&der_ndpf[m][d][n]);
    //     }
    //   }
    // }
    // del_rk3_tens(der_ndpf, dim_eta, max_num_contr);
    // return;

    if (mpc_lessthan_tol(roots[ep][0])) {
      mpc_set_ui(roots[ep][0], 0, MPFR_RNDN);
    }

    // needed for check with poly_frac
    // // uncomment to proceed with poly_frac check
    // wp2 = wp2_orig;
    // mpfr_clear(mpfr_tol);  mpfr_init2(mpfr_tol, wp2);
    // mpfr_set(mpfr_tol, mpfr_tol_orig, MPFR_RNDN);
    // exit(0);

    // update nroots globally across the matrix
    for (int i=0; i<dim_eta; i++) {
      for (int j=0; j<dim_eta; j++) {
        poly_frac_roots_append_zero(&pfmat[ep][i][j], nroots[ep]);
      }
    }

		// PRUNE
		// poly_frac_rk2_prune(pfmat[ep], dim_eta, dim_eta, wp2_rel_decr);
    // poly_frac_rk2_normal(pfmat[ep], dim_eta, dim_eta);
    poly_frac_rk2_trim_zero_p(pfmat[ep], dim_eta, dim_eta);
    int wp_bin = - mpfr_log2_int(mpfr_tol);
    mpc_rk1_prune_re_im(roots[ep], wp_bin, nroots[ep]);
    
    // // prune roots
    // int wp_bin = - mpfr_log2_int(mpfr_tol);
    // mpc_rk1_prune_re_im(roots[ep], nroots[ep], wp_bin);

    // zero_label
    zero_label[ep] = -1;
    for (int i=0; i<dim_eta; i++) {
      for (int j=0; j<dim_eta; j++) {
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

    //////
    // CACHE
    //////
    // MATRIX
    char tmp_filepath[MAX_PATH_LEN];
    snprintf(tmp_filepath, sizeof(tmp_filepath), "%s%s%d%s", filepath_cache, "pfmat", ep, ".txt");
    cout << endl; cout << "writing to " << tmp_filepath << endl;
    poly_frac_rk2_to_file(
      tmp_filepath,
      pfmat[ep], dim_eta, dim_eta
    );
    // poly_frac_rk2_free(pfmat[ep], dim_eta, dim_eta);
    // del_rk2_tens(pfmat[ep], dim_eta);

    // ROOTS
    snprintf(tmp_filepath, sizeof(tmp_filepath), "%s%s%d%s", filepath_cache, "roots", ep, ".txt");
    cout << "writing to " << tmp_filepath << endl;
    int_rk0_mpc_rk1_to_file(tmp_filepath, roots[ep], nroots[ep], zero_label[ep]);
    // mpc_rk1_clear(roots[ep], nroots[ep]);
    // delete[] roots[ep];

    // free memory
    mpc_rk1_clear(roots_bench, nroots_branch);
    delete[] roots_bench;
    mpfr_rk1_clear(tols_bench, nroots_branch);
    delete[] tols_bench;

    // for (int m=0; m<dim_eta; m++) {
    //   for (int d=0; d<MI_eta[m].nmass; d++) {
    //     for (int n=0; n<dim_eta; n++) {
    //       rpfrac_free(&rpf_der[m][d][n]);
    //     }
    //   }
    // }
    // del_rk3_tens(rpf_der, dim_eta, max_num_contr);
    // for (int m=0; m<dim_eta; m++) {
    //   for (int n=0; n<dim_eta; n++) {
    //     rpfrac_free(&rpf_mat[m][n]);
    //   }
    // }
    // del_rk2_tens(rpf_mat, dim_eta);
    pfaclist_free(&facs);
    // if (pf_der) {
    //   for (int m=0; m<dim_eta; m++) {
    //     for (int d=0; d<MI_eta[m].nmass; d++) {
    //       for (int n=0; n<dim_eta; n++) {
    //         poly_frac_free(&pf_der[m][d][n]);
    //       }
    //     }
    //   }
    //   del_rk3_tens(pf_der, dim_eta, max_num_contr);
    // }
    delete[] rpf_der;
    delete[] rpf_mat;
    if (pf_der) {delete[] pf_der;}
    poly_frac_free(&pf_check);

    mpz_rk1_clear(stdim_num_pows, stdim_npows);
    delete[] stdim_num_pows;
    mpz_rk1_clear(stdim_den_pows, stdim_npows);
    delete[] stdim_den_pows;
  }
  fprintf(terminal, "\033[22D\033[K"); fflush(terminal); usleep(sleep_time);
  fprintf(terminal, "\033[22D\033[K"); fflush(terminal); usleep(sleep_time);

  // FREE
  if (kira_str) {
    for (int m=0; m<dim_eta; m++) {
      for (int d=0; d<max_num_contr; d++) {
        for (int n=0; n<dim_eta; n++) {
          if (kira_str[m][d][n]) {
            free(kira_str[m][d][n]);
          }
        }
      }
    }
    del_rk3_tens(kira_str, dim_eta, max_num_contr);
  }

  for (int m=0; m<dim_eta; m++) {
    for (int d=0; d<MI_eta[m].nmass; d++) {
      for (int n=0; n<dim_eta; n++) {
        if (der_ndpf[m][d][n].info == -1) {
          continue;
        } else if (der_ndpf[m][d][n].info == 1) {
          continue;
        }

        if (der_ndpf[m][d][n].info == 2) {
          free(der_ndpf[m][d][n].num_pows);
          free(der_ndpf[m][d][n].den_pows);

          for (int nt=0; nt<der_ndpf[m][d][n].num_nterms; nt++) {
            mpz_rk1_clear(coeffs_num_den[m][d][n][0][nt], nterms_num_den[m][d][n][0][nt]);
            delete[] coeffs_num_den[m][d][n][0][nt];
            free(pows_num_den[m][d][n][0][nt]);
          }
          delete[] coeffs_num_den[m][d][n][0];
          delete[] pows_num_den[m][d][n][0];
          delete[] nterms_num_den[m][d][n][0];

          if (der_ndpf[m][d][n].den_is_one) {
            continue;
          }

          for (int nt=0; nt<der_ndpf[m][d][n].den_nterms; nt++) {
            mpz_rk1_clear(coeffs_num_den[m][d][n][1][nt], nterms_num_den[m][d][n][1][nt]);
            delete[] coeffs_num_den[m][d][n][1][nt];
            free(pows_num_den[m][d][n][1][nt]);
          }
          delete[] coeffs_num_den[m][d][n][1];
          delete[] pows_num_den[m][d][n][1];
          delete[] nterms_num_den[m][d][n][1];
        }
      }
    }
  }
  del_rk4_tens(coeffs_num_den, dim_eta, max_num_contr, dim_eta);
  del_rk4_tens(pows_num_den, dim_eta, max_num_contr, dim_eta);
  del_rk4_tens(nterms_num_den, dim_eta, max_num_contr, dim_eta);
  del_rk3_tens(der_ndpf, dim_eta, max_num_contr);

}


void kira_to_DE_str(
  // OUTPUT
  char ***mat,
  // INPUT
  LI* MI_eta, int dim_eta,
  char *kira_results_filepath, char *topo_name
) {
  char *tmp_str = NULL;
  size_t len = 0;
  char *line = NULL;
  FILE *fptr = fopen(kira_results_filepath, "r");
  int mi_idx, mi_count, c = 0, topo_len = strlen(topo_name);
  for (int m=0; m<dim_eta; m++) {
    // initialize to NULL
    for (int n=0; n<dim_eta; n++) {
      mat[m][n] = NULL;
    }

    for (int d=0; d<MI_eta[m].nmass; d++) {
      // cout << "m, d  = " << m << ", " << d << endl;

      if(tmp_str) {free(tmp_str);}
      tmp_str = (char*) malloc(((MAX_POW_DIGITS+1)*dim_eta+2)*sizeof(char));

      // check whether contribution is a master
      mi_idx = LI_rk1_get_idx(MI_eta[m].mder_contr[d].pows_str, MI_eta, dim_eta);
      if (mi_idx != -1) {
        // contribution is a master
        // cout << "contribution is a master" << endl;
        if (mat[m][mi_idx]) {
          // cumulate
          join_path(&mat[m][mi_idx], mat[m][mi_idx], (char*)"+(1)");
        } else {
          // initialize
          // multiply memory times number of propagators in order to lower
          // the number of allocations
          mat[m][mi_idx] = strdup("+(1)");
        }
        // multiply power factor and -1 since eta appears with minus sign
        snprintf(tmp_str, sizeof(tmp_str), "*(%d)", MI_eta[m].pows[MI_eta[m].mder_idx[d]]);
        join_path(&mat[m][mi_idx], mat[m][mi_idx], tmp_str);
        continue;
      }

      // skip first line
      getline(&line, &len, fptr);
      // cout << "first line: " << line << endl;
      
      // read coefficients
      mi_count = 0;
      while (getline(&line, &len, fptr)) {
        // cout << "mi_count = " << mi_count << endl;
        if (strcmp(line, "\n") == 0 || mi_count == dim_eta) {
          break;
        }
        // cout << "line: "<< line << endl;

        // skip topology name
        c = topo_len;
      
        // read power list
        while (line[c-1] != ']') {
          tmp_str[c-topo_len] = line[c];
          c++;
        }
        tmp_str[c-topo_len] = '\0';
        // cout << "power list: " << tmp_str << endl;

        // find corresponding MI
        // for (mi_idx=0; mi_idx<dim_eta; mi_idx++) {
        //   // cout << "candidate : " << MI_eta[mi_idx].pows_str << endl;
        //   if (strcmp(tmp_str, MI_eta[mi_idx].pows_str) == 0) {
        //     // cout << "found" << endl;
        //     break;
        //   }
        // }
        mi_idx = LI_rk1_get_idx(tmp_str, MI_eta, dim_eta);
        if (mi_idx == -1) {
          cout << "error for m, d  = " << m << ", " << d << endl;
          perror("master of kira results do not correspond to MIs");
          exit(1);
        }
        // if (m == 51) cout << "found mi_idx = " << mi_idx << endl;

        // mi_idx++; continue;

        // go to math expression
        while (line[c] != '*') {
          c++;
        }
        // store math expression
        line[c] = '+';  //change '*' with '+' to cumulate
        end_at_char(line+c, '\n');
        // cout << "cumulate: " << line+c << endl;
        if (mat[m][mi_idx]) {
          // cumulate
          join_path(&mat[m][mi_idx], mat[m][mi_idx], line+c);
        } else {
          // initialize
          // multiply memory times number of propagators in order to lower
          // the number of allocations
          mat[m][mi_idx] = (char*) malloc(MI_eta[m].nprop*(strlen(line+c)+1)*sizeof(char));
          strcpy(mat[m][mi_idx], line+c);
        }
        // multiply power factor and -1 since eta appears with minus sign
        snprintf(tmp_str, sizeof(tmp_str), "*(%d)", MI_eta[m].pows[MI_eta[m].mder_idx[d]]);
        join_path(&mat[m][mi_idx], mat[m][mi_idx], tmp_str);
        // if (m == 51) cout << "cumulated math expr:" << endl;
        // if (m == 51) cout << mat[m][mi_idx] << endl;
        mi_count++;
      }
    }
    
    // set to zero every other element
    for (int n=0; n<dim_eta; n++) {
      // if (m == 51) cout << "n = " << n << endl;
      if (!mat[m][n]) {
        // if (m == 51) cout << "set to zero" << endl;
        mat[m][n] = (char*) malloc(2*sizeof(char));
        mat[m][n][0] = '0';
        mat[m][n][1] = '\0';
      }
    }
  }
  fclose(fptr);

}


void call_kira(
  // OUTPUT
  LI **MI_eta, int *dim_eta, int **MI_idx, int *dim,
  // INPUT
  int redo, int opt_kira_parallel,
  char **kin, char *dir_parent, char *dir_amflow,
  FILE *terminal
) {
  char *dir_common = NULL;
  join_path(&dir_common, dir_parent, (char*)"common/");

  char *topo_filepath = NULL;
  join_path(&topo_filepath, dir_common, (char*)"topology.txt");

  nlist ext_mom_nl, loop_mom_nl, mass_nl, prop_nl, mom_cons_nl, sq_mom_nl, kin_inv_nl;
  FILE *fptr = fopen(topo_filepath, "r");

  cout << endl; cout << "read topology from " << topo_filepath << endl;
  nlist_read_file(&ext_mom_nl, (char*)"external-momenta", fptr);
  cout << ext_mom_nl.nitems << " external momenta: ";
  for (int p=0; p<ext_mom_nl.nitems; p++) {cout << ext_mom_nl.item[p].str << ", ";} cout << endl;
  fseek(fptr, 0, SEEK_SET);
  nlist_read_file(&loop_mom_nl, (char*)"loop-momenta", fptr);
  cout << loop_mom_nl.nitems << " loop momenta: ";
  for (int p=0; p<loop_mom_nl.nitems; p++) {cout << loop_mom_nl.item[p].str << ", ";} cout << endl;
  fseek(fptr, 0, SEEK_SET);
  nlist_read_file(&mass_nl, (char*)"masses", fptr);
  cout << mass_nl.nitems << " masses: ";
  for (int p=0; p<mass_nl.nitems; p++) {cout << mass_nl.item[p].str << ", ";} cout << endl;
  fseek(fptr, 0, SEEK_SET);
  nlist_read_file(&prop_nl, (char*)"propagators", fptr);
  cout << prop_nl.nitems << " propagators: ";
  for (int p=0; p<prop_nl.nitems; p++) {cout << prop_nl.item[p].str << ", ";} cout << endl;
  fseek(fptr, 0, SEEK_SET);
  nlist_read_file(&mom_cons_nl, (char*)"momentum-conservation", fptr);
  cout << "momentum conservation: " << mom_cons_nl.str << endl;
  fseek(fptr, 0, SEEK_SET);
  nlist_read_file(&kin_inv_nl, (char*)"kinematic-invariants", fptr);
  cout << "kinematic invariants: " << kin_inv_nl.str << endl;
  fseek(fptr, 0, SEEK_SET);
  nlist_read_file(&sq_mom_nl, (char*)"squared-momenta", fptr);
  cout << "squared momentum: " << sq_mom_nl.str << endl;
  fclose(fptr);

  // DEFINE TOPOLOGY SYMBOLS
  int ntopo_sym = ext_mom_nl.nitems + loop_mom_nl.nitems + mass_nl.nitems;
  char **topo_sym = new char*[ntopo_sym];
  int sym_idx = 0;
  for (int p=0; p<loop_mom_nl.nitems; p++) {
    topo_sym[sym_idx] = (char*) malloc(MAX_SYM_LEN*sizeof(char));
    strcpy(topo_sym[sym_idx], loop_mom_nl.item[p].str);
    sym_idx++;
  }
  for (int p=0; p<ext_mom_nl.nitems; p++) {
    topo_sym[sym_idx] = (char*) malloc(MAX_SYM_LEN*sizeof(char));
    strcpy(topo_sym[sym_idx], ext_mom_nl.item[p].str);
    sym_idx++;
  }
  for (int p=0; p<mass_nl.nitems; p++) {
    topo_sym[sym_idx] = (char*) malloc(MAX_SYM_LEN*sizeof(char));
    strcpy(topo_sym[sym_idx], mass_nl.item[p].str);
    sym_idx++;
  }
  cout << ntopo_sym << " symbols: ";
  for (int s=0; s<ntopo_sym; s++) { cout << topo_sym[s] << ", ";} cout << endl;  

  // DEFINE PROPAGATORS
  // eta-less propagators
  propagator *prop = new propagator[prop_nl.nitems];
  cout << "eta-less propagator representation:" << endl;
  for (int p=0; p<prop_nl.nitems; p++) {
    propagator_build(&prop[p]);
    propagator_set_nlist(&prop[p], &prop_nl.item[p], topo_sym, ntopo_sym, loop_mom_nl.nitems);
    cout << "p = " << p << ": "; propagator_print(&prop[p], topo_sym);  cout << endl;
  }

  // // propagators with eta
  // propagator *prop_eta = new propagator[prop_nl.nitems];
  // nlist prop_eta_nl;
  // cout << "propagator representation:" << endl;
  // for (int p=0; p<prop_nl.nitems; p++) {
  //   propagator_build(&prop_eta[p]);
  //   propagator_copy(&prop_eta[p], &prop[p]);
  //   prop_eta[p].is_massive = 1;
  //   if (prop[p].mass_tree) {
  //   } else {
  //   }
  // }

  // DEFINE KINEMATIC INVARIANTS SYMBOLS
  int ninvs = kin_inv_nl.nitems + mass_nl.nitems;
  cout << "ninvs = " << ninvs << endl;
  char **inv_sym = new char*[ninvs];
  for (int s=0; s<kin_inv_nl.nitems; s++) {
    inv_sym[s] = (char*) malloc((strlen(kin_inv_nl.item[s].str)+1)*sizeof(char));
    strcpy(inv_sym[s], kin_inv_nl.item[s].str);
    // cout << inv_sym[s] << endl;
  }
  for (int s=0; s<mass_nl.nitems; s++) {
    inv_sym[kin_inv_nl.nitems+s] = (char*) malloc((strlen(mass_nl.item[s].str)+1)*sizeof(char));
    strcpy(inv_sym[kin_inv_nl.nitems+s], mass_nl.item[s].str);
    // cout << inv_sym[kin_inv_nl.nitems+s] << endl;
  }

  // DEFINE MIs
  char *MIs_filepath = NULL;
  join_path(&MIs_filepath, dir_common, (char*)"MIs.txt");
  cout << endl; cout << "read MIs without eta from " << MIs_filepath << endl;

  LI *MI = NULL;
  LI_rk1_from_file(MIs_filepath, &MI, dim, (char*)"MI");
  cout << "eta-less MIs:" << endl;
  LI_rk1_pows_print(MI, *dim);

  // identify numerator ISP
  int *nISP = new int[prop_nl.nitems];
  LI_rk1_get_nISP(nISP, MI, *dim);
  cout << "numerator ISPs: [";
  for (int i=0; i<prop_nl.nitems; i++) {
    cout << nISP[i] << ", ";
  }
  cout << "]" << endl;

  // LOAD FIXED KINEMATICS
  mpq_t *PS = new mpq_t[ninvs];
  cout << endl; cout << "eta-less point: [" << endl;
  for (int s=0; s<ninvs; s++) {
    mpq_init(PS[s]);
    if (strchr(kin[s], (".")[0])) {
      // input is real

      // find_cf(3.141592653589793);
      // find_cf(1.41421356237309504880168872421, 1e-16);
      // mpfr_t realin;
      // mpfr_init2(realin, wp2);
      // mpfr_t delta; mpfr_init2(delta, wp2); mpfr_set_d(delta, 1e-50, MPFR_RNDN);
      // mpfr_sqrt_ui(realin, 2, MPFR_RNDN);
      // mpq_t rat; mpq_init(rat);
      // rationalize_cf(rat, realin, delta);
      // mpq_out_str(stdout, 10, rat); cout << endl;
      // rationalize_cf_from_str(PS[s], (char*)"1.41421356237309504880168872421");
      
      rationalize_cf_from_str(PS[s], kin[s]);

    } else {
      // input is rational
      if (mpq_set_str(PS[s], kin[s], 10)) {  // #rmp: change to custom function
        printf("error during conversion from user kinematic input to mpq_t rational\n");
        perror("error during conversion from user kinematic input to mpq_t rational");
        exit(1);
      }
    }

    cout << "  " << inv_sym[s] << " = "; mpq_out_str(stdout, 10, PS[s]); cout << endl;
  }
  cout << "]" << endl;

  //////
  // WRITE KIRA INPUT
  //////
  char topo_name[] = "MIeta";

  char *dir_kira = NULL;
  join_path(&dir_kira, dir_amflow, (char*)"kira/");
  make_dir(dir_kira);

  // write MIs as preferred for Kira
  char *preferred_filepath = NULL;
  join_path(&preferred_filepath, dir_kira, (char*)"preferred");
  cout << endl;
  cout << "writing eta-less MIs to " << preferred_filepath << endl;
  LI_rk1_to_file(preferred_filepath, MI, *dim, topo_name);

  // set config/ file path
  char *dir_config= NULL;
  join_path(&dir_config, dir_kira, (char*)"config/");
  make_dir(dir_config);
  char *kira_kin_filepath = NULL, *kira_int_fam_filepath = NULL;
  join_path(&kira_kin_filepath, dir_config, (char*)"kinematics.yaml");
  join_path(&kira_int_fam_filepath, dir_config, (char*)"integralfamilies.yaml");

  // config/integral_families.yaml
  cout << "writing to " << kira_int_fam_filepath << endl;
  fptr = fopen(kira_int_fam_filepath, "w");
  fprintf(fptr, "integralfamilies:\n");
  fprintf(fptr, "  - name: %s\n", topo_name);
  fprintf(fptr, "    loop_momenta: %s\n", loop_mom_nl.str);
  // fprintf(fptr, "    top_level_sectors: [%d]\n", LI_get_sector(&MI[*dim-1]));  // #review
  fprintf(fptr, "    top_level_sectors: [%d]\n", int_rk1_is_zero_to_decimal(nISP, prop_nl.nitems));  // #review
  fprintf(fptr, "    propagators:\n");
  for (int p=0; p<prop_nl.nitems; p++) {
    if (strcmp(prop_nl.item[p].item[1].str, "0") == 0) {
      if (nISP[p] == 0) {
        fprintf(fptr, "      - [%s,eta]\n", prop_nl.item[p].item[0].str);  // insert masses
      } else {
        fprintf(fptr, "      - [%s,0]\n", prop_nl.item[p].item[0].str);  // insert masses
      }
    } else {
      struct lnode nd;
      lnode_parse_expression(
        &nd,
        prop_nl.item[p].item[1].str,
        (const char **)inv_sym, kin_inv_nl.nitems, 0
      );
      char *nd_str = lnode_to_str(&nd, (char*)", ");
      mpq_t point; mpq_init(point);
      decode_tree_mpq(&point, PS, &nd_str, (char*)", ");
      free(nd_str);

      fprintf(fptr, "      - [%s,", prop_nl.item[p].item[0].str);
      mpq_out_str(fptr, 10, point);  // insert masses
      if (nISP[p] == 0) {
        fprintf(fptr, "+eta]\n");
      } else {
        fprintf(fptr, "]\n");
      }
    }    
  }
  fprintf(fptr, "    cut_propagators: []\n");
  fclose(fptr);

  // config/kinematics.yaml
  cout << "writing to " << kira_kin_filepath << endl;
  fptr = fopen(kira_kin_filepath, "w");
  fprintf(fptr, "kinematics:\n");
  fprintf(fptr, "  incoming_momenta: %s\n", ext_mom_nl.str);
  fprintf(fptr, "  outgoing_momenta: []\n");
  fprintf(fptr, "  momentum_conservation: %s\n", mom_cons_nl.str);
  fprintf(fptr, "  kinematic_invariants:\n");
  fprintf(fptr, "    - [eta, 2]\n");
  fprintf(fptr, "  scalarproduct_rules:\n");
  for (int i=0; i<sq_mom_nl.nitems; i++) {
    // process eta-less fixed kinematics
    struct lnode nd;
    lnode_parse_expression(
      &nd,
      sq_mom_nl.item[i].item[1].str,
      (const char **)inv_sym, kin_inv_nl.nitems, 0
    );
    char *nd_str = lnode_to_str(&nd, (char*)", ");
    mpq_t point; mpq_init(point);
    decode_tree_mpq(&point, PS, &nd_str, (char*)", ");
    free(nd_str);

    fprintf(fptr, "    - [[%s,%s], ",
      sq_mom_nl.item[i].item[0].str,
      sq_mom_nl.item[i].item[0].str
    );
    mpq_out_str(fptr, 10, point);
    fprintf(fptr, "]\n");
  }
  fclose(fptr);

  // write jobs for kira
  char *kira_jobs_MIs_eta_filepath = NULL, *kira_jobs_filepath = NULL;
  join_path(&kira_jobs_MIs_eta_filepath, dir_kira, (char*)"jobs_MIs-eta.yaml");
  join_path(&kira_jobs_filepath, dir_kira, (char*)"jobs.yaml");

  // jobs_MIs-eta.yaml
  cout << "writing to " << kira_jobs_MIs_eta_filepath << endl;
  fptr = fopen(kira_jobs_MIs_eta_filepath, "w");
  int r = LI_rk1_get_r(MI, *dim) + 1, s = LI_rk1_get_s(MI, *dim) + 1, d = 1;
  fprintf(fptr, "jobs:\n");
  fprintf(fptr, "  - reduce_sectors:\n");
  fprintf(fptr, "      reduce:\n");
  fprintf(fptr, "        - {r: %d, s: %d}\n", r, s+1);
  fprintf(fptr, "      select_integrals:\n");
  fprintf(fptr, "        select_mandatory_recursively:\n");
  fprintf(fptr, "          - {r: %d, s: %d, d: %d}\n", r, s, d);
  fprintf(fptr, "      preferred_masters: preferred\n");
  fprintf(fptr, "      integral_ordering: 5\n");
  fprintf(fptr, "      run_initiate: masters\n");
  fclose(fptr);

  // jobs.yaml
  cout << "writing to " << kira_jobs_filepath << endl;
  fptr = fopen(kira_jobs_filepath, "w");
  fprintf(fptr, "jobs:\n");
  fprintf(fptr, "  - reduce_sectors:\n");
  fprintf(fptr, "      reduce:\n");
  fprintf(fptr, "        - {r: %d, s: %d}\n", r, s+1);
  fprintf(fptr, "      select_integrals:\n");
  fprintf(fptr, "        select_mandatory_list:\n");
  fprintf(fptr, "          - [%s, target]\n", topo_name);
  fprintf(fptr, "      integral_ordering: 5\n");
  fprintf(fptr, "      run_initiate: true\n");
  fprintf(fptr, "      run_triangular: true\n");
  fprintf(fptr, "      run_back_substitution: true\n");
  fprintf(fptr, "  - kira2file:\n");
  fprintf(fptr, "      target:\n");
  fprintf(fptr, "        - [%s, target]\n", topo_name);
  fclose(fptr);

  //////
  // LAUNCH KIRA to find MIs with eta
  //////
  char *dir_common_eta = NULL;
  join_path(&dir_common_eta, dir_amflow, (char*)"common/");
  make_dir(dir_common_eta);
  char *MIs_eta_filepath = NULL;
  join_path(&MIs_eta_filepath, dir_common_eta, (char*)"MIs.txt");
  char *kira_MIs_eta_filepath = NULL;
  join_path(&kira_MIs_eta_filepath, dir_kira, (char*)"results/");
  join_path(&kira_MIs_eta_filepath, kira_MIs_eta_filepath, topo_name);
  join_path(&kira_MIs_eta_filepath, kira_MIs_eta_filepath, (char*)"/masters");

  char *command = (char*) malloc(MAX_VALUE_LEN*sizeof(char));

  int execute;
  if (redo == 0) {
    execute = 0;
    cout << "skip call to Kira" << endl;
  } else if (redo == 1) {
    execute = 1;
  } else if (redo == -1) {
    if (file_exists(MIs_eta_filepath)) {
      execute = 0;
      cout << "Kira output already exists" << endl;
    } else {
      execute = 1;
    }
  }

  if (execute) {
    fprintf(terminal, "call Kira to find MIs... "); fflush(terminal); usleep(sleep_time);
    snprintf(
      command, MAX_VALUE_LEN*sizeof(char),
      "cd %s && rm -rf results/ sectormappings/ tmp/; cd -", dir_kira
    );
    cout << endl; cout << "executing shell command: " << endl; cout << "  $ " << command << endl;
    cout << endl;
    system(command);
    snprintf(
      command, MAX_VALUE_LEN*sizeof(char),
      "cd %s && kira --parallel=%d jobs_MIs-eta.yaml; cd -", dir_kira, opt_kira_parallel
    );
    cout << endl; cout << "executing shell command: " << endl; cout << "  $ " << command << endl;
    cout << endl;
    system(command);
    snprintf(
      command, MAX_VALUE_LEN*sizeof(char),
      "cp %s %s", kira_MIs_eta_filepath, MIs_eta_filepath
    );
    cout << endl; cout << "executing shell command: " << endl; cout << "  $ " << command << endl;
    cout << endl;
    system(command);
    fprintf(terminal, "\033[25D\033[K"); fflush(terminal); usleep(sleep_time);
  }

  //////
  // READ MIs with eta
  //////
  cout << "reading MIs with eta from " << MIs_eta_filepath << endl;

  LI *MI_eta_tmp = NULL;
  LI_rk1_from_file(MIs_eta_filepath, &MI_eta_tmp, dim_eta, topo_name);
  // cout << "unsorted:" << endl;
  // LI_rk1_pows_print(MI_eta_tmp, *dim_eta);

  // SORT MIs
  if (*MI_eta) {
    delete[] *MI_eta;
  }
  *MI_eta = NULL;
  LI_rk1_qsort(MI_eta, MI_eta_tmp, *dim_eta);
  cout << endl; cout << "MIs with auxiliary mass" << endl;
  LI_rk1_pows_print(*MI_eta, *dim_eta);

  // STORE INDICES OF eta-LESS MIs
  *MI_idx = new int[*dim];
  LI_rk1_get_idx_rk1(*MI_idx, MI, *dim, *MI_eta, *dim_eta);

  //////
  // GENERATE BOUNDARY FILES
  //////
  char *bound_behav_filepath = NULL;
  join_path(&bound_behav_filepath, dir_common_eta, (char*)"bound_behav.txt");

  cout << endl;
  cout << "writing to " << bound_behav_filepath << endl;
  fptr = fopen(bound_behav_filepath, "w");
  fprintf(fptr, "num-eig: 1\n");
  for (int m=0; m<*dim_eta; m++) {
    fprintf(fptr, "MI%d: [[[eig, 0], [pow, %d]]]\n", m, LI_get_r(&(*MI_eta)[m])-LI_get_s(&(*MI_eta)[m]) - 2*loop_mom_nl.nitems);
  }
  fclose(fptr);

  // compute infinite mass limit structure
  char *bound_build_filepath = NULL;
  join_path(&bound_build_filepath, dir_common_eta, (char*)"bound_build.txt");

  cout << "writing to " << bound_build_filepath << endl;
  fptr = fopen(bound_build_filepath, "w");
  int ***loop_curr, curr_dim = int_pow_int(2, loop_mom_nl.nitems), factorized;
  malloc_rk3_tens(loop_curr, *dim_eta, 2, curr_dim);
  for (int m=0; m<*dim_eta; m++) {
    // cout << "m = " << m << endl;
    LI_get_loop_current(loop_curr[m], &(*MI_eta)[m], prop, nISP);

    // // check for loop factorization
    // factorized = 1;
    // for (int l=1; l<loop_mom_nl.nitems; l++) {
    //   for (int c=int_pow_int(2, l)+1; c<int_pow_int(2, l+1); c++) {
    //     if (loop_curr[m][c] > 0) {
    //       factorized = 0;
    //       break;
    //     }
    //   }
    // }

    if (loop_mom_nl.nitems == 1) {
      // tadpole with dots
      fprintf(
        fptr,
        "MI%d: [[[[loop, 1], [id, -1], [sym, eta], [point, fin], [intargs, [%d]]]]]\n",
        m, LI_get_r(&(*MI_eta)[m])-LI_get_s(&(*MI_eta)[m])
      );
    } else if (loop_mom_nl.nitems == 2) {
      fprintf(
        fptr,
        "MI%d: [[[[loop, 2], [id, 7], [sym, eta], [point, fin], [intargs, [%d, %d, %d, %d, %d, %d]]]]]\n",
        m, loop_curr[m][0][1], loop_curr[m][0][2], loop_curr[m][0][3], loop_curr[m][1][1], loop_curr[m][1][2], loop_curr[m][1][3]
      );
    }
  }
  fclose(fptr);

  //////
  // WRITE TARGETS FOR DEs
  //////
  char *target_filepath = NULL;
  join_path(&target_filepath, dir_kira, (char*)"target");
  cout << "writing targets to " << target_filepath << endl;
  fptr = fopen(target_filepath, "w");
  char *tmp_str = (char*) malloc(((MAX_POW_DIGITS+1)*(*dim_eta)+2+1)*sizeof(char));
  char *out_str = NULL;
  int *der_pow = new int[prop_nl.nitems], m = 0;
  for (m=0; m<*dim_eta; m++) {
    // cout << "m = " << m << endl;
    int_rk1_to_str(&tmp_str, (*MI_eta)[m].nprop, (*MI_eta)[m].pows);
    // cout << "power list member: " << (*MI_eta)[m].pows_str << endl;
    // cout << "power list: " << tmp_str << endl;

    (*MI_eta)[m].nmass = 0;
    for (int p=0; p<(*MI_eta)[m].nprop; p++) {
      der_pow[p] = (*MI_eta)[m].pows[p];
      if ((*MI_eta)[m].pows[p] != 0 && nISP[p] == 0) {
        (*MI_eta)[m].nmass++;
      }
    }
    LI_mder_alloc(&(*MI_eta)[m]);
    
    (*MI_eta)[m].nmass = 0;
    for (int p=0; p<(*MI_eta)[m].nprop; p++) {
      if (nISP[p] == 1) {continue;}
      if ((*MI_eta)[m].pows[p] != 0) {
        // cout << "derivative w.r.t. propagator " << p << ":" << endl;
        (*MI_eta)[m].mder_idx[(*MI_eta)[m].nmass] = p;  // store idx
        der_pow[p]++;
        // for (int p=0; p<(*MI_eta)[m].nprop; p++) {
        //   cout << der_pow[p] << ", ";
        // }
        // cout << endl;

        // convert to string
        int_rk1_to_str(&tmp_str, (*MI_eta)[m].nprop, der_pow);
        join_path(&tmp_str, topo_name, tmp_str);
        // cout << "string: " << tmp_str << endl;
        fprintf(fptr, "%s\n", tmp_str);
        // store
        LI_init(
          &(*MI_eta)[m].mder_contr[(*MI_eta)[m].nmass],
          (*MI_eta)[m].nprop, der_pow
        );
        // cout << "converted power list: ";
        // LI_pows_print(&(*MI_eta)[m].mder_contr[(*MI_eta)[m].nmass]); cout << endl;
        (*MI_eta)[m].nmass++;
        
        // reset
        der_pow[p]--;
      }
    }

  }
  fclose(fptr);

  //////
  // LAUNCH KIRA to reduce targets
  //////
  char *results_filepath = NULL;
  join_path(&results_filepath, dir_kira, (char*)"results/");
  join_path(&results_filepath, results_filepath, topo_name);
  join_path(&results_filepath, results_filepath, (char*)"/kira_target.kira");
  
  if (redo == 0) {
    execute = 0;
    cout << "skip call to Kira" << endl;
  } else if (redo == 1) {
    execute = 1;
  } else if (redo == -1) {
    if (file_exists(results_filepath)) {
      execute = 0;
      cout << "Kira output already exists" << endl;
    } else {
      execute = 1;
    }
  }

  if (execute) {
    fprintf(terminal, "call Kira to reduce derivatives... "); fflush(terminal); usleep(sleep_time);    
    snprintf(
      command, MAX_VALUE_LEN*sizeof(char),
      "cd %s && kira --parallel=%d jobs.yaml; cd -", dir_kira, opt_kira_parallel
    );
    cout << endl; cout << "executing shell command: " << endl; cout << "  $ " << command << endl;
    cout << endl;
    system(command);
    fprintf(terminal, "\033[45D\033[K"); fflush(terminal); usleep(sleep_time);
  }

  return;

  //////
  // ASSEMBLE DEs FROM KIRA RESULTS
  //////
  cout << "ASSEMBLE DEs FROM KIRA RESULTS" << endl;
  char ***mat;
  malloc_rk2_tens(mat, *dim_eta, *dim_eta);

  cout << "extracting DE from " << results_filepath << endl;
  kira_to_DE_str(mat, *MI_eta, *dim_eta, results_filepath, topo_name);

  // WRITE DE TO FILE
  char *DEs_filepath = NULL;
  join_path(&DEs_filepath, dir_common_eta, (char*)"0.txt");

  cout << "writing DEs to " << DEs_filepath << endl;
  fptr = fopen(DEs_filepath, "w");
  for (m=0; m<*dim_eta; m++) {
    for (int n=0; n<*dim_eta-1; n++) {
      fprintf(fptr, "%s\t", mat[m][n]);
    }
    fprintf(fptr, "%s\n", mat[m][*dim_eta-1]);
  }
  fclose(fptr);

}

