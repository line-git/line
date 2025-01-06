#include <iostream>
#include <cstring>
// #include <complex>
#include "mpc.h"

using namespace std;

#include "setup.h"
#include "utils.h"
#include "tensor_utils.h"
#include "topology.h"
extern "C" {
  #include "ex_tree.h"
}


void propagator_build(
  propagator *prop
) {
  prop->momentum_tree = (struct lnode*) malloc(sizeof(struct lnode));
  prop->mass_tree = (struct lnode*) malloc(sizeof(struct lnode));
  prop->nsym = -1;
  prop->nloop = -1;
  prop->is_massive = -1;
  prop->loop_current_dec = -1;
  prop->loop_current = NULL;
}


void propagator_set_nlist(
  // OUTPUT
  propagator *prop,
  // INPUT
  nlist *nl, char **symbols, int nsym, int nloop
) {
  // global info
  prop->nsym = nsym;
  prop->nloop = nloop;

  // momentum tree
  char *tmp_str = (char*) malloc((strlen(nl->str)+1)*sizeof(char));  // enough for both components
  strcpy(tmp_str, nl->item[0].str);  // backup copy
  lnode_parse_expression(
    prop->momentum_tree,
    tmp_str, (const char**) symbols, nsym, 0
  );

  // mass tree
  strcpy(tmp_str, nl->item[1].str);  // backup copy
  lnode_parse_expression(
    prop->mass_tree,
    tmp_str, (const char**) symbols, nsym, 0
  );

  // detect mass
  if (strcmp(nl->item[1].str, "0") == 0) {
    prop->is_massive = 0;
  } else {
    prop->is_massive = 1;
  }

  // loop current
  // int *loop_current_mirror = new int[nsym];
  prop->loop_current = new int[nsym];
  lnode_contains_sym(prop->loop_current, nsym, prop->momentum_tree);
  // for (int s=0; s<nloop; s++) {
  //   prop->loop_current[s] = loop_current_mirror[nloop-1-s];
  // }

  // cout << "loop current: ";
  // for (int l=0; l<nloop; l++) { cout << prop->loop_current[l] << ", ";} cout << endl;

  // decimal loop current
  prop->loop_current_dec = int_rk1_is_positive_to_decimal(prop->loop_current, nloop);

}


void propagator_copy(
  // OUTPUT
  propagator *out,
  // INPUT
  propagator *in
) {
  out->nsym = in->nsym;
  out->nloop = in->nloop;
  out->is_massive = in->is_massive;

	// COPY MOMENTUM TREE
	if (in->momentum_tree) {
		if (out->momentum_tree) {
			free(out->momentum_tree);
		}
		out->momentum_tree = (struct lnode*) malloc(sizeof(struct lnode));
		lnode_copy(&out->momentum_tree[0], &in->momentum_tree[0]);
	} else {
		out->momentum_tree = NULL;
	}

	// COPY MASS TREE
	if (in->mass_tree) {
		if (out->mass_tree) {
			free(out->mass_tree);
		}
		out->mass_tree = (struct lnode*) malloc(sizeof(struct lnode));
		lnode_copy(&out->mass_tree[0], &in->mass_tree[0]);
	} else {
		out->mass_tree = NULL;
	}

  // COPY LOOP CURRENT
  if (in->loop_current) {
    if (out->loop_current) {
      delete[] out->loop_current;
    }
    out->loop_current = new int[in->nloop];
    for (int l=0; l<in->nloop; l++) {
      out->loop_current[l] = in->loop_current[l];
    }
  } else {
    out->loop_current = NULL;
  }

}


void propagator_print(
  propagator *prop, char **symbols
) {
  printf("[mom: ");
  lnode_print(prop->momentum_tree, (const char**)symbols);
  printf(", mass: ");
  lnode_print(prop->mass_tree, (const char**)symbols);
  printf(", lcurr: %d]", prop->loop_current_dec);
  return;
}


void LI_pows_print(
  LI *li
) {
  printf("[");
  for (int p=0; p<li->nprop; p++) {
    printf("%d, ", li->pows[p]);
  }
  printf("]");
}


void LI_rk1_pows_print(
  LI *li, int dim
) {
  for (int i=0; i<dim; i++) {
    printf("i = %d: ", i);
    LI_pows_print(&li[i]);
    printf("\n");
  }
}


void LI_build(
  // IN-OUT
  LI *li
) {
  li->pows = NULL;
  li->pows_str = NULL;
  li->mder_contr = NULL;
  li->mder_idx = NULL;
  return;
}


void LI_rk1_build(
  // IN-OUT
  LI *li, int dim
) {
  for (int i=0; i<dim; i++) {
    LI_build(&li[i]);
  }
  return;
}


void LI_pows_str_generate(
  // IN-OUT
  LI *li
) {
  int_rk1_to_str(&li->pows_str, li->nprop, li->pows);  
}


void LI_pows_alloc(
  // IN-OUT
  LI *li
) {
  li->pows = new int[li->nprop];
}


void LI_mder_alloc(
  // IN-OUT
  LI *li
) {
  li->mder_contr = new LI[li->nmass];
  li->mder_idx = new int[li->nmass];
  for (int i=0; i<li->nmass; i++) {
    LI_build(&li->mder_contr[i]);
  }
}


void LI_init(
  // IN-OUT
  LI *li,
  // INPUT
  int nprop, int *pows
) {
  // set nprop
  li->nprop = nprop;

  // alloc power list
  if (li->pows) {
    delete[] li-> pows;
  }
  li->pows = new int[nprop];

  // set power list
  for (int p=0; p<nprop; p++) {
    li->pows[p] = pows[p];
  }

  // convert to string
  LI_pows_str_generate(li);

}


int pows_str_to_MI_idx(
  // INPUT
  char *pows_str, int dim, LI *MIs
) {
  int m;
  for (m=0; m<dim; m++) {
    if (strcmp(pows_str, MIs[m].pows_str) == 0) {
      break;
    }
  }

  if (m == dim) {
    return -1;
  } else {
    return m;
  }
}


void LI_set(
  // OUTPUT
  LI *liout,
  // INPUT
  LI *liin
) {
  // power list
  liout->nprop = liin->nprop;
  if (liout->pows) {
    delete[] liout-> pows;
  }
  liout->pows = new int[liin->nprop];
  for (int p=0; p<liin->nprop; p++) {
    liout->pows[p] = liin->pows[p];
    // cout << "liout->pows[" << p << "] = " << liout->pows[p] << endl;
  }

  // power list string
  if (liout->pows_str) {
    free(liout->pows_str);
  }
  liout->pows_str = (char*) malloc((strlen(liin->pows_str)+1)*sizeof(char));
  strcpy(liout->pows_str, liin->pows_str);

  // number of massive propagators
  liout->nmass = liin->nmass;

  // massive derivative contributions
  if (liout->mder_idx) {
    delete[] liout->mder_idx;
  }
  if (liin->mder_idx) {
    liout->mder_idx = new int[liin->nmass];
    for (int k=0; k<liin->nmass; k++) {
      liout->mder_idx[k] = liin->mder_idx[k];
    }
  } else {
    liout->mder_idx = NULL;
  }
  
  if (liout->mder_contr) {
    delete[] liout->mder_contr;
  }
  if (liin->mder_contr) {
    liout->mder_contr = new LI[liin->nmass];
    LI_rk1_build(liout->mder_contr, liin->nmass);
    for (int k=0; k<liin->nmass; k++) {
      LI_set(&liout->mder_contr[k], &liin->mder_contr[k]);
    }
  } else {
    liout->mder_contr = NULL;
  }

  return;
}


int LI_get_sector(
  LI *li
) {
  return int_rk1_is_positive_to_decimal(li->pows, li->nprop);
}


int LI_get_r(
  LI *li
) {
  return int_rk1_sum_postivie(li->pows, li->nprop);
}


int LI_get_s(
  LI *li
) {
  return int_rk1_sum_negative_abs(li->pows, li->nprop);
}


void LI_get_loop_current(
  // OUTPUT
  int **loop_current,
  // INPUT
  LI* li, propagator *prop, int *nISP
) {
  // initialize to zero
  for (int p=0; p<int_pow_int(2, prop[0].nloop); p++) {
    loop_current[0][p] = 0;
    loop_current[1][p] = 0;
  }

  for (int p=0; p<li->nprop; p++) {
    loop_current[nISP[p]][prop[p].loop_current_dec] += li->pows[p];
    // cout << "p = " << p << endl;
    // cout << "loop curr idx = " << prop[p].loop_current_dec << endl;
    // cout << "pow = " << li->pows[p] << endl;
    // cout << "loop curr cum = " << loop_current[prop[p].loop_current_dec] << endl;
  }
}


int LI_cmp(
  LI *li1, LI *li2
) {
  // cout << "comparing LIS" << endl;
  // cout << "li1: " << endl;
  // LI_pows_print(li1); cout << endl;
  // cout << "li2: " << endl;
  // LI_pows_print(li2); cout << endl;  
  int positive1 = int_rk1_count_postivie(li1->pows, li1->nprop);
  int positive2 = int_rk1_count_postivie(li2->pows, li2->nprop);
  // cout << "positive1 = " << positive1 << endl;
  // cout << "positive2 = " << positive2 << endl;
  if (positive1 < positive2) {
    return -1;
  } else if (positive1 > positive2) {
    return 1;
  } else if (positive1 == positive2) {
    int sec1 = LI_get_sector(li1);
    int sec2 = LI_get_sector(li2);
    // cout << "sec1 = " << sec1 << endl;
    // cout << "sec2 = " << sec2 << endl;
    if (sec1 < sec2) {
      return -1;
    } else if (sec1 > sec2) {
      return 1;
    } else if (sec1 == sec2) {
      return 0;
    } else {
      perror("error while comparing LIs sectors");
      exit(1);
    }
  } else {
    perror("error while comparing LIs");
    exit(1);
  }
}


void LI_rk1_qsort_labels(
  // IN-OUT
  int *label,
  // INPUT
  int first, int last, LI *li
) {
  int i, j, pivot, temp;
  int valuei, valuej, valuepivot;
  if(first<last){
    pivot = first;
    i = first;
    j = last;
    while(i<j){
      while(LI_cmp(&li[label[i]], &li[label[pivot]]) <= 0 && i<last) {
        i++;
      }
      while(LI_cmp(&li[label[j]], &li[label[pivot]]) > 0) {
        j--;
      }
      if(i<j){
        temp=label[i];
        label[i]=label[j];
        label[j]=temp;

      }
    }
    temp = label[pivot];
    label[pivot] = label[j];
    label[j] = temp;
    LI_rk1_qsort_labels(label, first, j-1, li);
    LI_rk1_qsort_labels(label, j+1, last, li);
  }
}


void LI_rk1_qsort(
  // OUTPUT
  LI **liout,
  // INPUT
  LI *liin, int dim
) {
  if (!*liout) {
    *liout = new LI[dim];
    LI_rk1_build(*liout, dim);
  }
  
  // starting order: 0, 1, 2, ...
  int *label = new int[dim];
  for (int k=0; k<dim; k++) {
    label[k] = k;
  }

  // sort labels
  LI_rk1_qsort_labels(label, 0, dim-1, liin);

  // cout << "sorted labels:" << endl;
  // for (int k=0; k<dim; k++) {
  //   cout << label[k] << ", ";
  // }
  // cout << endl;

  // write output in the right order
  for (int k=0; k<dim; k++) {
    // cout << "set " << k << " to " << label[k] << endl;
    // cout << "liin: "; LI_pows_print(&liin[label[k]]); cout << endl;
    LI_set(&(*liout)[k], &liin[label[k]]);
  }

}


void LI_rk1_get_nISP(
  // OUTPUT
  int *nISP,
  // INPUT
  LI *li, int dim
) {
  int nprop = li[0].nprop;
  for (int p=0; p<nprop; p++) {
    nISP[p] = 1;
    for (int m=0; m<dim; m++) {
      if (li[m].pows[p] > 0) {
        nISP[p] = 0;
        break;
      }
    }
  }
}


int LI_rk1_get_idx(
  char *pows_str,
  LI *li, int dim
) {
  int mi_idx = 0;
  for (mi_idx=0; mi_idx<dim; mi_idx++) {
    // cout << "candidate : " << li[mi_idx].pows_str << endl;
    if (strcmp(pows_str, li[mi_idx].pows_str) == 0) {
      // cout << "found" << endl;
      return mi_idx;
    }
  }

  return -1;
}


void LI_rk1_get_idx_rk1(
  // OUTPUT
  int *idx,
  // INPUT
  LI *li1, int dim1,
  LI *li2, int dim2
) {
  for (int m=0; m<dim1; m++) {
    idx[m] = LI_rk1_get_idx(li1[m].pows_str, li2, dim2);
    // cout << "m, pows_str, idx = " << m << ", " << li1[m].pows_str << ", " << idx[m] << endl;
  }
}


int LI_rk1_get_r(
  LI *li, int dim
) {
  int r = LI_get_r(&li[0]);
  int rp;
  for (int m=1; m<dim; m++) {
    rp = LI_get_r(&li[m]);
    if (rp > r) {
      r = rp;
    }
  }
  return r;
}


int LI_rk1_get_s(
  LI *li, int dim
) {
  int s = LI_get_s(&li[0]);
  int sp;
  for (int m=1; m<dim; m++) {
    sp = LI_get_s(&li[m]);
    if (sp > s) {
      s = sp;
    }
  }

  return s;
}


void LI_rk1_from_file(
  char *filepath,
  // OUTPUT
  LI **li, int *dim,
  // INPUT
  char *topo_name
) {
  FILE *fptr = fopen(filepath, "r");
  if (fptr == NULL) {
    printf("Could not open file %s\n", filepath);
    exit(1);
  }
  
  unsigned long int line_size = MAX_VALUE_LEN * sizeof(char);
  
  // COUNT MIs
  *dim = count_lines_fptr(fptr);  // #2BD: count only non-empty lines
  
  if (*li) {
    delete[] *li;
  }
  *li = new LI[*dim];

  // READ POWER LISTS
  fseek(fptr, 0, SEEK_SET);
  int nprop, c = 0, m = 0, topo_len = strlen(topo_name);
  char **pows_str = NULL;
  char *line = (char*) malloc(line_size);
  char *tmp_str = (char*) malloc(((MAX_POW_DIGITS+1)*(*dim)+2)*sizeof(char));
	while (fgets(line, MAX_VALUE_LEN, fptr)) {
    // // skip in-line comment
    // end_at_char(line, '#');

    // // go to left delimiter
    // while(line[c] != '[') {
    //   c++;
    // }

    // skip topology name
    c = topo_len;
  
    // read power list
    int cp = 0;
    while (line[c-1] != ']') {
      if (!isspace(line[c])) {
        tmp_str[cp++] = line[c];
      }
      c++;
    }
    tmp_str[cp] = '\0';
    c++;
    // cout << "power list: " << tmp_str << endl;

    // parse list
    parse_list(&pows_str, &nprop, tmp_str);
    LI_build(&(*li)[m]);
    (*li)[m].nprop = nprop;
    (*li)[m].pows = new int[nprop];
    // cout << "parsed power list: ";
    for (int p=0; p<nprop; p++) {
      // cout << pows_str[p] << ", ";
      (*li)[m].pows[p] = atoi(pows_str[p]);
    }
    // cout << endl;
    (*li)[m].pows_str = (char*) malloc((strlen(tmp_str)+1)*sizeof(char));
    strcpy((*li)[m].pows_str, tmp_str);
    // cout << "copied power list: " << (*li)[m].pows_str << endl;

    m++;
  }
  fclose(fptr);

}


void LI_rk1_to_file(
  char *filepath,
  LI *li, int dim,
  char *topo_name
) {
  FILE *fptr = fopen(filepath, "w");
  for (int m=0; m<dim; m++) {
    fprintf(fptr, "%s%s\n", topo_name, li[m].pows_str);
  }
}

