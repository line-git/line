#ifndef TOPOLOGY_H
#define TOPOLOGY_H

#include "setup.h"

typedef struct propagator {
  int nsym;  // global number of momenta
  int nloop;  // global number of masses
  struct lnode *momentum_tree;  // momentum math tree
  struct lnode *mass_tree;  // mass math tree
  int is_massive;  // 0: massless, 1: massive
  int *loop_current;  // 110: has k1 and k2, 101: has k1 and k3
  int loop_current_dec;  // (as decimal) 110: 3, 101: 5
} propagator;

void propagator_build(
  propagator *prop
);

void propagator_set_nlist(
  // OUTPUT
  propagator *prop,
  // INPUT
  nlist *nl, char **symbols, int nsym, int nloop
);

void propagator_copy(
  // OUTPUT
  propagator *out,
  // INPUT
  propagator *in
);

void propagator_print(
  propagator *prop, char **symbols
);

typedef struct LI {
  int nprop;  // global number of propagators
  int *pows;  // power list
  char *pows_str;  // power list as string
  int nmass;  // number of massive propagators
  struct LI *mder_contr;  // contributions to massive derivative
  int *mder_idx;  // propagator indices for contributions to massive derivative
} LI;

void LI_pows_print(
  LI *li
);

void LI_rk1_pows_print(
  LI *li, int dim
);

void LI_build(
  // IN-OUT
  LI *li
);

void LI_rk1_build(
  // IN-OUT
  LI *li, int dim
);

void LI_pows_str_generate(
  // IN-OUT
  LI *li
);

void LI_pows_alloc(
  // IN-OUT
  LI *li
);

void LI_mder_alloc(
  // IN-OUT
  LI *li
);

void LI_init(
  // IN-OUT
  LI *li,
  // INPUT
  int nprop, int *pows
);

int pows_str_to_MI_idx(
  // INPUT
  char *pows_str, int dim, LI *MIs
);

void LI_set(
  // OUTPUT
  LI *liout,
  // INPUT
  LI *liin
);

int LI_get_sector(
  LI *li
);

int LI_get_r(
  LI *li
);

int LI_get_s(
  LI *li
);

void LI_get_loop_current(
  // OUTPUT
  int **loop_current,
  // INPUT
  LI* li, propagator *prop, int *nISP
);

int LI_cmp(
  LI *li1, LI *li2
);

void LI_rk1_qsort_labels(
  // IN-OUT
  int *label,
  // INPUT
  int first, int last, LI *li
);

void LI_rk1_qsort(
  // OUTPUT
  LI **liout,
  // INPUT
  LI *liin, int dim
);

void LI_rk1_get_nISP(
  // OUTPUT
  int *nISP,
  // INPUT
  LI *li, int dim
);

int LI_rk1_get_idx(
  char *pows_str,
  LI *li, int dim
);

void LI_rk1_get_idx_rk1(
  // OUTPUT
  int *idx,
  // INPUT
  LI *li1, int dim1,
  LI *li2, int dim2
);

int LI_rk1_get_r(
  LI *li, int dim
);

int LI_rk1_get_s(
  LI *li, int dim
);

void LI_rk1_from_file(
  char *filepath,
  // OUTPUT
  LI **li, int *dim,
  // INPUT
  char *topo_name
);

void LI_rk1_to_file(
  char *filepath,
  LI *li, int dim,
  char *topo_name
);

#endif