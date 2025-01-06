#ifndef BC_GEN_H
#define BC_GEN_H

void mpc_parg(mpfr_t *out, mpc_t *in);

void mpc_plog(mpc_t *out, mpc_t *in);

void formulary(mpc_t *out, int loop, int id, mpc_t eps, mpc_t *scale, int *intargs);

void dispatcher(
  // OUTPUT
  mpc_t *bound,
  // INPUT
  // char *id_str, mpc_t eps,  int *is_mass, 
  char *eps_str, char ***kin, char **symbols, int ninvs,
  int dim, int *nadd, int **nfac, int ***loop, int ***id, char ****sym, char ****point, int ****intargs
);

void generate_boundaries(
  // OUTPUT
  mpc_t **bound,
  // INPUT
  char *filepath, int eps_num, int dim, int loops,
  char **eps_str, char ***kin, char **symbols, int ninvs
);

#endif
