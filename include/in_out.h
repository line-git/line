#ifndef IN_OUT_H
#define IN_OUT_H

void int_rk1_to_file(
  char *file_name, int *tens, int dim1
);

void int_rk1_from_file(
  char *file_name, int *tens, int dim1
);

void int_rk2_to_file(
  char *file_name, int **tens, int dim1, int dim2
);

void int_rk2_from_file(
  char *file_name, int **tens, int dim1, int dim2
);

void mpc_rk1_to_file(
  char *file_name, mpc_t *tens, int dim
);

void mpc_rk2_to_file(
  char *file_name, mpc_t **tens, int dim1, int dim2
);

void int_rk0_int_rk1_to_file(
  char *file_name, int *tens, int dim, int scalar
);

void int_rk0_mpc_rk1_to_file(
  char *file_name, mpc_t *tens, int dim, int scalar
);

void mpc_rk1_from_file(
  char *file_name, mpc_t *tens, int dim
);

void mpc_rk2_from_file(
  char *file_name, mpc_t **tens, int dim1, int dim2
);

void int_rk0_int_rk1_from_file(
  char *file_name, int *tens, int dim, int *scalar
);

void int_rk0_mpc_rk1_from_file(
  char *file_name, mpc_t *tens, int dim, int *scalar
);

int copy_file(
  const char *source_path, const char *destination_path
);

void print_result(
  FILE *resfptr, int precision,
  mpc_t **res, int dim, int order, int nloops
);

#endif