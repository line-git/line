#include <stdio.h>
#include <stdlib.h>
#include "mpc.h"
#include "global_vars.h"


void int_rk1_to_file(
  char *file_name, int *tens, int dim1
) {
  FILE *fptr = fopen(file_name, "w");
  for (int i1=0; i1<dim1; i1++) {
    fprintf(fptr, "%d\n", tens[i1]);
  }
  fclose(fptr);
}


void int_rk1_from_file(
  char *file_name, int *tens, int dim1
) {
  FILE *fptr = fopen(file_name, "r");
  
  // check if file exists
  if (fptr == NULL) {
    printf("Could not open file %s\n", file_name);
    exit(1);
  }

  int count;
  long int save_pos;
  char c, *line;
  for (int i1=0; i1<dim1; i1++) {
    fscanf(fptr, "%d\n", &tens[i1]);
  }

  fclose(fptr);
}


void int_rk2_to_file(
  char *file_name, int **tens, int dim1, int dim2
) {
  FILE *fptr = fopen(file_name, "w");
  for (int i1=0; i1<dim1; i1++) {
    for (int i2=0; i2<dim2; i2++) {
      fprintf(fptr, "%d", tens[i1][i2]);
      fprintf(fptr, "\n");
    }
  }
  fclose(fptr);
}


void int_rk2_from_file(
  char *file_name, int **tens, int dim1, int dim2
) {
  FILE *fptr = fopen(file_name, "r");
  
  // check if file exists
  if (fptr == NULL) {
    printf("Could not open file %s\n", file_name);
    exit(1);
  }

  int count;
  long int save_pos;
  char c, *line;
  for (int i1=0; i1<dim1; i1++) {
    for (int i2=0; i2<dim2; i2++) {
      fscanf(fptr, "%d\n", &tens[i1][i2]);
    }
  }

  fclose(fptr);
}


void mpc_rk1_to_file(
  char *file_name, mpc_t *tens, int dim
) {
  FILE *fptr = fopen(file_name, "w");
  for (int i=0; i<dim; i++) {
    mpc_out_str(fptr, 10, 0, tens[i], MPFR_RNDN);
    fprintf(fptr, "\n");
  }
  fclose(fptr);
}


void mpc_rk2_to_file(
  char *file_name, mpc_t **tens, int dim1, int dim2
) {
  FILE *fptr = fopen(file_name, "w");
  for (int i1=0; i1<dim1; i1++) {
    for (int i2=0; i2<dim2; i2++) {
      mpc_out_str(fptr, 10, 0, tens[i1][i2], MPFR_RNDN);
      fprintf(fptr, "\n");
    }
  }
  fclose(fptr);
}


void int_rk0_int_rk1_to_file(
  char *file_name, int *tens, int dim, int scalar
) {
  FILE *fptr = fopen(file_name, "w");
  fprintf(fptr, "%d\n", scalar);
  for (int i=0; i<dim; i++) {
    fprintf(fptr, "%d\n", tens[i]);
    fprintf(fptr, "\n");
  }
  fclose(fptr);
}


void int_rk0_mpc_rk1_to_file(
  char *file_name, mpc_t *tens, int dim, int scalar
) {
  FILE *fptr = fopen(file_name, "w");
  fprintf(fptr, "%d\n", scalar);
  for (int i=0; i<dim; i++) {
    mpc_out_str(fptr, 10, 0, tens[i], MPFR_RNDN);
    fprintf(fptr, "\n");
  }
  fclose(fptr);
}


void mpc_rk1_from_file(
  char *file_name, mpc_t *tens, int dim
) {
  FILE *fptr = fopen(file_name, "r");
  
  // check if file exists
  if (fptr == NULL) {
    printf("Could not open file %s\n", file_name);
    exit(1);
  }

  int count;
  long int save_pos;
  char c, *line;
  for (int i=0; i<dim; i++) {
    // // count characters in the line
    // save_pos = ftell(fptr);  // save position
    // count = 0;
    // for (c=getc(fptr); c!='\n'; c=getc(fptr)) {
    //   count++;
    // }
    // fseek(fptr, save_pos, SEEK_SET);  // restore position

    // // get line
    // line = (char*) malloc(count*sizeof(char));
    // fgets(line, count, fptr);

    // convert to mpc
    mpc_inp_str(tens[i], fptr, 0, 10, MPFR_RNDN);
    
  }

  fclose(fptr);
}


void mpc_rk2_from_file(
  char *file_name, mpc_t **tens, int dim1, int dim2
) {
  FILE *fptr = fopen(file_name, "r");
  
  // check if file exists
  if (fptr == NULL) {
    printf("Could not open file %s\n", file_name);
    exit(1);
  }

  int count;
  long int save_pos;
  char c, *line;
  for (int i1=0; i1<dim1; i1++) {
    for (int i2=0; i2<dim2; i2++) {
      mpc_inp_str(tens[i1][i2], fptr, 0, 10, MPFR_RNDN);
    }
  }

  fclose(fptr);
}


void int_rk0_int_rk1_from_file(
  char *file_name, int *tens, int dim, int *scalar
) {
  FILE *fptr = fopen(file_name, "r");
  
  // check if file exists
  if (fptr == NULL) {
    printf("Could not open file %s\n", file_name);
    exit(1);
  }

  fscanf(fptr, "%d\n", scalar);

  int count;
  long int save_pos;
  char c, *line;
  for (int i=0; i<dim; i++) {
    fscanf(fptr, "%d\n", &tens[i]);
  }

  fclose(fptr);
}


void int_rk0_mpc_rk1_from_file(
  char *file_name, mpc_t *tens, int dim, int *scalar
) {
  FILE *fptr = fopen(file_name, "r");
  
  // check if file exists
  if (fptr == NULL) {
    printf("Could not open file %s\n", file_name);
    exit(1);
  }

  fscanf(fptr, "%d\n", scalar);

  int count;
  long int save_pos;
  char c, *line;
  for (int i=0; i<dim; i++) {
    // // count characters in the line
    // save_pos = ftell(fptr);  // save position
    // count = 0;
    // for (c=getc(fptr); c!='\n'; c=getc(fptr)) {
    //   count++;
    // }
    // fseek(fptr, save_pos, SEEK_SET);  // restore position

    // // get line
    // line = (char*) malloc(count*sizeof(char));
    // fgets(line, count, fptr);

    // convert to mpc
    mpc_inp_str(tens[i], fptr, 0, 10, MPFR_RNDN);
  }

  fclose(fptr);
}


int copy_file(
  const char *source_path, const char *destination_path
) {
	FILE *source_file, *destination_file;
	char buffer[BUFSIZ];
	size_t bytes_read;

	// Open source file for reading
	source_file = fopen(source_path, "rb");
	if (source_file == NULL) {
		perror("Error opening source file");
		return 1;
	}

	// Open destination file for writing
	destination_file = fopen(destination_path, "wb");
	if (destination_file == NULL) {
		perror("Error opening destination file");
		fclose(source_file);
		return 1;
	}

	// Copy contents from source file to destination file
	while ((bytes_read = fread(buffer, 1, sizeof(buffer), source_file)) > 0) {
		fwrite(buffer, 1, bytes_read, destination_file);
	}

	// Close files
	fclose(source_file);
	fclose(destination_file);

	printf("File copied successfully.\n");

	return 0;
}


void print_result(
  FILE *resfptr, int precision,
  mpc_t **res, int dim, int order, int nloops
) {
  char res_format[20];
  snprintf(res_format, sizeof(res_format), "%%.%dRe", precision-1);
  int sign;
  for (int i=0; i<dim; i++) {
    fprintf(resfptr, "MI n. %d\n", i);
    for (int ep=0; ep<=order; ep++) {
      fprintf(resfptr, "eps^%d: ", ep - 2*nloops);
      if (ep - 2*nloops >= 0) fprintf(resfptr, " ");
      fprintf(resfptr, "(");

      // print real part
      sign = mpfr_cmp_ui(mpc_realref(res[i][ep]), 0);
      if (sign == 0) {
        fprintf(resfptr, "0");
      } else {
        if (sign > 0) fprintf(resfptr, "+");
        mpfr_fprintf(resfptr, res_format, mpc_realref(res[i][ep]));
        // mpfr_out_str(resfptr, 10, 0, mpc_realref(res[i][ep]), MPFR_RNDN);
      }

      fprintf(resfptr, " ");

      // print imag part
      sign = mpfr_cmp_ui(mpc_imagref(res[i][ep]), 0);
      if (sign == 0) {
        fprintf(resfptr, "0");
      } else {
        if (sign > 0) fprintf(resfptr, "+");
        mpfr_fprintf(resfptr, res_format, mpc_imagref(res[i][ep]));
        // mpfr_out_str(resfptr, 10, 0, mpc_imagref(res[i][ep]), MPFR_RNDN);
      }
      fprintf(resfptr, ")\n");
    }
    fprintf(resfptr, "\n");
  }
  fclose(resfptr);
}

