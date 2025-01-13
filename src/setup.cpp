#include <iostream>
#include <cstring>
// #include <fstream>
// #include <sstream>
// #include <string>
// #include <cmath>
// #include <complex>
#include <openssl/sha.h>
#include <openssl/evp.h>
#include <sys/stat.h>
#include "mpc.h"

using namespace std;

#include "global_vars.h"
#include "malloc_defs.h"
#include "utils.h"
#include "tensor_utils.h"
// #include "system_analyzer.h"
#include "conversions.h"
#include "setup.h"

extern "C" {
  #include "in_out.h"
}


#define MAX_KEY_LENGTH 50
#define MAX_VALUE_LENGTH 2000


void read_mat_to_str(
  // OUTPUT
  char ****matrix, int *dim,
  // INPUT
  char *filepath
) {
  // Open file
  FILE *file = fopen(filepath, "r");
  if (file == NULL) {
		printf("Error opening file.\n");
		exit(1);
	}

	// Read the first line to determine the dimension
	char tmp;
	*dim = 1;
	while ((tmp = fgetc(file)) != EOF && tmp != '\n') {
		if (tmp == '\t') {
			(*dim)++;
		}
	}
	// printf("dim = %d\n", *dim);
	
	// Reset file position to the beginning
	fseek(file, 0, SEEK_SET);
	
	// Allocate memory for matrix
  int max_len = 10; // Initialize max dimension
	(*matrix) = (char***) malloc((*dim)*sizeof(char**));
	for (int i = 0; i < *dim; i++) {
		(*matrix)[i] = (char**) malloc((*dim)*sizeof(char*));
	}
	
	// Read the matrix elements
  int len;
	for (int i = 0; i < *dim; i++) {
    for (int j = 0; j < *dim; j++) {
			// cout << "i, j = " << i << ", " << j << endl;
      (*matrix)[i][j] = (char*) malloc(max_len*sizeof(char*));
		  len = 0;
		  // while (temp = fgetc(file)) != '\t' && temp != '\n') {
      while (1) {
        (*matrix)[i][j][len] = fgetc(file);
        // if (i == 26 && j == 0) cout << (*matrix)[i][j][len];

        // exit conditions
        if ((*matrix)[i][j][len] == '\t' || (*matrix)[i][j][len] == '\n') {
          // element (and possibly line) is ended
					(*matrix)[i][j][len] = '\0';
          break;
        }

			  if ((*matrix)[i][j][len] == '0' && len == 0) {
					// if element starts with zero, then it's zero
          (*matrix)[i][j][len+1] = '\0'; // terminate the string
					// skip the rest of the string
          char tmp = fgetc(file);
					while (!(tmp == '\t' || tmp == '\n')) {
					  tmp = fgetc(file);
					}
          break;
					// char tmp = fgetc(file);
					// if (tmp != '.') {
					// 	(*matrix)[i][j][len+1] = '\0'; // terminate the string
					// 	// skip the rest of the string
					// 	while (!(tmp == '\t' || tmp == '\n')) {
					// 		tmp = fgetc(file);
					// 	}
					// 	break;
					// } else {
					// 	tmp = fgetc(file);
					// 	if (tmp == '\t' || tmp == '\n') {
					// 		// file might contain zeros as 0.
					// 		(*matrix)[i][j][len+1] = '\0';
					// 		break;
					// 	} else {
					// 		// file contains floating point numbers
					// 		printf("Error while reading matrix from file: found '.' in input file (not a 0.)");
					// 		exit(0);
					// 	}
					// }
        }

        len++;

        // update max len if needed
        if (len == max_len) {
          max_len *= 2;
          (*matrix)[i][j] = (char*) realloc((*matrix)[i][j], max_len*sizeof(char));
        }
			}
      // cout << (*matrix)[i][j] << endl;
		}
	}

	fclose(file);
}


void compute_setup(
  // OUTPUT
  int *workprec, int *xorder, int *number, char ***eps_list,
  // INPUT
  int order, int precision , int loops
) {
  int i;

  *number = ceil( 5*order/2. + 2*loops);
  double exponent = loops/2. + precision/(order+1.);

  int singlepr = max((int) ceil((*number+2*loops)*exponent), 30);
  *workprec = 2*singlepr;
  *xorder = 4*singlepr;

  // EPSILON LIST
  double eps0 = pow(10.,-exponent);
  mpq_t r_eps0; mpq_init(r_eps0);
  rationalize(r_eps0, eps0, eps0);
  
  *eps_list = new char*[*number];
  mpq_t tmp; mpq_init(tmp);
  for (int ep=0; ep<*number; ep++) {
    mpz_mul_ui(mpq_numref(tmp), mpq_numref(r_eps0), 101+ep);
    mpz_mul_ui(mpq_denref(tmp), mpq_denref(r_eps0), 100);
    mpq_canonicalize(tmp);
    (*eps_list)[ep] = mpq_get_str(NULL, 10, tmp);
  }

  // // #hard-coded
  // *eps_list = new char*[*number];
  // if (loops == 1) {
  //   (*eps_list)[0] = strdup("101/146700");
  //   (*eps_list)[1] = strdup("17/24450");
  //   (*eps_list)[2] = strdup("103/146700");
  //   (*eps_list)[3] = strdup("26/36675");
  //   (*eps_list)[4] = strdup("7/9780");
  //   (*eps_list)[5] = strdup("53/73350");
  //   (*eps_list)[6] = strdup("107/146700");
  //   (*eps_list)[7] = strdup("3/4075");
  //   (*eps_list)[8] = strdup("109/146700");
  //   (*eps_list)[9] = strdup("11/14670");
  //   (*eps_list)[10] = strdup("37/48900");
  //   (*eps_list)[11] = strdup("28/36675");
  //   (*eps_list)[12] = strdup("113/146700");
  //   (*eps_list)[13] = strdup("19/24450");
  //   (*eps_list)[14] = strdup("23/29340");
  // } else if (loops == 2) {
  //   (*eps_list)[0] = strdup("101/464100");
  //   (*eps_list)[1] = strdup("1/4550");
  //   (*eps_list)[2] = strdup("103/464100");
  //   (*eps_list)[3] = strdup("2/8925");
  //   (*eps_list)[4] = strdup("1/4420");
  //   (*eps_list)[5] = strdup("53/232050");
  //   (*eps_list)[6] = strdup("107/464100");
  //   (*eps_list)[7] = strdup("9/38675");
  //   (*eps_list)[8] = strdup("109/464100");
  //   (*eps_list)[9] = strdup("11/46410");
  //   (*eps_list)[10] = strdup("37/154700");
  //   (*eps_list)[11] = strdup("4/16575");
  //   (*eps_list)[12] = strdup("113/464100");
  //   (*eps_list)[13] = strdup("19/77350");
  //   (*eps_list)[14] = strdup("23/92820");
  //   (*eps_list)[15] = strdup("29/116025");
  //   (*eps_list)[16] = strdup("3/11900");
  // }

}


void mpfr_tol_set() {
  mpfr_init2(mpfr_tol, wp2);
  mpfr_set_d(mpfr_tol, 2, MPFR_RNDN);
	mpfr_pow_si(mpfr_tol, mpfr_tol, -0.90*((int) wp2), MPFR_RNDN);
  // mpfr_to_gnc(&gnc_tol, &mpfr_tol);  // #uncomment-for-ginac
}


void mpfr_tol_set_wp(int workprec) {
  mpfr_init2(mpfr_tol, wp2);
  mpfr_set_d(mpfr_tol, 2, MPFR_RNDN);
  int wp2p = (10*workprec)/3;
	mpfr_pow_si(mpfr_tol, mpfr_tol, -0.90*(wp2p), MPFR_RNDN);
  // mpfr_to_gnc(&gnc_tol, &mpfr_tol);  // #uncomment-for-ginac
}


void info_to_log(
  // INPUT
  FILE *logfptr,
  int loops, int order, int precision, int eta_ord, int workprec, int wp2
) {
  fprintf(logfptr, "\nINPUT PARAMETERS:\n");
  fprintf(logfptr, "loops = %d\n", loops);
  fprintf(logfptr, "order = %d\n", order);
  fprintf(logfptr, "precision = %d\n", precision);
  fprintf(logfptr, "\nINTERNAL PARAMETERS:\n");
  fprintf(logfptr, "eta_ord = %d\n", eta_ord);
  fprintf(logfptr, "workprec = %d\n", workprec);
  fprintf(logfptr, "binary workprec = %d\n", wp2);
  fprintf(logfptr, "mpfr tolerance = "); mpfr_out_str(logfptr, 10, 0, mpfr_tol, MPFR_RNDN); fprintf(logfptr, "\n");
}


void parse_list(
  // OUTPUT
  char ***out, int *dim,
  // INPUT
  char *list
) {
  // eliminate spaces
  char *ns_list = (char*) malloc(MAX_VALUE_LENGTH*sizeof(char*));
  int i = 0, c = 0;
  while (list[i] != '\0') {
    if (!isspace(list[i])) {
      ns_list[c++] = list[i];
    }
    i++;
  }
  ns_list[i] = '\0';
  // printf("no-space list: %s\n", ns_list);

  // read delimiter
  char rdel;
  switch (ns_list[0]) {
    case '(':
    rdel = ')';
    break;
    case '[':
    rdel = ']';
    break;
    case '{':
    rdel = '}';
    break;
    default:
    perror("left delimiter not allowed");
    exit(1);
    break;
  }
  // cout << "rdel = " << rdel << endl;

  // count number of elements
  *dim = 0, i = 0;
  while (ns_list[i] != rdel) {
    if (ns_list[i] == ',') {
      (*dim)++;
    }
    i++;
  }
  if (ns_list[i-1] != ',') {
  (*dim)++;
  }
  // cout << "dim = " << *dim << endl;

  // alloc memory for output
  if ((*out)) {
    delete[] (*out);
  }
  (*out) = new char*[*dim];
  i = 1; // skip left delimiter
  for (int k=0; k<*dim; k++) {
    (*out)[k] = (char*) malloc(MAX_VALUE_LENGTH*sizeof(char));
    c = 0;
    while (ns_list[i] != ',' && ns_list[i] != rdel) {
      (*out)[k][c++] = ns_list[i++];
    }
    (*out)[k][c] = '\0';
    i++;
  }

  // FREE
  free(ns_list);
}


int sym_to_idx(
  char *key, char **symbol, int ninvs
) {
  for (int s=0; s<ninvs; s++) {
    if (strcmp(key, symbol[s]) == 0) {
      return s;
    }
  }

  // if no match is found, return error
  perror("key not found in symbol array");
  printf("key: %s\n", key);
  exit(1);
}


int end_at_char(
  // IN-OUT
  char *line,
  // INPUT
  char ch
) {
  int i = 0;
  while(line[i] != '\0') {
    if (line[i] == ch) {
      // terminate string at target character
      line[i] = '\0';
      i = 0;
      return 1;
    }
    i++;
  }

  return 0;
}


int get_key(
  // OUTPUT
  char *key,
  char *line
) {
  // Find the colon and extract the key
  int i=0, key_length = 0;
  while (line[i] != ':' && line[i] != '\0') {
    if (!isspace(line[i])) {
      key[key_length++] = line[i];
    }
    i++;
  }
  key[key_length] = '\0';

  if (key_length) {
    return 1;
  } else {
    return 0;
  }
}


int goto_key(
  // OUTPUT
  char **line,
  // INPUT
  char *key, FILE *fptr
) {
  char tmp_key[5000];  // #2BD: set to lower value but change get_key function
	while (fgets(*line, MAX_VALUE_LENGTH, fptr)) {
    // cout << "line:" << endl; cout << line << endl;
    // skip comment lines
    if ((*line)[0] == '#') {
      continue;
    }

    end_at_char(*line, '#');

    if (get_key(tmp_key, *line)) {
      if (strcmp(tmp_key, key) == 0) {
        return 1;
        break;
      }
    } else {
      continue;
    }
  }

  return 0;
}


int read_param_str(
  // OUTPUT
  char **value,
  // INPUT
  char *key, FILE *fptr
) {
  const size_t line_size = MAX_VALUE_LENGTH * sizeof(char);
  char *line = (char*) malloc(line_size);
  if (line == NULL) {
    perror("error while allocating memory for line");
    exit(1);
  }
  if (!goto_key(&line, key, fptr)) {
    return 0;
  }

  // go to ':'
  int i = 0;
  while(line[i] != ':') {
    i++;
  };
  i++;

  // skip spaces
  while(isspace(line[i])) {
    i++;
  }

  // count length of value, skip spaces
  int len=0;
  while (!isspace(line[i+len])) {
    len++;
  }
  line[i+len] = '\0';

  *value = (char*) malloc((len+1)*sizeof(char));
  strcpy(*value, line+i);
  free(line);
  return 1;
}


int read_param_int(
  // OUTPUT
  int *value,
  // INPUT
  char *key, FILE *fptr
) {
  unsigned long int line_size = MAX_VALUE_LENGTH * sizeof(char);
  char *line = (char*) malloc(line_size);
  if (!goto_key(&line, key, fptr)) {
    free(line);
    return 0;
  }

  // go to ':'
  int i = 0;
  while(line[i] != ':') {
    i++;
  };
  i++;

  // skip spaces
  while(isspace(line[i])) {
    i++;
  }

  // extract value
  sscanf(line+i, "%d", value);
  free(line);
  return 1;
}


int read_list_str(
  // OUTPUT
  char ***list,
  // INPUT
  int *dim, char *key, FILE *fptr
) {
  unsigned long int line_size = MAX_VALUE_LENGTH * sizeof(char);
  char *line = (char*) malloc(line_size);
  if (!goto_key(&line, key, fptr)) {
    free(line);
    return 0;
  };
  
  // look for "[", ignore everything else
  int i = 0, security_count = 0;
  while (line[i] != '[') {
    security_count++;
    if (security_count == 500) {
      // ensure end of program if '[' is not found
      perror("error while parsing phase-space point: left delimiter missing or too many characters");
      fclose(fptr);
      free(line);
      exit(1);
    }
    i++;
  }
  i++; // skip '['

  int pos, num_len, security_count_comma;
  for (int s=0; s<*dim; s++) {
    // cout << "s = " << s << endl;

    // skip any space, tab or new line
    // cout << "skip any space, tab or newline" << endl;
    while (isspace(line[i])) {
      // read new line if necessary
      if (line[i] == '\n') {
        fgets(line, line_size, fptr);
        i = 0;
      } else {
        i++;
      }
    }

    // store list element
    num_len = 0;
    security_count_comma = 0;
    while(line[i+num_len] != ',' && line[i+num_len] != ']' && line[i+num_len] != '\n') {
      // if (isspace(line[i+num_len])) {
      //   perror("error while parsing phase-space point: found space, tab or end of line in symbol name\n");
      //   exit(1);
      // }
      security_count_comma++;
      if (security_count_comma == 500) {
        perror("error while parsing phase-space point: separator missing or too many characters");
        fclose(fptr);
        free(line);
        exit(1);
      }
      num_len++;
    }
    // cout << "numerical length = " << num_len << endl;
    (*list)[s] = (char*) malloc((num_len+1)*sizeof(char));
    for (int c=0; c<num_len; c++) {
      (*list)[s][c] = line[i];
      i++;
    }
    (*list)[s][num_len] = '\0';

    // skip any space, tab or new line
    // cout << "skip any space, tab or newline" << endl;
    while (isspace(line[i])) {
      // read new line if necessary
      if (line[i] == '\n') {
        fgets(line, line_size, fptr);
        i = 0;
      } else {
        i++;
      }
    }

    // skip spaces, ',' and ']'
    if (line[i] == ',') {
      i++;
      continue;
    } else if (line[i] == ']') {
      if (s == *dim - 1) {
        break;
      } else {
        perror("error while parsing phase-space point: input point has fewer invariants than symbols");
        fclose(fptr);
        free(line);
        exit(1);
      }
    }
  }

  free(line);
  return 1;
}


int read_PS_point(
  // OUTPUT
  char ***kin,
  // INPUT
  int *ninvs, char **symbols, char *key, FILE *fptr
) {
  unsigned long int line_size = MAX_VALUE_LENGTH * sizeof(char);
  char *line = (char*) malloc(line_size);
  if (!goto_key(&line, key, fptr)) {
    free(line);
    return 0;
  };
  
  // look for "[", ignore everything else
  int i = 0, security_count = 0;
  while (line[i] != '[') { 
    security_count++;
    if (security_count == 500) {
      // ensure end of program if '[' is not found
      perror("error while parsing phase-space point: left delimiter missing or too many characters");
      fclose(fptr);
      free(line);
      exit(1);
    }
    i++;
  }
  i++; // skip '['

  int pos, num_len, security_count_comma, sym_idx;
  for (int s=0; s<*ninvs; s++) {
    // cout << "s = " << s << endl;

    // skip any space, tab or new line
    // cout << "skip any space, tab or newline" << endl;
    while (isspace(line[i])) {
      // read new line if necessary
      if (line[i] == '\n') {
        fgets(line, line_size, fptr);
        i = 0;
      } else {
        i++;
      }
    }

    // get symbol name
    char *sym_key = (char*) malloc(50*sizeof(char));
    pos = 0;
    while(!isspace(line[i]) && line[i] != '=') {
      sym_key[pos++] = line[i];
      i++;
    }
    sym_key[pos] = '\0';
    // cout << "sym key: " << sym_key << endl;
    sym_idx = sym_to_idx(sym_key, symbols, *ninvs);
    // cout << "sym idx = " << sym_idx << endl;

    // skip spaces and '='
    while(isspace(line[i]) || line[i] == '=') {
      i++;
    }

    // store numerical value
    num_len = 0;
    security_count_comma = 0;
    while(line[i+num_len] != ',' && line[i+num_len] != ']' && line[i+num_len] != '\n') {
      // if (isspace(line[i+num_len])) {
      //   perror("error while parsing phase-space point: found space, tab or end of line in symbol name\n");
      //   exit(1);
      // }
      security_count_comma++;
      if (security_count_comma == 500) {
        perror("error while parsing phase-space point: separator missing or too many characters");
        fclose(fptr);
        free(line);
        exit(1);
      }
      num_len++;
    }
    // cout << "numerical length = " << num_len << endl;
    (*kin)[sym_idx] = (char*) malloc((num_len+1)*sizeof(char));
    for (int c=0; c<num_len; c++) {
      (*kin)[sym_idx][c] = line[i];
      i++;
    }
    (*kin)[sym_idx][num_len] = '\0';

    // skip any space, tab or new line
    // cout << "skip any space, tab or newline" << endl;
    while (isspace(line[i])) {
      // read new line if necessary
      if (line[i] == '\n') {
        fgets(line, line_size, fptr);
        i = 0;
      } else {
        i++;
      }
    }

    // skip spaces, ',' and ']'
    if (line[i] == ',') {
      i++;
      continue;
    } else if (line[i] == ']') {
      if (s == *ninvs - 1) {
        break;
      } else {
        perror("error while parsing phase-space point: input point has fewer invariants than symbols");
        fclose(fptr);
        free(line);
        exit(1);
      }
    }
  }

  free(line);
  return 1;
}


int read_bound_factors(
  // OUTPUT
  char ****out, int *nfac,
  // INPUT
  int nfac_max, int nparam, char ***symbols, char *key, FILE *fptr
) {
  unsigned long int line_size = MAX_VALUE_LENGTH * sizeof(char);
  char *line = (char*) malloc(line_size);
  if (!goto_key(&line, key, fptr)) {
    free(line);
    return 0;
  };

  // look for "{", ignore everything else
  int i = 0, security_count = 0;
  while (line[i] != '{') { 
    security_count++;
    if (security_count == 500) {
      // ensure end of program if '{' is not found
      perror("error while parsing phase-space point: left delimiter missing or too many characters");
      fclose(fptr);
      free(line);
      exit(1);
    }
    i++;
  }
  i++; // skip '{'
  // cout << "i = " << i << endl;

  for (*nfac=0; *nfac<nfac_max; (*nfac)++) {
    // cout << "nfac = " << *nfac << endl;
    // look for "[", ignore everything else
    security_count = 0;
    while (line[i] != '[') { 
      security_count++;
      if (security_count == 500) {
        // ensure end of program if '[' is not found
        perror("error while parsing phase-space point: left delimiter missing or too many characters");
        fclose(fptr);
        free(line);
        exit(1);
      } else if (line[i] == '}') {
        free(line);
        return 1;
      }
      i++;
    }
    i++; // skip '['
    // cout << "i = " << i << endl;

    int pos, num_len, security_count_comma, sym_idx;
    for (int s=0; s<nparam; s++) {
      // cout << "s = " << s << endl;

      // skip any space, tab or new line
      // cout << "skip any space, tab or newline" << endl;
      while (isspace(line[i])) {
        // read new line if necessary
        if (line[i] == '\n') {
          fgets(line, line_size, fptr);
          i = 0;
        } else {
          i++;
        }
      }

      // get symbol name
      char *sym_key = (char*) malloc(50*sizeof(char));
      pos = 0;
      while(!isspace(line[i]) && line[i] != '=') {
        sym_key[pos++] = line[i];
        i++;
      }
      sym_key[pos] = '\0';
      // cout << "sym key: " << sym_key << endl;
      sym_idx = sym_to_idx(sym_key, *symbols, nparam);
      // cout << "sym idx = " << sym_idx << endl;
      // exit(0);

      // skip spaces and '='
      while(isspace(line[i]) || line[i] == '=') {
        i++;
      }

      // store numerical value
      num_len = 0;
      security_count_comma = 0;
      int parenthesis_count = 0;
      while(
        parenthesis_count != 0 ||\
        (line[i+num_len] != ',' && line[i+num_len] != ']' && line[i+num_len] != '\n')
      ) {
        if (line[i+num_len] == '(') {
          parenthesis_count++;
        }
        if (line[i+num_len] == ')') {
          parenthesis_count--;
        }
        // if (isspace(line[i+num_len])) {
        //   perror("error while parsing phase-space point: found space, tab or end of line in symbol name\n");
        //   exit(1);
        // }
        security_count_comma++;
        if (security_count_comma == 500) {
          perror("error while parsing phase-space point: separator missing or too many characters");
          fclose(fptr);
          free(line);
          exit(1);
        }
        num_len++;
      }
      // cout << "numerical length = " << num_len << endl;
      (*out)[*nfac][sym_idx] = (char*) malloc((num_len+1)*sizeof(char));
      for (int c=0; c<num_len; c++) {
        (*out)[*nfac][sym_idx][c] = line[i];
        i++;
      }
      (*out)[*nfac][sym_idx][num_len] = '\0';

      // skip any space, tab or new line
      // cout << "skip any space, tab or newline" << endl;
      while (isspace(line[i])) {
        // read new line if necessary
        if (line[i] == '\n') {
          fgets(line, line_size, fptr);
          i = 0;
        } else {
          i++;
        }
      }

      // skip spaces, ',' and ']'
      if (line[i] == ',') {
        // cout << "found ',' with s = " << s << endl;
        i++;
        continue;
      } else if (line[i] == ']') {
        if (s == nparam - 1) {
          // cout << "found ']' with s = " << s << endl;
          break;
        } else {
          perror("error while parsing phase-space point: input point has fewer invariants than symbols");
          fclose(fptr);
          free(line);
          exit(1);
        }
      }
    }

    // skip any space, tab or new line
    // cout << "skip any space, tab or newline" << endl;
    while (isspace(line[i])) {
      // read new line if necessary
      if (line[i] == '\n') {
        fgets(line, line_size, fptr);
        i = 0;
      } else {
        i++;
      }
    }
    i++;

    // skip spaces, ',' and '}'
    if (line[i] == ',') {
      // cout << "found ',' with nfac = " << *nfac << endl;
      i++;
      continue;
    } else if (line[i] == '}') {
      if (*nfac < nfac_max) {
        // cout << "found '}' with nfac = " << *nfac << endl;
        break;
      } else if (*nfac >= nfac_max) {        
        perror("error while parsing boundary factors: input point has more factors than maximum allowed");
        fclose(fptr);
        free(line);
        exit(1);
      }
    }
  }

  (*nfac)++;
  free(line);
  return 1;
}


// <<load_parameters>>
void load_parameters(
  // OUTPUT
  int *loops, int *order, int *precision,
  char **prev_run,
  char ***kin,
  // INPUT
  int *ninvs, char ***symbols,
  char *filename
) {
	FILE *file = fopen(filename, "r");
	if (file == NULL) {
		printf("Error opening file.\n");
		return;
	}

  if (prev_run) {
    *prev_run = NULL;
  }

	char line[2000]; // Assuming each line has at most X characters

	while (fgets(line, sizeof(line), file)) {
    // cout << "line:" << endl; cout << line << endl;
    // skip comment lines
    if (line[0] == '#') {
      continue;
    }

    // skip in-line comments
    int i = 0;
    while(line[i] != '\0') {
      if (line[i] == '#') {
        line[i] = '\0';
        i = 0;
        break;
      }
      i++;
    }
    i = 0;
		
		char key[50]; // Assuming each key has at most 20 characters
		char value_str[50]; // Assuming each value has at most 20 characters
		int key_length = 0;

		// Find the colon and extract the key
		while (line[i] != ':' && line[i] != '\0') {
			if (!isspace(line[i])) {
				key[key_length++] = line[i];
			}
			i++;
		}
		key[key_length] = '\0';
    // cout << "key: " << key << endl;

		// Skip whitespace after the colon
    i++;
		while (isspace(line[i])) {
			i++;
		}

    //////
    // NON-NUMERICAL OPTIONS
    //////
    if (strcmp(key, "use-prev-run") == 0) {
      *prev_run = (char*) malloc(200*sizeof(char));
      int j=0;
      while (line[i] != '\0' && line[i] != '\n') {
        (*prev_run)[j++] = line[i++];
      }
      (*prev_run)[j] = '\0';
      continue;
    }

    //////
    // NON-NUMERICAL OPTIONS
    //////
    // SYMBOLS
    if (strcmp(key, "symbols") == 0) {
      // cout << "reading symbols..." << endl;
      char *sym_str = (char*) malloc(2000*sizeof(char));
      int sym_str_pos = 0;
      *ninvs = 0;

      int security_count = 0;
      // look for "[", ignore everything else
      while (line[i] != '[') {
        security_count++;
        if (security_count == 500) {
          // ensure end of program if '[' is not found
          perror("error while parsing symbols: left delimiter missing or too many characters");
          fclose(file);
          exit(1);
        }
        i++;
      }
      i++; // skip '['
      
      security_count = 0;
      while (line[i] != ']') {
        // cout << "next char: "; cout << line[i] << endl;

        // ensure end of program if ']' is not found
        security_count++;
        if (security_count == 500) {
          perror("error while parsing symbols: right delimiter missing or too many characters");
          fclose(file);
          exit(1);
        }

        // read new line if necessary
        if (line[i] == '\n') {
          fgets(line, sizeof(line), file);
          i = 0;
          continue;
        }

        // skip any space or tab
        while (isspace(line[i])) {
          i++;
        }
        // cout << "after spaces: " << line[i] << endl;

        // anything before "," is a symbol
        // cout << "found symbol" << endl;
        (*ninvs)++;
        // cout << "ninvs = " << *ninvs << endl;
        int security_count_comma = 0;
        // cout << "reading up to next separator..." << endl;
        while (line[i] != ',' && line[i] != ']') {
          if (isspace(line[i])) {
            perror("error while parsing symbols: found space, tab or end of line in symbol name");
            fclose(file);
            exit(1);
          }
          // cout << line[i] << endl;
          // ensure end of program if ',' or ']' is not found
          security_count_comma++;
          if (security_count_comma == 500) {
            perror("error while parsing symbols: separator missing or too many characters");
            fclose(file);
            exit(1);
          }

          sym_str[sym_str_pos++] = line[i];
          i++;
        }
        if (line[i] == ',') {i++;}
        sym_str[sym_str_pos++] = '$';

        // skip any space or tab
        while (isspace(line[i])) {
          i++;
        }

      }
      sym_str[sym_str_pos] = '\0';
      // cout << "sym_str: " << sym_str << endl;

      // convert symbol string into symbol array
      *symbols = (char**) malloc((*ninvs+1)*sizeof(char*));
	    (*symbols)[0] = (char*) malloc(sizeof(char)+1);
	    (*symbols)[0] = (char*) "d";
      sym_str_pos = 0;
      int sym_len, sym_pos;
      for (int s=1; s<=*ninvs; s++) {
        // cout << "s = " << s << endl;
        // get length of the symbol
        sym_len = 0;
        while (sym_str[sym_str_pos+sym_len] != '$') {
          sym_len++;
        }
        // cout << "sym len = " << sym_len << endl;

        // allocate memory for symbol
        (*symbols)[s] = (char*) malloc((sym_len+1)*sizeof(char));
        sym_pos = 0;
        while (sym_str[sym_str_pos] != '$') {
          (*symbols)[s][sym_pos++] = sym_str[sym_str_pos];
          sym_str_pos++;
        }
        sym_str_pos++;
        // cout << "s = " << s << ": " << (*symbols)[s] << endl;
      }
      continue;
    }

    //////
    // NUMERIC OPTIONS
    //////

    // Extract value
		int j = 0;
		while (line[i] != '\0' && line[i] != '\n') {
      if (!isspace(line[i])) {
        value_str[j++] = line[i++];
      } else {
        i++;
      }
		}
		value_str[j] = '\0';
    // cout << "value string: " << value_str << endl;

		// Convert the value string to an integer
		int value = atoi(value_str);
    // cout << "value = " << value << endl;

		// Assign the value based on the key
		if (strcmp(key, "loops") == 0) {
			*loops = value;
		} else if (strcmp(key, "precision") == 0) {
			*precision = value;
		} else if (strcmp(key, "order") == 0) {
			*order = value;
		}
	}

  //////
  // ON-DEMAND
  //////
  // PS POINT
  if (*ninvs == 0) {
    perror("error: parsing phase-space point before symbols");
    fclose(file);
    exit(1);
  }
  *kin = new char*[*ninvs];

  fseek(file, 0, SEEK_SET);
  if (!read_PS_point(
    kin,
    ninvs, *symbols+1, (char*)"point", file
  )) {
    perror("error while reading file: PS point not found");
    fclose(file);
    exit(1);
  }
  
	fclose(file);
}


void param_string(
  // OUTPUT
  char **out,
  // INPUT
  int num, char **keys, int *values
) {
  /*
  assemble string in the form:
    key1_value1_key2_val2  
  */

  // allocate memory for the output string
  int tot_len = 0, max_val_len = 5;
  for (int i=0; i<num; i++) {
    tot_len += strlen(keys[i]) + 1 + max_val_len + 1;
  }
  *out = (char*)malloc((tot_len + 1)*sizeof(char));

  // assemble string
  char sep = '_';
  // int pos = strlen(*out);
  int pos = 0;
  for (int k=0; k<num; k++) {
    if (k>0) {
      (*out)[pos++] = sep;
    }

    // attach key
    for (int i=0; keys[k][i] != '\0'; i++) {
      (*out)[pos++] = keys[k][i];
    }
    (*out)[pos++] = sep;

    // attacch value
    char *value_str = (char*) malloc((max_val_len + 1)*sizeof(char));
    snprintf(value_str, max_val_len+1, "%d", values[k]);
    for (int i=0; value_str[i] != '\0'; i++) {
      (*out)[pos++] = value_str[i];
    }
  }
  (*out)[pos] = '\0';
  
  return;
}


void point_string(
	// OUTPUT
	char **out,
	// INPUT
	// int order, int precision,
	int ninvs, char** symbols, char** PS
) {
	// allocate memory for the output string
	int total_length = 0;
	// char order_str[50], precision_str[50];
	// snprintf(order_str, sizeof(order_str), "%d", order);
	// total_length += strlen(order_str) + 1;
  // snprintf(precision_str, sizeof(precision_str), "%d", precision);
	// total_length += strlen(precision_str) + 1;
	for (int s=0; s<ninvs; s++) {
		total_length += strlen(symbols[s]) + 1 + strlen(PS[s]) + 1; // +1 for the separator "_"
	}

	*out = (char*)malloc(total_length + 1);
	if (*out == NULL) {
		printf("Error in memory allocation");
		exit(1);
	}

	// build string
	char sep = '_';

  // // insert parameters
	// strcat(*out, order_str);
	// strcat(*out, "_");
	// strcat(*out, precision_str);
	// strcat(*out, "_");
	
  // int pos = strlen(*out);
  int pos = 0;
	for (int s=0; s<ninvs; s++) {
		if (s > 0) {
			(*out)[pos++] = sep;
		}

		// insert symbol
		for (int i = 0; symbols[s][i] != '\0'; i++) {
			(*out)[pos++] = symbols[s][i];
		}
		(*out)[pos++] = sep;

		// insert value
		for (int i = 0; PS[s][i] != '\0' && PS[s][i] != '\n'; i++) {
			if (PS[s][i] == '/') {
				(*out)[pos++] = '%';
			} else if (PS[s][i] == '\n') {
				pos++;
			} else {
				(*out)[pos++] = PS[s][i];
			}
		}

	}
	(*out)[pos] = '\0';
	
	return;
}


void target_string(
	// OUTPUT
	char **out,
	// INPUT
	int nparam, char **param_keys, int *param_values,
	int ninvs, char** symbols, char** PS
) {

  // build parameter string
  char *param_str;
  param_string(
    &param_str,
    nparam, param_keys, param_values
  );

  // build point string
  char *pt_str;
	point_string(
		&pt_str,
		ninvs, symbols, PS
	);

  // concatenate parameters and point
  (*out) = (char*) malloc((strlen(param_str) + 1 + strlen(pt_str) + 1)*sizeof(char));
  strcpy(*out, param_str);
  strcat(*out, "_");
  strcat(*out, pt_str);
  
}


// Function to convert an unsigned char* hash to a char* string
char* hash_to_string(unsigned char *hash) {
	// Length of the SHA256 hash (in bytes)
	int hash_length = 32;

	// Allocate memory for the output string
	char *output = (char*)malloc((2 * hash_length + 1) * sizeof(char));
	if (output == NULL) {
		printf("Error allocating memory for output string");
		exit(1);
	}

	// Convert each byte of the hash to two hexadecimal characters
	for (int i = 0; i < hash_length; i++) {
		snprintf(output + 2 * i, sizeof(output + 2 * i), "%02x", hash[i]);
	}

	// Add the null terminator
	output[2 * hash_length] = '\0';

	return output;
}


void calculate_hash(
	// OUTPUT
	char **out,
	// INPUT
	char *in
) {

	unsigned char *hash_str = (unsigned char*) malloc(SHA256_DIGEST_LENGTH);
	
	// // #DEPRECATED IN OPENSSL 3.3.0
	// SHA256_CTX context;
	// SHA256_Init(&context);
	// SHA256_Update(&context, in, strlen(in));
	// SHA256_Final(hash_str, &context);

  EVP_MD_CTX *mdctx = EVP_MD_CTX_new();
  EVP_DigestInit_ex(mdctx, EVP_sha256(), NULL);
  EVP_DigestUpdate(mdctx, in, strlen(in));
  unsigned int out_len;
  EVP_DigestFinal_ex(mdctx, hash_str, &out_len);

	// printf("hash: ");
  // for(int i=0; i<SHA256_DIGEST_LENGTH; i++) {
	// 	printf("%02x", hash_str[i]);
  // }
  // printf("\n");

	*out = hash_to_string(hash_str);

}


int file_exists(const char *path) {
	struct stat buffer;
  return (stat(path, &buffer) == 0);
}


int directory_exists(const char *path) {
	struct stat buffer;
  return (stat(path, &buffer) == 0 && S_ISDIR(buffer.st_mode));
}


int make_dir(char* path) {
	if (directory_exists(path)) {
		return 1;
	} else {
		if (mkdir(path, 0777) == -1) {
			perror("Error while creating directory");
			exit(1);
		}
		return 0;
	}
}


void log_param(
  char *key, int value,
  char *filepath, const char* mode
) {
  FILE *fptr = fopen(filepath, mode);
  fprintf(fptr, "%s: %d\n", key, value);
  fclose(fptr);
}


void log_params(
  int num, char **keys, int *values,
  char *filepath, const char* mode
) {
  FILE *fptr = fopen(filepath, mode);
  for (int k=0; k<num; k++) {
    fprintf(fptr, "%s: %d\n", keys[k], values[k]);
  }
  fclose(fptr);
}


void log_point(
  int ninvs, char **symbols, char **kin,
  char *key, char *filepath, const char* mode
) {
  FILE *fptr = fopen(filepath, mode);
  fprintf(fptr, "%s: [", key);
  for (int k=0; k<ninvs; k++) {
    fprintf(fptr, "\n\t%s = %s", symbols[k], kin[k]);
    if (k != ninvs-1) {
      fprintf(fptr, ",");
    }
  }
  fprintf(fptr, "\n]\n");
  fclose(fptr);
}


void prepare_prev_run(
  // OUTPUT
  char **filepath_prev_run_PS, char **filepath_prev_run_sol,
  // INPUT
  char *dir_parent, char *prev_run
) {
  char path_common[200];
  strcpy(path_common, dir_parent);
  strcat(path_common, "runs/");
  strcat(path_common, prev_run);

  join_path(filepath_prev_run_PS, path_common, (char*)"/input/PS_points.txt");
  if (*filepath_prev_run_sol) {
    delete[] *filepath_prev_run_sol;
  }
  join_path(filepath_prev_run_sol, path_common, (char*)"/output/sol");
}


// <<load_symbols>>
void load_symbols(
	// OUTPUT
	int *ninvs, char ***symbols,
	// IN-OUT
	char *filepath_vars
) {
  FILE *file = fopen(filepath_vars, "r");
	if (file == NULL) {
		perror("Error opening file");
		exit(EXIT_FAILURE);
  }

	// count the number of lines in the file
	*ninvs = 0;
	char *line = NULL;
	size_t line_len = 0;
	while (getline(&line, &line_len, file) != -1) {
		(*ninvs)++;
	}
	rewind(file);

	// allocate memory for the symbols array
	*symbols = (char**) malloc((*ninvs+1)*sizeof(char*));
	(*symbols)[0] = (char*) malloc(sizeof(char)+1);
	(*symbols)[0] = (char*) "d";

	// read symbols from file
	for (int i=1; i<=*ninvs; i++) {
		(*symbols)[i] = NULL;
		getline(&(*symbols)[i], &line_len, file);
		(*symbols)[i][strlen((*symbols)[i])-1] = '\0';
		// printf("line address = %p\n", symbols[i]);
		// printf("line size = %zu\n", line_len);
		// printf("line len = %zu\n", strlen(symbols[i]));
		// printf("line = %s\n", symbols[i]);
		// printf("last char = %c\n", symbols[i][strlen(symbols[i])-1]);
	}

  fclose(file);
}


// <<detect_masses>>
void detect_masses(
  // OUTPUT
  int *nmass, int **mass, int **is_mass,
  // INPUT
  int ninvs, char **symbols
) {
  // count masses
  *nmass = 0;
  for (int s=1; s<=ninvs; s++) {
    if (symbols[s][0] == 'm') {
      (*nmass)++;
    }
  }
  *mass = new int[*nmass];
  *is_mass = new int[ninvs];

  // get mass indices
  int count = 0;
  for (int s=1; s<=ninvs; s++) {
    if (symbols[s][0] == 'm') {
      (*is_mass)[s-1] = 1;
      (*mass)[count] = s-1;
      (count)++;
    } else {
      (*is_mass)[s-1] = 0;
    }
  }
}


// <<load_DE_matrices>>
void load_DE_matrices(
  // OUTPUT
  char ****mats_str, int *dim,
  // INPUT
  int ninvs, char **symbols, int *skip_inv,
  char *filepath_mats, char *file_ext
) {
	for (int s=0; s<ninvs; s++) {
    if (skip_inv[s] == 1) {continue;}
		char tmp_filepath[200];
		snprintf(tmp_filepath, sizeof(tmp_filepath), "%s%d%s", filepath_mats, s, file_ext);
		read_mat_to_str(&mats_str[s], dim, tmp_filepath);
	}
}


void str_rk1_to_mpc_rk1(
  // OUTPUT
  mpc_t *out,
  // INPUT
  char **kin, int ninvs
) {
  for (int s=0; s<ninvs; s++) {
    mpc_init3(out[s], wp2, wp2);
    mpc_set_str_rat(&out[s], kin[s]);
  }
}


void detect_active_invs(
  // OUTPUT
  int *skip_inv,
  // INPUT
  int ninvs, mpc_t *PS_ini, mpc_t *PS_fin
) {
  for (int s=0; s<ninvs; s++) {
    if (mpc_equal_within_tol(PS_ini[s], PS_fin[s])) {
      skip_inv[s] = 1;
    } else {
      skip_inv[s] = 0;
    }
  }
}


// <<load_PS>>
void load_PS(
  // OUT
  mpc_t *PS_ini, mpc_t *PS_fin,
  int *skip_inv, char ***kin,
  // INPUT
  int ninvs, char *filepath_PS
) {
  FILE *fptr = fopen(filepath_PS, "r");
	if (fptr == NULL) {
		perror("Error opening file");
    printf("file: %s", filepath_PS);
		exit(EXIT_FAILURE);
  }

  size_t len = 0;
  // cout << "initial point:" << endl;
  for (int s=0; s<ninvs; s++) {
    kin[0][s] = NULL;
    getline(&kin[0][s], &len, fptr);
    mpc_init3(PS_ini[s], wp2, wp2);
    mpc_set_str_rat(&PS_ini[s], kin[0][s]);
  }

  // cout << "final point:" << endl;
  char *tmp_line = NULL;
  getline(&tmp_line, &len, fptr); // skip a line
	for (int s=0; s<ninvs; s++) {
    kin[1][s] = NULL;
    getline(&kin[1][s], &len, fptr);
    mpc_init3(PS_fin[s], wp2, wp2);
    mpc_set_str_rat(&PS_fin[s], kin[1][s]);
	}
  fclose(fptr);

  // identify active invariants
  // (i.e. invariants that are actually varied)
  // cout << "active invariants:" << endl;
  mpc_t tmpc;
  mpc_init3(tmpc, wp2, wp2);
  for (int s=0; s<ninvs; s++) {
    mpc_sub(tmpc, PS_fin[s], PS_ini[s], MPFR_RNDN);
    if (mpc_lessthan_tol(tmpc)) {
      skip_inv[s] = 1;
    } else {
      skip_inv[s] = 0;
    }
    // cout << "s: " << s << ", skip: " << skip_inv[s] << endl;
  }
}


// <<load_PS_prev_run>>
void load_PS_prev_run(
  // OUT
  mpc_t *PS_ini, mpc_t *PS_fin,
  int *skip_inv, char ***kin,
  // INPUT
  int ninvs, char *filepath_PS
) {
  FILE *fptr = fopen(filepath_PS, "r");
  size_t len = 0;

  // skip lines of previous run initial point
  char *skip_line = NULL;
  for (int s=0; s<=ninvs; s++) {
    getline(&skip_line, &len, fptr);
  }

  // load previous run final point as initial point
  for (int s=0; s<ninvs; s++) {
    if (kin[0][s]) {
      delete[] kin[0][s];
    }
    kin[0][s] = NULL;
    getline(&kin[0][s], &len, fptr);
    mpc_init3(PS_ini[s], wp2, wp2);
    mpc_set_str_rat(&PS_ini[s], kin[0][s]);
  }

  // identify active invariants
  // (i.e. invariants that are actually varied)
  // cout << "active invariants:" << endl;
  mpc_t tmpc;
  mpc_init3(tmpc, wp2, wp2);
  for (int s=0; s<ninvs; s++) {
    mpc_sub(tmpc, PS_fin[s], PS_ini[s], MPFR_RNDN);
    if (mpc_lessthan_tol(tmpc)) {
      skip_inv[s] = 1;
    } else {
      skip_inv[s] = 0;
    }
    // cout << "s: " << s << ", skip: " << skip_inv[s] << endl;
  }
}


// <<linearize_masses_mpc>>
void linearize_masses_mpc(
  // IN-OUT
  mpc_t *PS_ini, mpc_t *PS_fin,
  // INPUT
  int nmass, int *mass
) {
  for (int m=0; m<nmass; m++) {
    mpc_sqrt(PS_ini[mass[m]], PS_ini[mass[m]], MPFR_RNDN);
    mpc_sqrt(PS_fin[mass[m]], PS_fin[mass[m]], MPFR_RNDN);
  }
}


// <<linearize_masses_str>>
void linearize_masses_str(
  // OUTPUT
  char ***lin_symbols,
  // INPUT
  int ninvs, char **symbols, int *is_mass
) {
  (*lin_symbols) = new char*[ninvs+1];

  // copy space-time dim
  (*lin_symbols)[0] = strdup(symbols[0]);

  // change symbol of mass variables (quadratic to linear)
  int len;
  for (int s=1; s<=ninvs; s++) {
    if (is_mass[s-1]) {
      len = strlen(symbols[s]);
      (*lin_symbols)[s] = new char[len+1];
      strcpy((*lin_symbols)[s], symbols[s]);
      (*lin_symbols)[s][len-1] = '\0';
    } else {
      (*lin_symbols)[s] = strdup(symbols[s]);
    }
  }

}


// <<load_expr_str>>
void load_branch_expr_str(
	// OUTPUT
	int *nexpr, char ***expr,
	// IN-OUT
	int linearize_mass, char *filepath
) {
  FILE *file = fopen(filepath, "r");
	if (file == NULL) {
		perror("Error opening file");
		exit(EXIT_FAILURE);
  }

	// // count the number of lines in the file
	// *nexpr = 0;
	// char *line = NULL;
	// size_t line_len = 0;
	// while (getline(&line, &line_len, file) != -1) {
  //   if (line[0] != '\n' && line[0] != '\0')
	// 	(*nexpr)++;
	// }
	// rewind(file);

  if (linearize_mass) {
    read_param_int(nexpr, (char*)"tot-branch", file);
  } else {
    read_param_int(nexpr, (char*)"massless-branch", file);
  }
  // cout << "nexpr = " << *nexpr << endl;

	// allocate memory for the symbols array
	*expr = (char**) malloc((*nexpr)*sizeof(char*));

	// // read expr from file
	// char *line = NULL;
	// size_t line_len = 0;
  // // goto_key(&line, "")

	// for (int i=0; i<*nexpr; i++) {
	// 	(*expr)[i] = NULL;
	// 	getline(&(*expr)[i], &line_len, file);
	// 	(*expr)[i][strlen((*expr)[i])-1] = '\0';
	// }

  read_list_str(expr, nexpr, (char*)"branches", file);

  fclose(file);
}


void bound_behav_from_file(
  char *filepath,
  int ***bound_behav, int ***mi_eig, int **mi_eig_num, int *num_region_classes,
  int dim
) {
  FILE *fptr = fopen(filepath, "r");
  if (!read_param_int(num_region_classes, (char*)"num-eig", fptr)) {
    perror("option num-eig not specified in bound_behav file");
    exit(1);
  };
  cout << "number of regions = " << *num_region_classes << endl;
  malloc_rk2_tens(*bound_behav, dim, *num_region_classes);
  malloc_rk2_tens(*mi_eig, dim, *num_region_classes);
  *mi_eig_num = new int[dim];

  nlist *BB = new nlist[dim];
  char MI_key[MAX_SYM_LEN];
  for (int m=0; m<dim; m++) {    
    snprintf(MI_key, sizeof(MI_key), "MI%d", m);
    nlist_read_file(&BB[m], MI_key, fptr);
    (*mi_eig_num)[m] = BB[m].nitems;
    for (int n=0; n<BB[m].nitems; n++) {
      for (int k=0; k<BB[m].item[n].nitems; k++) {
        if (strcmp(BB[m].item[n].item[k].item[0].str, "eig") == 0) {
          (*mi_eig)[m][n] = atoi(BB[m].item[n].item[k].item[1].str);
        }
      }
      for (int k=0; k<BB[m].item[n].nitems; k++) {
        if (strcmp(BB[m].item[n].item[k].item[0].str, "pow") == 0) {
          (*bound_behav)[m][(*mi_eig)[m][n]] = atoi(BB[m].item[n].item[k].item[1].str);
        }
      }
      cout << "m, np, behav = " << m << ", " << (*mi_eig)[m][n] << ", " << (*bound_behav)[m][(*mi_eig)[m][n]] << endl;
    }
  }

  // FREE
  nlist_rk1_free(BB, dim);
  
}


void bound_behav_from_file_old(
  char *file_name,
  int ***pow_behav, int ***mi_eig, int **mi_eig_num, int *num_classes,
  int num_mi
) {
  FILE *fptr = fopen(file_name, "r");
  
  // check if file exists
  if (fptr == NULL) {
    printf("Could not open file %s\n", file_name);
    exit(1);
  }

  // read number of regions
  fscanf(fptr, "%d\n", num_classes);
  fscanf(fptr, "\n");

  // alloc memory for boundary behavior
  malloc_rk2_tens(*pow_behav, num_mi, *num_classes);
  malloc_rk2_tens(*mi_eig, num_mi, *num_classes);
  *mi_eig_num = new int[num_mi];

  int count;
  long int save_pos;
  char c, *line;
  int num_eig;
  for (int m=0; m<num_mi; m++) {
    fscanf(fptr, "%d\n", &(*mi_eig_num)[m]);
    for (int n=0; n<(*mi_eig_num)[m]; n++) {
      fscanf(fptr, "%d\n", &(*mi_eig)[m][n]);
    }
    for (int lam, n=0; n<(*mi_eig_num)[m]; n++) {
      lam = (*mi_eig)[m][n];
      fscanf(fptr, "%d\n", &(*pow_behav)[m][lam]);
    }
  }

  fclose(fptr);
}


void qsort_singular_labels(int *label, int first, int last, complex_d *poles, int *sort_also_these){
  int i, j, pivot, temp;
  if(first<last){
    pivot = first;
    i = first;
    j = last;
    while(i<j){
        while(poles[label[i]].real()<=poles[label[pivot]].real() && i<last) {
          i++;
        }
        while(poles[label[j]].real()>poles[label[pivot]].real()) {
          j--;
        }
        if(i<j){
          temp=label[i];
          label[i]=label[j];
          label[j]=temp;
          // sort also these
          temp=sort_also_these[i];
          sort_also_these[i]=sort_also_these[j];
          sort_also_these[j]=temp;          
        }
    }
    temp = label[pivot];
    label[pivot] = label[j];
    label[j] = temp;
    // sort also these
    temp = sort_also_these[pivot];
    sort_also_these[pivot] = sort_also_these[j];
    sort_also_these[j] = temp;
    qsort_singular_labels(label, first, j-1, poles, sort_also_these);
    qsort_singular_labels(label, j+1, last, poles, sort_also_these);
  }
}


void qsort_singular_labels_mp(int *label, int first, int last, mpc_t *poles, int *sort_also_these){
  int i, j, pivot, temp;
  if(first<last){
    pivot = first;
    i = first;
    j = last;
    while(i<j){
        while(
          (
            mpfr_equal_within_tol(mpc_realref(poles[label[i]]), mpc_realref(poles[label[pivot]])) ||\
            mpfr_cmp(mpc_realref(poles[label[i]]), mpc_realref(poles[label[pivot]])) < 0
          ) &&\
          i<last
        ) {
          i++;
        }
        while(
          mpfr_cmp(mpc_realref(poles[label[j]]), mpc_realref(poles[label[pivot]])) > 0 &&\
          !mpfr_equal_within_tol(mpc_realref(poles[label[j]]), mpc_realref(poles[label[pivot]]))
        ) {
          j--;
        }
        if(i<j){
          temp=label[i];
          label[i]=label[j];
          label[j]=temp;
          // sort also these
          temp=sort_also_these[i];
          sort_also_these[i]=sort_also_these[j];
          sort_also_these[j]=temp;          
        }
    }
    temp = label[pivot];
    label[pivot] = label[j];
    label[j] = temp;
    // sort also these
    temp = sort_also_these[pivot];
    sort_also_these[pivot] = sort_also_these[j];
    sort_also_these[j] = temp;
    qsort_singular_labels_mp(label, first, j-1, poles, sort_also_these);
    qsort_singular_labels_mp(label, j+1, last, poles, sort_also_these);
  }
}


void generate_coordinates(
  // OUTPUT
  double *pts,
  int *npts,
  // INPUT
  complex_d *poles, int npoles,
  int RunRadius, int len_max, int mall
) {
  /*
  Given a list of poles in the complex plain, generate a set of points
  lying in the interval [0, 1] according to the following procedure:
    
    1. the 1st point is 0;
    2. consider the closest pole and move a step forward of a size
      given by the distance from such a pole divided by RunRadius;
    3. if the arrival point falls inside [0, 1] repeat from step 2,
      otherwise take 1 as final point.
  
  OUTPUT:
    - pts: array of output points;
    - *npts: number of output points;
  INPUT:
    - poles: list of output poles;
    - npoles: list of input poles;
    - RunRadius: parameter;
    - len_max: max number of generated points;
    - mall: if true, allocate memory for the ouput array of points and
      considering the exact number of generated points.
      If false, then use the input vector without allocating memory but
      still generating at most `len_max` points.
  */

  double *int_pts;
  if (mall) {
    int_pts = new double[len_max];
  } else {
    int_pts = pts;
  }

  double current = 0;
  double dists[npoles];

  int_pts[0] = current;
  *npts = 1;
  while (current < 1 && *npts < len_max) {
    for (int k=0; k<npoles; k++) {
      dists[k] = abs(current - poles[k]);
    }
    current = min(1., current + min1(dists, npoles)/RunRadius);
    // cout << "current = " << current << endl;
    int_pts[*npts] = current;
    (*npts)++;
  }

  if (mall) {
    pts = new double[*npts];
    for (int l=0; l<*npts; l++) {
      pts[l] = int_pts[l];
    }
    delete[] int_pts;
  }
}


void generate_coordinates_mp(
  // OUTPUT
  mpfr_t *pts,
  int *npts,
  // INPUT
  mpc_t *poles, int npoles,
  int RunRadius, int len_max, int mall
) {
  mpfr_t mpfr_one; mpfr_init2(mpfr_one, wp2);
  mpfr_set_ui(mpfr_one, 1, MPFR_RNDN);
  mpc_t tmpc; mpc_init3(tmpc, wp2, wp2);
  mpfr_t min_dist; mpfr_init2(min_dist, wp2);

  mpfr_t *int_pts;
  if (mall) {
    int_pts = new mpfr_t[len_max];
    init_rk1_mpfr(int_pts, len_max);
  } else {
    int_pts = pts;
  }

  mpfr_t current;
  mpfr_init2(current, wp2);
  mpfr_set_ui(current, 0, MPFR_RNDN);

  mpfr_t dists[npoles];
  init_rk1_mpfr(dists, npoles);

  mpfr_set(int_pts[0], current, MPFR_RNDN);
  *npts = 1;
  while (
    mpfr_cmp_si(current, 1) < 0 && !mpfr_equal_within_tol(current, mpfr_one) &&\
    *npts < len_max
  ) {
    for (int k = 0; k < npoles; k++) {
      mpc_sub_fr(tmpc, poles[k], current, MPFR_RNDN);
      mpc_abs(dists[k], tmpc, MPFR_RNDN);
    }
    min1_mp(&min_dist, dists, npoles);
    mpfr_div_si(min_dist, min_dist, RunRadius, MPFR_RNDN);
    mpfr_add(current, current, min_dist, MPFR_RNDN);
    if (mpfr_cmp_ui(current, 1) < 0) {
      mpfr_set(int_pts[*npts], current, MPFR_RNDN);
    } else {
      mpfr_set_ui(int_pts[*npts], 1, MPFR_RNDN);
    }
    (*npts)++;
  }

  if (mall) {
    pts = new mpfr_t[*npts];
    for (int l = 0; l < *npts; l++) {
      mpfr_init2(pts[l], wp2);
      mpfr_set(pts[l], int_pts[l], MPFR_RNDN);
    }
    for (int i = 0; i < len_max; i++) {
      mpfr_clear(int_pts[i]);
    }
    delete[] int_pts;
  }

  mpfr_clear(mpfr_one);
  mpc_clear(tmpc);
  mpfr_clear(min_dist);
  mpfr_clear(current);
  for (int i = 0; i < npoles; i++) {
    mpfr_clear(dists[i]);
  }
}


void generate_regular_points(
  // OUTPUT
  double *pts, int *npts,
  // INPUT
  complex_d *poles, int npoles,
  int RunRadius, int len_max,
  complex_d pt1, complex_d pt2
) {
  /*
  Given a list of poles in the complex plain, generate a set of points
  lying on the segment between two input complex points,
  according to the following procedure:
    
    1. the 1st point is the 1st end of the segment;
    2. consider the closest pole and move a step forward of a size
      given by the distance from such a pole divided by RunRadius;
    3. if the arrival point falls inside the segment, repeat from
      step 2, otherwise discard it and take the 2nd end of the segment
      as final point.
  
  The output are the coordinates in [0, 1] of the generated points
  along the segment.
  
  OUTPUT:
    - pts: array of output points;
    - *npts: number of output points;
  INPUT:
    - poles: list of output poles;
    - npoles: list of input poles;
    - RunRadius: parameter;
    - len_max: max number of generated points;
    - pt1, pt2: 1st and 2nd end of the segment in the complex plane;
  */
 
  complex_d mod_poles[npoles];
  complex_d dist = pt2 - pt1;
  for (int i=0; i<npoles; i++) {
    mod_poles[i] = (poles[i]-pt1)/dist;
    // cout << "mod_poles[" << i << "] = " << mod_poles[i] << endl;
  }

  generate_coordinates(
    pts, npts,
    mod_poles, npoles,
    RunRadius, len_max, 0
  );
  // for (int i=0; i<*npts; i++) {
  //   cout << "pts[" << i << "] = " << pts[i] << endl;
  // }

}


void generate_regular_points_mp(
  // OUTPUT
  mpfr_t *pts, int *npts,
  // INPUT
  mpc_t *poles, int npoles,
  int RunRadius, int len_max,
  mpc_t pt1, mpc_t pt2
) { 
  mpc_t mod_poles[npoles];
  mpc_t dist;
  mpc_init3(dist, wp2, wp2);
  mpc_sub(dist, pt2, pt1, MPFR_RNDN);
  for (int i = 0; i < npoles; i++) {
    mpc_init3(mod_poles[i], wp2, wp2);
    mpc_sub(mod_poles[i], poles[i], pt1, MPFR_RNDN);
    mpc_div(mod_poles[i], mod_poles[i], dist, MPFR_RNDN);
    // cout << "mod_poles[" << i << "] = " << mod_poles[i] << endl;
  }

  generate_coordinates_mp(
    pts, npts,
    mod_poles, npoles,
    RunRadius, len_max, 0
  );
  // for (int i = 0; i < *npts; i++) {
  //   cout << "pts[" << i << "] = "; print_mpfr(&pts[i]); cout << endl;
  // }

  mpc_clear(dist);
  for (int i = 0; i < npoles; i++) {
    mpc_clear(mod_poles[i]);
  }
}


void get_path_PS(
  // OUTPUT
  mpc_t **path, int **path_tags, int *neta_vals, int **sing_lab, int *nsings,
  // INPUT
  mpc_t *roots, int nroots, int zero_label
) {
  /*
  (DESCRIPTION IS DEPRECATED)
  Given the list of poles in the complex plane, generate a path of
  points lying on the real interval [0, 1] including all the poles
  that are in [0, 1] (singular points) and a certain number of points
  between them (regular and matching points). The points adjacent to
  the singular ones are matching points, while all the other are regular
  points. Between two singular points, at least one matching point has
  to be generated: if the disposition of the poles do not require its
  presence, it is manually placed at half the distance between the
  singular points.

  The matching points around a singular point are decided before the
  regular points, by looking at the closest pole: for every singular
  point, two matching points are placed right before and right after it
  at a distance given by the distance of the closest pole divided
  by a parameter `RunRadius`. Then, regular points are generated
  between the matching points lying immediatly after and immediatly
  before two consecutive singular points. The procedure is summarised
  in the following scheme:

    1. s, ---------------, s, ---------------, s, ---------------, s

    2. s, m, ---------, m, s, m, ---------, m, s, m, ---------, m, s

    3. s, m, r, ..., r, m, s, m, r, ..., r, m, s, m, r, ..., r, m, s
  
  The algorithm that generates the regular points between two matching
  points is the following:

    1. start from the first matching point;
    2. consider the closest pole and move a step forward of a size
      given by the distance from such a pole divided by `RunRadius`;
    3. if the arrival point falls before the 2nd matching point repeat
      from step 2, otherwise discard it and stop.

  */

  int print = 0;

  if (zero_label == -1) {
    roots++;
    nroots--;
  }

  // define parameters
  int RunRadius = 2, len_max = 200;

  // CONVERT POLES TO DOUBLE AND IDENTIFY SINGULAR POINTS
  // int npoles = nroots;
  int npoles = nroots;
  // double realp;
  complex_d *poles = new complex_d[npoles];
  int segs[npoles+2], nsegs = 1;
  int *tmp_sing_lab = new int[npoles];
  int count = 0;
  *nsings = 0;
  if (print) cout << "poles:" << endl;
  for (int k=0; k<nroots; k++) {
    if (print) {
    cout << "k = " << k << endl;
    cout << "mpc root:" << endl; print_mpc(&roots[k]); cout << endl;
    }
    // cout << "real and imaginary doubles:" << endl;
    // cout << mpfr_get_d(mpc_realref(roots[k]), MPFR_RNDN) << endl;
    // cout << mpfr_get_d(mpc_imagref(roots[k]), MPFR_RNDN) << endl;
    // realp = mpfr_get_d(mpc_realref(roots[k]), MPFR_RNDN);
    // if (realp == 0 || realp == 1) {
    //   continue;
    // }
    if (mpc_lessthan_tol(roots[k])) {
      poles[count] = (complex_d) {0,0};
    } else {
      poles[count] = (complex_d) {
        mpfr_get_d(mpc_realref(roots[k]), MPFR_RNDN),
        mpfr_get_d(mpc_imagref(roots[k]), MPFR_RNDN)
      };
    }
    // cout << poles[count] << endl;

    // IDENTIFY SINGULAR POINTS
    if (
      poles[count].imag() == 0 &&\
      0 <= poles[count].real() && poles[count].real() <= 1
    ) {
      // cout << poles[count] << endl;
      segs[(*nsings)] = count;
      tmp_sing_lab[(*nsings)] = k;
      (*nsings)++;
      if (poles[count].real() > 0) {
        nsegs++;
      }
    }
    count++;
  }
  if (print) {
  cout << endl;
  cout << "nsegs = " << nsegs << endl;
  cout << "nsings = " << *nsings << endl;
  }

  // arrange output singular labels
  *sing_lab = new int[*nsings];
  // *nsings = nsegs-1;
  for (int k=0; k<*nsings; k++) {
    (*sing_lab)[k] = tmp_sing_lab[k];    
  }

  // cout << "unsorted singular points:" << endl;
  // for (int k=0; k<*nsings; k++) {
  //   cout << poles[segs[k]] << endl;
  // }
  // cout << endl;

  // sort singular points in ascending order of real part
  if (*nsings >= 2) {
    qsort_singular_labels(segs, 0, (*nsings)-1, poles, *sing_lab);
  }

  // cout << "sorted singular points:" << endl;
  // for (int k=0; k<*nsings; k++) {
  //   cout << poles[segs[k]] << endl;
  // }
  // cout << endl;
  // cout << "sorted output sin lab:" << endl;
  // for (int k=0; k<*nsings; k++) {
  //   cout << (*sing_lab)[k] << endl;
  // }

  // COMPUTE DISTANCE OF MATCHING POINTS
  if (print) cout << "COMPUTE DISTANCE OF MATCHING POINTS" << endl;
  double *match = new double[nsegs];
  double dists[npoles];
  for (int s=0; s<*nsings; s++) {
    // calculate distance of the poles
    for (int k=0; k<npoles; k++) {
      dists[k] = abs(poles[k] - poles[segs[s]]);
    }
    match[s] = min1(dists, npoles)/RunRadius;
    if (print) cout << "match[" << s << "] = " << match[s] << endl;
  }
  
  // GENERATE MATCHING AND REGULAR POINTS
  double m1, m2;
  complex_d **int_points;
  int_points = new complex_d*[nsegs+1];
  int **tags;
  tags = new int*[nsegs+1];
  int *nints;
  nints = new int[nsegs+1];
  double **coordinates;
  coordinates = new double*[nsegs+1];
  int ncoords;
  complex_d pt1, pt2;
  *neta_vals = nsegs-1;

  int s_start = 0, s_end = 0;
  if (*nsings == 0) {
    s_end = 0;
    int_points[0] = new complex_d[len_max];
    tags[0] = new int[len_max];
    coordinates[0] = new double[len_max];
    tags[0][0] = 1;
    int_points[0][ncoords-1] = (complex_d) {0, 0};
    generate_regular_points(
      coordinates[0], &ncoords,
      poles, npoles,
      RunRadius, len_max,
      (complex_d) {0, 0}, (complex_d) {1, 0}
    );
    if (print) cout << "ncoords = " << ncoords << endl;
    for (int i=1; i<ncoords-1; i++) {
      if (print) {cout << "i = " << i << ", coordinate = " << coordinates[0][i] << endl;}
      int_points[0][i] = (complex_d) {0, 0} +\
        coordinates[0][i]*((complex_d) {1, 0}-(complex_d) {0, 0});
      tags[0][i] = 2;
    }
    // count internal points
    nints[0] = ncoords;
    int_points[0][ncoords-1] = (complex_d) {1, 0};
    tags[0][ncoords-1] = 1;
    *neta_vals += nints[0];
    goto merge_segments;
  }

  // deal with 1st segment if necessary
  if (poles[segs[0]].real() > 0) {
    if (print) cout << "processing 1st segment" << endl;
    int_points[0] = new complex_d[len_max];
    tags[0] = new int[len_max];
    coordinates[0] = new double[len_max];
    //// add 1st matching point
    int_points[0][0] = (complex_d) {0, 0};
    tags[0][0] = 1;
    m2 = poles[segs[0]].real() - match[0];
    // cout << "m2 = " << m2 << endl;
    if (m2 <= 0) {
      nints[0] = 1;
    } else {
      // generate regular points
      generate_regular_points(
        coordinates[0], &ncoords,
        poles, npoles,
        RunRadius, len_max,
        (complex_d) {0, 0}, (complex_d) {m2, 0}
      );
      for (int i=1; i<ncoords-1; i++) {
        // cout << "coordinate = " << coordinates[0][i] << endl;
        int_points[0][i] = (complex_d) {0, 0} +\
          coordinates[0][i]*((complex_d) {m2, 0}-(complex_d) {0, 0});
        tags[0][i] = 2;
      }
      // add 2nd matching point
      int_points[0][ncoords-1] = (complex_d) {m2, 0};
      tags[0][ncoords-1] = 1;
      // count internal points
      nints[0] = ncoords;
    }
    *neta_vals += nints[0];
    s_start = 0;
    s_end = nsegs - 1;
  } else {
    if (print) cout << "1st segment is trivial" << endl;
    s_start = 1;
    s_end = nsegs;
    (*neta_vals)++;
  }
  // cout << "1st segment:" << endl;
  // cout << "nints = " << nints[0] << endl;
  // for (int i=0; i<nints[0]; i++) {
  //   cout << tags[0][i] << ", " << int_points[0][i] << endl;
  // }

  if (print) {
  cout << "s_start = " << s_start << endl;
  cout << "s_end = " << s_end << endl;
  }
  
  // deal with internal segments
  if (print) cout << "deal with internal segments" << endl;
  for (int s=1; s<s_end; s++) {
    if (print) cout << "segment " << s << endl;
    // getchar();
    int_points[s] = new complex_d[len_max];
    tags[s] = new int[len_max];
    coordinates[s] = new double[len_max];

    // get position of the matching points
    m1 = poles[segs[s-1]].real() + match[s-1];
    m2 = poles[segs[s]].real() - match[s];
    if (print) cout << "m1 = " << m1 << endl;
    if (print) cout << "m2 = " << m2 << endl;
    
    if (m1 != poles[segs[s]].real() && m2 != poles[segs[s-1]].real()) {
      if (print) cout << "case 1" << endl;
      if (m1 >= m2) {
        //////
        // s1, m2, m1, s2
        //////
        if (print) cout << "s1, m2, m1, s2" << endl;

        // // add only 2nd matching point
        // int_points[s][0] = (complex_d) {m2, 0};
        // tags[s][0] = 1;
        // nints[s] = 1;
        // add middle point as matching point
        // int_points[s][0] = (poles[segs[s-1]] + poles[segs[s]])/2.;
        int_points[s][0] = (complex_d) {(m1+m2)/2., 0};
        tags[s][0] = 1;
        nints[s] = 1;
      } else {
        //////
        // s1, m1, ..., m2, s2
        //////
        if (print) cout << "s1, m1, ..., m2, s2" << endl;

        // add 1st matching point
        int_points[s][0] = (complex_d) {m1, 0};
        tags[s][0] = 1;
        // generate regular points
        generate_regular_points(
          coordinates[s], &ncoords,
          poles, npoles,
          RunRadius, len_max,
          (complex_d) {m1, 0}, (complex_d) {m2, 0}
        );
        for (int i=1; i<ncoords-1; i++) {
          // cout << "coordinate = " << coordinates[s][i] << endl;
          int_points[s][i] = (complex_d) {m1, 0} +\
            coordinates[s][i]*((complex_d) {m2, 0}-(complex_d) {m1, 0});
          tags[s][i] = 2;
        }
        // add 2nd matching point
        int_points[s][ncoords-1] = (complex_d) {m2, 0};
        tags[s][ncoords-1] = 1;
        // count internal points
        nints[s] = ncoords;
      }
    } else if (m1 != poles[segs[s]].real() && m2 == poles[segs[s-1]].real()) {
      cout << "case 2" << endl;
      //////
      // s1 = m2, m1, s2
      //////
      if (print) cout << "s1 = m2, m1, s2" << endl;

      // add 1st matching point
      int_points[s][0] = (complex_d) {m1, 0};
      tags[s][0] = 1;
      nints[s] = 1;
      // // generate regular points
      // generate_regular_points(
      //   coordinates[s], &ncoords,
      //   poles, npoles,
      //   RunRadius, len_max,
      //   int_points[s][0], poles[sing[s+1]]
      // );
      // for (int i=1; i<ncoords-1; i++) {
      //   int_points[s][i] = int_points[s][0] +\
      //     coordinates[s][i]*(poles[sing[s+1]]-(complex_d) {m1, 0});
      //   if (i == ncoords-1) {
      //     tags[s][i] = 1;
      //   } else {
      //     tags[s][i] = 2;
      //   }
      // }
      // // count internal points
      // nints[s] = ncoords - 1;
    // } else if (m1 != poles[segs[s]].real() && m2 == poles[segs[s-1]].real()) {
    } else if (m1 == poles[segs[s]].real() && m2 != poles[segs[s-1]].real()) {
      //////
      // s1, m2, m1 = s2
      //////
      if (print) cout << "s1, m2, m1 = s2" << endl;
      
      // add 1nd matching point
      // int_points[s][0] = (complex_d) {m1, 0};
      int_points[s][0] = (complex_d) {m2, 0};
      tags[s][0] = 1;
      nints[s] = 1;
    } else {
      if (print) cout << "else" << endl;
      // add middle point as matching point
      int_points[s][0] = (poles[segs[s-1]] + poles[segs[s]])/2.;
      tags[s][0] = 1;
      nints[s] = 1;
    }
    *neta_vals += nints[s];
  }

  // deal with last segments if necessary
  // (i.e. only when last point is not singular)
  if (poles[segs[s_end-1]].real() < 1) {
    if (print) cout << "processing last segment" << endl;
    int_points[s_end] = new complex_d[len_max];
    tags[s_end] = new int[len_max];
    coordinates[s_end] = new double[len_max];
    
    m1 = poles[segs[s_end-1]].real() + match[s_end-1];
    if (print) {
    cout << "last segment" << endl;
    cout << "nsegs = " << nsegs << endl;
    cout << "index of last pole = " << segs[s_end-1] << endl;
    cout << "last pole = " << poles[segs[s_end-1]] << endl;
    cout << "m1 = " << m1 << endl;
    }
    // getchar();
    if (m1 >= 1) {
      nints[s_end] = 1;
      ncoords = 1;
    } else {
      // add 1st matching point
      int_points[s_end][0] = (complex_d) {m1, 0};
      tags[s_end][0] = 1;
      // generate regular points
      generate_regular_points(
        coordinates[s_end], &ncoords,
        poles, npoles,
        RunRadius, len_max,
        (complex_d) {m1, 0}, (complex_d) {1, 0}
      );
      for (int i=1; i<ncoords-1; i++) {
        // cout << "coordinate = " << coordinates[s_end][i] << endl;
        int_points[s_end][i] = (complex_d) {m1, 0} +\
          coordinates[s_end][i]*((complex_d) {1, 0}-(complex_d) {m1, 0});
        tags[s_end][i] = 2;
      }
      // count internal points
      nints[s_end] = ncoords;
    }
    int_points[s_end][ncoords-1] = (complex_d) {1, 0};
    tags[s_end][ncoords-1] = 1;
    *neta_vals += nints[s_end];

    if (print) {
    cout << "last segment:" << endl;
    cout << "nints = " << nints[s_end] << endl;
    for (int i=0; i<nints[s_end]; i++) {
      cout << tags[s_end][i] << ", " << int_points[s_end][i] << endl;
    }
    }
  } else {
    if (print) cout << "last segment is trivial" << endl;
    nints[s_end] = 0;
  }

  merge_segments:
  // MERGE SEGMENTS IN THE OUPUT
  if (print) cout << "MERGING SEGMENTS IN PATH" << endl;
  if (print) cout << *neta_vals << " eta values" << endl;
  *path = new mpc_t[*neta_vals];
  init_rk1_mpc(*path, *neta_vals);
  *path_tags = new int[*neta_vals];
  int offset = 0;
  if (print) cout << "n. segments = " << nsegs << endl;
  for (int s=s_start; s<s_end+1; s++) {
  // for (int s=0; s<s_end+1; s++) {
    if (print) cout << "s = " << s << endl;
    if (s != 0) {
      mpc_set_d_d((*path)[offset], poles[segs[s-1]].real(), poles[segs[s-1]].imag(), MPFR_RNDN);
      (*path_tags)[offset] = 0;
      if (print) {
      cout << "sing point:" << endl;
      cout << offset << ": ";
      cout << (*path_tags)[offset] << ", "; print_mpc(&((*path)[offset])); cout << endl;
      }
      offset++;
    }
    if (print) cout << nints[s] << " internal points:" << endl;
    for (int i=0; i<nints[s]; i++) {
      mpc_set_d_d(
        (*path)[offset+i], int_points[s][i].real(), int_points[s][i].imag(),
        MPFR_RNDN
      );
      (*path_tags)[offset+i] = tags[s][i];
      if (print) cout << offset+i << ": ";
      if (print) {cout << (*path_tags)[offset+i] << ", "; print_mpc(&((*path)[offset+i])); cout << endl;}
      if (print) printf(" double: (%.27g %.27g)\n", int_points[s][i].real(), int_points[s][i].imag());
    }

    offset += nints[s];
  }

  if (zero_label == -1) {
    roots--;
    nroots++;
    for (int i=0; i<*nsings; i++) {
      (*sing_lab)[i]++;
    }
  }
  
}


void get_path_PS_mp(
  // OUTPUT
  mpc_t **path, int **path_tags, int *neta_vals, int **sing_lab, int *nsings,
  // INPUT
  mpc_t *roots, int nroots, int zero_label
) {
  int print = 0;

  mpc_t tmpc; mpc_init3(tmpc, wp2, wp2);
  mpfr_t tmpfr; mpfr_init2(tmpfr, wp2);
  mpc_t mpc_zero; mpc_init3(mpc_zero, wp2, wp2);
  mpc_set_ui(mpc_zero, 0, MPFR_RNDN);
  mpc_t mpc_one; mpc_init3(mpc_one, wp2, wp2);
  mpc_set_ui(mpc_one, 1, MPFR_RNDN);
  mpfr_t mpfr_one; mpfr_init2(mpfr_one, wp2);
  mpfr_set_ui(mpfr_one, 1, MPFR_RNDN);

  if (zero_label == -1) {
    roots++;
    nroots--;
  }

  // define parameters
  int RunRadius = 2, len_max = 200;

  // CONVERT POLES TO DOUBLE AND IDENTIFY SINGULAR POINTS
  // int npoles = nroots;
  int npoles = nroots;
  // double realp;
  mpc_t *poles = new mpc_t[npoles];
  init_rk1_mpc(poles, npoles);
  int segs[npoles+2], nsegs = 1;
  int *tmp_sing_lab = new int[npoles];
  int count = 0;
  *nsings = 0;
  if (print) cout << "poles:" << endl;
  for (int k=0; k<nroots; k++) {
    if (print) {
    cout << "k = " << k << endl;
    cout << "mpc root:" << endl; print_mpc(&roots[k]); cout << endl;
    }
    // cout << "real and imaginary doubles:" << endl;
    // cout << mpfr_get_d(mpc_realref(roots[k]), MPFR_RNDN) << endl;
    // cout << mpfr_get_d(mpc_imagref(roots[k]), MPFR_RNDN) << endl;
    // realp = mpfr_get_d(mpc_realref(roots[k]), MPFR_RNDN);
    // if (realp == 0 || realp == 1) {
    //   continue;
    // }
    if (mpc_lessthan_tol(roots[k])) {
      mpc_set_ui(poles[count], 0, MPFR_RNDN);
    } else {
      mpc_set(poles[count], roots[k], MPFR_RNDN);
    }
    // cout << poles[count] << endl;

    // IDENTIFY SINGULAR POINTS
    if (
      mpfr_lessthan_tol(mpc_imagref(poles[count])) &&\
      (
        mpfr_lessthan_tol(mpc_realref(poles[count])) ||\
        (
          mpfr_cmp_si(mpc_realref(poles[count]), 0) > 0 &&\
          mpfr_cmp_si(mpc_realref(poles[count]), 1) < 0
        ) ||\
        mpfr_equal_within_tol(mpc_realref(poles[count]), mpfr_one)
      )
    ) {
      // cout << poles[count] << endl;
      segs[(*nsings)] = count;
      tmp_sing_lab[(*nsings)] = k;
      (*nsings)++;
      if (
        !mpfr_lessthan_tol(mpc_realref(poles[count])) &&
        mpfr_cmp_ui(mpc_realref(poles[count]), 0) > 0
      ) {
        nsegs++;
      }
    }
    count++;
  }
  if (print) {
  cout << endl;
  cout << "nsegs = " << nsegs << endl;
  cout << "nsings = " << *nsings << endl;
  }

  // arrange output singular labels
  *sing_lab = new int[*nsings];
  // *nsings = nsegs-1;
  for (int k=0; k<*nsings; k++) {
    (*sing_lab)[k] = tmp_sing_lab[k];    
  }

  // cout << "unsorted singular points:" << endl;
  // for (int k=0; k<*nsings; k++) {
  //   cout << poles[segs[k]] << endl;
  // }
  // cout << endl;

  // sort singular points in ascending order of real part
  if (*nsings >= 2) {
    qsort_singular_labels_mp(segs, 0, (*nsings)-1, poles, *sing_lab);
  }

  // cout << "sorted singular points:" << endl;
  // for (int k=0; k<*nsings; k++) {
  //   cout << poles[segs[k]] << endl;
  // }
  // cout << endl;
  // cout << "sorted output sin lab:" << endl;
  // for (int k=0; k<*nsings; k++) {
  //   cout << (*sing_lab)[k] << endl;
  // }

  // COMPUTE DISTANCE OF MATCHING POINTS
  if (print) cout << "COMPUTE DISTANCE OF MATCHING POINTS" << endl;
  mpfr_t *match = new mpfr_t[nsegs];
  mpfr_t dists[npoles];
  for (int s=0; s<*nsings; s++) {
    // calculate distance of the poles
    for (int k=0; k<npoles; k++) {
      mpc_sub(tmpc, poles[k], poles[segs[s]], MPFR_RNDN);
      mpfr_init2(dists[k], wp2);
      mpc_abs(dists[k], tmpc, MPFR_RNDN);
    }
    mpfr_init2(match[s], wp2);
    min1_mp(&match[s], dists, npoles);
    mpfr_div_ui(match[s], match[s], RunRadius, MPFR_RNDN);
    if (print) {cout << "match[" << s << "] = "; print_mpfr(&match[s]); cout << endl;}
  }


  // GENERATE MATCHING AND REGULAR POINTS
  mpfr_t m1, m2;
  mpfr_init2(m1, wp2);
  mpfr_init2(m2, wp2);
  mpc_t **int_points;
  int_points = new mpc_t*[nsegs+1];
  int **tags;
  tags = new int*[nsegs+1];
  int *nints;
  nints = new int[nsegs+1];
  mpfr_t **coordinates;
  coordinates = new mpfr_t*[nsegs+1];
  int ncoords;
  mpc_t pt1, pt2;
  mpc_init3(pt1, wp2, wp2);
  mpc_init3(pt2, wp2, wp2);
  *neta_vals = nsegs-1;

  int s_start = 0, s_end = 0;
  if (*nsings == 0) {
    s_end = 0;
    int_points[0] = new mpc_t[len_max];
    init_rk1_mpc(int_points[0], len_max);
    tags[0] = new int[len_max];
    coordinates[0] = new mpfr_t[len_max];
    init_rk1_mpfr(coordinates[0], len_max);
    tags[0][0] = 1;
    mpc_set_ui(int_points[0][0], 0, MPFR_RNDN);
    generate_regular_points_mp(
      coordinates[0], &ncoords,
      poles, npoles,
      RunRadius, len_max,
      mpc_zero, mpc_one
    );
    if (print) cout << "ncoords = " << ncoords << endl;
    for (int i=1; i<ncoords-1; i++) {
      if (print) {cout << "i = " << i << ", coordinate = "; print_mpfr(&coordinates[0][i]); cout << endl;}
      mpc_set_fr(int_points[0][i], coordinates[0][i], MPFR_RNDN);
      tags[0][i] = 2;
    }
    // count internal points
    nints[0] = ncoords;
    mpc_set_ui(int_points[0][ncoords-1], 1, MPFR_RNDN);
    tags[0][ncoords-1] = 1;
    *neta_vals += nints[0];
    goto merge_segments;
  }

  // deal with 1st segment if necessary
  if (
    !mpfr_lessthan_tol(mpc_realref(poles[segs[0]])) &&\
    mpfr_cmp_ui(mpc_realref(poles[segs[0]]), 0) > 0
  ) {
    if (print) cout << "processing 1st segment" << endl;
    int_points[0] = new mpc_t[len_max];
    init_rk1_mpc(int_points[0], len_max);
    tags[0] = new int[len_max];
    coordinates[0] = new mpfr_t[len_max];
    init_rk1_mpfr(coordinates[0], len_max);
    //// add 1st matching point
    mpc_set_ui(int_points[0][0], 0, MPFR_RNDN);
    tags[0][0] = 1;
    mpfr_sub(m2, mpc_realref(poles[segs[0]]), match[0], MPFR_RNDN);
    // cout << "m2 = " << m2 << endl;
    if (
      mpfr_lessthan_tol(m2) ||\
      mpfr_cmp_ui(m2, 0) < 0
    ) {
      nints[0] = 1;
    } else {
      // generate regular points
      mpc_set_fr(tmpc, m2, MPFR_RNDN);
      generate_regular_points_mp(
        coordinates[0], &ncoords,
        poles, npoles,
        RunRadius, len_max,
        mpc_zero, tmpc
      );
      for (int i=1; i<ncoords-1; i++) {
        // cout << "coordinate = " << coordinates[0][i] << endl;
        mpc_mul_fr(int_points[0][i], tmpc, coordinates[0][i], MPFR_RNDN);
        tags[0][i] = 2;
      }
      // add 2nd matching point
      mpc_set_fr(int_points[0][ncoords-1], m2, MPFR_RNDN);
      tags[0][ncoords-1] = 1;
      // count internal points
      nints[0] = ncoords;
    }
    *neta_vals += nints[0];
    s_start = 0;
    s_end = nsegs - 1;
  } else {
    if (print) cout << "1st segment is trivial" << endl;
    s_start = 1;
    s_end = nsegs;
    (*neta_vals)++;
  }
  // cout << "1st segment:" << endl;
  // cout << "nints = " << nints[0] << endl;
  // for (int i=0; i<nints[0]; i++) {
  //   cout << tags[0][i] << ", " << int_points[0][i] << endl;
  // }

  if (print) {
  cout << "s_start = " << s_start << endl;
  cout << "s_end = " << s_end << endl;
  }

  // deal with internal segments
  if (print) cout << "deal with internal segments" << endl;
  for (int s=1; s<s_end; s++) {
    if (print) cout << "segment " << s << endl;
    // getchar();
    int_points[s] = new mpc_t[len_max];
    init_rk1_mpc(int_points[s], len_max);
    tags[s] = new int[len_max];
    coordinates[s] = new mpfr_t[len_max];
    init_rk1_mpfr(coordinates[s], len_max);

    // get position of the matching points
    mpfr_add(m1, mpc_realref(poles[segs[s-1]]), match[s-1], MPFR_RNDN);
    mpfr_sub(m2, mpc_realref(poles[segs[s]]), match[s], MPFR_RNDN);
    if (print) cout << "m1 = " << m1 << endl;
    if (print) cout << "m2 = " << m2 << endl;
    
    if (
      !mpfr_equal_within_tol(m1, mpc_realref(poles[segs[s]])) &&\
      !mpfr_equal_within_tol(m2, mpc_realref(poles[segs[s-1]]))
    ) {
      if (print) cout << "case 1" << endl;
      if (
        mpfr_equal_within_tol(m1, m2) ||\
        mpfr_cmp(m1, m2) > 0
      ) {
        //////
        // s1, m2, m1, s2
        //////
        if (print) cout << "s1, m2, m1, s2" << endl;

        // // add only 2nd matching point
        // int_points[s][0] = (complex_d) {m2, 0};
        // tags[s][0] = 1;
        // nints[s] = 1;
        // add middle point as matching point
        mpc_set_fr(int_points[s][0], m1, MPFR_RNDN);
        mpc_add_fr(int_points[s][0], int_points[s][0], m2, MPFR_RNDN);
        mpc_div_ui(int_points[s][0], int_points[s][0], 2, MPFR_RNDN);
        tags[s][0] = 1;
        nints[s] = 1;
      } else {
        //////
        // s1, m1, ..., m2, s2
        //////
        if (print) cout << "s1, m1, ..., m2, s2" << endl;

        // add 1st matching point
        mpc_set_fr(int_points[s][0], m1, MPFR_RNDN);
        tags[s][0] = 1;
        // generate regular points
        mpc_set_fr(tmpc, m2, MPFR_RNDN);
        generate_regular_points_mp(
          coordinates[s], &ncoords,
          poles, npoles,
          RunRadius, len_max,
          int_points[s][0], tmpc
        );
        for (int i=1; i<ncoords-1; i++) {
          // cout << "coordinate = " << coordinates[s][i] << endl;
          mpfr_sub(tmpfr, m2, m1, MPFR_RNDN);
          mpfr_fma(
            mpc_realref(int_points[s][i]),
            tmpfr, coordinates[s][i], m1,
            MPFR_RNDN
          );
          mpfr_set_ui(mpc_imagref(int_points[s][i]), 0, MPFR_RNDN);
          tags[s][i] = 2;
        }
        // add 2nd matching point
        mpc_set_fr(int_points[s][ncoords-1], m2, MPFR_RNDN);
        tags[s][ncoords-1] = 1;
        // count internal points
        nints[s] = ncoords;
      }
    } else if (
      !mpfr_equal_within_tol(m1, mpc_realref(poles[segs[s]])) &&\
      mpfr_equal_within_tol(m2, mpc_realref(poles[segs[s-1]]))
    ) {
      cout << "case 2" << endl;
      //////
      // s1 = m2, m1, s2
      //////
      if (print) cout << "s1 = m2, m1, s2" << endl;

      // add 1st matching point
      mpc_set_fr(int_points[s][0], m1, MPFR_RNDN);
      tags[s][0] = 1;
      nints[s] = 1;
    } else if (
      mpfr_equal_within_tol(m1, mpc_realref(poles[segs[s]])) &&\
      !mpfr_equal_within_tol(m2, mpc_realref(poles[segs[s-1]]))
    ) {
      // #adjust
      //////
      // s1, m2, m1 = s2
      //////
      if (print) cout << "s1, m2, m1 = s2" << endl;
      
      // add 1nd matching point
      mpc_set_fr(int_points[s][0], m2, MPFR_RNDN);
      tags[s][0] = 1;
      nints[s] = 1;
    } else {
      if (print) cout << "else" << endl;
      // add middle point as matching point
      mpc_add(int_points[s][0], poles[segs[s-1]], poles[segs[s]], MPFR_RNDN);
      mpc_div_ui(int_points[s][0], int_points[s][0], 2, MPFR_RNDN);
      tags[s][0] = 1;
      nints[s] = 1;
    }
    *neta_vals += nints[s];
  }

  // deal with last segments if necessary
  // (i.e. only when last point is not singular)
  if (
    !mpfr_equal_within_tol(mpc_realref(poles[segs[s_end-1]]), mpfr_one) &&\
    mpfr_cmp_ui(mpc_realref(poles[segs[s_end-1]]), 1) < 0
  ) {
    if (print) cout << "processing last segment" << endl;
    int_points[s_end] = new mpc_t[len_max];
    init_rk1_mpc(int_points[s_end], len_max);
    tags[s_end] = new int[len_max];
    coordinates[s_end] = new mpfr_t[len_max];
    init_rk1_mpfr(coordinates[s_end], len_max);

    mpfr_add(m1, mpc_realref(poles[segs[s_end-1]]), match[s_end-1], MPFR_RNDN);
    if (print) {
    cout << "last segment" << endl;
    cout << "nsegs = " << nsegs << endl;
    cout << "index of last pole = " << segs[s_end-1] << endl;
    cout << "last pole = "; print_mpc(&poles[segs[s_end-1]]); cout << endl;
    cout << "m1 = " << m1 << endl;
    }
    // getchar();
    if (
      mpfr_equal_within_tol(m1, mpfr_one) ||\
      mpfr_cmp_ui(m1, 1) > 0
    ) {
      nints[s_end] = 1;
      ncoords = 1;
    } else {
      // add 1st matching point
      mpc_set_fr(int_points[s_end][0], m1, MPFR_RNDN);
      tags[s_end][0] = 1;
      // generate regular points
      generate_regular_points_mp(
        coordinates[s_end], &ncoords,
        poles, npoles,
        RunRadius, len_max,
        int_points[s_end][0], mpc_one
      );
      for (int i=1; i<ncoords-1; i++) {
        // cout << "coordinate = " << coordinates[s_end][i] << endl;
        mpfr_sub(tmpfr, mpfr_one, m1, MPFR_RNDN);
        mpfr_fma(
          mpc_realref(int_points[s_end][i]),
          tmpfr, coordinates[s_end][i], m1,
          MPFR_RNDN
        );
        mpfr_set_ui(mpc_imagref(int_points[s_end][i]), 0, MPFR_RNDN);
        tags[s_end][i] = 2;
      }
      // count internal points
      nints[s_end] = ncoords;
    }
    mpc_set_ui(int_points[s_end][ncoords-1], 1, MPFR_RNDN);
    tags[s_end][ncoords-1] = 1;
    *neta_vals += nints[s_end];

    if (print) {
    cout << "last segment:" << endl;
    cout << "nints = " << nints[s_end] << endl;
    for (int i=0; i<nints[s_end]; i++) {
      cout << tags[s_end][i] << ", "; print_mpc(&int_points[s_end][i]); cout << endl;
    }
    }
  } else {
    if (print) cout << "last segment is trivial" << endl;
    nints[s_end] = 0;
  }

  merge_segments:
  // MERGE SEGMENTS IN THE OUPUT
  if (print) cout << "MERGING SEGMENTS IN PATH" << endl;
  if (print) cout << *neta_vals << " eta values" << endl;
  *path = new mpc_t[*neta_vals];
  init_rk1_mpc(*path, *neta_vals);
  *path_tags = new int[*neta_vals];
  int offset = 0;
  if (print) cout << "n. segments = " << nsegs << endl;
  for (int s=s_start; s<s_end+1; s++) {
  // for (int s=0; s<s_end+1; s++) {
    if (print) cout << "s = " << s << endl;
    if (s != 0) {
      mpc_set((*path)[offset], poles[segs[s-1]], MPFR_RNDN);
      (*path_tags)[offset] = 0;
      if (print) {
      cout << "sing point:" << endl;
      cout << offset << ": ";
      cout << (*path_tags)[offset] << ", "; print_mpc(&((*path)[offset])); cout << endl;
      }
      offset++;
    }
    if (print) cout << nints[s] << " internal points:" << endl;
    for (int i=0; i<nints[s]; i++) {
      mpc_set((*path)[offset+i], int_points[s][i], MPFR_RNDN);
      (*path_tags)[offset+i] = tags[s][i];
      if (print) cout << offset+i << ": ";
      if (print) {cout << (*path_tags)[offset+i] << ", "; print_mpc(&((*path)[offset+i])); cout << endl;}
      if (print) {cout << "  int_points = "; print_mpc(&int_points[s][i]); cout << endl;}
    }

    offset += nints[s];
  }

  if (zero_label == -1) {
    roots--;
    nroots++;
    for (int i=0; i<*nsings; i++) {
      (*sing_lab)[i]++;
    }
  }

  // FREE
  mpc_clear(tmpc);
  mpfr_clear(tmpfr);
  mpc_clear(mpc_zero);
  mpc_clear(mpc_one);
  mpfr_clear(mpfr_one);
  mpc_clear(pt1);
  mpc_clear(pt2);
  
}


void get_path(mpc_t **eta_values, int* neta_vals, complex_d *poles, int npoles) {
  int print = 0;

  // cout << "poles:" << endl;
  // for (int k=0; k<npoles; k++) {
  //   cout << poles[k] << endl;
  // }

  // define parameters
  int RunCandidate = 10, RunRadius = 2, len_max = 200;

  // compute initial and final point in the [0, 1] interval
  double ini, fin;
  if (print) cout << "npoles = " << npoles << endl;
  if (npoles == 0 && poles[0] == (complex_d){0, 0}) {
    ini = 1;
    fin = 1;
  } else {
    double pole_mods[npoles];
    for(int i=0; i<npoles; i++) {
      pole_mods[i] = abs(poles[i]);
    }
    double min, max;
    minmax0(&min, &max, pole_mods, npoles);
    ini = max*RunRadius;
    fin = min/RunRadius;
  }
  if (print) cout << "ini = " << ini << endl;
  if (print) cout << "fin = " << fin << endl;

  // define candidate directions
  complex_d dirs[2*RunCandidate+1];
  if (print) cout << "candidates" << endl;
  for (int i=-RunCandidate; i<=RunCandidate; i++) {
    dirs[i+RunCandidate] = ((complex_d) {(double) i, -1.} )/ sqrt(1+i*i);
    if (print) cout << dirs[i+RunCandidate] << endl;
  }
  
  // loop that generates points in the [0, 1] interval
  int run_len = len_max, fails, count_run, count_run_sav;
  double run[len_max], current, dists[npoles];
  double* run_pts;
  complex_d cen_rot_poles[npoles], dir_sav;
  if (ini==fin)
    *run_pts = ini;
  else {
    // // #debug: force to select 1st direction
    // for (int i=0; i<1; i++) {
    // i = 10;
    for (int i=0; i<=2*RunCandidate; i++) {
      if (print) cout << "candidate " << i << ":" << endl;
      fails = 0;
      for (int j=0; j<npoles; j++) {
        cen_rot_poles[j] = (poles[j]/dirs[i] -  ini)/(fin - ini);
        if (print) cout << "cen_rot_poles[" << j << "] = " << cen_rot_poles[j] << endl;
        if (cen_rot_poles[j].imag() == 0 &&  0 <= cen_rot_poles[j].real() && cen_rot_poles[j].real() <= 1) {
          if (print) cout << "fails" << endl;
          fails = 1;
          break;
        }
      }
      if (fails) {
        run[1] = {-1};
      } else {
        if (print) cout << "pass" << endl;
        current = 0;
        run[0] = current;
        count_run = 1;
        while (current < 1 && count_run < len_max) {
          for (int k=0; k<npoles; k++) {
            dists[k] = abs(current - cen_rot_poles[k]);
          }
          current = min(1., current + min1(dists, npoles)/RunRadius);
          run[count_run] = current;
          count_run++;
        }
        if (print) cout << "count run = " << count_run << endl;
      }
      if (count_run < run_len && run[0] != -1) {
        run_len = count_run;
        run_pts = new double[count_run];
        for (int l=0; l<count_run; l++)
          run_pts[l] = run[l];
        dir_sav = dirs[i];
        count_run_sav = count_run;
      }
    }
  }
  if (print) cout << "dir = " << dir_sav << endl;
  if (print) cout << "count run = " << count_run << endl;

  // map back [0, 1] to complex half-line
  if (print) cout << "map back [0, 1] to complex half-line" << endl;
  *eta_values = new mpc_t[count_run_sav];
  complex_d tmp;
  for (int i=0; i<count_run_sav ; i++) {
    if (print) cout << "i = " << i << endl;
    if (print) cout << "run_pts = " << run_pts[i] << endl;
    tmp = (ini + run_pts[i]*(fin - ini))*dir_sav;
    mpc_init3((*eta_values)[i], wp2, wp2);
    mpc_set_d_d((*eta_values)[i], tmp.real(), tmp.imag(), MPFR_RNDN);
    // print_mpc(&(*eta_values[i])); cout << endl; fflush(stdout);
  }
  delete[] run_pts;
  
  // save number of path values
  *neta_vals = count_run_sav;

}


void get_path_mp(mpc_t **eta_values, int* neta_vals, mpc_t *poles, int npoles) {
  int print = 0;

  // define parameters
  int RunCandidate = 10, RunRadius = 2, len_max = 200;

  mpc_t tmpc; mpc_init3(tmpc, wp2, wp2);
  mpfr_t tmpfr; mpfr_init2(tmpfr, wp2);
  mpfr_t mpfr_one; mpfr_init2(mpfr_one, wp2);
  mpfr_set_ui(mpfr_one, 1, MPFR_RNDN);
  mpfr_t min_dist;
  mpfr_init2(min_dist, wp2);

  // compute initial and final point in the [0, 1] interval
  mpfr_t ini, fin;
  mpfr_init2(ini, wp2);
  mpfr_init2(fin, wp2);
  mpfr_t min, max;
  mpfr_init2(min, wp2);
  mpfr_init2(max, wp2);
  if (print) cout << "npoles = " << npoles << endl;
  if (npoles == 0 && mpc_cmp_si_si(poles[0], 0, 0) == 0) {
    mpfr_set_si(ini, 1, MPFR_RNDN);
    mpfr_set_si(fin, 1, MPFR_RNDN);
  } else {
    mpfr_t pole_mods[npoles];
    for(int i=0; i<npoles; i++) {
      mpfr_init2(pole_mods[i], wp2);
      mpc_abs(pole_mods[i], poles[i], MPFR_RNDN);
    }
    minmax0_mp(&min, &max, pole_mods, npoles);
    mpfr_mul_si(ini, max, RunRadius, MPFR_RNDN);
    mpfr_div_si(fin, min, RunRadius, MPFR_RNDN);
    for(int i=0; i<npoles; i++) {
      mpfr_clear(pole_mods[i]);
    }
  }
  if (print) {cout << "ini = "; print_mpfr(&ini); cout << endl;}
  if (print) {cout << "fin = "; print_mpfr(&fin); cout << endl;}

  // define candidate directions
  mpc_t dirs[2*RunCandidate+1];
  mpfr_t denom;
  mpfr_init2(denom, wp2);
  if (print) cout << "candidates" << endl;
  for (int i=-RunCandidate; i<=RunCandidate; i++) {
    mpc_init2(dirs[i+RunCandidate], wp2);
    mpc_set_d_d(dirs[i+RunCandidate], (double) i, -1., MPFR_RNDN);
    mpfr_set_si(denom, 1 + i * i, MPFR_RNDN);
    mpfr_sqrt(denom, denom, MPFR_RNDN);
    mpc_div_fr(dirs[i+RunCandidate], dirs[i+RunCandidate], denom, MPFR_RNDN);
    if (print) {print_mpc(&dirs[i+RunCandidate]); cout << endl;}
  }
  mpfr_clear(denom);
  
  // loop that generates points in the [0, 1] interval
  int run_len = len_max, fails, count_run, count_run_sav;
  mpfr_t run[len_max], current, dists[npoles], run_pts[len_max];
  mpc_t cen_rot_poles[npoles], dir_sav;
  for (int i=0; i<len_max; i++) {
    mpfr_init2(run[i], wp2);
    mpfr_init2(run_pts[i], wp2);
  }
  for (int i=0; i<npoles; i++) {
    mpfr_init2(dists[i], wp2);
    mpc_init3(cen_rot_poles[i], wp2, wp2);
  }
  mpc_init2(dir_sav, wp2);
  mpfr_init2(current, wp2);
  if (mpfr_equal_within_tol(ini, fin)) {
    mpfr_set(run_pts[0], ini, MPFR_RNDN);
  } else {
    for (int i=0; i<=2*RunCandidate; i++) {
      if (print) cout << "candidate " << i << ":" << endl;
      fails = 0;
      for (int j=0; j<npoles; j++) {
        mpc_div(cen_rot_poles[j], poles[j], dirs[i], MPFR_RNDN);
        mpc_sub_fr(cen_rot_poles[j], cen_rot_poles[j], ini, MPFR_RNDN);
        mpfr_sub(tmpfr, fin, ini, MPFR_RNDN);
        mpc_div_fr(cen_rot_poles[j], cen_rot_poles[j], tmpfr, MPFR_RNDN);
        if (print) {cout << "cen_rot_poles[" << j << "] = "; print_mpc(&cen_rot_poles[j]); cout << endl;}
        if (
          mpfr_lessthan_tol(mpc_imagref(cen_rot_poles[j])) &&\
          (
            mpfr_lessthan_tol(mpc_realref(cen_rot_poles[j])) ||\
            (
              mpfr_cmp_si(mpc_realref(cen_rot_poles[j]), 0) > 0 &&\
              mpfr_cmp_si(mpc_realref(cen_rot_poles[j]), 1) < 0
            ) ||\
            mpfr_equal_within_tol(mpc_realref(cen_rot_poles[j]), mpfr_one)
          )
        ) {
          if (print) cout << "fails" << endl;
          fails = 1;
          break;
        }
      }
      if (fails) {
        mpfr_set_si(run[1], -1, MPFR_RNDN);
      } else {
        if (print) cout << "pass" << endl;
        mpfr_set_si(current, 0, MPFR_RNDN);
        mpfr_set(run[0], current, MPFR_RNDN);
        count_run = 1;
        while (
          (mpfr_cmp_si(current, 1) < 0 && !mpfr_equal_within_tol(current, mpfr_one)) &&\
          count_run < len_max
        ) {
          for (int k=0; k<npoles; k++) {
            mpfr_sub(mpc_realref(tmpc), current, mpc_realref(cen_rot_poles[k]), MPFR_RNDN);
            mpfr_neg(mpc_imagref(tmpc), mpc_imagref(cen_rot_poles[k]), MPFR_RNDN);
            mpc_abs(dists[k], tmpc, MPFR_RNDN);
          }
          min1_mp(&min_dist, dists, npoles);
          mpfr_div_si(min_dist, min_dist, RunRadius, MPFR_RNDN);
          mpfr_add(current, current, min_dist, MPFR_RNDN);
          if (mpfr_cmp_ui(current, 1) < 0) {
            mpfr_set(run[count_run], current, MPFR_RNDN);
          } else {
            mpfr_set_ui(run[count_run], 1, MPFR_RNDN);
          }
          count_run++;
        }
        if (print) cout << "count run = " << count_run << endl;
      }
      if (count_run < run_len && mpfr_cmp_si(run[0], -1) != 0) {
        run_len = count_run;
        for (int l=0; l<count_run; l++) {
          mpfr_set(run_pts[l], run[l], MPFR_RNDN);
        }
        mpc_set(dir_sav, dirs[i], MPFR_RNDN);
        count_run_sav = count_run;
      }
    }
  }
  if (print) {cout << "dir = "; print_mpc(&dir_sav); cout << endl;}
  if (print) cout << "count run = " << count_run << endl;

  // map back [0, 1] to complex half-line
  if (print) cout << "map back [0, 1] to complex half-line" << endl;
  *eta_values = new mpc_t[count_run_sav];
  for (int i=0; i<count_run_sav ; i++) {
    if (print) cout << "i = " << i << endl;
    if (print) {cout << "run_pts = "; print_mpfr(&run_pts[i]); cout << endl;}
    mpfr_sub(tmpfr, fin, ini, MPFR_RNDN);
    mpfr_fma(tmpfr, tmpfr, run_pts[i], ini, MPFR_RNDN);
    mpc_init3((*eta_values)[i], wp2, wp2);
    mpc_mul_fr((*eta_values)[i], dir_sav, tmpfr, MPFR_RNDN);
    // print_mpc(&(*eta_values[i])); cout << endl;
  }
  
  // save number of path values
  *neta_vals = count_run_sav;

  // clear allocated memory
  mpfr_clear(ini);
  mpfr_clear(fin);
  mpfr_clear(min);
  mpfr_clear(max);
  for (int i=0; i<len_max; i++) {
    mpfr_clear(run[i]);
    mpfr_clear(run_pts[i]);
  }
  for (int i=0; i<npoles; i++) {
    mpfr_clear(dists[i]);
    mpc_clear(cen_rot_poles[i]);
  }
  for (int i=0; i<=2*RunCandidate; i++) {
    mpc_clear(dirs[i]);
  }
  mpc_clear(dir_sav);
  mpfr_clear(current);
  mpfr_clear(min_dist);
  mpc_clear(tmpc);
  mpfr_clear(tmpfr);
  mpfr_clear(mpfr_one);

}


void get_path_PS_infty(
  // OUTPUT
  mpc_t **path, int **path_tags, int *neta_vals, int **sing_lab, int *nsings,
  // INPUT
  mpc_t *roots, int nroots, int zero_label
) {
  // specific infty parameters
  *nsings = 1;
  *sing_lab = new int[*nsings];
  (*sing_lab)[0] = 0;

  // CONVERT POLES TO DOUBLE
  complex_d *poles = new complex_d[nroots];
  for (int k=0; k<nroots; k++) {
    if (mpc_lessthan_tol(roots[k])) {
      poles[k] = (complex_d) {0,0};
    } else {
      poles[k] = (complex_d) {
        mpfr_get_d(mpc_realref(roots[k]), MPFR_RNDN),
        mpfr_get_d(mpc_imagref(roots[k]), MPFR_RNDN)
      };
    }  
  }

  // CALL ROUTINE FOR PATH FROM INFTY
  get_path(path, neta_vals, poles, nroots);

  // add zero
  mpc_t *zero = new mpc_t[1];
  mpc_init3(zero[0], wp2, wp2);
  mpc_set_ui(zero[0], 0, MPFR_RNDN);

  mpc_rk1_append(neta_vals, path, 1, zero);
  (*neta_vals)++;

  // for (int k=0; k<*neta_vals; k++) {
  //   cout << "k = " << k << ": "; print_mpc(&(*path)[k]); cout << endl;
  // }
  // exit(0);

  // SET PATH TAGS
  *path_tags = new int[*neta_vals];
  (*path_tags)[0] = 1;
  for (int k=1; k<*neta_vals-1; k++) {
    (*path_tags)[k] = 1;
  }
  (*path_tags)[*neta_vals-1] = 0;

}


void get_path_PS_infty_mp(
  // OUTPUT
  mpc_t **path, int **path_tags, int *neta_vals, int **sing_lab, int *nsings,
  // INPUT
  mpc_t *roots, int nroots, int zero_label
) {
  // specific infty parameters
  *nsings = 1;
  *sing_lab = new int[*nsings];
  (*sing_lab)[0] = 0;

  // prune poles
  mpc_t *poles = new mpc_t[nroots];
  for (int k=0; k<nroots; k++) {
    mpc_init3(poles[k], wp2, wp2);
    if (mpc_lessthan_tol(roots[k])) {
      mpc_set_ui(poles[k], 0, MPFR_RNDN);
    } else {
      mpc_set(poles[k], roots[k], MPFR_RNDN);
    }  
  }

  // CALL ROUTINE FOR PATH FROM INFTY
  get_path_mp(path, neta_vals, poles, nroots);

  // add zero
  mpc_t *zero = new mpc_t[1];
  mpc_init3(zero[0], wp2, wp2);
  mpc_set_ui(zero[0], 0, MPFR_RNDN);

  mpc_rk1_append(neta_vals, path, 1, zero);
  (*neta_vals)++;

  // for (int k=0; k<*neta_vals; k++) {
  //   cout << "k = " << k << ": "; print_mpc(&(*path)[k]); cout << endl;
  // }
  // exit(0);

  // SET PATH TAGS
  *path_tags = new int[*neta_vals];
  (*path_tags)[0] = 1;
  for (int k=1; k<*neta_vals-1; k++) {
    (*path_tags)[k] = 1;
  }
  (*path_tags)[*neta_vals-1] = 0;

}


void path_to_math(
  mpc_t *path, int *path_tags, int neta_values,
  mpc_t *roots, int nroots
) {
  printf("{\n");
  for (int i=0; i<nroots; i++) {
    printf(
      "{%f, %f}",
      mpfr_get_d(mpc_realref(roots[i]), MPFR_RNDN),
      mpfr_get_d(mpc_imagref(roots[i]), MPFR_RNDN)
    );
    if (i < nroots-1) {
      printf(",");
    }
    printf("\n");   
  }
  printf("}\n");
  printf("{\n");
  printf("{\n");
  for (int i=0; i<neta_values; i++) {
    if (path_tags[i] != 2) {
      continue;
    }
    printf(
      "{%f, %f}",
      mpfr_get_d(mpc_realref(path[i]), MPFR_RNDN),
      mpfr_get_d(mpc_imagref(path[i]), MPFR_RNDN)
    );
    if (i < neta_values-1) {
      printf(",");
    }
    printf("\n");
  }
  printf("}\n");
  printf(",\n");
  printf("{\n");
  for (int i=0; i<neta_values; i++) {
    if (path_tags[i] != 1) {
      continue;
    }
    printf(
      "{%f, %f}",
      mpfr_get_d(mpc_realref(path[i]), MPFR_RNDN),
      mpfr_get_d(mpc_imagref(path[i]), MPFR_RNDN)
    );
    if (i < neta_values-1) {
      printf(",");
    }
    printf("\n");
  }
  printf("}\n");
  printf("}\n");
}


void wrt_cmp_path(
  mpc_t **path, int **path_tags, int eps_num, int *neta_values,
  int *nsings, int **sing_lab, int **perm,
  char *file_ext, char *filepath_path, char *filepath_path_tags, char* filepath_sing_lab,
  int opt_write
) {
  char tmp_filepath[200];
  for (int ep=0; ep<eps_num; ep++) {
    if (opt_write) {
      // WRITE TO FILE
      snprintf(tmp_filepath, sizeof(tmp_filepath), "%s%d%s", filepath_path, ep, file_ext);
      cout << endl; cout << "writing to " << tmp_filepath << endl;
      mpc_rk1_to_file(tmp_filepath, path[ep], neta_values[ep]);
      snprintf(tmp_filepath, sizeof(tmp_filepath), "%s%d%s", filepath_path_tags, ep, file_ext);
      cout << "writing to " << tmp_filepath << endl;
      int_rk1_to_file(tmp_filepath, path_tags[ep], neta_values[ep]);
      snprintf(tmp_filepath, sizeof(tmp_filepath), "%s%d%s", filepath_sing_lab, ep, file_ext);
      cout << "writing to " << tmp_filepath << endl;
      int_rk1_to_file(tmp_filepath, sing_lab[ep], nsings[ep]);
    } else {
      // COMPARE PATH POINTS      
      snprintf(tmp_filepath, sizeof(tmp_filepath), "%s%d%s", filepath_path, ep, file_ext);
      cout << endl; cout << "reading path points from file " << tmp_filepath << endl;
      mpc_t *bench_path;
      int bench_neta_values = count_lines(tmp_filepath);
      if (bench_neta_values != neta_values[ep]) {
        cout << "CHECK: FAIL" << endl;
        cout << "benchmark nsings = " << bench_neta_values << endl;
        cout << "nsings = " << neta_values[ep] << endl;
        exit(1);
      }
      bench_path = new mpc_t[bench_neta_values];
      init_rk1_mpc(bench_path, bench_neta_values);
      mpc_rk1_from_file(tmp_filepath, bench_path, bench_neta_values);
      cout << "CHECK path..." << endl;
      mpc_rk1_compare_double(bench_path, path[ep], neta_values[ep]);

      // COMPARE PATH TAGS
      snprintf(tmp_filepath, sizeof(tmp_filepath), "%s%d%s", filepath_path_tags, ep, file_ext);
      cout << "reading path tags from file " << tmp_filepath << endl;
      int *bench_path_tags = new int[neta_values[ep]];
      int_rk1_from_file(tmp_filepath, bench_path_tags, neta_values[ep]);
      cout << "CHECK path tags..." << endl;
      int_rk1_compare(bench_path_tags, path_tags[ep], neta_values[ep]);

      // COMPARE SINGULAR LABELS
      // permutation of roots is needed in order to compare singular labels
      snprintf(tmp_filepath, sizeof(tmp_filepath), "%s%d%s", filepath_sing_lab, ep, file_ext);
      cout << "reading sing labs from file " << tmp_filepath << endl;
      int bench_nsings, *bench_sing_lab;
      bench_nsings = count_lines(tmp_filepath);
      if (bench_nsings != nsings[ep]) {
        cout << "CHECK: FAIL" << endl;
        cout << "benchmark nsings = " << bench_nsings << endl;
        cout << "nsings = " << nsings[ep] << endl;
        exit(1);
      }
      bench_sing_lab = new int[bench_nsings];
      int_rk1_from_file(tmp_filepath, bench_sing_lab, bench_nsings);
      cout << "CHECK sing labs..." << endl;
      int_rk1_compare_perm(perm[ep], bench_sing_lab, sing_lab[ep], nsings[ep]);
    }
  }
}


//////
// NESTED LISTS
//////
void nlist_build(nlist *nl) {
  nl->nitems = 0;
  nl->item = NULL;
  nl->str = NULL;
}


void nlist_self_parse(nlist *nl) {
  if (nl->str[0] != '[') {
    // bottom of the tree
    nl->nitems = 0;
    nl->item = NULL;
    return;
  }

  if (nl->str[1] == ']') {
    // empty list
    nl->nitems = 0;
    nl->item = NULL;
    return;
  }

  // get number of elements
  int del_count = 1, c = 1;
  nl->nitems = 0;
  while (del_count > 0) {
    if (nl->str[c] == '[') {
      del_count++;
    } else if (nl->str[c] == ']') {
      del_count--;
      if (del_count == 0 && nl->str[c-1] != ',') {
        nl->nitems++;
      }
    } else if (nl->str[c] == ',' && del_count == 1) {
      nl->nitems++;
    }
    c++;
  }
  // cout << "nitems = " << nl->nitems << endl;

  // prepare child items
  nl->item = new nlist[nl->nitems];
  for (int i=0; i<nl->nitems; i++) {
    nlist_build(&nl->item[i]);
    nl->item[i].str = (char*) malloc(MAX_VALUE_LENGTH*sizeof(char));
  }

  // get string for child items
  del_count = 1, c = 1;
  int item_count = 0, cp = 0;
  while (item_count < nl->nitems) {
    // cout << "c, nl->str[c] = " << c << ", " << nl->str[c] << endl;
    // getchar();
    if (nl->str[c] == '[') {
      del_count++;
    } else if (nl->str[c] == ']') {
      del_count--;
      if (del_count == 0) {
        nl->item[item_count].str[cp] = '\0';
        // cout << "item " << item_count << ": " << nl->item[item_count].str << endl;
        item_count++;
        continue;
      }
    } else if (nl->str[c] == ',' && del_count == 1) {
      (nl->item[item_count]).str[cp] = '\0';
      // cout << "item " << item_count << ": " << nl->item[item_count].str << endl;
      item_count++;
      cp = 0;
      c++;
      continue;
    }

    // cout << "store in item " << item_count << endl;
    // getchar();
    nl->item[item_count].str[cp++] = nl->str[c++];
    // cout << "accepted" << endl;
    // getchar();
  }

  // recursively parse child items
  for (int i=0; i<nl->nitems; i++) {
    nlist_self_parse(&nl->item[i]);
  }
  
}


int nlist_read_file(
  // OUTPUT
  nlist *nl,
  // INPUT
  char *key, FILE *fptr
) {
  const size_t line_size = MAX_VALUE_LENGTH * sizeof(char);
  char *line = (char*) malloc(line_size);
  if (line == NULL) {
    perror("error while allocating memory for line");
    exit(1);
  }
  if (!goto_key(&line, key, fptr)) {
    free(line);
    return 0;
  }

  // go to ':'
  int i = 0;
  while(line[i] != ':') {
    i++;
  }
  i++;

  // goto '['
  while(line[i] != '[') {
    i++;
  }
  nl->str = (char*) malloc(MAX_VALUE_LENGTH*sizeof(char));
  int c = 0;
  nl->str[c++] = line[i++];

  // get no-space string
  int del_count = 1;
  while (del_count > 0) {
    if (line[i] == '\n') {
      fgets(line, MAX_VALUE_LENGTH, fptr);
      i = 0;
      continue;
    }

    if (isspace(line[i])) {
      i++;
      continue;
    }
        
    if (line[i] == '[') {
      del_count++;
    } else if (line[i] == ']') {
      del_count--;
    }

    nl->str[c++] = line[i++];
  }
  nl->str[c] = '\0';
  // cout << "string:" << endl;
  // cout << nl->str << endl;

  // read child items
  nlist_self_parse(nl);

  free(line);
  return 1;
}


void nlist_free(nlist *nl) {
  if (nl->item != NULL) {
    for (int i=0; i<nl->nitems; i++) {
      nlist_free(&nl->item[i]);
    }
    delete[] nl->item;
  }
  if (nl->str != NULL) {
    free(nl->str);  // Free the allocated memory for str
  }
}


void nlist_rk1_free(nlist *nl, int dim) {
  for (int i=0; i<dim; i++) {
    nlist_free(&nl[i]);
  }
}

