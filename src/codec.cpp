#include <iostream>
#include <cstring>
#include "mpc.h"

using namespace std;

#include "global_vars.h"
// #include "malloc_defs.h"
#include "utils.h"
#include "tensor_utils.h"
#include "poly_frac.h"
#include "conversions.h"
#include "system_analyzer.h"
extern "C" {
	#include "mp_roots_poly.h"
	#include "cpoly.h"
}


//////
// FORWARD DECLARATION
//////
void decode_tree_pf(
	// OUTPUT
	struct poly_frac *out,
	// IN-OUT
	int *nroots, mpc_t **roots, mpfr_t **tols,
	// INPUT
	struct poly_frac *in,
	char ***kin, int *skip_inv, int *is_mass,
	char **s_tree, char *sep,
	int wp2, mpfr_t mpfr_tol, int incr_prec
);

//////
// ROUTINES FOR POLYNOMIALS BY COEFFICIENTS
//////
struct cpoly{
	int deg;
	mpc_t *coeffs;
};


void cpoly_build(
	struct cpoly *pol
) {
	(*pol).coeffs = NULL;
}


void cpoly_print(
	struct cpoly *pol
) {
	if ((*pol).deg == -1) {
		cout << "0" << endl;
		return;
	}

	cout << "degree: " << (*pol).deg << endl;
	if ((*pol).coeffs) {
		for (int k=0; k<=(*pol).deg; k++) {
			cout << "pow " << k << ": ";
			mpc_out_str(stdout, 10, 0, (*pol).coeffs[k], MPFR_RNDN);
			cout << endl;
		}
	}
}


void cpoly_set_zero(
	struct cpoly *pol
) {
	(*pol).deg = -1;
	if ((*pol).coeffs) {
		delete[] (*pol).coeffs;
	}
}


void cpoly_set_ui(
	struct cpoly *pol,
	int ui_num
) {
	if (ui_num == 0) {
		cpoly_set_zero(pol);
		return;
	}

	(*pol).deg = 0;
	if ((*pol).coeffs) {
		delete[] (*pol).coeffs;
	}
	(*pol).coeffs = new mpc_t[1];
	mpc_init3((*pol).coeffs[0], wp2, wp2);
	mpc_set_ui((*pol).coeffs[0], ui_num, MPFR_RNDN);
}


void cpoly_set_str(
	struct cpoly *pol,
	char *str
) {
	(*pol).deg = 0;
	if ((*pol).coeffs) {
		delete[] (*pol).coeffs;
	}
	(*pol).coeffs = new mpc_t[1];
	mpc_init3((*pol).coeffs[0], wp2, wp2);
	mpc_set_str((*pol).coeffs[0], str, 10, MPFR_RNDN);
}


void cpoly_set_pol(
	struct cpoly *polout,
	struct cpoly *polin
) {
	// check if input is zero
	if ((*polin).deg == -1) {
		cpoly_set_zero(polout);
		return;
	}

	(*polout).deg = (*polin).deg;
	if ((*polout).coeffs) {
		delete[] (*polout).coeffs;
	}
	(*polout).coeffs = new mpc_t[(*polout).deg+1];
	for (int k=0; k<=(*polin).deg; k++) {
		mpc_init3((*polout).coeffs[k], wp2, wp2);
		mpc_set((*polout).coeffs[k], (*polin).coeffs[k], MPFR_RNDN);
	}
}


void cpoly_add_pol(
	struct cpoly *pol,
	struct cpoly *pol1,
	struct cpoly *pol2
) {
	// check if one of the input is zero
	if ((*pol1).deg == -1) {
		cpoly_set_pol(pol, pol2);
		return;
	} else if ((*pol2).deg == -1) {
		cpoly_set_pol(pol, pol1);
		return;
	}

	// deal with overwriting
	struct cpoly *in1, *in2;
	if (pol == pol1) {
		in1 = new struct cpoly[1];
		cpoly_build(in1);
		cpoly_set_pol(in1, pol1);
		in2 = pol2;
	} else if (pol == pol2) {
		in1 = pol1;
		in2 = new struct cpoly[1];
		cpoly_build(in2);
		cpoly_set_pol(in2, pol2);
	} else {
		in1 = pol1;
		in2 = pol2;
	}

	// determine max and min degree
	int deg, deg_min, first_greater;
	if ((*pol1).deg >= (*pol2).deg) {
		deg = (*pol1).deg;
		deg_min = (*pol2).deg;
		first_greater = 1;
	} else {
		deg = (*pol2).deg;
		deg_min = (*pol1).deg;
		first_greater = 0;
	}
	
	// alloc memory
	if ((*pol).coeffs && (*pol).deg != deg) {
		delete[] (*pol).coeffs;
		(*pol).deg = deg;
		(*pol).coeffs = new mpc_t[(*pol).deg+1];
		for (int k=0; k<=(*pol).deg; k++) {
			mpc_init3((*pol).coeffs[k], wp2, wp2);
		}
	}

	// add
	for (int k=0; k<=deg_min; k++) {
		mpc_add((*pol).coeffs[k], (*in1).coeffs[k], (*in2).coeffs[k], MPFR_RNDN);
	}
	for (int k=deg_min+1; k<=(*pol).deg; k++) {
		if (first_greater) {
			mpc_set((*pol).coeffs[k], (*in1).coeffs[k], MPFR_RNDN);
		} else {
			mpc_set((*pol).coeffs[k], (*in2).coeffs[k], MPFR_RNDN);
		}
	}
	return;
}


void mul_poly_no_trunc(mpc_t *out_pol, mpc_t *pol1, mpc_t *pol2, int deg1, int deg2) {
	/*
	Perform multiplication of two polynomials:

		pol1 = \sum_{i=0}^{deg1} a_i * x^i
		pol2 = \sum_{j=0}^{deg2} b_j * x^j

		pol1*pol2 = \sum_{k=0}^{deg1+deg2} \sum_{i=\max{0, k-deg2}}^{k} a_i * b_{k-i}

		Note that deg1 has to be greater than deg2.
	*/
	int i, k, i_min=0;
	mpc_t tmp;
	mpc_init3(tmp, wp2, wp2);
	for (k=0; k<=deg1+deg2; k++) {
		// mpc_init3(out_pol[k], wp2, wp2);
		mpc_set_d(out_pol[k], 0, MPFR_RNDN);
		if (k > deg2) {
			i_min = k - deg2;
		}
		if (k < deg1) {
			for (i=i_min; i<=k; i++) {
				mpc_fma(out_pol[k], pol1[i], pol2[k-i], out_pol[k], MPFR_RNDN);
			}			
		} else {
			for (i=i_min; i<=deg1; i++) {
				mpc_fma(out_pol[k], pol1[i], pol2[k-i], out_pol[k], MPFR_RNDN);
			}
		}
	}
}


void cpoly_mul_pol(
	struct cpoly *pol,
	struct cpoly *pol1,
	struct cpoly *pol2
) {
	// if input is zero, set output to zero
	if ((*pol1).deg == -1) {
		if (pol != pol1) {
			cpoly_set_zero(pol);
		}
		return;
	}
	if ((*pol2).deg == -1) {
		if (pol != pol2) {
			cpoly_set_zero(pol);
		}
		return;
	}

	// deal with overwriting
	struct cpoly *in1, *in2;
	if (pol == pol1) {
		in1 = new struct cpoly[1];
		cpoly_build(in1);
		cpoly_set_pol(in1, pol1);
		in2 = pol2;
	} else if (pol == pol2) {
		in1 = pol1;
		in2 = new struct cpoly[1];
		cpoly_build(in2);
		cpoly_set_pol(in2, pol2);
	} else {
		in1 = pol1;
		in2 = pol2;
	}

	// alloc memory
	int deg = (*pol1).deg+(*pol2).deg;
	if ((*pol).coeffs && (*pol).deg != deg) {
		delete[] (*pol).coeffs;
		(*pol).deg = deg;
		(*pol).coeffs = new mpc_t[(*pol).deg+1];
		for (int k=0; k<=(*pol).deg; k++) {
			mpc_init3((*pol).coeffs[k], wp2, wp2);
		}
	}

	// multiply
	if ((*in1).deg < (*in2).deg) {
		struct cpoly *tmp;	
		tmp = in1;
		in1 = in2;
		in2 = tmp;
	}
	mul_poly_no_trunc((*pol).coeffs, (*in1).coeffs, (*in2).coeffs, (*in1).deg, (*in2).deg);
	return;
}


void cpoly_pow_ui(
	struct cpoly *pol,
	struct cpoly *pol1,
	int pow
) {
	// if input is zero, set output to zero
	if ((*pol1).deg == -1) {
		if (pol != pol1) {
			cpoly_set_zero(pol);
		}
		return;
	}

	// if pow is zero, set output to one
	if (pow == 0) {
		if (pol != pol1) {
			cpoly_set_ui(pol, 1);
		}
		return;
	}

	// if pow is one, set output to input
	if (pow == 1) {
		if (pol != pol1) {
			cpoly_set_pol(pol, pol1);
		}
		return;
	}

	// deal with overwriting
	struct cpoly *in1;
	if (pol == pol1) {
		in1 = new struct cpoly[1];
		cpoly_build(in1);
		cpoly_set_pol(in1, pol1);
	} else {
		in1 = pol1;
	}

	// multiply
	cpoly_set_pol(pol, in1);
	for (int p=1; p<pow; p++) {
		cpoly_mul_pol(pol, pol, in1);
	}
}


//////
// UTILS
//////
// int sym_idx(ex sym, symbol *sym_vec, int num_syms) {
// 	int s;	
// 	for (s=0; s<num_syms; s++) {
// 		if (sym.is_equal(sym_vec[s])) {
// 			return s;
// 		}
// 	}
// 	return -1;
// }

// void get_kin_inv_coeffs(
//   // OUTPUT
//   int 
//   // INPUT
//   ex mpoly,
//   symbol *g_s,
//   int num_inv
// ) {

// }


//////
// ROUTINES FOR ENCODER-DECODER
//////
// static void my_print(const ex & e) {
// void my_print(const ex & e) {
// 	cout << ex_to<basic>(e).class_name();
// 	cout << "(";
// 	size_t n = e.nops();
// 	if (n) {
// 		for (size_t i=0; i<n; i++) {
// 			my_print(e.op(i));
// 			if (i != n-1) {
// 				cout << ",";
// 			}
// 		}
// 	} else {
// 		cout << e;
// 	}
// 	cout << ")";
// }


// static void denom_factors(const ex & gncpf) {
// 	size_t n = gncpf.nops();
// 	char* kind = (char*) ex_to<basic>(gncpf).class_name();
// 	// cout << "enter denom_factors" << endl;
// 	// cout << "input: " << gncpf << endl;
// 	// cout << "kind: " << kind << endl;
// 	// cout << "n = " << n << endl;
// 	if (strcmp(kind, "power") == 0) {
// 		cout << "factor:" << endl;
// 		cout << gncpf << endl;
// 		cout << "power:" << endl;
// 		cout << gncpf.op(1) << endl;
// 		cout << "poly:" << endl;
// 		cout << gncpf.op(0) << endl;
// 		return;
// 	}
// 	if (n) {
// 		for (size_t i=0; i<n; i++) {
// 			denom_factors(gncpf.op(i));
// 		}
// 	}
// }


// // static void encode_tree(
// void encode_tree(
// 	// OUTPUT
// 	stringstream *stringa,
// 	// INPUT
// 	const ex & e,
// 	symbol *sym_vec,
// 	int num_syms
// ) {
// 	char sep[] = ", ";
// 	size_t n = e.nops();
// 	char* kind = (char*) ex_to<basic>(e).class_name();
// 	int idx;
// 	// if (strcmp(kind, "power") == 0) {
// 	// 	(*stringa) << kind[0] << sep;
// 	// 	// (*stringa) << e.op(0) << sep;
// 	// 	// recognize symbol
// 	// 	idx = sym_idx(e.op(0), sym_vec, num_syms);
// 	// 	(*stringa) << idx << sep;
// 	// 	(*stringa) << e.op(1) << sep;
// 	// 	return;
// 	// }
// 	if (strcmp(kind, "symbol") == 0) {
// 		// (*stringa) << kind[0] << sep;
// 		(*stringa) << 's' << sep;
// 		// recognize symbol
// 		idx = sym_idx(e, sym_vec, num_syms);
// 		(*stringa) << idx << sep;
// 		// (*stringa) << e << sep;
// 		return;
// 	}
// 	if (strcmp(kind, "numeric") == 0) {
// 		// (*stringa) << kind[0] << sep;
// 		(*stringa) << 'n' << sep;
// 		(*stringa) << e << sep;
// 		return;
// 	}
// 	(*stringa) << "o" << sep;
// 	// (*stringa) << kind[0] << sep;
// 	if (n) {
// 		if (strcmp(kind, "power") == 0) {
// 			if (e.op(1) > 0) {
// 				(*stringa) << "p" << sep;
// 			} else {
// 				(*stringa) << "i" << sep;
// 			}

// 			(*stringa) << e.op(1) << sep;
// 			encode_tree(stringa, e.op(0), sym_vec, num_syms);
// 		} else {
// 			// (*stringa) << kind[0] << sep;
// 			switch (kind[0]) {
// 				case 'a':
// 				(*stringa) << '+' << sep;
// 				break;
// 				case 'm':
// 				(*stringa) << '*' << sep;
// 				break;
// 				default:
// 				(*stringa) << kind[0] << sep;
// 				break;
// 			}
// 			(*stringa) << n << sep;
// 			for (size_t i=0; i<n; i++) {
// 				encode_tree(stringa, e.op(i), sym_vec, num_syms);
// 				// if (i != n-1) {
// 				// 	cout << endl;
// 				// }
// 			}
// 		}
// 	} else {
// 		(*stringa) << kind[0] << sep;
// 		(*stringa) << e;
// 	}
// }


void decode_tree(
	// OUTPUT
	mpc_t *out,
	// INPUT
	mpc_t *in,
	char **s_tree,
	char *sep
) {
	// FILE *fp = fopen(tree_file, "r");
	// for (char c = getc(fp); c != EOF; c = getc(fp)) {
	// }
	cout << endl;
	cout << "ENTER decode_tree" << endl;
	cout << "input tree:" << endl;
	cout << (*s_tree) << endl;
	char *s_tree_copy = (char*)malloc(sizeof(char)*strlen((*s_tree)));
	char *s_tree_copy2 = (char*)malloc(sizeof(char)*strlen((*s_tree)));
	strcpy(s_tree_copy, (*s_tree));
	strcpy(s_tree_copy2, (*s_tree));
	int stride = 0;
	int n, pow, idx;
	char *operation;
	mpc_t *inner_mpc;
	char *token = strtok((*s_tree), (char*) sep);
	stride += strlen(token) + strlen(sep);
	while (token) {
		printf("token: %s\n", token);
		if (strcmp(token, "o") == 0) {
			cout << "START operation" << endl;
			token = strtok(NULL, (char*) sep);
			stride += strlen(token) + strlen(sep);

			operation = (char*) malloc(sizeof(char)*strlen(token));
			strcpy(operation, token);
			printf("token: %s\n", token);
			
			// for power and division, n=1, otherwise it is written in the string
			if (strcmp(operation, "p") == 0 && strcmp(operation, "i") == 0) {
				n = 1;
			} else {
				token = strtok(NULL, (char*) sep);
				stride += strlen(token) + strlen(sep);
				printf("token: %s\n", token);
				n = atoi(token);
			}
			
			// sscanf(token, "%04d", &n);
			cout << "n = " << n << endl;
			inner_mpc = new mpc_t[n];
			for (int i=0; i<n; i++) {
				cout << "i = " << i << endl;
				mpc_init3(inner_mpc[i], wp2, wp2);
				strcpy((*s_tree), s_tree_copy + stride);
				decode_tree(&inner_mpc[i], in, s_tree, sep);
				// cout << "recovered tree:" << endl;
				// cout << *s_tree << endl;
				// printf("%p\n", *s_tree);
				free(s_tree_copy);
				s_tree_copy = (char*)malloc(sizeof(char)*strlen((*s_tree)));
				strcpy(s_tree_copy, (*s_tree));
				stride = 0;
			}
			if (strcmp(operation, "a") == 0) {
				cout << "adding " << n << " output" << endl;
				mpc_set_ui(*out, 0, MPFR_RNDN);
				mpc_out_str(stdout, 10, 0, *out, MPFR_RNDN); cout << endl;
				for (int i=0; i<n; i++) {
					mpc_add(*out, *out, inner_mpc[i], MPFR_RNDN);
					mpc_out_str(stdout, 10, 0, *out, MPFR_RNDN); cout << endl;
				}
			} else if (strcmp(operation, "m") == 0) {
				cout << "multiplying " << n << " output" << endl;
				mpc_set_ui(*out, 1, MPFR_RNDN);
				mpc_out_str(stdout, 10, 0, *out, MPFR_RNDN); cout << endl;
				for (int i=0; i<n; i++) {
					mpc_mul(*out, *out, inner_mpc[i], MPFR_RNDN);
					mpc_out_str(stdout, 10, 0, *out, MPFR_RNDN); cout << endl;
				}
			}
			cout << "END operation" << endl;
			return;
		} else {
			cout << "START ELSE" << endl;
			if (strcmp(token, "p") == 0) {
				token = strtok(NULL, (char*) sep);
				stride += strlen(token) + strlen(sep);
				idx = atoi(token);
				token = strtok(NULL, (char*) sep);
				stride += strlen(token) + strlen(sep);
				pow = atoi(token);
				(*s_tree) = (char*)malloc(sizeof(char)*strlen(s_tree_copy + stride));
				strcpy((*s_tree), s_tree_copy + stride);
				cout << "SET power" << endl;
				mpc_out_str(stdout, 10, 0, in[idx], MPFR_RNDN); cout << endl;
				mpc_pow_si(*out, in[idx], pow, MPFR_RNDN);
				mpc_out_str(stdout, 10, 0, *out, MPFR_RNDN); cout << endl;
			} else if (strcmp(token, "s") == 0) {
				token = strtok(NULL, (char*) sep);
				stride += strlen(token) + strlen(sep);
				idx = atoi(token);
				(*s_tree) = (char*)malloc(sizeof(char)*strlen(s_tree_copy + stride));
				strcpy((*s_tree), s_tree_copy + stride);
				cout << "SET symbol" << endl;
				mpc_out_str(stdout, 10, 0, in[idx], MPFR_RNDN); cout << endl;
				mpc_set(*out, in[idx], MPFR_RNDN);
				mpc_out_str(stdout, 10, 0, *out, MPFR_RNDN); cout << endl;
			} else if (strcmp(token, "n") == 0) {
				token = strtok(NULL, (char*) sep);
				stride += strlen(token) + strlen(sep);
				(*s_tree) = (char*)malloc(sizeof(char)*strlen(s_tree_copy + stride));
				strcpy((*s_tree), s_tree_copy + stride);
				cout << "SET numeric" << endl;
				cout << token << endl;
				mpc_set_str(*out, token, 10, MPFR_RNDN);
				mpc_out_str(stdout, 10, 0, *out, MPFR_RNDN); cout << endl;
			}
			return;
		}
		token = strtok(NULL, (char*) sep);
		stride += strlen(token) + strlen(sep);
	}
}


void decode_tree_mpq(
	// OUTPUT
	mpq_t *out,
	// INPUT
	mpq_t *in,
	char **s_tree,
	char *sep
) {
	// FILE *fp = fopen(tree_file, "r");
	// for (char c = getc(fp); c != EOF; c = getc(fp)) {
	// }
	// cout << endl;
	// cout << "ENTER decode_tree" << endl;
	// cout << "input tree:" << endl;
	// cout << (*s_tree) << endl;
	char *s_tree_copy = (char*)malloc(sizeof(char)*(strlen((*s_tree))+1));
	// char *s_tree_copy2 = (char*)malloc(sizeof(char)*(strlen((*s_tree))+1));
	strcpy(s_tree_copy, (*s_tree));
	// strcpy(s_tree_copy2, (*s_tree));
	int stride = 0;
	int n, pow, idx;
	char *operation;
	mpq_t *inner_mpq;
	char *token = strtok((*s_tree), (char*) sep);
	stride += strlen(token) + strlen(sep);
	while (token) {
		if (dbg) printf("token: %s\n", token);
		if (strcmp(token, "o") == 0) {
			if (dbg) cout << "START operation" << endl;
			token = strtok(NULL, (char*) sep);
			if (dbg) printf("token: %s\n", token);
			stride += strlen(token) + strlen(sep);
			operation = (char*) malloc(sizeof(char)*(strlen(token)+1));
			strcpy(operation, token);
			
			// for power and division, n=1, otherwise it is written in the string
			if (strcmp(operation, "p") == 0 || strcmp(operation, "i") == 0) {
				n = 1;
				token = strtok(NULL, (char*) sep);
				// printf("token: %s\n", token);
				stride += strlen(token) + strlen(sep);
				pow = atoi(token);
			} else if (strcmp(operation, "/") == 0) {
				n = 1;
				pow = -1;
			} else if (strcmp(operation, "-") == 0) {
				n = 1;
			} else {
				token = strtok(NULL, (char*) sep);
				// printf("token: %s\n", token);
				stride += strlen(token) + strlen(sep);
				n = atoi(token);
			}
			if (dbg) cout << "n = " << n << endl;

			inner_mpq = new mpq_t[n];
			for (int i=0; i<n; i++) {
				if (dbg) cout << "i = " << i << endl;
				mpq_init(inner_mpq[i]);
				free(*s_tree);
				*s_tree = (char*)malloc(sizeof(char)*(strlen(s_tree_copy + stride)+1));
				strcpy((*s_tree), s_tree_copy + stride);
				decode_tree_mpq(&inner_mpq[i], in, s_tree, sep);
				free(s_tree_copy);
				s_tree_copy = (char*)malloc(sizeof(char)*(strlen((*s_tree))+1));
				strcpy(s_tree_copy, (*s_tree));
				stride = 0;
			}
			if (strcmp(operation, "+") == 0) {
				if (dbg) cout << "adding " << n << " output" << endl;
				mpq_set_ui(*out, 0, 1);
				for (int i=0; i<n; i++) {
					mpq_add(*out, *out, inner_mpq[i]);
				}
			} else if (strcmp(operation, "-") == 0) {
				mpq_set(*out, inner_mpq[0]);
				mpq_neg(*out, *out);
			} else if (strcmp(operation, "*") == 0) {
				if (dbg) cout << "multiplying " << n << " output" << endl;
				mpq_set_ui(*out, 1, 1);
				for (int i=0; i<n; i++) {
					mpq_mul(*out, *out, inner_mpq[i]);
				}
				// cpoly_print(out); cout << endl;
			} else if (strcmp(operation, "/") == 0) {
				mpq_inv(*out, inner_mpq[0]);
			} else if (strcmp(operation, "p") == 0) {
				if (dbg) cout << "raising to power" << endl;
				// int pow = mpfr_get_si(mpc_realref(inner_pf[1].coeffs[0]), MPFR_RNDN);
				if (dbg) cout << "pow = " << pow << endl;
				mpz_pow_ui(mpq_numref(*out), mpq_numref(inner_mpq[0]), pow);
				mpz_pow_ui(mpq_denref(*out), mpq_denref(inner_mpq[0]), pow);
				// cout << "output:" << endl;
				// poly_frac_print(out);
				// cout << *s_tree << endl;
			} else if (strcmp(operation, "i") == 0) {
				mpz_pow_ui(mpq_numref(*out), mpq_denref(inner_mpq[0]), pow);
				mpz_pow_ui(mpq_denref(*out), mpq_numref(inner_mpq[0]), pow);
			}
			if (dbg) cout << "END operation" << endl;
			mpq_rk1_clear(inner_mpq, n); delete[] inner_mpq;
			free(s_tree_copy);
			return;
		} else {
			if (dbg) cout << "START ELSE" << endl;
			if (strcmp(token, "s") == 0) {
				token = strtok(NULL, (char*) sep);
				stride += strlen(token) + strlen(sep);
				idx = atoi(token);
				free(*s_tree);
				(*s_tree) = (char*)malloc(sizeof(char)*(strlen(s_tree_copy + stride)+1));
				strcpy((*s_tree), s_tree_copy + stride);
				if (dbg) cout << "SET symbol" << endl;
				mpq_set(*out, in[idx]);
			} else if (strcmp(token, "n") == 0) {
				token = strtok(NULL, (char*) sep);
				if (dbg) cout << "SET numeric" << endl;
				if (dbg) cout << token << endl;
				mpq_set_str(*out, token, 10);
				stride += strlen(token) + strlen(sep);
				free(*s_tree);
				(*s_tree) = (char*)malloc(sizeof(char)*(strlen(s_tree_copy + stride)+1));
				strcpy((*s_tree), s_tree_copy + stride);
			}
			free(s_tree_copy);
			return;
		}
		token = strtok(NULL, (char*) sep);
		stride += strlen(token) + strlen(sep);
	}
	free(s_tree_copy);
}


void decode_tree_poly(
	// OUTPUT
	struct cpoly *out,
	// INPUT
	struct cpoly *in,
	char **s_tree,
	char *sep
) {
	// FILE *fp = fopen(tree_file, "r");
	// for (char c = getc(fp); c != EOF; c = getc(fp)) {
	// }
	cout << endl;
	cout << "ENTER decode_tree" << endl;
	cout << "input tree:" << endl;
	cout << (*s_tree) << endl;
	char *s_tree_copy = (char*)malloc(sizeof(char)*(strlen((*s_tree))+1));
	// char *s_tree_copy2 = (char*)malloc(sizeof(char)*(strlen((*s_tree))+1));
	strcpy(s_tree_copy, (*s_tree));
	// strcpy(s_tree_copy2, (*s_tree));
	int stride = 0;
	int n, pow, idx;
	char *operation;
	struct cpoly *inner_pol;
	char *token = strtok((*s_tree), (char*) sep);
	stride += strlen(token) + strlen(sep);
	while (token) {
		printf("token: %s\n", token);
		if (strcmp(token, "o") == 0) {
			cout << "START operation" << endl;
			token = strtok(NULL, (char*) sep);
			stride += strlen(token) + strlen(sep);
			printf("token: %s\n", token);
			n = atoi(token);
			// sscanf(token, "%04d", &n);
			cout << "n = " << n << endl;
			inner_pol = new struct cpoly[n];
			token = strtok(NULL, (char*) sep);
			stride += strlen(token) + strlen(sep);
			operation = (char*) malloc(sizeof(char)*(strlen(token)+1));
			strcpy(operation, token);
			printf("token: %s\n", token);
			for (int i=0; i<n; i++) {
				cout << "i = " << i << endl;
				cpoly_build(&inner_pol[i]);
				*s_tree = (char*)malloc(sizeof(char)*(strlen(s_tree_copy + stride)+1));
				strcpy((*s_tree), s_tree_copy + stride);
				decode_tree_poly(&inner_pol[i], in, s_tree, sep);
				// cout << "recovered tree:" << endl;
				// cout << *s_tree << endl;
				// printf("%p\n", *s_tree);
				free(s_tree_copy);
				s_tree_copy = (char*)malloc(sizeof(char)*(strlen(*s_tree)+1));
				strcpy(s_tree_copy, (*s_tree));
				stride = 0;
			}
			if (strcmp(operation, "a") == 0) {
				cout << "adding " << n << " output" << endl;
				cpoly_set_zero(out);
				for (int i=0; i<n; i++) {
					cpoly_add_pol(out, out, &inner_pol[i]);
				}
				// cpoly_print(out); cout << endl;
			} else if (strcmp(operation, "m") == 0) {
				cout << "multiplying " << n << " output" << endl;
				cpoly_set_ui(out, 1);
				for (int i=0; i<n; i++) {
					cpoly_mul_pol(out, out, &inner_pol[i]);
				}
				// cpoly_print(out); cout << endl;
			} else if (strcmp(operation, "p") == 0) {
				cout << "raising to power" << endl;
				cpoly_pow_ui(out, &inner_pol[0], mpfr_get_ui(mpc_realref(inner_pol[1].coeffs[0]), MPFR_RNDN));
			}
			cout << "END operation" << endl;
			return;
		} else {
			cout << "START ELSE" << endl;
			// if (strcmp(token, "p") == 0) {
			// 	token = strtok(NULL, (char*) sep);
			// 	stride += strlen(token) + strlen(sep);
			// 	idx = atoi(token);
			// 	token = strtok(NULL, (char*) sep);
			// 	stride += strlen(token) + strlen(sep);
			// 	pow = atoi(token);
			// 	(*s_tree) = (char*)malloc(sizeof(char)*(strlen(s_tree_copy + stride)+1));
			// 	strcpy((*s_tree), s_tree_copy + stride);
			// 	cout << "SET power" << endl;
			// 	cout << "pow = " << pow << endl;
			// 	// cpoly_print(&in[idx]);
			// 	cpoly_pow_ui(out, &in[idx], pow);
			// 	// cpoly_print(out);
			// 	cout << endl;
			if (strcmp(token, "s") == 0) {
				token = strtok(NULL, (char*) sep);
				stride += strlen(token) + strlen(sep);
				idx = atoi(token);
				(*s_tree) = (char*)malloc(sizeof(char)*(strlen(s_tree_copy + stride)+1));
				strcpy((*s_tree), s_tree_copy + stride);
				cout << "SET symbol" << endl;
				cpoly_set_pol(out, &in[idx]);
				// cpoly_print(out);
				cout << endl;
			} else if (strcmp(token, "n") == 0) {
				token = strtok(NULL, (char*) sep);
				stride += strlen(token) + strlen(sep);
				(*s_tree) = (char*)malloc(sizeof(char)*(strlen(s_tree_copy + stride)+1));
				strcpy((*s_tree), s_tree_copy + stride);
				cout << "SET numeric" << endl;
				cout << token << endl;
				cpoly_set_str(out, token);
				// cpoly_print(out);
				cout << endl;
			}
			return;
		}
		token = strtok(NULL, (char*) sep);
		stride += strlen(token) + strlen(sep);
	}
}


// #pair: poly_mpq_find_roots
void poly_frac_numer_roots_from_tree(
	// OUTPUT
	struct poly_frac *pf, mpc_t **roots, mpfr_t **tols,
	int *vdeg, int *ldeg,
	// IN-OUT
	int *in_nroots, mpc_t **in_roots, mpfr_t **in_tols,
	// INPUT
	char **s_tree, char *sep,
	struct poly_frac *in,
	char ***kin, int *skip_inv, int *is_mass,
	int wp2, mpfr_t mpfr_tol, int incr_prec, int curr_mul
) {
	// cout << endl; cout << "ENTER ROUTINE FOR NUMERATOR ROOTS" << endl;
	// cout << "INPUT:" << endl;
	// cout << "string: " << (*s_tree) << endl;
	// cout << "wp2 = " << wp2 << endl;
	// cout << "curr_mul = " << curr_mul << endl;

	mpc_t tmpc;
	mpc_init3(tmpc, wp2, wp2);
	mpfr_t tmp_mpfr;
	mpfr_init(tmp_mpfr);

	//////
	// GET POLY_FRAC TO BE PROCESSED
	//////
	// cout << "GET POLY_FRAC TO BE PROCESSED" << endl;
	// struct poly_frac pfin;
	// poly_frac_build(&pfin);
	// *s_tree = (char*)malloc(sizeof(char)*(strlen(s_tree_copy + stride)+1));
	// strcpy((*s_tree), s_tree_copy + stride);
	// copy tree
	char *s_tree_copy = (char*)malloc(sizeof(char)*(strlen((*s_tree))+1));
	strcpy(s_tree_copy, (*s_tree));
	// cout << "s_tree_copy = " << s_tree_copy << endl;
	decode_tree_pf(
		pf,
		in_nroots, in_roots, in_tols,
		in,
		kin, skip_inv, is_mass,
		&s_tree_copy, sep,
		wp2, mpfr_tol, incr_prec
	);
	// cout << "pf to be processed:" << endl;
	// poly_frac_print(pf);
	// cout << "original tree: " << (*s_tree) << endl;
	// cout << "residual tree: " << (s_tree_copy) << endl;

	//////
	// FIRST CALL TO ROOT FINDER
	//////
	// #2BD: per polinomi di rango 2 manipolare il tree per salvarsi
	// il tree delle radici in funzione dei punti di spazio fase,
	// così da salvare solo quello in libreria e non dover cercare le
	// radici numeriche al variare dei punti di spazio fase
	*vdeg = (*pf).num_vdeg;
	*ldeg = (*pf).num_deg - (*pf).num_vdeg;
	*roots = new mpc_t[*vdeg+1];
	// init_rk1_mpc(*roots, *vdeg+1);
	for (int k=0; k<=*vdeg; k++) {
		mpc_init3((*roots)[k], wp2, wp2);
	}
	*tols = new mpfr_t[*vdeg+1];
	// init_rk1_mpfr(*tols, *vdeg+1);
	for (int k=0; k<=*vdeg; k++) {
		mpfr_init2((*tols)[k], wp2);
	}
	int mul;
	// cout << "call polynomial root finder and read multiplicity" << endl;
	// cout << "coeffs:" << endl;
	// print_poly((*pf).coeffs, *vdeg);
	if (*vdeg > 0) {
		mp_zroots_impr_out_mul(
			*roots, *tols, &mul,
			(*pf).coeffs, *vdeg, NULL,
			1, wp2, mpfr_tol, curr_mul
		);
	} else {
		mul = 1;
	}
	// cout << "mul = " << mul << endl;
	// cout << "1st attempt roots:" << endl;
	// print_poly(*roots, *vdeg);
	// cout << "1st attempt tols:" << endl;
	// for (int k=1; k<=*vdeg; k++) {
	// 	mpfr_out_str(stdout, 10, 0, (*tols)[k], MPFR_RNDN); cout << endl;
	// }
	// cout << "fake root:" << endl;
	// print_mpc(&(*roots)[0]); cout << endl;

	//////
	// ITERATIVE CALLS
	//////
	// cout << "ITERATIVE CALLS" << endl;
	// iteratively call itself with increased precision
	// if detected multiplicity is greater than the current one,
	// return otherwise
	int p_wp2;
	mpfr_t p_mpfr_tol;
	mpfr_init(p_mpfr_tol);
	mpc_t *p_roots;
	mpfr_t *p_tols;
	mpc_t fake_root;
	mpc_init3(fake_root, wp2, wp2);
	if (mul <= curr_mul) {
		// for (int k=1; k<=ldeg; k++) {
		//   mpc_set_ui(inroots[k], 0, MPFR_RNDN);
		//   mpfr_set(intols[k], mpfr_tol, MPFR_RNDN);
		// }
		// cout << "mul is enough" << endl;
		free(*s_tree);
		*s_tree = (char*)malloc(sizeof(char)*(strlen(s_tree_copy)+1));
		strcpy(*s_tree, s_tree_copy);
		return;
	} else {
		// save fake root
		mpc_set(fake_root, (*roots)[0], MPFR_RNDN);
		// increase precision
		p_wp2 = ((double) mul * wp2)/curr_mul;
	}

	// cout << "wp2 = " << wp2 << endl;
	// cout << "p_wp2 = " << p_wp2 << endl;
	mpfr_set_prec(p_mpfr_tol, p_wp2);
	// cout << "mpfr_tol = "; mpfr_out_str(stdout, 10, 0, mpfr_tol, MPFR_RNDN); cout << endl;
	mpfr_rootn_ui(p_mpfr_tol, mpfr_tol, curr_mul, MPFR_RNDN);
	// cout << "p_mpfr_tol = "; mpfr_out_str(stdout, 10, 0, p_mpfr_tol, MPFR_RNDN); cout << endl;
	mpfr_pow_si(p_mpfr_tol, p_mpfr_tol, mul, MPFR_RNDN);
	// cout << "p_mpfr_tol = "; mpfr_out_str(stdout, 10, 0, p_mpfr_tol, MPFR_RNDN); cout << endl;

	// call itself with icreased precision
	// cout << "call itself with icreased precision" << endl;
	free(s_tree_copy);
	// if ((*pf).coeffs) {
	// 	delete[] (*pf).coeffs;
	// 	(*pf).coeffs = NULL;
	// }
	poly_frac_numer_roots_from_tree(
		pf, &p_roots, &p_tols, vdeg, ldeg,
		in_nroots, in_roots, in_tols,
		s_tree, sep,
		in,
		kin, skip_inv, is_mass,
		p_wp2, p_mpfr_tol, 1, mul
	);

	// lower back precision of the output
	// cout << "lower back precision of the output" << endl;
	mpfr_t EPSS_red;
	mpfr_init2(EPSS_red, wp2);
	mpfr_rootn_ui(EPSS_red, mpfr_tol, mul, MPFR_RNDN);
	mpfr_mul_d(EPSS_red, EPSS_red, 1e10, MPFR_RNDN);
	mpfr_set_prec(tmp_mpfr, p_wp2);
	mpc_set_prec(tmpc, p_wp2);
	for (int k=1; k<=*vdeg; k++) {
		mpc_set((*roots)[k], p_roots[k], MPFR_RNDN);
		mpc_sub(tmpc, (*roots)[k], fake_root, MPFR_RNDN);
		mpc_abs(tmp_mpfr, tmpc, MPFR_RNDN);
		// mpfr_sub(tmp_mpfr, tmp_mpfr, mpfr_tol, MPFR_RNDN);
		// mpfr_abs(tmp_mpfr, tmp_mpfr, MPFR_RNDN);
		if (0 && mpfr_less_p(tmp_mpfr, EPSS_red)) {
			mpfr_rootn_ui((*tols)[k], p_tols[k], mul, MPFR_RNDN);
		} else {
			mpfr_set((*tols)[k], p_tols[k], MPFR_RNDN);
		}
	}
	// for (int k=1; k<=ldeg; k++) {
	// 	mpc_set_ui((*roots)[k], 0, MPFR_RNDN);
	// 	mpfr_set((*tols)[k], mpfr_tol, MPFR_RNDN);
	// }
}


// #pair: poly_frac_numer_roots_from_tree
void poly_mpq_find_roots(
	// OUTPUT
	int *found_mul_roots,
	int *vdeg, mpc_t **roots, mpfr_t **tols,
	// IN-OUT
	// int *in_nroots, mpc_t **in_roots, mpfr_t **in_tols,
	// INPUT
	mpq_t *coeffs_mpq, int deg,
	int wp2, mpfr_t mpfr_tol, int incr_prec, int curr_mul,
	int defl_nroots, int *defl_mults, mpc_t *defl_roots
) {
	if (dbg) {
		cout << "internal mpfr_tol = ";
		mpfr_out_str(stdout, 10, 0, mpfr_tol, MPFR_RNDN); cout << endl;
	}

	//////
	// GET WP COEFFICIENTS
	//////
	mpc_t *coeffs = new mpc_t[deg+1];
	// find ldeg
	int ldeg = 0;
	for (int k=0; k<=deg; k++) {
		// cout << "sign of "; mpq_out_str(stdout, 10, coeffs_mpq[k]); cout << endl;
		if (mpq_sgn(coeffs_mpq[k]) == 0) {
			ldeg++;
		}	else {
			break;
		}
	}
	// cout << "ldeg = " << ldeg << endl;
	*vdeg = deg-ldeg;
	*roots = new mpc_t[*vdeg+1];
	for (int k=0; k<=*vdeg; k++) {
		mpc_init3((*roots)[k], wp2, wp2);
	}
	*tols = new mpfr_t[*vdeg+1];
	for (int k=0; k<=*vdeg; k++) {
		mpfr_init2((*tols)[k], wp2);
	}

	// convert coeffs
	for (int k=ldeg; k<=deg; k++) {
		mpc_init3(coeffs[k-ldeg], wp2, wp2);
		mpc_set_q(coeffs[k-ldeg], coeffs_mpq[k], MPFR_RNDN);
		// cout << "k = " << k << ": "; print_mpc(&coeffs[k-ldeg]); cout << endl;
	}

	// perform deflation if necessary
	// cout << "vdeg = " << *vdeg << endl; fflush(stdout);
	int *defl_root_lab = new int[*vdeg];
	int count = 0;
	for (int k=1; k<defl_nroots; k++) {
		for (int m=0; m<defl_mults[k]; m++) {
			// cout << "count = " << count << endl; fflush(stdout);
			deflation(coeffs, *vdeg-count, &defl_roots[k], wp2);
			// (*vdeg)--;
			defl_root_lab[count] = k;
			count++;
		}
		defl_mults[k] = 0;
	}

	if (dbg) {
		cout << "remnant coeffs after deflation:" << endl;
		print_poly(coeffs, *vdeg-count);
	}

	//////
	// FIRST CALL TO ROOT FINDER
	//////
	int mul;
	if (dbg) cout << "FIRST CALL TO ROOT FINDER" << endl;
	if (*vdeg > 0) {
		mp_zroots_impr_out_mul(
			*roots+count, *tols+count, &mul,
			coeffs, *vdeg-count, NULL,
			1, wp2, mpfr_tol, curr_mul
		);
	} else {
		mul = 1;
	}
	if (dbg) {
		cout << "mul = " << mul << endl;
		cout << "roots:" << endl;
		print_poly(*roots+count, *vdeg-count);
	}

	// // FILL FIRST POSITIONS WITH DEFLATION ROOTS
	// // (cannot do that before because a NaN is inserted by convention
	// // in first position, so it would override the last deflation root)
	// count = 0;
	// for (int k=1; k<defl_nroots; k++) {
	// 	for (int m=0; m<defl_mults[k]; m++) {
	// 		perror("I am where I am not supposed to");
	// 		printf("I am where I am not supposed to");
	// 		mpc_set((*roots)[1+count++], defl_roots[k], MPFR_RNDN);
	// 	}
	// }

	//////
	// ITERATIVE CALLS
	//////
	int p_wp2;
	if (mul <= curr_mul) {
		*found_mul_roots = 0;
		delete[] defl_root_lab;
		mpc_rk1_clear(coeffs, deg+1);
		delete[] coeffs;
		return;
	} else {
		// increase precision
		p_wp2 = ((double) mul * wp2)/curr_mul;
	}
	mpfr_t p_mpfr_tol;
	mpfr_init(p_mpfr_tol);
	mpc_t *p_roots;
	mpfr_t *p_tols;

	mpfr_set_prec(p_mpfr_tol, p_wp2);
	mpfr_rootn_ui(p_mpfr_tol, mpfr_tol, curr_mul, MPFR_RNDN);
	mpfr_pow_si(p_mpfr_tol, p_mpfr_tol, mul, MPFR_RNDN);

	poly_mpq_find_roots(
		found_mul_roots,
		vdeg, &p_roots, &p_tols,
		// in_nroots, in_roots, in_tols,
		coeffs_mpq, deg,
		p_wp2, p_mpfr_tol, 1, mul,
		defl_nroots, defl_mults, defl_roots
	);

	// lower back precision of the output
	int recount = 0, idx = 0, found_defl = 0;
	if (dbg ) cout << "count = " << count << endl;
	for (int k=1; k<=*vdeg; k++) {
		if (dbg) {
		cout << "k = " << k << endl;
		cout << "recount = " << recount << endl;
		cout << "idx = " << idx << endl;
		}
		if (recount < count) {
			// first check deflation roots
			found_defl = 0;
			for (int kp=0; kp<count; kp++) {
				if (mpc_equal_within_tol(defl_roots[defl_root_lab[kp]], p_roots[k])) {
					if (dbg) cout << "found defl at kp = " << kp << endl;
					mpc_set((*roots)[1+recount], p_roots[k], MPFR_RNDN);
					mpfr_set((*tols)[1+recount], p_tols[k], MPFR_RNDN);
					recount++;
					found_defl = 1;
					break;
				}
			}
			if (found_defl) {
				continue;
			}
		}
		mpc_set((*roots)[1+count+idx], p_roots[k], MPFR_RNDN);
		mpfr_set((*tols)[1+count+idx], p_tols[k], MPFR_RNDN);
		idx++;
		// mpc_set((*roots)[k], p_roots[k], MPFR_RNDN);
		// mpfr_set((*tols)[k], p_tols[k], MPFR_RNDN);
	}

	*found_mul_roots = 1;
	// delete[] defl_root_lab;

	// FREE
	mpfr_clear(p_mpfr_tol);
	mpc_rk1_clear(p_roots, *vdeg+1);
	delete[] p_roots;
	mpfr_rk1_clear(p_tols, *vdeg+1);
	delete[] p_tols;
}


void generate_PS_pf(
	// OUTPUT
	struct poly_frac *out,
	// INPUT
	char ***kin, int skip_inv, int is_mass, int s, int nroots,
	// int ninvs, ***char kin, int *skip_inv, int *is_mass
	int wp2
) {
	(*out).nroots = nroots;
	if ((*out).mults) {
		delete [] (*out).mults;
	}
	(*out).mults = new int[nroots];
	for (int k=0; k<nroots; k++) {
		(*out).mults[k] = 0;
	}

	if (s == 0) {
		//////
		// EPSILON FROM SPACE-TIME DIMENSION
		//////
		poly_frac_set_ui_wp2(wp2, out, 4, nroots);

		mpc_t eps;
		mpc_init3(eps, wp2, wp2);
		mpc_set_str_rat(&eps, kin[0][s]);
		mpc_mul_ui(eps, eps, 2, MPFR_RNDN);
		mpc_sub((*out).coeffs[0], (*out).coeffs[0], eps, MPFR_RNDN);
		mpc_clear(eps);
		return;
	}

	//////
	// INVARIANTS
	//////
	if (skip_inv) {
		// only store initial point
		mpc_t ini;
		mpc_init3(ini, wp2, wp2);
		mpc_set_str_rat(&ini, kin[0][s]);
		poly_frac_set_mpc_wp2(wp2, out, &ini, nroots);
		mpc_clear(ini);
		if (out->num_vdeg == -1) {
			out->coeffs = new mpc_t[2];
			mpc_init3(out->coeffs[0], wp2, wp2);
			mpc_set_ui(out->coeffs[0], 0, MPFR_RNDN);
			mpc_init3(out->coeffs[1], wp2, wp2);
			mpc_set_ui(out->coeffs[1], 0, MPFR_RNDN);
		} else if (is_mass) {
			// linearize squared masses
			mpc_sqrt(out->coeffs[0], out->coeffs[0], MPFR_RNDN);
		}
	} else {
		// // AVOID COPY TO INTERMEDIATE ARRAY IN MOST CASES
		// if (kin[0][s][0] == "0" && strlen(kin[0][s][0]) == 1) {
		// 	// only store final point
		// 	mpc_t *fin;
		// 	mpc_init3(fin, wp2, wp2);
		// 	mpc_set_str_rat(&fin, kin[1][s]);
		// 	poly_frac_set_mpc_wp2(wp2, out, &fin, nroots);	
		// } else {
		// 	// store both initial and final points
		// 	(*out).num_deg = 1;
		// 	(*out).num_vdeg = 1;
		// 	(*out).den_deg = 0;
		// if ((*out).coeffs) {
    // 	delete[] (*out).coeffs;
  	// }
		// 	(*out).coeffs = new mpc_t[2];

		// 	mpc_init3((*out).coeffs[0], wp2, wp2);
		// 	mpc_set_str_rat(&(*out).coeffs[0], kin[1][s]);
		// 	mpc_((*out).coeffs[0], &fin, nroots);	
		// }
		mpc_t *PS_if = new mpc_t[2];
		mpc_init3(PS_if[0], wp2, wp2);
		mpc_set_str_rat(&PS_if[0], kin[0][s]);
		mpc_init3(PS_if[1], wp2, wp2);
		mpc_set_str_rat(&PS_if[1], kin[1][s]);
		
		if (is_mass) {
			// linearize squared masses
			mpc_sqrt(PS_if[0], PS_if[0], MPFR_RNDN);
			mpc_sqrt(PS_if[1], PS_if[1], MPFR_RNDN);
		}
		mpc_sub(PS_if[1], PS_if[1], PS_if[0], MPFR_RNDN);
		(*out).nroots = nroots;
		// if ((*out).mults) {
		// 	delete [] (*out).mults;
		// }
		// (*out).mults = new int[nroots];
		// for (int k=0; k<nroots; k++) {
		// 	(*out).mults[k] = 0;
		// }
		poly_frac_set_coeffs_wp2(out, PS_if, 1, wp2);
		mpc_clear(PS_if[0]); mpc_clear(PS_if[1]);
		delete[] PS_if;
		// poly_frac_print(out);
		// if (is_mass) {
		// 	poly_frac_pow_ui_wp2(wp2, out, out, 2);
		// 	// poly_frac_print(out);
		// }
	}

}


// # forward declared above
void decode_tree_pf(
	// OUTPUT
	struct poly_frac *out,
	// IN-OUT
	int *nroots, mpc_t **roots, mpfr_t **tols,
	// INPUT
	struct poly_frac *in,
	char ***kin, int *skip_inv, int *is_mass,
	char **s_tree, char *sep,
	int wp2, mpfr_t mpfr_tol, int incr_prec
)  {
	// cout << endl; cout << "ENTER decode_tree" << endl;
	// cout << "input tree:" << endl;
	// cout << (*s_tree) << endl;
	// cout << "input address: " << out << endl;
	char *s_tree_copy = (char*)malloc(sizeof(char)*(strlen((*s_tree))+1));
	// char *s_tree_copy2 = (char*)malloc(sizeof(char)*(strlen((*s_tree))+1));
	strcpy(s_tree_copy, (*s_tree));
	// strcpy(s_tree_copy2, (*s_tree));
	int stride = 0;
	int n, pow, idx;
	char *operation;
	struct poly_frac *inner_pf;
	char *token = strtok((*s_tree), (char*) sep);
	stride += strlen(token) + strlen(sep);
	while (token) {
		// printf("token: %s\n", token);
		if (strcmp(token, "o") == 0) {
			// cout << "START operation" << endl;
			token = strtok(NULL, (char*) sep);
			// printf("token: %s\n", token);
			stride += strlen(token) + strlen(sep);
			operation = (char*) malloc(sizeof(char)*(strlen(token)+1));
			strcpy(operation, token);
			// cout << "operation = " << operation << endl;

			// for power and division, n=1, otherwise it is written in the string
			if (strcmp(operation, "p") == 0 || strcmp(operation, "i") == 0) {
				n = 1;
				token = strtok(NULL, (char*) sep);
				// printf("token: %s\n", token);
				stride += strlen(token) + strlen(sep);
				pow = atoi(token);
			} else if (strcmp(operation, "/") == 0) {
				n = 1;
				pow = -1;
			} else if (strcmp(operation, "-") == 0) {
				n = 1;
			} else {
				token = strtok(NULL, (char*) sep);
				// printf("token: %s\n", token);
				stride += strlen(token) + strlen(sep);
				n = atoi(token);
			}
			// cout << "n = " << n << endl;

			inner_pf = new struct poly_frac[n];
			if (strcmp(operation, "i") == 0 || strcmp(operation, "/") == 0) {
				//////
				// INVERSE OF POLY_FRAC
				//////
				// cout << "COMPUTE INVERSE OF POLY_FRAC" << endl;

				// OUTPUT DENOMINATOR FROM INPUT NUMERATOR
				// cout << "OUTPUT DENOMINATOR FROM INPUT NUMERATOR" << endl;
				// cout << "find numerator roots from tree" << endl;
				free(*s_tree);
				*s_tree = (char*)malloc(sizeof(char)*(strlen(s_tree_copy + stride)+1));
				strcpy((*s_tree), s_tree_copy + stride);
				poly_frac_build(&inner_pf[0]);
				mpc_t *num_roots;
				mpfr_t *num_tols;
				int num_deg, num_vdeg, num_ldeg;
				poly_frac_numer_roots_from_tree(
					&inner_pf[0], &num_roots, &num_tols,
					&num_vdeg, &num_ldeg,
					nroots, roots, tols,
					s_tree, sep,
					in,
					kin, skip_inv, is_mass,
					wp2, mpfr_tol, incr_prec, 1
				);
				num_deg = inner_pf[0].num_deg;
				// cout << "inner poly_frac:" << endl;
				// poly_frac_print(&inner_pf[0]);
				// cout << "inner nroots = " << inner_pf[0].nroots << endl;
				// cout << "num vdeg, num ldeg = " << num_vdeg << ", " << num_ldeg << endl;				
				// cout << "numerator roots:" << endl;
				// print_poly(num_roots, inner_pf[0].num_vdeg);

				// update global roots and assign labels
				// cout << "update global roots and assign labels" << endl;
				int *root_prof = new int[num_deg];
				int num_new_roots = 0;
				mpc_t *new_roots = new mpc_t[num_vdeg];
				mpfr_t *new_tols = new mpfr_t[num_vdeg];
				for (int k=0; k<num_ldeg; k++) {
					root_prof[k] = 0;
				}
				// cout << "old nroots = " << *nroots << endl;
				// cout << "old roots:" << endl;
				// print_poly(*roots, *nroots-1);

				// #2BD: replace update_new_roots, and append function in order to
				// initialize new mpc_t variables at local precision
				// (needed for denominators in denominators)

				for (int k=1; k<=num_vdeg; k++) {
				// cout << "analyze k = " << k << endl;
				// cout << "tmp root:" << endl;
				
				// cout << "tmp tol:" << endl;
				// mpfr_out_str(stdout, 10, 0, tmp_tols[k], MPFR_RNDN); cout << endl;
					update_new_roots(
						&root_prof[num_ldeg + k-1], &num_new_roots, new_roots, new_tols,
						nroots, *roots, *tols,
						&num_roots[k] , &num_tols[k]
					);
				}
				// cout << "num new roots = " << num_new_roots << endl;
				// if (num_new_roots > 0) {
				// 	cout << "found new roots:" << endl;
				// 	print_poly(new_roots, num_new_roots-1);
				// }

				mpc_rk1_append(
					nroots, roots,
					num_new_roots, new_roots
				);
				mpfr_rk1_append(
					nroots, tols,
					num_new_roots, new_tols
				);
				*nroots += num_new_roots;
				(*out).nroots = *nroots;
				(*out).den_deg = num_deg;

				// assign multiplicities
				// cout << "assign multiplicities" << endl;
				root_prof_to_poly_frac_den(
					out,
					root_prof, num_deg,
					*nroots
				);

				// OUTPUT NUMERATOR FROM INPUT DENOMINATOR
				// cout << "OUTPUT NUMERATOR FROM INPUT DENOMINATOR" << endl;
				poly_frac_den_to_num(out, &inner_pf[0], *roots, wp2);
				// make polynomial monic
				poly_frac_div_mpc(out, out, &inner_pf[0].coeffs[inner_pf[0].num_vdeg]);
				// cout << "inverted poly_frac:" << endl;
				// poly_frac_print(out);
				// exit(0);
				free(s_tree_copy);
				s_tree_copy = (char*)malloc(sizeof(char)*(strlen(*s_tree)+1));
				strcpy(s_tree_copy, (*s_tree));
				stride = 0;
				// cout << "RAISE TO POWER" << endl;
				// *s_tree = (char*)malloc(sizeof(char)*(strlen(s_tree_copy + stride)+1));
				// strcpy((*s_tree), s_tree_copy + stride);
				// cout << *s_tree << endl;
				// exit(0);
				// token = strtok(NULL, (char*) sep);
				// stride += strlen(token) + strlen(sep);
				// token = strtok(NULL, (char*) sep);
				// stride += strlen(token) + strlen(sep);
				// int pow = -atoi(token);
				// cout << "pow = " << pow << endl;
				// cout << *s_tree << endl;
				free(*s_tree);
				*s_tree = (char*)malloc(sizeof(char)*(strlen(s_tree_copy + stride)+1));
				strcpy((*s_tree), s_tree_copy + stride);
				// cout << "s_tree = " << *s_tree << endl;
				poly_frac_pow_ui_wp2(wp2, out, out, -pow);
				// cout << "poly_frac raised to power:" << endl;
				// poly_frac_print(out);
			} else {
				for (int i=0; i<n; i++) {
					if (dbg) cout << "i = " << i << endl;
					poly_frac_build(&inner_pf[i]);
					free(*s_tree);
					*s_tree = (char*)malloc(sizeof(char)*(strlen(s_tree_copy + stride)+1));
					// cout << "size = " << strlen(s_tree_copy + stride)+1 << endl;
					strcpy((*s_tree), s_tree_copy + stride);
					decode_tree_pf(
						&inner_pf[i],
						nroots, roots, tols,
						in,
						kin, skip_inv, is_mass,
						s_tree,
						sep, wp2, mpfr_tol, incr_prec
					);
					free(s_tree_copy);
					s_tree_copy = (char*)malloc(sizeof(char)*(strlen(*s_tree)+1));
					strcpy(s_tree_copy, (*s_tree));
					stride = 0;
				}
				// update nroots of inner pf
				// (nroots could have been changed in the process of getting
				// inner polynomials, therefore it has to be changed here,
				// AFTER every inner element is obtained)
				for (int i=0; i<n; i++) {
					poly_frac_roots_append_zero(&inner_pf[i], *nroots);
				}
				if (strcmp(operation, "+") == 0) {
					// cout << "adding " << n << " output" << endl;
					poly_frac_set_zero(out);
					for (int i=0; i<n; i++) {
						poly_frac_add_pf_wp2(wp2, out, out, &inner_pf[i], *roots, NULL);
					}
					// cpoly_print(out); cout << endl;
				} else if (strcmp(operation, "-") == 0) {
					poly_frac_set_pf_wp2(wp2, out, &inner_pf[0]);
					poly_frac_neg(out);
				} else if (strcmp(operation, "*") == 0) {
					// cout << "multiplying " << n << " output" << endl;
					poly_frac_set_ui_wp2(wp2, out, 1, *nroots);
					for (int i=0; i<n; i++) {
						// cout << "mult i = " << i << endl;
						// cout << "factor 1:" << endl;
						// poly_frac_print(out);
						// cout << "nroots = " << (*out).nroots << endl;
						// cout << "factor 2:" << endl;
						// poly_frac_print(&inner_pf[i]);
						// cout << "nroots = " << inner_pf[i].nroots << endl;
						poly_frac_mul_pf_wp2(wp2, out, out, &inner_pf[i]);
						// cout << "product:" << endl;
						// poly_frac_print(out);
						// cout << "nroots = " << (*out).nroots << endl;
					}
					// cpoly_print(out); cout << endl;
				} else if (strcmp(operation, "p") == 0) {
					// cout << "raising to power" << endl;
					// int pow = mpfr_get_si(mpc_realref(inner_pf[1].coeffs[0]), MPFR_RNDN);
					// cout << "pow = " << pow << endl;
					poly_frac_pow_ui_wp2(wp2, out, &inner_pf[0], pow);
					// cout << "output:" << endl;
					// poly_frac_print(out);
					// cout << *s_tree << endl;
				}
			}
			// cout << "END operation" << endl;
			// cout << "output:" << endl;
			// poly_frac_print(out);
			// cout << "nroots = " << (*out).nroots << endl;
			// cout << "address = " << out << endl;
			poly_frac_rk1_free(inner_pf, n); delete[] inner_pf;
			free(s_tree_copy);
			return;
		} else {
			// cout << "START ELSE" << endl;
			// if (strcmp(token, "p") == 0) {
			// 	token = strtok(NULL, (char*) sep);
			// 	stride += strlen(token) + strlen(sep);
			// 	idx = atoi(token);
			// 	token = strtok(NULL, (char*) sep);
			// 	stride += strlen(token) + strlen(sep);
			// 	pow = atoi(token);
			// 	(*s_tree) = (char*)malloc(sizeof(char)*(strlen(s_tree_copy + stride)+1));
			// 	strcpy((*s_tree), s_tree_copy + stride);
			// 	cout << "SET power" << endl;
			// 	cout << "pow = " << pow << endl;
			// 	// cpoly_print(&in[idx]);
			// 	cpoly_pow_ui(out, &in[idx], pow);
			// 	// cpoly_print(out);
			// 	cout << endl;
			if (strcmp(token, "s") == 0) {
				token = strtok(NULL, (char*) sep);
				stride += strlen(token) + strlen(sep);
				idx = atoi(token);
				free(*s_tree);
				(*s_tree) = (char*)malloc(sizeof(char)*(strlen(s_tree_copy + stride)+1));
				strcpy((*s_tree), s_tree_copy + stride);
				// cout << "SET symbol" << endl;
				if (!incr_prec) {
					// cout << "standard prec" << endl;
					poly_frac_set_pf_wp2(wp2, out, &in[idx]);
				} else {
					// read phase-space from string with increased precision
					// cout << "increased prec PS" << endl;
					// cout << "increased wp2 = " << wp2 << endl;
					// generate_PS_pf(out, kin, skip_inv[idx-1], is_mass[idx-1], idx, *nroots, wp2);
					// cout << "idx = " << idx << endl;
					// cout << "skip_inv[idx] = " << skip_inv[idx] << endl;
					// cout << "is_mass[idx] = " << is_mass[idx] << endl;
					if (idx == 0) {
						generate_PS_pf(out, kin, 0, 0, idx, *nroots, wp2);
					} else if (idx > 0) {
						generate_PS_pf(out, kin, skip_inv[idx], is_mass[idx], idx, *nroots, wp2);
					}
					// if (idx == 0 || skip_inv[idx-1] == 1)  {
					// 	mpc_t PS_ini;
					// 	mpc_init3(PS_ini, wp2, wp2);
					// 	mpc_set_str_rat(&PS_ini, kin[0][idx]);
					// 	if (is_mass[idx-1]) {
					// 		mpc_sqrt(PS_ini, PS_ini, MPFR_RNDN);
					// 	}
					// 	poly_frac_set_mpc_wp2(wp2, out, &PS_ini, *nroots);
					// 	if (is_mass[idx-1]) {
					// 		poly_frac_pow_ui_wp2(wp2, out, out, 2);
					// 	}
					// } else {
					// 	mpc_t *PS_if = new mpc_t[2];
					// 	mpc_init3(PS_if[0], wp2, wp2);
					// 	mpc_set_str_rat(&PS_if[0], kin[0][idx]);
					// 	mpc_init3(PS_if[1], wp2, wp2);
					// 	mpc_set_str_rat(&PS_if[1], kin[1][idx]);
					// 	if (is_mass[idx-1]) {
					// 		cout << "sqrt of mass" << endl;
					// 		mpc_sqrt(PS_if[0], PS_if[0], MPFR_RNDN);
					// 		mpc_sqrt(PS_if[1], PS_if[1], MPFR_RNDN);
					// 	}
					// 	mpc_sub(PS_if[1], PS_if[1], PS_if[0], MPFR_RNDN);
					// 	(*out).nroots = *nroots;
					// 	if ((*out).mults) {
					// 		delete [] (*out).mults;
					// 	}
					// 	(*out).mults = new int[*nroots];
					// 	for (int k=0; k<*nroots; k++) {
					// 		(*out).mults[k] = 0;
					// 	}
					// 	poly_frac_set_coeffs(out, PS_if, 1);
					// 	if (is_mass[idx-1]) {
					// 		poly_frac_pow_ui_wp2(wp2, out, out, 2);
					// 		poly_frac_print(out);
					// 	}
					// }
				}
				// poly_frac_print(out);
			} else if (strcmp(token, "n") == 0) {
				token = strtok(NULL, (char*) sep);
				// cout << "SET numeric" << endl;
				// cout << token << endl;
				mpc_t number;
				mpc_init3(number, wp2, wp2);
				mpc_set_str_rat(&number, token);
				poly_frac_set_mpc_wp2(wp2, out, &number, *nroots);
				mpc_clear(number);
				stride += strlen(token) + strlen(sep);
				free(*s_tree);
				(*s_tree) = (char*)malloc(sizeof(char)*(strlen(s_tree_copy + stride)+1));
				strcpy((*s_tree), s_tree_copy + stride);
				// cpoly_print(out); cout << endl;
			}
			// cout << "output:" << endl;
			// poly_frac_print(out);
			// cout << "nroots = " << (*out).nroots << endl;
			// cout << "address = " << out << endl;
			free(s_tree_copy);
			return;
		}
		token = strtok(NULL, (char*) sep);
		stride += strlen(token) + strlen(sep);
	}
	free(s_tree_copy);
}
