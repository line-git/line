#include <stdio.h>
#include <stdlib.h>
#include <string.h>
// #include <ctype.h>
#include "mpc.h"
#include "global_vars.h"

//////
// 1ST IMPLEMENTATION
//////

// Structure to represent an expression tree node
struct node {
	int n; // number of args
	char op; // operation ('+' for addition, '-' for subtraction, '*' for multiplication, '/' for division, '^' for power raising)
	struct node *args; // pointer to the operand node
	int number;
};


// Function to print the expression tree
void print_expression_tree(struct node* root, const char **symbols) {
	if (root->op == 's') {
		if (symbols == 	NULL) {
			printf("(%c %d)", root->op, (*root).number);
		} else {
			printf("(%c %s)", root->op, symbols[(*root).number]);
		}
		return;
	} else if(root->op == 'n') {
		printf("(%c %d)", root->op, (*root).number);
		return;
	} else	if (root->op == 'I') {
		printf("(%c)", root->op);
		return;
	} else {
		printf("(%c ", root->op);
		for (int i = 0; i < root->n; i++) {
			print_expression_tree(&(root->args[i]), symbols);
			if (i < root->n - 1)
				printf(" ");
		}
		printf(")");
	}
}


// Function to create a new node
void create_node(
	// OUTPUT
	struct node *new_node,
	// INPUT
	int n, char op
) {
	// printf("creating node with n = %d, op = %c\n", n, op);
	new_node->n = n;
	new_node->op = op;
	if (op == 's' || op == 'n') {
		new_node->args = NULL;
	} else {
		new_node->args = (struct node*)malloc(n*sizeof(struct node));
	}

	return;
}


void find_precedence(
	// OUTPUT
	int *num_lowest_precedence, int *lowest_precedence_index, int *lowest_precedence,
	// INPUT
	char *expression
) {
	*num_lowest_precedence = 0; // Number of operators with the lowest precedence
	int parenthesis_count = 0;
	for (int i = 0; expression[i] != '\0'; i++) {
		// printf("%c", expression[i]);
		if (expression[i] == '(') {
			parenthesis_count++;
		} else if (expression[i] == ')') {
			parenthesis_count--;
		} else if (parenthesis_count == 0 && strchr("^*/+-", expression[i])) {
			int precedence;
			switch (expression[i]) {
				case '^':
					precedence = 4;
					break;
				case '*':
				case '/':
					precedence = 3;
					break;
				case '+':
					precedence = 2;
					break;
				case '-':
					if (i > 0) {
						if (expression[i-1] == '^') {
							precedence = 5;
							break;
						}
					}
					precedence = 2;
					break;
				default:
					precedence = 1; // Symbols have the lowest precedence
			}
			if (precedence <= *lowest_precedence) {
				// printf("\n");
				if (precedence < *lowest_precedence) {
					// printf("new lowest precedence\n");
					*num_lowest_precedence = 0;
					*lowest_precedence = precedence;
				}
				// printf("add lowest precedence operation index: %d\n", i);
				lowest_precedence_index[*num_lowest_precedence] = i;
				(*num_lowest_precedence)++;
			}			
		}
	}
	// printf("\n");
}


// Function to parse the expression recursively
void parse_expression(
	// OUTPUT
	struct node *root,
	// INPUT
	char* expression,
	const char** symbols, int num_symbols
) {
	printf("\ninput expression: %s\n", expression);

	// Remove white spaces from expression
	char* p = expression;
	char* q = expression;
	while (*q != '\0') {
		if (*q != ' ') {
			*p = *q;
			p++;
		}
		q++;
	}
	*p = '\0';

	// Remove outer parentheses if present
	int n = strlen(expression);
	int parenthesis_count = 0;

	// Check if the expression starts with '('
	while (expression[0] == '(') {
		// Scan the string starting from the character right after the '('
		int i;
		for (i = 1; i < n; i++) {
			if (expression[i] == '(') {
				parenthesis_count++;
			} else if (expression[i] == ')') {
				parenthesis_count--;
			}

			if (parenthesis_count == -1) {
				break;
			}
		}

		// printf("parenthesis_count = %d\n", parenthesis_count);
		if (i == n-1) {
			if (expression[n-1] != ')') {
				fprintf(stderr, "Error: Unmatched parenthesis.\n");
				exit(1);
			} else {
				expression[n-1] = '\0';
				expression++;
				printf("unnecessary outer parenthesis removed\n");
				continue;
			}
		} else {
			break;
		}
	}

	// Base case: expression is a symbol
	for (int i = 0; i < num_symbols; i++) {
		if (strcmp(expression, symbols[i]) == 0) {
			printf("expression is a symbol\n");
			create_node(root, 0, 's');
			(*root).args = NULL;
			(*root).number = i;
			return;
		}
	}

	int lowest_precedence = 1000; // Initialize to a large value
	int lowest_precedence_index[100]; // Maximum number of operators with the same precedence
	int num_lowest_precedence;
	find_precedence(
		&num_lowest_precedence, lowest_precedence_index, &lowest_precedence,
		expression
	);
	// int num_lowest_precedence = 0; // Number of operators with the lowest precedence
	// parenthesis_count = 0;
	// for (int i = 0; expression[i] != '\0'; i++) {
	// printf("%c", expression[i]);
	//     if (expression[i] == '(') {
	//         parenthesis_count++;
	//     } else if (expression[i] == ')') {
	//         parenthesis_count--;
	//     } else if (parenthesis_count == 0 && strchr("^*/+-", expression[i])) {
	//         int precedence;
	//         switch (expression[i]) {
	//             case '^':
	//                 precedence = 4;
	//                 break;
	//             case '*':
	//             case '/':
	//                 precedence = 3;
	//                 break;
	//             case '+':
	//             case '-':
	//                 precedence = 2;
	//                 break;
	//             default:
	//                 precedence = 1; // Symbols have the lowest precedence
	//         }
	//         if (precedence <= lowest_precedence) {
	// 		printf("\n");
	//             if (precedence < lowest_precedence) {
	// 			printf("new lowest precedence\n");
	//                 num_lowest_precedence = 0;
	//                 lowest_precedence = precedence;
	//             }
	// 		printf("add lowest precedence operation index: %d\n", i);
	//             lowest_precedence_index[num_lowest_precedence++] = i;
	//         }			
	//     }
	// }
	// printf("\n");

	// If no operator is found, then the operand is a number
	if (lowest_precedence == 1000) {
		create_node(root, 0, 'n');
		(*root).args = NULL;
		(*root).number = atoi(expression);
		printf("create numeric root with number = %d:\n", (*root).number);
		return;
	}

	// Split the expression into parts based on the lowest precedence operators
	printf("num_lowest_precedence = %d\n", num_lowest_precedence);
	if (num_lowest_precedence > 1) {
		// Parse each sub-string between lowest precedence operators
		// starting from the right-most one, so that we can read it's operation
		// character, use that info, and then replace that char with '\0'
		// when processing next sub-string

		struct node *root_args;
		if (lowest_precedence_index[0] != 0) {
			// if the first lowest precedence operation is not at the beginning
			// of the expression, the number of operation is the number of
			// sub-strings plus one
			create_node(root, num_lowest_precedence+1, '\0');
			root_args = (*root).args + 1;
		} else {
			create_node(root, num_lowest_precedence, '\0');
			root_args = (*root).args;
		}

		if (lowest_precedence == 2) {
			// create_node(root, num_lowest_precedence, '+');
			(*root).op = '+';

			// right-most sub-string
			char op = expression[lowest_precedence_index[num_lowest_precedence-1]];
			if (op == '+') {
				parse_expression(
					&root_args[num_lowest_precedence-1],
					expression + lowest_precedence_index[num_lowest_precedence-1] + 1,
					symbols, num_symbols
				);
			} else if (op == '-') {
				create_node(&root_args[num_lowest_precedence-1], 1, '-');
				parse_expression(
					&(*root).args[num_lowest_precedence-1].args[0],
					expression + lowest_precedence_index[num_lowest_precedence-1] + 1,
					symbols, num_symbols
				);			
			}

			// all the other sub-strings
			for (int i=num_lowest_precedence-2; i>=0; i--) {
				expression[lowest_precedence_index[i+1]] = '\0';
				op = expression[lowest_precedence_index[i]];
				if (op == '+') {
					// Apply parse_expression to the segment and assign it to the corresponding sub-string
					parse_expression(
						&root_args[i],
						expression + lowest_precedence_index[i] + 1,
						symbols, num_symbols
					);
				} else if (op == '-') {
					// Create a node with op='-', n=1, and arg[0] as the result of parse_expression for the sub-string
					create_node(&root_args[i], 1, '-');
					parse_expression(
						&root_args[i].args[0],
						expression + lowest_precedence_index[i] + 1,
						symbols, num_symbols
					);
				}
			}
		} else if (lowest_precedence == 3) {
			// create_node(root, num_lowest_precedence, '*');
			(*root).op = '*';

			// right-most sub-string
			char op = expression[lowest_precedence_index[num_lowest_precedence-1]];
			if (op == '*') {
				// printf("top ex: %s\n", expression + lowest_precedence_index[num_lowest_precedence-1] + 1);
				parse_expression(
					&root_args[num_lowest_precedence-1],
					expression + lowest_precedence_index[num_lowest_precedence-1] + 1,
					symbols, num_symbols
				);
			} else if (op == '/') {
				create_node(&root_args[num_lowest_precedence-1], 1, '/');
				parse_expression(
					&root_args[num_lowest_precedence-1].args[0],
					expression + lowest_precedence_index[num_lowest_precedence-1] + 1,
					symbols, num_symbols
				);			
			}

			// all the other sub-strings
			for (int i=num_lowest_precedence-2; i>=0; i--) {
				expression[lowest_precedence_index[i+1]] = '\0';
				printf("ex: %s\n", expression);
				op = expression[lowest_precedence_index[i]];
				if (op == '*') {
					// Apply parse_expression to the segment and assign it to the corresponding sub-string
					parse_expression(
						&root_args[i],
						expression + lowest_precedence_index[i] + 1,
						symbols, num_symbols
					);
				} else if (op == '/') {
					// Create a node with op='-', n=1, and arg[0] as the result of parse_expression for the sub-string
					create_node(&root_args[i], 1, '/');
					parse_expression(
						&root_args[i].args[0],
						expression + lowest_precedence_index[i] + 1,
						symbols, num_symbols
					);
				}				
			}
		}

		if (lowest_precedence_index[0] != 0) {
			// add the root corresponding to the first sub-string
			// printf("add the root corresponding to the first sub-string\n");
			// printf("ex: %s\n", expression);
			// printf("low prec idx = %d\n", lowest_precedence_index[0]);
			expression[lowest_precedence_index[0]] = '\0';
			parse_expression(
				&(*root).args[0],
				expression,
				symbols, num_symbols
			);
		}
		return;
	} else {
		int idx = lowest_precedence_index[0];
		char op = expression[idx];
		printf("op = %c, idx = %d\n", op, idx);
		if (idx == 0) {
			if (op == '+' || op == '-') {
				create_node(root, 1, op);
				parse_expression(
					&root->args[0],
					expression + 1,
					symbols, num_symbols
				);
			} else {
				fprintf(stderr, "Error: expression starting with an operation different from '+' or '-'.\n");
				exit(1);
			}
		} else {
			if (op == '-') {
				create_node(root, 2, '+');
				create_node(&root->args[1], 1, op);
				parse_expression(
					&(*root).args[1].args[0],
					expression + idx + 1,
					symbols, num_symbols
				);
				expression[idx] = '\0';
				parse_expression(
					&root->args[0],
					expression,
					symbols, num_symbols
				);
			} else if (op == '/'){
				create_node(root, 2, '*');
				create_node(&root->args[1], 1, op);
				parse_expression(
					&(*root).args[1].args[0],
					expression + idx + 1,
					symbols, num_symbols
				);
				expression[idx] = '\0';
				// printf("residual ex: %s\n", expression);
				parse_expression(
					&root->args[0],
					expression,
					symbols, num_symbols
				);
			} else {
				create_node(root, 2, op);
				parse_expression(
					&root->args[1],
					expression + idx + 1,
					symbols, num_symbols
				);
				expression[idx] = '\0';
				parse_expression(
					&root->args[0],
					expression,
					symbols, num_symbols
				);
			}
		}
		return;
	}

}


// Function to expand multiplications and divisions
void node_expand_mul_core(
	// IN-OUT
	struct node *nd
) {
	if ((*nd).op != '*') {
		return;
	}
	int count = 0;
	int *op_to_raise = (int*) malloc((*nd).n*sizeof(int));

	// find operation to raise up
	for (int i=0; i<(*nd).n; i++) {
		if ((*nd).args[i].op == '*') {
			op_to_raise[i] = (*nd).args[i].n;
			count += (*nd).args[i].n;
		} else {
			op_to_raise[i] = 0;
			count++;
		}
	}
	printf("count = %d\n", count);

	// raise up operations if necessary
	if (count>0) {
		struct node *holder = (*nd).args;
		(*nd).args = (struct node*) malloc(count*sizeof(struct node));
		count = 0;
		for (int i=0; i<(*nd).n; i++) {
			if (holder[i].op == '*') {
				for (int j=0; j<op_to_raise[i]; j++) {
					(*nd).args[count] = holder[i].args[j];
					count++;
				}
			} else {
				(*nd).args[count] = holder[i];
				count++;
			}
		}

		// update number of operations
		(*nd).n = count;
		
		// delete holder
		free(holder);
	}
	printf("Expression tree:\n");
	print_expression_tree(nd, NULL); printf("\n");
}


// Function to expand multiplications and divisions
// with bottom-up approach
void node_expand(
	// IN-OUT
	struct node *nd
) {
	if ((*nd).n == 0) {
		return;
	}

	// recursevely expand product terms in children nodes
	for (int i=0; i<(*nd).n; i++) {
		node_expand(&(*nd).args[i]);
	}

	node_expand_mul_core(nd);
}


// Function to expand multiplications and division
// with top-bottom approach
void node_expand_mul(
	// IN-OUT
	struct node *nd
) {
	int count;
	while (1) {
		// count operation to raise up
		count = 0;
		for (int i=0; i<(*nd).n; i++) {
			if ((*nd).args[i].op == '*') {
				count += (*nd).args[i].n;
			} else {
				count++;
			}
		}
		printf("count = %d\n", count);

		if (count == (*nd).n) {
			break;
		} else {
			node_expand_mul_core(nd);
		}
	}
}


void tree_mul_count_expanded_args(
	// IN-OUT
	int *count,
	// INPUT
	struct node *nd
) {
	if ((*nd).op != '*') {
		(*count)++;
		return;
	}
	
	for (int i=0; i<(*nd).n; i++) {
		tree_mul_count_expanded_args(count, &(*nd).args[i]);
	}

}


void tree_mul_fill_expanded_args(
	// IN-OUT
	struct node *args, int *count, 
	// INPUT
	struct node *nd
) {
	if ((*nd).op != '*') {
		args[*count] = *nd;
		(*count)++;
		return;
	}
	
	for (int i=0; i<(*nd).n; i++) {
		tree_mul_fill_expanded_args(args, count, &(*nd).args[i]);
	}
}


void tree_mul_expand_mul(
	// IN-OUT
	struct node *nd
) {
	struct node inner;
	inner.n = (*nd).n;
	inner.op = (*nd).op;
	inner.args = (*nd).args;

	int count = 0;
	tree_mul_count_expanded_args(&count, nd);
	(*nd).n = count;
	(*nd).args = (struct node*) malloc(count*sizeof(struct node));
	count = 0;
	for (int i=0; i<inner.n; i++) {
		tree_mul_fill_expanded_args((*nd).args, &count, &inner.args[i]);
	}
}


void tree_expand_mul(
	// IN-OUT
	struct node *nd,
	// INPUT
	char parent_op
) {
	printf("input:\n");
	print_expression_tree(nd, NULL); printf("\n");

	if ((*nd).n == 0) return;

	printf("checkpoint 0:\n");
	print_expression_tree(nd, NULL); printf("\n");

	for (int i=0; i<(*nd).n; i++) {
		tree_expand_mul(&(*nd).args[i], (*nd).op);
	}

	printf("checkpoint 1:\n");
	print_expression_tree(nd, NULL); printf("\n");
	printf("parent op: %c\n", parent_op);

	if ((*nd).op == '*' && parent_op != '*') {
		printf("checkpoint 2:\n");
		print_expression_tree(nd, NULL); printf("\n");

		tree_mul_expand_mul(nd);

		printf("expanded:\n");
		print_expression_tree(nd, NULL); printf("\n");
	}
}


// Function to free memory allocated for the expression tree
// #never-tried
void free_expression_tree(struct node* root) {
	if (root == NULL)
		return;
	if (root->n > 0)
		free(root->args);
	free(root);
}


void tree_mul_count_expanded_args_div(
	// IN-OUT
	int *count,
	// INPUT
	struct node *nd
) {
	if ((*nd).op != '*' && (*nd).op != '/') {
		// print_expression_tree(nd, NULL); printf("\n");
		(*count)++;
		return;
	}
	
	for (int i=0; i<(*nd).n; i++) {
		tree_mul_count_expanded_args_div(count, &(*nd).args[i]);
	}

}


void tree_mul_fill_expanded_args_div(
	// IN-OUT
	struct node *args, int *count, int sgn,
	// INPUT
	struct node *nd
) {
	printf("input of fill:\n");
	print_expression_tree(nd, NULL); printf("\n");

	for (int i=0; i<(*nd).n; i++) {
		if (!((*nd).args[i].op == '/' && (*nd).args[i].args[0].op == '*')) {
			printf("i = %d, tree: ", i); print_expression_tree((*nd).args, NULL); printf("\n");
			if (sgn == 1) {
				args[*count] = (*nd).args[i];
				(*count)++;
			} else {
				create_node(&args[*count], 1, '/');
				args[*count].args[0] = (*nd).args[i];
				(*count)++;
			}
			continue;
		}

		tree_mul_fill_expanded_args_div(args, count, -sgn, &(*nd).args[i].args[0]);
	}

}


void tree_mul_expand_div(
	// IN-OUT
	struct node *nd
) {
	struct node inner;
	inner.n = (*nd).n;
	inner.op = (*nd).op;
	inner.args = (*nd).args;

	int count = 0, sgn = 1;
	tree_mul_count_expanded_args_div(&count, nd);
	printf("count = %d\n", count);
	(*nd).n = count;
	(*nd).args = (struct node*) malloc(count*sizeof(struct node));
	
	count = 0;
	tree_mul_fill_expanded_args_div((*nd).args, &count, sgn, &inner);
}


void tree_expand_div(
	// IN-OUT
	struct node *nd,
	// INPUT
	char parent_op
) {
	printf("input:\n");
	print_expression_tree(nd, NULL); printf("\n");

	if ((*nd).n == 0) return;

	printf("checkpoint 0:\n");
	print_expression_tree(nd, NULL); printf("\n");

	for (int i=0; i<(*nd).n; i++) {
		tree_expand_div(&(*nd).args[i], (*nd).op);
	}

	printf("checkpoint 1:\n");
	print_expression_tree(nd, NULL); printf("\n");
	printf("parent op: %c\n", parent_op);

	if ((*nd).op == '*' && parent_op != '*' && parent_op != '/') {
		printf("checkpoint 2:\n");
		print_expression_tree(nd, NULL); printf("\n");

		tree_mul_expand_div(nd);

		printf("expanded:\n");
		print_expression_tree(nd, NULL); printf("\n");
	}

}


//////
// LINKED LIST IMPLEMENTATION
//////

// Structure to represent an expression tree node
struct lnode {
	int n; // number of operands
	char op; // operation ('+' for addition, '-' for subtraction, '*' for multiplication, '/' for division, '^' for power raising)
	struct lnode *son; // pointer to the first operand node
	struct lnode *bro; // pointer to neighbour brother node
	int number;
	mpz_t mpz;
};


void lnode_build(struct lnode *nd) {
	nd->son = NULL;
	nd->bro = NULL;
	nd->op = '$';
	nd->n = -1;
	nd->number = -1;
}


void lnode_free(struct lnode *nd) {
	if (nd->op == 'n') {
		mpz_clear(nd->mpz);
	}
	if (nd->son) {
		lnode_free(nd->son);
		free(nd->son);
	}
	if (nd->bro) {
		lnode_free(nd->bro);
		free(nd->bro);
	}
}


void lnode_rk2_free(
	struct lnode **nd,
	int dim1, int dim2
) {
	for (int i1=0; i1<dim1; i1++) {
		for (int i2=0; i2<dim2; i2++) {
			lnode_free(&nd[i1][i2]);
		}
	}
}


// Function to print the expression tree
void lnode_print(struct lnode* root, const char **symbols) {
	if (root->op == 's') {
		if (symbols == 	NULL) {
			printf("(%c %d)", root->op, (*root).number);
		} else {
			printf("(%c %s)", root->op, symbols[(*root).number]);
		}
		return;
	} else if(root->op == 'n') {
		// printf("(%c %d)", root->op, (*root).number);
		gmp_printf("(%c %Zd)", root->op, (*root).mpz);
		return;
	} else if(root->op == 'I') {
		printf("(%c)", root->op);
		return;
	} else if (root->op == '^') {
		printf("(%c%d ", root->op, (*root).number);
		lnode_print(root->son, symbols);
		printf(")");
		return;
	} else {
		printf("(%c ", root->op);
		struct lnode *nd = root->son;
		// printf("\nn operand = %d:\n", root->n);
		while (1) {
			// printf("\nbro = %p:\n", nd->bro);
			lnode_print(nd, symbols);
			if ((*nd).bro) {
				printf(" ");
				nd = (*nd).bro;
			} else {
				// printf("\nno more bros\n");
				break;
			}
		}
		printf(")");
	}
}


void lnode_copy(
	// OUTPUT
	struct lnode *out,
	// INPUT
	struct lnode *in
) {
	// printf("enter lnode_copy\n");
	// printf("in: "); lnode_print(in, NULL); printf("\n");
	// printf("copy info\n");
	out->n = in->n;
	out->op = in->op;
	out->number= in->number;
	if (in->op == 'n') {
		// printf("copy mpz\n");
		mpz_init(out->mpz);
		mpz_set(out->mpz, in->mpz);
	}

	// COPY SON
	if (in->son) {
		// printf("copy son\n");	
		if (out->son) {
			lnode_free(out->son);
			free(out->son);
		}
		out->son = (struct lnode*) malloc(sizeof(struct lnode));
		lnode_build(out->son);
		lnode_copy(out->son, in->son);
	} else {
		out->son = NULL;
	}

	// COPY BROTHER
	if (in->bro) {
		// printf("copy bro\n");
		if (out->bro) {
			lnode_free(out->bro);
			free(out->bro);
		}
		out->bro = (struct lnode*) malloc(sizeof(struct lnode));
		lnode_build(out->bro);
		lnode_copy(&out->bro[0], &in->bro[0]);
	} else {
		out->bro = NULL;
	}

}


void lnode_to_str_helper(
	// IN-OUT
	char **output, size_t *output_size,
	// INPUT
	struct lnode *root, char *sep
) {
	if (root == NULL) {
		return;
	}

	// Process the current node
	// #hard-coded: max dim of single leaf-node (to be changed using wp2)
	char *tmp = NULL;
	size_t mpz_buf_size;
	// char tmp[256];
	switch (root->op) {
		case 'I':
			tmp = (char*) malloc(MAX_LEAF_LEN*sizeof(char));
			// snprintf(tmp, sizeof(tmp), "I%s", sep);
			snprintf(tmp, MAX_LEAF_LEN, "I%s", sep);
			break;
		case 'n':
			mpz_buf_size= mpz_sizeinbase(root->mpz, 10) + 2 + MAX_LEAF_LEN;
			tmp = (char*) malloc(mpz_buf_size);
			// snprintf(tmp, sizeof(tmp), "n%s%d%s", sep, root->number, sep);
			gmp_snprintf(tmp, mpz_buf_size, "n%s%Zd%s", sep, root->mpz, sep);
			// gmp_printf("mpz: %Zd\n", root->mpz);
			// printf("tmp: %s\n", tmp);
			// printf("tmp_len: %ld\n", strlen(tmp));
			break;
		case 's':
			tmp = (char*) malloc(MAX_LEAF_LEN*sizeof(char));
			// snprintf(tmp, sizeof(tmp), "s%s%d%s", sep, root->number, sep);
			snprintf(tmp, MAX_LEAF_LEN, "s%s%d%s", sep, root->number, sep);
			break;
		case '-':
		case '/':
			tmp = (char*) malloc(MAX_LEAF_LEN*sizeof(char));
			// snprintf(tmp, sizeof(tmp), "o%s%c%s", sep, root->op, sep);
			snprintf(tmp, MAX_LEAF_LEN, "o%s%c%s", sep, root->op, sep);
			break;
		case '^':
			tmp = (char*) malloc(MAX_LEAF_LEN*sizeof(char));
			if (root->number < 0) {
				// snprintf(tmp, sizeof(tmp), "i%s%d%s", sep, root->number, sep);
				snprintf(tmp, MAX_LEAF_LEN, "o%si%s%d%s", sep, sep, root->number, sep);
			} else {
				// snprintf(tmp, sizeof(tmp), "p%s%d%s", sep, root->number, sep);
				snprintf(tmp, MAX_LEAF_LEN, "o%sp%s%d%s", sep, sep, root->number, sep);
			}
			break;
		default:
		  tmp = (char*) malloc(MAX_LEAF_LEN*sizeof(char));
			// snprintf(tmp, sizeof(tmp), "%c%s%d%s", root->op, sep, root->n, sep);
			// snprintf(tmp, sizeof(tmp), "o%s%c%s%d%s", sep, root->op, sep, root->n, sep);
			snprintf(tmp, MAX_LEAF_LEN, "o%s%c%s%d%s", sep, root->op, sep, root->n, sep);
			break;
	}

	size_t tmp_len = strlen(tmp);
	// printf("tmp: %s\n", tmp);
	// printf("tmp_len: %ld\n", tmp_len);

	// printf("(*output): %s\n", *output);
	size_t output_len = strlen(*output);
	// printf("output_len: %ld\n", output_len);
	// printf("output_size = : %ld\n", *output_size);

	while (output_len + tmp_len >= *output_size) {
		(*output_size) *= 2;
		*output = (char*) realloc(*output, (*output_size)*sizeof(char));
		// printf("new output_size = : %ld\n", *output_size);
		// printf("realloc output: %s\n", output);
		// printf("realloc output_len: %ld\n", strlen(output));
	}

	// Concatenate the current node's data to the output string
	strcat(*output, tmp);
	if (tmp) {free(tmp);}
	// printf("new output: %s\n", output);
	// printf("len: %ld\n", strlen(output));
	

	// Process the children nodes
	struct lnode *child = root->son;
	while (child != NULL) {
		// char child_temp[*output_size];
		lnode_to_str_helper(output, output_size, child, sep);
		child = child->bro;
	}
}


char *lnode_to_str(
	// INPUT
	struct lnode *root, char* sep
) {
	// #hard-coded
	size_t output_size = 1024;
	char *output = (char*) malloc((output_size+1)*sizeof(char));
	output[0] = '\0';
	lnode_to_str_helper(&output, &output_size, root, sep);
	return output;
}


// Function to create a new node
void create_lnode(
	// OUTPUT
	struct lnode *new_node,
	// INPUT
	int n, char op
) {
	// printf("creating lnode with n = %d, op = %c\n", n, op);
	new_node->n = n;
	new_node->op = op;
	if (op == 's' || op == 'n' || op == 'I') {
		new_node->son = NULL;
	} else {
		new_node->son = (struct lnode*) malloc(sizeof(struct lnode));
		lnode_build(new_node->son);
	}
	return;
}


void free_tree(struct lnode *root) {
	if (root == NULL) {
		return;
	}
	
	// Free son nodes
	free_tree(root->son);
	
	// Free brother nodes
	free_tree(root->bro);

	// Free mpz
	if (root->op == 'n') {
		mpz_clear(root->mpz);
	}

	// Free the current node
	free(root);
}

// <<lnode_remove_son>>
void lnode_remove_son(
	// IN-OUT
	struct lnode *parent,
	// INPUT
	struct lnode *son_to_remove
) {
	if (parent == NULL || son_to_remove == NULL) {
		return;
	}

	// If the son to remove is the first son
	if (parent->son == son_to_remove) {
		parent->son = son_to_remove->bro;
		free_tree(son_to_remove->son);  // Free memory of the subtree rooted at the son
		free(son_to_remove);
		parent->n -= 1;
		return;
	}

	// If the son to remove is not the first son, find it in the list of brothers
	struct lnode *current = parent->son;
	while (current != NULL && current->bro != son_to_remove) {
		current = current->bro;
	}

	// If the son to remove was found in the list of brothers
	if (current != NULL) {
		current->bro = son_to_remove->bro;
		free_tree(son_to_remove->son);  // Free memory of the subtree rooted at the son
		free(son_to_remove);
		parent->n -= 1;
	}
}


void lnode_attach_son_first(
	// IN-OUT
	struct lnode *nd
	// // INPUT
	// int n, char op
) {
	// printf("son address = %p\n", nd->son);
	struct lnode *holder = (*nd).son;
	nd->son = (struct lnode*) malloc(sizeof(struct lnode));
	lnode_build(nd->son);
	nd->son->bro = holder;
	// printf("son address = %p\n", nd->son);
	// printf("son bro address = %p\n", nd->son->bro);
	// create_lnode(&(*nd).son, n, op);
	// (*nd).son.bro = holder;
}


// Function to parse the expression recursively
void lnode_parse_expression(
	// OUTPUT
	struct lnode *root,
	// INPUT
	char* expression,
	const char** symbols, int num_symbols,
	int bro
) {
	// printf("\ninput expression: %s\n", expression);

	// Remove white spaces from expression
	char* p = expression;
	char* q = expression;
	while (*q != '\0') {
		if (*q != ' ') {
			*p = *q;
			p++;
		}
		q++;
	}
	if (*(p-1) == '\n') {
		*(p-1) = '\0';
	} else {
		*p = '\0';
	}
	p = NULL;
	q = NULL;

	// Remove outer parentheses if present
	int n = strlen(expression);
	int parenthesis_count = 0;

	// Check if the expression starts with '('
	while (expression[0] == '(') {
		// Scan the string starting from the character right after the '('
		int i;
		for (i = 1; i < n; i++) {
			if (expression[i] == '(') {
				parenthesis_count++;
			} else if (expression[i] == ')') {
				parenthesis_count--;
			}

			if (parenthesis_count == -1) {
				break;
			}
		}

		// printf("parenthesis_count = %d\n", parenthesis_count);
		if (i == n-1) {
			if (expression[n-1] != ')') {
				fprintf(stderr, "Error: Unmatched parenthesis.\n");
				exit(1);
			} else {
				expression[n-1] = '\0';
				expression++;
				// printf("unnecessary outer parenthesis removed\n");
				continue;
			}
		} else {
			break;
		}
	}

	// Base case: expression is a symbol
	for (int i = 0; i < num_symbols; i++) {
		if (strcmp(expression, symbols[i]) == 0) {
			// printf("expression is a symbol\n");
			create_lnode(root, 0, 's');
			(*root).number = i;
			if (bro == 0) root->bro = NULL;
			expression = NULL;
			return;
		}
	}

	// check whether expression is the imaginary unit
	if (strcmp(expression, "I") == 0) {
		// printf("expression is the imaginary unit\n");
		create_lnode(root, 0, 'I');
		if (bro == 0) root->bro = NULL;
		expression = NULL;
		return;
	}

	// int lowest_precedence = 1000; // Initialize to a large value
	int lowest_precedence = 1000000; // Initialize to a large value
	// int lowest_precedence_index[100]; // Maximum number of operators with the same precedence
	// int lowest_precedence_index[10000]; // Maximum number of operators with the same precedence
	int lowest_precedence_index[100000]; // Maximum number of operators with the same precedence
	int num_lowest_precedence;
	find_precedence(
		&num_lowest_precedence, lowest_precedence_index, &lowest_precedence,
		expression
	);

	// If no operator is found, then the operand is a number
	// if (lowest_precedence == 1000) {
	if (lowest_precedence == 1000000) {
		// printf("expression is a number\n");
		create_lnode(root, 0, 'n');
		mpz_init(root->mpz);
		mpz_set_str(root->mpz, expression, 10);
		// if (mpz_set_str(root->mpz, expression, 10)) {
		// 	perror("error while setting mpz from string during parsing");
		// 	printf("expression:\n%s\n", expression);
		// 	exit(1);
		// };
		// gmp_printf("created number with mpz = %Zd:\n", (*root).mpz);
		// (*root).number = atoi(expression);
		// printf("created numeric root with number = %d:\n", (*root).number);
		if (bro == 0) root->bro = NULL;
		expression = NULL;
		return;
	}

	// Split the expression into parts based on the lowest precedence operators
	// printf("num_lowest_precedence = %d\n", num_lowest_precedence);
	if (num_lowest_precedence > 1) {
		// Parse each sub-string between lowest precedence operators
		// starting from the right-most one, so that we can read it's operation
		// character, use that info, and then replace that char with '\0'
		// when processing next sub-string

		// // struct lnode *root_args;
		// if (lowest_precedence_index[0] != 0) {
		// 	// if the first lowest precedence operation is not at the beginning
		// 	// of the expression, the number of operation is the number of
		// 	// sub-strings plus one
		// 	create_node(root, num_lowest_precedence+1, '\0');
		// 	root_args = (*root).args + 1;
		// } else {
		// 	create_node(root, num_lowest_precedence, '\0');
		// 	root_args = (*root).args;
		// }

		if (lowest_precedence == 2) {
			create_lnode(root, num_lowest_precedence, '+');
			root->son->bro = NULL;
			(*root).op = '+';

			// right-most sub-string
			char op = expression[lowest_precedence_index[num_lowest_precedence-1]];
			if (op == '+') {
				lnode_parse_expression(
					(*root).son,
					expression + lowest_precedence_index[num_lowest_precedence-1] + 1,
					symbols, num_symbols,
					0
				);
			} else if (op == '-') {
				create_lnode((*root).son, 1, '-');
				root->son->son->bro = NULL;
				lnode_parse_expression(
					root->son->son,
					expression + lowest_precedence_index[num_lowest_precedence-1] + 1,
					symbols, num_symbols,
					0
				);			
			}

			// all the other sub-strings
			for (int i=num_lowest_precedence-2; i>=0; i--) {
				expression[lowest_precedence_index[i+1]] = '\0';
				op = expression[lowest_precedence_index[i]];
				lnode_attach_son_first(root);
				if (op == '+') {
					// Apply parse_expression to the segment and assign it to the corresponding sub-string
					lnode_parse_expression(
						(*root).son,
						expression + lowest_precedence_index[i] + 1,
						symbols, num_symbols,
						1
					);
				} else if (op == '-') {
					// Create a node with op='-', n=1, and arg[0] as the result of parse_expression for the sub-string
					create_lnode((*root).son, 1, '-');
					root->son->son->bro = NULL;
					lnode_parse_expression(
						root->son->son,
						expression + lowest_precedence_index[i] + 1,
						symbols, num_symbols,
						1
					);
				}
			}
		} else if (lowest_precedence == 3) {
			create_lnode(root, num_lowest_precedence, '*');
			root->son->bro = NULL;
			(*root).op = '*';

			// right-most sub-string
			char op = expression[lowest_precedence_index[num_lowest_precedence-1]];
			if (op == '*') {
				// printf("top ex: %s\n", expression + lowest_precedence_index[num_lowest_precedence-1] + 1);
				lnode_parse_expression(
					(*root).son,
					expression + lowest_precedence_index[num_lowest_precedence-1] + 1,
					symbols, num_symbols,
					0
				);
			} else if (op == '/') {
				create_lnode((*root).son, 1, '/');
				root->son->son->bro = NULL;
				lnode_parse_expression(
					root->son->son,
					expression + lowest_precedence_index[num_lowest_precedence-1] + 1,
					symbols, num_symbols,
					0
				);			
			}

			// all the other sub-strings
			for (int i=num_lowest_precedence-2; i>=0; i--) {
				expression[lowest_precedence_index[i+1]] = '\0';
				// printf("ex: %s\n", expression);
				op = expression[lowest_precedence_index[i]];
				lnode_attach_son_first(root);
				if (op == '*') {
					// Apply parse_expression to the segment and assign it to the corresponding sub-string
					lnode_parse_expression(
						(*root).son,
						expression + lowest_precedence_index[i] + 1,
						symbols, num_symbols,
						1
					);
				} else if (op == '/') {
					// Create a node with op='-', n=1, and arg[0] as the result of parse_expression for the sub-string
					create_lnode((*root).son, 1, '/');
					root->son->son->bro = NULL;
					lnode_parse_expression(
						root->son->son,
						expression + lowest_precedence_index[i] + 1,
						symbols, num_symbols,
						1
					);
				}				
			}
		}

		if (lowest_precedence_index[0] != 0) {
			// add the root corresponding to the first sub-string
			// printf("add the root corresponding to the first sub-string\n");
			// printf("ex: %s\n", expression);
			// printf("low prec idx = %d\n", lowest_precedence_index[0]);
			root->n++;
			expression[lowest_precedence_index[0]] = '\0';
			lnode_attach_son_first(root);
			lnode_parse_expression(
				(*root).son,
				expression,
				symbols, num_symbols,
				1
			);
		}
		expression = NULL;
		return;
	} else {
		int idx = lowest_precedence_index[0];
		char op = expression[idx];
		// printf("op = %c, idx = %d\n", op, idx);
		if (idx == 0) {
			if (op == '+' || op == '-') {
				create_lnode(root, 1, op);
				root->son->bro = NULL;
				lnode_parse_expression(
					(*root).son,
					expression + 1,
					symbols, num_symbols,
					0
				);
			} else {
				fprintf(stderr, "Error: expression starting with an operation different from '+' or '-'.\n");
				exit(1);
			}
		} else {
			if (op == '^') {
				create_lnode(root, 1, '^');
				root->son->bro = NULL;
				if (expression[idx+1] == '(') {
					expression[strlen(expression)] = '\0';
					root->number = atoi(expression+idx+2);
				} else {
					root->number = atoi(expression+idx+1);
				}
				// root->number = atoi(expression+idx+1);
				expression[idx] = '\0';
				lnode_parse_expression(
					root->son,
					expression,
					symbols, num_symbols,
					0
				);
			} else if (op == '-') {
				create_lnode(root, 2, '+');
				root->son->bro = NULL;
				create_lnode(root->son, 1, op);
				root->son->son->bro = NULL;
				lnode_parse_expression(
					root->son->son,
					expression + idx + 1,
					symbols, num_symbols,
					0
				);
				expression[idx] = '\0';
				lnode_attach_son_first(root);
				lnode_parse_expression(
					root->son,
					expression,
					symbols, num_symbols,
					1
				);
			} else if (op == '/') {
				create_lnode(root, 2, '*');
				root->son->bro = NULL;
				create_lnode(root->son, 1, op);
				root->son->son->bro = NULL;
				lnode_parse_expression(
					root->son->son,
					expression + idx + 1,
					symbols, num_symbols,
					0
				);
				expression[idx] = '\0';
				lnode_attach_son_first(root);
				// printf("residual ex: %s\n", expression);
				lnode_parse_expression(
					root->son,
					expression,
					symbols, num_symbols,
					1
				);
			} else {
				create_lnode(root, 2, op);
				root->son->bro = NULL;
				lnode_parse_expression(
					root->son,
					expression + idx + 1,
					symbols, num_symbols,
					0
				);
				expression[idx] = '\0';
				lnode_attach_son_first(root);
				lnode_parse_expression(
					root->son,
					expression,
					symbols, num_symbols,
					1
				);
			}
		}
		expression = NULL;
		return;
	}

}


void lnode_expand(
	// IN-OUT
	struct lnode *nd
) {
	// printf("lnode_expand input:\n"); lnode_print(nd, NULL); printf("\n");
	if (nd->son == NULL) {
		return;
	}

	struct lnode *son;

	// expand each child
	son = nd->son;
	while (son) {
		lnode_expand(son);
		son = son->bro;
	}

	if (nd->op == '*') {
		// deal with first son
		son = nd->son;
		if (son->op == '*') {
			nd->n += son->n - 1;

			struct lnode *grandson = son->son;

			// attach last grandson to next bro
			while (grandson->bro) {
				grandson = grandson->bro;
			}
			grandson->bro = son->bro;

			// attach first grandson to parent as first son
			nd->son = son->son;

			if (grandson->bro) {
				son = grandson;
			} else {
				return;
			}
		} else if (son->op == '/' && son->son->op == '*') {
			nd->n += son->son->n - 1;
			
			struct lnode *grandson = son->son;
			
			// change grandgrandson into divisions and attach last one
			while (1) {
				grandson->op = '/';
				grandson->n = 1;
				if (grandson->son->bro) {
					grandson->bro = (struct lnode*) malloc(sizeof(struct lnode));
					lnode_build(grandson->bro);
					grandson->bro->son = grandson->son->bro;
					grandson->son->bro = NULL;
					grandson = grandson->bro;
				} else {
					break;
				}
			}	
			grandson->bro = son->bro;

			// attach first grandgrandson to parent as first bro
			nd->son = son->son;

			// if (grandson->bro) {
			// 	son = grandson;
			// } else {
			// 	return;
			// }
		} else if (son->bro->op == '/' && son->bro->son->op == '/') {
			nd->n += son->son->n - 1;
			
			struct lnode *grandson = son->son->son;

			// attach grandson to next bro
			grandson->bro = son->bro;

			// attach first grandson to current bro
			son->bro = son->bro->son->son;

			// if (grandson->bro) {
			// 	son = grandson;
			// } else {
			// 	return;
			// }
		} else if (son->op == '^' && son->son->op == '*') {
			nd->n += son->son->n - 1;
			
			struct lnode *grandson = son->son;
			
			// change grandgrandson into powers and attach last one
			while (1) {
				grandson->op = '^';
				grandson->number = son->number;
				grandson->n = 1;
				if (grandson->son->bro) {
					grandson->bro = (struct lnode*) malloc(sizeof(struct lnode));
					lnode_build(grandson->bro);
					grandson->bro->son = grandson->son->bro;
					grandson->son->bro = NULL;
					grandson = grandson->bro;
				} else {
					break;
				}
			}
			grandson->bro = son->bro;

			// attach first grandgrandson to parent as first bro
			nd->son = son->son;

			// if (grandson->bro) {
			// 	son = grandson;
			// } else {
			// 	return;
			// }
		}

		// deal with other siblings
		while (son->bro) {
			if (son->bro->op == '*') {
				// lnode_print(son, NULL); printf("\n");
				nd->n += son->bro->n - 1;
				
				struct lnode *grandson = son->bro->son;				

				// attach last grandson to next bro
				while (grandson->bro) {
					grandson = grandson->bro;
				}
				grandson->bro = son->bro->bro;

				// attach first grandson to current bro
				son->bro = son->bro->son;

				if (grandson->bro) {
					son = grandson;
					continue;
				} else {
					break;
				}
			} else if (son->bro->op == '/' && son->bro->son->op == '*') {
				// printf("here1\n");
				// printf("son:\n"); lnode_print(son, NULL); printf("\n");
				// printf("son->bro:\n"); lnode_print(son->bro, NULL); printf("\n");
				// printf("son->bro->bro:\n"); lnode_print(son->bro->bro, NULL); printf("\n");
				nd->n += son->bro->son->n - 1;

				struct lnode *grandson = son->bro->son;

				// change grandson into divisions and attach last one
				// printf("start while loop\n");
				while (1) {
					// printf("son:\n"); lnode_print(son, NULL); printf("\n");
					// if (son->bro) {printf("son->bro:\n"); lnode_print(son->bro, NULL); printf("\n");}
					// if (son->bro->bro) {printf("son->bro->bro:\n"); lnode_print(son->bro->bro, NULL); printf("\n");}
					// printf("grandson:\n"); lnode_print(grandson, NULL); printf("\n");
					grandson->op = '/';
					grandson->n = 1;
					// printf("grandson:\n"); lnode_print(grandson, NULL); printf("\n");
					if (grandson->son->bro) {
						grandson->bro = (struct lnode*) malloc(sizeof(struct lnode));
						lnode_build(grandson->bro);
						grandson->bro->son = grandson->son->bro;
						grandson->son->bro = NULL;
						// printf("grandson:\n"); lnode_print(grandson, NULL); printf("\n");
						grandson = grandson->bro;
					} else {
						break;
					}
				}
				// printf("end while loop\n");
				// printf("son:\n"); lnode_print(son, NULL); printf("\n");
				// if (son->bro) {printf("son->bro:\n"); lnode_print(son->bro, NULL); printf("\n");}
				// if (son->bro->bro) {printf("son->bro->bro:\n"); lnode_print(son->bro->bro, NULL); printf("\n");}
				// printf("grandson:\n"); lnode_print(grandson, NULL); printf("\n");
				grandson->bro = son->bro->bro;

				// attach first grandson to current bro
				son->bro = son->bro->son;
				// printf("nd:\n"); lnode_print(nd, NULL); printf("\n");
				
				// if (grandson->bro) {
				// 	son = grandson;
				// 	printf("son:\n"); lnode_print(son, NULL); printf("\n");
				// 	continue;
				// } else {
				// 	break;
				// }
			} else if (son->bro->op == '/' && son->bro->son->op == '/') {
				// printf("enter / / with:\n"); lnode_print(son, NULL); printf("\n");
				nd->n += son->bro->n - 1;
				
				struct lnode *grandson = son->bro->son->son;

				// attach grandson to next bro
				grandson->bro = son->bro->bro;

				// attach grandson to current bro
				son->bro = son->bro->son->son;

				// if (grandson->bro) {
				// 	son = grandson;
				// 	continue;
				// } else {
				// 	break;
				// }
			} else if (son->op == '^' && son->son->op == '*') {
				nd->n += son->bro->son->n - 1;

				struct lnode *grandson = son->bro->son;

				// change grandson into powers and attach last one
				while (1) {
					grandson->op = '^';
					grandson->number = son->number;
					grandson->n = 1;
					if (grandson->son->bro) {
						grandson->bro = (struct lnode*) malloc(sizeof(struct lnode));
						lnode_build(grandson->bro);
						grandson->bro->son = grandson->son->bro;
						grandson->son->bro = NULL;
						grandson = grandson->bro;
					} else {
						break;
					}
				}
				grandson->bro = son->bro->bro;

				// attach first grandson to current bro
				son->bro = son->bro->son;

			} else {
				// Move to the next sibling
				son = son->bro;
			}
		}

	} else if (nd->op == '^') {
		if (nd->son->op == '^') {
			nd->number *= nd->son->number;
			nd->son = nd->son->son;
			// lnode_expand(nd);
		} else if (nd->son->op == '/') {
			nd->number = -nd->number;
			nd->son = nd->son->son;
			// lnode_expand(nd);
		} else if (nd->son->op == '*') {
			// printf("changing ^* into *^\n");
			// printf("before:\n"); lnode_print(nd, NULL); printf("\n");
			struct lnode *grandson = nd->son;
			
			// change grandgrandson into brother powers
			while (1) {
				grandson->op = '^';
				grandson->number = nd->number;
				grandson->n = 1;
				if (grandson->son->bro) {
					grandson->bro = (struct lnode*) malloc(sizeof(struct lnode));
					lnode_build(grandson->bro);
					grandson->bro->son = grandson->son->bro;
					grandson->son->bro = NULL;
					grandson = grandson->bro;
				} else {
					break;
				}
			}

			// change original power node into product node
			nd->number = -1;
			nd->op = '*';
			// printf("after 1:\n"); lnode_print(nd, NULL); printf("\n");

			// // expand new product
			// lnode_expand(nd);
			// printf("after 2:\n"); lnode_print(nd, NULL); printf("\n");

		}
	} else if (nd->op == '/') {
		if (nd->son->op == '^') {
			nd->op = '^';
			nd->number = -nd->son->number;
			nd->son = nd->son->son;
			// lnode_expand(nd);
		}
	}

	// printf("lnode_expand checkpoint:\n"); lnode_print(nd, NULL); printf("\n");
	
	// expand each child
	son = nd->son;
	while (son) {
		lnode_expand(son);
		son = son->bro;
	}
}


void lnode_div_to_pow(
	// IN-OUT
	struct lnode *root
) {
	if (root == NULL) {
		return;
	}

	// Recursively process the children nodes
	struct lnode *child = root->son;
	while (child != NULL) {
		lnode_div_to_pow(child);
		child = child->bro;
	}

	// Process the current node
	if (root->op == '/') {
		root->op = '^';
		root->number = -1;
	}
}


void lnode_contains_sym_inner(
	// OUTPUT
	int *sym_found, int nsym,
	// INPUT
	struct lnode *root
) {
	if (root->op == 's') {
		sym_found[root->number] = 1;
		return;
	} else if(root->op == 'n') {
		return;
	} else {
		struct lnode *nd = root->son;
		while (1) {
			lnode_contains_sym_inner(sym_found, nsym, nd);
			if ((*nd).bro) {
				nd = (*nd).bro;
			} else {
				break;
			}
		}
	}
}


void lnode_contains_sym(
	// OUTPUT
	int *sym_found, int nsym,
	// INPUT
	struct lnode *root
) {
	// initialize to zero
	for (int s=0; s<nsym; s++) {
		sym_found[s] = 0;
	}

	lnode_contains_sym_inner(sym_found, nsym, root);
	return;
}


void lnode_coefficient_list_helper(
	// OUTPUT
	struct lnode **coeffs, int **pows, int *nterms,
	// INPUT
	struct lnode *pol, int sym_idx
) {
	// printf("input: "); lnode_print(pol, NULL); printf("\n");
	struct lnode *tmp, *tmp1, *coeff_nd;
	int dim = 0;
	struct lnode *nd;
	int nsons, sign;

	// remove '-' is input is a negative
	if (pol->op == '-') {
		tmp = pol->son;
		sign = -1;
	} else {
		tmp = pol;
		sign = 1;
	}

	// select children if what remains is a sum
	if (tmp->op == '+') {
		nd = tmp->son;
		nsons = tmp->n;
	} else {
		nd = tmp;
		nsons = 1;
	}
	// printf("nsons = %d\n", nsons);

	*coeffs = (struct lnode*) malloc(nsons*sizeof(struct lnode));
	for (int s=0; s<nsons; s++) {
		lnode_build(&(*coeffs)[s]);
	}
	*pows = (int*) malloc(nsons*sizeof(int));

	int coeff_sign;
	for (int s=0; s<nsons; s++) {
		// printf("term s = %d: ", s); lnode_print(nd, NULL); printf("\n");
		if (nd->op == '-') {
			coeff_sign = -1;
			tmp1 = nd->son;
		} else {
			coeff_sign = 1;
			tmp1 = nd;
		}
		if (sign*coeff_sign == -1) {
			create_lnode(&(*coeffs)[dim], 1, '-');
			coeff_nd = (*coeffs)[dim].son;
			lnode_build(coeff_nd);
		} else {
			coeff_nd = &(*coeffs)[dim];
		}
		if (tmp1->op == '*') {
			lnode_copy(coeff_nd, tmp1);
			// lnode_print(&(*coeffs)[dim], NULL); printf("\n");
			// term is tree*eta or tree*eta^pow
			// search eta or eta^pow through its factors
			tmp = coeff_nd->son;
			while (tmp) {
				// printf("factor: "); lnode_print(tmp, NULL); printf("\n");
				if (tmp->op == '^') {
					if (tmp->son->op == 's') {
						if (tmp->son->number == sym_idx) {
							// found eta^pow
							(*pows)[dim] = tmp->number;
							break;
						}
					}
				} else if (tmp->op == 's') {
					if (tmp->number == sym_idx) {
						// found eta
						(*pows)[dim] = 1;
						break;
					}
				}
				tmp = tmp->bro;
			}
			if (tmp) {
				// printf("remove node: "); lnode_print(tmp, NULL); printf("\n");
				lnode_remove_son(coeff_nd, tmp);
				// printf("coeffs[%d]: ", dim); lnode_print(&(*coeffs)[dim], NULL); printf("\n");
			} else {
				(*pows)[dim] = 0;
			}
		} else {
			// term is tree, eta or eta^pow
			int etaless_tree = 0;
			// check for eta or eta^pow
			if (tmp1->op == '^') {
				if (tmp1->son->op == 's') {
					if (tmp1->son->number == sym_idx) {
						// found eta^pow
						(*pows)[dim] = tmp1->number;
						create_lnode(coeff_nd, 0, 'n');
						mpz_init(coeff_nd->mpz);
						mpz_set_ui(coeff_nd->mpz, 1);
					} else {
						// found sym^pow
						etaless_tree = 1;
					}
				} else {
					// fount tree^pow
					etaless_tree = 1;
				}
			} else if (tmp1->op == 's') {
				if (tmp1->number == sym_idx) {
					// found eta
					(*pows)[dim] = 1;
					create_lnode(coeff_nd, 0, 'n');
					mpz_init(coeff_nd->mpz);
					mpz_set_ui(coeff_nd->mpz, 1);
				} else {
					// found symbol different from eta
					etaless_tree = 1;
				}
			} else {
				// everything else
				etaless_tree = 1;
			}

			if (etaless_tree) {
				// found eta-independent tree
				(*pows)[dim] = 0;
				lnode_copy(coeff_nd, tmp1);
			}

		}

		dim++;
		nd = nd->bro;
	}

	*nterms = dim;
	return;
}


void lnode_neg(
	// IN-OUT
	struct lnode *nd
) {

	if (nd->op == '-') {
		printf("neg of a neg\n");
		perror("neg of a neg");
		//////
		// REMOVE NEG
		//////
		// #to-be-reviewed
		// transfer son to parent
		nd->op = nd->son->op;
		nd->n = nd->son->n;
		nd->number = nd->son->number;
		if (nd->op == 'n') {
			mpz_init(nd->son->mpz);
			mpz_set(nd->son->mpz, nd->mpz);
			mpz_clear(nd->mpz);
		}
		nd->bro = nd->son->bro;
		struct lnode *holder = nd->son;
		nd->son = nd->son->son;
		free(holder);
		
	} else {
		//////
		// ADD NEG
		//////
		// hold son
		struct lnode *holder = nd->son;

		// transfer parent to son
		nd->son = (struct lnode*) malloc(sizeof(struct lnode));
		nd->son->op = nd->op;
		nd->son->n = nd->n;
		nd->son->number = nd->number;
		if (nd->op == 'n') {
			mpz_init(nd->son->mpz);
			mpz_set(nd->son->mpz, nd->mpz);
			mpz_clear(nd->mpz);
		}
		nd->son->bro = nd->bro;
		nd->bro = NULL;
		nd->son->son = holder;
			
		// change lnode to neg
		nd->op = '-';
		nd->n = 1;
	}

	return;
}


//////
// IN/OUT
//////
void lnode_to_file(FILE *file, struct lnode *nd) {
	fwrite(&nd->n, sizeof(int), 1, file);
	fwrite(&nd->op, sizeof(char), 1, file);
	fwrite(&nd->number, sizeof(int), 1, file);

	if (nd->op == 'n') {
		// Calcola correttamente la dimensione in byte
		size_t mpz_size = (mpz_sizeinbase(nd->mpz, 2) + 7) / 8; // Arrotonda verso l'alto
		fwrite(&mpz_size, sizeof(size_t), 1, file);

		void *mpz_data = calloc(1, mpz_size); // Usa calloc per inizializzare a zero
		mpz_export(mpz_data, NULL, 1, 1, 0, 0, nd->mpz);
		fwrite(mpz_data, 1, mpz_size, file);
		free(mpz_data);
	}

	// Scrivi ricorsivamente son e bro
	int has_son = (nd->son != NULL);
	fwrite(&has_son, sizeof(int), 1, file);
	if (has_son) {
		lnode_to_file(file, nd->son);
	}

	int has_bro = (nd->bro != NULL);
	fwrite(&has_bro, sizeof(int), 1, file);
	if (has_bro) {
		lnode_to_file(file, nd->bro);
	}
}


void lnode_from_file(FILE *file, struct lnode *nd) {
	fread(&nd->n, sizeof(int), 1, file);
	fread(&nd->op, sizeof(char), 1, file);
	fread(&nd->number, sizeof(int), 1, file);

	// Leggi mpz_t
	if (nd->op == 'n') {
		size_t mpz_size;
		fread(&mpz_size, sizeof(size_t), 1, file);
		void *mpz_data = malloc(mpz_size);
		fread(mpz_data, 1, mpz_size, file);
		mpz_init(nd->mpz);
		mpz_import(nd->mpz, mpz_size, 1, 1, 0, 0, mpz_data);
		free(mpz_data);
	}

	// Leggi ricorsivamente son e bro
	int has_son;
	fread(&has_son, sizeof(int), 1, file);
	if (has_son) {
		nd->son = (struct lnode *)malloc(sizeof(struct lnode));
		lnode_from_file(file, nd->son);
	} else {
		nd->son = NULL;
	}

	int has_bro;
	fread(&has_bro, sizeof(int), 1, file);
	if (has_bro) {
		nd->bro = (struct lnode *)malloc(sizeof(struct lnode));
		lnode_from_file(file, nd->bro);
	} else {
		nd->bro = NULL;
	}
}


void lnode_rk2_to_file(
	FILE *file, struct lnode **nd,
	int dim1, int dim2
) {
	for (int i1=0; i1<dim1; i1++) {
		for (int i2=0; i2<dim2; i2++) {
			lnode_to_file(file, &nd[i1][i2]);
		}
	}
}


void lnode_rk2_from_file(
	FILE *file, struct lnode **nd,
	int dim1, int dim2
) {
	for (int i1=0; i1<dim1; i1++) {
		for (int i2=0; i2<dim2; i2++) {
			lnode_from_file(file, &nd[i1][i2]);
		}
	}
}

