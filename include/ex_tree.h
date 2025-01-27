#ifndef EX_TREE_H
#define EX_TREE_H

struct node {
	int n;
	char op;
	struct node *args;
	int number;
};

void print_expression_tree(struct node* root, const char **symbols);

void create_node(struct node *new_node, int n, char op);

void parse_expression(
	struct node *root,
	char *expression,
	const char** symbols, int num_symbols
);

void node_expand_mul(struct node *nd);

void node_expand(struct node *nd);

void free_expression_tree(struct node* root);

void tree_mul_count_expanded_args(
	// IN-OUT
	int *count,
	// INPUT
	struct node *nd
);

void tree_mul_fill_expanded_args(
	// IN-OUT
	struct node *args, int *count, 
	// INPUT
	struct node *nd
);

void tree_mul_expand_mul(
	// IN-OUT
	struct node *nd
);

void tree_expand_mul(
	// IN-OUT
	struct node *nd,
	// INPUT
	char parent_op
);

void tree_mul_count_expanded_args_div(
	// IN-OUT
	int *count,
	// INPUT
	struct node *nd
);

void tree_mul_fill_expanded_args_div(
	// IN-OUT
	struct node *args, int *count, int sgn,
	// INPUT
	struct node *nd
);

void tree_mul_expand_div(
	// IN-OUT
	struct node *nd
);

void tree_expand_div(
	// IN-OUT
	struct node *nd,
	// INPUT
	char parent_op
);

//////
// LINKED LIST IMPLEMENTATION
//////
struct lnode {
	int n;
	char op;
	struct lnode *son;
	struct lnode *bro;
	int number;
	mpz_t mpz;
};

void lnode_build(struct lnode *nd);

void lnode_free(struct lnode *nd);

void lnode_rk2_free(
	struct lnode **nd,
	int dim1, int dim2
);

void lnode_copy(
	// OUTPUT
	struct lnode *out,
	// INPUT
	struct lnode *in
);

void lnode_print(struct lnode* root, const char **symbols);

char *lnode_to_str(
	// INPUT
	struct lnode *root, char* sep
);

void create_lnode(
	// OUTPUT
	struct lnode *new_node,
	// INPUT
	int n, char op
);

void free_tree(struct lnode *root);

void lnode_remove_son(
	// IN-OUT
	struct lnode *parent,
	// INPUT
	struct lnode *son_to_remove
);

void lnode_attach_son_first(
	// IN-OUT
	struct lnode *nd
	// // INPUT
	// int n, char op
);

void lnode_parse_expression(
	// OUTPUT
	struct lnode *root,
	// INPUT
	char* expression,
	const char** symbols, int num_symbols,
	int bro
);

void lnode_expand(
	// IN-OUT
	struct lnode *nd
);

void lnode_div_to_pow(
	// IN-OUT
	struct lnode *root
);

void lnode_contains_sym(
	// OUTPUT
	int *sym_found, int nsym,
	// INPUT
	struct lnode *root
);

void lnode_coefficient_list_helper(
	// OUTPUT
	struct lnode **coeffs, int **pows, int *nterms,
	// INPUT
	struct lnode *pol, int sym_idx
);

void lnode_neg(
	// IN-OUT
	struct lnode *nd
);

void lnode_to_file(
	FILE *file, struct lnode *nd
);

void lnode_from_file(
	FILE *file, struct lnode *nd
);

void lnode_rk2_to_file(
	FILE *file, struct lnode **nd,
	int dim1, int dim2
);

void lnode_rk2_from_file(
	FILE *file, struct lnode **nd,
	int dim1, int dim2
);

#endif