#include "LAGraph.h"
#include "symbol_list.h"
#include "terminal/convert_graph.h"
#include "terminal/terminal_format.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static char msg[LAGRAPH_MSG_LEN];

// Utils

static GrB_Matrix make_empty_matrix(GrB_Index n) {
	GrB_Matrix m;
	GrB_Matrix_new(&m, GrB_BOOL, n, n);
	return m;
}

static GrB_Index matrix_nvals(GrB_Matrix m) {
	GrB_Index nvals = 0;
	GrB_Matrix_nvals(&nvals, m);
	return nvals;
}

static bool matrix_has_edge(GrB_Matrix m, GrB_Index i, GrB_Index j) {
	bool val = false;
	GrB_Info info = GrB_Matrix_extractElement_BOOL(&val, m, i, j);
	return info == GrB_SUCCESS && val;
}

typedef struct {
	GrB_Matrix *matrices;
	SymbolList symbol_list;
	size_t count;
} MatrixSet;

static MatrixSet make_matrix_set(GrB_Index n, char **labels, size_t count) {
	MatrixSet ms;
	ms.count = count;
	ms.matrices = malloc(count * sizeof(GrB_Matrix));
	ms.symbol_list = symbol_list_create();

	for (size_t i = 0; i < count; i++) {
		ms.matrices[i] = make_empty_matrix(n);
		symbol_list_add_str(&ms.symbol_list, labels[i], false);
	}

	return ms;
}

static void free_matrix_set(MatrixSet *ms) {
	for (size_t i = 0; i < ms->count; i++) {
		GrB_Matrix_free(&ms->matrices[i]);
	}
	free(ms->matrices);
	symbol_list_free(&ms->symbol_list);
}

static void free_mr_graph_arrays(MRGraph *g) {
	free(g->open_par);
	free(g->close_par);
	free(g->open_bra);
	free(g->close_bra);
}

// Tests

static void test_null_arguments(void) {
	MRGraph out = {0};
	GrB_Matrix matrices[1];
	SymbolList sl = symbol_list_create();

	assert(build_mr_graph(NULL, &sl, 3, &DefaultTerminalFormat, &out) ==
		   GrB_INVALID_VALUE);
	assert(build_mr_graph(matrices, NULL, 3, &DefaultTerminalFormat, &out) ==
		   GrB_INVALID_VALUE);
	assert(build_mr_graph(matrices, &sl, 3, NULL, &out) == GrB_INVALID_VALUE);
	assert(build_mr_graph(matrices, &sl, 3, &DefaultTerminalFormat, NULL) ==
		   GrB_INVALID_VALUE);

	symbol_list_free(&sl);
}

static void test_empty_symbol_list(void) {
	SymbolList sl = symbol_list_create();
	MRGraph out = {0};

	GrB_Info info = build_mr_graph(NULL, &sl, 3, &DefaultTerminalFormat, &out);
	assert(info == GrB_INVALID_VALUE); // matrices is NULL

	GrB_Matrix dummy[1] = {NULL};
	info = build_mr_graph(dummy, &sl, 3, &DefaultTerminalFormat, &out);
	assert(info == GrB_SUCCESS);
	assert(out.n_par == 0);
	assert(out.n_bra == 0);
	assert(out.normal == NULL);
	assert(out.open_par == NULL);
	assert(out.close_par == NULL);

	symbol_list_free(&sl);
}

// One open_par and one close_par are correctly assigned
static void test_single_parenthesis_pair(void) {
	GrB_Index n = 3;
	char *labels[] = {"op_0", "cp_0"};
	MatrixSet ms = make_matrix_set(n, labels, 2);

	GrB_Matrix_setElement_BOOL(ms.matrices[0], true, 0, 1);
	GrB_Matrix_setElement_BOOL(ms.matrices[1], true, 1, 0);

	MRGraph out = {0};
	GrB_Info info = build_mr_graph(ms.matrices, &ms.symbol_list, n,
								   &DefaultTerminalFormat, &out);
	assert(info == GrB_SUCCESS);

	assert(out.n_par == 1);
	assert(out.n_bra == 0);
	assert(out.normal == NULL);
	assert(out.open_par[0] == ms.matrices[0]);
	assert(out.close_par[0] == ms.matrices[1]);
	assert(out.open_bra == NULL);
	assert(out.close_bra == NULL);
	assert(matrix_has_edge(out.open_par[0], 0, 1));
	assert(matrix_has_edge(out.close_par[0], 1, 0));

	free_mr_graph_arrays(&out);
	free_matrix_set(&ms);
}

// One open_bra and one close_bra are correctly assigned
static void test_single_bracket_pair(void) {
	GrB_Index n = 3;
	char *labels[] = {"ob_0", "cb_0"};
	MatrixSet ms = make_matrix_set(n, labels, 2);

	GrB_Matrix_setElement_BOOL(ms.matrices[0], true, 0, 2);
	GrB_Matrix_setElement_BOOL(ms.matrices[1], true, 2, 0);

	MRGraph out = {0};
	GrB_Info info = build_mr_graph(ms.matrices, &ms.symbol_list, n,
								   &DefaultTerminalFormat, &out);
	assert(info == GrB_SUCCESS);

	assert(out.n_par == 0);
	assert(out.n_bra == 1);
	assert(out.normal == NULL);
	assert(out.open_par == NULL);
	assert(out.close_par == NULL);
	assert(out.open_bra[0] == ms.matrices[0]);
	assert(out.close_bra[0] == ms.matrices[1]);
	assert(matrix_has_edge(out.open_bra[0], 0, 2));
	assert(matrix_has_edge(out.close_bra[0], 2, 0));

	free_mr_graph_arrays(&out);
	free_matrix_set(&ms);
}

// Multiple parenthesis pairs are assigned to correct slots
static void test_multiple_parenthesis_pairs(void) {
	GrB_Index n = 5;
	char *labels[] = {"op_0", "op_1", "cp_0", "cp_1"};
	MatrixSet ms = make_matrix_set(n, labels, 4);

	GrB_Matrix_setElement_BOOL(ms.matrices[0], true, 0, 1);
	GrB_Matrix_setElement_BOOL(ms.matrices[1], true, 2, 3);
	GrB_Matrix_setElement_BOOL(ms.matrices[2], true, 1, 0);
	GrB_Matrix_setElement_BOOL(ms.matrices[3], true, 3, 2);

	MRGraph out = {0};
	GrB_Info info = build_mr_graph(ms.matrices, &ms.symbol_list, n,
								   &DefaultTerminalFormat, &out);
	assert(info == GrB_SUCCESS);

	assert(out.n_par == 2);
	assert(out.n_bra == 0);

	// Find which slot has each edge
	int slot0 = -1, slot1 = -1;
	for (int64_t i = 0; i < out.n_par; i++) {
		if (matrix_has_edge(out.open_par[i], 0, 1)) {
			slot0 = i;
		}
		if (matrix_has_edge(out.open_par[i], 2, 3)) {
			slot1 = i;
		}
	}
	assert(slot0 != -1 && slot1 != -1 && slot0 != slot1);
	assert(matrix_has_edge(out.close_par[slot0], 1, 0));
	assert(matrix_has_edge(out.close_par[slot1], 3, 2));

	assert(out.normal == NULL);
	assert(out.open_bra == NULL);
	assert(out.close_bra == NULL);

	free_mr_graph_arrays(&out);
	free_matrix_set(&ms);
}

// Normal matrix is correctly assigned
static void test_normal_matrix(void) {
	GrB_Index n = 3;
	char *labels[] = {"normal"};
	MatrixSet ms = make_matrix_set(n, labels, 1);

	GrB_Matrix_setElement_BOOL(ms.matrices[0], true, 0, 2);

	MRGraph out = {0};
	GrB_Info info = build_mr_graph(ms.matrices, &ms.symbol_list, n,
								   &DefaultTerminalFormat, &out);
	assert(info == GrB_SUCCESS);

	assert(out.n_par == 0);
	assert(out.n_bra == 0);
	assert(out.open_par == NULL);
	assert(out.close_par == NULL);
	assert(out.open_bra == NULL);
	assert(out.close_bra == NULL);
	assert(out.normal == ms.matrices[0]);
	assert(matrix_has_edge(out.normal, 0, 2));

	free_mr_graph_arrays(&out);
	free_matrix_set(&ms);
}

// Empty closure matrix does not register a pair
static void test_empty_close_registers_pair(void) {
	GrB_Index n = 3;
	char *labels[] = {"op_0", "cp_0"};
	MatrixSet ms = make_matrix_set(n, labels, 2);

	// cp_0 is empty, op_0 has an edge
	GrB_Matrix_setElement_BOOL(ms.matrices[0], true, 1, 0);

	MRGraph out = {0};
	GrB_Info info = build_mr_graph(ms.matrices, &ms.symbol_list, n,
								   &DefaultTerminalFormat, &out);
	assert(info == GrB_SUCCESS);

	assert(out.n_par == 0);
	assert(out.open_bra == NULL);
	assert(out.close_bra == NULL);
	assert(out.open_par == NULL);
	assert(out.close_bra == NULL);

	free_mr_graph_arrays(&out);
	free_matrix_set(&ms);
}

// Empty open matrix does not register the pair
static void test_empty_open_does_not_register_pair(void) {
	GrB_Index n = 3;
	char *labels[] = {"op_0", "cp_0"};
	MatrixSet ms = make_matrix_set(n, labels, 2);

	// op_0 is empty, cp_0 has an edge
	GrB_Matrix_setElement_BOOL(ms.matrices[1], true, 1, 0);

	MRGraph out = {0};
	GrB_Info info = build_mr_graph(ms.matrices, &ms.symbol_list, n,
								   &DefaultTerminalFormat, &out);
	assert(info == GrB_SUCCESS);

	assert(out.n_par == 0);
	assert(out.open_par == NULL);
	assert(out.close_par == NULL);
	assert(out.open_bra == NULL);
	assert(out.close_bra == NULL);
	assert(out.normal == NULL);

	free_mr_graph_arrays(&out);
	free_matrix_set(&ms);
}

// Mixed parentheses, brackets and normal all together
static void test_mixed_all_types(void) {
	GrB_Index n = 4;
	char *labels[] = {"op_0", "cp_0", "ob_0", "cb_0", "normal"};
	MatrixSet ms = make_matrix_set(n, labels, 5);

	GrB_Matrix_setElement_BOOL(ms.matrices[0], true, 0, 1);
	GrB_Matrix_setElement_BOOL(ms.matrices[1], true, 1, 0);
	GrB_Matrix_setElement_BOOL(ms.matrices[2], true, 0, 2);
	GrB_Matrix_setElement_BOOL(ms.matrices[3], true, 2, 0);
	GrB_Matrix_setElement_BOOL(ms.matrices[4], true, 0, 3);

	MRGraph out = {0};
	GrB_Info info = build_mr_graph(ms.matrices, &ms.symbol_list, n,
								   &DefaultTerminalFormat, &out);
	assert(info == GrB_SUCCESS);

	assert(out.n_par == 1);
	assert(out.n_bra == 1);
	assert(out.normal == ms.matrices[4]);
	assert(matrix_has_edge(out.open_par[0], 0, 1));
	assert(matrix_has_edge(out.close_par[0], 1, 0));
	assert(matrix_has_edge(out.open_bra[0], 0, 2));
	assert(matrix_has_edge(out.close_bra[0], 2, 0));
	assert(matrix_has_edge(out.normal, 0, 3));

	free_mr_graph_arrays(&out);
	free_matrix_set(&ms);
}

int main(void) {
	LAGraph_Init(msg);

	test_null_arguments();
	test_empty_symbol_list();
	test_single_parenthesis_pair();
	test_single_bracket_pair();
	test_multiple_parenthesis_pairs();
	test_normal_matrix();
	test_empty_close_registers_pair();
	test_empty_open_does_not_register_pair();
	test_mixed_all_types();

	LAGraph_Finalize(msg);
	return 0;
}
