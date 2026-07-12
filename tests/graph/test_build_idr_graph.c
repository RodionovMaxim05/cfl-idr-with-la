#include <LAGraph.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cfl_idr.h"
#include "cfl_idr_graph_builder.h"
#include "symbol_list.h"

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
	GraphMatrices gm;
	SymbolList symbol_list;
} MatrixSet;

static MatrixSet make_matrix_set(GrB_Index n, char **labels, size_t count) {
	MatrixSet ms;
	ms.symbol_list = symbol_list_create();

	GrB_Matrix *matrices = malloc(count * sizeof(GrB_Matrix));
	MatrixSymbolInfo *symbols = malloc(count * sizeof(MatrixSymbolInfo));

	for (size_t i = 0; i < count; i++) {
		matrices[i] = make_empty_matrix(n);
		size_t sym_idx =
			(size_t)symbol_list_add_str(&ms.symbol_list, labels[i], false);

		size_t block_index = 0;
		const char *underscore = strrchr(labels[i], '_');
		if (underscore != NULL && *(underscore + 1) != '\0') {
			block_index = (size_t)atoi(underscore + 1);
		}

		symbols[i] =
			(MatrixSymbolInfo){.symbol_index = sym_idx, .block_index = block_index};
	}

	ms.gm = (GraphMatrices){
		.matrices = matrices,
		.matrix_symbols = symbols,
		.count = count,
	};

	return ms;
}

static void free_matrix_set(MatrixSet *ms) {
	for (size_t i = 0; i < ms->gm.count; i++) {
		GrB_Matrix_free(&ms->gm.matrices[i]);
	}
	free(ms->gm.matrices);
	free(ms->gm.matrix_symbols);
	symbol_list_free(&ms->symbol_list);
}

static void free_idr_graph_arrays(IdrGraph *g) {
	free(g->open_par);
	free(g->close_par);
	free(g->open_bra);
	free(g->close_bra);
}

// Tests

static void test_null_arguments(void) {
	IdrGraph out = {0};
	GraphMatrices gm = {0};
	SymbolList sl = symbol_list_create();

	assert(get_idr_graph(&out, NULL, &sl, 3, &DefaultTerminalFormat) ==
		   GrB_INVALID_VALUE);
	assert(get_idr_graph(&out, &gm, NULL, 3, &DefaultTerminalFormat) ==
		   GrB_INVALID_VALUE);
	assert(get_idr_graph(&out, &gm, &sl, 3, NULL) == GrB_INVALID_VALUE);
	assert(get_idr_graph(NULL, &gm, &sl, 3, &DefaultTerminalFormat) ==
		   GrB_INVALID_VALUE);

	symbol_list_free(&sl);
}

static void test_empty_symbol_list(void) {
	SymbolList sl = symbol_list_create();
	GraphMatrices gm = {.matrices = NULL, .matrix_symbols = NULL, .count = 0};
	IdrGraph out = {0};

	GrB_Info info = get_idr_graph(&out, &gm, &sl, 3, &DefaultTerminalFormat);
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

	GrB_Matrix_setElement_BOOL(ms.gm.matrices[0], true, 0, 1);
	GrB_Matrix_setElement_BOOL(ms.gm.matrices[1], true, 1, 0);

	IdrGraph out = {0};
	GrB_Info info =
		get_idr_graph(&out, &ms.gm, &ms.symbol_list, n, &DefaultTerminalFormat);

	assert(info == GrB_SUCCESS);

	assert(out.n_par == 1);
	assert(out.n_bra == 0);
	assert(out.normal == NULL);
	assert(out.open_bra == NULL);
	assert(out.close_bra == NULL);
	assert(matrix_has_edge(out.open_par[0], 0, 1));
	assert(matrix_has_edge(out.close_par[0], 1, 0));

	free_idr_graph_arrays(&out);
	free_matrix_set(&ms);
}

// One open_bra and one close_bra are correctly assigned
static void test_single_bracket_pair(void) {
	GrB_Index n = 3;
	char *labels[] = {"ob_0", "cb_0"};
	MatrixSet ms = make_matrix_set(n, labels, 2);

	GrB_Matrix_setElement_BOOL(ms.gm.matrices[0], true, 0, 2);
	GrB_Matrix_setElement_BOOL(ms.gm.matrices[1], true, 2, 0);

	IdrGraph out = {0};
	GrB_Info info =
		get_idr_graph(&out, &ms.gm, &ms.symbol_list, n, &DefaultTerminalFormat);

	assert(info == GrB_SUCCESS);

	assert(out.n_par == 0);
	assert(out.n_bra == 1);
	assert(out.normal == NULL);
	assert(out.open_par == NULL);
	assert(out.close_par == NULL);
	assert(matrix_has_edge(out.open_bra[0], 0, 2));
	assert(matrix_has_edge(out.close_bra[0], 2, 0));

	free_idr_graph_arrays(&out);
	free_matrix_set(&ms);
}

// Multiple parenthesis pairs are assigned to correct slots
static void test_multiple_parenthesis_pairs(void) {
	GrB_Index n = 5;
	char *labels[] = {"op_0", "op_1", "cp_0", "cp_1"};
	MatrixSet ms = make_matrix_set(n, labels, 4);

	GrB_Matrix_setElement_BOOL(ms.gm.matrices[0], true, 0, 1);
	GrB_Matrix_setElement_BOOL(ms.gm.matrices[1], true, 2, 3);
	GrB_Matrix_setElement_BOOL(ms.gm.matrices[2], true, 1, 0);
	GrB_Matrix_setElement_BOOL(ms.gm.matrices[3], true, 3, 2);

	IdrGraph out = {0};
	GrB_Info info =
		get_idr_graph(&out, &ms.gm, &ms.symbol_list, n, &DefaultTerminalFormat);
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

	free_idr_graph_arrays(&out);
	free_matrix_set(&ms);
}

// Normal matrix is correctly assigned
static void test_normal_matrix(void) {
	GrB_Index n = 3;
	char *labels[] = {"normal"};
	MatrixSet ms = make_matrix_set(n, labels, 1);

	GrB_Matrix_setElement_BOOL(ms.gm.matrices[0], true, 0, 2);

	IdrGraph out = {0};
	GrB_Info info =
		get_idr_graph(&out, &ms.gm, &ms.symbol_list, n, &DefaultTerminalFormat);
	assert(info == GrB_SUCCESS);

	assert(out.n_par == 0);
	assert(out.n_bra == 0);
	assert(out.open_par == NULL);
	assert(out.close_par == NULL);
	assert(out.open_bra == NULL);
	assert(out.close_bra == NULL);
	assert(matrix_has_edge(out.normal, 0, 2));

	free_idr_graph_arrays(&out);
	free_matrix_set(&ms);
}

// Empty closure matrix does not register a pair
static void test_empty_close_registers_pair(void) {
	GrB_Index n = 3;
	char *labels[] = {"op_0", "cp_0"};
	MatrixSet ms = make_matrix_set(n, labels, 2);

	// cp_0 is empty, op_0 has an edge
	GrB_Matrix_setElement_BOOL(ms.gm.matrices[0], true, 1, 0);

	IdrGraph out = {0};
	GrB_Info info =
		get_idr_graph(&out, &ms.gm, &ms.symbol_list, n, &DefaultTerminalFormat);
	assert(info == GrB_SUCCESS);

	assert(out.n_par == 0);
	assert(out.open_bra == NULL);
	assert(out.close_bra == NULL);
	assert(out.open_par == NULL);
	assert(out.close_bra == NULL);

	free_idr_graph_arrays(&out);
	free_matrix_set(&ms);
}

// Empty open matrix does not register the pair
static void test_empty_open_does_not_register_pair(void) {
	GrB_Index n = 3;
	char *labels[] = {"op_0", "cp_0"};
	MatrixSet ms = make_matrix_set(n, labels, 2);

	// op_0 is empty, cp_0 has an edge
	GrB_Matrix_setElement_BOOL(ms.gm.matrices[1], true, 1, 0);

	IdrGraph out = {0};
	GrB_Info info =
		get_idr_graph(&out, &ms.gm, &ms.symbol_list, n, &DefaultTerminalFormat);
	assert(info == GrB_SUCCESS);

	assert(out.n_par == 0);
	assert(out.open_par == NULL);
	assert(out.close_par == NULL);
	assert(out.open_bra == NULL);
	assert(out.close_bra == NULL);
	assert(out.normal == NULL);

	free_idr_graph_arrays(&out);
	free_matrix_set(&ms);
}

// Mixed parentheses, brackets and normal all together
static void test_mixed_all_types(void) {
	GrB_Index n = 4;
	char *labels[] = {"op_0", "cp_0", "ob_0", "cb_0", "normal"};
	MatrixSet ms = make_matrix_set(n, labels, 5);

	GrB_Matrix_setElement_BOOL(ms.gm.matrices[0], true, 0, 1);
	GrB_Matrix_setElement_BOOL(ms.gm.matrices[1], true, 1, 0);
	GrB_Matrix_setElement_BOOL(ms.gm.matrices[2], true, 0, 2);
	GrB_Matrix_setElement_BOOL(ms.gm.matrices[3], true, 2, 0);
	GrB_Matrix_setElement_BOOL(ms.gm.matrices[4], true, 0, 3);

	IdrGraph out = {0};
	GrB_Info info =
		get_idr_graph(&out, &ms.gm, &ms.symbol_list, n, &DefaultTerminalFormat);
	assert(info == GrB_SUCCESS);

	assert(out.n_par == 1);
	assert(out.n_bra == 1);
	assert(matrix_has_edge(out.open_par[0], 0, 1));
	assert(matrix_has_edge(out.close_par[0], 1, 0));
	assert(matrix_has_edge(out.open_bra[0], 0, 2));
	assert(matrix_has_edge(out.close_bra[0], 2, 0));
	assert(matrix_has_edge(out.normal, 0, 3));

	free_idr_graph_arrays(&out);
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
