#include <LAGraph.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "cfl_idr.h"
#include "graph/valueflow_extensions.h"

static char msg[LAGRAPH_MSG_LEN];

// Utils

static IdrGraph init_graph(GrB_Index n, int64_t n_bra) {
	IdrGraph g;
	g.n = n;
	g.n_par = 0;
	g.open_par = NULL;
	g.close_par = NULL;
	g.n_bra = n_bra;
	g.open_bra = calloc(n_bra, sizeof(GrB_Matrix));
	g.close_bra = calloc(n_bra, sizeof(GrB_Matrix));
	for (int64_t i = 0; i < n_bra; i++) {
		GrB_Matrix_new(&g.open_bra[i], GrB_BOOL, n, n);
		GrB_Matrix_new(&g.close_bra[i], GrB_BOOL, n, n);
	}
	g.normal = NULL;
	return g;
}
static GrB_Matrix make_paths(GrB_Index n, const GrB_Index *rows,
							 const GrB_Index *cols, GrB_Index count) {
	GrB_Matrix paths;
	GrB_Matrix_new(&paths, GrB_BOOL, n, n);
	for (GrB_Index i = 0; i < count; i++) {
		GrB_Matrix_setElement_BOOL(paths, true, rows[i], cols[i]);
	}
	GrB_Matrix_wait(paths, GrB_MATERIALIZE);
	return paths;
}

static GrB_Index count_matrix_entries(GrB_Matrix m) {
	GrB_Index nvals = 0;
	GrB_Matrix_nvals(&nvals, m);
	return nvals;
}

static bool matrix_has_entry(GrB_Matrix m, GrB_Index row, GrB_Index col) {
	bool val = false;
	GrB_Info info = GrB_Matrix_extractElement_BOOL(&val, m, row, col);
	return info == GrB_SUCCESS && val;
}

// Tests

// Empty paths returns empty result
static void test_empty_paths(void) {
	IdrGraph g = init_graph(2, 1);

	GrB_Matrix paths = make_paths(2, NULL, NULL, 0);
	GrB_Matrix filtered = NULL;

	filter_bracket_paths(&g, paths, &filtered);

	assert(count_matrix_entries(filtered) == 0);

	GrB_Matrix_free(&paths);
	GrB_Matrix_free(&filtered);
	idr_graph_free(&g);
}

// No matching brackets removes all paths
static void test_no_matching_brackets(void) {
	// A(0) -normal-> B(1), path A->B
	IdrGraph g = init_graph(2, 1);
	GrB_Matrix_new(&g.normal, GrB_BOOL, 2, 2);
	GrB_Matrix_setElement_BOOL(g.normal, true, 0, 1);

	GrB_Index rows[] = {0}, cols[] = {1};
	GrB_Matrix paths = make_paths(2, rows, cols, 1);
	GrB_Matrix filtered = NULL;

	filter_bracket_paths(&g, paths, &filtered);

	assert(count_matrix_entries(filtered) == 0);

	GrB_Matrix_free(&paths);
	GrB_Matrix_free(&filtered);
	idr_graph_free(&g);
}

// Simple matching bracket path A -ob-> B -cb-> C, path A->C preserved
static void test_simple_matching_bracket(void) {
	// A(0) -ob-> B(1) -cb-> C(2)
	IdrGraph g = init_graph(3, 1);
	GrB_Matrix_setElement_BOOL(g.open_bra[0], true, 0, 1);
	GrB_Matrix_setElement_BOOL(g.close_bra[0], true, 1, 2);

	GrB_Index rows[] = {0}, cols[] = {2};
	GrB_Matrix paths = make_paths(3, rows, cols, 1);
	GrB_Matrix filtered = NULL;

	filter_bracket_paths(&g, paths, &filtered);

	assert(count_matrix_entries(filtered) == 1);
	assert(matrix_has_entry(filtered, 0, 2));

	GrB_Matrix_free(&paths);
	GrB_Matrix_free(&filtered);
	idr_graph_free(&g);
}

// Mismatched bracket ids filtered out
static void test_mismatched_bracket_ids(void) {
	// A(0) -ob[0]-> B(1) -cb[1]-> C(2)
	IdrGraph g = init_graph(3, 2);
	GrB_Matrix_setElement_BOOL(g.open_bra[0], true, 0, 1);
	GrB_Matrix_setElement_BOOL(g.close_bra[1], true, 1, 2);

	GrB_Index rows[] = {0}, cols[] = {2};
	GrB_Matrix paths = make_paths(3, rows, cols, 1);
	GrB_Matrix filtered = NULL;

	filter_bracket_paths(&g, paths, &filtered);

	assert(count_matrix_entries(filtered) == 0);

	GrB_Matrix_free(&paths);
	GrB_Matrix_free(&filtered);
	idr_graph_free(&g);
}

// Open bracket exists but no close bracket - filtered out
static void test_no_close_bracket(void) {
	// A(0) -ob-> B(1)
	IdrGraph g = init_graph(2, 1);
	GrB_Matrix_setElement_BOOL(g.open_bra[0], true, 0, 1);

	GrB_Index rows[] = {0}, cols[] = {1};
	GrB_Matrix paths = make_paths(2, rows, cols, 1);
	GrB_Matrix filtered = NULL;

	filter_bracket_paths(&g, paths, &filtered);

	assert(count_matrix_entries(filtered) == 0);

	GrB_Matrix_free(&paths);
	GrB_Matrix_free(&filtered);
	idr_graph_free(&g);
}

// Brackets exist but no reachability between them - filtered out
static void test_no_reachability_between_brackets(void) {
	// A(0) -ob-> B(1), C(2) -cb-> D(3)
	IdrGraph g = init_graph(4, 1);
	GrB_Matrix_setElement_BOOL(g.open_bra[0], true, 0, 1);
	GrB_Matrix_setElement_BOOL(g.close_bra[0], true, 2, 3);

	GrB_Index rows[] = {0}, cols[] = {3};
	GrB_Matrix paths = make_paths(4, rows, cols, 1);
	GrB_Matrix filtered = NULL;

	filter_bracket_paths(&g, paths, &filtered);

	assert(count_matrix_entries(filtered) == 0);

	GrB_Matrix_free(&paths);
	GrB_Matrix_free(&filtered);
	idr_graph_free(&g);
}

// Bracket path through intermediate normal edges accepted
static void test_path_through_normal_edges(void) {
	// A(0) -ob-> B(1) -normal-> C(2) -cb-> D(3)
	IdrGraph g = init_graph(4, 1);
	GrB_Matrix_new(&g.normal, GrB_BOOL, 4, 4);
	GrB_Matrix_setElement_BOOL(g.open_bra[0], true, 0, 1);
	GrB_Matrix_setElement_BOOL(g.normal, true, 1, 2);
	GrB_Matrix_setElement_BOOL(g.close_bra[0], true, 2, 3);

	GrB_Index rows[] = {0}, cols[] = {3};
	GrB_Matrix paths = make_paths(4, rows, cols, 1);
	GrB_Matrix filtered = NULL;

	filter_bracket_paths(&g, paths, &filtered);

	assert(count_matrix_entries(filtered) == 1);
	assert(matrix_has_entry(filtered, 0, 3));

	GrB_Matrix_free(&paths);
	GrB_Matrix_free(&filtered);
	idr_graph_free(&g);
}

// Multiple valid bracket paths all preserved
static void test_multiple_valid_paths(void) {
	// A(0) -ob[0]-> B(1) -cb[0]-> C(2)
	// A(0) -ob[1]-> D(3) -cb[1]-> E(4)
	IdrGraph g = init_graph(5, 2);
	GrB_Matrix_setElement_BOOL(g.open_bra[0], true, 0, 1);
	GrB_Matrix_setElement_BOOL(g.close_bra[0], true, 1, 2);
	GrB_Matrix_setElement_BOOL(g.open_bra[1], true, 0, 3);
	GrB_Matrix_setElement_BOOL(g.close_bra[1], true, 3, 4);

	GrB_Index rows[] = {0, 0}, cols[] = {2, 4};
	GrB_Matrix paths = make_paths(5, rows, cols, 2);
	GrB_Matrix filtered = NULL;

	filter_bracket_paths(&g, paths, &filtered);

	assert(count_matrix_entries(filtered) == 2);
	assert(matrix_has_entry(filtered, 0, 2));
	assert(matrix_has_entry(filtered, 0, 4));

	GrB_Matrix_free(&paths);
	GrB_Matrix_free(&filtered);
	idr_graph_free(&g);
}

// Multiple open/close edges with same bracket id, path accepted
static void test_multiple_open_close_same_id(void) {
	// A(0) -ob-> B(1), A(0) -ob-> C(2)
	// B(1) -cb-> D(3), C(2) -cb-> D(3)
	IdrGraph g = init_graph(4, 1);
	GrB_Matrix_setElement_BOOL(g.open_bra[0], true, 0, 1);
	GrB_Matrix_setElement_BOOL(g.open_bra[0], true, 0, 2);
	GrB_Matrix_setElement_BOOL(g.close_bra[0], true, 1, 3);
	GrB_Matrix_setElement_BOOL(g.close_bra[0], true, 2, 3);

	GrB_Index rows[] = {0}, cols[] = {3};
	GrB_Matrix paths = make_paths(4, rows, cols, 1);
	GrB_Matrix filtered = NULL;

	filter_bracket_paths(&g, paths, &filtered);

	assert(count_matrix_entries(filtered) == 1);
	assert(matrix_has_entry(filtered, 0, 3));

	GrB_Matrix_free(&paths);
	GrB_Matrix_free(&filtered);
	idr_graph_free(&g);
}

// Nested brackets: A->C, C->E, A->E all preserved
static void test_nested_brackets(void) {
	// A(0) -ob-> B(1) -cb-> C(2) -ob-> D(3) -cb-> E(4)
	IdrGraph g = init_graph(5, 1);
	GrB_Matrix_setElement_BOOL(g.open_bra[0], true, 0, 1);
	GrB_Matrix_setElement_BOOL(g.close_bra[0], true, 1, 2);
	GrB_Matrix_setElement_BOOL(g.open_bra[0], true, 2, 3);
	GrB_Matrix_setElement_BOOL(g.close_bra[0], true, 3, 4);

	GrB_Index rows[] = {0, 2, 0}, cols[] = {2, 4, 4};
	GrB_Matrix paths = make_paths(5, rows, cols, 3);
	GrB_Matrix filtered = NULL;

	filter_bracket_paths(&g, paths, &filtered);

	assert(count_matrix_entries(filtered) == 3);
	assert(matrix_has_entry(filtered, 0, 2));
	assert(matrix_has_entry(filtered, 2, 4));
	assert(matrix_has_entry(filtered, 0, 4));

	GrB_Matrix_free(&paths);
	GrB_Matrix_free(&filtered);
	idr_graph_free(&g);
}

int main(void) {
	LAGraph_Init(msg);

	test_empty_paths();
	test_no_matching_brackets();
	test_simple_matching_bracket();
	test_mismatched_bracket_ids();
	test_no_close_bracket();
	test_no_reachability_between_brackets();
	test_path_through_normal_edges();
	test_multiple_valid_paths();
	test_multiple_open_close_same_id();
	test_nested_brackets();

	LAGraph_Finalize(msg);
	return 0;
}
