
#include <LAGraph.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cfl_idr.h"
#include "graph/condensate_graph.h"

static char msg[LAGRAPH_MSG_LEN];

// Utils

static GrB_Matrix make_empty_matrix(GrB_Index n) {
	GrB_Matrix m;
	GrB_Matrix_new(&m, GrB_BOOL, n, n);
	return m;
}

static bool matrix_has_edge(GrB_Matrix m, GrB_Index i, GrB_Index j) {
	bool val = false;
	GrB_Info info = GrB_Matrix_extractElement_BOOL(&val, m, i, j);
	return info == GrB_SUCCESS && val;
}

static GrB_Index matrix_nvals(GrB_Matrix m) {
	GrB_Index nvals = 0;
	GrB_Matrix_nvals(&nvals, m);
	return nvals;
}

static GrB_Index idr_graph_nvals(const IdrGraph *g) {
	GrB_Index total = 0;
	for (int64_t i = 0; i < g->n_par; i++) {
		total += matrix_nvals(g->open_par[i]);
		total += matrix_nvals(g->close_par[i]);
	}
	for (int64_t i = 0; i < g->n_bra; i++) {
		total += matrix_nvals(g->open_bra[i]);
		total += matrix_nvals(g->close_bra[i]);
	}
	if (g->normal != NULL)
		total += matrix_nvals(g->normal);
	return total;
}

static GrB_Index get_rep(GrB_Vector components, GrB_Index i) {
	uint64_t rep = 0;
	GrB_Vector_extractElement_UINT64(&rep, components, i);
	return (GrB_Index)rep;
}

static IdrGraph make_simple_graph(GrB_Index n, GrB_Matrix open_par,
								  GrB_Matrix close_par) {
	GrB_Matrix *op = malloc(sizeof(GrB_Matrix));
	GrB_Matrix *cp = malloc(sizeof(GrB_Matrix));
	op[0] = open_par;
	cp[0] = close_par;

	IdrGraph g = {.open_par = op,
				  .close_par = cp,
				  .n_par = 1,
				  .open_bra = NULL,
				  .close_bra = NULL,
				  .n_bra = 0,
				  .normal = NULL,
				  .n = n};
	return g;
}

static void free_simple_graph_arrays(IdrGraph *g) {
	free(g->open_par);
	free(g->close_par);
}

// Tests

static void test_empty_graph(void) {
	GrB_Matrix op = make_empty_matrix(3);
	GrB_Matrix cp = make_empty_matrix(3);
	IdrGraph graph = make_simple_graph(3, op, cp);

	CondensationResult res = {0};
	GrB_Info info = condensate_from_under_approx(&graph, NULL, &res, msg);
	assert(info == GrB_SUCCESS);

	// With an empty under_approx, its own vertex is its own representative
	assert(get_rep(res.components, 0) == 0);
	assert(get_rep(res.components, 1) == 1);
	assert(get_rep(res.components, 2) == 2);

	assert(idr_graph_nvals(&res.condensed_graph) == 0);

	condensation_result_free(&res, msg);
	free_simple_graph_arrays(&graph);
	GrB_Matrix_free(&op);
	GrB_Matrix_free(&cp);
}

// Mutual paths cause vertices to be merged
static void test_mutual_paths_cause_merging(void) {
	GrB_Matrix op = make_empty_matrix(2);
	GrB_Matrix cp = make_empty_matrix(2);
	GrB_Matrix_setElement_BOOL(op, true, 0, 1);
	GrB_Matrix_setElement_BOOL(cp, true, 1, 0);

	// Graph: 0->1 (open), 1->0 (close)
	IdrGraph graph = make_simple_graph(2, op, cp);

	// under_approx: (0->1) and (1->0)
	GrB_Matrix under = make_empty_matrix(2);
	GrB_Matrix_setElement_BOOL(under, true, 0, 1);
	GrB_Matrix_setElement_BOOL(under, true, 1, 0);

	CondensationResult cr = {0};
	GrB_Info info = condensate_from_under_approx(&graph, under, &cr, msg);
	assert(info == GrB_SUCCESS);

	GrB_Index rep0 = get_rep(cr.components, 0);
	GrB_Index rep1 = get_rep(cr.components, 1);
	assert(rep0 == rep1);

	assert(idr_graph_nvals(&cr.condensed_graph) == 2);

	condensation_result_free(&cr, msg);
	free_simple_graph_arrays(&graph);
	GrB_Matrix_free(&op);
	GrB_Matrix_free(&cp);
	GrB_Matrix_free(&under);
}

// One-way path does not cause merging
static void test_one_way_path_no_merge(void) {
	GrB_Matrix op = make_empty_matrix(2);
	GrB_Matrix cp = make_empty_matrix(2);
	GrB_Matrix_setElement_BOOL(op, true, 0, 1);
	GrB_Matrix_setElement_BOOL(cp, true, 1, 0);

	// Graph: 0->1 (open), 1->0 (close)
	IdrGraph graph = make_simple_graph(2, op, cp);

	// under_approx: only (0→1)
	GrB_Matrix under = make_empty_matrix(2);
	GrB_Matrix_setElement_BOOL(under, true, 0, 1);

	CondensationResult cr = {0};
	GrB_Info info = condensate_from_under_approx(&graph, under, &cr, msg);
	assert(info == GrB_SUCCESS);

	GrB_Index rep0 = get_rep(cr.components, 0);
	GrB_Index rep1 = get_rep(cr.components, 1);
	assert(rep0 != rep1);

	assert(matrix_has_edge(cr.condensed_graph.open_par[0], rep0, rep1));

	condensation_result_free(&cr, msg);
	free_simple_graph_arrays(&graph);
	GrB_Matrix_free(&op);
	GrB_Matrix_free(&cp);
	GrB_Matrix_free(&under);
}

// Three vertices with mutual pairs merged correctly
static void test_three_vertices_all_merged(void) {
	GrB_Matrix op = make_empty_matrix(3);
	GrB_Matrix cp = make_empty_matrix(3);
	GrB_Matrix_setElement_BOOL(op, true, 0, 1);
	GrB_Matrix_setElement_BOOL(op, true, 1, 2);
	GrB_Matrix_setElement_BOOL(cp, true, 1, 0);
	GrB_Matrix_setElement_BOOL(cp, true, 2, 1);

	// Graph: 0->1, 1->2 (open), 1->0, 2->1 (close)
	IdrGraph graph = make_simple_graph(3, op, cp);

	// under_approx: A<->B and B<->C
	GrB_Matrix under = make_empty_matrix(3);
	GrB_Matrix_setElement_BOOL(under, true, 0, 1);
	GrB_Matrix_setElement_BOOL(under, true, 1, 0);
	GrB_Matrix_setElement_BOOL(under, true, 1, 2);
	GrB_Matrix_setElement_BOOL(under, true, 2, 1);

	CondensationResult cr = {0};
	GrB_Info info = condensate_from_under_approx(&graph, under, &cr, msg);
	assert(info == GrB_SUCCESS);

	GrB_Index rep0 = get_rep(cr.components, 0);
	GrB_Index rep1 = get_rep(cr.components, 1);
	GrB_Index rep2 = get_rep(cr.components, 2);
	assert(rep0 == rep1);
	assert(rep1 == rep2);

	assert(idr_graph_nvals(&cr.condensed_graph) == 2);

	condensation_result_free(&cr, msg);
	free_simple_graph_arrays(&graph);
	GrB_Matrix_free(&op);
	GrB_Matrix_free(&cp);
	GrB_Matrix_free(&under);
}

// Partial mutual paths - only fully mutual groups are merged
static void test_partial_mutual_paths(void) {
	GrB_Matrix op = make_empty_matrix(3);
	GrB_Matrix cp = make_empty_matrix(3);
	GrB_Matrix_setElement_BOOL(op, true, 0, 1);
	GrB_Matrix_setElement_BOOL(op, true, 1, 2);
	GrB_Matrix_setElement_BOOL(cp, true, 1, 0);

	// Graph: 0->1, 1->2 (open), 1->0 (close)
	IdrGraph graph = make_simple_graph(3, op, cp);

	// under_approx: 0<->1, but 1 and 2 not mutually reachable
	GrB_Matrix under = make_empty_matrix(3);
	GrB_Matrix_setElement_BOOL(under, true, 0, 1);
	GrB_Matrix_setElement_BOOL(under, true, 1, 0);
	GrB_Matrix_setElement_BOOL(under, true, 1, 2);

	CondensationResult cr = {0};
	GrB_Info info = condensate_from_under_approx(&graph, under, &cr, msg);
	assert(info == GrB_SUCCESS);

	GrB_Index rep0 = get_rep(cr.components, 0);
	GrB_Index rep1 = get_rep(cr.components, 1);
	GrB_Index rep2 = get_rep(cr.components, 2);
	assert(rep0 == rep1);
	assert(rep1 != rep2);

	assert(idr_graph_nvals(&cr.condensed_graph) == 3);

	condensation_result_free(&cr, msg);
	free_simple_graph_arrays(&graph);
	GrB_Matrix_free(&op);
	GrB_Matrix_free(&cp);
	GrB_Matrix_free(&under);
}

// Duplicate edges are not added to condensed graph
static void test_edges_between_components_deduplicated(void) {
	GrB_Matrix op = make_empty_matrix(4);
	GrB_Matrix cp = make_empty_matrix(4);
	GrB_Matrix_setElement_BOOL(op, true, 0, 2);
	GrB_Matrix_setElement_BOOL(op, true, 1, 2);
	GrB_Matrix_setElement_BOOL(cp, true, 0, 3);

	// Graph: 0->2, 1->2 (open), 0->3 (close)
	IdrGraph graph = make_simple_graph(4, op, cp);

	// under_approx: 0<->1
	GrB_Matrix under = make_empty_matrix(4);
	GrB_Matrix_setElement_BOOL(under, true, 0, 1);
	GrB_Matrix_setElement_BOOL(under, true, 1, 0);

	CondensationResult cr = {0};
	GrB_Info info = condensate_from_under_approx(&graph, under, &cr, msg);
	assert(info == GrB_SUCCESS);

	GrB_Index rep0 = get_rep(cr.components, 0);
	GrB_Index rep1 = get_rep(cr.components, 1);
	GrB_Index rep2 = get_rep(cr.components, 2);
	GrB_Index rep3 = get_rep(cr.components, 3);
	assert(rep0 == rep1);
	assert(rep0 != rep2);
	assert(rep0 != rep3);

	assert(idr_graph_nvals(&cr.condensed_graph) == 2);

	condensation_result_free(&cr, msg);
	free_simple_graph_arrays(&graph);
	GrB_Matrix_free(&op);
	GrB_Matrix_free(&cp);
	GrB_Matrix_free(&under);
}

// Edges between different components are preserved once
static void test_null_under_approx(void) {
	GrB_Matrix op = make_empty_matrix(3);
	GrB_Matrix cp = make_empty_matrix(3);
	GrB_Matrix_setElement_BOOL(op, true, 0, 1);
	GrB_Matrix_setElement_BOOL(cp, true, 1, 2);

	// Graph: 0->1, 1->2 (open)
	IdrGraph graph = make_simple_graph(3, op, cp);

	CondensationResult cr = {0};
	GrB_Info info = condensate_from_under_approx(&graph, NULL, &cr, msg);
	assert(info == GrB_SUCCESS);

	assert(get_rep(cr.components, 0) == 0);
	assert(get_rep(cr.components, 1) == 1);
	assert(get_rep(cr.components, 2) == 2);

	assert(idr_graph_nvals(&cr.condensed_graph) == 2);
	assert(matrix_has_edge(cr.condensed_graph.open_par[0], 0, 1));
	assert(matrix_has_edge(cr.condensed_graph.close_par[0], 1, 2));

	condensation_result_free(&cr, msg);
	free_simple_graph_arrays(&graph);
	GrB_Matrix_free(&op);
	GrB_Matrix_free(&cp);
}

int main(void) {
	LAGraph_Init(msg);

	test_empty_graph();
	test_mutual_paths_cause_merging();
	test_one_way_path_no_merge();
	test_three_vertices_all_merged();
	test_partial_mutual_paths();
	test_edges_between_components_deduplicated();
	test_null_under_approx();

	LAGraph_Finalize(msg);
	return 0;
}
