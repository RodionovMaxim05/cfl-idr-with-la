#include "LAGraph.h"
#include "LAGraphX.h"

#include "graph/remove_not_path.h"

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

static GrB_Index mrgraph_nvals(const MRGraph *g) {
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

static bool scc_reach_has(GrB_Matrix scc_reach, GrB_Index i, GrB_Index j) {
	bool val = false;
	GrB_Info info = GrB_Matrix_extractElement_BOOL(&val, scc_reach, i, j);
	return info == GrB_SUCCESS && val;
}

static MRGraph make_simple_graph(GrB_Index n, GrB_Matrix open_par,
								 GrB_Matrix close_par) {
	GrB_Matrix *op = malloc(sizeof(GrB_Matrix));
	GrB_Matrix *cp = malloc(sizeof(GrB_Matrix));
	op[0] = open_par;
	cp[0] = close_par;

	MRGraph g = {
		.n = n,
		.n_par = 1,
		.open_par = op,
		.close_par = cp,
		.n_bra = 0,
		.open_bra = NULL,
		.close_bra = NULL,
		.normal = NULL,
	};
	return g;
}

static void free_simple_graph_arrays(MRGraph *g) {
	free(g->open_par);
	free(g->close_par);
}

static GrB_Matrix make_over_approx(GrB_Index n, const GrB_Index *srcs,
								   const GrB_Index *dsts, GrB_Index count) {
	GrB_Matrix m;
	GrB_Matrix_new(&m, GrB_BOOL, n, n);
	for (GrB_Index k = 0; k < count; k++)
		GrB_Matrix_setElement_BOOL(m, true, srcs[k], dsts[k]);
	return m;
}

// =================================================================================
// Tests for compute_sccs
// =================================================================================

// single vertex with no edges -> single SCC, reflexive pair
static void test_computeSccs_single_vertex(void) {
	GrB_Matrix op = make_empty_matrix(1);
	GrB_Matrix cp = make_empty_matrix(1);
	MRGraph graph = make_simple_graph(1, op, cp);

	SccResult sr = {0};
	GrB_Info info = compute_sccs(&graph, &sr, msg);
	assert(info == GrB_SUCCESS);

	assert(sr.n_scc == 1);
	assert(sr.scc_ids[0] == 0);
	assert(scc_reach_has(sr.scc_reach, 0, 0));
	assert(matrix_nvals(sr.scc_reach) == 1);

	scc_result_free(&sr);
	free_simple_graph_arrays(&graph);
	GrB_Matrix_free(&op);
	GrB_Matrix_free(&cp);
}

// Two isolated vertices -> 2 SCC, reflexive pairs only
static void test_computeSccs_two_disconnected_vertices(void) {
	GrB_Matrix op = make_empty_matrix(2);
	GrB_Matrix cp = make_empty_matrix(2);
	MRGraph graph = make_simple_graph(2, op, cp);

	SccResult sr = {0};
	GrB_Info info = compute_sccs(&graph, &sr, msg);
	assert(info == GrB_SUCCESS);

	assert(sr.n_scc == 2);

	GrB_Index c0 = sr.scc_ids[0];
	GrB_Index c1 = sr.scc_ids[1];
	assert(c0 != c1);

	assert(scc_reach_has(sr.scc_reach, c0, c0));
	assert(scc_reach_has(sr.scc_reach, c1, c1));
	assert(!scc_reach_has(sr.scc_reach, c0, c1));
	assert(!scc_reach_has(sr.scc_reach, c1, c0));
	assert(matrix_nvals(sr.scc_reach) == 2);

	scc_result_free(&sr);
	free_simple_graph_arrays(&graph);
	GrB_Matrix_free(&op);
	GrB_Matrix_free(&cp);
}

// Directed cycle A->B->C->A -> all three vertices in one SCC
static void test_computeSccs_directed_cycle_is_one_scc(void) {
	GrB_Matrix op = make_empty_matrix(3);
	GrB_Matrix cp = make_empty_matrix(3);
	GrB_Matrix_setElement_BOOL(op, true, 0, 1); // A->B
	GrB_Matrix_setElement_BOOL(op, true, 1, 2); // B->C
	GrB_Matrix_setElement_BOOL(cp, true, 2, 0); // C->A
	MRGraph graph = make_simple_graph(3, op, cp);

	SccResult sr = {0};
	GrB_Info info = compute_sccs(&graph, &sr, msg);
	assert(info == GrB_SUCCESS);

	assert(sr.n_scc == 1);
	assert(sr.scc_ids[0] == sr.scc_ids[1]);
	assert(sr.scc_ids[1] == sr.scc_ids[2]);

	assert(matrix_nvals(sr.scc_reach) == 1);
	assert(scc_reach_has(sr.scc_reach, 0, 0));

	scc_result_free(&sr);
	free_simple_graph_arrays(&graph);
	GrB_Matrix_free(&op);
	GrB_Matrix_free(&cp);
}

// DAG: A->B (open), C->D (open), B->C (close)
static void test_computeSccs_dag_transitive_closure(void) {
	GrB_Matrix op = make_empty_matrix(4);
	GrB_Matrix cp = make_empty_matrix(4);
	GrB_Matrix_setElement_BOOL(op, true, 0, 1); // A->B
	GrB_Matrix_setElement_BOOL(op, true, 2, 3); // C->D
	GrB_Matrix_setElement_BOOL(cp, true, 1, 2); // B->C
	MRGraph graph = make_simple_graph(4, op, cp);

	SccResult sr = {0};
	GrB_Info info = compute_sccs(&graph, &sr, msg);
	assert(info == GrB_SUCCESS);

	assert(sr.n_scc == 4);

	GrB_Index ca = sr.scc_ids[0];
	GrB_Index cb = sr.scc_ids[1];
	GrB_Index cc = sr.scc_ids[2];
	GrB_Index cd = sr.scc_ids[3];

	// All in different components
	assert(ca != cb && ca != cc && ca != cd);
	assert(cb != cc && cb != cd);
	assert(cc != cd);

	// Reflexive
	assert(scc_reach_has(sr.scc_reach, ca, ca));
	assert(scc_reach_has(sr.scc_reach, cb, cb));
	assert(scc_reach_has(sr.scc_reach, cc, cc));
	assert(scc_reach_has(sr.scc_reach, cd, cd));

	// Direct edges
	assert(scc_reach_has(sr.scc_reach, ca, cb));
	assert(scc_reach_has(sr.scc_reach, cb, cc));
	assert(scc_reach_has(sr.scc_reach, cc, cd));

	// Transitive
	assert(scc_reach_has(sr.scc_reach, ca, cc));
	assert(scc_reach_has(sr.scc_reach, ca, cd));
	assert(scc_reach_has(sr.scc_reach, cb, cd));

	// No inverse
	assert(!scc_reach_has(sr.scc_reach, cd, ca));
	assert(!scc_reach_has(sr.scc_reach, cc, ca));

	assert(matrix_nvals(sr.scc_reach) == 10);

	scc_result_free(&sr);
	free_simple_graph_arrays(&graph);
	GrB_Matrix_free(&op);
	GrB_Matrix_free(&cp);
}

// Parallel edges A->B (three labels) -> 2 SCC, one edge in condensation
static void test_computeSccs_parallel_edges(void) {
	GrB_Matrix op = make_empty_matrix(2);
	GrB_Matrix cp = make_empty_matrix(2);
	GrB_Matrix_setElement_BOOL(op, true, 0, 1);
	GrB_Matrix_setElement_BOOL(op, true, 0, 1);
	GrB_Matrix_setElement_BOOL(cp, true, 0, 1);
	MRGraph graph = make_simple_graph(2, op, cp);

	SccResult sr = {0};
	GrB_Info info = compute_sccs(&graph, &sr, msg);
	assert(info == GrB_SUCCESS);

	assert(sr.n_scc == 2);
	GrB_Index c0 = sr.scc_ids[0];
	GrB_Index c1 = sr.scc_ids[1];
	assert(c0 != c1);

	assert(scc_reach_has(sr.scc_reach, c0, c0));
	assert(scc_reach_has(sr.scc_reach, c1, c1));
	assert(scc_reach_has(sr.scc_reach, c0, c1));
	assert(!scc_reach_has(sr.scc_reach, c1, c0));
	assert(matrix_nvals(sr.scc_reach) == 3);

	scc_result_free(&sr);
	free_simple_graph_arrays(&graph);
	GrB_Matrix_free(&op);
	GrB_Matrix_free(&cp);
}

// =================================================================================
// Tests for remove_not_path
// =================================================================================

// Empty graph + empty over_approx -> empty result
static void test_removeNotPath_empty_graph_empty_approx(void) {
	GrB_Matrix op = make_empty_matrix(0);
	GrB_Matrix cp = make_empty_matrix(0);
	MRGraph graph = make_simple_graph(0, op, cp);

	GrB_Matrix over_approx = make_empty_matrix(0);

	MRGraph result = {0};
	GrB_Info info = remove_not_path(&graph, over_approx, &result, msg);
	assert(info == GrB_SUCCESS);

	assert(mrgraph_nvals(&result) == 0);

	mr_graph_free(&result);
	GrB_Matrix_free(&op);
	GrB_Matrix_free(&cp);
	GrB_Matrix_free(&over_approx);
}

// Graph: A->B (open), B->C (close)
// over_approx = {(A,B)}
static void test_removeNotPath_removes_edge_not_in_approx(void) {
	GrB_Matrix op = make_empty_matrix(3);
	GrB_Matrix cp = make_empty_matrix(3);
	GrB_Matrix_setElement_BOOL(op, true, 0, 1); // A->B
	GrB_Matrix_setElement_BOOL(cp, true, 1, 2); // B->C
	MRGraph graph = make_simple_graph(3, op, cp);

	GrB_Index srcs[] = {0};
	GrB_Index dsts[] = {1};
	GrB_Matrix over_approx = make_over_approx(3, srcs, dsts, 1);

	MRGraph result = {0};
	GrB_Info info = remove_not_path(&graph, over_approx, &result, msg);
	assert(info == GrB_SUCCESS);

	assert(mrgraph_nvals(&result) == 1);
	assert(matrix_has_edge(result.open_par[0], 0, 1));
	assert(!matrix_has_edge(result.close_par[0], 1, 2));

	mr_graph_free(&result);
	free_simple_graph_arrays(&graph);
	GrB_Matrix_free(&op);
	GrB_Matrix_free(&cp);
	GrB_Matrix_free(&over_approx);
}

// Cycle A->B->C->A, all vertices in one SCC
// over_approx = {(A,C)}
static void test_removeNotPath_keeps_edge_via_indirect_scc(void) {
	GrB_Matrix op = make_empty_matrix(3);
	GrB_Matrix cp = make_empty_matrix(3);
	GrB_Matrix_setElement_BOOL(op, true, 0, 1); // A->B
	GrB_Matrix_setElement_BOOL(op, true, 1, 2); // B->C
	GrB_Matrix_setElement_BOOL(cp, true, 2, 0); // C->A
	MRGraph graph = make_simple_graph(3, op, cp);

	GrB_Index srcs[] = {0};
	GrB_Index dsts[] = {2};
	GrB_Matrix over_approx = make_over_approx(3, srcs, dsts, 1);

	MRGraph result = {0};
	GrB_Info info = remove_not_path(&graph, over_approx, &result, msg);
	assert(info == GrB_SUCCESS);

	assert(mrgraph_nvals(&result) == 3);
	assert(matrix_has_edge(result.open_par[0], 0, 1));
	assert(matrix_has_edge(result.open_par[0], 1, 2));
	assert(matrix_has_edge(result.close_par[0], 2, 0));

	mr_graph_free(&result);
	free_simple_graph_arrays(&graph);
	GrB_Matrix_free(&op);
	GrB_Matrix_free(&cp);
	GrB_Matrix_free(&over_approx);
}

// Two SCCs: {A,B} with cycle A->B->A, {C,D} with cycle C->D->C, inter-component B->C
// over_approx = {(A,D)}
static void test_removeNotPath_multiple_sccs_allowed_path(void) {
	GrB_Matrix op = make_empty_matrix(4);
	GrB_Matrix cp = make_empty_matrix(4);
	GrB_Matrix_setElement_BOOL(op, true, 0, 1); // A->B
	GrB_Matrix_setElement_BOOL(op, true, 1, 0); // B->A
	GrB_Matrix_setElement_BOOL(op, true, 2, 3); // C->D
	GrB_Matrix_setElement_BOOL(op, true, 3, 2); // D->C
	GrB_Matrix_setElement_BOOL(cp, true, 1, 2); // B->C
	MRGraph graph = make_simple_graph(4, op, cp);

	GrB_Index srcs[] = {0};
	GrB_Index dsts[] = {3};
	GrB_Matrix over_approx = make_over_approx(4, srcs, dsts, 1);

	MRGraph result = {0};
	GrB_Info info = remove_not_path(&graph, over_approx, &result, msg);
	assert(info == GrB_SUCCESS);

	assert(mrgraph_nvals(&result) == 5);
	assert(matrix_has_edge(result.open_par[0], 0, 1));
	assert(matrix_has_edge(result.open_par[0], 1, 0));
	assert(matrix_has_edge(result.open_par[0], 2, 3));
	assert(matrix_has_edge(result.open_par[0], 3, 2));
	assert(matrix_has_edge(result.close_par[0], 1, 2));

	mr_graph_free(&result);
	free_simple_graph_arrays(&graph);
	GrB_Matrix_free(&op);
	GrB_Matrix_free(&cp);
	GrB_Matrix_free(&over_approx);
}

// Graph: A→B (open), B→C (close), C→C (open, self-loop)
// over_approx = {(A,B)}
static void test_removeNotPath_removes_inter_scc_edge_not_covered(void) {
	GrB_Matrix op = make_empty_matrix(3);
	GrB_Matrix cp = make_empty_matrix(3);
	GrB_Matrix_setElement_BOOL(op, true, 0, 1); // A->B
	GrB_Matrix_setElement_BOOL(op, true, 2, 2); // C->C
	GrB_Matrix_setElement_BOOL(cp, true, 1, 2); // B->C
	MRGraph graph = make_simple_graph(3, op, cp);

	GrB_Index srcs[] = {0};
	GrB_Index dsts[] = {1};
	GrB_Matrix over_approx = make_over_approx(3, srcs, dsts, 1);

	MRGraph result = {0};
	GrB_Info info = remove_not_path(&graph, over_approx, &result, msg);
	assert(info == GrB_SUCCESS);

	assert(mrgraph_nvals(&result) == 1);
	assert(matrix_has_edge(result.open_par[0], 0, 1));
	assert(!matrix_has_edge(result.close_par[0], 1, 2));
	assert(!matrix_has_edge(result.open_par[0], 2, 2));

	mr_graph_free(&result);
	free_simple_graph_arrays(&graph);
	GrB_Matrix_free(&op);
	GrB_Matrix_free(&cp);
	GrB_Matrix_free(&over_approx);
}

int main(void) {
	LAGraph_Init(msg);

	test_computeSccs_single_vertex();
	test_computeSccs_two_disconnected_vertices();
	test_computeSccs_directed_cycle_is_one_scc();
	test_computeSccs_dag_transitive_closure();
	test_computeSccs_parallel_edges();

	test_removeNotPath_empty_graph_empty_approx();
	test_removeNotPath_removes_edge_not_in_approx();
	test_removeNotPath_keeps_edge_via_indirect_scc();
	test_removeNotPath_multiple_sccs_allowed_path();
	test_removeNotPath_removes_inter_scc_edge_not_covered();

	LAGraph_Finalize(msg);
	return 0;
}
