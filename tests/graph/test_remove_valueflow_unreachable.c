#include <LAGraph.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "cfl_idr.h"

static char msg[LAGRAPH_MSG_LEN];

// Utils

static IdrGraph init_graph(GrB_Index n) {
	IdrGraph g;
	g.n = n;
	g.n_par = 0;
	g.open_par = NULL;
	g.close_par = NULL;
	g.n_bra = 1;
	g.open_bra = calloc(1, sizeof(GrB_Matrix));
	g.close_bra = calloc(1, sizeof(GrB_Matrix));
	GrB_Matrix_new(&g.open_bra[0], GrB_BOOL, n, n);
	GrB_Matrix_new(&g.close_bra[0], GrB_BOOL, n, n);
	g.normal = NULL;
	return g;
}

// Tests

static void test_no_ob_or_cb_removes_all(void) {
	IdrGraph g = init_graph(2);
	GrB_Matrix_new(&g.normal, GrB_BOOL, 2, 2);
	GrB_Matrix_setElement_BOOL(g.normal, true, 0, 1); // Edge A->B

	IdrGraph out = {0};
	idr_remove_valueflow_unreachable(&out, &g);

	assert(idr_graph_count_edges(&out) == 0);

	idr_graph_free(&g);
	idr_graph_free(&out);
}

static void test_simple_path_preserved(void) {
	// A(0) -ob-> B(1) -normal-> C(2) -cb-> C(2)
	IdrGraph g = init_graph(3);
	GrB_Matrix_new(&g.normal, GrB_BOOL, 3, 3);
	GrB_Matrix_setElement_BOOL(g.open_bra[0], true, 0, 1);
	GrB_Matrix_setElement_BOOL(g.normal, true, 1, 2);
	GrB_Matrix_setElement_BOOL(g.close_bra[0], true, 2, 2);

	IdrGraph out = {0};
	idr_remove_valueflow_unreachable(&out, &g);

	assert(idr_graph_count_edges(&out) == 3);

	idr_graph_free(&g);
	idr_graph_free(&out);
}

static void test_reachable_but_no_cb_removed(void) {
	// A(0) -ob-> B(1) -normal-> C(2)
	IdrGraph g = init_graph(3);
	GrB_Matrix_new(&g.normal, GrB_BOOL, 3, 3);
	GrB_Matrix_setElement_BOOL(g.open_bra[0], true, 0, 1);
	GrB_Matrix_setElement_BOOL(g.normal, true, 1, 2);

	IdrGraph out = {0};
	idr_remove_valueflow_unreachable(&out, &g);

	assert(idr_graph_count_edges(&out) == 0);

	idr_graph_free(&g);
	idr_graph_free(&out);
}

static void test_cycle_preserved(void) {
	// A(0)-ob->B(1), B-normal->C(2), C-normal->B, C-cb->C
	IdrGraph g = init_graph(3);
	GrB_Matrix_new(&g.normal, GrB_BOOL, 3, 3);
	GrB_Matrix_setElement_BOOL(g.open_bra[0], true, 0, 1);
	GrB_Matrix_setElement_BOOL(g.normal, true, 1, 2);
	GrB_Matrix_setElement_BOOL(g.normal, true, 2, 1);
	GrB_Matrix_setElement_BOOL(g.close_bra[0], true, 2, 2);

	IdrGraph out = {0};
	idr_remove_valueflow_unreachable(&out, &g);

	assert(idr_graph_count_edges(&out) == 4);

	idr_graph_free(&g);
	idr_graph_free(&out);
}

static void test_invalid_cb_to_ob_path(void) {
	// A(0)-cb->B(1)-normal->C(2)-ob->D(3)
	IdrGraph g = init_graph(4);
	GrB_Matrix_new(&g.normal, GrB_BOOL, 4, 4);
	GrB_Matrix_setElement_BOOL(g.close_bra[0], true, 0, 1);
	GrB_Matrix_setElement_BOOL(g.normal, true, 1, 2);
	GrB_Matrix_setElement_BOOL(g.open_bra[0], true, 2, 3);

	IdrGraph out = {0};
	idr_remove_valueflow_unreachable(&out, &g);

	assert(idr_graph_count_edges(&out) == 0);

	idr_graph_free(&g);
	idr_graph_free(&out);
}

int main(void) {
	LAGraph_Init(msg);

	test_no_ob_or_cb_removes_all();
	test_simple_path_preserved();
	test_reachable_but_no_cb_removed();
	test_cycle_preserved();
	test_invalid_cb_to_ob_path();

	LAGraph_Finalize(msg);
	return 0;
}
