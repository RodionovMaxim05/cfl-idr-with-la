#include <LAGraph.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cfl_idr.h"
#include "graph/split_into_components.h"

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

static void free_components(IdrGraph *components, GrB_Index **vertex_maps,
							GrB_Index count) {
	for (GrB_Index c = 0; c < count; c++) {
		idr_graph_free(&components[c]);
		free(vertex_maps[c]);
	}
	free(components);
	free(vertex_maps);
}

// Finds local index of global vertex v in component c
static GrB_Index find_local(GrB_Index *vmap, GrB_Index n, GrB_Index v) {
	for (GrB_Index i = 0; i < n; i++) {
		if (vmap[i] == v)
			return i;
	}
	return GrB_INDEX_MAX;
}

// Finds component index that contains global vertex v
static GrB_Index find_comp(GrB_Index **vertex_maps, IdrGraph *components,
						   GrB_Index count, GrB_Index v) {
	for (GrB_Index c = 0; c < count; c++) {
		for (GrB_Index i = 0; i < components[c].n; i++) {
			if (vertex_maps[c][i] == v)
				return c;
		}
	}
	return GrB_INDEX_MAX;
}

// Tests

// empty graph returns empty list
static void test_empty_graph(void) {
	GrB_Matrix op = make_empty_matrix(3);
	GrB_Matrix cp = make_empty_matrix(3);
	IdrGraph graph = make_simple_graph(3, op, cp);

	IdrGraph *components = NULL;
	GrB_Index **vertex_maps = NULL;
	GrB_Index count = 0;

	GrB_Info info =
		split_IdrGraph_into_components(&graph, &components, &vertex_maps, &count);
	assert(info == GrB_SUCCESS);
	assert(count == 0);

	free_components(components, vertex_maps, count);
	free_simple_graph_arrays(&graph);
	GrB_Matrix_free(&op);
	GrB_Matrix_free(&cp);
}

// single edge forms one component
static void test_single_edge_one_component(void) {
	GrB_Matrix op = make_empty_matrix(2);
	GrB_Matrix cp = make_empty_matrix(2);
	GrB_Matrix_setElement_BOOL(op, true, 0, 1);

	IdrGraph graph = make_simple_graph(2, op, cp);

	IdrGraph *components = NULL;
	GrB_Index **vertex_maps = NULL;
	GrB_Index count = 0;

	GrB_Info info =
		split_IdrGraph_into_components(&graph, &components, &vertex_maps, &count);
	assert(info == GrB_SUCCESS);
	assert(count == 1);
	assert(components[0].n == 2);
	assert(matrix_nvals(components[0].open_par[0]) == 1);

	free_components(components, vertex_maps, count);
	free_simple_graph_arrays(&graph);
	GrB_Matrix_free(&op);
	GrB_Matrix_free(&cp);
}

// two disconnected edges form two components
static void test_two_disconnected_edges_two_components(void) {
	GrB_Matrix op = make_empty_matrix(4);
	GrB_Matrix cp = make_empty_matrix(4);
	GrB_Matrix_setElement_BOOL(op, true, 0, 1);
	GrB_Matrix_setElement_BOOL(op, true, 2, 3);

	IdrGraph graph = make_simple_graph(4, op, cp);

	IdrGraph *components = NULL;
	GrB_Index **vertex_maps = NULL;
	GrB_Index count = 0;

	GrB_Info info =
		split_IdrGraph_into_components(&graph, &components, &vertex_maps, &count);
	assert(info == GrB_SUCCESS);
	assert(count == 2);
	for (GrB_Index c = 0; c < count; c++)
		assert(components[c].n == 2);

	free_components(components, vertex_maps, count);
	free_simple_graph_arrays(&graph);
	GrB_Matrix_free(&op);
	GrB_Matrix_free(&cp);
}

// directed cycle is one weakly connected component
static void test_directed_cycle_one_component(void) {
	GrB_Matrix op = make_empty_matrix(3);
	GrB_Matrix cp = make_empty_matrix(3);
	GrB_Matrix_setElement_BOOL(op, true, 0, 1);
	GrB_Matrix_setElement_BOOL(op, true, 1, 2);
	GrB_Matrix_setElement_BOOL(op, true, 2, 0);

	IdrGraph graph = make_simple_graph(3, op, cp);

	IdrGraph *components = NULL;
	GrB_Index **vertex_maps = NULL;
	GrB_Index count = 0;

	GrB_Info info =
		split_IdrGraph_into_components(&graph, &components, &vertex_maps, &count);
	assert(info == GrB_SUCCESS);
	assert(count == 1);
	assert(components[0].n == 3);

	free_components(components, vertex_maps, count);
	free_simple_graph_arrays(&graph);
	GrB_Matrix_free(&op);
	GrB_Matrix_free(&cp);
}

// two separate cycles form two components
static void test_two_separate_cycles_two_components(void) {
	GrB_Matrix op = make_empty_matrix(4);
	GrB_Matrix cp = make_empty_matrix(4);
	GrB_Matrix_setElement_BOOL(op, true, 0, 1);
	GrB_Matrix_setElement_BOOL(cp, true, 1, 0);
	GrB_Matrix_setElement_BOOL(op, true, 2, 3);
	GrB_Matrix_setElement_BOOL(cp, true, 3, 2);

	IdrGraph graph = make_simple_graph(4, op, cp);

	IdrGraph *components = NULL;
	GrB_Index **vertex_maps = NULL;
	GrB_Index count = 0;

	GrB_Info info =
		split_IdrGraph_into_components(&graph, &components, &vertex_maps, &count);
	assert(info == GrB_SUCCESS);
	assert(count == 2);
	for (GrB_Index c = 0; c < count; c++)
		assert(components[c].n == 2);

	free_components(components, vertex_maps, count);
	free_simple_graph_arrays(&graph);
	GrB_Matrix_free(&op);
	GrB_Matrix_free(&cp);
}

// vertex with self-loop is zero component (isolated vertex is deleted)
static void test_self_loop_one_component(void) {
	GrB_Matrix op = make_empty_matrix(1);
	GrB_Matrix cp = make_empty_matrix(1);
	GrB_Matrix_setElement_BOOL(op, true, 0, 0);

	IdrGraph graph = make_simple_graph(1, op, cp);

	IdrGraph *components = NULL;
	GrB_Index **vertex_maps = NULL;
	GrB_Index count = 0;

	GrB_Info info =
		split_IdrGraph_into_components(&graph, &components, &vertex_maps, &count);
	assert(info == GrB_SUCCESS);
	assert(count == 0);

	free_simple_graph_arrays(&graph);
	GrB_Matrix_free(&op);
	GrB_Matrix_free(&cp);
}

// isolated vertex plus connected pair -> one components (isolated vertex is deleted)
static void test_isolated_vertex_plus_connected_pair(void) {
	// 0 isolated, 1->2
	GrB_Matrix op = make_empty_matrix(3);
	GrB_Matrix cp = make_empty_matrix(3);
	GrB_Matrix_setElement_BOOL(op, true, 1, 2);

	IdrGraph graph = make_simple_graph(3, op, cp);

	IdrGraph *components = NULL;
	GrB_Index **vertex_maps = NULL;
	GrB_Index count = 0;

	GrB_Info info =
		split_IdrGraph_into_components(&graph, &components, &vertex_maps, &count);
	assert(info == GrB_SUCCESS);
	assert(count == 1);

	GrB_Index sizes[1] = {components[0].n};
	assert(sizes[0] == 2);

	free_components(components, vertex_maps, count);
	free_simple_graph_arrays(&graph);
	GrB_Matrix_free(&op);
	GrB_Matrix_free(&cp);
}

// complex graph with three components
static void test_complex_graph_three_components(void) {
	// Comp 1: 0<->1->2 (3 edges total: op[0,1], cp[1,0], op[1,2])
	// Comp 2: 3->4<-5  (2 edges total: op[3,4], op[5,4])
	GrB_Matrix op = make_empty_matrix(7);
	GrB_Matrix cp = make_empty_matrix(7);

	GrB_Matrix_setElement_BOOL(op, true, 0, 1);
	GrB_Matrix_setElement_BOOL(cp, true, 1, 0);
	GrB_Matrix_setElement_BOOL(op, true, 1, 2);

	GrB_Matrix_setElement_BOOL(op, true, 3, 4);
	GrB_Matrix_setElement_BOOL(op, true, 5, 4);

	IdrGraph graph = make_simple_graph(7, op, cp);

	IdrGraph *components = NULL;
	GrB_Index **vertex_maps = NULL;
	GrB_Index count = 0;

	GrB_Info info =
		split_IdrGraph_into_components(&graph, &components, &vertex_maps, &count);
	assert(info == GrB_SUCCESS);
	assert(count == 2);

	GrB_Index comp_size3_edges3 = GrB_INDEX_MAX;
	GrB_Index comp_size3_edges2 = GrB_INDEX_MAX;

	for (GrB_Index c = 0; c < count; c++) {
		GrB_Index n = components[c].n;
		GrB_Index e = idr_graph_count_edges(&components[c]);
		if (n == 3 && e == 3)
			comp_size3_edges3 = c;
		else if (n == 3 && e == 2)
			comp_size3_edges2 = c;
	}

	assert(comp_size3_edges3 != GrB_INDEX_MAX);
	assert(comp_size3_edges2 != GrB_INDEX_MAX);

	free_components(components, vertex_maps, count);
	free_simple_graph_arrays(&graph);
	GrB_Matrix_free(&op);
	GrB_Matrix_free(&cp);
}

// edges are remapped correctly to local indices
static void test_components_preserve_edges(void) {
	// 0->1->2, single component
	GrB_Matrix op = make_empty_matrix(3);
	GrB_Matrix cp = make_empty_matrix(3);
	GrB_Matrix_setElement_BOOL(op, true, 0, 1);
	GrB_Matrix_setElement_BOOL(op, true, 1, 2);

	IdrGraph graph = make_simple_graph(3, op, cp);

	IdrGraph *components = NULL;
	GrB_Index **vertex_maps = NULL;
	GrB_Index count = 0;

	GrB_Info info =
		split_IdrGraph_into_components(&graph, &components, &vertex_maps, &count);
	assert(info == GrB_SUCCESS);
	assert(count == 1);
	assert(components[0].n == 3);
	assert(matrix_nvals(components[0].open_par[0]) == 2);

	GrB_Index l0 = find_local(vertex_maps[0], components[0].n, 0);
	GrB_Index l1 = find_local(vertex_maps[0], components[0].n, 1);
	GrB_Index l2 = find_local(vertex_maps[0], components[0].n, 2);
	assert(l0 != GrB_INDEX_MAX);
	assert(l1 != GrB_INDEX_MAX);
	assert(l2 != GrB_INDEX_MAX);
	assert(matrix_has_edge(components[0].open_par[0], l0, l1));
	assert(matrix_has_edge(components[0].open_par[0], l1, l2));

	free_components(components, vertex_maps, count);
	free_simple_graph_arrays(&graph);
	GrB_Matrix_free(&op);
	GrB_Matrix_free(&cp);
}

int main(void) {
	LAGraph_Init(msg);

	test_empty_graph();
	test_single_edge_one_component();
	test_two_disconnected_edges_two_components();
	test_directed_cycle_one_component();
	test_two_separate_cycles_two_components();
	test_self_loop_one_component();
	test_isolated_vertex_plus_connected_pair();
	test_complex_graph_three_components();
	test_components_preserve_edges();

	LAGraph_Finalize(msg);
	return 0;
}
