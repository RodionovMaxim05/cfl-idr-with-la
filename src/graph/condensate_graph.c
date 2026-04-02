#include "condensate_graph.h"
#include "LAGraph.h"
#include "LAGraphX.h"
#include "approximation/approximation.h"
#include "approximation/mutual_refinement.h"

GrB_Info condensate_from_under_approx(const MRGraph *graph, GrB_Matrix under_approx,
									  CondensationResult *out, char *msg) {
	GrB_Info info = GrB_SUCCESS;

	GrB_Vector components = NULL;
	GrB_Matrix under_T = NULL;
	GrB_Matrix mutual = NULL;
	GrB_Matrix_new(&under_T, GrB_BOOL, graph->n, graph->n);
	GrB_Matrix_new(&mutual, GrB_BOOL, graph->n, graph->n);

	GrB_transpose(under_T, NULL, NULL, under_approx, NULL);
	GrB_Matrix_eWiseMult_BinaryOp(mutual, NULL, NULL, GrB_LAND, under_approx,
								  under_T, NULL);
	GrB_Matrix_free(&under_T);

	LAGraph_Graph G = NULL;
	LAGraph_New(&G, &mutual, LAGraph_ADJACENCY_UNDIRECTED, msg);
	G->is_symmetric_structure = LAGraph_TRUE;

	info = LAGr_ConnectedComponents(&components, G, msg);
	LAGraph_Delete(&G, msg);
	if (info != GrB_SUCCESS) {
		return info;
	}

	// Construct a permutation matrix P and P_T
	GrB_Matrix P = NULL;
	GrB_Matrix P_T = NULL;
	GrB_Matrix_new(&P, GrB_BOOL, graph->n, graph->n);
	GrB_Matrix_new(&P_T, GrB_BOOL, graph->n, graph->n);

	for (GrB_Index i = 0; i < graph->n; i++) {
		uint64_t rep;
		GrB_Vector_extractElement_UINT64(&rep, components, i);
		GrB_Matrix_setElement_BOOL(P, true, rep, i);
	}
	GrB_transpose(P_T, NULL, NULL, P, NULL);

	// Condense each matrix: condensed = P * adj * P_T
	int64_t terms_count = get_terms_count(graph->n_par, graph->n_bra, graph->normal);
	GrB_Matrix *condensed_matrices =
		(GrB_Matrix *)malloc(terms_count * sizeof(GrB_Matrix));
	GrB_Matrix *adj = assemble_adj_matrices(graph);

	for (int64_t t = 0; t < terms_count; t++) {
		GrB_Matrix tmp = NULL;
		GrB_Matrix_new(&condensed_matrices[t], GrB_BOOL, graph->n, graph->n);
		GrB_Matrix_new(&tmp, GrB_BOOL, graph->n, graph->n);

		GrB_mxm(tmp, NULL, NULL, GrB_LOR_LAND_SEMIRING_BOOL, P, adj[t], NULL);
		GrB_mxm(condensed_matrices[t], NULL, NULL, GrB_LOR_LAND_SEMIRING_BOOL, tmp,
				P_T, NULL);

		GrB_Matrix_free(&tmp);
	}

	free((void *)adj);
	GrB_Matrix_free(&P);
	GrB_Matrix_free(&P_T);

	build_refined_graph(&out->condensed_graph, condensed_matrices, graph->n_par,
						graph->n_bra, graph->normal != NULL, graph->n);
	out->components = components;

	return GrB_SUCCESS;
}

void condensation_result_free(CondensationResult *cr, char *msg) {
	free_refined_graph(&cr->condensed_graph);
	GrB_Vector_free(&cr->components);
}

GrB_Info expand_result(GrB_Matrix mr_result, GrB_Vector components, GrB_Index n,
					   GrB_Matrix *result, char *msg) {
	GrB_Info info = GrB_SUCCESS;

	// Construct a permutation matrix P: P[rep(i)][i] = true
	GrB_Matrix P = NULL;
	GrB_Matrix_new(&P, GrB_BOOL, n, n);
	for (GrB_Index i = 0; i < n; i++) {
		uint64_t rep;
		GrB_Info vi = GrB_Vector_extractElement_UINT64(&rep, components, i);
		if (vi == GrB_SUCCESS) {
			GrB_Matrix_setElement_BOOL(P, true, rep, i);
		}
	}

	// P_T: transposed P[i][rep(i)] = true
	GrB_Matrix P_T = NULL;
	GrB_Matrix_new(&P_T, GrB_BOOL, n, n);
	GrB_transpose(P_T, NULL, NULL, P, NULL);

	// result = P_T * mr_result * P
	GrB_Matrix tmp = NULL;
	GrB_Matrix_new(&tmp, GrB_BOOL, n, n);
	info =
		GrB_mxm(tmp, NULL, NULL, GrB_LOR_LAND_SEMIRING_BOOL, P_T, mr_result, NULL);
	if (info != GrB_SUCCESS) {
		goto cleanup;
	}

	GrB_Matrix_new(result, GrB_BOOL, n, n);
	info = GrB_mxm(*result, NULL, NULL, GrB_LOR_LAND_SEMIRING_BOOL, tmp, P, NULL);
	if (info != GrB_SUCCESS) {
		goto cleanup;
	}

	// Add intra-cluster connections
	// To do this, construct a cluster matrix:
	// cluster[i][j] = true if rep(i) == rep(j) cluster = P_T * P
	GrB_Matrix cluster = NULL;
	GrB_Matrix_new(&cluster, GrB_BOOL, n, n);
	info = GrB_mxm(cluster, NULL, NULL, GrB_LOR_LAND_SEMIRING_BOOL, P_T, P, NULL);
	if (info != GrB_SUCCESS) {
		goto cleanup;
	}

	// Remove the diagonal from cluster (i != j)
	GrB_Matrix diag = NULL;
	GrB_Matrix_new(&diag, GrB_BOOL, n, n);
	for (GrB_Index i = 0; i < n; i++) {
		GrB_Matrix_setElement_BOOL(diag, true, i, i);
	}
	// cluster = cluster - diag (remove i==j)
	GrB_Matrix_eWiseAdd_BinaryOp(cluster, diag, NULL, GrB_LOR, cluster, cluster,
								 GrB_DESC_RC);

	// result = result OR cluster
	GrB_Matrix_eWiseAdd_BinaryOp(*result, NULL, NULL, GrB_LOR, *result, cluster,
								 NULL);

cleanup:
	GrB_Matrix_free(&P);
	GrB_Matrix_free(&P_T);
	GrB_Matrix_free(&tmp);
	GrB_Matrix_free(&cluster);
	GrB_Matrix_free(&diag);

	return info;
}
