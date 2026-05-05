#include "condensate_graph.h"

#include <LAGraph.h>
#include <LAGraphX.h>

#include "approximation/idr_graph.h"
#include "approximation/mutual_refinement.h"
#include "internal/grb_utils.h"

GrB_Info condensate_from_under_approx(const IdrGraph *graph, GrB_Matrix under_approx,
									  CondensationResult *out) {
	GrB_Info info = GrB_SUCCESS;
	char msg[LAGRAPH_MSG_LEN];

	GrB_Vector components = NULL;
	GrB_Matrix under_T = NULL;
	GrB_Matrix mutual = NULL;
	LAGraph_Graph G = NULL;
	GrB_Matrix P = NULL;
	GrB_Matrix P_T = NULL;
	GrB_Matrix *condensed_matrices = NULL;
	GrB_Matrix *adj = NULL;

	GRB_TRY(GrB_Matrix_new(&under_T, GrB_BOOL, graph->n, graph->n));
	GRB_TRY(GrB_Matrix_new(&mutual, GrB_BOOL, graph->n, graph->n));
	GRB_TRY(GrB_transpose(under_T, NULL, NULL, under_approx, NULL));
	GRB_TRY(GrB_Matrix_eWiseMult_BinaryOp(mutual, NULL, NULL, GrB_LAND, under_approx,
										  under_T, NULL));
	GrB_Matrix_free(&under_T);

	GRB_TRY(LAGraph_New(&G, &mutual, LAGraph_ADJACENCY_UNDIRECTED, msg));
	G->is_symmetric_structure = LAGraph_TRUE;
	GRB_TRY(LAGr_ConnectedComponents(&components, G, msg));
	GRB_TRY(LAGraph_Delete(&G, msg));

	// Construct a permutation matrix P and P_T
	GRB_TRY(GrB_Matrix_new(&P, GrB_BOOL, graph->n, graph->n));
	GRB_TRY(GrB_Matrix_new(&P_T, GrB_BOOL, graph->n, graph->n));

	for (GrB_Index i = 0; i < graph->n; i++) {
		uint64_t rep;
		GRB_TRY(GrB_Vector_extractElement_UINT64(&rep, components, i));
		GRB_TRY(GrB_Matrix_setElement_BOOL(P, true, rep, i));
	}
	GRB_TRY(GrB_transpose(P_T, NULL, NULL, P, NULL));

	// Condense each matrix: condensed = P * adj * P_T
	int64_t terms_count = get_terms_count(graph->n_par, graph->n_bra, graph->normal);
	condensed_matrices = (GrB_Matrix *)malloc(terms_count * sizeof(GrB_Matrix));
	if (!condensed_matrices) {
		info = GrB_OUT_OF_MEMORY;
		goto cleanup;
	}

	GRB_TRY(idr_graph_get_adj_matrices(graph, &adj));

	for (int64_t t = 0; t < terms_count; t++) {
		GrB_Matrix tmp = NULL;
		GRB_TRY(
			GrB_Matrix_new(&condensed_matrices[t], GrB_BOOL, graph->n, graph->n));
		GRB_TRY(GrB_Matrix_new(&tmp, GrB_BOOL, graph->n, graph->n));

		GrB_Info tmp_info =
			GrB_mxm(tmp, NULL, NULL, GrB_LOR_LAND_SEMIRING_BOOL, P, adj[t], NULL);
		if (tmp_info >= GrB_SUCCESS) {
			tmp_info = GrB_mxm(condensed_matrices[t], NULL, NULL,
							   GrB_LOR_LAND_SEMIRING_BOOL, tmp, P_T, NULL);
		}
		GrB_Matrix_free(&tmp);
		if (tmp_info < GrB_SUCCESS) {
			info = tmp_info;
			goto cleanup;
		}
	}

	free((void *)adj);
	adj = NULL;

	build_idr_graph(&out->condensed_graph, condensed_matrices, graph->n_par,
					graph->n_bra, graph->normal != NULL, graph->n, true);
	out->components = components;
	components = NULL;
	free(condensed_matrices);
	condensed_matrices = NULL;

cleanup:
	GrB_Matrix_free(&under_T);
	GrB_Matrix_free(&mutual);
	if (G) {
		LAGraph_Delete(&G, msg);
	}
	GrB_Vector_free(&components);
	GrB_Matrix_free(&P);
	GrB_Matrix_free(&P_T);
	if (condensed_matrices) {
		for (int64_t t = 0; t < terms_count; t++) {
			GrB_Matrix_free(&condensed_matrices[t]);
		}
		free((void *)condensed_matrices);
	}
	free((void *)adj);
	return info;
}

void condensation_result_free(CondensationResult *cr) {
	idr_graph_free(&cr->condensed_graph);
	GrB_Vector_free(&cr->components);
}

GrB_Info expand_result(GrB_Matrix mr_result, GrB_Vector components, GrB_Index n,
					   GrB_Matrix *result) {
	GrB_Info info = GrB_SUCCESS;

	GrB_Matrix P = NULL;
	GrB_Matrix P_T = NULL;
	GrB_Matrix tmp = NULL;
	GrB_Matrix cluster = NULL;
	GrB_Matrix diag = NULL;

	// Construct a permutation matrix P: P[rep(i)][i] = true
	GRB_TRY(GrB_Matrix_new(&P, GrB_BOOL, n, n));
	for (GrB_Index i = 0; i < n; i++) {
		uint64_t rep;
		GrB_Info vi = GrB_Vector_extractElement_UINT64(&rep, components, i);
		if (vi == GrB_SUCCESS) {
			GRB_TRY(GrB_Matrix_setElement_BOOL(P, true, rep, i));
		}
	}

	// P_T: transposed P[i][rep(i)] = true
	GRB_TRY(GrB_Matrix_new(&P_T, GrB_BOOL, n, n));
	GRB_TRY(GrB_transpose(P_T, NULL, NULL, P, NULL));

	// result = P_T * mr_result * P
	GRB_TRY(GrB_Matrix_new(&tmp, GrB_BOOL, n, n));
	GRB_TRY(
		GrB_mxm(tmp, NULL, NULL, GrB_LOR_LAND_SEMIRING_BOOL, P_T, mr_result, NULL));

	GRB_TRY(GrB_Matrix_new(result, GrB_BOOL, n, n));
	GRB_TRY(GrB_mxm(*result, NULL, NULL, GrB_LOR_LAND_SEMIRING_BOOL, tmp, P, NULL));

	// Add intra-cluster connections
	// To do this, construct a cluster matrix:
	// cluster[i][j] = true if rep(i) == rep(j) cluster = P_T * P
	GRB_TRY(GrB_Matrix_new(&cluster, GrB_BOOL, n, n));
	GRB_TRY(GrB_mxm(cluster, NULL, NULL, GrB_LOR_LAND_SEMIRING_BOOL, P_T, P, NULL));

	// Remove the diagonal from cluster (i != j)
	GRB_TRY(GrB_Matrix_new(&diag, GrB_BOOL, n, n));
	for (GrB_Index i = 0; i < n; i++) {
		GRB_TRY(GrB_Matrix_setElement_BOOL(diag, true, i, i));
	}
	// cluster = cluster - diag (remove i==j)
	GRB_TRY(GrB_Matrix_eWiseAdd_BinaryOp(cluster, diag, NULL, GrB_LOR, cluster,
										 cluster, GrB_DESC_RC));

	// result = result OR cluster
	GRB_TRY(GrB_Matrix_eWiseAdd_BinaryOp(*result, NULL, NULL, GrB_LOR, *result,
										 cluster, NULL));

cleanup:
	GrB_Matrix_free(&P);
	GrB_Matrix_free(&P_T);
	GrB_Matrix_free(&tmp);
	GrB_Matrix_free(&cluster);
	GrB_Matrix_free(&diag);

	return info;
}
