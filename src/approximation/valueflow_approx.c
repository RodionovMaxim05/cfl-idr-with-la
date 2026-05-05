#include "valueflow_approx.h"

#include "graph/valueflow_extensions.h"
#include "idr_graph.h"
#include "internal/grb_utils.h"
#include "utils/extract_edges.h"

GrB_Info apply_valueflow_under_approx(GrB_Matrix *paths, GrB_Matrix *adj_matrices,
									  MRGrammar_t grammar, const IdrGraph *comp,
									  GrB_Matrix *comp_result) {
	GrB_Info info = GrB_SUCCESS;
	GrB_Matrix *out_edges = NULL;
	IdrGraph updated_graph = {0};
	GrB_Matrix filtered = NULL;

	out_edges = (GrB_Matrix *)malloc(grammar.terms_count * sizeof(GrB_Matrix));
	if (!out_edges) {
		info = GrB_OUT_OF_MEMORY;
		goto cleanup;
	}

	GRB_TRY(extract_edges_from_outputs(paths, adj_matrices, grammar, comp->n, NULL,
									   out_edges));

	GRB_TRY(build_idr_graph(&updated_graph, out_edges, comp->n_par, comp->n_bra,
							comp->normal != NULL, comp->n, false));

	GRB_TRY(GrB_Matrix_new(&filtered, GrB_BOOL, comp->n, comp->n));
	GRB_TRY(filter_bracket_paths(&updated_graph, *comp_result, &filtered));

	GrB_Matrix_free(comp_result);
	*comp_result = filtered;
	filtered = NULL;

cleanup:
	free((void *)out_edges);
	idr_graph_free(&updated_graph);
	GrB_Matrix_free(&filtered);
	return info;
}

GrB_Info apply_valueflow_over_approx(const IdrGraph *graph, GrB_Matrix *beta_reach,
									 IdrGraph *filtered_graph) {
	GrB_Info info = GrB_SUCCESS;
	GrB_Matrix filtered_paths = NULL;

	GRB_TRY(GrB_Matrix_new(&filtered_paths, GrB_BOOL, graph->n, graph->n));
	GRB_TRY(filter_bracket_paths(graph, *beta_reach, &filtered_paths));

	GrB_Matrix_free(beta_reach);
	*beta_reach = filtered_paths;
	filtered_paths = NULL;

	GRB_TRY(idr_remove_valueflow_unreachable(graph, filtered_graph));

cleanup:
	GrB_Matrix_free(&filtered_paths);
	return info;
}
