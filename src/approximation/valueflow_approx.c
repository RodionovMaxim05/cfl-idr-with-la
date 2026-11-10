#include "valueflow_approx.h"

#include "graph/valueflow_extensions.h"
#include "idr_graph.h"
#include "internal/grb_utils.h"
#include "utils/extract_edges.h"

GrB_Info apply_valueflow_under_approx(GrB_Matrix *comp_result, GrB_Matrix *edges,
									  MRGrammar_t grammar, const IdrGraph *comp) {
	GrB_Info info = GrB_SUCCESS;
	IdrGraph updated_graph = {0};
	GrB_Matrix filtered = NULL;

	GRB_TRY(build_idr_graph(&updated_graph, &edges[grammar.nonterms_count],
							comp->n_par, comp->n_bra, comp->normal != NULL, comp->n,
							false));

	GRB_TRY(filter_bracket_paths(&filtered, &updated_graph, *comp_result));

	GrB_Matrix_free(comp_result);
	*comp_result = filtered;
	filtered = NULL;

cleanup:
	idr_graph_free(&updated_graph);
	return info;
}

GrB_Info apply_valueflow_over_approx(IdrGraph *out_graph, const IdrGraph *graph,
									 GrB_Matrix *reach) {
	GrB_Info info = GrB_SUCCESS;
	GrB_Matrix filtered_paths = NULL;

	GRB_TRY(filter_bracket_paths(&filtered_paths, graph, *reach));

	GrB_Matrix_free(reach);
	*reach = filtered_paths;
	filtered_paths = NULL;

	GRB_TRY(idr_remove_valueflow_unreachable(out_graph, graph));

cleanup:
	return info;
}
