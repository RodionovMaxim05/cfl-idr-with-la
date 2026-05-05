#include <GraphBLAS.h>
#include <LAGraph.h>
#include <LAGraphX.h>

#include "cfl_idr.h"
#include "grammar/grammar.h"
#include "graph/condensate_graph.h"
#include "graph/remove_not_path.h"
#include "graph/split_into_components.h"
#include "graph/valueflow_extensions.h"
#include "idr_graph.h"
#include "internal/grb_utils.h"
#include "mr_cache.h"
#include "mutual_refinement.h"
#include "utils/extract_edges.h"
#include "utils/extract_paths.h"
#include "valueflow_approx.h"

static GrB_Info process_under_approx_component(IdrGraph *comp, GrB_Index *vmap,
											   GrB_Matrix result, bool valueflow) {
	GrB_Info info = GrB_SUCCESS;
	char msg[LAGRAPH_MSG_LEN];

	GrB_Type all_paths_t = NULL;
	GrB_Matrix *paths = NULL;
	GrB_Matrix *adj_matrices = NULL;
	GrB_Matrix comp_result = NULL;
	GrB_Index *rows = NULL;
	GrB_Index *cols = NULL;
	MRGrammar_t grammar = {0};

	grammar = dyck_grammar(comp->n_par, comp->n_bra, comp->normal != NULL);
	if (!grammar.rules) {
		info = GrB_OUT_OF_MEMORY;
		goto cleanup;
	}

	GRB_TRY(idr_graph_get_adj_matrices(comp, &adj_matrices));

	paths = (GrB_Matrix *)calloc(grammar.nonterms_count, sizeof(GrB_Matrix));
	if (!paths) {
		info = GrB_OUT_OF_MEMORY;
		goto cleanup;
	}

	GRB_TRY(LAGraph_CFL_AllPaths(paths, &all_paths_t, adj_matrices,
								 grammar.terms_count, grammar.nonterms_count,
								 grammar.rules, grammar.rules_count, msg, 0));

	GRB_TRY(GrB_Matrix_new(&comp_result, GrB_BOOL, comp->n, comp->n));
	GRB_TRY(extract_non_trivial_paths(paths[0], &comp_result));

	if (valueflow) {
		GRB_TRY(apply_valueflow_under_approx(paths, adj_matrices, grammar, comp,
											 &comp_result));
	}

	GrB_Index nnz = 0;
	GRB_TRY(GrB_Matrix_nvals(&nnz, comp_result));

	if (nnz > 0) {
		rows = malloc(nnz * sizeof(GrB_Index));
		cols = malloc(nnz * sizeof(GrB_Index));
		if (!rows || !cols) {
			info = GrB_OUT_OF_MEMORY;
			goto cleanup;
		}

		GRB_TRY(GrB_Matrix_extractTuples_BOOL(rows, cols, NULL, &nnz, comp_result));

		for (GrB_Index k = 0; k < nnz; k++) {
			GRB_TRY(GrB_Matrix_setElement_BOOL(result, true, vmap[rows[k]],
											   vmap[cols[k]]));
		}
	}

cleanup:
	free(rows);
	free(cols);
	GrB_Matrix_free(&comp_result);
	LAGraph_CFL_AllPaths_free_outputs(paths, grammar.nonterms_count, &all_paths_t);
	GrB_free(&all_paths_t);
	grammar_free(&grammar);
	free((void *)adj_matrices);
	return info;
}

GrB_Info idr_get_under_approx(const IdrGraph *graph, bool valueflow,
							  GrB_Matrix *result) {
	GrB_Info info = GrB_SUCCESS;

	IdrGraph *components = NULL;
	GrB_Index **vertex_maps = NULL;
	GrB_Index comp_count = 0;
	GrB_Index *rows = NULL;
	GrB_Index *cols = NULL;

	GRB_TRY(GrB_Matrix_new(result, GrB_BOOL, graph->n, graph->n));
	GRB_TRY(split_IdrGraph_into_components(graph, &components, &vertex_maps,
										   &comp_count));

	for (GrB_Index c = 0; c < comp_count; c++) {
		info = process_under_approx_component(&components[c], vertex_maps[c],
											  *result, valueflow);
		idr_graph_free(&components[c]);
		free(vertex_maps[c]);
		if (info != GrB_SUCCESS) {
			break;
		}
	}

cleanup:
	free(rows);
	free(cols);
	free(components);
	free((void *)vertex_maps);
	return info;
}

GrB_Info idr_get_over_approx(const IdrGraph *graph, IdrGrammarType grammar_type,
							 GrB_Matrix under_approx, GrB_Matrix *result,
							 bool valueflow, bool filter_empty) {
	GrB_Info info = GrB_SUCCESS;

	MRCache cache = {0};
	mr_cache_init(&cache);

	CondensationResult cr = {0};
	GrB_Matrix mr_result = NULL;

	if (under_approx == NULL) {
		GRB_TRY(mutual_refinement(graph, grammar_type, result, valueflow,
								  filter_empty, &cache));
		goto cleanup;
	}

	// Collapse mutually reachable vertices
	GRB_TRY(condensate_from_under_approx(graph, under_approx, &cr));

	GRB_TRY(mutual_refinement(&cr.condensed_graph, grammar_type, &mr_result,
							  valueflow, filter_empty, &cache));

	GRB_TRY(expand_result(mr_result, cr.components, graph->n, result));

cleanup:
	GrB_Matrix_free(&mr_result);
	condensation_result_free(&cr);
	mr_cache_free(&cache);
	return info;
}
