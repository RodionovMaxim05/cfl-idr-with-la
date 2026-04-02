#include "approximation.h"
#include "LAGraph.h"
#include "LAGraphX.h"
#include "grammar/grammar.h"
#include "graph/condensate_graph.h"
#include "utils/extract_edges.h"
#include "utils/extract_paths.h"
#include "utils_LAGraph.h"

GrB_Matrix *assemble_adj_matrices(const MRGraph *graph) {
	int64_t terms_count = get_terms_count(graph->n_par, graph->n_bra, graph->normal);
	GrB_Matrix *adj = (GrB_Matrix *)malloc(terms_count * sizeof(GrB_Matrix));

	for (int64_t i = 0; i < graph->n_par; i++) {
		adj[2 * i] = graph->open_par[i];
		adj[2 * i + 1] = graph->close_par[i];
	}
	for (int64_t i = 0; i < graph->n_bra; i++) {
		adj[2 * graph->n_par + 2 * i] = graph->open_bra[i];
		adj[2 * graph->n_par + 2 * i + 1] = graph->close_bra[i];
	}
	if (graph->normal != NULL) {
		adj[2 * graph->n_par + 2 * graph->n_bra] = graph->normal;
	}

	return adj;
}

GrB_Info get_under_approx(const MRGraph *graph, GrB_Matrix *result) {
	char msg[LAGRAPH_MSG_LEN];
	GrB_Info info = GrB_SUCCESS;

	MRGrammar_t grammar =
		dyck_grammar(graph->n_par, graph->n_bra, graph->normal != NULL);
	GrB_Matrix *adj_matrices = assemble_adj_matrices(graph);

	GrB_Type all_paths_t = NULL;
	GrB_Matrix *paths = NULL;
	LAGraph_Calloc((void **)&paths, grammar.nonterms_count, sizeof(GrB_Matrix), msg);

	info = LAGraph_CFL_AllPaths(paths, adj_matrices, &all_paths_t,
								grammar.terms_count, grammar.nonterms_count,
								grammar.rules, grammar.rules_count, msg);
	if (info != GrB_SUCCESS) {
		goto cleanup;
	}

	GrB_Matrix_new(result, GrB_BOOL, graph->n, graph->n);
	extractNonTrivialPaths(paths[0], graph, result);

cleanup:
	for (int64_t a = 0; a < grammar.nonterms_count; a++) {
		if (paths[a] != NULL) {
			free_AllPaths_matrix(&paths[a]);
		}
	}
	LAGraph_Free((void **)&paths, msg);

	GrB_free(&all_paths_t);
	grammar_free(&grammar);
	free((void *)adj_matrices);

	return info;
}

GrB_Info get_over_approx(const MRGraph *graph, MRGrammarType grammar_type,
						 GrB_Matrix under_approx, GrB_Matrix *result) {
	GrB_Info info = GrB_SUCCESS;
	char msg[LAGRAPH_MSG_LEN];

	if (under_approx == NULL) {
		return mutual_refinement(graph, grammar_type, result);
	}

	CondensationResult cr = {0};

	// Collapse mutually reachable vertices
	info = condensate_from_under_approx(graph, under_approx, &cr, msg);
	if (info != GrB_SUCCESS) {
		return info;
	}

	GrB_Matrix mr_result = NULL;
	info = mutual_refinement(&cr.condensed_graph, grammar_type, &mr_result);
	if (info != GrB_SUCCESS) {
		goto cleanup;
	}

	info = expand_result(mr_result, cr.components, graph->n, result, msg);

cleanup:
	GrB_Matrix_free(&mr_result);
	condensation_result_free(&cr, msg);
	return info;
}
