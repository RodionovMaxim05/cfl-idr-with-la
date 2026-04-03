#include "mutual_refinement.h"
#include "LAGraph.h"
#include "approximation.h"
#include "grammar/grammar_analysis_utils.h"
#include "utils/extract_edges.h"
#include "utils/extract_paths.h"
#include "utils_LAGraph.h"
#include <stdlib.h>
#include <string.h>

static GrB_Index count_edges(const MRGraph *graph) {
	GrB_Index total = 0, nvals = 0;
	for (int64_t i = 0; i < graph->n_par; i++) {
		GrB_Matrix_nvals(&nvals, graph->open_par[i]);
		total += nvals;
		GrB_Matrix_nvals(&nvals, graph->close_par[i]);
		total += nvals;
	}
	for (int64_t i = 0; i < graph->n_bra; i++) {
		GrB_Matrix_nvals(&nvals, graph->open_bra[i]);
		total += nvals;
		GrB_Matrix_nvals(&nvals, graph->close_bra[i]);
		total += nvals;
	}
	if (graph->normal != NULL) {
		GrB_Matrix_nvals(&nvals, graph->normal);
		total += nvals;
	}
	return total;
}

static GrB_Info matrix_intersect(GrB_Matrix *result, GrB_Matrix A, GrB_Matrix B,
								 GrB_Index n, char *msg) {
	GrB_Matrix_new(result, GrB_BOOL, n, n);
	return GrB_Matrix_eWiseMult_BinaryOp(*result, NULL, NULL, GrB_LAND, A, B, NULL);
}

GrB_Info build_refined_graph(MRGraph *out, GrB_Matrix *result_matrices,
							 int64_t n_par, int64_t n_bra, bool has_normal,
							 GrB_Index n) {
	out->n = n;
	out->n_par = n_par;
	out->n_bra = n_bra;

	out->open_par = (GrB_Matrix *)malloc(n_par * sizeof(GrB_Matrix));
	out->close_par = (GrB_Matrix *)malloc(n_par * sizeof(GrB_Matrix));
	out->open_bra = (GrB_Matrix *)malloc(n_bra * sizeof(GrB_Matrix));
	out->close_bra = (GrB_Matrix *)malloc(n_bra * sizeof(GrB_Matrix));

	for (int64_t i = 0; i < n_par; i++) {
		out->open_par[i] = result_matrices[2 * i];
		out->close_par[i] = result_matrices[2 * i + 1];
	}
	for (int64_t i = 0; i < n_bra; i++) {
		out->open_bra[i] = result_matrices[2 * n_par + 2 * i];
		out->close_bra[i] = result_matrices[2 * n_par + 2 * i + 1];
	}
	out->normal = has_normal ? result_matrices[2 * n_par + 2 * n_bra] : NULL;

	return GrB_SUCCESS;
}

void free_refined_graph(MRGraph *graph) {
	for (int64_t i = 0; i < graph->n_par; i++) {
		GrB_Matrix_free(&graph->open_par[i]);
		GrB_Matrix_free(&graph->close_par[i]);
	}
	for (int64_t i = 0; i < graph->n_bra; i++) {
		GrB_Matrix_free(&graph->open_bra[i]);
		GrB_Matrix_free(&graph->close_bra[i]);
	}
	if (graph->normal != NULL) {
		GrB_Matrix_free(&graph->normal);
	}
	free((void *)graph->open_par);
	free((void *)graph->close_par);
	free((void *)graph->open_bra);
	free((void *)graph->close_bra);
}

static GrB_Info run_cfl_step(const MRGraph *graph, MRGrammar_t grammar,
							 GrB_Matrix *out_reachability, GrB_Matrix *out_edges,
							 char *msg) {
	GrB_Info info = GrB_SUCCESS;

	GrB_Matrix *adj = assemble_adj_matrices(graph);

	GrB_Type all_paths_t = NULL;
	GrB_Matrix *paths = NULL;
	LAGraph_Calloc((void **)&paths, grammar.nonterms_count, sizeof(GrB_Matrix), msg);

	info = LAGraph_CFL_AllPaths(paths, adj, &all_paths_t, grammar.terms_count,
								grammar.nonterms_count, grammar.rules,
								grammar.rules_count, msg);
	if (info != GrB_SUCCESS) {
		goto cleanup;
	}

	info = extractEdgesFromOutputs(paths, adj, grammar, graph->n, GrB_INDEX_MAX,
								   GrB_INDEX_MAX, out_edges, msg);
	if (info != GrB_SUCCESS) {
		goto cleanup;
	}

	GrB_Matrix_new(out_reachability, GrB_BOOL, graph->n, graph->n);
	extractNonTrivialPaths(paths[0], graph, out_reachability);

cleanup:
	for (int64_t a = 0; a < grammar.nonterms_count; a++) {
		if (paths[a] != NULL) {
			free_AllPaths_matrix(&paths[a]);
		}
	}
	LAGraph_Free((void **)&paths, msg);
	GrB_free(&all_paths_t);
	free((void *)adj);

	return info;
}

GrB_Info mutual_refinement(const MRGraph *graph, MRGrammarType grammar_type,
						   GrB_Matrix *result) {
	GrB_Info info = GrB_SUCCESS;
	char msg[LAGRAPH_MSG_LEN];
	bool has_normal = (graph->normal != NULL);

	GrB_Index initial_edge_count = count_edges(graph);
	if (initial_edge_count == 0) {
		GrB_Matrix_new(result, GrB_BOOL, graph->n, graph->n);
		return GrB_SUCCESS;
	}

	int64_t terms_count = get_terms_count(graph->n_par, graph->n_bra, graph->normal);

	// Alpha phase

	MRGrammar_t alpha_grammar =
		get_alpha_grammar(grammar_type, graph->n_par, graph->n_bra, has_normal);

	GrB_Matrix alpha_reach = NULL;
	GrB_Matrix *alpha_edges = (GrB_Matrix *)malloc(terms_count * sizeof(GrB_Matrix));

	info = run_cfl_step(graph, alpha_grammar, &alpha_reach, alpha_edges, msg);
	grammar_free(&alpha_grammar);
	if (info != GrB_SUCCESS) {
		goto cleanup_alpha;
	}

	MRGraph alpha_graph;
	build_refined_graph(&alpha_graph, alpha_edges, graph->n_par, graph->n_bra,
						has_normal, graph->n);

	// Beta phase

	MRGrammar_t beta_grammar = get_beta_grammar(grammar_type, alpha_graph.n_par,
												alpha_graph.n_bra, has_normal);

	GrB_Matrix beta_reach = NULL;
	GrB_Matrix *beta_edges = (GrB_Matrix *)malloc(terms_count * sizeof(GrB_Matrix));

	info = run_cfl_step(&alpha_graph, beta_grammar, &beta_reach, beta_edges, msg);
	grammar_free(&beta_grammar);
	if (info != GrB_SUCCESS) {
		goto cleanup_beta;
	}

	MRGraph beta_graph;
	build_refined_graph(&beta_graph, beta_edges, alpha_graph.n_par,
						alpha_graph.n_bra, has_normal, graph->n);

	GrB_Index final_edge_count = count_edges(&beta_graph);

	// Check convergence
	if (final_edge_count == 0 || initial_edge_count == final_edge_count) {
		// Stability has been achieved
		info = matrix_intersect(result, alpha_reach, beta_reach, graph->n, msg);
	} else {
		// Recursive refinement
		info = mutual_refinement(&beta_graph, grammar_type, result);
	}

cleanup_beta:
	GrB_Matrix_free(&beta_reach);
	free((void *)beta_edges);
	free_refined_graph(&beta_graph);

cleanup_alpha:
	GrB_Matrix_free(&alpha_reach);
	free((void *)alpha_edges);
	free_refined_graph(&alpha_graph);

	return info;
}
