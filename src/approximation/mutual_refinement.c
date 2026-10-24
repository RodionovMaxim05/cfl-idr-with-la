#include "mutual_refinement.h"

#include <stdlib.h>
#include <string.h>

#include "grammar/grammar_analysis_utils.h"
#include "graph/remove_not_path.h"
#include "graph/split_into_components.h"
#include "idr_graph.h"
#include "internal/grb_utils.h"
#include "mutual_refinement_step.h"
#include "valueflow_approx.h"

/**
 * @brief Computes the element-wise logical intersection of three boolean matrices.
 *
 * Performs `result = A ∧ B ∧ C` using GraphBLAS `GrB_Matrix_eWiseMult_BinaryOp`
 * with the `GrB_LAND` operator. If `C` is `NULL`, computes `A ∧ B`.
 *
 * @param[out] result  Output matrix (newly allocated, n × n, `GrB_BOOL`).
 * @param[in]  A       First input matrix.
 * @param[in]  B       Second input matrix.
 * @param[in]  C       Optional third input matrix; may be `NULL`.
 * @param[in]  n       Dimension of the square matrices.
 *
 * @return `GrB_SUCCESS` on success, or a GraphBLAS error code on failure.
 *         On error, `*result` is set to `NULL` and any temporary allocations
 *         are freed.
 */
static GrB_Info matrix_intersect3(GrB_Matrix *result, GrB_Matrix A, GrB_Matrix B,
								  GrB_Matrix C, GrB_Index n) {
	*result = NULL;

	GrB_Info info = GrB_Matrix_new(result, GrB_BOOL, n, n);
	if (info != GrB_SUCCESS) {
		return info;
	}

	info = GrB_Matrix_eWiseMult_BinaryOp(*result, C, NULL, GrB_LAND, A, B, NULL);
	if (info != GrB_SUCCESS) {
		GrB_Matrix_free(result);
		*result = NULL;
	}

	return info;
}

#define RETURN_EMPTY                                                                \
	do {                                                                            \
		GRB_TRY(GrB_Matrix_new(result, GrB_BOOL, graph->n, graph->n));              \
		goto cleanup;                                                               \
	} while (0)

/**
 * @brief Core mutual-refinement algorithm for a single (non-decomposed) graph.
 *
 * Implements the iterative grammar-based refinement loop:
 * - Executes Alpha, Beta, and optional Project/Exclude phases.
 * - Intersects results from different phases to tighten the approximation.
 * - Checks for convergence (edge count stability) and recurses if needed.
 * - Supports on-demand analysis via `target_path` for early termination.
 * - Applies value-flow post-processing if requested.
 *
 * Uses `run_cfl_step` for each grammar phase and `build_idr_graph` to construct
 * refined graphs from extracted edges.
 *
 * @param[out] result        Output reachability matrix (n × n, `GrB_BOOL`).
 * @param[in]  graph         Input graph to analyze.
 * @param[in]  grammar_type  Grammar variant controlling refinement behavior.
 * @param[in]  valueflow     If `true`, apply value-flow post-processing.
 * @param[in]  filter_empty  If `true`, filter empty parenthesis/bracket pairs.
 * @param[in]  target_path   Optional target path for on-demand analysis.
 * @param[in]  cache         Cache for intermediate CFL results.
 *
 * @return `GrB_SUCCESS` on success, or a GraphBLAS/LAGraph error code on failure.
 */
static GrB_Info mutual_refinement_single(GrB_Matrix *result, const IdrGraph *graph,
										 IdrGrammarType grammar_type, bool valueflow,
										 bool filter_empty,
										 const TargetPath *target_path,
										 MRCache *cache) {
	GrB_Info info = GrB_SUCCESS;
	bool has_normal = (graph->normal != NULL);
	bool target_found = false;
	const GrB_Index n = graph->n;

	GrB_Matrix alpha_reach = NULL;
	GrB_Matrix *alpha_edges = NULL;
	IdrGraph alpha_graph = {0};

	GrB_Matrix beta_reach = NULL;
	GrB_Matrix *beta_edges = NULL;
	IdrGraph beta_graph = {0};

	GrB_Matrix project_reach = NULL;
	GrB_Matrix *project_edges = NULL;
	IdrGraph project_graph = {0};

	GrB_Matrix *exclude_edges = NULL;
	IdrGraph exclude_graph = {0};

	IdrGraph filtered_graph = {0};

	GrB_Index initial_edge_count = idr_graph_count_edges(graph);
	if (initial_edge_count == 0) {
		RETURN_EMPTY;
	}

	int64_t terms_count = get_terms_count(graph->n_par, graph->n_bra, graph->normal);

	int64_t alpha_terms_count = 0;
	int64_t beta_terms_count = 0;
	int64_t cur_terms_count = 0;

	// Alpha phase
	int64_t graph_n_par = graph->n_par;
	int64_t graph_n_bra = graph->n_bra;

	MRGrammar_t alpha_grammar =
		get_alpha_grammar(grammar_type, graph->n_par, graph->n_bra, has_normal);
	if (!alpha_grammar.rules) {
		info = GrB_OUT_OF_MEMORY;
		goto cleanup;
	}

	alpha_edges = (GrB_Matrix *)malloc(terms_count * sizeof(GrB_Matrix));
	if (!alpha_edges) {
		info = GrB_OUT_OF_MEMORY;
		goto cleanup;
	}

	GRB_TRY(run_cfl_step(graph, alpha_grammar, GRAMMAR_TAG_ALPHA, target_path, cache,
						 &alpha_reach, &alpha_edges, &target_found));
	grammar_free(&alpha_grammar);
	idr_graph_free(graph);

	if (target_path != NULL && !target_found) {
		RETURN_EMPTY;
	}

	GRB_TRY(build_idr_graph(&alpha_graph, alpha_edges, graph_n_par, graph_n_bra,
							has_normal, n, filter_empty));
	free((void *)alpha_edges);
	alpha_edges = NULL;

	// Beta phase
	int64_t alpha_n_par = alpha_graph.n_par;
	int64_t alpha_n_bra = alpha_graph.n_bra;

	MRGrammar_t beta_grammar = get_beta_grammar(grammar_type, alpha_graph.n_par,
												alpha_graph.n_bra, has_normal);
	if (!beta_grammar.rules) {
		info = GrB_OUT_OF_MEMORY;
		goto cleanup;
	}

	alpha_terms_count =
		get_terms_count(alpha_graph.n_par, alpha_graph.n_bra, alpha_graph.normal);
	beta_edges = (GrB_Matrix *)malloc(alpha_terms_count * sizeof(GrB_Matrix));
	if (!beta_edges) {
		info = GrB_OUT_OF_MEMORY;
		goto cleanup;
	}

	GRB_TRY(run_cfl_step(&alpha_graph, beta_grammar, GRAMMAR_TAG_BETA, target_path,
						 cache, &beta_reach, &beta_edges, &target_found));
	grammar_free(&beta_grammar);
	idr_graph_free(&alpha_graph);

	if (target_path != NULL && !target_found) {
		RETURN_EMPTY;
	}

	GRB_TRY(build_idr_graph(&beta_graph, beta_edges, alpha_n_par, alpha_n_bra,
							has_normal, n, filter_empty));
	free((void *)beta_edges);
	beta_edges = NULL;

	// Optional: Project phase

	if (grammar_type == IDR_PROJECT || grammar_type == IDR_ALL) {
		int64_t beta_n_par = beta_graph.n_par;
		int64_t beta_n_bra = beta_graph.n_bra;

		MRGrammar_t project_grammar =
			get_project_grammar(beta_graph.n_par, beta_graph.n_bra, has_normal);
		if (!project_grammar.rules) {
			info = GrB_OUT_OF_MEMORY;
			goto cleanup;
		}

		beta_terms_count =
			get_terms_count(beta_graph.n_par, beta_graph.n_bra, beta_graph.normal);
		project_edges = (GrB_Matrix *)malloc(beta_terms_count * sizeof(GrB_Matrix));
		if (!project_edges) {
			info = GrB_OUT_OF_MEMORY;
			goto cleanup;
		}

		GRB_TRY(run_cfl_step(&beta_graph, project_grammar, GRAMMAR_TAG_PROJECT,
							 target_path, cache, &project_reach, &project_edges,
							 &target_found));
		grammar_free(&project_grammar);
		idr_graph_free(&beta_graph);

		if (target_path != NULL && !target_found) {
			RETURN_EMPTY;
		}

		GRB_TRY(build_idr_graph(&project_graph, project_edges, beta_n_par,
								beta_n_bra, has_normal, n, filter_empty));

		free((void *)project_edges);
		project_edges = NULL;
	}

	// Optional: Exclude phase

	const IdrGraph *final_graph =
		(grammar_type == IDR_PROJECT || grammar_type == IDR_ALL) ? &project_graph
																 : &beta_graph;

	if (grammar_type == IDR_EXCLUDE || grammar_type == IDR_ALL) {
		const IdrGraph *cur_graph = final_graph;

		for (int64_t ex_bra = 0; ex_bra < cur_graph->n_bra; ex_bra++) {
			int64_t cur_n_par = cur_graph->n_par;
			int64_t cur_n_bra = cur_graph->n_bra;

			MRGrammar_t exclude_grammar =
				get_exclude_grammar(cur_n_par, cur_n_bra, has_normal, ex_bra);
			if (!exclude_grammar.rules) {
				info = GrB_OUT_OF_MEMORY;
				goto cleanup;
			}

			cur_terms_count =
				get_terms_count(cur_n_par, cur_n_bra, cur_graph->normal);
			exclude_edges =
				(GrB_Matrix *)malloc(cur_terms_count * sizeof(GrB_Matrix));
			if (!exclude_edges) {
				info = GrB_OUT_OF_MEMORY;
				goto cleanup;
			}

			GrB_Matrix cur_reach = NULL;

			GRB_TRY(run_cfl_step(cur_graph, exclude_grammar,
								 GRAMMAR_TAG_EXCLUDE_BASE + (uint32_t)ex_bra,
								 target_path, cache, &cur_reach, &exclude_edges,
								 &target_found));
			grammar_free(&exclude_grammar);
			GrB_Matrix_free(&cur_reach);

			if (target_path != NULL && !target_found) {
				RETURN_EMPTY;
			}

			idr_graph_free(&exclude_graph);
			GRB_TRY(build_idr_graph(&exclude_graph, exclude_edges, cur_n_par,
									cur_n_bra, has_normal, n, filter_empty));

			free((void *)exclude_edges);
			exclude_edges = NULL;
			cur_graph = &exclude_graph;
		}
		final_graph = cur_graph;
	}

	if (valueflow) {
		GRB_TRY(
			apply_valueflow_over_approx(&filtered_graph, final_graph, &beta_reach));
		final_graph = &filtered_graph;
	}

	GrB_Index final_edge_count = idr_graph_count_edges(final_graph);

	// Check convergence
	if (final_edge_count == 0 || initial_edge_count == final_edge_count) {
		// Stability has been achieved

		if (target_path) {
			*result = alpha_reach;
			alpha_reach = NULL;
		} else {
			info =
				matrix_intersect3(result, alpha_reach, beta_reach, project_reach, n);
		}
	} else {
		GrB_Matrix_free(&alpha_reach);
		GrB_Matrix_free(&beta_reach);
		GrB_Matrix_free(&project_reach);

		// Recursive refinement
		GRB_TRY(mutual_refinement(result, final_graph, grammar_type, valueflow,
								  filter_empty, target_path, cache));
	}

cleanup:
	idr_graph_free(&filtered_graph);
	if (exclude_edges) {
		for (int64_t i = 0; i < cur_terms_count; i++) {
			GrB_Matrix_free(&exclude_edges[i]);
		}
		free((void *)exclude_edges);
	}
	idr_graph_free(&exclude_graph);
	GrB_Matrix_free(&project_reach);
	if (project_edges) {
		for (int64_t i = 0; i < beta_terms_count; i++) {
			GrB_Matrix_free(&project_edges[i]);
		}
		free((void *)project_edges);
	}
	idr_graph_free(&project_graph);
	GrB_Matrix_free(&beta_reach);
	if (beta_edges) {
		for (int64_t i = 0; i < alpha_terms_count; i++) {
			GrB_Matrix_free(&beta_edges[i]);
		}
		free((void *)beta_edges);
	}
	idr_graph_free(&beta_graph);
	GrB_Matrix_free(&alpha_reach);
	if (alpha_edges) {
		for (int64_t i = 0; i < terms_count; i++) {
			GrB_Matrix_free(&alpha_edges[i]);
		}
		free((void *)alpha_edges);
	}
	idr_graph_free(&alpha_graph);

	return info;
}

GrB_Info mutual_refinement_with_components(GrB_Matrix *result, IdrGraph *components,
										   GrB_Index **vertex_maps,
										   GrB_Index comp_count, GrB_Index global_n,
										   IdrGrammarType grammar_type,
										   bool valueflow, bool filter_empty,
										   const TargetPath *target_path,
										   MRCache *cache) {
	GrB_Info info = GrB_SUCCESS;

	GrB_Matrix single_path = NULL;
	GrB_Matrix comp_result = NULL;
	IdrGraph reduced_comp = {0};
	bool reduced_owned = false;

	GRB_TRY(GrB_Matrix_new(result, GrB_BOOL, global_n, global_n));

	for (GrB_Index c = 0; c < comp_count; c++) {
		IdrGraph *comp = &components[c];
		GrB_Index *vmap = vertex_maps[c];
		const GrB_Index comp_n = comp->n;

		if (target_path != NULL) {
			// Search for local src/tgt indices in this component

			GrB_Index local_src = GrB_INDEX_MAX;
			GrB_Index local_tgt = GrB_INDEX_MAX;
			for (GrB_Index i = 0; i < comp_n; i++) {
				if (vmap[i] == target_path->src) {
					local_src = i;
				}
				if (vmap[i] == target_path->tgt) {
					local_tgt = i;
				}
				if (local_src != GrB_INDEX_MAX && local_tgt != GrB_INDEX_MAX) {
					break;
				}
			}
			if (local_src == GrB_INDEX_MAX || local_tgt == GrB_INDEX_MAX) {
				continue;
			}

			GRB_TRY(GrB_Matrix_new(&single_path, GrB_BOOL, comp_n, comp_n));
			GRB_TRY(
				GrB_Matrix_setElement_BOOL(single_path, true, local_src, local_tgt));

			const IdrGraph *to_process = comp;
			if (!is_all_pairs(single_path, comp_n)) {
				GRB_TRY(remove_not_path(&reduced_comp, comp, single_path));
				reduced_owned = true;
				to_process = &reduced_comp;
			}
			GrB_Matrix_free(&single_path);

			TargetPath local_target = {.src = local_src, .tgt = local_tgt};
			GRB_TRY(mutual_refinement_single(&comp_result, to_process, grammar_type,
											 valueflow, filter_empty, &local_target,
											 cache));
			if (reduced_owned) {
				idr_graph_free(&reduced_comp);
				reduced_owned = false;
			}
		} else {
			GRB_TRY(mutual_refinement_single(&comp_result, comp, grammar_type,
											 valueflow, filter_empty, NULL, cache));
			idr_graph_free(comp);
		}

		GrB_Index nnz = 0;
		GRB_TRY(GrB_Matrix_nvals(&nnz, comp_result));
		if (nnz > 0) {
			GRB_TRY(GrB_Matrix_assign(*result, NULL, GrB_LOR, comp_result, vmap,
									  comp_n, vmap, comp_n, NULL));
		}
		GrB_Matrix_free(&comp_result);

		if (target_path == NULL) {
			free(vmap);
			vertex_maps[c] = NULL;
		}
	}

cleanup:
	GrB_Matrix_free(&single_path);
	GrB_Matrix_free(&comp_result);
	if (reduced_owned) {
		idr_graph_free(&reduced_comp);
	}
	return info;
}

GrB_Info mutual_refinement(GrB_Matrix *result, const IdrGraph *graph,
						   IdrGrammarType grammar_type, bool valueflow,
						   bool filter_empty, const TargetPath *target_path,
						   MRCache *cache) {
	GrB_Info info = GrB_SUCCESS;

	IdrGraph *components = NULL;
	GrB_Index **vertex_maps = NULL;
	GrB_Index comp_count = 0;

	GRB_TRY(split_IdrGraph_into_components(graph, &components, &vertex_maps,
										   &comp_count));

	GRB_TRY(mutual_refinement_with_components(
		result, components, vertex_maps, comp_count, graph->n, grammar_type,
		valueflow, filter_empty, target_path, cache));

cleanup:
	for (GrB_Index c = 0; c < comp_count; c++) {
		idr_graph_free(&components[c]);
		free(vertex_maps[c]);
	}
	free(components);
	free((void *)vertex_maps);
	return info;
}
