#include "mutual_refinement.h"

#include <LAGraph.h>
#include <LAGraphX.h>
#include <stdlib.h>
#include <string.h>

#include "grammar/grammar_analysis_utils.h"
#include "graph/remove_not_path.h"
#include "graph/split_into_components.h"
#include "idr_graph.h"
#include "internal/grb_utils.h"
#include "utils/extract_paths.h"
#include "valueflow_approx.h"

/**
 * @brief Checks whether a CFL-reachability result matrix contains a specific
 * source-to-target path.
 *
 * @param[in] m            Matrix to query. Must be valid and of compatible type.
 * @param[in] target_path  Target path specification (source and target indices).
 *
 * @return `true` if the path exists in the matrix, `false` otherwise.
 */
static bool matrix_has_path_udt(GrB_Matrix m, const TargetPath *target_path) {
	AllPathsElem val = {0};
	GrB_Info info =
		GrB_Matrix_extractElement_UDT(&val, m, target_path->src, target_path->tgt);
	return info == GrB_SUCCESS;
}

/**
 * @brief Computes the element-wise logical intersection of three boolean matrices.
 *
 * Performs `result = A ∧ B ∧ C` using GraphBLAS `GrB_Matrix_eWiseMult_BinaryOp`
 * with the `GrB_LAND` operator. If `C` is `NULL`, computes only `A ∧ B`.
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

	if (C == NULL) {
		return GrB_Matrix_eWiseMult_BinaryOp(*result, NULL, NULL, GrB_LAND, A, B,
											 NULL);
	}

	GrB_Matrix temp;
	info = GrB_Matrix_new(&temp, GrB_BOOL, n, n);
	if (info != GrB_SUCCESS) {
		goto cleanup;
	}

	info = GrB_Matrix_eWiseMult_BinaryOp(temp, NULL, NULL, GrB_LAND, A, B, NULL);
	if (info != GrB_SUCCESS) {
		goto cleanup;
	}

	info =
		GrB_Matrix_eWiseMult_BinaryOp(*result, NULL, NULL, GrB_LAND, temp, C, NULL);
	if (info != GrB_SUCCESS) {
		goto cleanup;
	}

	GrB_Matrix_free(&temp);
	return info;

cleanup:
	GrB_Matrix_free(result);
	GrB_Matrix_free(&temp);
	*result = NULL;
	return info;
}

/**
 * @brief Executes a single CFL-reachability computation step with caching support.
 *
 * This function orchestrates one phase of the mutual-refinement algorithm:
 * 1. Computes a cache key via `get_graph_cache_hash`.
 * 2. Checks the cache for a precomputed result matching `(graph_key, grammar_tag)`.
 * 3. If cached: extracts edges and reachability from stored results.
 * 4. If not cached: invokes `LAGraph_CFL_AllPaths`, extracts results, and inserts
 *    them into the cache via `mr_cache_insert`.
 * 5. Optionally checks for the presence of a target path.
 *
 * @param[in]  graph             Input graph for CFL analysis.
 * @param[in]  grammar           Grammar specification (rules, nonterminals, etc.).
 * @param[in]  grammar_tag       Tag identifying the grammar variant for caching.
 * @param[in]  target_path       Optional target path for early termination checks.
 * @param[in]  cache             Cache for storing/retrieving intermediate results.
 * @param[out] out_reachability  Output matrix: reachability via start symbol.
 * @param[out] out_edges         Array of edge matrices extracted for each terminal.
 * @param[out] target_found      Output flag: `true` if target path was found.
 *
 * @return `GrB_SUCCESS` on success, or a GraphBLAS/LAGraph error code on failure.
 *
 * @note `out_edges` must point to a pre-allocated array of sufficient size
 *       (at least `grammar.terms_count` elements).
 */
static GrB_Info run_cfl_step(const IdrGraph *graph, MRGrammar_t grammar,
							 uint32_t grammar_tag, const TargetPath *target_path,
							 MRCache *cache, GrB_Matrix *out_reachability,
							 GrB_Matrix **out_edges, bool *target_found) {
	GrB_Info info = GrB_SUCCESS;
	char msg[LAGRAPH_MSG_LEN];
	GrB_Matrix *adj = NULL;
	GrB_Type all_paths_t = NULL;
	GrB_Matrix *paths = NULL;

	uint64_t key = get_graph_cache_hash(graph);
	const MRStepResult *hit = mr_cache_lookup(cache, key, grammar_tag);

	GRB_TRY(idr_graph_collect_matrices(&adj, graph, grammar.nonterms_count,
									   &grammar.config, grammar.k));

	if (hit) {
		GRB_TRY(extract_edges_from_outputs(*out_edges, hit->matrices, adj, grammar,
										   graph->n, graph->n_par, graph->n_bra,
										   &grammar.config, target_path));

		GRB_TRY(GrB_Matrix_new(out_reachability, GrB_BOOL, graph->n, graph->n));

		if (target_path) {
			*target_found = matrix_has_path_udt(hit->matrices[0], target_path);
			if (*target_found) {
				GRB_TRY(GrB_Matrix_setElement_BOOL(
					*out_reachability, true, target_path->src, target_path->tgt));
			} else {
				GrB_Matrix_free(out_reachability);
			}
		} else {
			GRB_TRY(extract_non_trivial_paths(out_reachability, hit->matrices[0]));
		}

		goto cleanup;
	}

	paths = (GrB_Matrix *)calloc(grammar.nonterms_count + grammar.terms_count,
								 sizeof(GrB_Matrix));
	if (!paths) {
		info = GrB_OUT_OF_MEMORY;
		goto cleanup;
	}

	GRB_TRY(LAGraph_CFL_AllPaths_adv(
		paths, &all_paths_t, adj, grammar.terms_count + grammar.nonterms_count,
		grammar.rules, grammar.rules_count, msg, OPT_EMPTY | OPT_BLOCK));

	GRB_TRY(extract_edges_from_outputs(*out_edges, paths, adj, grammar, graph->n,
									   graph->n_par, graph->n_bra, &grammar.config,
									   target_path));

	mr_cache_insert(cache, key, grammar_tag, paths,
					grammar.nonterms_count + grammar.terms_count, all_paths_t);

	GRB_TRY(GrB_Matrix_new(out_reachability, GrB_BOOL, graph->n, graph->n));

	if (target_path) {
		*target_found = matrix_has_path_udt(paths[0], target_path);
		if (*target_found) {
			GRB_TRY(GrB_Matrix_setElement_BOOL(*out_reachability, true,
											   target_path->src, target_path->tgt));
		} else {
			GrB_Matrix_free(out_reachability);
		}
	} else {
		GRB_TRY(extract_non_trivial_paths(out_reachability, paths[0]));
	}

cleanup:
	if (adj) {
		for (int i = 0; i < grammar.nonterms_count; i++) {
			GrB_free(&adj[i]);
		}
		free((void *)adj);
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

	if (target_path != NULL && !target_found) {
		RETURN_EMPTY;
	}

	GRB_TRY(build_idr_graph(&alpha_graph, alpha_edges, graph->n_par, graph->n_bra,
							has_normal, graph->n, filter_empty));
	free((void *)alpha_edges);
	alpha_edges = NULL;

	// Beta phase

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

	if (target_path != NULL && !target_found) {
		RETURN_EMPTY;
	}

	GRB_TRY(build_idr_graph(&beta_graph, beta_edges, alpha_graph.n_par,
							alpha_graph.n_bra, has_normal, alpha_graph.n,
							filter_empty));
	free((void *)beta_edges);
	beta_edges = NULL;

	// Optional: Project phase

	if (grammar_type == IDR_PROJECT || grammar_type == IDR_ALL) {
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

		if (target_path != NULL && !target_found) {
			RETURN_EMPTY;
		}

		GRB_TRY(build_idr_graph(&project_graph, project_edges, beta_graph.n_par,
								beta_graph.n_bra, has_normal, beta_graph.n,
								filter_empty));

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
			MRGrammar_t exclude_grammar = get_exclude_grammar(
				cur_graph->n_par, cur_graph->n_bra, has_normal, ex_bra);
			if (!exclude_grammar.rules) {
				info = GrB_OUT_OF_MEMORY;
				goto cleanup;
			}

			cur_terms_count = get_terms_count(cur_graph->n_par, cur_graph->n_bra,
											  cur_graph->normal);
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
			GRB_TRY(build_idr_graph(&exclude_graph, exclude_edges, cur_graph->n_par,
									cur_graph->n_bra, has_normal, cur_graph->n,
									filter_empty));

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
			info = matrix_intersect3(result, alpha_reach, beta_reach, project_reach,
									 graph->n);
		}
	} else {
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

	GrB_Index *rows = NULL;
	GrB_Index *cols = NULL;

	GRB_TRY(GrB_Matrix_new(result, GrB_BOOL, global_n, global_n));

	for (GrB_Index c = 0; c < comp_count; c++) {
		IdrGraph *comp = &components[c];
		GrB_Index *vmap = vertex_maps[c];
		GrB_Matrix comp_result = NULL;

		if (target_path != NULL) {
			// Search for local src/tgt indices in this component

			GrB_Index local_src = GrB_INVALID_VALUE;
			GrB_Index local_tgt = GrB_INVALID_VALUE;
			for (GrB_Index i = 0; i < comp->n; i++) {
				if (vmap[i] == target_path->src) {
					local_src = i;
				}
				if (vmap[i] == target_path->tgt) {
					local_tgt = i;
				}
			}
			if (local_src == (GrB_Index)GrB_INVALID_VALUE ||
				local_tgt == (GrB_Index)GrB_INVALID_VALUE) {
				continue;
			}

			GrB_Matrix single_path = NULL;
			GRB_TRY(GrB_Matrix_new(&single_path, GrB_BOOL, comp->n, comp->n));
			GRB_TRY(
				GrB_Matrix_setElement_BOOL(single_path, true, local_src, local_tgt));

			IdrGraph reduced_comp = {0};
			bool reduced_owned = false;
			if (!is_all_pairs(single_path, comp->n)) {
				info = remove_not_path(&reduced_comp, comp, single_path);
				GrB_Matrix_free(&single_path);
				if (info != GrB_SUCCESS) {
					goto cleanup;
				}
				reduced_owned = true;
			} else {
				GrB_Matrix_free(&single_path);
				reduced_comp = *comp;
			}

			TargetPath local_target = {.src = local_src, .tgt = local_tgt};
			info = mutual_refinement_single(&comp_result, &reduced_comp,
											grammar_type, valueflow, filter_empty,
											&local_target, cache);
			if (reduced_owned) {
				idr_graph_free(&reduced_comp);
			}
		} else {
			info = mutual_refinement_single(&comp_result, comp, grammar_type,
											valueflow, filter_empty, NULL, cache);
		}

		if (info != GrB_SUCCESS) {
			GrB_Matrix_free(&comp_result);
			goto cleanup;
		}

		GrB_Index nnz = 0;
		GrB_Matrix_nvals(&nnz, comp_result);
		if (nnz > 0) {
			GrB_Index *rows = malloc(nnz * sizeof(GrB_Index));
			GrB_Index *cols = malloc(nnz * sizeof(GrB_Index));
			if (!rows || !cols) {
				GrB_Matrix_free(&comp_result);
				info = GrB_OUT_OF_MEMORY;
				goto cleanup;
			}
			GRB_TRY(
				GrB_Matrix_extractTuples_BOOL(rows, cols, NULL, &nnz, comp_result));

			for (GrB_Index k = 0; k < nnz; k++) {
				GRB_TRY(GrB_Matrix_setElement_BOOL(*result, true, vmap[rows[k]],
												   vmap[cols[k]]));
			}
			free(rows);
			free(cols);
			rows = NULL;
			cols = NULL;
		}
		GrB_Matrix_free(&comp_result);
	}

cleanup:
	free(rows);
	free(cols);
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
