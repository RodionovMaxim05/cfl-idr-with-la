#include "mutual_refinement_step.h"

#include <LAGraphX.h>

#include "idr_graph.h"
#include "internal/grb_utils.h"
#include "internal/system_memory.h"
#include "mutual_refinement.h"
#include "utils/extract_paths.h"

/**
 * @brief Checks whether a CFL-reachability result matrix contains a specific
 * source-to-target path.
 *
 * @param[in] m            Matrix to query. Must be valid and of compatible type.
 * @param[in] target_path  Target path specification (source and target indices).
 *
 * @return `true` if the path exists in the matrix, `false` otherwise.
 */
static bool matrix_has_path_bool(GrB_Matrix m, const TargetPath *target_path) {
	if (!m || !target_path) {
		return false;
	}
	bool val = false;
	GrB_Info info =
		GrB_Matrix_extractElement_BOOL(&val, m, target_path->src, target_path->tgt);
	return (info == GrB_SUCCESS && val);
}

/**
 * @brief Duplicates an array of matrices.
 *
 * Allocates and copies `count` matrices from `src` into `dst`. On failure,
 * all previously allocated matrices in `dst` are safely freed.
 *
 * @param[out] dst    Destination matrix array. Must be preallocated for `count`
 *                    elements.
 * @param[in]  src    Source matrix array.
 * @param[in]  count  Number of matrices to duplicate.
 *
 * @return `GrB_SUCCESS` on success, or GraphBLAS error code on failure.
 */
static GrB_Info dup_matrices_into(GrB_Matrix *dst, GrB_Matrix *src, int64_t count) {
	for (int64_t i = 0; i < count; i++) {
		dst[i] = NULL;
	}
	for (int64_t i = 0; i < count; i++) {
		if (src[i] != NULL) {
			GrB_Info info = GrB_Matrix_dup(&dst[i], src[i]);
			if (info != GrB_SUCCESS) {
				for (int64_t j = 0; j <= i; j++) {
					GrB_free(&dst[j]);
				}
				return info;
			}
		}
	}
	return GrB_SUCCESS;
}

/**
 * @brief Releases per-symbol reachability matrices produced during CFL evaluation.
 *
 * @param[in,out] reach_matrices Array of reachability matrices to free. May be
 *                               `NULL`.
 * @param[in]     count          Number of matrices in `reach_matrices`.
 */
static void free_reach_matrices(GrB_Matrix *reach_matrices, int64_t count) {
	if (!reach_matrices) {
		return;
	}
	for (int64_t i = 0; i < count; i++) {
		GrB_free(&reach_matrices[i]);
	}
	LAGraph_Free((void **)&reach_matrices, NULL);
}

GrB_Info run_cfl_step(const IdrGraph *graph, MRGrammar_t grammar,
					  uint32_t grammar_tag, const TargetPath *target_path,
					  MRCache *cache, GrB_Matrix *out_reachability,
					  GrB_Matrix **out_edges, bool *target_found) {
	GrB_Info info = GrB_SUCCESS;
	char msg[LAGRAPH_MSG_LEN];
	GrB_Matrix *adj_matrices = NULL;
	GrB_Matrix *reach_matrices = NULL;
	GrB_Matrix *edges = NULL;
	int32_t *grammar_to_straight = NULL;
	bool reach_owned = false;
	GrB_Matrix *cached_edges = NULL;
	GrB_Matrix cached_reach = NULL;

	if (out_reachability) {
		*out_reachability = NULL;
	}
	if (out_edges && *out_edges) {
		for (int64_t t = 0; t < grammar.terms_count; t++) {
			(*out_edges)[t] = NULL;
		}
	}
	if (target_found) {
		*target_found = false;
	}

	const uint64_t key = get_graph_cache_hash(graph);

	if (!target_path && cache) {
		const MRExtractedResult *extracted_hit =
			mr_cache_lookup_extracted(cache, key, grammar_tag);

		if (extracted_hit) {
			GRB_TRY(dup_matrices_into(*out_edges, extracted_hit->edges,
									  extracted_hit->edges_count));
			info = GrB_Matrix_dup(out_reachability, extracted_hit->reachability);
			if (info != GrB_SUCCESS) {
				goto cleanup;
			}
			return GrB_SUCCESS;
		}
	}

	GrB_Matrix *reach_res = NULL;
	const MRReachabilityResult *reach_hit =
		(target_path && cache)
			? mr_cache_lookup_reachability(cache, key, grammar_tag)
			: NULL;

	if (reach_hit) {
		reach_res = reach_hit->matrices;
		*target_found = matrix_has_path_bool(reach_res[NT_START], target_path);
		if (!*target_found) {
			return GrB_SUCCESS;
		}
	}

	GRB_TRY(idr_graph_collect_matrices(&adj_matrices, graph, grammar.nonterms_count,
									   &grammar.config, grammar.k));

	size_t symbols_amount = (size_t)(grammar.nonterms_count + grammar.terms_count);

	if (!reach_res) {
		if (cache && (cache->reachability_count + cache->extracted_count > 0)) {
			if (is_system_memory_critical(RAM_BUDGET_PERCENT) == 1) {
				mr_cache_free(cache);
				release_unused_memory_to_os();
				mr_cache_init(cache);
			}
		}

		reach_matrices = (GrB_Matrix *)calloc(symbols_amount, sizeof(GrB_Matrix));
		if (!reach_matrices) {
			info = GrB_OUT_OF_MEMORY;
			goto cleanup;
		}

		GRB_TRY(LAGraph_CFL_reachability_adv(
			reach_matrices, adj_matrices, symbols_amount, grammar.rules,
			grammar.rules_count, msg, OPT_EMPTY | OPT_BLOCK | OPT_PAT));

		reach_owned = true;
		reach_res = reach_matrices;

		if (target_path) {
			*target_found = matrix_has_path_bool(reach_res[NT_START], target_path);

			if (cache && mr_cache_insert_reachability(
							 cache, key, grammar_tag, reach_matrices,
							 (int64_t)symbols_amount) == GrB_SUCCESS) {
				reach_owned = false;
			}

			if (!*target_found) {
				goto cleanup;
			}
		}
	}

	GRB_TRY(GrB_Matrix_new(out_reachability, GrB_BOOL, graph->n, graph->n));

	if (target_path) {
		GRB_TRY(GrB_Matrix_setElement_BOOL(*out_reachability, true, target_path->src,
										   target_path->tgt));
	} else if (reach_res[NT_START] != NULL) {
		GRB_TRY(extract_non_trivial_paths(out_reachability, reach_res[NT_START]));
	}

	GrB_Index reach_nnz = 0;
	GRB_TRY(GrB_Matrix_nvals(&reach_nnz, *out_reachability));

	if (reach_nnz > 0) {
		edges = (GrB_Matrix *)calloc(symbols_amount, sizeof(GrB_Matrix));
		if (!edges) {
			info = GrB_OUT_OF_MEMORY;
			goto cleanup;
		}

		grammar_to_straight =
			build_grammar_to_straight_map(&grammar.config, grammar.terms_count,
										  graph->n_par, graph->n_bra, grammar.k);
		if (!grammar_to_straight) {
			info = GrB_OUT_OF_MEMORY;
			goto cleanup;
		}

		GrB_Matrix pairs = target_path ? *out_reachability : NULL;
		GRB_TRY(LAGraph_CFL_extract_edges(
			edges, reach_res, adj_matrices, symbols_amount, grammar.rules,
			grammar.rules_count, NT_START, pairs, msg));

		for (int64_t i = 0; i < grammar.terms_count; i++) {
			int32_t straight_idx = grammar_to_straight[i];
			(*out_edges)[straight_idx] = edges[grammar.nonterms_count + i];
			edges[grammar.nonterms_count + i] = NULL;
		}
	}

	if (!target_path) {
		if (reach_owned) {
			free_reach_matrices(reach_matrices, (int64_t)symbols_amount);
			reach_matrices = NULL;
			reach_owned = false;
		}

		if (cache) {
			cached_edges =
				(GrB_Matrix *)calloc(grammar.terms_count, sizeof(GrB_Matrix));

			if (cached_edges &&
				dup_matrices_into(cached_edges, *out_edges, grammar.terms_count) ==
					GrB_SUCCESS &&
				GrB_Matrix_dup(&cached_reach, *out_reachability) == GrB_SUCCESS) {

				if (mr_cache_insert_extracted(cache, key, grammar_tag, cached_edges,
											  grammar.terms_count,
											  cached_reach) == GrB_SUCCESS) {
					cached_edges = NULL;
					cached_reach = NULL;
				}
			}
		}
	}

cleanup:
	if (reach_owned) {
		free_reach_matrices(reach_matrices, (int64_t)symbols_amount);
	}
	free(grammar_to_straight);
	if (edges) {
		for (size_t i = 0; i < symbols_amount; i++) {
			GrB_free(&edges[i]);
		}
		free((void *)edges);
	}
	if (cached_edges) {
		for (int64_t t = 0; t < grammar.terms_count; t++) {
			GrB_free(&cached_edges[t]);
		}
		free((void *)cached_edges);
	}
	GrB_free(&cached_reach);
	if (adj_matrices) {
		for (int i = 0; i < grammar.nonterms_count; i++) {
			GrB_free(&adj_matrices[i]);
		}
		free((void *)adj_matrices);
	}
	return info;
}
