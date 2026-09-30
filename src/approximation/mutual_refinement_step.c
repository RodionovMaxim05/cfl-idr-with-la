#include "mutual_refinement_step.h"

#include <LAGraphX.h>

#include "idr_graph.h"
#include "internal/grb_utils.h"
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
static bool matrix_has_path_udt(GrB_Matrix m, const TargetPath *target_path) {
	AllPathsElem val = {0};
	GrB_Info info =
		GrB_Matrix_extractElement_UDT(&val, m, target_path->src, target_path->tgt);
	return info == GrB_SUCCESS;
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
		GrB_Info info = GrB_Matrix_dup(&dst[i], src[i]);
		if (info != GrB_SUCCESS) {
			for (int64_t j = 0; j <= i; j++) {
				GrB_Matrix_free(&dst[j]);
			}
			return info;
		}
	}
	return GrB_SUCCESS;
}

/**
 * @brief Releases intermediate outputs produced during CFL path evaluation.
 *
 * Invokes `LAGraph_CFL_AllPaths_adv_free_outputs` on `paths`, frees `rule_table`,
 * and releases `all_paths_t`.
 *
 * @param[in,out] paths       Array of output matrices to free. May be `NULL`.
 * @param[in]     count       Number of matrices in `paths`.
 * @param[in,out] rule_table  Binary rule metadata array. May be `NULL`.
 * @param[in,out] all_paths_t GraphBLAS type descriptor to release.
 */
static void free_raw(GrB_Matrix *paths, int64_t count, BinaryRuleInfo *rule_table,
					 GrB_Type *all_paths_t) {
	if (!paths) {
		return;
	}
	LAGraph_CFL_AllPaths_adv_free_outputs(paths, count, all_paths_t);
	free(rule_table);
	GrB_free(all_paths_t);
}

GrB_Info run_cfl_step(const IdrGraph *graph, MRGrammar_t grammar,
					  uint32_t grammar_tag, const TargetPath *target_path,
					  MRCache *cache, GrB_Matrix *out_reachability,
					  GrB_Matrix **out_edges, bool *target_found) {
	GrB_Info info = GrB_SUCCESS;
	char msg[LAGRAPH_MSG_LEN];
	GrB_Matrix *adj = NULL;
	GrB_Type all_paths_t = NULL;
	GrB_Matrix *paths = NULL;
	BinaryRuleInfo *rule_table = NULL;
	bool raw_ready = false;
	bool raw_cached = false;
	GrB_Matrix *cached_edges = NULL;
	GrB_Matrix cached_reach = NULL;

	const uint64_t key = get_graph_cache_hash(graph);

	if (!target_path) {
		const MRRefinedResult *refined_hit =
			mr_cache_lookup_refined(cache, key, grammar_tag);

		if (refined_hit) {
			GRB_TRY(dup_matrices_into(*out_edges, refined_hit->edges,
									  refined_hit->edges_count));
			GRB_TRY(GrB_Matrix_dup(out_reachability, refined_hit->reachability));
			goto cleanup;
		}
	}

	GRB_TRY(idr_graph_collect_matrices(&adj, graph, grammar.nonterms_count,
									   &grammar.config, grammar.k));

	GrB_Matrix *raw_res = NULL;
	const MRStepResult *raw_hit =
		target_path ? mr_cache_lookup_raw(cache, key, grammar_tag) : NULL;

	if (raw_hit) {
		raw_res = raw_hit->matrices;
		rule_table = raw_hit->rule_table;
	} else {
		paths = (GrB_Matrix *)calloc(grammar.nonterms_count + grammar.terms_count,
									 sizeof(GrB_Matrix));
		if (!paths) {
			info = GrB_OUT_OF_MEMORY;
			goto cleanup;
		}

		GRB_TRY(LAGraph_CFL_AllPaths_adv(
			paths, &all_paths_t, &rule_table, adj,
			grammar.nonterms_count + grammar.terms_count, grammar.rules,
			grammar.rules_count, msg, OPT_EMPTY | OPT_BLOCK));

		raw_ready = true;
		raw_res = paths;
		rule_table = rule_table;
	}

	GRB_TRY(extract_edges_from_outputs(*out_edges, raw_res, rule_table, adj, grammar,
									   graph->n, graph->n_par, graph->n_bra,
									   &grammar.config, target_path));

	GRB_TRY(GrB_Matrix_new(out_reachability, GrB_BOOL, graph->n, graph->n));

	if (target_path) {
		*target_found = matrix_has_path_udt(raw_res[0], target_path);
		if (*target_found) {
			GRB_TRY(GrB_Matrix_setElement_BOOL(*out_reachability, true,
											   target_path->src, target_path->tgt));
		} else {
			GrB_Matrix_free(out_reachability);
		}

		if (!raw_hit && cache) {
			GrB_Info ins = mr_cache_insert_raw(
				cache, key, grammar_tag, paths, rule_table,
				grammar.nonterms_count + grammar.terms_count, all_paths_t);
			raw_cached = (ins == GrB_SUCCESS);
		}
	} else {
		GRB_TRY(extract_non_trivial_paths(out_reachability, raw_res[0]));

		if (cache) {
			cached_edges =
				(GrB_Matrix *)calloc(grammar.terms_count, sizeof(GrB_Matrix));

			if (cached_edges &&
				dup_matrices_into(cached_edges, *out_edges, grammar.terms_count) ==
					GrB_SUCCESS &&
				GrB_Matrix_dup(&cached_reach, *out_reachability) == GrB_SUCCESS) {

				GrB_Info ins =
					mr_cache_insert_refined(cache, key, grammar_tag, cached_edges,
											grammar.terms_count, cached_reach);
				if (ins == GrB_SUCCESS) {
					cached_edges = NULL;
					cached_reach = NULL;
				}
			}
		}
	}

cleanup:
	if (raw_ready && !raw_cached) {
		free_raw(paths, grammar.nonterms_count + grammar.terms_count, rule_table,
				 &all_paths_t);
	}
	if (cached_edges) {
		for (int64_t t = 0; t < grammar.terms_count; t++) {
			GrB_Matrix_free(&cached_edges[t]);
		}
		free((void *)cached_edges);
	}
	GrB_Matrix_free(&cached_reach);
	if (adj) {
		for (int i = 0; i < grammar.nonterms_count; i++) {
			GrB_free(&adj[i]);
		}
		free((void *)adj);
	}
	return info;
}
