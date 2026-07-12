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

GrB_Info run_cfl_step(const IdrGraph *graph, MRGrammar_t grammar,
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
