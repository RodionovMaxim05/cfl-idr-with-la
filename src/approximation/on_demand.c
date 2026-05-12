#include <GraphBLAS.h>
#include <stdlib.h>
#include <string.h>

#include "cfl_idr.h"
#include "grammar/grammar_analysis_utils.h"
#include "graph/condensate_graph.h"
#include "graph/remove_not_path.h"
#include "graph/split_into_components.h"
#include "internal/grb_utils.h"
#include "mr_cache.h"
#include "mutual_refinement.h"
#include "utils/extract_edges.h"
#include "utils/extract_paths.h"

/**
 * @brief Computes the set of paths that are in the over-approximation but not
 * yet confirmed by the under-approximation.
 *
 * @param[out] out           Newly allocated `GrB_Matrix` (n × n, `GrB_BOOL`)
 *                           holding the unknown paths.
 * @param[in]  over_approx   Matrix representing the current over-approximation.
 * @param[in]  under_approx  Matrix representing the current under-approximation.
 * @param[in]  n             Dimension of the square matrices.
 *
 * @return `GrB_SUCCESS` on success, or a GraphBLAS error code on failure.
 */
static GrB_Info compute_unknown_paths(GrB_Matrix *out, GrB_Matrix over_approx,
									  GrB_Matrix under_approx, GrB_Index n) {
	GrB_Info info = GrB_SUCCESS;

	GRB_TRY(GrB_Matrix_new(out, GrB_BOOL, n, n));
	GRB_TRY(GrB_Matrix_assign(*out, under_approx, NULL, over_approx, GrB_ALL, n,
							  GrB_ALL, n, GrB_DESC_RSC));
cleanup:
	return info;
}

/**
 * @brief Refines mutual-refinement results using grammar-based analysis on
 * condensed graph components.
 *
 * This function implements the core refinement loop for on-demand CFL-reachability:
 * 1. Condenses the input graph using the under-approximation to merge mutually
 *    reachable vertices (`condensate_from_under_approx`).
 * 2. Identifies "unknown" vertex pairs: edges present in `over_approx` but absent
 *    from `under_approx`.
 * 3. Projects unknown pairs onto the condensed SCC graph to find candidate SCC
 *    pairs requiring refinement.
 * 4. For each candidate SCC pair, invokes `mutual_refinement_with_components`
 *    with a specific grammar to determine if the path is truly reachable.
 * 5. Confirmed paths are merged back into the result via matrix multiplication
 *    with the component mapping matrix `S`.
 *
 * @param[out] result        On success, a newly allocated `GrB_Matrix` (n × n,
 *                           `GrB_BOOL`) with the refined reachability result.
 * @param[in]  graph         Original input graph. Must not be `NULL`.
 * @param[in]  under_approx  Current under-approximation matrix. Must not be `NULL`.
 * @param[in]  over_approx   Current over-approximation matrix. Must not be `NULL`.
 * @param[in]  grammar_type  Grammar variant for refinement (see `IdrGrammarType`).
 * @param[in]  valueflow     If `true`, apply value-flow specific constraints.
 * @param[in]  filter_empty  If `true`, filter empty parenthesis/bracket pairs
 *                           during graph construction.
 *
 * @return `GrB_SUCCESS` on success, or a GraphBLAS/LAGraph error code on failure.
 */
static GrB_Info refine_mr_with_grammar(GrB_Matrix *result, const IdrGraph *graph,
									   GrB_Matrix under_approx,
									   GrB_Matrix over_approx,
									   IdrGrammarType grammar_type, bool valueflow,
									   bool filter_empty) {
	GrB_Info info = GrB_SUCCESS;

	MRCache cache = {0};
	mr_cache_init(&cache);

	GrB_Index n = graph->n;
	CondensationResult cr = {0};
	GrB_Matrix unknown_paths = NULL;
	GrB_Index *comp_map = NULL;
	GrB_Matrix S = NULL;
	GrB_Scalar s_true = NULL;
	GrB_Index *v_ids = NULL;
	GrB_Matrix root_candidates = NULL;
	GrB_Matrix tmp_scc = NULL;
	GrB_Matrix confirmed_roots = NULL;
	GrB_Matrix tmp_n = NULL;
	GrB_Matrix res = NULL;
	GrB_Index *r_rows = NULL;
	GrB_Index *r_cols = NULL;
	IdrGraph *components = NULL;
	GrB_Index **vertex_maps = NULL;
	GrB_Index comp_count = 0;

	GRB_TRY(condensate_from_under_approx(&cr, graph, under_approx));

	// Identify unknown resultPaths: candidates not yet confirmed by under-approx
	GRB_TRY(compute_unknown_paths(&unknown_paths, over_approx, under_approx, n));

	comp_map = malloc(n * sizeof(GrB_Index));
	if (!comp_map) {
		info = GrB_OUT_OF_MEMORY;
		goto cleanup;
	}
	GRB_TRY(GrB_Vector_extractTuples_UINT64(NULL, comp_map, &n, cr.components));

	GrB_Index n_scc = cr.condensed_graph.n;

	GRB_TRY(GrB_Matrix_new(&S, GrB_BOOL, n, n_scc));

	v_ids = malloc(n * sizeof(GrB_Index));
	if (!v_ids) {
		info = GrB_OUT_OF_MEMORY;
		goto cleanup;
	}
	for (GrB_Index i = 0; i < n; i++) {
		v_ids[i] = i;
	}

	GRB_TRY(GrB_Scalar_new(&s_true, GrB_BOOL));
	GRB_TRY(GrB_Scalar_setElement_BOOL(s_true, true));
	GRB_TRY(GxB_Matrix_build_Scalar(S, v_ids, comp_map, s_true, n));
	free(v_ids);
	v_ids = NULL;
	GrB_Scalar_free(&s_true);

	// root_candidates = S^T * unknown_paths * S
	GRB_TRY(GrB_Matrix_new(&tmp_scc, GrB_BOOL, n_scc, n));
	GRB_TRY(GrB_Matrix_new(&root_candidates, GrB_BOOL, n_scc, n_scc));
	GRB_TRY(GrB_mxm(tmp_scc, NULL, NULL, GrB_LOR_LAND_SEMIRING_BOOL, S,
					unknown_paths, GrB_DESC_T0));
	GRB_TRY(GrB_mxm(root_candidates, NULL, NULL, GrB_LOR_LAND_SEMIRING_BOOL, tmp_scc,
					S, NULL));
	GrB_Matrix_free(&tmp_scc);

	GRB_TRY(GrB_Matrix_new(&confirmed_roots, GrB_BOOL, n_scc, n_scc));
	GRB_TRY(split_IdrGraph_into_components(&cr.condensed_graph, &components,
										   &vertex_maps, &comp_count));

	GrB_Index n_root_pairs = 0;
	GRB_TRY(GrB_Matrix_nvals(&n_root_pairs, root_candidates));

	if (n_root_pairs > 0) {
		r_rows = malloc(n_root_pairs * sizeof(GrB_Index));
		r_cols = malloc(n_root_pairs * sizeof(GrB_Index));
		if (!r_rows || !r_cols) {
			info = GrB_OUT_OF_MEMORY;
			goto cleanup;
		}

		GRB_TRY(GrB_Matrix_extractTuples_BOOL(r_rows, r_cols, NULL, &n_root_pairs,
											  root_candidates));

		for (GrB_Index i = 0; i < n_root_pairs; i++) {
			if (i % 10 == 0) {
				mr_cache_free(&cache);
			}

			GrB_Index src = r_rows[i];
			GrB_Index tgt = r_cols[i];
			GrB_Matrix mr_result = NULL;
			TargetPath target_path = {.src = src, .tgt = tgt};

			info = mutual_refinement_with_components(
				&mr_result, components, vertex_maps, comp_count, n_scc, grammar_type,
				valueflow, filter_empty, &target_path, &cache);

			if (info == GrB_SUCCESS && mr_result != NULL) {
				bool confirmed = false;
				GrB_Matrix_extractElement_BOOL(&confirmed, mr_result, src, tgt);
				if (confirmed) {
					GRB_TRY(
						GrB_Matrix_setElement_BOOL(confirmed_roots, true, src, tgt));
				}
			}
			GrB_Matrix_free(&mr_result);
		}
	}

	// res = under_approx | (S * confirmed_roots * S^T)
	GRB_TRY(GrB_Matrix_dup(&res, under_approx));
	GRB_TRY(GrB_Matrix_new(&tmp_n, GrB_BOOL, n, n_scc));

	// tmp_n = S * confirmed_roots
	GRB_TRY(GrB_mxm(tmp_n, NULL, NULL, GrB_LOR_LAND_SEMIRING_BOOL, S,
					confirmed_roots, NULL));
	// res = res | (tmp_n * S^T)
	GRB_TRY(GrB_mxm(res, NULL, GrB_LOR, GrB_LOR_LAND_SEMIRING_BOOL, tmp_n, S,
					GrB_DESC_T1));

	*result = res;
	res = NULL;

cleanup:
	for (GrB_Index c = 0; c < comp_count; c++) {
		idr_graph_free(&components[c]);
		free(vertex_maps[c]);
	}
	free(components);
	free((void *)vertex_maps);
	free(r_rows);
	free(r_cols);
	free(comp_map);
	free(v_ids);
	GrB_Matrix_free(&tmp_n);
	GrB_Matrix_free(&res);
	GrB_Matrix_free(&confirmed_roots);
	GrB_Matrix_free(&root_candidates);
	GrB_Matrix_free(&tmp_scc);
	GrB_Matrix_free(&S);
	GrB_Scalar_free(&s_true);
	GrB_Matrix_free(&unknown_paths);
	condensation_result_free(&cr);
	mr_cache_free(&cache);
	return info;
}

GrB_Info idr_get_on_demand(GrB_Matrix *result, const IdrGraph *graph,
						   GrB_Matrix under_approx, GrB_Matrix over_approx,
						   bool parity_d, bool valueflow, bool filter_empty) {
	GrB_Info info = GrB_SUCCESS;

	IdrGraph reduced1 = {0};
	bool reduced1_owned = false;
	IdrGraph reduced2 = {0};
	bool reduced2_owned = false;
	GrB_Matrix default_paths = NULL;
	GrB_Matrix final_paths = NULL;

	// Step 1: Apply "default" grammar refinement

	if (is_all_pairs(over_approx, graph->n)) {
		reduced1 = *graph;
	} else {
		GRB_TRY(remove_not_path(&reduced1, graph, over_approx));
		reduced1_owned = true;
	}

	GRB_TRY(refine_mr_with_grammar(&default_paths, &reduced1, under_approx,
								   over_approx, IDR_DEFAULT, valueflow,
								   filter_empty));

	if (parity_d) {
		*result = default_paths;
		default_paths = NULL;
		goto cleanup;
	}

	// Step 2: Apply "all" grammar refinement on the result of step 1

	if (is_all_pairs(default_paths, graph->n)) {
		reduced2 = reduced1;
	} else {
		GRB_TRY(remove_not_path(&reduced2, &reduced1, default_paths));
		reduced2_owned = true;
	}

	if (reduced1_owned && reduced2_owned) {
		idr_graph_free(&reduced1);
		reduced1_owned = false;
	}

	GRB_TRY(refine_mr_with_grammar(&final_paths, &reduced2, under_approx,
								   default_paths, IDR_ALL, valueflow, filter_empty));

	*result = final_paths;
	final_paths = NULL;

cleanup:
	if (reduced1_owned) {
		idr_graph_free(&reduced1);
	}
	if (reduced2_owned) {
		idr_graph_free(&reduced2);
	}
	GrB_Matrix_free(&default_paths);
	GrB_Matrix_free(&final_paths);
	return info;
}
