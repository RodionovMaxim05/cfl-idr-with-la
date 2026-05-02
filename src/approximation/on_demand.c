#include "on_demand.h"

#include <stdlib.h>
#include <string.h>

#include "graph/condensate_graph.h"
#include "graph/remove_not_path.h"
#include "graph/split_into_components.h"
#include "mr_cache.h"
#include "mutual_refinement.h"
#include "utils/extract_edges.h"
#include "utils/extract_paths.h"

static GrB_Info compute_unknown_paths(GrB_Matrix over_approx,
									  GrB_Matrix under_approx, GrB_Index n,
									  GrB_Matrix *unknown, char *msg) {
	GrB_Info info;
	GrB_Matrix_new(unknown, GrB_BOOL, n, n);

	return GrB_Matrix_assign(*unknown, under_approx, NULL, over_approx, GrB_ALL, n,
							 GrB_ALL, n, GrB_DESC_RSC);
}

static GrB_Info refine_mr_with_grammar(const MRGraph *graph, GrB_Matrix under_approx,
									   GrB_Matrix over_approx,
									   MRGrammarType grammar_type,
									   GrB_Matrix *result, bool filter_empty,
									   char *msg) {
	GrB_Info info = GrB_SUCCESS;

	MRCache cache = {0};
	mr_cache_init(&cache);
	GrB_Index n = graph->n;

	CondensationResult cr = {0};
	info = condensate_from_under_approx(graph, under_approx, &cr, msg);

	// Identify unknown resultPaths: candidates not yet confirmed by under-approx
	GrB_Matrix unknown_paths = NULL;
	info = compute_unknown_paths(over_approx, under_approx, n, &unknown_paths, msg);

	GrB_Index *comp_map = malloc(n * sizeof(GrB_Index));
	if (!comp_map) {
		info = GrB_OUT_OF_MEMORY;
		goto cleanup_base;
	}
	GrB_Vector_extractTuples_UINT64(NULL, comp_map, &n, cr.components);

	GrB_Index n_scc = cr.condensed_graph.n;
	GrB_Matrix S = NULL;
	GrB_Matrix_new(&S, GrB_BOOL, n, n_scc);

	GrB_Index *v_ids = malloc(n * sizeof(GrB_Index));
	for (GrB_Index i = 0; i < n; i++) {
		v_ids[i] = i;
	}

	GrB_Scalar s_true = NULL;
	GrB_Scalar_new(&s_true, GrB_BOOL);
	GrB_Scalar_setElement_BOOL(s_true, true);
	GxB_Matrix_build_Scalar(S, v_ids, comp_map, s_true, n);
	free(v_ids);

	// root_candidates = S^T * unknown_paths * S
	GrB_Matrix root_candidates = NULL;
	GrB_Matrix tmp_scc = NULL;
	GrB_Matrix_new(&tmp_scc, GrB_BOOL, n_scc, n);
	GrB_Matrix_new(&root_candidates, GrB_BOOL, n_scc, n_scc);
	GrB_mxm(tmp_scc, NULL, NULL, GrB_LOR_LAND_SEMIRING_BOOL, S, unknown_paths,
			GrB_DESC_T0);
	GrB_mxm(root_candidates, NULL, NULL, GrB_LOR_LAND_SEMIRING_BOOL, tmp_scc, S,
			NULL);
	GrB_Matrix_free(&tmp_scc);

	GrB_Matrix confirmed_roots = NULL;
	GrB_Matrix_new(&confirmed_roots, GrB_BOOL, n_scc, n_scc);

	MRGraph *components = NULL;
	GrB_Index **vertex_maps = NULL;
	GrB_Index comp_count = 0;
	split_MRGraph_into_components(&cr.condensed_graph, &components, &vertex_maps,
								  &comp_count, msg);

	GrB_Index n_root_pairs = 0;
	GrB_Matrix_nvals(&n_root_pairs, root_candidates);

	if (n_root_pairs > 0) {
		GrB_Index *r_rows = malloc(n_root_pairs * sizeof(GrB_Index));
		GrB_Index *r_cols = malloc(n_root_pairs * sizeof(GrB_Index));
		GrB_Matrix_extractTuples_BOOL(r_rows, r_cols, NULL, &n_root_pairs,
									  root_candidates);

		for (GrB_Index i = 0; i < n_root_pairs; i++) {
			if (i % 10 == 0) {
				mr_cache_free(&cache);
			}
			GrB_Index sr = r_rows[i];
			GrB_Index tr = r_cols[i];

			GrB_Matrix mr_result = NULL;
			TargetPath target_path = {.src = sr, .tgt = tr};

			info = mutual_refinement_with_components(
				components, vertex_maps, comp_count, n_scc, grammar_type, &mr_result,
				filter_empty, &target_path, &cache);

			if (info == GrB_SUCCESS && mr_result != NULL) {
				bool confirmed = false;
				GrB_Matrix_extractElement_BOOL(&confirmed, mr_result, sr, tr);
				if (confirmed) {
					GrB_Matrix_setElement_BOOL(confirmed_roots, true, sr, tr);
				}
			}
			GrB_Matrix_free(&mr_result);
		}
		free(r_rows);
		free(r_cols);
	}

	// res = under_approx | (S * confirmed_roots * S^T)
	GrB_Matrix res = NULL;
	GrB_Matrix_dup(&res, under_approx);

	GrB_Matrix tmp_n = NULL;
	GrB_Matrix_new(&tmp_n, GrB_BOOL, n, n_scc);

	// tmp_n = S * confirmed_roots
	GrB_mxm(tmp_n, NULL, NULL, GrB_LOR_LAND_SEMIRING_BOOL, S, confirmed_roots, NULL);
	// res = res | (tmp_n * S^T)
	GrB_mxm(res, NULL, GrB_LOR, GrB_LOR_LAND_SEMIRING_BOOL, tmp_n, S, GrB_DESC_T1);

	*result = res;

cleanup_split:
	for (GrB_Index c = 0; c < comp_count; c++) {
		mr_graph_free(&components[c]);
		free(vertex_maps[c]);
	}
	free(components);
	free((void *)vertex_maps);
	GrB_Matrix_free(&tmp_n);
	GrB_Matrix_free(&confirmed_roots);
	GrB_Matrix_free(&root_candidates);
	GrB_Matrix_free(&S);
	GrB_Scalar_free(&s_true);
	free(comp_map);
cleanup_base:
	GrB_Matrix_free(&unknown_paths);
	condensation_result_free(&cr, msg);
	mr_cache_free(&cache);
	return info;
}

GrB_Info get_on_demand(const MRGraph *graph, GrB_Matrix under_approx,
					   GrB_Matrix over_approx, bool parityD, GrB_Matrix *result,
					   bool filter_empty, char *msg) {
	GrB_Info info = GrB_SUCCESS;

	// Step 1: Apply "default" grammar refinement

	MRGraph reduced1 = {0};
	bool reduced1_owned = false;
	if (is_all_pairs(over_approx, graph->n)) {
		reduced1 = *graph;
	} else {
		info = remove_not_path(graph, over_approx, &reduced1, msg);
		if (info != GrB_SUCCESS) {
			return info;
		}
		reduced1_owned = true;
	}

	GrB_Matrix default_paths = NULL;
	info = refine_mr_with_grammar(&reduced1, under_approx, over_approx, DEFAULT,
								  &default_paths, filter_empty, msg);
	if (info != GrB_SUCCESS) {
		if (reduced1_owned) {
			mr_graph_free(&reduced1);
		}
		return info;
	}

	if (parityD) {
		if (reduced1_owned) {
			mr_graph_free(&reduced1);
		}
		*result = default_paths;
		return GrB_SUCCESS;
	}

	// Step 2: Apply "all" grammar refinement on the result of step 1

	MRGraph reduced2 = {0};
	bool reduced2_owned = true;
	if (is_all_pairs(default_paths, graph->n)) {
		reduced2 = reduced1;
	} else {
		info = remove_not_path(&reduced1, default_paths, &reduced2, msg);
		if (info != GrB_SUCCESS) {
			if (reduced1_owned) {
				mr_graph_free(&reduced1);
			}
			return info;
		}
		reduced2_owned = true;
	}

	if (reduced1_owned && reduced1_owned != reduced2_owned) {
		mr_graph_free(&reduced1);
	}

	GrB_Matrix final_paths = NULL;
	info = refine_mr_with_grammar(&reduced2, under_approx, default_paths, ALL,
								  &final_paths, filter_empty, msg);
	if (reduced2_owned) {
		mr_graph_free(&reduced2);
	}
	GrB_Matrix_free(&default_paths);
	if (info != GrB_SUCCESS) {
		return info;
	}

	*result = final_paths;
	return GrB_SUCCESS;
}
