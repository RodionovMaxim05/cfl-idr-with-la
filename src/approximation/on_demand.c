#include "on_demand.h"

#include <stdlib.h>
#include <string.h>

#include "graph/condensate_graph.h"
#include "graph/remove_not_path.h"
#include "graph/split_into_components.h"
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

	GrB_Index n = graph->n;

	CondensationResult cr = {0};
	info = condensate_from_under_approx(graph, under_approx, &cr, msg);

	// Identify unknown resultPaths: candidates not yet confirmed by under-approx
	GrB_Matrix unknown_paths = NULL;
	info = compute_unknown_paths(over_approx, under_approx, n, &unknown_paths, msg);

	GrB_Index unknown_nnz = 0;
	GrB_Matrix_nvals(&unknown_nnz, unknown_paths);

	GrB_Index *u_rows = NULL, *u_cols = NULL;
	if (unknown_nnz > 0) {
		u_rows = malloc(unknown_nnz * sizeof(GrB_Index));
		u_cols = malloc(unknown_nnz * sizeof(GrB_Index));
		if (!u_rows || !u_cols) {
			info = GrB_OUT_OF_MEMORY;
			goto cleanup_unknown;
		}
		GrB_Matrix_extractTuples_BOOL(u_rows, u_cols, NULL, &unknown_nnz,
									  unknown_paths);
	}

	bool *is_root = calloc(unknown_nnz, sizeof(bool));
	GrB_Index *src_root = malloc(unknown_nnz * sizeof(GrB_Index));
	GrB_Index *tgt_root = malloc(unknown_nnz * sizeof(GrB_Index));
	if (!is_root || !src_root || !tgt_root) {
		info = GrB_OUT_OF_MEMORY;
		goto cleanup_classify;
	}

	GrB_Index root_count = 0;
	GrB_Index derived_count = 0;
	for (GrB_Index k = 0; k < unknown_nnz; k++) {
		GrB_Vector_extractElement_UINT64(&src_root[k], cr.components, u_rows[k]);
		GrB_Vector_extractElement_UINT64(&tgt_root[k], cr.components, u_cols[k]);

		bool is_root_path = (u_rows[k] == src_root[k] && u_cols[k] == tgt_root[k]);
		is_root[k] = is_root_path;
		if (is_root_path) {
			root_count++;
		} else {
			derived_count++;
		}
	}

	GrB_Index *ordered = malloc(unknown_nnz * sizeof(GrB_Index));
	if (!ordered) {
		info = GrB_OUT_OF_MEMORY;
		goto cleanup_classify;
	}

	GrB_Index ri = 0, di = root_count;
	for (GrB_Index k = 0; k < unknown_nnz; k++) {
		if (is_root[k]) {
			ordered[ri++] = k;
		} else {
			ordered[di++] = k;
		}
	}

	GrB_Matrix confirmed_pairs = NULL;
	GrB_Matrix_new(&confirmed_pairs, GrB_BOOL, cr.condensed_graph.n,
				   cr.condensed_graph.n);

	GrB_Matrix res = NULL;
	GrB_Matrix_dup(&res, under_approx);

	for (GrB_Index idx = 0; idx < unknown_nnz; idx++) {
		GrB_Index k = ordered[idx];
		GrB_Index sr = src_root[k];
		GrB_Index tr = tgt_root[k];

		if (!is_root[k]) {
			// Derived path: inherit decision of its root pair
			bool confirmed = false;
			GrB_Info qi =
				GrB_Matrix_extractElement_BOOL(&confirmed, confirmed_pairs, sr, tr);
			if (qi == GrB_SUCCESS && confirmed) {
				GrB_Matrix_setElement_BOOL(res, true, u_rows[k], u_cols[k]);
			}
			continue;
		}

		// Root path: run mutual refinement on the condensed graph
		GrB_Matrix mr_result = NULL;
		TargetPath target_path = {.src = sr, .tgt = tr};
		info = mutual_refinement(&cr.condensed_graph, grammar_type, &mr_result,
								 filter_empty, &target_path);
		if (info != GrB_SUCCESS) {
			goto cleanup_loop;
		}

		// Check if the specific (sr, tr) pair was confirmed
		GrB_Index mr_nnz = 0;
		GrB_Matrix_nvals(&mr_nnz, mr_result);

		bool pair_found = false;
		if (mr_nnz > 0) {
			bool val = false;
			GrB_Info qi = GrB_Matrix_extractElement_BOOL(&val, mr_result, sr, tr);
			pair_found = (qi == GrB_SUCCESS && val);
		}

		if (pair_found) {
			GrB_Matrix_setElement_BOOL(confirmed_pairs, true, sr, tr);
			GrB_Matrix_setElement_BOOL(res, true, u_rows[k], u_cols[k]);
		}

	cleanup_loop:
		GrB_Matrix_free(&mr_result);
		if (info != GrB_SUCCESS) {
			goto cleanup_pairs;
		}
	}

	*result = res;
	res = NULL;

cleanup_pairs:
	GrB_Matrix_free(&confirmed_pairs);
	GrB_Matrix_free(&res);
	free(ordered);
cleanup_classify:
	free(is_root);
	free(src_root);
	free(tgt_root);
cleanup_unknown:
	free(u_rows);
	free(u_cols);
	GrB_Matrix_free(&unknown_paths);
	condensation_result_free(&cr, msg);
	return info;
}

GrB_Info get_on_demand(const MRGraph *graph, GrB_Matrix under_approx,
					   GrB_Matrix over_approx, bool parityD, GrB_Matrix *result,
					   bool filter_empty, char *msg) {
	GrB_Info info = GrB_SUCCESS;

	// Step 1: Apply "classic" grammar refinement

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

	GrB_Matrix classic_paths = NULL;
	info = refine_mr_with_grammar(&reduced1, under_approx, over_approx, DEFAULT,
								  &classic_paths, filter_empty, msg);
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
		*result = classic_paths;
		return GrB_SUCCESS;
	}

	// Step 2: Apply "all" grammar refinement on the result of step 1

	MRGraph reduced2 = {0};
	bool reduced2_owned = true;
	if (is_all_pairs(classic_paths, graph->n)) {
		reduced2 = reduced1;
	} else {
		info = remove_not_path(&reduced1, classic_paths, &reduced2, msg);
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
	info = refine_mr_with_grammar(&reduced2, under_approx, classic_paths, ALL,
								  &final_paths, filter_empty, msg);
	if (reduced2_owned) {
		mr_graph_free(&reduced2);
	}
	GrB_Matrix_free(&classic_paths);
	if (info != GrB_SUCCESS) {
		return info;
	}

	*result = final_paths;
	return GrB_SUCCESS;
}
