#include "mutual_refinement.h"

#include <LAGraph.h>
#include <LAGraphX.h>
#include <stdlib.h>
#include <string.h>

#include "approximation.h"
#include "grammar/grammar.h"
#include "grammar/grammar_analysis_utils.h"
#include "graph/remove_not_path.h"
#include "graph/split_into_components.h"
#include "utils/extract_paths.h"
#include "valueflow_approx.h"

static GrB_Index count_edges(const IdrGraph *graph) {
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

static bool matrix_has_path_udt(GrB_Matrix m, TargetPath *target_path) {
	AllPathsElem val = {0};
	GrB_Info info =
		GrB_Matrix_extractElement_UDT(&val, m, target_path->src, target_path->tgt);
	return info == GrB_SUCCESS;
}

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

GrB_Info build_refined_graph(IdrGraph *out, GrB_Matrix *result_matrices,
							 int64_t n_par, int64_t n_bra, bool has_normal,
							 GrB_Index n, bool filter_empty) {
	out->n = n;

	int64_t actual_n_par = 0;
	int64_t actual_n_bra = 0;

	if (!filter_empty) {
		// Quick build without filtration

		actual_n_par = n_par;
		actual_n_bra = n_bra;

		out->open_par =
			n_par ? (GrB_Matrix *)malloc(n_par * sizeof(GrB_Matrix)) : NULL;
		out->close_par =
			n_par ? (GrB_Matrix *)malloc(n_par * sizeof(GrB_Matrix)) : NULL;
		out->open_bra =
			n_bra ? (GrB_Matrix *)malloc(n_bra * sizeof(GrB_Matrix)) : NULL;
		out->close_bra =
			n_bra ? (GrB_Matrix *)malloc(n_bra * sizeof(GrB_Matrix)) : NULL;

		for (int64_t i = 0; i < n_par; i++) {
			out->open_par[i] = result_matrices[2 * i];
			out->close_par[i] = result_matrices[2 * i + 1];
		}

		int64_t bra_base = 2 * n_par;
		for (int64_t i = 0; i < n_bra; i++) {
			out->open_bra[i] = result_matrices[bra_base + 2 * i];
			out->close_bra[i] = result_matrices[bra_base + 2 * i + 1];
		}
	} else {
		// Build with empty pair filtering

		for (int64_t i = 0; i < n_par; i++) {
			GrB_Index n_open = 0, n_close = 0;
			GrB_Matrix_nvals(&n_open, result_matrices[2 * i]);
			GrB_Matrix_nvals(&n_close, result_matrices[2 * i + 1]);
			if (n_open > 0 && n_close > 0) {
				actual_n_par++;
			}
		}

		int64_t bra_base = 2 * n_par;
		for (int64_t i = 0; i < n_bra; i++) {
			GrB_Index n_open = 0, n_close = 0;
			GrB_Matrix_nvals(&n_open, result_matrices[bra_base + 2 * i]);
			GrB_Matrix_nvals(&n_close, result_matrices[bra_base + 2 * i + 1]);
			if (n_open > 0 && n_close > 0) {
				actual_n_bra++;
			}
		}

		out->open_par = actual_n_par
							? (GrB_Matrix *)malloc(actual_n_par * sizeof(GrB_Matrix))
							: NULL;
		out->close_par =
			actual_n_par ? (GrB_Matrix *)malloc(actual_n_par * sizeof(GrB_Matrix))
						 : NULL;
		out->open_bra = actual_n_bra
							? (GrB_Matrix *)malloc(actual_n_bra * sizeof(GrB_Matrix))
							: NULL;
		out->close_bra =
			actual_n_bra ? (GrB_Matrix *)malloc(actual_n_bra * sizeof(GrB_Matrix))
						 : NULL;

		int64_t p_idx = 0, b_idx = 0;
		for (int64_t i = 0; i < n_par; i++) {
			GrB_Index n_valid_open = 0, n_valid_close = 0;
			GrB_Matrix_nvals(&n_valid_open, result_matrices[2 * i]);
			GrB_Matrix_nvals(&n_valid_close, result_matrices[2 * i + 1]);
			if (n_valid_open > 0 && n_valid_close > 0) {
				out->open_par[p_idx] = result_matrices[2 * i];
				out->close_par[p_idx] = result_matrices[2 * i + 1];
				p_idx++;
			}
		}
		for (int64_t i = 0; i < n_bra; i++) {
			GrB_Index n_valid_open = 0, n_valid_close = 0;
			GrB_Matrix_nvals(&n_valid_open, result_matrices[bra_base + 2 * i]);
			GrB_Matrix_nvals(&n_valid_close, result_matrices[bra_base + 2 * i + 1]);
			if (n_valid_open > 0 && n_valid_close > 0) {
				out->open_bra[b_idx] = result_matrices[bra_base + 2 * i];
				out->close_bra[b_idx] = result_matrices[bra_base + 2 * i + 1];
				b_idx++;
			}
		}
	}

	out->n_par = actual_n_par;
	out->n_bra = actual_n_bra;
	out->normal = has_normal ? result_matrices[2 * n_par + 2 * n_bra] : NULL;

	return GrB_SUCCESS;
}

void free_refined_graph(IdrGraph *graph) {
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

static GrB_Info run_cfl_step(const IdrGraph *graph, MRGrammar_t grammar,
							 uint32_t grammar_tag, GrB_Matrix *out_reachability,
							 GrB_Matrix **out_edges, const TargetPath *target_path,
							 bool *target_found, MRCache *cache, char *msg) {
	GrB_Info info = GrB_SUCCESS;

	uint64_t key = get_graph_cache_hash(graph);
	const MRStepResult *hit = mr_cache_lookup(cache, key, grammar_tag);

	GrB_Matrix *adj = assemble_adj_matrices(graph);

	if (hit) {
		info = extractEdgesFromOutputs(hit->matrices, adj, grammar, graph->n,
									   target_path, *out_edges, msg);
		if (info != GrB_SUCCESS) {
			return info;
		}
		free((void *)adj);

		GrB_Matrix_new(out_reachability, GrB_BOOL, graph->n, graph->n);
		extractNonTrivialPaths(hit->matrices[0], graph, out_reachability);

		if (target_path) {
			*target_found = matrix_has_path_udt(hit->matrices[0], target_path);
			if (!*target_found) {
				GrB_Matrix_free(out_reachability);
			}
		}
		return GrB_SUCCESS;
	}

	GrB_Type all_paths_t = NULL;
	GrB_Matrix *paths = NULL;
	LAGraph_Calloc((void **)&paths, grammar.nonterms_count, sizeof(GrB_Matrix), msg);

	info = LAGraph_CFL_AllPaths(paths, &all_paths_t, adj, grammar.terms_count,
								grammar.nonterms_count, grammar.rules,
								grammar.rules_count, msg, 0);
	if (info != GrB_SUCCESS) {
		goto cleanup;
	}

	info = extractEdgesFromOutputs(paths, adj, grammar, graph->n, target_path,
								   *out_edges, msg);
	if (info != GrB_SUCCESS) {
		goto cleanup;
	}

	GrB_Matrix_new(out_reachability, GrB_BOOL, graph->n, graph->n);
	extractNonTrivialPaths(paths[0], graph, out_reachability);

	mr_cache_insert(cache, key, grammar_tag, *out_reachability, paths,
					grammar.nonterms_count, all_paths_t);

	if (target_path) {
		*target_found = matrix_has_path_udt(paths[0], target_path);
	}

cleanup:
	free((void *)adj);
	return info;
}

GrB_Info mutual_refinement_single(const IdrGraph *graph, IdrGrammarType grammar_type,
								  GrB_Matrix *result, bool valueflow,
								  bool filter_empty, const TargetPath *target_path,
								  MRCache *cache) {
	GrB_Info info = GrB_SUCCESS;
	char msg[LAGRAPH_MSG_LEN];
	bool has_normal = (graph->normal != NULL);
	bool target_found = false;

	GrB_Index initial_edge_count = count_edges(graph);
	if (initial_edge_count == 0) {
		GrB_Matrix_new(result, GrB_BOOL, graph->n, graph->n);
		return info;
	}

	int64_t terms_count = get_terms_count(graph->n_par, graph->n_bra, graph->normal);

	// Alpha phase

	MRGrammar_t alpha_grammar =
		get_alpha_grammar(grammar_type, graph->n_par, graph->n_bra, has_normal);

	GrB_Matrix alpha_reach = NULL;
	GrB_Matrix *alpha_edges = (GrB_Matrix *)malloc(terms_count * sizeof(GrB_Matrix));

	info = run_cfl_step(graph, alpha_grammar, GRAMMAR_TAG_ALPHA, &alpha_reach,
						&alpha_edges, target_path, &target_found, cache, msg);
	grammar_free(&alpha_grammar);
	if (info != GrB_SUCCESS) {
		goto cleanup_alpha;
	}

	IdrGraph alpha_graph = {0};

	if (target_path != NULL && !target_found) {
		GrB_Matrix_new(result, GrB_BOOL, graph->n, graph->n);
		goto cleanup_alpha;
	}

	build_refined_graph(&alpha_graph, alpha_edges, graph->n_par, graph->n_bra,
						has_normal, graph->n, filter_empty);

	GrB_Index alpha_edge_count = count_edges(&alpha_graph);
	if (alpha_edge_count == 0) {
		GrB_Matrix_new(result, GrB_BOOL, graph->n, graph->n);
		goto cleanup_alpha;
	}

	// Beta phase

	MRGrammar_t beta_grammar = get_beta_grammar(grammar_type, alpha_graph.n_par,
												alpha_graph.n_bra, has_normal);

	GrB_Matrix beta_reach = NULL;
	int64_t alpha_terms_count =
		get_terms_count(alpha_graph.n_par, alpha_graph.n_bra, graph->normal);
	GrB_Matrix *beta_edges =
		(GrB_Matrix *)malloc(alpha_terms_count * sizeof(GrB_Matrix));

	info = run_cfl_step(&alpha_graph, beta_grammar, GRAMMAR_TAG_BETA, &beta_reach,
						&beta_edges, target_path, &target_found, cache, msg);
	grammar_free(&beta_grammar);
	if (info != GrB_SUCCESS) {
		goto cleanup_beta;
	}

	IdrGraph beta_graph = {0};

	if (target_path != NULL && !target_found) {
		GrB_Matrix_new(result, GrB_BOOL, graph->n, graph->n);
		goto cleanup_beta;
	}

	build_refined_graph(&beta_graph, beta_edges, alpha_graph.n_par,
						alpha_graph.n_bra, has_normal, graph->n, filter_empty);

	GrB_Index beta_edge_count = count_edges(&beta_graph);
	if (beta_edge_count == 0) {
		GrB_Matrix_new(result, GrB_BOOL, graph->n, graph->n);
		goto cleanup_beta;
	}

	// Optional: Project phase

	GrB_Matrix project_reach = NULL;
	GrB_Matrix *project_edges = NULL;
	IdrGraph project_graph = {0};
	if (grammar_type == IDR_PROJECT || grammar_type == IDR_ALL) {
		MRGrammar_t project_grammar =
			dyck_project_grammar(beta_graph.n_par, beta_graph.n_bra, has_normal);

		int64_t beta_terms_count =
			get_terms_count(beta_graph.n_par, beta_graph.n_bra, graph->normal);
		project_edges = (GrB_Matrix *)malloc(beta_terms_count * sizeof(GrB_Matrix));

		info = run_cfl_step(&beta_graph, project_grammar, GRAMMAR_TAG_PROJECT,
							&project_reach, &project_edges, target_path,
							&target_found, cache, msg);
		grammar_free(&project_grammar);
		if (info != GrB_SUCCESS) {
			goto cleanup_project;
		}

		if (target_path != NULL && !target_found) {
			GrB_Matrix_new(result, GrB_BOOL, graph->n, graph->n);
			goto cleanup_project;
		}

		build_refined_graph(&project_graph, project_edges, beta_graph.n_par,
							beta_graph.n_bra, has_normal, graph->n, filter_empty);
	}

	// Optional: Exclude phase

	const IdrGraph *final_graph =
		(grammar_type == IDR_PROJECT || grammar_type == IDR_ALL) ? &project_graph
																 : &beta_graph;
	GrB_Matrix *exclude_edges = NULL;
	IdrGraph exclude_graph = {0};
	if (grammar_type == IDR_EXCLUDE || grammar_type == IDR_ALL) {
		const IdrGraph *cur_graph = final_graph;

		for (int64_t ex_bra = 0; ex_bra < cur_graph->n_bra; ex_bra++) {
			MRGrammar_t exclude_grammar = dyck_alpha_grammar_k_parity_exclude(
				cur_graph->n_par, cur_graph->n_bra, has_normal,
				/*k=*/2, ex_bra);

			int64_t cur_terms_count =
				get_terms_count(cur_graph->n_par, cur_graph->n_bra, graph->normal);
			exclude_edges =
				(GrB_Matrix *)malloc(cur_terms_count * sizeof(GrB_Matrix));

			GrB_Matrix cur_reach = NULL;

			info =
				run_cfl_step(cur_graph, exclude_grammar,
							 GRAMMAR_TAG_EXCLUDE_BASE + (uint32_t)ex_bra, &cur_reach,
							 &exclude_edges, target_path, &target_found, cache, msg);
			grammar_free(&exclude_grammar);
			GrB_Matrix_free(&cur_reach);
			if (info != GrB_SUCCESS) {
				goto cleanup_exclude;
			}

			if (target_path != NULL && !target_found) {
				GrB_Matrix_new(result, GrB_BOOL, graph->n, graph->n);
				goto cleanup_exclude;
			}

			free_refined_graph(&exclude_graph);
			build_refined_graph(&exclude_graph, exclude_edges, cur_graph->n_par,
								cur_graph->n_bra, has_normal, cur_graph->n,
								filter_empty);

			free((void *)exclude_edges);
			exclude_edges = NULL;
			cur_graph = &exclude_graph;
		}
		final_graph = cur_graph;
	}

	IdrGraph filtered_graph = {0};
	if (valueflow) {
		info = apply_valueflow_over_approx(final_graph, &beta_reach, &filtered_graph,
										   msg);
		if (info != GrB_SUCCESS) {
			free_refined_graph(&filtered_graph);
			goto cleanup_exclude;
		}
	}

	GrB_Index final_edge_count = count_edges(final_graph);

	// Check convergence
	if (final_edge_count == 0 || initial_edge_count == final_edge_count) {
		// Stability has been achieved

		if (target_path != NULL) {
			*result = alpha_reach;
			alpha_reach = NULL;
			goto cleanup_exclude;
		}

		info = matrix_intersect3(result, alpha_reach, beta_reach, project_reach,
								 graph->n);
	} else {
		// Recursive refinement
		info = mutual_refinement_single(final_graph, grammar_type, result, valueflow,
										filter_empty, target_path, cache);
	}

	free_refined_graph(&filtered_graph);

cleanup_exclude:
	free((void *)exclude_edges);
	free_refined_graph(&exclude_graph);

cleanup_project:
	GrB_Matrix_free(&project_reach);
	free((void *)project_edges);
	free_refined_graph(&project_graph);

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

GrB_Info
mutual_refinement_with_components(IdrGraph *components, GrB_Index **vertex_maps,
								  GrB_Index comp_count, GrB_Index global_n,
								  IdrGrammarType grammar_type, GrB_Matrix *result,
								  bool valueflow, bool filter_empty,
								  const TargetPath *target_path, MRCache *cache) {
	GrB_Info info = GrB_SUCCESS;
	char msg[LAGRAPH_MSG_LEN];

	GrB_Matrix_new(result, GrB_BOOL, global_n, global_n);

	for (GrB_Index c = 0; c < comp_count; c++) {
		IdrGraph *comp = &components[c];
		GrB_Index *vmap = vertex_maps[c];

		GrB_Matrix comp_result = NULL;

		if (target_path != NULL) {
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
			if (local_src == GrB_INVALID_VALUE || local_tgt == GrB_INVALID_VALUE) {
				continue;
			}

			GrB_Matrix single_path = NULL;
			GrB_Matrix_new(&single_path, GrB_BOOL, comp->n, comp->n);
			GrB_Matrix_setElement_BOOL(single_path, true, local_src, local_tgt);

			IdrGraph reduced_comp = {0};
			bool reduced_owned = false;
			if (!is_all_pairs(single_path, comp->n)) {
				info = remove_not_path(comp, single_path, &reduced_comp, msg);
				if (info != GrB_SUCCESS) {
					GrB_Matrix_free(&single_path);
					goto cleanup;
				}
				reduced_owned = true;
			} else {
				reduced_comp = *comp;
			}
			GrB_Matrix_free(&single_path);

			TargetPath local_target = {.src = local_src, .tgt = local_tgt};
			info = mutual_refinement_single(&reduced_comp, grammar_type,
											&comp_result, valueflow, filter_empty,
											&local_target, cache);
			if (reduced_owned) {
				idr_graph_free(&reduced_comp);
			}
		} else {
			info = mutual_refinement_single(comp, grammar_type, &comp_result,
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
			GrB_Matrix_extractTuples_BOOL(rows, cols, NULL, &nnz, comp_result);

			for (GrB_Index k = 0; k < nnz; k++) {
				GrB_Matrix_setElement_BOOL(*result, true, vmap[rows[k]],
										   vmap[cols[k]]);
			}
			free(rows);
			free(cols);
		}
		GrB_Matrix_free(&comp_result);
	}

cleanup:
	return info;
}

GrB_Info mutual_refinement(const IdrGraph *graph, IdrGrammarType grammar_type,
						   GrB_Matrix *result, bool valueflow, bool filter_empty,
						   MRCache *cache) {
	GrB_Info info = GrB_SUCCESS;
	char msg[LAGRAPH_MSG_LEN];

	IdrGraph *components = NULL;
	GrB_Index **vertex_maps = NULL;
	GrB_Index comp_count = 0;

	info = split_IdrGraph_into_components(graph, &components, &vertex_maps,
										  &comp_count, msg);
	if (info != GrB_SUCCESS) {
		return info;
	}

	info = mutual_refinement_with_components(components, vertex_maps, comp_count,
											 graph->n, grammar_type, result,
											 valueflow, filter_empty, NULL, cache);

	for (GrB_Index c = 0; c < comp_count; c++) {
		idr_graph_free(&components[c]);
		free(vertex_maps[c]);
	}
	free(components);
	free((void *)vertex_maps);
	return info;
}
