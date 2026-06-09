#include "idr_graph.h"

#include "internal/grb_utils.h"

void idr_graph_free(IdrGraph *g) {
	if (!g) {
		return;
	}

	for (int64_t i = 0; i < g->n_par; i++) {
		GrB_Matrix_free(&g->open_par[i]);
		GrB_Matrix_free(&g->close_par[i]);
	}
	for (int64_t i = 0; i < g->n_bra; i++) {
		GrB_Matrix_free(&g->open_bra[i]);
		GrB_Matrix_free(&g->close_bra[i]);
	}
	if (g->normal != NULL) {
		GrB_Matrix_free(&g->normal);
	}

	free((void *)g->open_par);
	free((void *)g->close_par);
	free((void *)g->open_bra);
	free((void *)g->close_bra);
}

void idr_graph_print(const IdrGraph *g) {
	if (g == NULL) {
		printf("IdrGraph: NULL\n");
		return;
	}

	printf("IdrGraph: n=%ld, n_par=%ld, n_bra=%ld\n", g->n, g->n_par, g->n_bra);

	for (int64_t i = 0; i < g->n_par; i++) {
		printf("  open_par[%ld]:", i);
		GxB_print(g->open_par[i], 2);
		printf("  close_par[%ld]:", i);
		GxB_print(g->close_par[i], 2);
	}
	for (int64_t i = 0; i < g->n_bra; i++) {
		printf("  open_bra[%ld]:", i);
		GxB_print(g->open_bra[i], 2);
		printf("  close_bra[%ld]:", i);
		GxB_print(g->close_bra[i], 2);
	}

	if (g->normal != NULL) {
		printf("  normal:");
		GxB_print(g->normal, 2);
	} else {
		printf("  normal:\n  NULL\n");
	}
}

GrB_Index idr_graph_count_edges(const IdrGraph *graph) {
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

static inline bool pair_nonempty(GrB_Matrix open, GrB_Matrix close) {
	GrB_Index o = 0, c = 0;
	GrB_Matrix_nvals(&o, open);
	GrB_Matrix_nvals(&c, close);
	return o > 0 && c > 0;
}

GrB_Info build_idr_graph(IdrGraph *out, GrB_Matrix *input_matrices, int64_t n_par,
						 int64_t n_bra, bool has_normal, GrB_Index n,
						 bool filter_empty) {
	GrB_Info info = GrB_SUCCESS;

	out->n = n;
	out->n_par = 0;
	out->n_bra = 0;
	out->normal = NULL;
	out->open_par = NULL;
	out->close_par = NULL;
	out->open_bra = NULL;
	out->close_bra = NULL;

	int64_t open_par_base = 0;
	int64_t close_par_base = n_par;
	int64_t open_bra_base = 2 * n_par;
	int64_t close_bra_base = 2 * n_par + n_bra;

	int64_t actual_n_par = 0;
	int64_t actual_n_bra = 0;

	if (!filter_empty) {
		// Quick build without filtration

		actual_n_par = n_par;
		actual_n_bra = n_bra;
	} else {
		// Build with empty pair filtering

		for (int64_t i = 0; i < n_par; i++) {
			if (pair_nonempty(input_matrices[open_par_base + i],
							  input_matrices[close_par_base + i])) {
				actual_n_par++;
			}
		}
		for (int64_t i = 0; i < n_bra; i++) {
			if (pair_nonempty(input_matrices[open_bra_base + i],
							  input_matrices[close_bra_base + i])) {
				actual_n_bra++;
			}
		}
	}

	out->open_par = actual_n_par
						? (GrB_Matrix *)malloc(actual_n_par * sizeof(GrB_Matrix))
						: NULL;
	out->close_par = actual_n_par
						 ? (GrB_Matrix *)malloc(actual_n_par * sizeof(GrB_Matrix))
						 : NULL;
	out->open_bra = actual_n_bra
						? (GrB_Matrix *)malloc(actual_n_bra * sizeof(GrB_Matrix))
						: NULL;
	out->close_bra = actual_n_bra
						 ? (GrB_Matrix *)malloc(actual_n_bra * sizeof(GrB_Matrix))
						 : NULL;

	if ((actual_n_par && (!out->open_par || !out->close_par)) ||
		(actual_n_bra && (!out->open_bra || !out->close_bra))) {
		info = GrB_OUT_OF_MEMORY;
		goto cleanup;
	}

	if (!filter_empty) {
		// Quick build without filtration

		for (int64_t i = 0; i < n_par; i++) {
			out->open_par[i] = input_matrices[open_par_base + i];
			out->close_par[i] = input_matrices[close_par_base + i];
		}
		for (int64_t i = 0; i < n_bra; i++) {
			out->open_bra[i] = input_matrices[open_bra_base + i];
			out->close_bra[i] = input_matrices[close_bra_base + i];
		}
	} else {
		// Build with empty pair filtering

		int64_t p_idx = 0, b_idx = 0;
		for (int64_t i = 0; i < n_par; i++) {
			if (pair_nonempty(input_matrices[open_par_base + i],
							  input_matrices[close_par_base + i])) {
				out->open_par[p_idx] = input_matrices[open_par_base + i];
				out->close_par[p_idx] = input_matrices[close_par_base + i];
				p_idx++;
			} else {
				GrB_Matrix_free(&input_matrices[open_par_base + i]);
				GrB_Matrix_free(&input_matrices[close_par_base + i]);
			}
		}
		for (int64_t i = 0; i < n_bra; i++) {
			if (pair_nonempty(input_matrices[open_bra_base + i],
							  input_matrices[close_bra_base + i])) {
				out->open_bra[b_idx] = input_matrices[open_bra_base + i];
				out->close_bra[b_idx] = input_matrices[close_bra_base + i];
				b_idx++;
			} else {
				GrB_Matrix_free(&input_matrices[open_bra_base + i]);
				GrB_Matrix_free(&input_matrices[close_bra_base + i]);
			}
		}
	}

	out->n_par = actual_n_par;
	out->n_bra = actual_n_bra;
	out->normal = has_normal ? input_matrices[2 * n_par + 2 * n_bra] : NULL;
	return info;

cleanup:
	free((void *)out->open_par);
	free((void *)out->close_par);
	free((void *)out->open_bra);
	free((void *)out->close_bra);
	out->open_par = out->close_par = NULL;
	out->open_bra = out->close_bra = NULL;
	return info;
}

static void fill_grouped(GrB_Matrix *adj, GrB_Matrix *open, GrB_Matrix *close,
						 int64_t n, int64_t k, int64_t base_offset) {
	int64_t group_open_base[64] = {0};
	int64_t group_close_base[64] = {0};
	int64_t t_idx = base_offset;

	for (int64_t b = 0; b < k; b++) {
		group_open_base[b] = t_idx;
		int64_t size = (b < n) ? (n - b + k - 1) / k : 0;
		t_idx += size;
	}
	for (int64_t b = 0; b < k; b++) {
		group_close_base[b] = t_idx;
		int64_t size = (b < n) ? (n - b + k - 1) / k : 0;
		t_idx += size;
	}

	for (int64_t i = 0; i < n; i++) {
		int64_t b = i % k;
		int64_t j = i / k;
		adj[group_open_base[b] + j] = open[i];
		adj[group_close_base[b] + j] = close[i];
	}
}

GrB_Info idr_graph_collect_matrices(GrB_Matrix **out, const IdrGraph *graph,
									int64_t nonterms_count,
									bool is_beta_parity_group, int64_t k) {
	int64_t n_par = graph->n_par;
	int64_t n_bra = graph->n_bra;
	bool has_normal = (graph->normal != NULL);
	int64_t terms_count = get_terms_count(graph->n_par, graph->n_bra, graph->normal);
	size_t symbols_amount = (size_t)(terms_count + nonterms_count);
	GrB_Matrix *adj = (GrB_Matrix *)malloc(symbols_amount * sizeof(GrB_Matrix));
	if (!adj) {
		return GrB_OUT_OF_MEMORY;
	}

	for (int64_t i = 0; i < nonterms_count; i++) {
		GrB_Info info = GrB_Matrix_new(&adj[i], GrB_BOOL, graph->n, graph->n);
		if (info != GrB_SUCCESS) {
			for (int64_t j = 0; j < i; j++) {
				GrB_free(&adj[j]);
			}
			free((void *)adj);
			return info;
		}
	}

	int64_t t_offset = nonterms_count;

	if (!is_beta_parity_group) {
		// Parentheses (n_par) go linearly
		for (int64_t i = 0; i < n_par; i++) {
			adj[t_offset + i] = graph->open_par[i];
			adj[t_offset + n_par + i] = graph->close_par[i];
		}

		// Brackets (n_bra) are grouped by mask k
		fill_grouped(adj, graph->open_bra, graph->close_bra, n_bra, k,
					 t_offset + 2 * n_par);
	} else {
		// Parentheses (n_par) are grouped by mask k
		fill_grouped(adj, graph->open_par, graph->close_par, n_par, k, t_offset + 0);

		// Brackets (n_bra) go linearly
		for (int64_t i = 0; i < n_bra; i++) {
			adj[t_offset + 2 * n_par + i] = graph->open_bra[i];
			adj[t_offset + 2 * n_par + n_bra + i] = graph->close_bra[i];
		}
	}

	if (has_normal) {
		adj[t_offset + 2 * n_par + 2 * n_bra] = graph->normal;
	}

	*out = adj;
	return GrB_SUCCESS;
}

GrB_Info idr_graph_to_adjacency(GrB_Matrix *out, const IdrGraph *graph) {
	GrB_Info info = GrB_SUCCESS;
	GrB_Matrix adj = NULL;
	GrB_Index *row_indices = NULL;
	GrB_Index *col_indices = NULL;
	bool *edge_values = NULL;

	GrB_Index total_nvals = 0;

	// Count total number of non-zero elements in all matrices
	for (int64_t i = 0; i < graph->n_par; i++) {
		if (graph->open_par[i]) {
			GrB_Index nv;
			GRB_TRY(GrB_Matrix_nvals(&nv, graph->open_par[i]));
			total_nvals += nv;
		}
		if (graph->close_par[i]) {
			GrB_Index nv;
			GRB_TRY(GrB_Matrix_nvals(&nv, graph->close_par[i]));
			total_nvals += nv;
		}
	}
	for (int64_t i = 0; i < graph->n_bra; i++) {
		if (graph->open_bra[i]) {
			GrB_Index nv;
			GRB_TRY(GrB_Matrix_nvals(&nv, graph->open_bra[i]));
			total_nvals += nv;
		}
		if (graph->close_bra[i]) {
			GrB_Index nv;
			GRB_TRY(GrB_Matrix_nvals(&nv, graph->close_bra[i]));
			total_nvals += nv;
		}
	}
	if (graph->normal) {
		GrB_Index nv;
		GRB_TRY(GrB_Matrix_nvals(&nv, graph->normal));
		total_nvals += nv;
	}

	GRB_TRY(GrB_Matrix_new(&adj, GrB_BOOL, graph->n, graph->n));

	if (total_nvals == 0) {
		*out = adj;
		return GrB_SUCCESS;
	}

	row_indices = malloc(total_nvals * sizeof(GrB_Index));
	col_indices = malloc(total_nvals * sizeof(GrB_Index));
	edge_values = malloc(total_nvals * sizeof(bool));
	if (!row_indices || !col_indices || !edge_values) {
		info = GrB_OUT_OF_MEMORY;
		goto cleanup;
	}

	// Extract tuples from all submatrices into one buffer in a row
	GrB_Index offset = 0;

#define EXTRACT_TUPLES(matrix)                                                      \
	do {                                                                            \
		if (matrix) {                                                               \
			GrB_Index nv;                                                           \
			GRB_TRY(GrB_Matrix_nvals(&nv, matrix));                                 \
			if (nv > 0) {                                                           \
				GRB_TRY(GrB_Matrix_extractTuples_BOOL(                              \
					&row_indices[offset], &col_indices[offset],                     \
					&edge_values[offset], &nv, matrix));                            \
				offset += nv;                                                       \
			}                                                                       \
		}                                                                           \
	} while (0)

	for (int64_t i = 0; i < graph->n_par; i++) {
		EXTRACT_TUPLES(graph->open_par[i]);
		EXTRACT_TUPLES(graph->close_par[i]);
	}
	for (int64_t i = 0; i < graph->n_bra; i++) {
		EXTRACT_TUPLES(graph->open_bra[i]);
		EXTRACT_TUPLES(graph->close_bra[i]);
	}
	EXTRACT_TUPLES(graph->normal);

#undef EXTRACT_TUPLES

	// Construct the final matrix
	GRB_TRY(GrB_Matrix_build_BOOL(adj, row_indices, col_indices, edge_values,
								  total_nvals, GrB_LOR));

	*out = adj;
	adj = NULL;

cleanup:
	free(row_indices);
	free(col_indices);
	free(edge_values);
	GrB_Matrix_free(&adj);
	return info;
}
