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

	int64_t bra_base = 2 * n_par;
	int64_t actual_n_par = 0;
	int64_t actual_n_bra = 0;

	if (!filter_empty) {
		// Quick build without filtration

		actual_n_par = n_par;
		actual_n_bra = n_bra;
	} else {
		// Build with empty pair filtering

		for (int64_t i = 0; i < n_par; i++) {
			if (pair_nonempty(input_matrices[2 * i], input_matrices[2 * i + 1])) {
				actual_n_par++;
			}
		}
		for (int64_t i = 0; i < n_bra; i++) {
			if (pair_nonempty(input_matrices[bra_base + 2 * i],
							  input_matrices[bra_base + 2 * i + 1])) {
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
			out->open_par[i] = input_matrices[2 * i];
			out->close_par[i] = input_matrices[2 * i + 1];
		}
		for (int64_t i = 0; i < n_bra; i++) {
			out->open_bra[i] = input_matrices[bra_base + 2 * i];
			out->close_bra[i] = input_matrices[bra_base + 2 * i + 1];
		}
	} else {
		// Build with empty pair filtering

		int64_t p_idx = 0, b_idx = 0;
		for (int64_t i = 0; i < n_par; i++) {
			if (pair_nonempty(input_matrices[2 * i], input_matrices[2 * i + 1])) {
				out->open_par[p_idx] = input_matrices[2 * i];
				out->close_par[p_idx] = input_matrices[2 * i + 1];
				p_idx++;
			} else {
				GrB_Matrix_free(&input_matrices[2 * i]);
				GrB_Matrix_free(&input_matrices[2 * i + 1]);
			}
		}
		for (int64_t i = 0; i < n_bra; i++) {
			if (pair_nonempty(input_matrices[bra_base + 2 * i],
							  input_matrices[bra_base + 2 * i + 1])) {
				out->open_bra[b_idx] = input_matrices[bra_base + 2 * i];
				out->close_bra[b_idx] = input_matrices[bra_base + 2 * i + 1];
				b_idx++;
			} else {
				GrB_Matrix_free(&input_matrices[bra_base + 2 * i]);
				GrB_Matrix_free(&input_matrices[bra_base + 2 * i + 1]);
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

GrB_Info idr_graph_get_adj_matrices(GrB_Matrix **out, const IdrGraph *graph) {
	int64_t terms_count = get_terms_count(graph->n_par, graph->n_bra, graph->normal);
	GrB_Matrix *adj = (GrB_Matrix *)malloc(terms_count * sizeof(GrB_Matrix));
	if (!adj) {
		return GrB_OUT_OF_MEMORY;
	}

	for (int64_t i = 0; i < graph->n_par; i++) {
		adj[2 * i] = graph->open_par[i];
		adj[2 * i + 1] = graph->close_par[i];
	}
	for (int64_t i = 0; i < graph->n_bra; i++) {
		adj[2 * graph->n_par + 2 * i] = graph->open_bra[i];
		adj[2 * graph->n_par + 2 * i + 1] = graph->close_bra[i];
	}
	if (graph->normal != NULL) {
		adj[2 * graph->n_par + 2 * graph->n_bra] = graph->normal;
	}

	*out = adj;
	return GrB_SUCCESS;
}
