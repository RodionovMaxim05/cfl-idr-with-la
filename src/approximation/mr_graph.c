#include "mr_graph.h"

void mr_graph_free(const MRGraph *g) {
	if (g == NULL) {
		printf("MRGraph: NULL\n");
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

void mr_graph_print(const MRGraph *g) {
	if (g == NULL) {
		printf("MRGraph: NULL\n");
		return;
	}

	printf("MRGraph: n=%llu, n_par=%ld, n_bra=%ld\n", (unsigned long long)g->n,
		   g->n_par, g->n_bra);

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
