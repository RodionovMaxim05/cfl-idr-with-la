#pragma once

#include <GraphBLAS.h>
#include <stdbool.h>
#include <stdint.h>

#include "cfl_idr.h"

static inline int64_t get_terms_count(int64_t n_par, int64_t n_bra,
									  GrB_Matrix normal) {
	return 2 * n_par + 2 * n_bra + (normal != NULL ? 1 : 0);
}

void idr_graph_print(const IdrGraph *g);

GrB_Index idr_graph_count_edges(const IdrGraph *graph);

GrB_Info build_idr_graph(IdrGraph *out, GrB_Matrix *input_matrices, int64_t n_par,
						 int64_t n_bra, bool has_normal, GrB_Index n,
						 bool filter_empty);

GrB_Info idr_graph_get_adj_matrices(GrB_Matrix **out, const IdrGraph *graph);
