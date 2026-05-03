#pragma once

#include <GraphBLAS.h>

#include "cfl_idr.h"
#include "mutual_refinement.h"

static inline int64_t get_terms_count(int64_t n_par, int64_t n_bra,
									  GrB_Matrix normal) {
	return 2 * n_par + 2 * n_bra + (normal != NULL ? 1 : 0);
}

GrB_Matrix *assemble_adj_matrices(const IdrGraph *graph);
