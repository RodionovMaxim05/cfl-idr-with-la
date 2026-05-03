#pragma once

#include <GraphBLAS.h>

#include "mutual_refinement.h"

static inline int64_t get_terms_count(int64_t n_par, int64_t n_bra,
									  GrB_Matrix normal) {
	return 2 * n_par + 2 * n_bra + (normal != NULL ? 1 : 0);
}

GrB_Matrix *assemble_adj_matrices(const MRGraph *graph);

GrB_Info get_under_approx(const MRGraph *graph, bool valueflow, GrB_Matrix *result);

GrB_Info get_over_approx(const MRGraph *graph, MRGrammarType grammar_type,
						 GrB_Matrix under_approx, // NULL if not
						 GrB_Matrix *result, bool valueflow,
						 bool filter_empty); // whether to filter out empty paths
