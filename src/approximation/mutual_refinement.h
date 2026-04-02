#pragma once

#include "GraphBLAS.h"
#include "grammar/grammar.h"
#include "mr_graph.h"
#include <stdbool.h>
#include <stdint.h>

GrB_Info mutual_refinement(const MRGraph *graph, MRGrammarType grammar_type,
						   GrB_Matrix *result);

GrB_Info build_refined_graph(MRGraph *out, GrB_Matrix *result_matrices,
							 int64_t n_par, int64_t n_bra, bool has_normal,
							 GrB_Index n);

void free_refined_graph(MRGraph *g);
