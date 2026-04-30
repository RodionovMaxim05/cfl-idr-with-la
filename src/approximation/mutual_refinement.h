#pragma once

#include <stdbool.h>
#include <stdint.h>

#include "GraphBLAS.h"

#include "grammar/grammar.h"
#include "mr_graph.h"
#include "utils/extract_edges.h"

GrB_Info mutual_refinement(const MRGraph *graph, MRGrammarType grammar_type,
						   GrB_Matrix *result, bool filter_empty,
						   const TargetPath *target_path);

GrB_Info build_refined_graph(MRGraph *out, GrB_Matrix *result_matrices,
							 int64_t n_par, int64_t n_bra, bool has_normal,
							 GrB_Index n, bool filter_empty);

void free_refined_graph(MRGraph *g);
