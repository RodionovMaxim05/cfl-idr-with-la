#pragma once

#include <GraphBLAS.h>

#include "approximation/mr_graph.h"
#include "grammar/grammar.h"

GrB_Info apply_valueflow_under_approx(GrB_Matrix *paths, GrB_Matrix *adj_matrices,
									  MRGrammar_t grammar, const MRGraph *comp,
									  GrB_Matrix *comp_result, char *msg);

GrB_Info apply_valueflow_over_approx(const MRGraph *graph, GrB_Matrix *beta_reach,
									 MRGraph *filtered_graph, char *msg);
