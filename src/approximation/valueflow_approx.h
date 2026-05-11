#pragma once

#include <GraphBLAS.h>

#include "cfl_idr.h"
#include "grammar/grammar.h"

GrB_Info apply_valueflow_under_approx(GrB_Matrix *comp_result, GrB_Matrix *paths,
									  GrB_Matrix *adj_matrices, MRGrammar_t grammar,
									  const IdrGraph *comp);

GrB_Info apply_valueflow_over_approx(IdrGraph *out_graph, const IdrGraph *graph,
									 GrB_Matrix *reach);
