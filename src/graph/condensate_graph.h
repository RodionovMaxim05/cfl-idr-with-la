#pragma once

#include <GraphBLAS.h>

#include "cfl_idr.h"

typedef struct {
	IdrGraph condensed_graph;
	GrB_Vector components;
} CondensationResult;

GrB_Info condensate_from_under_approx(CondensationResult *out, const IdrGraph *graph,
									  GrB_Matrix under_approx);

void condensation_result_free(CondensationResult *cr);

GrB_Info expand_result(GrB_Matrix *result, GrB_Matrix mr_result,
					   GrB_Vector components, GrB_Index n);
