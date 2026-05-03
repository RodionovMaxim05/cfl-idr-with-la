#pragma once

#include <GraphBLAS.h>

#include "cfl_idr.h"

typedef struct {
	IdrGraph condensed_graph;
	GrB_Vector components;
} CondensationResult;

GrB_Info condensate_from_under_approx(const IdrGraph *graph, GrB_Matrix under_approx,
									  CondensationResult *out, char *msg);

void condensation_result_free(CondensationResult *cr, char *msg);

GrB_Info expand_result(GrB_Matrix mr_result, GrB_Vector components, GrB_Index n,
					   GrB_Matrix *result, char *msg);
