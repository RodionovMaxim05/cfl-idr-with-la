#pragma once

#include <GraphBLAS.h>

#include "approximation/mr_graph.h"

typedef struct {
	MRGraph condensed_graph;
	GrB_Vector components;
} CondensationResult;

GrB_Info condensate_from_under_approx(const MRGraph *graph, GrB_Matrix under_approx,
									  CondensationResult *out, char *msg);

void condensation_result_free(CondensationResult *cr, char *msg);

GrB_Info expand_result(GrB_Matrix mr_result, GrB_Vector components, GrB_Index n,
					   GrB_Matrix *result, char *msg);
