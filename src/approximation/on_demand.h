#pragma once

#include "GraphBLAS.h"

#include "grammar/grammar.h"
#include "mr_graph.h"

GrB_Info get_on_demand(const MRGraph *graph, GrB_Matrix under_approx,
					   GrB_Matrix over_approx, bool parityD, GrB_Matrix *result,
					   bool filter_empty, char *msg);
