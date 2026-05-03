#pragma once

#include <GraphBLAS.h>

#include "approximation/mr_graph.h"
#include "remove_not_path.h"

GrB_Info remove_valueflow_unreachable(const MRGraph *graph, MRGraph *out, char *msg);

GrB_Info filter_bracket_paths(const MRGraph *graph, GrB_Matrix paths,
							  GrB_Matrix *filtered_out, char *msg);
