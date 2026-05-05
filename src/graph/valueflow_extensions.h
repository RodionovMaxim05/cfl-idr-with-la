#pragma once

#include <GraphBLAS.h>

#include "cfl_idr.h"

GrB_Info filter_bracket_paths(const IdrGraph *graph, GrB_Matrix paths,
							  GrB_Matrix *filtered_out);
