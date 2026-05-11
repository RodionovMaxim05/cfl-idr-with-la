#pragma once

#include <GraphBLAS.h>

#include "cfl_idr.h"

GrB_Info filter_bracket_paths(GrB_Matrix *out, const IdrGraph *graph,
							  GrB_Matrix paths);
