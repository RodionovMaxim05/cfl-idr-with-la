#pragma once

#include <GraphBLAS.h>

#include "cfl_idr.h"

void extractNonTrivialPaths(GrB_Matrix paths, const IdrGraph *graph,
							GrB_Matrix *result);
