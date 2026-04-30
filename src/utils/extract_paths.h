#pragma once

#include <GraphBLAS.h>

#include "approximation/mr_graph.h"

void extractNonTrivialPaths(GrB_Matrix paths, const MRGraph *graph,
							GrB_Matrix *result);
