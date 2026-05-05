#pragma once

#include <GraphBLAS.h>

#include "cfl_idr.h"

GrB_Info split_IdrGraph_into_components(IdrGraph *graph, IdrGraph **out_components,
										GrB_Index ***out_vertex_maps,
										GrB_Index *out_count);
