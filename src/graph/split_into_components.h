#pragma once

#include <GraphBLAS.h>

#include "approximation/mr_graph.h"

GrB_Info split_MRGraph_into_components(MRGraph *graph, MRGraph **out_components,
									   GrB_Index ***out_vertex_maps,
									   GrB_Index *out_count, char *msg);
