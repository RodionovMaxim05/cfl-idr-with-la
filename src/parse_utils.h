#pragma once

#include "GraphBLAS.h"
#include "approximation/mr_graph.h"

GrB_Info parse_graph(const char *filename, MRGraph *out);
