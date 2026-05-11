#pragma once

#include <GraphBLAS.h>

#include "cfl_idr.h"

GrB_Info parse_graph(IdrGraph *out, const char *filename);
