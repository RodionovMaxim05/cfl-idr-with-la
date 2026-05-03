#pragma once

#include <GraphBLAS.h>

#include "cfl_idr.h"

GrB_Info parse_graph(const char *filename, IdrGraph *out);
