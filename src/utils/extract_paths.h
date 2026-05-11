#pragma once

#include <GraphBLAS.h>

#include "cfl_idr.h"

GrB_Info extract_non_trivial_paths(GrB_Matrix *out, GrB_Matrix paths);
