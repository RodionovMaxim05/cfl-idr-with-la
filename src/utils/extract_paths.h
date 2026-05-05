#pragma once

#include <GraphBLAS.h>

#include "cfl_idr.h"

GrB_Info extract_non_trivial_paths(GrB_Matrix paths, GrB_Matrix *result);
