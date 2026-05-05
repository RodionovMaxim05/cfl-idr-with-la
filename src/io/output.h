#pragma once

#include <GraphBLAS.h>

#include "io/cli.h"

GrB_Info write_matrix_pairs(FILE *out, GrB_Matrix matrix, GrB_Index nvals);

void resolve_output_path(const Args *args, char *output_file, size_t size);
