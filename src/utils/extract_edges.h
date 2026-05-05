#pragma once

#include "grammar/grammar.h"

typedef struct {
	GrB_Index src;
	GrB_Index tgt;
} TargetPath;

GrB_Info extract_edges_from_outputs(GrB_Matrix *paths, GrB_Matrix *adj_matrices,
									MRGrammar_t grammar, GrB_Index n,
									TargetPath *target_path,
									GrB_Matrix *result_matrices);
