#pragma once

#include "grammar/grammar.h"

GrB_Info extractEdgesFromOutputs(GrB_Matrix *paths, GrB_Matrix *adj_matrices,
								 MRGrammar_t grammar, GrB_Index n, GrB_Index start_i,
								 GrB_Index start_j, GrB_Matrix *result_matrices,
								 char *msg);
