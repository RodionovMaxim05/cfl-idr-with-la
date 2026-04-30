#pragma once

#include <GraphBLAS.h>
#include <stdint.h>

typedef struct {
	GrB_Index n;
	int64_t n_par;
	GrB_Matrix *open_par;
	GrB_Matrix *close_par;
	int64_t n_bra;
	GrB_Matrix *open_bra;
	GrB_Matrix *close_bra;
	GrB_Matrix normal;
} MRGraph;

void mr_graph_free(const MRGraph *g);

void mr_graph_print(const MRGraph *g);
