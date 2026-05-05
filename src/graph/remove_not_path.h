#pragma once

#include <GraphBLAS.h>

#include "cfl_idr.h"

typedef struct {
	GrB_Index n_scc;
	GrB_Index *scc_ids;
	GrB_Matrix scc_reach;
} SccResult;

GrB_Info build_adjacency(const IdrGraph *graph, GrB_Matrix *adj_out);

GrB_Info compute_sccs(const IdrGraph *graph, GrB_Matrix adj, SccResult *out);

void scc_result_free(SccResult *r);

bool is_all_pairs(GrB_Matrix over_approx, GrB_Index n);

GrB_Info remove_not_path(const IdrGraph *graph, GrB_Matrix over_approx,
						 IdrGraph *out);
