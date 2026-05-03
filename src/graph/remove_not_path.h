#pragma once

#include <GraphBLAS.h>

#include "approximation/mr_graph.h"

typedef struct {
	GrB_Index n_scc;
	GrB_Index *scc_ids;
	GrB_Matrix scc_reach;
} SccResult;

GrB_Info build_adjacency(const MRGraph *graph, GrB_Matrix *adj_out);

GrB_Info compute_sccs(const MRGraph *graph, GrB_Matrix adj, SccResult *out,
					  char *msg);

void scc_result_free(SccResult *r);

bool is_all_pairs(GrB_Matrix over_approx, GrB_Index n);

GrB_Info remove_not_path(const MRGraph *graph, GrB_Matrix over_approx, MRGraph *out,
						 char *msg);
