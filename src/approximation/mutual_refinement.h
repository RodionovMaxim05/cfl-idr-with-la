#pragma once

#include <GraphBLAS.h>
#include <stdbool.h>
#include <stdint.h>

#include "cfl_idr.h"
#include "grammar/grammar.h"
#include "mr_cache.h"
#include "utils/extract_edges.h"

GrB_Info build_refined_graph(IdrGraph *out, GrB_Matrix *result_matrices,
							 int64_t n_par, int64_t n_bra, bool has_normal,
							 GrB_Index n, bool filter_empty);

GrB_Info mutual_refinement(const IdrGraph *graph, IdrGrammarType grammar_type,
						   GrB_Matrix *result, bool valueflow, bool filter_empty,
						   MRCache *cache);

GrB_Info
mutual_refinement_with_components(IdrGraph *components, GrB_Index **vertex_maps,
								  GrB_Index comp_count, GrB_Index global_n,
								  IdrGrammarType grammar_type, GrB_Matrix *result,
								  bool valueflow, bool filter_empty,
								  const TargetPath *target_path, MRCache *cache);

void free_refined_graph(IdrGraph *g);
