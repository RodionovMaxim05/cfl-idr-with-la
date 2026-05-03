#pragma once

#include <GraphBLAS.h>
#include <stdbool.h>
#include <stdint.h>

#include "grammar/grammar.h"
#include "mr_cache.h"
#include "mr_graph.h"
#include "utils/extract_edges.h"

GrB_Info build_refined_graph(MRGraph *out, GrB_Matrix *result_matrices,
							 int64_t n_par, int64_t n_bra, bool has_normal,
							 GrB_Index n, bool filter_empty);

GrB_Info mutual_refinement(const MRGraph *graph, MRGrammarType grammar_type,
						   GrB_Matrix *result, bool valueflow, bool filter_empty,
						   MRCache *cache);

GrB_Info
mutual_refinement_with_components(MRGraph *components, GrB_Index **vertex_maps,
								  GrB_Index comp_count, GrB_Index global_n,
								  MRGrammarType grammar_type, GrB_Matrix *result,
								  bool valueflow, bool filter_empty,
								  const TargetPath *target_path, MRCache *cache);

void free_refined_graph(MRGraph *g);
