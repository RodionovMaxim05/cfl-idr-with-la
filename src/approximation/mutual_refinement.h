#pragma once

#include <GraphBLAS.h>
#include <stdbool.h>
#include <stdint.h>

#include "cfl_idr.h"
#include "grammar/grammar.h"
#include "mr_cache.h"
#include "utils/extract_edges.h"

GrB_Info mutual_refinement(GrB_Matrix *result, const IdrGraph *graph,
						   IdrGrammarType grammar_type, bool valueflow,
						   bool filter_empty, MRCache *cache);

GrB_Info mutual_refinement_with_components(GrB_Matrix *result, IdrGraph *components,
										   GrB_Index **vertex_maps,
										   GrB_Index comp_count, GrB_Index global_n,
										   IdrGrammarType grammar_type,
										   bool valueflow, bool filter_empty,
										   const TargetPath *target_path,
										   MRCache *cache);
