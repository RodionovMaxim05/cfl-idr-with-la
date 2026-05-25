#pragma once

#include <GraphBLAS.h>
#include <stdbool.h>
#include <stdint.h>

#include "cfl_idr.h"
#include "grammar/grammar_analysis_utils.h"
#include "mr_cache.h"
#include "utils/extract_edges.h"

/**
 * @brief Computes a refined over-approximation of reachable paths using the
 * mutual-refinement (MR) algorithm.
 *
 * This function implements the core mutual-refinement algorithm for CFL-reachability
 * analysis. It iteratively applies grammar-based transformations to converge on a
 * precise over-approximation of Dyck-reachable vertex pairs.
 *
 * The algorithm proceeds as follows for each connected component:
 * 1. **Alpha phase**: Analyze with Alpha grammar to get initial path set
 * 2. **Beta phase**: Analyze filtered graph with Beta grammar
 * 3. **Optional Project phase**: Apply projection grammar if specified
 * 4. **Optional Exclude phase**: Apply exclusion grammar if specified
 * 5. **Check convergence**: If graph unchanged, return intersection of all phases
 * 6. **Recursive refinement**: Otherwise, recurse on the refined graph
 *
 * @param[out] result        On success, a newly allocated `GrB_Matrix` (n × n,
 *                           `GrB_BOOL`) holding the over-approximation.
 * @param[in]  graph         Input graph to analyze. Must not be `NULL`.
 * @param[in]  grammar_type  Grammar variant controlling refinement behavior
 *                           (see `IdrGrammarType`).
 * @param[in]  valueflow     If `true`, apply value-flow specific post-processing
 *                           via `apply_valueflow_over_approx`.
 * @param[in]  filter_empty  If `true`, remove parenthesis/bracket types where
 *                           either the opening or closing matrix is empty.
 * @param[in]  target_path   Optional target path for on-demand analysis.
 *                           Pass `NULL` for full all-pairs analysis.
 * @param[in]  cache         Cache for intermediate CFL-reachability results.
 *
 * @return `GrB_SUCCESS` on success, or a GraphBLAS/LAGraph error code on failure.
 */
GrB_Info mutual_refinement(GrB_Matrix *result, const IdrGraph *graph,
						   IdrGrammarType grammar_type, bool valueflow,
						   bool filter_empty, const TargetPath *target_path,
						   MRCache *cache);

/**
 * @brief Runs mutual-refinement analysis on a pre-split set of graph components.
 *
 * This variant is used when the input graph has already been decomposed into
 * strongly-connected components (via `split_IdrGraph_into_components`).
 * It processes each component independently and merges results into a global
 * reachability matrix using the provided vertex mappings.
 *
 * When `target_path` is non-`NULL`, the function performs on-demand analysis:
 * - Only components containing both source and target vertices are processed.
 * - Each component is optionally reduced via `remove_not_path` to eliminate
 *   vertices not lying on any path between the local source/target.
 *
 * @param[out] result        Global output matrix (global_n × global_n, `GrB_BOOL`).
 * @param[in]  components    Array of `comp_count` component graphs.
 * @param[in]  vertex_maps   Array of vertex mappings: `vertex_maps[c][local_idx]`
 *                           maps to global vertex index.
 * @param[in]  comp_count    Number of components in the arrays.
 * @param[in]  global_n      Total number of vertices in the original graph.
 * @param[in]  grammar_type  Grammar variant for refinement (see `IdrGrammarType`).
 * @param[in]  valueflow     If `true`, apply value-flow post-processing.
 * @param[in]  filter_empty  If `true`, filter empty parenthesis/bracket pairs.
 * @param[in]  target_path   Optional target path for on-demand analysis.
 *                           Pass `NULL` for full all-pairs analysis.
 * @param[in]  cache         Cache for intermediate results.
 *
 * @return `GrB_SUCCESS` on success, or a GraphBLAS/LAGraph error code on failure.
 */
GrB_Info mutual_refinement_with_components(GrB_Matrix *result, IdrGraph *components,
										   GrB_Index **vertex_maps,
										   GrB_Index comp_count, GrB_Index global_n,
										   IdrGrammarType grammar_type,
										   bool valueflow, bool filter_empty,
										   const TargetPath *target_path,
										   MRCache *cache);
