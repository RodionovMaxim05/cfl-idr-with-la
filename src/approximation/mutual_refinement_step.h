#pragma once

#include <GraphBLAS.h>
#include <stdbool.h>

#include "cfl_idr.h"
#include "mr_cache.h"
#include "utils/extract_edges.h"

/**
 * @brief Executes a single CFL-reachability computation step with caching
 * support.
 *
 * This function orchestrates one phase of the mutual-refinement algorithm:
 * 1. Computes a cache key via `get_graph_cache_hash`.
 * 2. Checks the cache for a precomputed result matching `(graph_key,
 * grammar_tag)`.
 * 3. If cached: extracts edges and reachability from stored results.
 * 4. If not cached: invokes `LAGraph_CFL_AllPaths`, extracts results, and
 * inserts them into the cache via `mr_cache_insert`.
 * 5. Optionally checks for the presence of a target path.
 *
 * @param[in]  graph             Input graph for CFL analysis.
 * @param[in]  grammar           Grammar specification (rules, nonterminals,
 * etc.).
 * @param[in]  grammar_tag       Tag identifying the grammar variant for caching.
 * @param[in]  target_path       Optional target path for early termination
 * checks.
 * @param[in]  cache             Cache for storing/retrieving intermediate
 * results.
 * @param[out] out_reachability  Output matrix: reachability via start symbol.
 * @param[out] out_edges         Array of edge matrices extracted for each
 * terminal.
 * @param[out] target_found      Output flag: `true` if target path was found.
 *
 * @return `GrB_SUCCESS` on success, or a GraphBLAS/LAGraph error code on
 * failure.
 *
 * @note `out_edges` must point to a pre-allocated array of sufficient size
 *       (at least `grammar.terms_count` elements).
 */
GrB_Info run_cfl_step(const IdrGraph *graph, MRGrammar_t grammar,
					  uint32_t grammar_tag, const TargetPath *target_path,
					  MRCache *cache, GrB_Matrix *out_reachability,
					  GrB_Matrix **out_edges, bool *target_found);
