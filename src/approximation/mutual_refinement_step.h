#pragma once

#include <GraphBLAS.h>
#include <stdbool.h>

#include "cfl_idr.h"
#include "mr_cache.h"
#include "utils/extract_edges.h"

/**
 * @brief Executes a single CFL-reachability and edge-extraction step with caching
 * support.
 *
 * Orchestrates one phase of the mutual-refinement algorithm using composite
 * cache keys `(graph_key, grammar_tag)` computed via `get_graph_cache_hash`:
 * - **Full-graph extraction mode (`target_path == NULL`)**:
 *   1. Checks the `extracted` cache via `mr_cache_lookup_extracted`; on a hit,
 *      duplicates the cached terminal edge matrices and non-trivial start-symbol
 *      reachability matrix.
 *   2. On a miss, computes per-symbol reachability via
 *      `LAGraph_CFL_reachability_adv`, extracts witness terminal edges across
 *      the entire graph via `LAGraph_CFL_extract_edges`, filters non-trivial
 *      start-symbol paths via `extract_non_trivial_paths`, and stores the
 *      results via `mr_cache_insert_extracted`.
 * - **Targeted pair mode (`target_path != NULL`)**:
 *   1. Checks the `reachability` cache via `mr_cache_lookup_reachability`; on a
 *      miss, computes per-symbol reachability via `LAGraph_CFL_reachability_adv`
 *      and caches the resulting matrix array via `mr_cache_insert_reachability`.
 *   2. Verifies whether `target_path` is reachable under the start nonterminal
 *      and sets `*target_found`.
 *   3. If reachable, extracts witness terminal edges restricted to
 *      `target_path` via `LAGraph_CFL_extract_edges` and outputs a single-entry
 *      reachability matrix; otherwise, sets `*out_reachability` and `*out_edges`
 *      entries to `NULL`.
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
