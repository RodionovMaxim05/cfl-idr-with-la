#pragma once

#include <GraphBLAS.h>

#include "cfl_idr.h"
#include "grammar/grammar.h"

/**
 * @brief Applies value-flow specific constraints to an under-approximation result.
 *
 * This function refines the computed under-approximation by enforcing store/load
 * consistency rules required for value-flow analysis.
 *
 * The refinement pipeline:
 * 1. Extracts terminal-labeled edges from CFL-reachability outputs via
 *    `extract_edges_from_outputs`.
 * 2. Reconstructs an `IdrGraph` from the extracted edges.
 * 3. Applies `filter_bracket_paths` to remove paths inconsistent with value-flow
 *    semantics.
 * 4. Replaces the input `comp_result` with the filtered matrix.
 *
 * @param[in,out] comp_result  On entry: under-approximation matrix to refine.
 *                             On exit: replaced with the filtered result.
 *                             Must not be `NULL`; ownership of the matrix is
 *                             transferred (old value is freed internally).
 * @param[in]     paths        Array of CFL-reachability output matrices (one per
 *                             grammar nonterminal). Must have at least
 *                             `grammar.nonterms_count` elements.
 * @param[in]     adj_matrices Input adjacency matrices for edge extraction.
 * @param[in]     grammar      Grammar specification used in the CFL analysis.
 * @param[in]     comp         Component graph being analyzed. Used for metadata
 *                             (parenthesis/bracket counts, vertex count).
 *
 * @return `GrB_SUCCESS` on success, or a GraphBLAS error code on failure.
 */
GrB_Info apply_valueflow_under_approx(GrB_Matrix *comp_result, GrB_Matrix *paths,
									  GrB_Matrix *adj_matrices, MRGrammar_t grammar,
									  const IdrGraph *comp);

/**
 * @brief Applies value-flow specific constraints to an over-approximation result.
 *
 * This function post-processes an over-approximation to enforce value-flow
 * discipline and prune unreachable vertices:
 * 1. Filters the reachability matrix via `filter_bracket_paths` to retain only
 *    paths consistent with store/load semantics.
 * 2. Removes vertices from the graph that cannot participate in any valid
 *    value-flow path via `idr_remove_valueflow_unreachable`.
 *
 * @param[out] out_graph  Output graph with unreachable vertices removed.
 *                        Must be a valid, uninitialized `IdrGraph` structure.
 * @param[in]  graph      Input graph to filter. Must not be `NULL`.
 * @param[in,out] reach   On entry: over-approximation matrix to refine.
 *                        On exit: replaced with the filtered reachability matrix.
 *                        Ownership is transferred (old value is freed internally).
 *
 * @return `GrB_SUCCESS` on success, or a GraphBLAS error code on failure.
 */
GrB_Info apply_valueflow_over_approx(IdrGraph *out_graph, const IdrGraph *graph,
									 GrB_Matrix *reach);
