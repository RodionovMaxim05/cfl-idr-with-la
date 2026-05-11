#pragma once

#include <GraphBLAS.h>

#include "cfl_idr.h"

/**
 * @brief Result of graph condensation via under-approximation analysis.
 *
 * Contains:
 * - `condensed_graph`: The graph with vertices merged according to
 *   strongly-connected components (SCCs) identified in the under-approximation.
 * - `components`: A `GrB_Vector` mapping each original vertex index to its
 *   representative SCC index (component ID). Used for expanding results back
 *   to the original graph space.
 *
 * @note The `components` vector has length `graph->n` and contains `uint64_t`
 *       values representing component IDs in the range `[0, num_components)`.
 */
typedef struct {
	IdrGraph condensed_graph; // Condensed graph with merged SCC vertices
	GrB_Vector components;	  // Vertex-to-component mapping vector
} CondensationResult;

/**
 * @brief Condenses a graph by merging vertices that are mutually reachable
 * in the under-approximation.
 *
 * This function performs graph condensation to improve the efficiency of
 * subsequent over-approximation analysis:
 * 1. Computes the symmetric closure of `under_approx` to identify mutually
 *    reachable vertex pairs: `mutual = under_approx ∧ under_approxᵀ`.
 * 2. Treats `mutual` as an undirected graph and computes connected components
 *    via `LAGr_ConnectedComponents`.
 * 3. Constructs a permutation matrix `P` where `P[rep(i)][i] = true` maps each
 *    vertex to its component representative.
 * 4. Condenses each adjacency matrix via matrix multiplication:
 *    `condensed[t] = P * adj[t] * Pᵀ`, effectively merging vertices within
 *    the same SCC.
 * 5. Builds a new `IdrGraph` from the condensed matrices with empty-pair
 *    filtering enabled.
 *
 * @param[out] out            Output structure to populate. Must be valid and
 *                            uninitialized.
 * @param[in]  graph          Original input graph. Must not be `NULL`.
 * @param[in]  under_approx   Under-approximation matrix (`GrB_BOOL`, n × n)
 *                            identifying confirmed reachable pairs.
 *
 * @return `GrB_SUCCESS` on success, or a GraphBLAS/LAGraph error code on failure.
 *
 * @note On success, `out->components` and `out->condensed_graph` are owned by
 *       the caller and must be freed via `condensation_result_free`.
 *       On failure, all internal allocations are cleaned up automatically.
 */
GrB_Info condensate_from_under_approx(CondensationResult *out, const IdrGraph *graph,
									  GrB_Matrix under_approx);

/**
 * @brief Releases all resources owned by a `CondensationResult`.
 *
 * After calling this function, the structure is zero-initialized and may be
 * safely reused or discarded.
 *
 * @param[in,out] cr  Result structure to free. May be `NULL` (no-op).
 */
void condensation_result_free(CondensationResult *cr);

/**
 * @brief Expands a reachability result from condensed space back to the original
 * vertex space.
 *
 * Given a result matrix computed on a condensed graph, this function maps it
 * back to the original graph by:
 * 1. Reconstructing the permutation matrix `P` from the `components` vector
 *    (same construction as in `condensate_from_under_approx`).
 * 2. Expanding via matrix multiplication: `expanded = Pᵀ * mr_result * P`.
 * 3. Adding intra-cluster edges: for any two vertices `i`, `j` belonging to
 *    the same SCC (`rep(i) == rep(j)` and `i != j`), sets `result[i][j] = true`.
 *    This ensures that mutual reachability within components is preserved.
 *
 * @param[out] result      Newly allocated output matrix (n × n, `GrB_BOOL`)
 *                         with expanded reachability information.
 * @param[in]  mr_result   Reachability matrix computed on the condensed graph
 *                         (size: num_components × num_components).
 * @param[in]  components  Vertex-to-component mapping vector from condensation.
 * @param[in]  n           Number of vertices in the original graph.
 *
 * @return `GrB_SUCCESS` on success, or a GraphBLAS error code on failure.
 */
GrB_Info expand_result(GrB_Matrix *result, GrB_Matrix mr_result,
					   GrB_Vector components, GrB_Index n);
