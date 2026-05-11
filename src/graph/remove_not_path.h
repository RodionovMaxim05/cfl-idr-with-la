#pragma once

#include <GraphBLAS.h>

#include "cfl_idr.h"

/**
 * @brief Result of strongly-connected component (SCC) analysis on a graph.
 */
typedef struct {
	GrB_Index n_scc;	  // Number of strongly-connected components
	GrB_Index *scc_ids;	  // Vertex-to-SCC mapping array (length = graph->n)
	GrB_Matrix scc_reach; // SCC-level transitive closure matrix
} SccResult;

/**
 * @brief Computes strongly-connected components and their condensation reachability.
 *
 * This function performs SCC decomposition and builds the SCC-level transitive
 * closure:
 * 1. Uses `LAGraph_scc` to assign each vertex to an SCC (raw component IDs).
 * 2. Remaps raw IDs to a dense range `[0, n_scc)`.
 * 3. Constructs a selection matrix `S` (`n × n_scc`) where `S[i][scc_id[i]] = true`.
 * 4. Projects the input adjacency to SCC space: `reach_scc = Sᵀ * adj * S`.
 * 5. Computes the transitive closure of `reach_scc` via iterative squaring
 *    (`reach = reach ∨ reach² ∨ reach⁴ ∨ …`) until convergence.
 * 6. Adds reflexive edges to ensure each SCC reaches itself.
 *
 * Special cases:
 * - `n == 0`: Returns empty result with `n_scc = 0`.
 * - `n == 1`: Returns single-component result without invoking `LAGraph_scc`.
 *
 * @param[out] out    Output structure to populate. Must be valid and uninitialized.
 * @param[in]  graph  Input graph. Must not be `NULL`.
 * @param[in]  adj    Unified adjacency matrix (from `idr_graph_to_adjacency`).
 *                    Ownership remains with caller; not modified.
 *
 * @return `GrB_SUCCESS` on success, or a GraphBLAS/LAGraph error code on failure.
 */
GrB_Info compute_sccs(SccResult *out, const IdrGraph *graph, GrB_Matrix adj);

/**
 * @brief Releases all resources owned by an `SccResult`.
 *
 * @param[in,out] r  Result structure to free. May be `NULL` (no-op).
 */
void scc_result_free(SccResult *r);

/**
 * @brief Checks whether a boolean matrix represents all-pairs reachability.
 *
 * @param[in] over_approx  Boolean matrix to test (`GrB_BOOL`, n × n).
 * @param[in] n            Dimension of the square matrix.
 *
 * @return `true` if `nnz(over_approx) ≥ n²`, `false` otherwise.
 */
bool is_all_pairs(GrB_Matrix over_approx, GrB_Index n);

/**
 * @brief Removes vertices and edges from a graph that cannot participate in
 * any path represented by an over-approximation matrix.
 *
 * This optimization prunes the input graph to retain only vertices and edges
 * that lie on at least one path from a source to a target in `over_approx`.
 *
 * The algorithm:
 * 1. Builds a unified adjacency matrix via `idr_graph_to_adjacency`.
 * 2. Computes SCCs and their condensation reachability via `compute_sccs`.
 * 3. Projects `over_approx` to SCC space: `allowed_scc = Sᵀ * over_approx * S`.
 * 4. Computes indirect SCC reachability:
 *    `indirect = scc_reachᵀ * allowed_scc * scc_reachᵀ`.
 * 5. Expands back to vertex space: `M = S * indirect * Sᵀ`.
 * 6. Masks the original adjacency: `keep_mask = adj ∧ M`.
 * 7. Filters each labeled matrix in the `IdrGraph` via element-wise AND with
 *    `keep_mask`, producing the pruned output graph.
 *
 * @param[out] out           Output pruned graph. Must be valid and uninitialized.
 * @param[in]  graph         Input graph to filter. Must not be `NULL`.
 * @param[in]  over_approx   Over-approximation matrix (`GrB_BOOL`, n × n)
 *                           defining valid paths. Must not be `NULL`.
 *
 * @return `GrB_SUCCESS` on success, or a GraphBLAS error code on failure.
 */
GrB_Info remove_not_path(IdrGraph *out, const IdrGraph *graph,
						 GrB_Matrix over_approx);
