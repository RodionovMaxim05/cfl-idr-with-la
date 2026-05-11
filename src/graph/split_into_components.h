#pragma once

#include <GraphBLAS.h>

#include "cfl_idr.h"

/**
 * @brief Splits an `IdrGraph` into its weakly-connected components for
 * independent analysis.
 *
 * This function decomposes the input graph into disjoint subgraphs, each
 * corresponding to a weakly-connected component (WCC) of the underlying
 * undirected graph.
 *
 * Algorithm:
 * 1. Builds a unified undirected adjacency matrix via `idr_graph_to_adjacency`,
 *    symmetrizing directed edges: `A_undir = A ∨ Aᵀ`.
 * 2. Computes connected components using `LAGr_ConnectedComponents`.
 * 3. Groups vertices by component ID and sorts for contiguous extraction.
 * 4. For each component with ≥2 vertices:
 *    - Extracts submatrices for all parenthesis/bracket/normal matrices via
 *      `GrB_Matrix_extract`, reindexing vertices to `[0, component_size)`.
 *    - Filters out empty parenthesis/bracket pairs (both open and close empty).
 *    - Stores the component graph and a vertex mapping array.
 * 5. Returns arrays of component graphs and their corresponding vertex maps.
 *
 * Vertex mapping semantics:
 * - `vertex_maps[c][local_idx]` gives the original vertex index in the input graph.
 * - Local indices within a component are contiguous and start at 0.
 *
 *
 * @param[in]  graph           Input graph to split. Must not be `NULL`.
 * @param[out] out_components  Array of `IdrGraph` structures, one per component
 *                             (length = `*out_count`). Ownership transferred to
 *                             caller.
 * @param[out] out_vertex_maps Array of vertex mapping arrays. Each
 *                             `out_vertex_maps[c]` has length equal to the number
 *                             of vertices in component `c`. Ownership transferred.
 * @param[out] out_count       Number of components produced (excluding isolated
 *                             single-vertex components).
 *
 * @return `GrB_SUCCESS` on success, or a GraphBLAS/LAGraph error code on failure.
 *
 * @note Single-vertex components (isolated vertices) are excluded from the output,
 *       as they cannot contribute to non-trivial CFL-reachable paths.
 */
GrB_Info split_IdrGraph_into_components(const IdrGraph *graph,
										IdrGraph **out_components,
										GrB_Index ***out_vertex_maps,
										GrB_Index *out_count);
