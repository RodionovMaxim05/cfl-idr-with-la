#pragma once

#include <GraphBLAS.h>

#include "cfl_idr.h"

/**
 * @brief Filters a reachability matrix to retain only paths consistent with
 * value-flow semantics for bracket-labeled edges.
 *
 * This function enforces the core value-flow constraint: a valid path must contain
 * at least one store operation (`[i`) followed later by a matching load operation
 * (`]i`) for some bracket type `i`. Paths that violate this discipline are removed.
 *
 * Algorithm:
 * 1. Computes SCCs of the input graph and their transitive closure via
 *    `compute_sccs`.
 * 2. Expands SCC-level reachability to vertex space: `R_node = S * R_scc * Sᵀ`.
 * 3. For each bracket type `i`, computes paths of the form:
 *    `open_bra[i] * R_node * close_bra[i]`, i.e., edges where an opening bracket
 *    can reach a closing bracket through the given reachability relation.
 * 4. Accumulates results across all bracket types via logical OR.
 * 5. Intersects the accumulated "bracket-valid" paths with the input `paths`
 *    matrix to produce the final filtered result.
 *
 * The output satisfies: `out[u][v] = true` iff there exists a path from `u` to `v`
 * in `paths` that includes at least one properly ordered store-load pair.
 *
 * @param[out] out     Newly allocated `GrB_Matrix` (n × n, `GrB_BOOL`) holding
 *                     the filtered reachability result.
 * @param[in]  graph   Input `IdrGraph` containing bracket-labeled edges.
 *                     Must not be `NULL`.
 * @param[in]  paths   Input reachability matrix to filter (`GrB_BOOL`, n × n).
 *                     Must not be `NULL`.
 *
 * @return `GrB_SUCCESS` on success, or a GraphBLAS error code on failure.
 */
GrB_Info filter_bracket_paths(GrB_Matrix *out, const IdrGraph *graph,
							  GrB_Matrix paths);
