#pragma once

#include <GraphBLAS.h>

#include "cfl_idr.h"

/**
 * @brief Extracts non-trivial paths from a CFL-reachability result matrix.
 *
 * Filters the input `paths` matrix to produce a boolean reachability matrix
 * containing only pairs `(i, j)` where `i ≠ j`.
 *
 * @param[out] out    Pointer to the output matrix (`GrB_BOOL`, n × n) to populate
 *                    with non-trivial reachability pairs. `*out` must be
 *                    pre-allocated to match the dimensions of `paths`.
 * @param[in]  paths  Input CFL-reachability result matrix (`GrB_BOOL`, n × n)
 *                    produced by `LAGraph_CFL_reachability_adv`. Must not be `NULL`.
 *
 * @return `GrB_SUCCESS` on success, or a GraphBLAS error code on failure.
 */
GrB_Info extract_non_trivial_paths(GrB_Matrix *out, GrB_Matrix paths);
