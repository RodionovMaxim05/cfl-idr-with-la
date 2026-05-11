#pragma once

#include <GraphBLAS.h>

#include "cfl_idr.h"

/**
 * @brief Extracts non-trivial paths from a CFL-reachability result matrix.
 *
 * Filters the input `paths` matrix to produce a boolean reachability matrix
 * containing only pairs `(i, j)` where `i ≠ j`.
 *
 * @param[out] out     Output matrix (`GrB_BOOL`, n × n) to populate with non-trivial
 *                     reachability pairs. Must be pre-allocated by the caller with
 *                     the same dimensions as `paths`.
 * @param[in]  paths   Input CFL-reachability result matrix (`AllPathsElem` UDT type)
 *                     from `LAGraph_CFL_AllPaths`. Must not be `NULL`.
 *
 * @return `GrB_SUCCESS` on success, or a GraphBLAS error code on failure.
 */
GrB_Info extract_non_trivial_paths(GrB_Matrix *out, GrB_Matrix paths);
