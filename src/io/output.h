#pragma once

#include <GraphBLAS.h>

#include "io/cli.h"

/**
 * @brief Writes all non-zero entries of a boolean matrix as vertex pairs to a file.
 *
 * @param[out] out     Open file stream for writing. Must not be `NULL`.
 * @param[in]  matrix  Input `GrB_Matrix` (`GrB_BOOL`) to extract pairs from.
 * @param[in]  nvals   Number of non-zero entries in `matrix` (from
 * `GrB_Matrix_nvals`).
 *
 * @return `GrB_SUCCESS` on success, or `GrB_OUT_OF_MEMORY` if allocation fails.
 */
GrB_Info write_matrix_pairs(FILE *out, GrB_Matrix matrix, GrB_Index nvals);

/**
 * @brief Constructs the output file path for a benchmark result based on CLI
 * arguments.
 *
 * Example:
 * ```
 * graph_file_path = "/path/to/example.g"
 * output_path     = NULL (use default)
 * Result: "output/example.out"
 * ```
 *
 * @param[in]  args        Parsed CLI arguments containing input/output paths.
 * @param[out] output_file Buffer to receive the constructed path. Must have at
 *                         least `size` bytes available.
 * @param[in]  size        Size of the `output_file` buffer in bytes.
 */
void resolve_output_path(const Args *args, char *output_file, size_t size);
