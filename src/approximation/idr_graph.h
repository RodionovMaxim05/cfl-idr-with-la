#pragma once

#include <GraphBLAS.h>
#include <stdbool.h>
#include <stdint.h>

#include "cfl_idr.h"

/**
 * @brief Computes the total number of terminal symbols (edge labels) in the graph.
 *
 * @param[in] n_par   Number of parenthesis types.
 * @param[in] n_bra   Number of bracket types.
 * @param[in] normal  Pointer to the normal (epsilon) matrix, or `NULL` if absent.
 *
 * @return Total number of terminals: `2 * n_par + 2 * n_bra + (normal ? 1 : 0)`.
 */
static inline int64_t get_terms_count(int64_t n_par, int64_t n_bra,
									  GrB_Matrix normal) {
	return 2 * n_par + 2 * n_bra + (normal != NULL ? 1 : 0);
}

/**
 * @brief Prints a human-readable debug representation of an `IdrGraph`.
 *
 * Outputs graph metadata (vertex count, parenthesis/bracket counts) followed by
 * a summary of each adjacency matrix using `GxB_print` with verbosity level 2.
 * If the graph pointer is `NULL`, prints a placeholder message.
 *
 * @param[in] g  Graph to print.
 */
void idr_graph_print(const IdrGraph *g);

/**
 * @brief Counts the total number of edges across all labeled and unlabeled matrices.
 *
 * Iterates over every matrix stored in the graph and sums their non-zero entry
 * counts (`nvals`). The result represents the total edge count in the multi-labeled
 * graph representation.
 *
 * @param[in] graph  Input graph. Must not be `NULL`.
 *
 * @return Total number of edges (non-zero entries) in the graph.
 */
GrB_Index idr_graph_count_edges(const IdrGraph *graph);

/**
 * @brief Constructs an `IdrGraph` from an array of input adjacency matrices.
 *
 * The input array `input_matrices` must contain matrices in the following order:
 * 1. Parenthesis pairs: `[open_par[0], close_par[0], ..., open_par[n_par-1],
 * close_par[n_par-1]]`
 * 2. Bracket pairs:     `[open_bra[0], close_bra[0], ..., open_bra[n_bra-1],
 * close_bra[n_bra-1]]`
 * 3. Normal matrix:     `[normal]` (only if `has_normal == true`)
 *
 * When `filter_empty` is `true`, parenthesis/bracket types where either the opening
 * or closing matrix is empty (`nvals == 0`) are excluded from the output graph, and
 * their matrices are freed via `GrB_Matrix_free`. When `false`, all matrices are
 * adopted without filtering.
 *
 * @note The function adopts ownership of the input matrices: successful calls
 * transfer matrix handles to the output graph. On error or when filtering removes a
 * pair, the corresponding matrices are freed.
 *
 * @param[out] out            Output graph structure to populate. Must point to valid
 *                            memory.
 * @param[in]  input_matrices Array of input `GrB_Matrix` objects in the order
 *                            described above.
 * @param[in]  n_par          Number of parenthesis types in the input.
 * @param[in]  n_bra          Number of bracket types in the input.
 * @param[in]  has_normal     `true` if the input array includes a normal (epsilon)
 *                            matrix.
 * @param[in]  n              Number of vertices in the graph.
 * @param[in]  filter_empty   If `true`, exclude empty parenthesis/bracket pairs.
 *
 * @return `GrB_SUCCESS` on success, or `GrB_OUT_OF_MEMORY` if allocation fails.
 */
GrB_Info build_idr_graph(IdrGraph *out, GrB_Matrix *input_matrices, int64_t n_par,
						 int64_t n_bra, bool has_normal, GrB_Index n,
						 bool filter_empty);

/**
 * @brief Extracts all adjacency matrices from an `IdrGraph` into a flat array.
 *
 * Allocates a new array of `GrB_Matrix` pointers containing references to the
 * graph's internal matrices in the canonical order:
 * - `open_par[0]`, `close_par[0]`, ..., `open_par[n_par-1]`, `close_par[n_par-1]`
 * - `open_bra[0]`, `close_bra[0]`, ..., `open_bra[n_bra-1]`, `close_bra[n_bra-1]`
 * - `normal` (if present)
 *
 * The caller is responsible for freeing the returned array with `free()`, but
 * must not free the individual matrices, as they remain owned by the `IdrGraph`.
 *
 * @param[out] out    Pointer to receive the allocated array of matrices.
 * @param[in]  graph  Source graph. Must not be `NULL`.
 *
 * @return `GrB_SUCCESS` on success, or `GrB_OUT_OF_MEMORY` if allocation fails.
 */
GrB_Info idr_graph_collect_matrices(GrB_Matrix **out, const IdrGraph *graph);

/**
 * @brief Constructs a unified boolean adjacency matrix from an `IdrGraph`.
 *
 * Combines all labeled edge matrices (parentheses, brackets, and optional `normal`)
 * into a single unweighted adjacency matrix via element-wise logical OR.
 *
 * @param[out] out    Newly allocated `GrB_Matrix` (n × n, `GrB_BOOL`) holding the
 *                    unified adjacency.
 * @param[in]  graph  Input `IdrGraph`. Must not be `NULL`.
 *
 * @return `GrB_SUCCESS` on success, or a GraphBLAS error code on failure.
 */
GrB_Info idr_graph_to_adjacency(GrB_Matrix *out, const IdrGraph *graph);
