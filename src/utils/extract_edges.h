#pragma once

#include "grammar/grammar_analysis_utils.h"

/**
 * @brief Specifies a source-target vertex pair for on-demand CFL-reachability
 * analysis.
 *
 * Used to constrain path extraction to a single query pair `(src, tgt)` instead of
 * computing all-pairs results. This enables efficient incremental analysis when
 * only specific reachability questions need to be answered.
 */
typedef struct {
	GrB_Index src; // Source vertex index
	GrB_Index tgt; // Target vertex index
} TargetPath;

/**
 * @brief Extracts terminal-labeled edges from CFL-reachability output matrices.
 *
 * This function performs grammar-guided backtracking to recover the concrete edge
 * labels that witness CFL-reachable paths. Given:
 * - `paths[A]`: Matrices indicating which vertex pairs are reachable via nonterminal
 * `A`
 * - `adj_matrices`: Input graph adjacency matrices where terminal `t` is located at
 * index `grammar.nonterms_count + t`
 * - `grammar`: CFG rules in WCNF format defining the Dyck-language constraints
 *
 * The algorithm:
 * 1. Maps interleaved grammar terminal indices to straight output indices based on
 *   layout parity.
 * 2. Builds fast lookup indices:
 *    - Terminal grammar rules index (`A → a`).
 *    - Grouped binary rules index (`BinaryIndex` for `A → B C`).
 *    - Dense row presence bitset matrix for O(1) checks of nonterminal `B` in row
 *      `i`.
 * 3. Initializes an explicit DFS work stack with:
 *    - A single `(src, tgt, NT_START)` tuple for on-demand analysis, or
 *    - All non-reflexive pairs from `paths[NT_START]` for all-pairs extraction.
 * 4. Performs DFS over the derivation forest:
 *    - For terminal rules `A → t`: checks if edge `(i, j)` exists in
 *      `adj_matrices[nonterms_count + t]` and sets `out[straight_idx](i, j) = true`.
 *    - For binary rules `A → B C`:checks subpaths `(i, mid)` via `B` and `(mid, j)`
 *      via `C`, then pushes unvisited subproblems onto the stack.
 *
 * @param[out] out           Array of `grammar.terms_count` output matrices
 *                           (`GrB_BOOL`, n × n). Each `out[t]` receives edges
 *                           labeled with terminal `t`. Must be pre-allocated by
 *                           caller or set to `NULL` pointers.
 * @param[in]  paths         Array of `grammar.nonterms_count` CFL-reachability
 *                           result matrices from `LAGraph_CFL_AllPaths`.
 * @param[in]  adj_matrices  Input adjacency matrices for each terminal symbol.
 * @param[in]  grammar       Grammar specification.
 * @param[in]  n             Number of vertices (matrix dimension).
 * @param[in]  n_par         Number of parenthesis types.
 * @param[in]  n_bra         Number of bracket types.
 * @param[in]  config        Grammar configuration specifying layout strategy and
 *                           exclusions.
 * @param[in]  target_path   Optional query pair for on-demand analysis. Pass `NULL`
 *                           for all-pairs extraction.
 *
 * @return `GrB_SUCCESS` on success, or a GraphBLAS error code on failure.
 */
GrB_Info extract_edges_from_outputs(GrB_Matrix *out, GrB_Matrix *paths,
									GrB_Matrix *adj_matrices, MRGrammar_t grammar,
									GrB_Index n, int64_t n_par, int64_t n_bra,
									const MRGrammarConfig *config,
									const TargetPath *target_path);
