#pragma once

#include <GraphBLAS.h>
#include <stdbool.h>
#include <stdint.h>

/**
 * @brief Sparse-matrix representation of a Dyck-reachability graph.
 *
 * The graph encodes parenthesis-labeled edges as separate GraphBLAS matrices.
 * Each pair of opening / closing matrices corresponds to one parenthesis *type*
 * (indexed 0 … n_par-1 for parentheses, 0 … n_bra-1 for brackets).
 * Unlabeled edges are stored in `normal`.
 */
typedef struct {
	GrB_Index n;		   // Number of vertices in the graph
	int64_t n_par;		   // Number of parenthesis types ( and )
	GrB_Matrix *open_par;  // open_par[i]  - edges labeled `(i`
	GrB_Matrix *close_par; // close_par[i] - edges labeled `)i`
	int64_t n_bra;		   // Number of bracket types [ and ]
	GrB_Matrix *open_bra;  // open_bra[i]  - edges labeled `[i`
	GrB_Matrix *close_bra; // close_bra[i] - edges labeled `]i`
	GrB_Matrix normal;	   // Unlabeled (epsilon) edges
} IdrGraph;

/**
 * @brief Grammar types used to control the over-approximation algorithm.
 *
 * Each value selects a different Dyck-CFL grammar, which determines the
 * trade-off between precision and computational cost.
 */
typedef enum {
	IDR_UNKNOWN = -1, // Uninitialized / error sentinel
	IDR_DEFAULT = 0,  // Default mutual refinement grammar
	IDR_PARITY,		  // PAR: Parity grammar with k=1 - Parity condition
	IDR_PARITY2,	  // PAR2: Parity grammar with k=2 - Extended parity condition
	IDR_SE,			  // PAR2E: Structured equality grammar - Valid endpoints
	IDR_PROJECT,	  // PARUnl: Projection grammar - Projection to an unlabeled Dyck
					  // grammar
	IDR_EXCLUDE,	  // PARErase: Exclusion grammar - Erasing labels
	IDR_ALL			  // COM: Comprehensive grammar
} IdrGrammarType;

/**
 * @brief Releases all GraphBLAS matrices contained in an `IdrGraph`.
 *
 * Calls `GrB_free` on every matrix stored in @p g (including all elements of the
 * `open_par`, `close_par`, `open_bra`, and `close_bra` arrays and the `normal`
 * matrix).
 *
 * @param[in] g  Graph whose matrices are to be freed.
 */
void idr_graph_free(IdrGraph *g);

/**
 * @brief Computes an under-approximation f reachable paths in a graph using Dyck
 * language analysis.
 *
 * This function analyzes the input graph to identify paths that are guaranteed to be
 * present (a conservative under-approximation).
 *
 * @param[out] result     On success, a newly allocated `GrB_Matrix` (n × n,
 *                        `GrB_BOOL`) holding the under-approximation.
 * @param[in]  graph      Input graph to analyze.  Must not be NULL.
 * @param[in]  valueflow  If `true`, apply value-flow specific constraints.
 *
 * @return `GrB_SUCCESS` on success, or a GraphBLAS error code otherwise.
 */
GrB_Info idr_get_under_approx(GrB_Matrix *result, const IdrGraph *graph,
							  bool valueflow);

/**
 * @brief Computes a refined over-approximation of reachable paths using mutual
 * refinement (MR) algorithm.
 *
 * This function implements the mutual refinement algorithm to compute a precise
 * over-approximation of reachable paths. It combines graph condensation with
 * iterative grammar-based analysis to refine path approximations.
 *
 * @param[out] result        On success, a newly allocated `GrB_Matrix` (n × n,
 *                           `GrB_BOOL`).
 * @param[in]  graph         Input graph to analyze.  Must not be NULL.
 * @param[in]  grammar_type  Grammar used for refinement (see `IdrGrammarType`).
 * @param[in]  under_approx  Optional under-approximation matrix produced by
 *                           `idr_get_under_approx`.  Pass `NULL` to skip
 *                           condensation.
 * @param[in]  valueflow     If `true`, apply value-flow specific optimizations
 *                           and constraints.
 * @param[in]  filter_empty  If `true`, remove bracket/square bracket types
 *                           where one of the matrices (opening/closing matrix)
 *                           is empty (nvals == 0 for one of the matrix matrices).
 *
 * @return `GrB_SUCCESS` on success, or a GraphBLAS error code otherwise.
 */
GrB_Info idr_get_over_approx(GrB_Matrix *result, const IdrGraph *graph,
							 IdrGrammarType grammar_type, GrB_Matrix under_approx,
							 bool valueflow, bool filter_empty);

/**
 * @brief Filters vertices from the graph that cannot participate in value-flow
 * analysis paths.
 *
 * This function implements a value-flow specific optimization that removes vertices
 * which cannot be part of any valid store-to-load path. In value-flow analysis:
 * - **Store operations** are represented by opening brackets `[i`
 * - **Load operations** are represented by closing brackets `]i`
 *
 * A vertex is retained only if there exists:
 * 1. A reachable path from some store operation to the vertex
 * 2. A reachable path from the vertex to some load operation
 *
 * @param[out] out    Output graph with unreachable vertices removed.
 * @param[in]  graph  Input graph to filter.  Must not be NULL.
 *
 * @return `GrB_SUCCESS` on success, or a GraphBLAS error code otherwise.
 */
GrB_Info idr_remove_valueflow_unreachable(IdrGraph *out, const IdrGraph *graph);

/**
 * @brief Performs on-demand mutual refinement using two-phase grammar analysis.
 *
 * This approach is more efficient than full mutual refinement when only incremental
 * updates are needed.
 *
 * @param[out] result        On success, a newly allocated `GrB_Matrix`
 *                           (n × n, `GrB_BOOL`) with the refined result.
 * @param[in]  graph         Input graph to analyze.  Must not be NULL.
 * @param[in]  under_approx  Current under-approximation (e.g., from
 *                           `idr_get_under_approx`).  Must not be NULL.
 * @param[in]  over_approx   Current over-approximation (e.g., from
 *                           `idr_get_over_approx`).  Must not be NULL.
 * @param[in]  parity_d      If `true`, execute only the first refinement phase
 *                           (PARD).  If `false`, execute both phases (COMD).
 * @param[in]  valueflow     If `true`, apply value-flow specific optimizations
 *                           and constraints.
 * @param[in]  filter_empty  If `true`, remove bracket/square bracket types where
 *                           one of the matrices (opening/closing matrix) is empty
 *                           (nvals == 0 for one of the matrix matrices).
 *
 * @return `GrB_SUCCESS` on success, or a GraphBLAS error code otherwise.
 */
GrB_Info idr_get_on_demand(GrB_Matrix *result, const IdrGraph *graph,
						   GrB_Matrix under_approx, GrB_Matrix over_approx,
						   bool parity_d, bool valueflow, bool filter_empty);
