#pragma once

/**
 * @file cfl_idr_graph_builder.h
 * @brief Utilities for constructing `IdrGraph` instances from parsed graph matrices.
 */

#include <stdbool.h>
#include <stddef.h>

#include "cfl_idr.h"
#include "parser.h"
#include "symbol_list.h"

/**
 * @brief Classification of bracket/parenthesis types for CFL grammar construction.
 */
typedef enum {
	BRACKET_TYPE_PARENTHESES, ///< Labels representing parentheses: `(` and `)`
	BRACKET_TYPE_BRACKETS,	  ///< Labels representing brackets: `[` and `]`
	BRACKET_TYPE_UNKNOWN	  ///< Unrecognized or non-bracket label
} BracketType;

/**
 * @brief Interface for parsing terminal labels in CFL-reachability graphs.
 *
 * Provides a pluggable strategy to classify edge labels as opening/closing
 * parentheses or brackets. This abstraction allows the same graph parsing
 * pipeline to support different labeling conventions (e.g., `op_0`/`cp_0`
 * vs `(`/`)`, `[`/`]`).
 *
 * A `TerminalFormat` instance must implement both function pointers to be
 * usable with `get_idr_graph` and related parsing utilities.
 */
typedef struct {
	/**
	 * @brief Determines the bracket category of a given label.
	 *
	 * @param[in] label  Null-terminated string representing the edge label.
	 * @return `BRACKET_TYPE_PARENTHESES`, `BRACKET_TYPE_BRACKETS`, or
	 *         `BRACKET_TYPE_UNKNOWN`.
	 */
	BracketType (*get_type)(const char *label);

	/**
	 * @brief Checks whether a label represents an opening or closing symbol.
	 *
	 * @param[in] label  Null-terminated string representing the edge label.
	 * @return `1` if opening, `0` if closing, `-1` if unrecognized.
	 */
	int (*is_open)(const char *label);
} TerminalFormat;

/**
 * @brief Default terminal format implementation.
 *
 * Uses a simple character-based convention:
 * - First character: `'o'` = opening, `'c'` = closing
 * - Second character: `'p'` = parentheses, `'b'` = brackets
 * - Example labels: `"op_0"`, `"cp_1"`, `"ob_2"`, `"cb_3"`
 *
 * This format is used by default in `get_idr_graph` when no custom format
 * is provided.
 */
extern const TerminalFormat DefaultTerminalFormat;

/**
 * @brief Constructs an `IdrGraph` from parsed graph matrices and symbol metadata.
 *
 * This function transforms a flat list of labeled adjacency matrices into the
 * structured `IdrGraph` representation required for CFL-reachability analysis.
 *
 * It performs the following steps:
 * 1. **Classification**: Uses `TerminalFormat` to classify each matrix label as:
 *    - Opening/closing parenthesis (`(`/`)`)
 *    - Opening/closing bracket (`[`/`]`)
 *    - Special `"normal"` label for unlabeled (epsilon) edges
 *    - Unknown (ignored)
 *
 * 2. **Grouping**: Groups matrices by `block_index` to pair opening and closing
 *    matrices belonging to the same bracket type.
 *
 * 3. **Filtering**: Excludes:
 *    - Empty matrices (`nvals == 0`)
 *    - Incomplete pairs (types with only opening OR only closing matrix)
 *
 * 4. **Assembly**: Allocates and populates the `IdrGraph` structure with valid
 *    parenthesis/bracket pairs and the optional `normal` matrix.
 *
 * Matrix ownership is transferred from `gm->matrices` to the output `IdrGraph`:
 * successfully adopted matrices are set to `NULL` in the input array to prevent
 * double-free.
 *
 * @param[out] out           Output `IdrGraph` structure to populate. Must be valid
 *                           and uninitialized.
 * @param[in]  gm            Parsed graph matrices with associated symbol metadata.
 *                           Must not be `NULL`.
 * @param[in]  symbol_list   Symbol name lookup table for resolving label strings.
 *                           Must not be `NULL`.
 * @param[in]  n             Number of vertices in the graph (matrix dimension).
 * @param[in]  fmt           Terminal format parser for classifying bracket labels.
 *                           Must not be `NULL`.
 *
 * @return `GrB_SUCCESS` on success, or a GraphBLAS error code on failure.
 */
GrB_Info get_idr_graph(IdrGraph *out, GraphMatrices *gm,
					   const SymbolList *symbol_list, GrB_Index n,
					   const TerminalFormat *fmt);
