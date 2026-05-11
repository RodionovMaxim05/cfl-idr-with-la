#pragma once

#include <GraphBLAS.h>

#include "cfl_idr.h"
#include "parser.h"
#include "symbol_list.h"
#include "terminal/terminal_format.h"

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
