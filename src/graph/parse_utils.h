#pragma once

#include <GraphBLAS.h>

#include "cfl_idr.h"

/**
 * @brief Parses a graph definition file and constructs an `IdrGraph` for
 * CFL-reachability analysis.
 *
 * This is the primary high-level entry point for loading graphs from disk. It
 * performs the complete pipeline:
 * 1. **File I/O**: Opens and reads the graph file in the internal text format.
 * 2. **Parsing**: Invokes `process_graph` to parse edges, labels, and symbol
 * metadata.
 * 3. **Matrix conversion**: Transforms parsed edges into GraphBLAS adjacency
 * matrices via `get_grb_matrices_from_graph`.
 * 4. **Graph assembly**: Groups matrices by bracket type and builds the structured
 *    `IdrGraph` via `get_idr_graph`, using `DefaultTerminalFormat` for label
 * classification.
 *
 * @param[out] out       Output `IdrGraph` structure to populate. Must be valid and
 *                       uninitialized. On success, contains the loaded graph.
 * @param[in]  filename  Path to the graph definition file. Must be a readable text
 *                       file in the expected format.
 *
 * @return `GrB_SUCCESS` on success, or a GraphBLAS error code on failure.
 */
GrB_Info parse_graph(IdrGraph *out, const char *filename);
