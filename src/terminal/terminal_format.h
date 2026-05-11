#pragma once

#include <stdbool.h>
#include <stddef.h>

/**
 * @brief Classification of bracket/parenthesis types for CFL grammar construction.
 */
typedef enum {
	BRACKET_TYPE_PARENTHESES, // Labels representing parentheses: `(` and `)`
	BRACKET_TYPE_BRACKETS,	  // Labels representing brackets: `[` and `]`
	BRACKET_TYPE_UNKNOWN	  // Unrecognized or non-bracket label
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
