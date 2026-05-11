#pragma once

#include "grammar/grammar.h"

/**
 * @brief Default output directory for result files.
 *
 * Used when the `-o` option is not specified on the command line.
 */
#define DEFAULT_OUTPUT_DIR "output"

/**
 * @brief Parsed command-line arguments for the CFL-reachability tool.
 *
 * This structure holds all configuration options after argument parsing,
 * providing a unified interface for the main program logic.
 *
 * @note String fields (`graph_file_path`, `grammar_str`, `output_path`) point
 *       to argv[] storage and must not be freed. The `output_path` field may
 *       be `NULL` to indicate the default directory should be used.
 */
typedef struct {
	const char *graph_file_path; // Path to input graph file (required)
	const char *grammar_str;	 // Original grammar name string from CLI
	const char *output_path;	 // Optional output directory, or `NULL` for default
	IdrGrammarType grammar_type; // Parsed grammar type enum
	int quiet;					 // Non-zero to suppress pair output (flag: `-q`)
	int valueflow; // Non-zero to enable value-flow constraints (flag: `-valueflow`)
	int on_demand; // Non-zero for on-demand two-phase refinement
	int parity_d;  // Non-zero for parity-only on-demand mode (PARD)
} Args;

/**
 * @brief Prints usage information and available options to stderr.
 *
 * @param[in] prog  Program name (typically `argv[0]`). Must not be `NULL`.
 */
void print_usage(const char *prog);

/**
 * @brief Parses command-line arguments into a structured `Args` configuration.
 *
 * Supported syntax:
 * ```
 * program [options] <input.g> <grammar>
 * ```
 *
 * Options:
 * - `-h`: Print help and exit immediately.
 * - `-o <path>`: Set custom output directory (default: `DEFAULT_OUTPUT_DIR`).
 * - `-q`: Enable quiet mode (suppress detailed pair output).
 * - `-valueflow`: Enable value-flow specific constraints and optimizations.
 *
 * Grammar types:
 * - `parity`: (PAR) Parity grammar with k=1 - Parity condition
 * - `parity2`: (PAR2) Parity grammar with k=2 - Extended parity condition
 * - `se`: (PAR2E) Structured equality grammar - Valid endpoints
 * - `project`: (PARUnl) Projection grammar - Projection to an unlabeled Dyck grammar
 * - `exclude`: (PARErase) Exclusion grammar - Erasing labels
 * - `all`: (COM) Comprehensive grammar
 * - `parityD`: (PARD) On-demand parity grammar
 * - `on-demand`: (COMD) On-demand combined method
 *
 * Special handling:
 * - `parityD` and `on-demand` are meta-modes that adjust both `grammar_type`
 *   and the corresponding mode flags (`parity_d` or `on_demand`), while
 *   preserving the base grammar string for internal use.
 *
 * @param[out] out   Structure to populate with parsed arguments. Must not be `NULL`.
 * @param[in]  argc  Argument count from `main()`.
 * @param[in]  argv  Argument vector from `main()`.
 *
 * @return Non-zero (true) on successful parsing, zero (false) on error.
 *         On error, a descriptive message is printed to stderr.
 */
int parse_args(Args *out, int argc, char *argv[]);
