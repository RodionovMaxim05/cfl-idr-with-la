#pragma once

#include <LAGraphX.h>
#include <stdbool.h>

#include "cfl_idr.h"

#define NT_START 0

/**
 * @brief Basic type of context-free grammar (transition diagram).
 */
typedef enum {
	MR_GRAMMAR_ALPHA = 0, // Standard alpha grammar (direct transitions).
	MR_GRAMMAR_BETA		  // Beta grammar (inverted logic for paired groups).
} MRGrammarKind;

/**
 * @brief Configuration of specific grammar generation parameters.
 *
 * Encapsulates the rule generation type and additional context data
 * required for grammar modifications.
 */
typedef struct {
	MRGrammarKind kind;	   // Basic type of transition diagram (ALPHA or BETA).
	int64_t exclude_index; // The index of the parenthesis to be isolated. If < 0
						   // (e.g., -1), then there is no exception.
} MRGrammarConfig;

/**
 * @brief Grammar specification for CFL-reachability analysis in mutual refinement.
 *
 * Encodes a weighted context-free grammar in EWCNF format suitable for
 * `LAGraph_CFL_AllPaths`. Contains mode configuration, group size, counts for
 * nonterminals, terminals, and rules, plus a pointer to the rule array.
 */
typedef struct {
	MRGrammarConfig config;	   // Configuration of specific grammar parameters
	int64_t k;				   // Group size for terminal matrices
	int64_t nonterms_count;	   // Number of nonterminal symbols
	int64_t terms_count;	   // Number of terminal symbols (edge labels)
	int64_t rules_count;	   // Number of production rules
	LAGraph_rule_EWCNF *rules; // Array of grammar rules in EWCNF format
} MRGrammar_t;

/**
 * @brief Creates an basic grammar for Dyck language analysis.
 *
 * Basic Dyck grammar is used to analyze Dyck languages.
 *
 * @param[in] n_par       Number of parenthesis types (pairs of `(`/`)`).
 * @param[in] n_bra       Number of bracket types (pairs of `[`/`]`).
 * @param[in] has_normal  Whether to include an epsilon/unlabeled edge terminal.
 *
 * @return A fully initialized `MRGrammar_t` accepting the Dyck language over
 *         the specified alphabet, or a zeroed struct on allocation failure.
 */
MRGrammar_t get_dyck_grammar(int64_t n_par, int64_t n_bra, bool has_normal);

/**
 * @brief Creates an project grammar for Dyck language analysis.
 *
 * Project grammar is used to analyze Dyck languages while projecting
 * specific bracket identifiers.
 *
 * @param[in] n_par       Number of parenthesis types (pairs of `(`/`)`).
 * @param[in] n_bra       Number of bracket types (pairs of `[`/`]`).
 * @param[in] has_normal  Whether to include an epsilon/unlabeled edge terminal.
 *
 * @return A fully initialized `MRGrammar_t` accepting unlabeled Dyck sequences,
 *         or a zeroed struct on allocation failure.
 */
MRGrammar_t get_project_grammar(int64_t n_par, int64_t n_bra, bool has_normal);

/**
 * @brief Creates an exclusion grammar for Dyck language analysis.
 *
 * Exclusion grammars are used to analyze Dyck languages while excluding
 * specific bracket identifiers from consideration.
 *
 * @param[in] n_par       Number of parenthesis types (pairs of `(`/`)`).
 * @param[in] n_bra       Number of bracket types (pairs of `[`/`]`).
 * @param[in] has_normal  Whether to include an epsilon/unlabeled edge terminal.
 * @param[in] ex_bra      Index of the bracket type to exclude from parity tracking
 *                        (0 ≤ ex_bra < n_bra).
 *
 * @return A fully initialized Alpha `MRGrammar_t` configured for the specified
 *         analysis, or a zeroed struct on allocation failure.
 */
MRGrammar_t get_exclude_grammar(int64_t n_par, int64_t n_bra, bool has_normal,
								int64_t ex_bra);

/**
 * @brief Creates an Alpha grammar based on the specified grammar type and
 * parameters.
 *
 * This function serves as a factory for Alpha grammars used in Dyck language
 * analysis. It selects the appropriate grammar generator based on the `curGrammar`
 * parameter, which determines the specific language class to analyze.
 *
 * The grammar type keys correspond to the article designations as follows:
 * - `IDR_PARITY`  → PAR (`D_par(Σ_α, Σ_β)`, k=1)
 * - `IDR_PARITY2` → PAR2 (`D_par²(Σ_α, Σ_β)`, k=2)
 * - `IDR_PARITYK` → PARk (`D_parᵏ(Σ_α, Σ_β)`, parameterized)
 * - `IDR_SE`      → PAR2E (`D⁺_par²(Σ_α, Σ_β)`, k=2 + valid endpoints)
 * - `IDR_PROJECT` → PARUnl — uses PAR grammar; the unlabeled-graph transformation
 *                   is applied externally in `mutualRefinement` via `getProjectSppf`
 * - `IDR_EXCLUDE` → PARErase — uses PAR grammar; label-erasing iteration
 *                   is applied externally in `mutualRefinement` via `getExcludeSppf`
 * - `IDR_ALL`     → COM (PAR2E grammar)
 * - `IDR_DEFAULT` → Standard mutual-refinement base grammar
 *
 * Note: for `IDR_PROJECT`, `IDR_EXCLUDE` and `IDR_ALL`, the grammar itself is the
 * same as a simpler variant — the additional algorithmic logic that defines these
 * methods is handled in `mutualRefinement`, not here.
 *
 * @param[in] grammar_type  Grammar variant selector (see `IdrGrammarType`).
 * @param[in] n_par         Number of parenthesis types (pairs of `(`/`)`).
 * @param[in] n_bra         Number of bracket types (pairs of `[`/`]`).
 * @param[in] has_normal    Whether to include an epsilon/unlabeled edge terminal.
 *
 * @return A fully initialized Alpha `MRGrammar_t` configured for the specified
 *         analysis, or a zeroed struct on allocation failure.
 */
MRGrammar_t get_alpha_grammar(IdrGrammarType grammar_type, int64_t n_par,
							  int64_t n_bra, bool has_normal);

/**
 * Creates a Beta grammar based on the specified grammar type and parameters.
 *
 * This function serves as a factory for Beta grammars used in Dyck language
 * analysis. Beta grammars typically complement Alpha grammars in the analysis
 * pipeline.
 *
 * The grammar type keys correspond to the article designations as follows:
 * - `IDR_PARITY`  → PAR (`D_par(Σ_β, Σ_α)`, k=1)
 * - `IDR_PARITY2` → PAR2 (`D_par²(Σ_β, Σ_α)`, k=2)
 * - `IDR_PARITYK` → PARk (`D_parᵏ(Σ_β, Σ_α)`, parameterized)
 * - `IDR_SE`      → PAR2E (`D⁺_par²(Σ_β, Σ_α)`, k=2 + valid endpoints)
 * - `IDR_PROJECT` → PARUnl — uses PAR grammar; the unlabeled-graph transformation
 *                   is applied externally in `mutualRefinement` via `getProjectSppf`
 * - `IDR_EXCLUDE` → PARErase — uses PAR grammar; label-erasing iteration
 *                   is applied externally in `mutualRefinement` via `getExcludeSppf`
 * - `IDR_ALL`     → COM (PAR2E grammar)
 * - `IDR_DEFAULT` → Standard mutual-refinement base grammar
 *
 * Note: for `IDR_PROJECT`, `IDR_EXCLUDE`, and `IDR_ALL`, the grammar itself is the
 * same as a simpler variant — the additional algorithmic logic that defines these
 * methods is handled in [mutualRefinement], not here.
 *
 * @param[in] grammar_type  Grammar variant selector (see `IdrGrammarType`).
 * @param[in] n_par         Number of parenthesis types (pairs of `(`/`)`).
 * @param[in] n_bra         Number of bracket types (pairs of `[`/`]`).
 * @param[in] has_normal    Whether to include an epsilon/unlabeled edge terminal.
 *
 * @return A fully initialized Beta `MRGrammar_t` configured for the specified
 *         analysis, or a zeroed struct on allocation failure.
 */
MRGrammar_t get_beta_grammar(IdrGrammarType grammar_type, int64_t n_par,
							 int64_t n_bra, bool has_normal);

/**
 * @brief Frees all resources owned by an `MRGrammar_t`.
 *
 * Deallocates the `rules` array via `free()` and zero-initializes the structure.
 *
 * @param[in] gr  Grammar to free. May be `NULL` (no-op).
 */
void grammar_free(MRGrammar_t *gr);

/**
 * @brief Calculates the size of a specific bracket group under the given grammar
 * configuration.
 *
 * @param[in] config  Pointer to the grammar configuration.
 * @param[in] b       Group index (0 to k).
 * @param[in] n       Total number of bracket matrices to distribute.
 * @param[in] k       The base number of target groups.
 *
 * @return The number of elements allocated to group b.
 */
int64_t MR_grammar_get_group_size(const MRGrammarConfig *config, int64_t b,
								  int64_t n, int64_t k);

/**
 * @brief Determines the target group bucket and sub-index offset for the i-th
 * bracket element.
 *
 * @param[in]  config       Pointer to the grammar configuration.
 * @param[in]  i            The original sequential index of the bracket element.
 * @param[in]  n            Total number of bracket matrices.
 * @param[in]  k            The base number of target groups.
 * @param[out] out_group    Pointer to receive the assigned target group index.
 * @param[out] out_sub_idx  Pointer to receive the internal zero-based offset within
 *                          the group.
 */
void MR_grammar_get_bracket_layout(const MRGrammarConfig *config, int64_t i,
								   int64_t n, int64_t k, int64_t *out_group,
								   int64_t *out_sub_idx);
