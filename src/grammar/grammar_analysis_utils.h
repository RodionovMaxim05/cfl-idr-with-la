#pragma once

#include <LAGraphX.h>
#include <stdbool.h>

#include "cfl_idr.h"

#define NT_START 0

/**
 * @brief Grammar specification for CFL-reachability analysis in mutual refinement.
 *
 * Encodes a weighted context-free grammar in WCNF format suitable for
 * `LAGraph_CFL_AllPaths`. Contains counts for nonterminals, terminals, and rules,
 * plus a pointer to the rule array.
 */
typedef struct {
	int64_t nonterms_count;	  // Number of nonterminal symbols
	int64_t terms_count;	  // Number of terminal symbols (edge labels)
	int64_t rules_count;	  // Number of production rules
	LAGraph_rule_WCNF *rules; // Array of grammar rules in WCNF format
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
