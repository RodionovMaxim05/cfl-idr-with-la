#pragma once

#include "grammar.h"

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
