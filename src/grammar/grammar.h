#pragma once

#include <GraphBLAS.h>
#include <LAGraph.h>
#include <LAGraphX.h>
#include <stdbool.h>
#include <stdlib.h>

#include "cfl_idr.h"
#include "grammar_analysis_utils.h"

/**
 * Represents the states for the "valid endpoints" condition in structured equality
 * grammars.
 *
 * The automaton tracks whether brackets have been properly opened and closed:
 * - [QE]: Start/end state (q0) - no unmatched brackets
 * - [QO]: Open state (q1) - an opening bracket has been seen
 * - [QC]: Close state (qf) - a closing bracket has been seen after opening
 */
#define RSTATE_QE 0 // Start/end state: balanced or empty bracket sequence
#define RSTATE_QO 1 // Open state: unmatched opening bracket seen
#define RSTATE_QC 2 // Close state: closing bracket seen after opening

#define RSTATE_COUNT (int64_t)3

/**
 * @brief Allocates and initializes a new `MRGrammar_t` structure.
 *
 * @param[in] rules_count     Number of production rules to allocate space for.
 * @param[in] terms_count     Number of terminal symbols in the grammar.
 * @param[in] nonterms_count  Number of nonterminal symbols in the grammar.
 *
 * @return Initialized `MRGrammar_t` with allocated `rules` array, or a zeroed
 *         struct if allocation fails.
 */
static inline MRGrammar_t make_grammar(int64_t rules_count, int64_t terms_count,
									   int64_t nonterms_count) {
	MRGrammar_t gr = {0};
	gr.rules = malloc(rules_count * sizeof(LAGraph_rule_WCNF));
	if (!gr.rules) {
		return gr;
	}
	gr.rules_count = rules_count;
	gr.terms_count = terms_count;
	gr.nonterms_count = nonterms_count;
	return gr;
}

/**
 * Creates a basic Dyck grammar that accepts properly nested parenthesis and bracket
 * sequences.
 *
 * This grammar corresponds to the standard Dyck language without interleaving
 * constraints. It generates strings where parentheses and brackets are properly
 * matched and nested.
 *
 * @param[in] n_par       Number of parenthesis types (pairs of `(`/`)`).
 * @param[in] n_bra       Number of bracket types (pairs of `[`/`]`).
 * @param[in] has_normal  Whether to include an epsilon/unlabeled edge terminal.
 *
 * @return A fully initialized `MRGrammar_t` accepting the Dyck language over
 *         the specified alphabet, or a zeroed struct on allocation failure.
 */
MRGrammar_t dyck_grammar(int64_t n_par, int64_t n_bra, bool has_normal);

/**
 * Creates a projection grammar that abstracts away specific bracket labels.
 *
 * Corresponds to the unlabeled Dyck grammar `D(Σ_u)` from the article (Sect. 4.2),
 * used in the **PARUnl** method.
 *
 * This grammar projects all opening symbols to a single abstract '(' and all closing
 * symbols to a single abstract ')', then accepts properly nested sequences. It's
 * used for unlabeled Dyck language analysis.
 *
 * @param[in] n_par       Number of parenthesis types (pairs of `(`/`)`).
 * @param[in] n_bra       Number of bracket types (pairs of `[`/`]`).
 * @param[in] has_normal  Whether to include an epsilon/unlabeled edge terminal.
 *
 * @return A fully initialized `MRGrammar_t` accepting unlabeled Dyck sequences,
 *         or a zeroed struct on allocation failure.
 */
MRGrammar_t dyck_project_grammar(int64_t n_par, int64_t n_bra, bool has_normal);

/**
 * Creates an Alpha grammar that precisely tracks parentheses but approximates
 * brackets.
 *
 * Corresponds to the base grammar `D'(Σ_α, Σ_β)` from the article,
 * used as the Alpha grammar in the baseline mutual refinement method.
 *
 * Alpha grammars are used in the mutual refinement algorithm to:
 * 1. Precisely match parentheses (require exact open/close pairs)
 * 2. Approximate brackets (accept any bracket symbol individually)
 *
 * @param[in] n_par       Number of parenthesis types (pairs of `(`/`)`).
 * @param[in] n_bra       Number of bracket types (pairs of `[`/`]`).
 * @param[in] has_normal  Whether to include an epsilon/unlabeled edge terminal.
 *
 * @return A fully initialized Alpha `MRGrammar_t` for mutual refinement,
 *         or a zeroed struct on allocation failure.
 */
MRGrammar_t dyck_alpha_grammar(int64_t n_par, int64_t n_bra, bool has_normal);

/**
 * Creates an Alpha grammar with k-parity conditions for enhanced precision.
 *
 * Corresponds to the grammar `D_parᵏ(Σ_α, Σ_β)` from the article (Sect. 4.1).
 * Special cases by k:
 * - k=1 → `D_par(Σ_α, Σ_β)`, used in method **PAR**
 * - k=2 → `D_par²(Σ_α, Σ_β)`, used in method **PAR2**
 * - k=n → `D_parⁿ(Σ_α, Σ_β)`, used in method **PARk** (parameterized)
 *
 * This grammar extends the basic Alpha grammar by tracking parity of bracket counts
 * across k groups.
 *
 * @param[in] n_par       Number of parenthesis types (pairs of `(`/`)`).
 * @param[in] n_bra       Number of bracket types (pairs of `[`/`]`).
 * @param[in] has_normal  Whether to include an epsilon/unlabeled edge terminal.
 * @param[in] k           Number of parity groups (1 ≤ k ≤ n_bra).
 *
 * @return A fully initialized Alpha `MRGrammar_t` with parity tracking,
 *         or a zeroed struct on allocation failure.
 */
MRGrammar_t dyck_alpha_grammar_k_parity(int64_t n_par, int64_t n_bra,
										bool has_normal, int64_t k);

/**
 * Creates an Alpha grammar with k-parity conditions and label exclusion.
 *
 * Corresponds to the grammar used in the **PARErase** method (Sect. 4.2).
 *
 * This grammar extends [dyckAlphaGrammarKParity] by excluding a specific bracket
 * label from parity tracking. The excluded label does not affect parity states.
 *
 * @param[in] n_par       Number of parenthesis types (pairs of `(`/`)`).
 * @param[in] n_bra       Number of bracket types (pairs of `[`/`]`).
 * @param[in] has_normal  Whether to include an epsilon/unlabeled edge terminal.
 * @param[in] ex_bra      Index of the bracket type to exclude from parity tracking
 *                        (0 ≤ ex_bra < n_bra).
 *
 * @return A fully initialized Alpha `MRGrammar_t` with parity tracking and
 *         label exclusion, or a zeroed struct on allocation failure.
 */
MRGrammar_t dyck_alpha_grammar_k_parity_exclude(int64_t n_par, int64_t n_bra,
												bool has_normal, int64_t k,
												int64_t ex_bra);

/**
 * Creates an Alpha grammar with k-parity conditions and structured equality (valid
 * endpoints).
 *
 * Corresponds to the grammar `D⁺_parᵏ(Σ_α, Σ_β)` from the article (Sect. 4.1),
 * defined as `D_parᵏ(Σ_α, Σ_β) ∩ R'(Σ_β)`, where `R'(Σ_β)` is the regular language
 * accepting only strings that start with an opening bracket and end with a closing
 * bracket. Special cases by k:
 * - k=2 → `D⁺_par²(Σ_α, Σ_β)`, used in methods **PAR2E** (`se`) and **COM** (`all`)
 *
 * This grammar combines parity tracking with the "valid endpoints" condition.
 * It ensures strings start in state QE and end in either QE or QC, improving
 * precision by tracking bracket opening/closing patterns.
 *
 * @param[in] n_par       Number of parenthesis types (pairs of `(`/`)`).
 * @param[in] n_bra       Number of bracket types (pairs of `[`/`]`).
 * @param[in] has_normal  Whether to include an epsilon/unlabeled edge terminal.
 * @param[in] k           Number of parity groups (1 ≤ k ≤ n_bra).
 *
 * @return A fully initialized Alpha `MRGrammar_t` with parity and structured
 *         equality constraints, or a zeroed struct on allocation failure.
 */
MRGrammar_t dyck_alpha_grammar_k_parity_se(int64_t n_par, int64_t n_bra,
										   bool has_normal, int64_t k);

/**
 * Creates a Beta grammar that precisely tracks brackets but approximates
 * parentheses.
 *
 * Corresponds to the base grammar `D'(Σ_β, Σ_α)` from the article,
 * used as the Beta grammar in the baseline mutual refinement method.
 *
 * Beta grammars are used in the mutual refinement algorithm to:
 * 1. Precisely match brackets (require exact open/close pairs)
 * 2. Approximate parentheses (accept any parenthesis symbol individually)
 *
 * @param[in] n_par       Number of parenthesis types (pairs of `(`/`)`).
 * @param[in] n_bra       Number of bracket types (pairs of `[`/`]`).
 * @param[in] has_normal  Whether to include an epsilon/unlabeled edge terminal.
 *
 * @return A fully initialized Beta `MRGrammar_t` for mutual refinement,
 *         or a zeroed struct on allocation failure.
 */
MRGrammar_t dyck_beta_grammar(int64_t n_par, int64_t n_bra, bool has_normal);

/**
 * Creates a Beta grammar with k-parity conditions for enhanced precision.
 *
 * Corresponds to the grammar `D_parᵏ(Σ_β, Σ_α)` from the article (Sect. 4.1).
 * Special cases by k:
 * - k=1 → `D_par(Σ_β, Σ_α)`, used in method **PAR**
 * - k=2 → `D_par²(Σ_β, Σ_α)`, used in method **PAR2**
 * - k=n → `D_parⁿ(Σ_β, Σ_α)`, used in method **PARk** (parameterized)
 *
 * This grammar extends the basic Beta grammar by tracking parity of parenthesis
 * counts across k groups.
 *
 * @param[in] n_par       Number of parenthesis types (pairs of `(`/`)`).
 * @param[in] n_bra       Number of bracket types (pairs of `[`/`]`).
 * @param[in] has_normal  Whether to include an epsilon/unlabeled edge terminal.
 * @param[in] k           Number of parity groups (1 ≤ k ≤ n_par).
 *
 * @return A fully initialized Beta `MRGrammar_t` with parity tracking,
 *         or a zeroed struct on allocation failure.
 */
MRGrammar_t dyck_beta_grammar_k_parity(int64_t n_par, int64_t n_bra, bool has_normal,
									   int64_t k);

/**
 * Creates a Beta grammar with k-parity conditions and structured equality (valid
 * endpoints).
 *
 * Corresponds to the grammar `D⁺_parᵏ(Σ_β, Σ_α)` from the article (Sect. 4.1),
 * defined as `D_parᵏ(Σ_β, Σ_α) ∩ R'(Σ_α)`, where `R'(Σ_α)` is the regular language
 * accepting only strings that start with an opening parenthesis and end with a
 * closing parenthesis. Special cases by k:
 * - k=2 → `D⁺_par²(Σ_β, Σ_α)`, used in methods **PAR2E** (`se`) and **COM** (`all`)
 *
 * This is the Beta counterpart to [dyckAlphaGrammarKParitySe], combining parity
 * tracking with the "valid endpoints" condition for parentheses.
 *
 * @param[in] n_par       Number of parenthesis types (pairs of `(`/`)`).
 * @param[in] n_bra       Number of bracket types (pairs of `[`/`]`).
 * @param[in] has_normal  Whether to include an epsilon/unlabeled edge terminal.
 * @param[in] k           Number of parity groups (1 ≤ k ≤ n_par).
 *
 * @return A fully initialized Beta `MRGrammar_t` with parity and structured
 *         equality constraints, or a zeroed struct on allocation failure.
 */
MRGrammar_t dyck_beta_grammar_k_parity_se(int64_t n_par, int64_t n_bra,
										  bool has_normal, int64_t k);
