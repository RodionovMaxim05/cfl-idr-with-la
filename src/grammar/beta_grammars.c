#include "grammar.h"

/**
 * S -> `(_i` S `)_i` S | `[_i` S | `]_i` S | normal S | eps in WCNF :
 *
 * S -> eps
 * S -> N S
 * N -> normal
 *
 * S -> A_i S
 * S -> B_i S
 * A_i -> `(_i`
 * B_i -> `)_i`
 *
 * S -> C_i F_i
 * F_i -> S H_i
 * H_i -> D_i S
 * C_i -> `[_i`
 * D_i -> `]_i`
 */
MRGrammar_t dyck_beta_grammar(int64_t n_par, int64_t n_bra, bool has_normal) {
	int64_t nonterms_count = 1 + 2 * n_par + 4 * n_bra + (has_normal ? 1 : 0);
	int64_t terms_count = 2 * n_par + 2 * n_bra + (has_normal ? 1 : 0);
	int64_t rules_count = 1 + 4 * n_par + 5 * n_bra + (has_normal ? 2 : 0);

	LAGraph_rule_WCNF *rules = malloc(rules_count * sizeof(LAGraph_rule_WCNF));
	int64_t r = 0;

	// S -> eps
	rules[r++] = (LAGraph_rule_WCNF){NT_S, -1, -1, 0};

	if (has_normal) {
		// N - last nonterm
		int32_t N = (int32_t)(nonterms_count - 1);
		// `normal` - last term
		int32_t term_normal = (int32_t)(terms_count - 1);

		// S -> N S
		rules[r++] = (LAGraph_rule_WCNF){NT_S, N, NT_S, 0};
		// N -> `normal`
		rules[r++] = (LAGraph_rule_WCNF){N, term_normal, -1, 0};
	}

	// Adding rules for each parenthese
	for (int64_t i = 0; i < n_par; i++) {
		int32_t A_i = (int32_t)(1 + 2 * i);
		int32_t B_i = (int32_t)(1 + 2 * i + 1);

		int32_t par_open = (int32_t)(2 * i);
		int32_t par_close = (int32_t)(2 * i + 1);

		// S -> A_i S
		rules[r++] = (LAGraph_rule_WCNF){NT_S, A_i, NT_S, 0};
		// S -> B_i S
		rules[r++] = (LAGraph_rule_WCNF){NT_S, B_i, NT_S, 0};
		// A_i -> `(_i`
		rules[r++] = (LAGraph_rule_WCNF){A_i, par_open, -1, 0};
		// B_i -> `)_i`
		rules[r++] = (LAGraph_rule_WCNF){B_i, par_close, -1, 0};
	}

	// Adding rules for each pair of square brackets
	for (int64_t i = 0; i < n_bra; i++) {
		int32_t C_i = (int32_t)(1 + 2 * n_par + 4 * i);
		int32_t D_i = (int32_t)(1 + 2 * n_par + 4 * i + 1);
		int32_t F_i = (int32_t)(1 + 2 * n_par + 4 * i + 2);
		int32_t H_i = (int32_t)(1 + 2 * n_par + 4 * i + 3);

		int32_t bra_open = (int32_t)(2 * n_par + 2 * i);
		int32_t bra_close = (int32_t)(2 * n_par + 2 * i + 1);

		// S -> C_i F_i
		rules[r++] = (LAGraph_rule_WCNF){NT_S, C_i, F_i, 0};
		// F_i -> S H_i
		rules[r++] = (LAGraph_rule_WCNF){F_i, NT_S, H_i, 0};
		// H_i -> D_i S
		rules[r++] = (LAGraph_rule_WCNF){H_i, D_i, NT_S, 0};
		// C_i -> `[_i`
		rules[r++] = (LAGraph_rule_WCNF){C_i, bra_open, -1, 0};
		// D_i -> `]_i`
		rules[r++] = (LAGraph_rule_WCNF){D_i, bra_close, -1, 0};
	}

	MRGrammar_t gr;
	gr.rules = rules;
	gr.rules_count = rules_count;
	gr.terms_count = terms_count;
	gr.nonterms_count = nonterms_count;
	return gr;
}

MRGrammar_t dyck_beta_grammar_k_parity(int64_t n_par, int64_t n_bra, bool has_normal,
									   int64_t k) {
	int64_t num_states = (int64_t)1 << k; // 2^k

	// Nonterminal counts (mirror of alpha: par and bra roles swapped):
	// S_mask:                    num_states
	// A_i, B_i (par auxiliary):  2 * n_par
	// C_i, D_i (bra auxiliary):  2 * n_bra
	// E_{i,im,m}, G_{i,im,m}:    2 * n_bra * num_states * num_states
	// N (if has_normal):         1
	int64_t nonterms_count = num_states + 2 * n_par + 2 * n_bra +
							 2 * n_bra * num_states * num_states +
							 (has_normal ? 1 : 0);

	// Terminal layout same as alpha
	int64_t terms_count = 2 * n_par + 2 * n_bra + (has_normal ? 1 : 0);

	// Rules count:
	// S_0 -> eps:                             1
	// N -> `normal`, S_m -> N S_m:            1 + num_states   (if has_normal)
	// A_i -> `(_i`, B_i -> `)_i`:             2 * n_par
	// S_m -> A_i S_next, S_m -> B_i S_next:   2 * n_par * num_states
	// C_i -> `[_i`, D_i -> `]_i`:             2 * n_bra
	// S_m -> C_i E, E -> S_im G, G -> D_i S:  3 * n_bra * num_states * num_states
	int64_t rules_count = 1 + (has_normal ? num_states + 1 : 0) + 2 * n_par +
						  2 * n_par * num_states + 2 * n_bra +
						  3 * n_bra * num_states * num_states;

	LAGraph_rule_WCNF *rules = malloc(rules_count * sizeof(LAGraph_rule_WCNF));
	int64_t r = 0;

	// === Nonterminal index helpers ===

#define NT_S_MASK(mask) ((int32_t)(mask))
#define NT_A(i) ((int32_t)(num_states + 2 * (i)))
#define NT_B(i) ((int32_t)(num_states + 2 * (i) + 1))
#define NT_C(i) ((int32_t)(num_states + 2 * n_par + 2 * (i)))
#define NT_D(i) ((int32_t)(num_states + 2 * n_par + 2 * (i) + 1))

	int64_t eg_base = num_states + 2 * n_par + 2 * n_bra;
#define NT_E(i, im, cm)                                                             \
	((int32_t)(eg_base +                                                            \
			   ((i) * num_states * num_states + (im) * num_states + (cm)) * 2))
#define NT_G(i, im, cm)                                                             \
	((int32_t)(eg_base +                                                            \
			   ((i) * num_states * num_states + (im) * num_states + (cm)) * 2 + 1))

	int32_t NT_NORMAL = (int32_t)(nonterms_count - 1);

	// === Terminal index helpers ===

#define TERM_PAR_OPEN(i) ((int32_t)(2 * (i)))
#define TERM_PAR_CLOSE(i) ((int32_t)(2 * (i) + 1))
#define TERM_BRA_OPEN(i) ((int32_t)(2 * n_par + 2 * (i)))
#define TERM_BRA_CLOSE(i) ((int32_t)(2 * n_par + 2 * (i) + 1))

	int32_t TERM_NORMAL = (int32_t)(terms_count - 1);

	// --- S_0 -> eps ---
	rules[r++] = (LAGraph_rule_WCNF){NT_S_MASK(0), -1, -1, 0};

	// --- Normal rules ---
	if (has_normal) {
		rules[r++] = (LAGraph_rule_WCNF){NT_NORMAL, TERM_NORMAL, -1, 0};
		for (int64_t m = 0; m < num_states; m++) {
			rules[r++] =
				(LAGraph_rule_WCNF){NT_S_MASK(m), NT_NORMAL, NT_S_MASK(m), 0};
		}
	}

	// --- Parenthesis rules ---
	// A_i -> `(_i`,  B_i -> `)_i`
	for (int64_t i = 0; i < n_par; i++) {
		rules[r++] = (LAGraph_rule_WCNF){NT_A(i), TERM_PAR_OPEN(i), -1, 0};
		rules[r++] = (LAGraph_rule_WCNF){NT_B(i), TERM_PAR_CLOSE(i), -1, 0};
	}

	// S_m -> A_i S_{m ^ bit}
	// S_m -> B_i S_{m ^ bit}
	// bit for parenthesis i = 1 << (i % k)
	for (int64_t i = 0; i < n_par; i++) {
		int64_t bit = (int64_t)1 << (i % k);
		for (int64_t m = 0; m < num_states; m++) {
			int64_t next = m ^ bit;
			rules[r++] =
				(LAGraph_rule_WCNF){NT_S_MASK(m), NT_A(i), NT_S_MASK(next), 0};
			rules[r++] =
				(LAGraph_rule_WCNF){NT_S_MASK(m), NT_B(i), NT_S_MASK(next), 0};
		}
	}

	// --- Bracket rules ---
	// C_i -> `[_i`,  D_i -> `]_i`
	for (int64_t i = 0; i < n_bra; i++) {
		rules[r++] = (LAGraph_rule_WCNF){NT_C(i), TERM_BRA_OPEN(i), -1, 0};
		rules[r++] = (LAGraph_rule_WCNF){NT_D(i), TERM_BRA_CLOSE(i), -1, 0};
	}

	// For each bra i, currentMask m, innerMask im:
	//   nextMask = m ^ im
	//   S_m        -> C_i E_{i,im,m}
	//   E_{i,im,m} -> S_im G_{i,im,m}
	//   G_{i,im,m} -> D_i S_next
	for (int64_t i = 0; i < n_bra; i++) {
		for (int64_t m = 0; m < num_states; m++) {
			for (int64_t im = 0; im < num_states; im++) {
				int64_t next = m ^ im;
				rules[r++] =
					(LAGraph_rule_WCNF){NT_S_MASK(m), NT_C(i), NT_E(i, im, m), 0};
				rules[r++] = (LAGraph_rule_WCNF){NT_E(i, im, m), NT_S_MASK(im),
												 NT_G(i, im, m), 0};
				rules[r++] =
					(LAGraph_rule_WCNF){NT_G(i, im, m), NT_D(i), NT_S_MASK(next), 0};
			}
		}
	}

#undef NT_S_MASK
#undef NT_A
#undef NT_B
#undef NT_C
#undef NT_D
#undef NT_E
#undef NT_G
#undef TERM_PAR_OPEN
#undef TERM_PAR_CLOSE
#undef TERM_BRA_OPEN
#undef TERM_BRA_CLOSE

	MRGrammar_t gr;
	gr.rules = rules;
	gr.rules_count = rules_count;
	gr.terms_count = terms_count;
	gr.nonterms_count = nonterms_count;
	return gr;
}
