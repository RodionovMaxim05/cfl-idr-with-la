#include "grammar.h"

/**
 * S -> `(_i` S `)_i` S | `[_i` S | `]_i` S | normal S | eps in WCNF :
 *
 * S -> eps
 * S -> N S
 * N -> normal
 *
 * S -> A_i E_i
 * E_i -> S G_i
 * G_i -> B_i S
 * A_i -> `(_i`
 * B_i -> `)_i`
 *
 * S -> C_i S
 * S -> D_i S
 * C_i -> `[_i`
 * D_i -> `]_i`
 */
MRGrammar_t dyck_alpha_grammar(int64_t n_par, int64_t n_bra, bool has_normal) {
	int64_t nonterms_count = 1 + 4 * n_par + 2 * n_bra + (has_normal ? 1 : 0);
	int64_t terms_count = 2 * n_par + 2 * n_bra + (has_normal ? 1 : 0);
	int64_t rules_count = 1 + 5 * n_par + 4 * n_bra + (has_normal ? 2 : 0);

	LAGraph_rule_WCNF *rules = malloc(rules_count * sizeof(LAGraph_rule_WCNF));
	int64_t r = 0;

	// S -> eps
	rules[r++] = (LAGraph_rule_WCNF){NT_START, -1, -1, 0};

	if (has_normal) {
		// N - last nonterm
		int32_t N = (int32_t)(nonterms_count - 1);
		// `normal` - last term
		int32_t term_normal = (int32_t)(terms_count - 1);

		// S -> N S
		rules[r++] = (LAGraph_rule_WCNF){NT_START, N, NT_START, 0};
		// N -> `normal`
		rules[r++] = (LAGraph_rule_WCNF){N, term_normal, -1, 0};
	}

	// Adding rules for each pair of parentheses
	for (int64_t i = 0; i < n_par; i++) {
		int32_t A_i = (int32_t)(1 + 4 * i);
		int32_t B_i = (int32_t)(1 + 4 * i + 1);
		int32_t E_i = (int32_t)(1 + 4 * i + 2);
		int32_t G_i = (int32_t)(1 + 4 * i + 3);

		int32_t par_open = (int32_t)(2 * i);
		int32_t par_close = (int32_t)(2 * i + 1);

		// S -> A_i E_i
		rules[r++] = (LAGraph_rule_WCNF){NT_START, A_i, E_i, 0};
		// E_i -> S G_i
		rules[r++] = (LAGraph_rule_WCNF){E_i, NT_START, G_i, 0};
		// G_i -> B_i S
		rules[r++] = (LAGraph_rule_WCNF){G_i, B_i, NT_START, 0};
		// A_i -> `(_i`
		rules[r++] = (LAGraph_rule_WCNF){A_i, par_open, -1, 0};
		// B_i -> `)_i`
		rules[r++] = (LAGraph_rule_WCNF){B_i, par_close, -1, 0};
	}

	// Adding rules for each square bracket
	for (int64_t i = 0; i < n_bra; i++) {
		int32_t C_i = (int32_t)(1 + 4 * n_par + 2 * i);
		int32_t D_i = (int32_t)(1 + 4 * n_par + 2 * i + 1);

		int32_t bra_open = (int32_t)(2 * n_par + 2 * i);
		int32_t bra_close = (int32_t)(2 * n_par + 2 * i + 1);

		// S -> C_i S
		rules[r++] = (LAGraph_rule_WCNF){NT_START, C_i, NT_START, 0};
		// S -> D_i S
		rules[r++] = (LAGraph_rule_WCNF){NT_START, D_i, NT_START, 0};
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
MRGrammar_t dyck_alpha_grammar_k_parity(int64_t n_par, int64_t n_bra,
										bool has_normal, int64_t k) {
	int64_t num_states = (int64_t)1 << k; // 2^k

	// Nonterminal counts:
	// S_mask:          num_states
	// A_i, B_i:        2 * n_par
	// C_i, D_i:        2 * n_bra
	// E_{i,im,m}, G_{i,im,m}: for each par i, for each (innerMask, currentMask) pair
	//   we need E and G.
	//   Total: 2 * n_par * num_states * num_states (E and G for each combination)
	// N (if has_normal): 1
	int64_t nonterms_count = num_states + 2 * n_par + 2 * n_bra +
							 2 * n_par * num_states * num_states +
							 (has_normal ? 1 : 0);

	// Terminal counts:
	// par_open_i, par_close_i: 2 * n_par
	// bra_open_i, bra_close_i: 2 * n_bra
	// `normal` (if has_normal): 1
	int64_t terms_count = 2 * n_par + 2 * n_bra + (has_normal ? 1 : 0);

	// Rules count:
	// S_0 -> eps:                          1
	// S_m -> N S_m (for each m):           num_states (if has_normal)
	// N -> `normal`:                       1          (if has_normal)
	// For each brId, each mask:
	//   S_m -> C_i S_next:                 n_bra * num_states
	//   S_m -> D_i S_next:                 n_bra * num_states
	// C_i -> `[_i`, D_i -> `]_i`:          2 * n_bra
	// For each parId, each (innerMask, currentMask):
	//   S_m -> A_i E_{i,im,m}:             n_par * num_states * num_states
	//   E_{i,im,m} -> S_im G_{i,im,m}:     n_par * num_states * num_states
	//   G_{i,im,m} -> B_i S_next:          n_par * num_states * num_states
	// 	A_i -> `(_i`, B_i -> `)_i`:         2 * n_par
	int64_t rules_count = 1 + (has_normal ? num_states + 1 : 0) +
						  2 * n_bra * num_states + 2 * n_bra +
						  3 * n_par * num_states * num_states + 2 * n_par;

	LAGraph_rule_WCNF *rules = malloc(rules_count * sizeof(LAGraph_rule_WCNF));
	int64_t r = 0;

// === Nonterminal index helpers ===

// S_mask: index = mask (0 .. num_states-1)
#define NT_S_MASK(mask) ((int32_t)(mask))

// A_i (open par):  num_states + 2*i
#define NT_A(i) ((int32_t)(num_states + 2 * (i)))
// B_i (close par): num_states + 2*i + 1
#define NT_B(i) ((int32_t)(num_states + 2 * (i) + 1))

// C_i (open bra):  num_states + 2*n_par + 2*i
#define NT_C(i) ((int32_t)(num_states + 2 * n_par + 2 * (i)))
// D_i (close bra): num_states + 2*n_par + 2*i + 1
#define NT_D(i) ((int32_t)(num_states + 2 * n_par + 2 * (i) + 1))

	int64_t eg_base = num_states + 2 * n_par + 2 * n_bra;

// E_{i, innerMask, currentMask}: base + (i * num_states * num_states + innerMask
//                                   * num_states + currentMask) * 2
#define NT_E(i, im, cm)                                                             \
	((int32_t)(eg_base +                                                            \
			   ((i) * num_states * num_states + (im) * num_states + (cm)) * 2))
// G_{i, innerMask, currentMask}: base + (i * num_states * num_states + innerMask
//                                   * num_states + currentMask) * 2 + 1
#define NT_G(i, im, cm)                                                             \
	((int32_t)(eg_base +                                                            \
			   ((i) * num_states * num_states + (im) * num_states + (cm)) * 2 + 1))

	// N (normal auxiliary, if has_normal): last nonterminal
	int32_t NT_NORMAL = (int32_t)(nonterms_count - 1);

// === Terminal index helpers ===

// par_open_i:  2*i
#define TERM_PAR_OPEN(i) ((int32_t)(2 * (i)))
// par_close_i: 2*i + 1
#define TERM_PAR_CLOSE(i) ((int32_t)(2 * (i) + 1))

// bra_open_i:  2*n_par + 2*i
#define TERM_BRA_OPEN(i) ((int32_t)(2 * n_par + 2 * (i)))
// bra_close_i: 2*n_par + 2*i + 1
#define TERM_BRA_CLOSE(i) ((int32_t)(2 * n_par + 2 * (i) + 1))

	// `normal`: last terminal
	int32_t TERM_NORMAL = (int32_t)(terms_count - 1);

	// --- S_0 -> eps ---
	rules[r++] = (LAGraph_rule_WCNF){NT_S_MASK(0), -1, -1, 0};

	// --- Normal rules ---
	if (has_normal) {
		// N -> `normal`
		rules[r++] = (LAGraph_rule_WCNF){NT_NORMAL, TERM_NORMAL, -1, 0};
		// S_m -> N S_m  (for each mask)
		for (int64_t m = 0; m < num_states; m++) {
			rules[r++] =
				(LAGraph_rule_WCNF){NT_S_MASK(m), NT_NORMAL, NT_S_MASK(m), 0};
		}
	}

	// --- Bracket rules (beta) ---

	// C_i -> `[_i`,  D_i -> `]_i`
	for (int64_t i = 0; i < n_bra; i++) {
		rules[r++] = (LAGraph_rule_WCNF){NT_C(i), TERM_BRA_OPEN(i), -1, 0};
		rules[r++] = (LAGraph_rule_WCNF){NT_D(i), TERM_BRA_CLOSE(i), -1, 0};
	}

	// S_m -> C_i S_{m ^ bit}
	// S_m -> D_i S_{m ^ bit}
	// bit for bracket i = 1 << (i % k)
	for (int64_t i = 0; i < n_bra; i++) {
		int64_t bit = (int64_t)1 << (i % k);
		for (int64_t m = 0; m < num_states; m++) {
			int64_t next = m ^ bit;
			rules[r++] =
				(LAGraph_rule_WCNF){NT_S_MASK(m), NT_C(i), NT_S_MASK(next), 0};
			rules[r++] =
				(LAGraph_rule_WCNF){NT_S_MASK(m), NT_D(i), NT_S_MASK(next), 0};
		}
	}

	// --- Parenthesis rules (alpha) ---

	// A_i -> `(_i`,  B_i -> `)_i`
	for (int64_t i = 0; i < n_par; i++) {
		rules[r++] = (LAGraph_rule_WCNF){NT_A(i), TERM_PAR_OPEN(i), -1, 0};
		rules[r++] = (LAGraph_rule_WCNF){NT_B(i), TERM_PAR_CLOSE(i), -1, 0};
	}

	// For each par i, currentMask m, innerMask im:
	//   nextMask = m ^ im
	//   S_m        -> A_i E_{i,im,m}
	//   E_{i,im,m} -> S_im G_{i,im,m}
	//   G_{i,im,m} -> B_i S_next
	for (int64_t i = 0; i < n_par; i++) {
		for (int64_t m = 0; m < num_states; m++) {
			for (int64_t im = 0; im < num_states; im++) {
				int64_t next = m ^ im;
				rules[r++] =
					(LAGraph_rule_WCNF){NT_S_MASK(m), NT_A(i), NT_E(i, im, m), 0};
				rules[r++] = (LAGraph_rule_WCNF){NT_E(i, im, m), NT_S_MASK(im),
												 NT_G(i, im, m), 0};
				rules[r++] =
					(LAGraph_rule_WCNF){NT_G(i, im, m), NT_B(i), NT_S_MASK(next), 0};
			}
		}
	}

// Cleanup macros
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

MRGrammar_t dyck_alpha_grammar_k_parity_se(int64_t n_par, int64_t n_bra,
										   bool has_normal, int64_t k) {
	int64_t num_parity_states = (int64_t)1 << k;

	// Nonterminal counts:
	// S (start):                1
	// Eps:                      1
	// S_{mask, qStart, qEnd}:   num_parity_states * RSTATE_COUNT * RSTATE_COUNT
	// A_i, B_i (par aux):       2 * n_par
	// C_i, D_i (bra aux):       2 * n_bra
	// E_{i,im,m,qs,qmid,qe},
	// G_{i,im,m,qs,qmid,qe}:    2 * n_par * num_parity_states * num_parity_states
	//                              * RSTATE_COUNT^3 (i, im, m, qs, qmid, qe)
	// N (if has_normal):        1
	int64_t s_count = num_parity_states * RSTATE_COUNT * RSTATE_COUNT;
	int64_t eg_count = 2 * n_par * num_parity_states * num_parity_states *
					   RSTATE_COUNT * RSTATE_COUNT * RSTATE_COUNT;
	int64_t nonterms_count =
		2 + s_count + 2 * n_par + 2 * n_bra + eg_count + (has_normal ? 1 : 0);

	int64_t terms_count = 2 * n_par + 2 * n_bra + (has_normal ? 1 : 0);

	// Rules count:
	// Eps -> eps:                                        1
	// S -> S_{0,QE,QE} Eps, S -> S_{0,QE,QC} Eps:        2
	// S_{0,q,q} -> eps (for each q):                     RSTATE_COUNT
	// N -> `normal`, S_{m,qs,qe} -> N S_{m,qs,qe}:       1 + s_count (if has_normal)
	// C_i -> `[_i`, D_i -> `]_i`:                        2 * n_bra
	// S_{m,qs,qe} -> C_i S_{next,QO,qe}:                 n_bra * s_count
	// S_{m,qs,qe} -> D_i S_{next,QC,qe} (qs!=QE):        n_bra * num_parity_states *
	//                                                       RSTATE_COUNT *
	//                                                       (RSTATE_COUNT - 1)
	// A_i -> `(_i`, B_i -> `)_i`:                        2 * n_par
	// S_{m,qs,qe} -> A_i E_{...}:                        n_par * s_count *
	//                                                       num_parity_states *
	//                                                       RSTATE_COUNT
	// E_{i,im,m,qs,qmid,qe} -> S_{im,qs,qmid} G_{...}:   n_par * s_count *
	//                                                       num_parity_states *
	//                                                       RSTATE_COUNT
	// G_{i,im,m,qs,qmid,qe} -> B_i S_{next,qmid,qe}:     n_par * s_count *
	//                                                       num_parity_states *
	//                                                       RSTATE_COUNT
	int64_t bra_close_count =
		n_bra * num_parity_states * RSTATE_COUNT * (RSTATE_COUNT - 1);
	int64_t par_eg_count = n_par * s_count * num_parity_states * RSTATE_COUNT;
	int64_t rules_count = 3 + RSTATE_COUNT + (has_normal ? s_count + 1 : 0) +
						  2 * n_bra + n_bra * s_count + bra_close_count + 2 * n_par +
						  3 * par_eg_count;

	LAGraph_rule_WCNF *rules = malloc(rules_count * sizeof(LAGraph_rule_WCNF));
	int64_t r = 0;

	// === Nonterminal index helpers ===

	int32_t NT_EPS = 1;

#define NT_S(mask, qs, qe)                                                          \
	((int32_t)(2 + (mask) * RSTATE_COUNT * RSTATE_COUNT + (qs) * RSTATE_COUNT +     \
			   (qe)))

#define NT_A(i) ((int32_t)(2 + s_count + 2 * (i)))
#define NT_B(i) ((int32_t)(2 + s_count + 2 * (i) + 1))
#define NT_C(i) ((int32_t)(2 + s_count + 2 * n_par + 2 * (i)))
#define NT_D(i) ((int32_t)(2 + s_count + 2 * n_par + 2 * (i) + 1))

	int64_t eg_base = 2 + s_count + 2 * n_par + 2 * n_bra;
	int64_t eg_lin_i = num_parity_states * num_parity_states * RSTATE_COUNT *
					   RSTATE_COUNT * RSTATE_COUNT;
	int64_t eg_lin_im =
		num_parity_states * RSTATE_COUNT * RSTATE_COUNT * RSTATE_COUNT;
	int64_t eg_lin_m = RSTATE_COUNT * RSTATE_COUNT * RSTATE_COUNT;
	int64_t eg_lin_qs = RSTATE_COUNT * RSTATE_COUNT;
	int64_t eg_lin_qmid = RSTATE_COUNT;

#define NT_E(i, im, m, qs, qmid, qe)                                                \
	((int32_t)(eg_base + ((i) * eg_lin_i + (im) * eg_lin_im + (m) * eg_lin_m +      \
						  (qs) * eg_lin_qs + (qmid) * eg_lin_qmid + (qe)) *         \
							 2))

#define NT_G(i, im, m, qs, qmid, qe) ((int32_t)(NT_E(i, im, m, qs, qmid, qe) + 1))

	// N (if has_normal): last nonterminal
	int32_t NT_NORMAL = (int32_t)(nonterms_count - 1);

#define TERM_PAR_OPEN(i) ((int32_t)(2 * (i)))
#define TERM_PAR_CLOSE(i) ((int32_t)(2 * (i) + 1))
#define TERM_BRA_OPEN(i) ((int32_t)(2 * n_par + 2 * (i)))
#define TERM_BRA_CLOSE(i) ((int32_t)(2 * n_par + 2 * (i) + 1))

	int32_t TERM_NORMAL = (int32_t)(terms_count - 1);

	// --- Start rules: S -> S_{0,QE,QE} Eps | S_{0,QE,QC} Eps ---
	rules[r++] = (LAGraph_rule_WCNF){NT_EPS, -1, -1, 0};
	rules[r++] =
		(LAGraph_rule_WCNF){NT_START, NT_S(0, RSTATE_QE, RSTATE_QE), NT_EPS, 0};
	rules[r++] =
		(LAGraph_rule_WCNF){NT_START, NT_S(0, RSTATE_QE, RSTATE_QC), NT_EPS, 0};

	// --- S_{0,q,q} -> eps  (for each q) ---
	for (int64_t q = 0; q < RSTATE_COUNT; q++) {
		rules[r++] = (LAGraph_rule_WCNF){NT_S(0, q, q), -1, -1, 0};
	}

	// --- Normal rules ---
	if (has_normal) {
		// N -> `normal`
		rules[r++] = (LAGraph_rule_WCNF){NT_NORMAL, TERM_NORMAL, -1, 0};
		// S_{m,qs,qe} -> N S_{m,qs,qe}  (for each m, qs, qe)
		for (int64_t m = 0; m < num_parity_states; m++) {
			for (int64_t qs = 0; qs < RSTATE_COUNT; qs++) {
				for (int64_t qe = 0; qe < RSTATE_COUNT; qe++) {
					rules[r++] = (LAGraph_rule_WCNF){NT_S(m, qs, qe), NT_NORMAL,
													 NT_S(m, qs, qe), 0};
				}
			}
		}
	}

	// --- Bracket rules ---

	// C_i -> `[_i`,  D_i -> `]_i`
	for (int64_t i = 0; i < n_bra; i++) {
		rules[r++] = (LAGraph_rule_WCNF){NT_C(i), TERM_BRA_OPEN(i), -1, 0};
		rules[r++] = (LAGraph_rule_WCNF){NT_D(i), TERM_BRA_CLOSE(i), -1, 0};
	}

	for (int64_t i = 0; i < n_bra; i++) {
		int64_t bit = (int64_t)1 << (i % k);
		for (int64_t m = 0; m < num_parity_states; m++) {
			int64_t next = m ^ bit;
			for (int64_t qs = 0; qs < RSTATE_COUNT; qs++) {
				for (int64_t qe = 0; qe < RSTATE_COUNT; qe++) {
					// S_{m,qs,qe} -> C_i S_{next,QO,qe}  (open bra, always allowed)
					rules[r++] = (LAGraph_rule_WCNF){NT_S(m, qs, qe), NT_C(i),
													 NT_S(next, RSTATE_QO, qe), 0};
					// S_{m,qs,qe} -> D_i S_{next,QC,qe}  (close bra, forbidden from
					//   QE)
					if (qs != RSTATE_QE) {
						rules[r++] = (LAGraph_rule_WCNF){
							NT_S(m, qs, qe), NT_D(i), NT_S(next, RSTATE_QC, qe), 0};
					}
				}
			}
		}
	}

	// --- Parenthesis rules ---

	// A_i -> `(_i`,  B_i -> `)_i`
	for (int64_t i = 0; i < n_par; i++) {
		rules[r++] = (LAGraph_rule_WCNF){NT_A(i), TERM_PAR_OPEN(i), -1, 0};
		rules[r++] = (LAGraph_rule_WCNF){NT_B(i), TERM_PAR_CLOSE(i), -1, 0};
	}

	// For each par i, currentMask m, innerMask im, qStart qs, qMid qmid:
	//   nextMask = m ^ im
	//   S_{m,qs,qe}        -> A_i E_{i,im,m,qs,qmid}
	//   E_{i,im,m,qs,qmid} -> S_{im,qs,qmid} G_{i,im,m,qs,qmid}
	//   G_{i,im,m,qs,qmid} -> B_i S_{next,qmid,qe}
	for (int64_t i = 0; i < n_par; i++) {
		for (int64_t m = 0; m < num_parity_states; m++) {
			for (int64_t im = 0; im < num_parity_states; im++) {
				int64_t next = m ^ im;
				for (int64_t qs = 0; qs < RSTATE_COUNT; qs++) {
					for (int64_t qmid = 0; qmid < RSTATE_COUNT; qmid++) {
						for (int64_t qe = 0; qe < RSTATE_COUNT; qe++) {
							rules[r++] =
								(LAGraph_rule_WCNF){NT_S(m, qs, qe), NT_A(i),
													NT_E(i, im, m, qs, qmid, qe), 0};

							rules[r++] = (LAGraph_rule_WCNF){
								NT_E(i, im, m, qs, qmid, qe), NT_S(im, qs, qmid),
								NT_G(i, im, m, qs, qmid, qe), 0};

							rules[r++] = (LAGraph_rule_WCNF){
								NT_G(i, im, m, qs, qmid, qe), NT_B(i),
								NT_S(next, qmid, qe), 0};
						}
					}
				}
			}
		}
	}

// Cleanup macros
#undef NT_S
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

MRGrammar_t dyck_alpha_grammar_k_parity_exclude(int64_t n_par, int64_t n_bra,
												bool has_normal, int64_t k,
												int64_t ex_bra) {
	int64_t num_states = (int64_t)1 << k;

	// Same nonterminal layout as dyck_alpha_grammar_k_parity
	int64_t nonterms_count = num_states + 2 * n_par + 2 * n_bra +
							 2 * n_par * num_states * num_states +
							 (has_normal ? 1 : 0);
	int64_t terms_count = 2 * n_par + 2 * n_bra + (has_normal ? 1 : 0);

	// Same rule count as dyck_alpha_grammar_k_parity
	int64_t rules_count = 1 + (has_normal ? num_states + 1 : 0) +
						  2 * n_bra * num_states + 2 * n_bra +
						  3 * n_par * num_states * num_states + 2 * n_par;

	LAGraph_rule_WCNF *rules = malloc(rules_count * sizeof(LAGraph_rule_WCNF));
	int64_t r = 0;

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

	// --- Bracket rules ---

	// C_i -> `[_i`,  D_i -> `]_i`
	for (int64_t i = 0; i < n_bra; i++) {
		rules[r++] = (LAGraph_rule_WCNF){NT_C(i), TERM_BRA_OPEN(i), -1, 0};
		rules[r++] = (LAGraph_rule_WCNF){NT_D(i), TERM_BRA_CLOSE(i), -1, 0};
	}

	// S_m -> C_i S_{next},  S_m -> D_i S_{next}
	// next = m ^ bit  where bit = 1 << (i % k), except for the excluded label
	// (ex_bra) where next = m
	for (int64_t i = 0; i < n_bra; i++) {
		int64_t bit = (i == ex_bra) ? 0 : (int64_t)1 << (i % k);
		for (int64_t m = 0; m < num_states; m++) {
			int64_t next = m ^ bit;
			rules[r++] =
				(LAGraph_rule_WCNF){NT_S_MASK(m), NT_C(i), NT_S_MASK(next), 0};
			rules[r++] =
				(LAGraph_rule_WCNF){NT_S_MASK(m), NT_D(i), NT_S_MASK(next), 0};
		}
	}

	// --- Parenthesis rules ---

	// A_i -> `(_i`,  B_i -> `)_i`
	for (int64_t i = 0; i < n_par; i++) {
		rules[r++] = (LAGraph_rule_WCNF){NT_A(i), TERM_PAR_OPEN(i), -1, 0};
		rules[r++] = (LAGraph_rule_WCNF){NT_B(i), TERM_PAR_CLOSE(i), -1, 0};
	}

	for (int64_t i = 0; i < n_par; i++) {
		for (int64_t m = 0; m < num_states; m++) {
			for (int64_t im = 0; im < num_states; im++) {
				int64_t next = m ^ im;
				rules[r++] =
					(LAGraph_rule_WCNF){NT_S_MASK(m), NT_A(i), NT_E(i, im, m), 0};
				rules[r++] = (LAGraph_rule_WCNF){NT_E(i, im, m), NT_S_MASK(im),
												 NT_G(i, im, m), 0};
				rules[r++] =
					(LAGraph_rule_WCNF){NT_G(i, im, m), NT_B(i), NT_S_MASK(next), 0};
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
