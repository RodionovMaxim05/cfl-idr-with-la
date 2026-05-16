#include "grammar.h"

/**
 * S -> `(_i` S `)_i` S | `[_i` S | `]_i` S | `normal` S | eps in WCNF :
 *
 * S -> eps
 * S -> N S
 * N -> `normal`
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

	int64_t rules_count = 1;
	if (has_normal) {
		rules_count += 2;
	}
	if (n_par > 0) {
		rules_count += 4;
	}
	if (n_bra > 0) {
		rules_count += 5;
	}

	MRGrammar_t gr = make_grammar(rules_count, terms_count, nonterms_count);
	if (!gr.rules) {
		return gr;
	}

	LAGraph_rule_EWCNF *rules = gr.rules;
	int64_t r = 0;

	// --- Nonterminal IDs ---

	int32_t offset = 1;

	// Blocks for parenthesis
	int32_t A_0 = (n_par > 0) ? offset : -1;
	offset += (int32_t)n_par;
	int32_t B_0 = (n_par > 0) ? offset : -1;
	offset += (int32_t)n_par;

	// Blocks for brackets
	int32_t C_0 = (n_bra > 0) ? offset : -1;
	offset += (int32_t)n_bra;
	int32_t D_0 = (n_bra > 0) ? offset : -1;
	offset += (int32_t)n_bra;
	int32_t F_0 = (n_bra > 0) ? offset : -1;
	offset += (int32_t)n_bra;
	int32_t H_0 = (n_bra > 0) ? offset : -1;
	offset += (int32_t)n_bra;

	int32_t NT_N = has_normal ? (int32_t)(nonterms_count - 1) : -1;

	// --- Terminal IDs ---

	int32_t t_idx = 0;
	int32_t t_par_open_0 = 0;
	t_idx += (int32_t)n_par;
	int32_t t_par_close_0 = t_idx;
	t_idx += (int32_t)n_par;
	int32_t t_bra_open_0 = t_idx;
	t_idx += (int32_t)n_bra;
	int32_t t_bra_close_0 = t_idx;
	t_idx += (int32_t)n_bra;

	int32_t t_normal_idx = has_normal ? t_idx : -1;

	// --- Rules ---

	// S -> eps
	rules[r++] = (LAGraph_rule_EWCNF){NT_START, -1, -1, 0, 0};

	if (has_normal) {
		// S -> N S
		rules[r++] = (LAGraph_rule_EWCNF){NT_START, NT_N, NT_START, 0, 0};
		// N -> `normal`
		rules[r++] = (LAGraph_rule_EWCNF){NT_N, TERM(t_normal_idx), -1, 0, 0};
	}

	// Rules for parentheses
	if (n_par > 0) {
		uint32_t count = (uint32_t)n_par;

		// S -> A_i S
		rules[r++] = (LAGraph_rule_EWCNF){NT_START, A_0, NT_START, count,
										  LAGraph_EWNCF_INDEX_PROD_A};
		// S -> B_i S
		rules[r++] = (LAGraph_rule_EWCNF){NT_START, B_0, NT_START, count,
										  LAGraph_EWNCF_INDEX_PROD_A};
		// A_i -> `(_i`
		rules[r++] = (LAGraph_rule_EWCNF){A_0, TERM(t_par_open_0), -1, count,
										  LAGraph_EWNCF_INDEX_NONTERM |
											  LAGraph_EWNCF_INDEX_PROD_A};
		// B_i -> `)_i`
		rules[r++] = (LAGraph_rule_EWCNF){B_0, TERM(t_par_close_0), -1, count,
										  LAGraph_EWNCF_INDEX_NONTERM |
											  LAGraph_EWNCF_INDEX_PROD_A};
	}

	// Rules for brackets
	if (n_bra > 0) {
		uint32_t count = (uint32_t)n_bra;

		// S -> C_i F_i
		rules[r++] = (LAGraph_rule_EWCNF){NT_START, C_0, F_0, count,
										  LAGraph_EWNCF_INDEX_PROD_A |
											  LAGraph_EWNCF_INDEX_PROD_B};
		// F_i -> S H_i
		rules[r++] = (LAGraph_rule_EWCNF){F_0, NT_START, H_0, count,
										  LAGraph_EWNCF_INDEX_NONTERM |
											  LAGraph_EWNCF_INDEX_PROD_B};
		// H_i -> D_i S
		rules[r++] = (LAGraph_rule_EWCNF){H_0, D_0, NT_START, count,
										  LAGraph_EWNCF_INDEX_NONTERM |
											  LAGraph_EWNCF_INDEX_PROD_A};
		// C_i -> `[_i`
		rules[r++] = (LAGraph_rule_EWCNF){C_0, TERM(t_bra_open_0), -1, count,
										  LAGraph_EWNCF_INDEX_NONTERM |
											  LAGraph_EWNCF_INDEX_PROD_A};
		// D_i -> `]_i`
		rules[r++] = (LAGraph_rule_EWCNF){D_0, TERM(t_bra_close_0), -1, count,
										  LAGraph_EWNCF_INDEX_NONTERM |
											  LAGraph_EWNCF_INDEX_PROD_A};
	}

	return gr;
}

MRGrammar_t dyck_beta_grammar_k_parity(int64_t n_par, int64_t n_bra, bool has_normal,
									   int64_t k) {
	int64_t num_states = (int64_t)1 << k; // 2^k

	// Nonterminal counts:
	// S_mask:          num_states
	// A_i, B_i:        2 * n_par
	// C_i, D_i:        2 * n_bra
	// E_{i,im,m}, G_{i,im,m}: for each par i, for each (innerMask, currentMask) pair
	//   we need E and G.
	//   Total: 2 * n_bra * num_states * num_states (E and G for each combination)
	// N (if has_normal): 1
	int64_t nonterms_count = num_states + 2 * n_par + 2 * n_bra +
							 2 * (n_bra * num_states * num_states) +
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
	//   S_m -> A_i S_next:                 n_par * num_states
	//   S_m -> B_i S_next:                 n_par * num_states
	// A_i -> `(_i`, B_i -> `)_i`:          2 * n_par
	// For each parId, each (innerMask, currentMask):
	//   S_m -> C_i E_{i,im,m}:             num_states * num_states
	//   E_{i,im,m} -> S_im G_{i,im,m}:     num_states * num_states
	//   G_{i,im,m} -> D_i S_next:          num_states * num_states
	// C_i -> `[_i`, D_i -> `]_i`:          2 * n_bra
	int64_t rules_to_alloc = 1 + (has_normal ? num_states + 1 : 0) +
							 2 * (n_par * num_states) + 2 * n_par +
							 (n_bra > 0 ? (3 * num_states * num_states + 2) : 0);

	MRGrammar_t gr = make_grammar(rules_to_alloc, terms_count, nonterms_count);
	if (!gr.rules) {
		return gr;
	}

	LAGraph_rule_EWCNF *rules = gr.rules;
	int64_t r = 0;

	// --- Nonterminal IDs ---

	int32_t offset = (int32_t)num_states;

	int32_t A_0 = (n_par > 0) ? offset : -1;
	offset += (int32_t)n_par;
	int32_t B_0 = (n_par > 0) ? offset : -1;
	offset += (int32_t)n_par;
	int32_t C_0 = (n_bra > 0) ? offset : -1;
	offset += (int32_t)n_bra;
	int32_t D_0 = (n_bra > 0) ? offset : -1;
	offset += (int32_t)n_bra;

	int32_t E_base = offset;
	offset += (int32_t)(num_states * num_states * n_bra);
	int32_t G_base = offset;
	offset += (int32_t)(num_states * num_states * n_bra);
	int32_t NT_N = (has_normal) ? (int32_t)(nonterms_count - 1) : -1;

#define E_ID(m, im, i) (E_base + ((m) * num_states * n_bra + (im) * n_bra + (i)))
#define G_ID(m, im, i) (G_base + ((m) * num_states * n_bra + (im) * n_bra + (i)))

	// --- Terminal IDs ---

	int32_t t_idx = 0;
	int32_t t_par_open_0 = 0;
	t_idx += (int32_t)n_par;
	int32_t t_par_close_0 = t_idx;
	t_idx += (int32_t)n_par;
	int32_t t_bra_open_0 = t_idx;
	t_idx += (int32_t)n_bra;
	int32_t t_bra_close_0 = t_idx;
	t_idx += (int32_t)n_bra;

	int32_t t_normal_idx = has_normal ? t_idx : -1;

	// --- Rules ---

	// S -> eps
	rules[r++] = (LAGraph_rule_EWCNF){NT_START, -1, -1, 0, 0};

	if (has_normal) {
		// N -> `normal`
		rules[r++] = (LAGraph_rule_EWCNF){NT_N, TERM(t_normal_idx), -1, 0, 0};
		// S_m -> N S_m  (for each mask)
		for (int64_t m = 0; m < num_states; m++) {
			rules[r++] = (LAGraph_rule_EWCNF){(int32_t)m, NT_N, (int32_t)m, 0, 0};
		}
	}

	// --- Terminal rules ---
	for (int64_t i = 0; i < n_par; i++) {
		// A_i -> `(_i`
		rules[r++] = (LAGraph_rule_EWCNF){A_0 + (int32_t)i,
										  TERM(t_par_open_0) + (int32_t)i, -1, 0, 0};
	}
	for (int64_t i = 0; i < n_par; i++) {
		// B_i -> `)_i`
		rules[r++] = (LAGraph_rule_EWCNF){
			B_0 + (int32_t)i, TERM(t_par_close_0) + (int32_t)i, -1, 0, 0};
	}
	if (n_bra > 0) {
		// C_i -> `[_i`
		rules[r++] = (LAGraph_rule_EWCNF){
			C_0, TERM(t_bra_open_0), -1, (uint32_t)n_bra,
			LAGraph_EWNCF_INDEX_NONTERM | LAGraph_EWNCF_INDEX_PROD_A};
		// D_i -> `]_i`
		rules[r++] = (LAGraph_rule_EWCNF){
			D_0, TERM(t_bra_close_0), -1, (uint32_t)n_bra,
			LAGraph_EWNCF_INDEX_NONTERM | LAGraph_EWNCF_INDEX_PROD_A};
	}

	// --- Alpha transition rules ---
	for (int64_t i = 0; i < n_par; i++) {
		int64_t bit = (int64_t)1 << (i % k);
		for (int64_t m = 0; m < num_states; m++) {
			int32_t next = (int32_t)(m ^ bit);

			// S_m -> A_i S_{m ^ bit}
			rules[r++] = (LAGraph_rule_EWCNF){
				NT_START + (int32_t)m, A_0 + (int32_t)i, NT_START + next, 0, 0};
			// S_m -> B_i S_{m ^ bit}
			rules[r++] = (LAGraph_rule_EWCNF){
				NT_START + (int32_t)m, B_0 + (int32_t)i, NT_START + next, 0, 0};
		}
	}

	// --- Beta transition rules ---
	if (n_bra > 0) {
		for (int64_t m = 0; m < num_states; m++) {
			for (int64_t im = 0; im < num_states; im++) {
				int64_t next = m ^ im;

				// S_m -> C_i E_{m,im,i}
				rules[r++] = (LAGraph_rule_EWCNF){
					NT_START + (int32_t)m, C_0, E_ID(m, im, 0), (uint32_t)n_bra,
					LAGraph_EWNCF_INDEX_PROD_A | LAGraph_EWNCF_INDEX_PROD_B};

				// E_{m,im,i} -> S_im G_{m,im,i}
				rules[r++] = (LAGraph_rule_EWCNF){
					E_ID(m, im, 0), NT_START + (int32_t)im, G_ID(m, im, 0),
					(uint32_t)n_bra,
					LAGraph_EWNCF_INDEX_NONTERM | LAGraph_EWNCF_INDEX_PROD_B};

				// G_{m,im,i} -> D_i S_next
				rules[r++] = (LAGraph_rule_EWCNF){
					G_ID(m, im, 0), D_0, NT_START + next, (uint32_t)n_bra,
					LAGraph_EWNCF_INDEX_NONTERM | LAGraph_EWNCF_INDEX_PROD_A};
			}
		}
	}

// Cleanup macros
#undef E_ID
#undef G_ID

	return gr;
}

MRGrammar_t dyck_beta_grammar_k_parity_se(int64_t n_par, int64_t n_bra,
										  bool has_normal, int64_t k) {
	int64_t num_parity_states = (int64_t)1 << k;
	int64_t s_count = num_parity_states * RSTATE_COUNT * RSTATE_COUNT;

	// Nonterminal counts (beta-SE: bra roles swapped with alpha-SE):
	// S (start):                1
	// Eps:                      1
	// S_{mask, qStart, qEnd}:   s_count
	// A_i, B_i (par aux):       2 * n_par
	// C_i, D_i (bra aux):       2 * n_bra
	// E_{im,m,qs,qmid,qe,i},
	// G_{im,m,qs,qmid,qe,i}:    2 * n_bra * num_parity_states * num_parity_states
	//                              * RSTATE_COUNT^3 (im, m, qs, qmid, qe, i)
	// N (if has_normal):        1
	int64_t eg_count = n_bra * num_parity_states * num_parity_states * RSTATE_COUNT *
					   RSTATE_COUNT * RSTATE_COUNT;

	int64_t nonterms_count =
		2 + s_count + 2 * n_par + 2 * n_bra + 2 * eg_count + (has_normal ? 1 : 0);

	int64_t terms_count = 2 * n_par + 2 * n_bra + (has_normal ? 1 : 0);

	// Rules count:
	// Eps -> eps:                                        1
	// S -> S_{0,QE,QE} Eps, S -> S_{0,QE,QC} Eps:        2
	// S_{0,q,q} -> eps (for each q):                     RSTATE_COUNT
	// N -> `normal`, S_{m,qs,qe} -> N S_{m,qs,qe}:       1 + s_count (if has_normal)
	// A_i -> `(_i`, B_i -> `)_i`:                        2 * n_par
	// S_{m,qs,qe} -> A_i S_{next,QO,qe}:                 n_par * s_count
	// S_{m,qs,qe} -> B_i S_{next,QC,qe} (qs!=QE):        n_par * num_parity_states *
	//                                                       RSTATE_COUNT *
	//                                                       (RSTATE_COUNT - 1)
	// C_i -> `[_i`, D_i -> `]_i`:                        2
	// S_{m,qs,qe} -> C_i E_{...}:                        s_count *
	//                                                       num_parity_states *
	//                                                       RSTATE_COUNT
	// E_{i,im,m,qs,qmid,qe} -> S_{im,qs,qmid} G_{...}:   s_count *
	//                                                       num_parity_states *
	//                                                       RSTATE_COUNT
	// G_{i,im,m,qs,qmid,qe} -> D_i S_{next,qmid,qe}:     s_count *
	//                                                       num_parity_states *
	//                                                       RSTATE_COUNT
	int64_t par_close_count =
		n_par * num_parity_states * RSTATE_COUNT * (RSTATE_COUNT - 1);
	int64_t bra_eg_count = s_count * num_parity_states * RSTATE_COUNT;
	int64_t rules_to_alloc = 3 + RSTATE_COUNT + (has_normal ? s_count + 1 : 0) +
							 2 * n_par + n_par * s_count + par_close_count +
							 ((n_bra > 0) ? (2 + 3 * bra_eg_count) : 0);

	MRGrammar_t gr = make_grammar(rules_to_alloc, terms_count, nonterms_count);
	if (!gr.rules) {
		return gr;
	}

	LAGraph_rule_EWCNF *rules = gr.rules;
	int64_t r = 0;

	// --- Nonterminal IDs ---

	int32_t NT_EPS = 1;
	int32_t S_base = 2;
	int32_t offset = S_base + (int32_t)s_count;

	int32_t A_0 = (n_par > 0) ? offset : -1;
	offset += (int32_t)n_par;
	int32_t B_0 = (n_par > 0) ? offset : -1;
	offset += (int32_t)n_par;
	int32_t C_0 = (n_bra > 0) ? offset : -1;
	offset += (int32_t)n_bra;
	int32_t D_0 = (n_bra > 0) ? offset : -1;
	offset += (int32_t)n_bra;

	int32_t E_base = offset;
	int32_t G_base = offset + (int32_t)eg_count;
	int32_t NT_N = (has_normal) ? (int32_t)(nonterms_count - 1) : -1;

	// --- Nonterminal IDs helpers ---

#define NT_S(mask, qs, qe)                                                          \
	((int32_t)(S_base + (mask) * RSTATE_COUNT * RSTATE_COUNT +                      \
			   (qs) * RSTATE_COUNT + (qe)))

#define E_IDX(im, m, qs, qmid, qe, i)                                               \
	((int64_t)(im) * num_parity_states * RSTATE_COUNT * RSTATE_COUNT *              \
		 RSTATE_COUNT * n_bra +                                                     \
	 (int64_t)(m) * RSTATE_COUNT * RSTATE_COUNT * RSTATE_COUNT * n_bra +            \
	 (int64_t)(qs) * RSTATE_COUNT * RSTATE_COUNT * n_bra +                          \
	 (int64_t)(qmid) * RSTATE_COUNT * n_bra + (int64_t)(qe) * n_bra + (int64_t)(i))

#define E_ID(im, m, qs, qmid, qe, i)                                                \
	(E_base + (int32_t)E_IDX(im, m, qs, qmid, qe, i))
#define G_ID(im, m, qs, qmid, qe, i)                                                \
	(G_base + (int32_t)E_IDX(im, m, qs, qmid, qe, i))

	// --- Terminal IDs ---

	int32_t t_idx = 0;
	int32_t t_par_open_0 = 0;
	t_idx += (int32_t)n_par;
	int32_t t_par_close_0 = t_idx;
	t_idx += (int32_t)n_par;
	int32_t t_bra_open_0 = t_idx;
	t_idx += (int32_t)n_bra;
	int32_t t_bra_close_0 = t_idx;
	t_idx += (int32_t)n_bra;

	int32_t t_normal_idx = has_normal ? t_idx : -1;

	// --- Rules ---

	rules[r++] = (LAGraph_rule_EWCNF){NT_EPS, -1, -1, 0, 0};
	rules[r++] =
		(LAGraph_rule_EWCNF){NT_START, NT_S(0, RSTATE_QE, RSTATE_QE), NT_EPS, 0, 0};
	rules[r++] =
		(LAGraph_rule_EWCNF){NT_START, NT_S(0, RSTATE_QE, RSTATE_QC), NT_EPS, 0, 0};

	for (int q = 0; q < RSTATE_COUNT; q++) {
		rules[r++] = (LAGraph_rule_EWCNF){NT_S(0, q, q), -1, -1, 0, 0};
	}

	// --- Normal rules ---
	if (has_normal) {
		// N -> `normal`
		rules[r++] = (LAGraph_rule_EWCNF){NT_N, TERM(t_normal_idx), -1, 0, 0};
		// S_{m,qs,qe} -> N S_{m,qs,qe}  (for each m, qs, qe)
		for (int64_t m = 0; m < num_parity_states; m++) {
			for (int64_t qs = 0; qs < RSTATE_COUNT; qs++) {
				for (int64_t qe = 0; qe < RSTATE_COUNT; qe++) {
					rules[r++] = (LAGraph_rule_EWCNF){NT_S(m, qs, qe), NT_N,
													  NT_S(m, qs, qe), 0};
				}
			}
		}
	}

	// --- Terminal rules ---
	for (int64_t i = 0; i < n_par; i++) {
		// A_i -> `(_i`
		rules[r++] = (LAGraph_rule_EWCNF){A_0 + (int32_t)i,
										  TERM(t_par_open_0) + (int32_t)i, -1, 0, 0};
	}
	for (int64_t i = 0; i < n_par; i++) {
		// B_i -> `)_i`
		rules[r++] = (LAGraph_rule_EWCNF){
			B_0 + (int32_t)i, TERM(t_par_close_0) + (int32_t)i, -1, 0, 0};
	}
	if (n_bra > 0) {
		// C_i -> `[_i`
		rules[r++] = (LAGraph_rule_EWCNF){
			C_0, TERM(t_bra_open_0), -1, (uint32_t)n_bra,
			LAGraph_EWNCF_INDEX_NONTERM | LAGraph_EWNCF_INDEX_PROD_A};
		// D_i -> `]_i`
		rules[r++] = (LAGraph_rule_EWCNF){
			D_0, TERM(t_bra_close_0), -1, (uint32_t)n_bra,
			LAGraph_EWNCF_INDEX_NONTERM | LAGraph_EWNCF_INDEX_PROD_A};
	}

	// --- Alpha transitions rules ---
	for (int64_t i = 0; i < n_par; i++) {
		int64_t bit = (int64_t)1 << (i % k);
		for (int64_t m = 0; m < num_parity_states; m++) {
			int64_t next = m ^ bit;
			for (int qs = 0; qs < RSTATE_COUNT; qs++) {
				for (int64_t qe = 0; qe < RSTATE_COUNT; qe++) {
					// S_{m,qs,qe} -> A_i S_{next,QO,qe}  (open par, always allowed)
					rules[r++] =
						(LAGraph_rule_EWCNF){NT_S(m, qs, qe), A_0 + (int32_t)i,
											 NT_S(next, RSTATE_QO, qe), 0, 0};
					// S_{m,qs,qe} -> B_i S_{next,QC,qe}  (close par, forbidden from
					//   QE)
					if (qs != RSTATE_QE) {
						rules[r++] =
							(LAGraph_rule_EWCNF){NT_S(m, qs, qe), B_0 + (int32_t)i,
												 NT_S(next, RSTATE_QC, qe), 0, 0};
					}
				}
			}
		}
	}

	// --- Beta transitions rules ---
	if (n_bra > 0) {
		for (int64_t m = 0; m < num_parity_states; m++) {
			for (int64_t im = 0; im < num_parity_states; im++) {
				int64_t next = m ^ im;
				for (int qs = 0; qs < RSTATE_COUNT; qs++) {
					for (int qmid = 0; qmid < RSTATE_COUNT; qmid++) {
						for (int64_t qe = 0; qe < RSTATE_COUNT; qe++) {
							int32_t cur_e = E_ID(im, m, qs, qmid, qe, 0);
							int32_t cur_g = G_ID(im, m, qs, qmid, qe, 0);

							// S_{m,qs,qe} -> C_i E_{im,m,qs,qmid,qe,i}
							rules[r++] = (LAGraph_rule_EWCNF){
								NT_S(m, qs, qe), C_0, cur_e, (uint32_t)n_bra,
								LAGraph_EWNCF_INDEX_PROD_A |
									LAGraph_EWNCF_INDEX_PROD_B};

							// E_{im,m,qs,qmid,qe,i} -> S_{im,qs,qmid}
							//  G_{im,m,qs,qmid,qe,i}
							rules[r++] = (LAGraph_rule_EWCNF){
								cur_e, NT_S(im, qs, qmid), cur_g, (uint32_t)n_bra,
								LAGraph_EWNCF_INDEX_NONTERM |
									LAGraph_EWNCF_INDEX_PROD_B};

							// G_{im,m,qs,qmid,qe,i} -> D_i S_{next,qmid,qe}
							rules[r++] = (LAGraph_rule_EWCNF){
								cur_g, D_0, NT_S(next, qmid, qe), (uint32_t)n_bra,
								LAGraph_EWNCF_INDEX_NONTERM |
									LAGraph_EWNCF_INDEX_PROD_A};
						}
					}
				}
			}
		}
	}

#undef NT_S
#undef EG_IDX
#undef E_ID
#undef G_ID

	return gr;
}
