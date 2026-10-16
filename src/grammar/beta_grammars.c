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

	MRGrammarConfig config = {.kind = MR_GRAMMAR_BETA, .exclude_index = -1};
	MRGrammar_t gr =
		make_grammar(rules_count, terms_count, nonterms_count, config, /*k=*/1);
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
	uint32_t R = 1;						  // without a state machine

	// Arrays for storing the sizes of parentheses and their base IDs
	int32_t par_group_size[64] = {0};
	int32_t M_group_base[64] = {0};
	int32_t t_par_open_group_base[64] = {0};
	int32_t t_par_close_group_base[64] = {0};

	// Distribute n_par parentheses into k bit groups
	int64_t active_groups = 0;
	for (int64_t b = 0; b < k; b++) {
		if (b < n_par) {
			par_group_size[b] = (int32_t)((n_par - b + k - 1) / k);
			active_groups++;
		} else {
			par_group_size[b] = 0;
		}
	}

	// Nonterminal counts:
	// S_mask:             num_states
	// A_i, B_i:           2 * n_par
	// C_i, D_i:           2 * n_bra
	// E_{im,i}, G_{im,i}: for each bra i, for each innerMask im,
	//   we need E and G.
	//   Total: 2 * n_bra * num_states
	// N:                  1 (if has_normal)
	int64_t nonterms_count = num_states + n_par + 2 * n_bra +
							 2 * (n_bra * num_states) + (has_normal ? 1 : 0);

	// Terminal counts:
	// par_open_i, par_close_i: 2 * n_par
	// bra_open_i, bra_close_i: 2 * n_bra
	// `normal`               : 1 (if has_normal)
	int64_t terms_count = 2 * n_par + 2 * n_bra + (has_normal ? 1 : 0);

	// Rules count:
	// S_0 -> eps:                       1
	// S_m -> N S_m:                     1 (if has_normal)
	// N -> `normal`:                    1 (if has_normal)
	// Terminal rules for C_i, D_i:      2 (if n_bra > 0)
	// M_{b,j} -> open_{b,j} by groups:  active_groups
	// M_{b,j} -> close_{b,j} by groups: active_groups
	// S_m -> M_{b,j} S_{1 << b}:        active_groups
	// E_{im,i} -> C_i G_{im,i}:         1 (if n_bra > 0)
	// G_{im,i} -> S_im D_i:             1 (if n_bra > 0)
	// S_m -> E_{im,i} S_{m ^ im}:       1 (if n_bra > 0)
	int64_t rules_to_alloc =
		1 + (has_normal ? 2 : 0) + 3 * active_groups + ((n_bra > 0) ? 5 : 0);

	MRGrammarConfig config = {.kind = MR_GRAMMAR_BETA, .exclude_index = -1};
	MRGrammar_t gr =
		make_grammar(rules_to_alloc, terms_count, nonterms_count, config, k);
	if (!gr.rules) {
		return gr;
	}

	LAGraph_rule_EWCNF *rules = gr.rules;
	int64_t r = 0;

	// --- Nonterminal IDs ---

	int32_t offset = (int32_t)num_states;

	// Label nonterminal M with continuous groups
	for (int64_t b = 0; b < k; b++) {
		M_group_base[b] = (par_group_size[b] > 0) ? offset : -1;
		offset += par_group_size[b];
	}

	int32_t C_0 = (n_bra > 0) ? offset : -1;
	offset += (int32_t)n_bra;
	int32_t D_0 = (n_bra > 0) ? offset : -1;
	offset += (int32_t)n_bra;

	int32_t E_base = offset;
	offset += (int32_t)(n_bra * num_states);
	int32_t G_base = offset;
	offset += (int32_t)(n_bra * num_states);

	int32_t NT_N = (has_normal) ? (int32_t)(nonterms_count - 1) : -1;

	// --- Terminal IDs ---

	int32_t t_idx = 0;

	// Label terminals for opening and closing parentheses by groups
	for (int64_t b = 0; b < k; b++) {
		t_par_open_group_base[b] = (par_group_size[b] > 0) ? t_idx : -1;
		t_idx += par_group_size[b];
	}
	for (int64_t b = 0; b < k; b++) {
		t_par_close_group_base[b] = (par_group_size[b] > 0) ? t_idx : -1;
		t_idx += par_group_size[b];
	}

	int32_t t_bra_open_0 = (n_bra > 0) ? t_idx : -1;
	t_idx += (int32_t)n_bra;
	int32_t t_bra_close_0 = (n_bra > 0) ? t_idx : -1;
	t_idx += (int32_t)n_bra;

	int32_t t_normal_idx = has_normal ? t_idx : -1;

	// --- PAT roles ---
	LAGraph_Pat_role S_role = {LAGraph_PAT_SYM_S, NT_START, 1, 1, 1};
	LAGraph_Pat_role E_role = {LAGraph_PAT_SYM_PATILY, E_base, n_bra, n_bra, n_bra};
	LAGraph_Pat_role G_role = {LAGraph_PAT_SYM_PATILY, G_base, n_bra, n_bra, n_bra};
	LAGraph_Pat_role plain = {LAGraph_PAT_SYM_PLAIN, 0, 0, 0, 0};

	// --- Rules ---

	const LAGraph_rule_Pat *d_eps =
		pat_new((LAGraph_rule_Pat){.kind = LAGraph_PAT_EPS,
								   .group = LAGraph_PAT_GROUP_XOR,
								   .masks = (uint32_t)num_states,
								   .states = 1,
								   .idx_count = 1,
								   .role = {S_role, plain, plain},
								   .block_symbol = NT_START});

	// S_0 -> eps
	rules[r++] = (LAGraph_rule_EWCNF){NT_START, -1, -1, 0, 0, d_eps};

	if (has_normal) {
		const LAGraph_rule_Pat *d_scale =
			pat_new((LAGraph_rule_Pat){.kind = LAGraph_PAT_SCALE,
									   .group = LAGraph_PAT_GROUP_XOR,
									   .masks = (uint32_t)num_states,
									   .states = 1,
									   .idx_count = 1,
									   .role = {S_role, plain, S_role}});

		// N -> `normal`
		rules[r++] = (LAGraph_rule_EWCNF){NT_N, TERM(t_normal_idx), -1, 0, 0, NULL};
		// S_m -> N S_m
		rules[r++] = (LAGraph_rule_EWCNF){NT_START, NT_N, NT_START, 0, 0, d_scale};
	}

	// Terminal rules for C_i, D_i
	if (n_bra > 0) {
		// C_i -> `[_i`
		rules[r++] = (LAGraph_rule_EWCNF){C_0,
										  TERM(t_bra_open_0),
										  -1,
										  (uint32_t)n_bra,
										  LAGraph_EWNCF_INDEX_NONTERM |
											  LAGraph_EWNCF_INDEX_PROD_A,
										  NULL};
		// D_i -> `]_i`
		rules[r++] = (LAGraph_rule_EWCNF){D_0,
										  TERM(t_bra_close_0),
										  -1,
										  (uint32_t)n_bra,
										  LAGraph_EWNCF_INDEX_NONTERM |
											  LAGraph_EWNCF_INDEX_PROD_A,
										  NULL};
	}

	// --- M-patily transition rules (merged open/close for group b) ---
	for (int64_t b = 0; b < k; b++) {
		if (par_group_size[b] == 0) {
			continue;
		}

		// M_{b,j} -> `(_{b,j}`
		rules[r++] = (LAGraph_rule_EWCNF){M_group_base[b],
										  TERM(t_par_open_group_base[b]),
										  -1,
										  (uint32_t)par_group_size[b],
										  LAGraph_EWNCF_INDEX_NONTERM |
											  LAGraph_EWNCF_INDEX_PROD_A,
										  NULL};

		// M_{b,j} -> `)_{b,j}`
		rules[r++] = (LAGraph_rule_EWCNF){M_group_base[b],
										  TERM(t_par_close_group_base[b]),
										  -1,
										  (uint32_t)par_group_size[b],
										  LAGraph_EWNCF_INDEX_NONTERM |
											  LAGraph_EWNCF_INDEX_PROD_A,
										  NULL};

		const LAGraph_rule_Pat *d_conv =
			pat_new((LAGraph_rule_Pat){.kind = LAGraph_PAT_CONV,
									   .group = LAGraph_PAT_GROUP_XOR,
									   .masks = (uint32_t)num_states,
									   .states = 1,
									   .idx_count = (uint32_t)par_group_size[b],
									   .role = {S_role, plain, S_role},
									   .element = (int64_t)1 << b});

		// S_m -> M_{b,j} S_{m ^ (1 << b)}
		rules[r++] = (LAGraph_rule_EWCNF){NT_START,
										  M_group_base[b],
										  NT_START,
										  (uint32_t)par_group_size[b],
										  LAGraph_EWNCF_INDEX_PROD_A,
										  d_conv};
	}

	// --- Beta transition rules ---
	if (n_bra > 0) {
		const LAGraph_rule_Pat *d_left =
			pat_new((LAGraph_rule_Pat){.kind = LAGraph_PAT_LEFT,
									   .group = LAGraph_PAT_GROUP_XOR,
									   .masks = (uint32_t)num_states,
									   .states = 1,
									   .idx_count = (uint32_t)n_bra,
									   .role = {E_role, plain, G_role}});

		const LAGraph_rule_Pat *d_right =
			pat_new((LAGraph_rule_Pat){.kind = LAGraph_PAT_RIGHT,
									   .group = LAGraph_PAT_GROUP_XOR,
									   .masks = (uint32_t)num_states,
									   .states = 1,
									   .idx_count = (uint32_t)n_bra,
									   .role = {G_role, S_role, plain}});

		const LAGraph_rule_Pat *d_alpha =
			pat_new((LAGraph_rule_Pat){.kind = LAGraph_PAT_CONV,
									   .group = LAGraph_PAT_GROUP_XOR,
									   .masks = (uint32_t)num_states,
									   .states = 1,
									   .idx_count = (uint32_t)n_bra,
									   .role = {S_role, E_role, S_role},
									   .element = LAGraph_PAT_ALL_ELEMENTS});

		// LEFT: E_{im,i} -> A_i G_{im,i}
		rules[r++] = (LAGraph_rule_EWCNF){E_base,
										  C_0,
										  G_base,
										  (uint32_t)n_bra,
										  LAGraph_EWNCF_INDEX_NONTERM |
											  LAGraph_EWNCF_INDEX_PROD_A |
											  LAGraph_EWNCF_INDEX_PROD_B,
										  d_left};

		// RIGHT: G_{im,i} -> S_im B_i
		rules[r++] = (LAGraph_rule_EWCNF){G_base,
										  NT_START,
										  D_0,
										  (uint32_t)n_bra,
										  LAGraph_EWNCF_INDEX_NONTERM |
											  LAGraph_EWNCF_INDEX_PROD_B,
										  d_right};

		// CONV ALL: S_m -> E_{im, i} S_{m ^ im}
		rules[r++] = (LAGraph_rule_EWCNF){NT_START, E_base, NT_START, 0, 0, d_alpha};
	}

	return gr;
}

MRGrammar_t dyck_beta_grammar_k_parity_se(int64_t n_par, int64_t n_bra,
										  bool has_normal, int64_t k) {
	int64_t num_parity_states = (int64_t)1 << k;
	uint32_t R = RSTATE_COUNT;
	int64_t s_count = num_parity_states * RSTATE_COUNT * RSTATE_COUNT;

	// Arrays for storing the sizes of brackets and their base IDs
	int32_t par_group_size[64] = {0};
	int32_t A_group_base[64] = {0};
	int32_t B_group_base[64] = {0};
	int32_t t_par_open_group_base[64] = {0};
	int32_t t_par_close_group_base[64] = {0};

	// Distribute n_bra brackets into k bit groups
	int64_t active_groups = 0;
	for (int64_t b = 0; b < k; b++) {
		if (b < n_par) {
			par_group_size[b] = (int32_t)((n_par - b + k - 1) / k);
			active_groups++;
		} else {
			par_group_size[b] = 0;
		}
	}

	// Nonterminal counts (beta-SE: bra roles swapped with alpha-SE):
	// S (start):                          1
	// Eps:                                1
	// S_{mask, qStart, qEnd}:             s_count
	// A_i, B_i (par aux):                 2 * n_par
	// C_i, D_i (bra aux):                 2 * n_bra
	// E_{im,qs,qmid,i}, G_{im,qs,qmid,i}: for each bra i, for each (im,qs,qmid)
	//   triple, we need E and G.
	//   Total: 2 * n_bra * num_parity_states * RSTATE_COUNT^2
	// N:                                  1 (if has_normal)
	int64_t eg_count = n_bra * num_parity_states * RSTATE_COUNT * RSTATE_COUNT;
	int64_t nonterms_count =
		2 + s_count + 2 * n_par + 2 * n_bra + 2 * eg_count + (has_normal ? 1 : 0);

	int64_t terms_count = 2 * n_par + 2 * n_bra + (has_normal ? 1 : 0);

	// Rules count:
	// Eps -> eps:                                       1
	// S -> S_{0,QE,QE} Eps, S -> S_{0,QE,QC} Eps:       2
	// S_{0,q,q} -> eps (for each q):                    RSTATE_COUNT
	// N -> `normal`:                                    1 (if has_normal)
	// S_base -> N S_base:                               1 (if has_normal)
	// Terminal rules C_i, D_i:                          2 (if n_bra > 0)
	// Terminal rules C, D by groups:                    2 * active_groups
	// S_{m,qs,qe} -> A_i S_{next,QO,qe}:                active_groups
	// S_{m,qs,qe} -> B_i S_{next,QC,qe} (qs!=QE):       active_groups
	// E_{im,qs,qmid,i} -> C_i G_{...}:                  1 (if n_bra > 0)
	// G_{im,qs,qmid,i} -> S_{im,qs,qmid} D_i:           1 (if n_bra > 0)
	// S_{m,qs,qe} -> E_{im,qs,qmid,i} S_{next,qmid,qe}: 1 (if n_bra > 0)
	int64_t rules_to_alloc = 3 + RSTATE_COUNT + (has_normal ? 2 : 0) +
							 2 * active_groups + 2 * active_groups +
							 ((n_bra > 0) ? 5 : 0);

	MRGrammarConfig config = {.kind = MR_GRAMMAR_BETA, .exclude_index = -1};
	MRGrammar_t gr =
		make_grammar(rules_to_alloc, terms_count, nonterms_count, config, k);
	if (!gr.rules) {
		return gr;
	}

	LAGraph_rule_EWCNF *rules = gr.rules;
	int64_t r = 0;

	// --- Nonterminal IDs ---

	int32_t NT_EPS = 1;
	int32_t S_base = 2;
	int32_t offset = S_base + (int32_t)s_count;

	// Label nonterminals A and B with continuous groups
	for (int64_t b = 0; b < k; b++) {
		A_group_base[b] = (par_group_size[b] > 0) ? offset : -1;
		offset += par_group_size[b];
	}
	for (int64_t b = 0; b < k; b++) {
		B_group_base[b] = (par_group_size[b] > 0) ? offset : -1;
		offset += par_group_size[b];
	}

	int32_t C_0 = (n_bra > 0) ? offset : -1;
	offset += (int32_t)n_bra;
	int32_t D_0 = (n_bra > 0) ? offset : -1;
	offset += (int32_t)n_bra;

	int32_t E_base = offset;
	int32_t G_base = offset + (int32_t)eg_count;
	int32_t NT_N = (has_normal) ? (int32_t)(nonterms_count - 1) : -1;

#define NT_S(mask, qs, qe)                                                          \
	((int32_t)(S_base + (mask) * RSTATE_COUNT * RSTATE_COUNT +                      \
			   (qs) * RSTATE_COUNT + (qe)))

	// --- Terminal IDs ---

	int32_t t_idx = 0;

	// Label terminals for opening and closing parenthesis by groups
	for (int64_t b = 0; b < k; b++) {
		t_par_open_group_base[b] = (par_group_size[b] > 0) ? t_idx : -1;
		t_idx += par_group_size[b];
	}
	for (int64_t b = 0; b < k; b++) {
		t_par_close_group_base[b] = (par_group_size[b] > 0) ? t_idx : -1;
		t_idx += par_group_size[b];
	}

	int32_t t_bra_open_0 = (n_bra > 0) ? t_idx : -1;
	t_idx += (int32_t)n_bra;
	int32_t t_bra_close_0 = (n_bra > 0) ? t_idx : -1;
	t_idx += (int32_t)n_bra;

	int32_t t_normal_idx = has_normal ? t_idx : -1;

	// --- PAT roles ---
	LAGraph_Pat_role S_role = {LAGraph_PAT_SYM_S, S_base, (int64_t)R * (int64_t)R,
							   (int64_t)R, 1};
	LAGraph_Pat_role E_role = {LAGraph_PAT_SYM_PATILY, E_base,
							   (int64_t)R * (int64_t)R * n_bra, (int64_t)R * n_bra,
							   n_bra};
	LAGraph_Pat_role G_role = {LAGraph_PAT_SYM_PATILY, G_base,
							   (int64_t)R * (int64_t)R * n_bra, (int64_t)R * n_bra,
							   n_bra};
	LAGraph_Pat_role plain = {LAGraph_PAT_SYM_PLAIN, 0, 0, 0, 0};

	// --- Rules ---

	// Eps -> eps
	rules[r++] = (LAGraph_rule_EWCNF){NT_EPS, -1, -1, 0, 0, NULL};

	// S -> S_{0,QE,QE} Eps, S -> S_{0,QE,QC} Eps

	int32_t proj_blocks[2] = {NT_S(0, RSTATE_QE, RSTATE_QE),
							  NT_S(0, RSTATE_QE, RSTATE_QC)};
	for (int j = 0; j < 2; j++) {
		const LAGraph_rule_Pat *d_proj =
			pat_new((LAGraph_rule_Pat){.kind = LAGraph_PAT_PROJECT,
									   .group = LAGraph_PAT_GROUP_XOR,
									   .masks = (uint32_t)num_parity_states,
									   .states = R,
									   .idx_count = 1,
									   .role = {plain, S_role, plain},
									   .block_symbol = proj_blocks[j]});

		rules[r++] =
			(LAGraph_rule_EWCNF){NT_START, proj_blocks[j], NT_EPS, 0, 0, d_proj};
	}

	// S_{0,q,q} -> eps
	for (int q = 0; q < RSTATE_COUNT; q++) {
		const LAGraph_rule_Pat *d_eps =
			pat_new((LAGraph_rule_Pat){.kind = LAGraph_PAT_EPS,
									   .group = LAGraph_PAT_GROUP_XOR,
									   .masks = (uint32_t)num_parity_states,
									   .states = R,
									   .idx_count = 1,
									   .role = {S_role, plain, plain},
									   .block_symbol = NT_S(0, q, q)});
		rules[r++] = (LAGraph_rule_EWCNF){NT_S(0, q, q), -1, -1, 0, 0, d_eps};
	}

	// --- Normal rules ---
	if (has_normal) {
		const LAGraph_rule_Pat *d_scale =
			pat_new((LAGraph_rule_Pat){.kind = LAGraph_PAT_SCALE,
									   .group = LAGraph_PAT_GROUP_XOR,
									   .masks = (uint32_t)num_parity_states,
									   .states = R,
									   .idx_count = 1,
									   .role = {S_role, plain, S_role}});

		// N -> `normal`
		rules[r++] = (LAGraph_rule_EWCNF){NT_N, TERM(t_normal_idx), -1, 0, 0, NULL};
		// S_{m,qs,qe} -> N S_{m,qs,qe}
		rules[r++] = (LAGraph_rule_EWCNF){S_base, NT_N, S_base, 0, 0, d_scale};
	}

	// --- Terminal rules ---

	// Terminal rules for parenthesis by groups
	for (int64_t b = 0; b < k; b++) {
		if (par_group_size[b] == 0) {
			continue;
		}

		// A_{b, j} -> `(_{b, j}`
		rules[r++] = (LAGraph_rule_EWCNF){A_group_base[b],
										  TERM(t_par_open_group_base[b]),
										  -1,
										  (uint32_t)par_group_size[b],
										  LAGraph_EWNCF_INDEX_NONTERM |
											  LAGraph_EWNCF_INDEX_PROD_A,
										  NULL};
		// B_{b, j} -> `)_{b, j}`
		rules[r++] = (LAGraph_rule_EWCNF){B_group_base[b],
										  TERM(t_par_close_group_base[b]),
										  -1,
										  (uint32_t)par_group_size[b],
										  LAGraph_EWNCF_INDEX_NONTERM |
											  LAGraph_EWNCF_INDEX_PROD_A,
										  NULL};
	}

	// Terminal rules for brackets
	if (n_bra > 0) {
		// C_i -> `[_i`
		rules[r++] = (LAGraph_rule_EWCNF){C_0,
										  TERM(t_bra_open_0),
										  -1,
										  (uint32_t)n_bra,
										  LAGraph_EWNCF_INDEX_NONTERM |
											  LAGraph_EWNCF_INDEX_PROD_A,
										  NULL};
		// D_i -> `]_i`
		rules[r++] = (LAGraph_rule_EWCNF){D_0,
										  TERM(t_bra_close_0),
										  -1,
										  (uint32_t)n_bra,
										  LAGraph_EWNCF_INDEX_NONTERM |
											  LAGraph_EWNCF_INDEX_PROD_A,
										  NULL};
	}

	// --- Alpha transitions rules ---
	for (int64_t b = 0; b < k; b++) {
		if (par_group_size[b] == 0) {
			continue;
		}

		LAGraph_Pat_state_entry open_ents[RSTATE_COUNT];
		LAGraph_Pat_state_entry close_ents[RSTATE_COUNT];
		uint32_t n_open = 0, n_close = 0;
		for (int qs = 0; qs < RSTATE_COUNT; qs++) {
			open_ents[n_open++] = (LAGraph_Pat_state_entry){qs, RSTATE_QO};
			if (qs != RSTATE_QE) {
				close_ents[n_close++] = (LAGraph_Pat_state_entry){qs, RSTATE_QC};
			}
		}
		const LAGraph_Pat_state_map *open_map = pat_map_new(n_open, open_ents);
		const LAGraph_Pat_state_map *close_map = pat_map_new(n_close, close_ents);

		const LAGraph_rule_Pat *d_open =
			pat_new((LAGraph_rule_Pat){.kind = LAGraph_PAT_CONV,
									   .group = LAGraph_PAT_GROUP_XOR,
									   .masks = (uint32_t)num_parity_states,
									   .states = R,
									   .idx_count = (uint32_t)par_group_size[b],
									   .role = {S_role, plain, S_role},
									   .element = (int64_t)1 << b,
									   .state = open_map});

		const LAGraph_rule_Pat *d_close =
			pat_new((LAGraph_rule_Pat){.kind = LAGraph_PAT_CONV,
									   .group = LAGraph_PAT_GROUP_XOR,
									   .masks = (uint32_t)num_parity_states,
									   .states = R,
									   .idx_count = (uint32_t)par_group_size[b],
									   .role = {S_role, plain, S_role},
									   .element = (int64_t)1 << b,
									   .state = close_map});

		// S_{m,qs,qe} -> A_{b,j} S_{1<<b,QO,qe}  (open par, always allowed)
		rules[r++] = (LAGraph_rule_EWCNF){S_base,
										  A_group_base[b],
										  S_base,
										  (uint32_t)par_group_size[b],
										  LAGraph_EWNCF_INDEX_PROD_A,
										  d_open};

		// S_{m,qs,qe} -> B_{b,j} S_{1<<b,QC,qe}  (close par, forbidden from QE)
		rules[r++] = (LAGraph_rule_EWCNF){S_base,
										  B_group_base[b],
										  S_base,
										  (uint32_t)par_group_size[b],
										  LAGraph_EWNCF_INDEX_PROD_A,
										  d_close};
	}

	// --- Beta transitions rules ---
	if (n_bra > 0) {
		const LAGraph_rule_Pat *d_left =
			pat_new((LAGraph_rule_Pat){.kind = LAGraph_PAT_LEFT,
									   .group = LAGraph_PAT_GROUP_XOR,
									   .masks = (uint32_t)num_parity_states,
									   .states = R,
									   .idx_count = (uint32_t)n_bra,
									   .role = {E_role, plain, G_role}});

		const LAGraph_rule_Pat *d_right =
			pat_new((LAGraph_rule_Pat){.kind = LAGraph_PAT_RIGHT,
									   .group = LAGraph_PAT_GROUP_XOR,
									   .masks = (uint32_t)num_parity_states,
									   .states = R,
									   .idx_count = (uint32_t)n_bra,
									   .role = {G_role, S_role, plain}});

		const LAGraph_rule_Pat *d_alpha =
			pat_new((LAGraph_rule_Pat){.kind = LAGraph_PAT_CONV,
									   .group = LAGraph_PAT_GROUP_XOR,
									   .masks = (uint32_t)num_parity_states,
									   .states = R,
									   .idx_count = (uint32_t)n_bra,
									   .role = {S_role, E_role, S_role},
									   .element = LAGraph_PAT_ALL_ELEMENTS});

		// E_{im,qs,qmid,i} -> C_i G_{im,qs,qmid,i}
		rules[r++] = (LAGraph_rule_EWCNF){E_base,
										  C_0,
										  G_base,
										  (uint32_t)n_bra,
										  LAGraph_EWNCF_INDEX_NONTERM |
											  LAGraph_EWNCF_INDEX_PROD_A |
											  LAGraph_EWNCF_INDEX_PROD_B,
										  d_left};

		// G_{im,qs,qmid,i} -> S_{im,qs,qmid} D_i
		rules[r++] = (LAGraph_rule_EWCNF){G_base,
										  S_base,
										  D_0,
										  (uint32_t)n_bra,
										  LAGraph_EWNCF_INDEX_NONTERM |
											  LAGraph_EWNCF_INDEX_PROD_B,
										  d_right};

		// S_{m,qs,qe} -> E_{im,qs,qmid,i} S_{m^im,qmid,qe}
		rules[r++] = (LAGraph_rule_EWCNF){S_base, E_base, S_base, 0, 0, d_alpha};
	}

#undef NT_S

	return gr;
}
