#include "grammar.h"

/**
 * Start: S
 *
 * Basic balanced language:
 * B -> eps | N B | C_i B | D_i B | A_i F_i
 * F_i -> B H_i
 * H_i -> B_i B
 *
 * Ascent along stack U (unmatched returns `)_i` above B):
 * U -> eps | N B | C_i B | D_i B | A_i F_i | U H_i
 *
 * A complete reachable path S (an ascent U, followed by a descent `(_i` + any S):
 * S -> eps | N B | C_i B | D_i B | A_i F_i | U H_i | S K_i
 * K_i -> A_i B
 *
 * Terminals:
 * N -> `normal`
 * A_i -> `(_i`
 * B_i -> `)_i`
 * C_i -> `[_i`
 * D_i -> `]_i`
 */
MRGrammar_t dyck_alpha_grammar(int64_t n_par, int64_t n_bra, bool has_normal) {
	int64_t nonterms_count = 3 + 5 * n_par + 2 * n_bra + (has_normal ? 1 : 0);
	int64_t terms_count = 2 * n_par + 2 * n_bra + (has_normal ? 1 : 0);

	int64_t rules_count = 3;
	if (has_normal) {
		rules_count += 4;
	}
	if (n_par > 0) {
		rules_count += 11;
	}
	if (n_bra > 0) {
		rules_count += 8;
	}

	MRGrammarConfig config = {.kind = MR_GRAMMAR_ALPHA, .exclude_index = -1};
	MRGrammar_t gr =
		make_grammar(rules_count, terms_count, nonterms_count, config, /*k=*/1);
	if (!gr.rules) {
		return gr;
	}

	LAGraph_rule_EWCNF *rules = gr.rules;
	int64_t r = 0;

	// --- Nonterminal IDs ---

	int32_t NT_U = 1;
	int32_t NT_B = 2;
	int32_t offset = 3;

	// Blocks for parenthesis
	int32_t A_0 = (n_par > 0) ? offset : -1;
	offset += (int32_t)n_par;
	int32_t B_0 = (n_par > 0) ? offset : -1;
	offset += (int32_t)n_par;
	int32_t F_0 = (n_par > 0) ? offset : -1;
	offset += (int32_t)n_par;
	int32_t H_0 = (n_par > 0) ? offset : -1;
	offset += (int32_t)n_par;
	int32_t K_0 = (n_par > 0) ? offset : -1;
	offset += (int32_t)n_par;

	// Blocks for brackets
	int32_t C_0 = (n_bra > 0) ? offset : -1;
	offset += (int32_t)n_bra;
	int32_t D_0 = (n_bra > 0) ? offset : -1;
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

	// Epsilon rules
	rules[r++] = (LAGraph_rule_EWCNF){NT_B, -1, -1, 0, 0, NULL};
	rules[r++] = (LAGraph_rule_EWCNF){NT_U, -1, -1, 0, 0, NULL};
	rules[r++] = (LAGraph_rule_EWCNF){NT_START, -1, -1, 0, 0, NULL};

	if (has_normal) {
		// B -> N B, U -> N B, S -> N B
		rules[r++] = (LAGraph_rule_EWCNF){NT_B, NT_N, NT_B, 0, 0, NULL};
		rules[r++] = (LAGraph_rule_EWCNF){NT_U, NT_N, NT_B, 0, 0, NULL};
		rules[r++] = (LAGraph_rule_EWCNF){NT_START, NT_N, NT_B, 0, 0, NULL};
		// N -> `normal`
		rules[r++] = (LAGraph_rule_EWCNF){NT_N, TERM(t_normal_idx), -1, 0, 0, NULL};
	}

	// Rules for parentheses
	if (n_par > 0) {
		uint32_t count = (uint32_t)n_par;

		// A_i F_i across B, U, S
		rules[r++] = (LAGraph_rule_EWCNF){
			NT_B,
			A_0,
			F_0,
			count,
			LAGraph_EWNCF_INDEX_PROD_A | LAGraph_EWNCF_INDEX_PROD_B,
			NULL};
		rules[r++] = (LAGraph_rule_EWCNF){
			NT_U,
			A_0,
			F_0,
			count,
			LAGraph_EWNCF_INDEX_PROD_A | LAGraph_EWNCF_INDEX_PROD_B,
			NULL};
		rules[r++] = (LAGraph_rule_EWCNF){
			NT_START,
			A_0,
			F_0,
			count,
			LAGraph_EWNCF_INDEX_PROD_A | LAGraph_EWNCF_INDEX_PROD_B,
			NULL};

		// U -> U H_i, S -> U H_i
		rules[r++] = (LAGraph_rule_EWCNF){
			NT_U, NT_U, H_0, count, LAGraph_EWNCF_INDEX_PROD_B, NULL};
		rules[r++] = (LAGraph_rule_EWCNF){
			NT_START, NT_U, H_0, count, LAGraph_EWNCF_INDEX_PROD_B, NULL};

		// S -> S K_i, K_i -> A_i B
		rules[r++] = (LAGraph_rule_EWCNF){
			NT_START, NT_START, K_0, count, LAGraph_EWNCF_INDEX_PROD_B, NULL};
		rules[r++] = (LAGraph_rule_EWCNF){
			K_0,
			A_0,
			NT_B,
			count,
			LAGraph_EWNCF_INDEX_NONTERM | LAGraph_EWNCF_INDEX_PROD_A,
			NULL};

		// F_i -> B H_i, H_i -> B_i B
		rules[r++] = (LAGraph_rule_EWCNF){
			F_0,
			NT_B,
			H_0,
			count,
			LAGraph_EWNCF_INDEX_NONTERM | LAGraph_EWNCF_INDEX_PROD_B,
			NULL};
		rules[r++] = (LAGraph_rule_EWCNF){
			H_0,
			B_0,
			NT_B,
			count,
			LAGraph_EWNCF_INDEX_NONTERM | LAGraph_EWNCF_INDEX_PROD_A,
			NULL};

		// A_i -> `(_i`, B_i -> `)_i`
		rules[r++] = (LAGraph_rule_EWCNF){A_0,
										  TERM(t_par_open_0),
										  -1,
										  count,
										  LAGraph_EWNCF_INDEX_NONTERM |
											  LAGraph_EWNCF_INDEX_PROD_A,
										  NULL};
		rules[r++] = (LAGraph_rule_EWCNF){B_0,
										  TERM(t_par_close_0),
										  -1,
										  count,
										  LAGraph_EWNCF_INDEX_NONTERM |
											  LAGraph_EWNCF_INDEX_PROD_A,
										  NULL};
	}

	// Rules for brackets
	if (n_bra > 0) {
		uint32_t count = (uint32_t)n_bra;

		// C_i B and D_i B across B, U, S
		rules[r++] = (LAGraph_rule_EWCNF){
			NT_B, C_0, NT_B, count, LAGraph_EWNCF_INDEX_PROD_A, NULL};
		rules[r++] = (LAGraph_rule_EWCNF){
			NT_B, D_0, NT_B, count, LAGraph_EWNCF_INDEX_PROD_A, NULL};
		rules[r++] = (LAGraph_rule_EWCNF){
			NT_U, C_0, NT_B, count, LAGraph_EWNCF_INDEX_PROD_A, NULL};
		rules[r++] = (LAGraph_rule_EWCNF){
			NT_U, D_0, NT_B, count, LAGraph_EWNCF_INDEX_PROD_A, NULL};
		rules[r++] = (LAGraph_rule_EWCNF){
			NT_START, C_0, NT_B, count, LAGraph_EWNCF_INDEX_PROD_A, NULL};
		rules[r++] = (LAGraph_rule_EWCNF){
			NT_START, D_0, NT_B, count, LAGraph_EWNCF_INDEX_PROD_A, NULL};

		// C_i -> `[_i`, D_i -> `]_i`
		rules[r++] = (LAGraph_rule_EWCNF){C_0,
										  TERM(t_bra_open_0),
										  -1,
										  count,
										  LAGraph_EWNCF_INDEX_NONTERM |
											  LAGraph_EWNCF_INDEX_PROD_A,
										  NULL};
		rules[r++] = (LAGraph_rule_EWCNF){D_0,
										  TERM(t_bra_close_0),
										  -1,
										  count,
										  LAGraph_EWNCF_INDEX_NONTERM |
											  LAGraph_EWNCF_INDEX_PROD_A,
										  NULL};
	}

	return gr;
}

/**
 * Start: S_0
 *
 * Basic balanced language with bracket parity mask m:
 * B_m -> eps (m = 0) | N B_m | M_{b,j} B_{m ^ (1 << b)} | E_{im,i} B_{m ^ im}
 * E_{im,i} -> A_i G_{im,i}
 * G_{im,i} -> B_{im} B_i
 *
 * Ascent along stack U (unmatched returns `)_i` above B):
 * U_m -> eps (m = 0) | N B_m | M_{b,j} B_{m ^ (1 << b)} | E_{im,i} B_{m ^ im} |
 * H_{im,i} B_{m ^ im} H_{im,i} -> U_{im} B_i
 *
 * A complete reachable path S (an ascent U, followed by a descent `(_i` + any S):
 * S_m -> eps (m = 0) | N B_m | M_{b,j} B_{m ^ (1 << b)} | E_{im,i} B_{m ^ im} |
 * H_{im,i} B_{m ^ im} | K_{im,i} B_{m ^ im} K_{im,i} -> S_{im} A_i
 *
 * Terminals:
 * N -> `normal`
 * A_i -> `(_i`
 * B_i -> `)_i`
 * M_{b,j} -> `[_{b,j}` | `]_{b,j}`
 */
MRGrammar_t dyck_alpha_grammar_k_parity(int64_t n_par, int64_t n_bra,
										bool has_normal, int64_t k) {
	int64_t num_states = (int64_t)1 << k; // N = 2^k
	uint32_t R = 1;						  // without a state machine

	// Arrays for storing the sizes of brackets and their base IDs
	int32_t bra_group_size[64] = {0};
	int32_t M_group_base[64] = {0};
	int32_t t_bra_open_group_base[64] = {0};
	int32_t t_bra_close_group_base[64] = {0};

	// Distribute n_bra brackets into k bit groups
	int64_t active_groups = 0;
	for (int64_t b = 0; b < k; b++) {
		if (b < n_bra) {
			bra_group_size[b] = (int32_t)((n_bra - b + k - 1) / k);
			active_groups++;
		} else {
			bra_group_size[b] = 0;
		}
	}

	// Nonterminal counts:
	// S_mask, U_mask, B_mask:                 3 * num_states
	// A_i, B_i:                               2 * n_par
	// M_{b,j}:                                n_bra
	// E_{im,i}, G_{im,i}, H_{im,i}, K_{im,i}: for each par i, for each innerMask im,
	//   we need E, G, H, and K.
	//   Total: 4 * n_par * num_states
	// N:                                      1 (if has_normal)
	int64_t nonterms_count = 3 * num_states + 2 * n_par + n_bra +
							 4 * (n_par * num_states) + (has_normal ? 1 : 0);

	// Terminal counts:
	// par_open_i, par_close_i: 2 * n_par
	// bra_open_i, bra_close_i: 2 * n_bra
	// `normal`               : 1 (if has_normal)
	int64_t terms_count = 2 * n_par + 2 * n_bra + (has_normal ? 1 : 0);

	// Rules count:
	// B_0 -> eps, U_0 -> eps, S_0 -> eps:            3
	// B_m -> N B_m, U_m -> N B_m, S_m -> N B_m:      3 (if has_normal)
	// N -> `normal`:                                 1 (if has_normal)
	// Terminal rules for A_i, B_i:                   2 (if n_par > 0)
	// M_{b,j} -> open_{b,j} by groups:               active_groups
	// M_{b,j} -> close_{b,j} by groups:              active_groups
	// B_m, U_m, S_m -> M_{b,j} B_{m ^ (1 << b)}:     3 * active_groups
	// E_{im,i} -> A_i G_{im,i}:                      1 (if n_par > 0)
	// G_{im,i} -> B_im B_i:                          1 (if n_par > 0)
	// H_{im,i} -> U_im B_i:                          1 (if n_par > 0)
	// K_{im,i} -> S_im A_i:                          1 (if n_par > 0)
	// B_m, U_m, S_m -> E_{im,i} B_{m ^ im}:          3 (if n_par > 0)
	// U_m, S_m -> H_{im,i} B_{m ^ im}:               2 (if n_par > 0)
	// S_m -> K_{im,i} B_{m ^ im}:                    1 (if n_par > 0)
	int64_t rules_to_alloc =
		3 + (has_normal ? 4 : 0) + 5 * active_groups + (n_par > 0 ? 12 : 0);

	MRGrammarConfig config = {.kind = MR_GRAMMAR_ALPHA, .exclude_index = -1};
	MRGrammar_t gr =
		make_grammar(rules_to_alloc, terms_count, nonterms_count, config, k);
	if (!gr.rules) {
		return gr;
	}

	LAGraph_rule_EWCNF *rules = gr.rules;
	int64_t r = 0;

	// --- Nonterminal IDs ---

	int32_t NT_U = (int32_t)num_states;
	int32_t NT_B = (int32_t)(2 * num_states);
	int32_t offset = (int32_t)(3 * num_states);

	int32_t A_0 = (n_par > 0) ? offset : -1;
	offset += (int32_t)n_par;
	int32_t B_0 = (n_par > 0) ? offset : -1;
	offset += (int32_t)n_par;

	// Label nonterminal M with continuous groups
	for (int64_t b = 0; b < k; b++) {
		M_group_base[b] = (bra_group_size[b] > 0) ? offset : -1;
		offset += bra_group_size[b];
	}

	int32_t E_base = offset;
	offset += (int32_t)(n_par * num_states);
	int32_t G_base = offset;
	offset += (int32_t)(n_par * num_states);
	int32_t H_base = offset;
	offset += (int32_t)(n_par * num_states);
	int32_t K_base = offset;
	offset += (int32_t)(n_par * num_states);
	int32_t NT_N = has_normal ? (int32_t)(nonterms_count - 1) : -1;

	// --- Terminal IDs ---

	int32_t t_idx = 0;
	int32_t t_par_open_0 = 0;
	t_idx += (int32_t)n_par;
	int32_t t_par_close_0 = t_idx;
	t_idx += (int32_t)n_par;

	// Label terminals for opening and closing brackets by groups
	for (int64_t b = 0; b < k; b++) {
		t_bra_open_group_base[b] = (bra_group_size[b] > 0) ? t_idx : -1;
		t_idx += bra_group_size[b];
	}
	for (int64_t b = 0; b < k; b++) {
		t_bra_close_group_base[b] = (bra_group_size[b] > 0) ? t_idx : -1;
		t_idx += bra_group_size[b];
	}

	int32_t t_normal_idx = has_normal ? t_idx : -1;

	// --- PAT Roles ---
	LAGraph_Pat_role S_role = {LAGraph_PAT_SYM_S, NT_START, /*mask_stride=*/1, 0, 0};
	LAGraph_Pat_role U_role = {LAGraph_PAT_SYM_S, NT_U, /*mask_stride=*/1, 0, 0};
	LAGraph_Pat_role B_role = {LAGraph_PAT_SYM_S, NT_B, /*mask_stride=*/1, 0, 0};
	LAGraph_Pat_role E_role = {LAGraph_PAT_SYM_PATILY, E_base, /*mask_stride=*/n_par,
							   0, n_par};
	LAGraph_Pat_role G_role = {LAGraph_PAT_SYM_PATILY, G_base, /*mask_stride=*/n_par,
							   0, n_par};
	LAGraph_Pat_role H_role = {LAGraph_PAT_SYM_PATILY, H_base, /*mask_stride=*/n_par,
							   0, n_par};
	LAGraph_Pat_role K_role = {LAGraph_PAT_SYM_PATILY, K_base, /*mask_stride=*/n_par,
							   0, n_par};
	LAGraph_Pat_role plain = {LAGraph_PAT_SYM_PLAIN, 0, 0, 0, 0};

	// --- Rules Generation ---

	const LAGraph_rule_Pat *d_eps_b =
		pat_new((LAGraph_rule_Pat){.kind = LAGraph_PAT_EPS,
								   .group = LAGraph_PAT_GROUP_XOR,
								   .masks = (uint32_t)num_states,
								   .states = R,
								   .role = {B_role, plain, plain},
								   .block_symbol = NT_B});

	const LAGraph_rule_Pat *d_eps_u =
		pat_new((LAGraph_rule_Pat){.kind = LAGraph_PAT_EPS,
								   .group = LAGraph_PAT_GROUP_XOR,
								   .masks = (uint32_t)num_states,
								   .states = R,
								   .role = {U_role, plain, plain},
								   .block_symbol = NT_U});

	const LAGraph_rule_Pat *d_eps_s =
		pat_new((LAGraph_rule_Pat){.kind = LAGraph_PAT_EPS,
								   .group = LAGraph_PAT_GROUP_XOR,
								   .masks = (uint32_t)num_states,
								   .states = R,
								   .role = {S_role, plain, plain},
								   .block_symbol = NT_START});

	// B_0 -> eps, U_0 -> eps, S_0 -> eps
	rules[r++] = (LAGraph_rule_EWCNF){NT_B, -1, -1, 0, 0, d_eps_b};
	rules[r++] = (LAGraph_rule_EWCNF){NT_U, -1, -1, 0, 0, d_eps_u};
	rules[r++] = (LAGraph_rule_EWCNF){NT_START, -1, -1, 0, 0, d_eps_s};

	if (has_normal) {
		const LAGraph_rule_Pat *d_scale_b =
			pat_new((LAGraph_rule_Pat){.kind = LAGraph_PAT_SCALE,
									   .group = LAGraph_PAT_GROUP_XOR,
									   .masks = (uint32_t)num_states,
									   .states = R,
									   .idx_count = 1,
									   .role = {B_role, plain, B_role}});

		const LAGraph_rule_Pat *d_scale_u =
			pat_new((LAGraph_rule_Pat){.kind = LAGraph_PAT_SCALE,
									   .group = LAGraph_PAT_GROUP_XOR,
									   .masks = (uint32_t)num_states,
									   .states = R,
									   .idx_count = 1,
									   .role = {U_role, plain, B_role}});

		const LAGraph_rule_Pat *d_scale_s =
			pat_new((LAGraph_rule_Pat){.kind = LAGraph_PAT_SCALE,
									   .group = LAGraph_PAT_GROUP_XOR,
									   .masks = (uint32_t)num_states,
									   .states = R,
									   .idx_count = 1,
									   .role = {S_role, plain, B_role}});

		// B_m -> N B_m, U_m -> N B_m, S_m -> N B_m
		rules[r++] = (LAGraph_rule_EWCNF){NT_B, NT_N, NT_B, 0, 0, d_scale_b};
		rules[r++] = (LAGraph_rule_EWCNF){NT_U, NT_N, NT_B, 0, 0, d_scale_u};
		rules[r++] = (LAGraph_rule_EWCNF){NT_START, NT_N, NT_B, 0, 0, d_scale_s};
		// N -> `normal`
		rules[r++] = (LAGraph_rule_EWCNF){NT_N, TERM(t_normal_idx), -1, 0, 0, NULL};
	}

	// Terminal rules for A_i, B_i
	if (n_par > 0) {
		// A_i -> `(_i`
		rules[r++] = (LAGraph_rule_EWCNF){A_0,
										  TERM(t_par_open_0),
										  -1,
										  (uint32_t)n_par,
										  LAGraph_EWNCF_INDEX_NONTERM |
											  LAGraph_EWNCF_INDEX_PROD_A,
										  NULL};
		// B_i -> `)_i`
		rules[r++] = (LAGraph_rule_EWCNF){B_0,
										  TERM(t_par_close_0),
										  -1,
										  (uint32_t)n_par,
										  LAGraph_EWNCF_INDEX_NONTERM |
											  LAGraph_EWNCF_INDEX_PROD_A,
										  NULL};
	}

	// --- M-patily transition rules (merged open/close for group b) ---
	for (int64_t b = 0; b < k; b++) {
		if (bra_group_size[b] == 0) {
			continue;
		}

		// M_{b,j} -> `[_{b, j}`
		rules[r++] = (LAGraph_rule_EWCNF){M_group_base[b],
										  TERM(t_bra_open_group_base[b]),
										  -1,
										  (uint32_t)bra_group_size[b],
										  LAGraph_EWNCF_INDEX_NONTERM |
											  LAGraph_EWNCF_INDEX_PROD_A,
										  NULL};

		// M_{b,j} -> `]_{b, j}`
		rules[r++] = (LAGraph_rule_EWCNF){M_group_base[b],
										  TERM(t_bra_close_group_base[b]),
										  -1,
										  (uint32_t)bra_group_size[b],
										  LAGraph_EWNCF_INDEX_NONTERM |
											  LAGraph_EWNCF_INDEX_PROD_A,
										  NULL};

		const LAGraph_rule_Pat *d_conv_b =
			pat_new((LAGraph_rule_Pat){.kind = LAGraph_PAT_CONV,
									   .group = LAGraph_PAT_GROUP_XOR,
									   .masks = (uint32_t)num_states,
									   .states = R,
									   .idx_count = (uint32_t)bra_group_size[b],
									   .role = {B_role, plain, B_role},
									   .element = (int64_t)1 << b});

		const LAGraph_rule_Pat *d_conv_u =
			pat_new((LAGraph_rule_Pat){.kind = LAGraph_PAT_CONV,
									   .group = LAGraph_PAT_GROUP_XOR,
									   .masks = (uint32_t)num_states,
									   .states = R,
									   .idx_count = (uint32_t)bra_group_size[b],
									   .role = {U_role, plain, B_role},
									   .element = (int64_t)1 << b});

		const LAGraph_rule_Pat *d_conv_s =
			pat_new((LAGraph_rule_Pat){.kind = LAGraph_PAT_CONV,
									   .group = LAGraph_PAT_GROUP_XOR,
									   .masks = (uint32_t)num_states,
									   .states = R,
									   .idx_count = (uint32_t)bra_group_size[b],
									   .role = {S_role, plain, B_role},
									   .element = (int64_t)1 << b});

		// B_m -> M_{b, j} B_{m ^ (1 << b)}
		rules[r++] = (LAGraph_rule_EWCNF){NT_B,
										  M_group_base[b],
										  NT_B,
										  (uint32_t)bra_group_size[b],
										  LAGraph_EWNCF_INDEX_PROD_A,
										  d_conv_b};

		// U_m -> M_{b, j} B_{m ^ (1 << b)}
		rules[r++] = (LAGraph_rule_EWCNF){NT_U,
										  M_group_base[b],
										  NT_B,
										  (uint32_t)bra_group_size[b],
										  LAGraph_EWNCF_INDEX_PROD_A,
										  d_conv_u};

		// S_m -> M_{b, j} B_{m ^ (1 << b)}
		rules[r++] = (LAGraph_rule_EWCNF){NT_START,
										  M_group_base[b],
										  NT_B,
										  (uint32_t)bra_group_size[b],
										  LAGraph_EWNCF_INDEX_PROD_A,
										  d_conv_s};
	}

	// --- Alpha transition rules ---
	if (n_par > 0) {
		const LAGraph_rule_Pat *d_left =
			pat_new((LAGraph_rule_Pat){.kind = LAGraph_PAT_LEFT,
									   .group = LAGraph_PAT_GROUP_XOR,
									   .masks = (uint32_t)num_states,
									   .states = R,
									   .idx_count = (uint32_t)n_par,
									   .role = {E_role, plain, G_role}});

		const LAGraph_rule_Pat *d_right_gb =
			pat_new((LAGraph_rule_Pat){.kind = LAGraph_PAT_RIGHT,
									   .group = LAGraph_PAT_GROUP_XOR,
									   .masks = (uint32_t)num_states,
									   .states = R,
									   .idx_count = (uint32_t)n_par,
									   .role = {G_role, B_role, plain}});

		const LAGraph_rule_Pat *d_right_hb =
			pat_new((LAGraph_rule_Pat){.kind = LAGraph_PAT_RIGHT,
									   .group = LAGraph_PAT_GROUP_XOR,
									   .masks = (uint32_t)num_states,
									   .states = R,
									   .idx_count = (uint32_t)n_par,
									   .role = {H_role, U_role, plain}});

		const LAGraph_rule_Pat *d_right_kb =
			pat_new((LAGraph_rule_Pat){.kind = LAGraph_PAT_RIGHT,
									   .group = LAGraph_PAT_GROUP_XOR,
									   .masks = (uint32_t)num_states,
									   .states = R,
									   .idx_count = (uint32_t)n_par,
									   .role = {K_role, S_role, plain}});

		const LAGraph_rule_Pat *d_alpha_beb =
			pat_new((LAGraph_rule_Pat){.kind = LAGraph_PAT_CONV,
									   .group = LAGraph_PAT_GROUP_XOR,
									   .masks = (uint32_t)num_states,
									   .states = R,
									   .idx_count = (uint32_t)n_par,
									   .role = {B_role, E_role, B_role},
									   .element = LAGraph_PAT_ALL_ELEMENTS});

		const LAGraph_rule_Pat *d_alpha_ueb =
			pat_new((LAGraph_rule_Pat){.kind = LAGraph_PAT_CONV,
									   .group = LAGraph_PAT_GROUP_XOR,
									   .masks = (uint32_t)num_states,
									   .states = R,
									   .idx_count = (uint32_t)n_par,
									   .role = {U_role, E_role, B_role},
									   .element = LAGraph_PAT_ALL_ELEMENTS});

		const LAGraph_rule_Pat *d_alpha_seb =
			pat_new((LAGraph_rule_Pat){.kind = LAGraph_PAT_CONV,
									   .group = LAGraph_PAT_GROUP_XOR,
									   .masks = (uint32_t)num_states,
									   .states = R,
									   .idx_count = (uint32_t)n_par,
									   .role = {S_role, E_role, B_role},
									   .element = LAGraph_PAT_ALL_ELEMENTS});

		const LAGraph_rule_Pat *d_alpha_uhb =
			pat_new((LAGraph_rule_Pat){.kind = LAGraph_PAT_CONV,
									   .group = LAGraph_PAT_GROUP_XOR,
									   .masks = (uint32_t)num_states,
									   .states = R,
									   .idx_count = (uint32_t)n_par,
									   .role = {U_role, H_role, B_role},
									   .element = LAGraph_PAT_ALL_ELEMENTS});

		const LAGraph_rule_Pat *d_alpha_shb =
			pat_new((LAGraph_rule_Pat){.kind = LAGraph_PAT_CONV,
									   .group = LAGraph_PAT_GROUP_XOR,
									   .masks = (uint32_t)num_states,
									   .states = R,
									   .idx_count = (uint32_t)n_par,
									   .role = {S_role, H_role, B_role},
									   .element = LAGraph_PAT_ALL_ELEMENTS});

		const LAGraph_rule_Pat *d_alpha_skb =
			pat_new((LAGraph_rule_Pat){.kind = LAGraph_PAT_CONV,
									   .group = LAGraph_PAT_GROUP_XOR,
									   .masks = (uint32_t)num_states,
									   .states = R,
									   .idx_count = (uint32_t)n_par,
									   .role = {S_role, K_role, B_role},
									   .element = LAGraph_PAT_ALL_ELEMENTS});

		// LEFT: E_{im,i} -> A_i G_{im,i}
		rules[r++] = (LAGraph_rule_EWCNF){E_base,
										  A_0,
										  G_base,
										  (uint32_t)n_par,
										  LAGraph_EWNCF_INDEX_NONTERM |
											  LAGraph_EWNCF_INDEX_PROD_A |
											  LAGraph_EWNCF_INDEX_PROD_B,
										  d_left};

		// RIGHT: G_{im,i} -> B_im B_i
		rules[r++] = (LAGraph_rule_EWCNF){G_base,
										  NT_B,
										  B_0,
										  (uint32_t)n_par,
										  LAGraph_EWNCF_INDEX_NONTERM |
											  LAGraph_EWNCF_INDEX_PROD_B,
										  d_right_gb};

		// RIGHT: H_{im,i} -> U_im B_i
		rules[r++] = (LAGraph_rule_EWCNF){H_base,
										  NT_U,
										  B_0,
										  (uint32_t)n_par,
										  LAGraph_EWNCF_INDEX_NONTERM |
											  LAGraph_EWNCF_INDEX_PROD_B,
										  d_right_hb};

		// RIGHT: K_{im,i} -> S_im A_i
		rules[r++] = (LAGraph_rule_EWCNF){K_base,
										  NT_START,
										  A_0,
										  (uint32_t)n_par,
										  LAGraph_EWNCF_INDEX_NONTERM |
											  LAGraph_EWNCF_INDEX_PROD_B,
										  d_right_kb};

		// CONV ALL: B_m -> E_{im, i} B_{m ^ im}
		rules[r++] = (LAGraph_rule_EWCNF){NT_B, E_base, NT_B, 0, 0, d_alpha_beb};

		// CONV ALL: U_m -> E_{im, i} B_{m ^ im}
		rules[r++] = (LAGraph_rule_EWCNF){NT_U, E_base, NT_B, 0, 0, d_alpha_ueb};

		// CONV ALL: S_m -> E_{im, i} B_{m ^ im}
		rules[r++] = (LAGraph_rule_EWCNF){NT_START, E_base, NT_B, 0, 0, d_alpha_seb};

		// CONV ALL: U_m -> H_{im, i} B_{m ^ im}
		rules[r++] = (LAGraph_rule_EWCNF){NT_U, H_base, NT_B, 0, 0, d_alpha_uhb};

		// CONV ALL: S_m -> H_{im, i} B_{m ^ im}
		rules[r++] = (LAGraph_rule_EWCNF){NT_START, H_base, NT_B, 0, 0, d_alpha_shb};

		// CONV ALL: S_m -> K_{im, i} B_{m ^ im}
		rules[r++] = (LAGraph_rule_EWCNF){NT_START, K_base, NT_B, 0, 0, d_alpha_skb};
	}

	return gr;
}

MRGrammar_t dyck_alpha_grammar_k_parity_se(int64_t n_par, int64_t n_bra,
										   bool has_normal, int64_t k) {
	int64_t num_parity_states = (int64_t)1 << k;
	uint32_t R = RSTATE_COUNT;
	int64_t s_count = num_parity_states * RSTATE_COUNT * RSTATE_COUNT;

	// Arrays for storing the sizes of brackets and their base IDs
	int32_t bra_group_size[64] = {0};
	int32_t C_group_base[64] = {0};
	int32_t D_group_base[64] = {0};
	int32_t t_bra_open_group_base[64] = {0};
	int32_t t_bra_close_group_base[64] = {0};

	// Distribute n_bra brackets into k bit groups
	int64_t active_groups = 0;
	for (int64_t b = 0; b < k; b++) {
		if (b < n_bra) {
			bra_group_size[b] = (int32_t)((n_bra - b + k - 1) / k);
			active_groups++;
		} else {
			bra_group_size[b] = 0;
		}
	}

	// Nonterminal counts:
	// S (start):                          1
	// Eps:                                1
	// S_{mask, qStart, qEnd}:             s_count
	// A_i, B_i (par aux):                 2 * n_par
	// C_i, D_i (bra aux):                 2 * n_bra
	// E_{im,qs,qmid,i}, G_{im,qs,qmid,i}: for each par i, for each (im,qs,qmid)
	//   triple, we need E and G.
	//   Total: 2 * n_par * num_parity_states * RSTATE_COUNT^2
	// N:                                  1 (if has_normal)
	int64_t eg_count = n_par * num_parity_states * RSTATE_COUNT * RSTATE_COUNT;
	int64_t nonterms_count =
		2 + s_count + 2 * n_par + 2 * n_bra + 2 * eg_count + (has_normal ? 1 : 0);

	int64_t terms_count = 2 * n_par + 2 * n_bra + (has_normal ? 1 : 0);

	// Rules count:
	// Eps -> eps:                                       1
	// S -> S_{0,QE,QE} Eps, S -> S_{0,QE,QC} Eps:       2
	// S_{0,q,q} -> eps (for each q):                    RSTATE_COUNT
	// N -> `normal`:                                    1 (if has_normal)
	// S_base -> N S_base:                               1 (if has_normal)
	// Terminal rules A_i, B_i:                          2 (if n_par > 0)
	// Terminal rules C, D by groups:                    2 * active_groups
	// S_{m,qs,qe} -> C_i S_{next,QO,qe}:                active_groups
	// S_{m,qs,qe} -> D_i S_{next,QC,qe} (qs!=QE):       active_groups
	// E_{im,qs,qmid,i} -> A_i G_{...}:                  1 (if n_par > 0)
	// G_{im,qs,qmid,i} -> S_{im,qs,qmid} B_i:           1 (if n_par > 0)
	// S_{m,qs,qe} -> E_{im,qs,qmid,i} S_{m^im,qmid,qe}: 1 (if n_par > 0)
	int64_t rules_to_alloc = 3 + RSTATE_COUNT + (has_normal ? 2 : 0) +
							 2 * active_groups + 2 * active_groups +
							 ((n_par > 0) ? 5 : 0);

	MRGrammarConfig config = {.kind = MR_GRAMMAR_ALPHA, .exclude_index = -1};
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

	int32_t A_0 = (n_par > 0) ? offset : -1;
	offset += (int32_t)n_par;
	int32_t B_0 = (n_par > 0) ? offset : -1;
	offset += (int32_t)n_par;

	// Label nonterminals C and D with continuous groups
	for (int64_t b = 0; b < k; b++) {
		C_group_base[b] = (bra_group_size[b] > 0) ? offset : -1;
		offset += bra_group_size[b];
	}
	for (int64_t b = 0; b < k; b++) {
		D_group_base[b] = (bra_group_size[b] > 0) ? offset : -1;
		offset += bra_group_size[b];
	}

	int32_t E_base = offset;
	int32_t G_base = offset + (int32_t)eg_count;
	int32_t NT_N = has_normal ? (int32_t)(nonterms_count - 1) : -1;

#define NT_S(mask, qs, qe)                                                          \
	((int32_t)(S_base + (mask) * RSTATE_COUNT * RSTATE_COUNT +                      \
			   (qs) * RSTATE_COUNT + (qe)))

	// --- Terminal IDs ---

	int32_t t_idx = 0;
	int32_t t_par_open_0 = 0;
	t_idx += (int32_t)n_par;
	int32_t t_par_close_0 = t_idx;
	t_idx += (int32_t)n_par;

	// Label terminals for opening and closing brackets by groups
	for (int64_t b = 0; b < k; b++) {
		t_bra_open_group_base[b] = (bra_group_size[b] > 0) ? t_idx : -1;
		t_idx += bra_group_size[b];
	}
	for (int64_t b = 0; b < k; b++) {
		t_bra_close_group_base[b] = (bra_group_size[b] > 0) ? t_idx : -1;
		t_idx += bra_group_size[b];
	}

	int32_t t_normal_idx = has_normal ? t_idx : -1;

	// --- PAT roles ---
	LAGraph_Pat_role S_role = {LAGraph_PAT_SYM_S, S_base, (int64_t)R * (int64_t)R,
							   (int64_t)R, 1};
	LAGraph_Pat_role E_role = {LAGraph_PAT_SYM_PATILY, E_base,
							   (int64_t)R * (int64_t)R * n_par, (int64_t)R * n_par,
							   n_par};
	LAGraph_Pat_role G_role = {LAGraph_PAT_SYM_PATILY, G_base,
							   (int64_t)R * (int64_t)R * n_par, (int64_t)R * n_par,
							   n_par};
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

	if (n_par > 0) {
		// A_i -> `(_i`
		rules[r++] = (LAGraph_rule_EWCNF){A_0,
										  TERM(t_par_open_0),
										  -1,
										  (uint32_t)n_par,
										  LAGraph_EWNCF_INDEX_NONTERM |
											  LAGraph_EWNCF_INDEX_PROD_A,
										  NULL};
		// B_i -> `)_i`
		rules[r++] = (LAGraph_rule_EWCNF){B_0,
										  TERM(t_par_close_0),
										  -1,
										  (uint32_t)n_par,
										  LAGraph_EWNCF_INDEX_NONTERM |
											  LAGraph_EWNCF_INDEX_PROD_A,
										  NULL};
	}

	// Terminal rules for brackets by groups
	for (int64_t b = 0; b < k; b++) {
		if (bra_group_size[b] == 0) {
			continue;
		}

		// C_{b, j} -> `[_{b, j}`
		rules[r++] = (LAGraph_rule_EWCNF){C_group_base[b],
										  TERM(t_bra_open_group_base[b]),
										  -1,
										  (uint32_t)bra_group_size[b],
										  LAGraph_EWNCF_INDEX_NONTERM |
											  LAGraph_EWNCF_INDEX_PROD_A,
										  NULL};
		// D_{b, j} -> `]_{b, j}`
		rules[r++] = (LAGraph_rule_EWCNF){D_group_base[b],
										  TERM(t_bra_close_group_base[b]),
										  -1,
										  (uint32_t)bra_group_size[b],
										  LAGraph_EWNCF_INDEX_NONTERM |
											  LAGraph_EWNCF_INDEX_PROD_A,
										  NULL};
	}

	// --- Beta transitions rules ---
	for (int64_t b = 0; b < k; b++) {
		if (bra_group_size[b] == 0) {
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
									   .idx_count = (uint32_t)bra_group_size[b],
									   .role = {S_role, plain, S_role},
									   .element = (int64_t)1 << b,
									   .state = open_map});

		const LAGraph_rule_Pat *d_close =
			pat_new((LAGraph_rule_Pat){.kind = LAGraph_PAT_CONV,
									   .group = LAGraph_PAT_GROUP_XOR,
									   .masks = (uint32_t)num_parity_states,
									   .states = R,
									   .idx_count = (uint32_t)bra_group_size[b],
									   .role = {S_role, plain, S_role},
									   .element = (int64_t)1 << b,
									   .state = close_map});

		// S_{m,qs,qe} -> C_{b,j} S_{1<<b,QO,qe}  (open bra, always allowed)
		rules[r++] = (LAGraph_rule_EWCNF){S_base,
										  C_group_base[b],
										  S_base,
										  (uint32_t)bra_group_size[b],
										  LAGraph_EWNCF_INDEX_PROD_A,
										  d_open};

		// S_{m,qs,qe} -> D_{b,j} S_{1<<b,QC,qe}  (close bra, forbidden from QE)
		rules[r++] = (LAGraph_rule_EWCNF){S_base,
										  D_group_base[b],
										  S_base,
										  (uint32_t)bra_group_size[b],
										  LAGraph_EWNCF_INDEX_PROD_A,
										  d_close};
	}

	// --- Alpha transitions rules ---
	if (n_par > 0) {
		const LAGraph_rule_Pat *d_left =
			pat_new((LAGraph_rule_Pat){.kind = LAGraph_PAT_LEFT,
									   .group = LAGraph_PAT_GROUP_XOR,
									   .masks = (uint32_t)num_parity_states,
									   .states = R,
									   .idx_count = (uint32_t)n_par,
									   .role = {E_role, plain, G_role}});

		const LAGraph_rule_Pat *d_right =
			pat_new((LAGraph_rule_Pat){.kind = LAGraph_PAT_RIGHT,
									   .group = LAGraph_PAT_GROUP_XOR,
									   .masks = (uint32_t)num_parity_states,
									   .states = R,
									   .idx_count = (uint32_t)n_par,
									   .role = {G_role, S_role, plain}});

		const LAGraph_rule_Pat *d_alpha =
			pat_new((LAGraph_rule_Pat){.kind = LAGraph_PAT_CONV,
									   .group = LAGraph_PAT_GROUP_XOR,
									   .masks = (uint32_t)num_parity_states,
									   .states = R,
									   .idx_count = (uint32_t)n_par,
									   .role = {S_role, E_role, S_role},
									   .element = LAGraph_PAT_ALL_ELEMENTS});

		// E_{im,qs,qmid,i} -> A_i G_{im,qs,qmid,i}
		rules[r++] = (LAGraph_rule_EWCNF){E_base,
										  A_0,
										  G_base,
										  (uint32_t)n_par,
										  LAGraph_EWNCF_INDEX_NONTERM |
											  LAGraph_EWNCF_INDEX_PROD_A |
											  LAGraph_EWNCF_INDEX_PROD_B,
										  d_left};

		// G_{im,qs,qmid,i} -> S_{im,qs,qmid} B_i
		rules[r++] = (LAGraph_rule_EWCNF){G_base,
										  S_base,
										  B_0,
										  (uint32_t)n_par,
										  LAGraph_EWNCF_INDEX_NONTERM |
											  LAGraph_EWNCF_INDEX_PROD_B,
										  d_right};

		// S_{m,qs,qe} -> E_{im,qs,qmid,i} S_{m^im,qmid,qe}
		rules[r++] = (LAGraph_rule_EWCNF){S_base, E_base, S_base, 0, 0, d_alpha};
	}

#undef NT_S

	return gr;
}

MRGrammar_t dyck_alpha_grammar_k_parity_exclude(int64_t n_par, int64_t n_bra,
												bool has_normal, int64_t k,
												int64_t ex_bra) {
	int64_t num_states = (int64_t)1 << k;
	uint32_t R = 1;

	// Arrays for storing the sizes of brackets and their base IDs.
	// Size is 65 to safely accommodate k normal groups + 1 exception group
	int32_t bra_group_size[65] = {0};
	int32_t M_group_base[65] = {0};
	int32_t t_bra_open_group_base[65] = {0};
	int32_t t_bra_close_group_base[65] = {0};

	// Distribute n_bra brackets into k bit groups + 1 exception group
	for (int64_t i = 0; i < n_bra; i++) {
		if (i == ex_bra) {
			bra_group_size[k]++; // Special group for ex_bra (bit shift = 0)
		} else {
			bra_group_size[i % k]++; // Normal bit groups
		}
	}

	int64_t active_groups = 0;
	for (int b = 0; b <= k; b++) {
		if (bra_group_size[b] > 0) {
			active_groups++;
		}
	}

	// Same nonterminal layout as dyck_alpha_grammar_k_parity
	int64_t nonterms_count = num_states + 2 * n_par + n_bra +
							 2 * (n_par * num_states) + (has_normal ? 1 : 0);

	int64_t terms_count = 2 * n_par + 2 * n_bra + (has_normal ? 1 : 0);

	// Rules count:
	// S_0 -> eps:                       1
	// N -> `normal`:                    1 (if has_normal)
	// S_m -> N S_m:                     1 (if has_normal)
	// Terminal rules for A_i, B_i:      2 (if n_par > 0)
	// M_{b,j} -> open_{b,j} by groups:  active_groups
	// M_{b,j} -> close_{b,j} by groups: active_groups
	// S_m -> M_{b,j} S_{1<<b}:          active_groups
	// E_{im,i} -> A_i G_{im,i}:         1 (if n_par > 0)
	// G_{im,i} -> S_im B_i:             1 (if n_par > 0)
	// S_m -> E_{im,i} S_{m ^ im}:       1 (if n_par > 0)
	int64_t rules_to_alloc =
		1 + (has_normal ? 2 : 0) + 3 * active_groups + ((n_par > 0) ? 5 : 0);

	MRGrammarConfig config = {.kind = MR_GRAMMAR_ALPHA, .exclude_index = ex_bra};
	MRGrammar_t gr =
		make_grammar(rules_to_alloc, terms_count, nonterms_count, config, k);
	if (!gr.rules) {
		return gr;
	}

	LAGraph_rule_EWCNF *rules = gr.rules;
	int64_t r = 0;

	// --- Nonterminal IDs ---

	int32_t S_base = 0;
	int32_t offset = (int32_t)num_states;

	int32_t A_0 = (n_par > 0) ? offset : -1;
	offset += (int32_t)n_par;
	int32_t B_0 = (n_par > 0) ? offset : -1;
	offset += (int32_t)n_par;

	// Label nonterminal M with continuous groups (including exception group)
	for (int b = 0; b <= k; b++) {
		M_group_base[b] = (bra_group_size[b] > 0) ? offset : -1;
		offset += bra_group_size[b];
	}

	int32_t E_base = offset;
	offset += (int32_t)(n_par * num_states);
	int32_t G_base = offset;
	offset += (int32_t)(n_par * num_states);

	int32_t NT_N = (has_normal) ? (int32_t)(nonterms_count - 1) : -1;

	// --- Terminal IDs ---

	int32_t t_idx = 0;
	int32_t t_par_open_0 = 0;
	t_idx += (int32_t)n_par;
	int32_t t_par_close_0 = t_idx;
	t_idx += (int32_t)n_par;

	// Label terminals for opening and closing brackets by groups
	for (int b = 0; b <= k; b++) {
		t_bra_open_group_base[b] = (bra_group_size[b] > 0) ? t_idx : -1;
		t_idx += bra_group_size[b];
	}
	for (int b = 0; b <= k; b++) {
		t_bra_close_group_base[b] = (bra_group_size[b] > 0) ? t_idx : -1;
		t_idx += bra_group_size[b];
	}

	int32_t t_normal_idx = has_normal ? t_idx : -1;

	// --- PAT roles ---
	LAGraph_Pat_role S_role = {LAGraph_PAT_SYM_S, S_base, 1, 1, 1};
	LAGraph_Pat_role E_role = {LAGraph_PAT_SYM_PATILY, E_base, n_par, n_par, n_par};
	LAGraph_Pat_role G_role = {LAGraph_PAT_SYM_PATILY, G_base, n_par, n_par, n_par};
	LAGraph_Pat_role plain = {LAGraph_PAT_SYM_PLAIN, 0, 0, 0, 0};

	// --- Rules ---

	const LAGraph_rule_Pat *d_eps =
		pat_new((LAGraph_rule_Pat){.kind = LAGraph_PAT_EPS,
								   .group = LAGraph_PAT_GROUP_XOR,
								   .masks = (uint32_t)num_states,
								   .states = R,
								   .idx_count = 1,
								   .role = {S_role, plain, plain},
								   .block_symbol = S_base});

	// S -> eps
	rules[r++] = (LAGraph_rule_EWCNF){S_base, -1, -1, 0, 0, d_eps};

	if (has_normal) {
		const LAGraph_rule_Pat *d_scale =
			pat_new((LAGraph_rule_Pat){.kind = LAGraph_PAT_SCALE,
									   .group = LAGraph_PAT_GROUP_XOR,
									   .masks = (uint32_t)num_states,
									   .states = R,
									   .idx_count = 1,
									   .role = {S_role, plain, S_role}});

		// N -> `normal`
		rules[r++] = (LAGraph_rule_EWCNF){NT_N, TERM(t_normal_idx), -1, 0, 0, NULL};
		// S_m -> N S_m
		rules[r++] = (LAGraph_rule_EWCNF){S_base, NT_N, S_base, 0, 0, d_scale};
	}

	// --- Terminal rules ---

	if (n_par > 0) {
		// A_i -> `(_i`
		rules[r++] = (LAGraph_rule_EWCNF){A_0,
										  TERM(t_par_open_0),
										  -1,
										  (uint32_t)n_par,
										  LAGraph_EWNCF_INDEX_NONTERM |
											  LAGraph_EWNCF_INDEX_PROD_A,
										  NULL};
		// B_i -> `)_i`
		rules[r++] = (LAGraph_rule_EWCNF){B_0,
										  TERM(t_par_close_0),
										  -1,
										  (uint32_t)n_par,
										  LAGraph_EWNCF_INDEX_NONTERM |
											  LAGraph_EWNCF_INDEX_PROD_A,
										  NULL};
	}

	// Terminal and Beta rules for M by groups
	for (int b = 0; b <= k; b++) {
		if (bra_group_size[b] == 0) {
			continue;
		}

		// M_{b,j} -> `[_{b, j}`
		rules[r++] = (LAGraph_rule_EWCNF){M_group_base[b],
										  TERM(t_bra_open_group_base[b]),
										  -1,
										  (uint32_t)bra_group_size[b],
										  LAGraph_EWNCF_INDEX_NONTERM |
											  LAGraph_EWNCF_INDEX_PROD_A,
										  NULL};

		// M_{b,j} -> `]_{b, j}`
		rules[r++] = (LAGraph_rule_EWCNF){M_group_base[b],
										  TERM(t_bra_close_group_base[b]),
										  -1,
										  (uint32_t)bra_group_size[b],
										  LAGraph_EWNCF_INDEX_NONTERM |
											  LAGraph_EWNCF_INDEX_PROD_A,
										  NULL};

		int64_t bit = (b == k) ? 0 : (int64_t)1 << b;

		const LAGraph_rule_Pat *d_conv =
			pat_new((LAGraph_rule_Pat){.kind = LAGraph_PAT_CONV,
									   .group = LAGraph_PAT_GROUP_XOR,
									   .masks = (uint32_t)num_states,
									   .states = R,
									   .idx_count = (uint32_t)bra_group_size[b],
									   .role = {S_role, plain, S_role},
									   .element = bit,
									   .state = NULL});

		// S_m -> M_{b, j} S_{m ^ bit}
		rules[r++] = (LAGraph_rule_EWCNF){S_base,
										  M_group_base[b],
										  S_base,
										  (uint32_t)bra_group_size[b],
										  LAGraph_EWNCF_INDEX_PROD_A,
										  d_conv};
	}

	// --- Alpha transition rules ---
	if (n_par > 0) {
		const LAGraph_rule_Pat *d_left =
			pat_new((LAGraph_rule_Pat){.kind = LAGraph_PAT_LEFT,
									   .group = LAGraph_PAT_GROUP_XOR,
									   .masks = (uint32_t)num_states,
									   .states = R,
									   .idx_count = (uint32_t)n_par,
									   .role = {E_role, plain, G_role}});

		const LAGraph_rule_Pat *d_right =
			pat_new((LAGraph_rule_Pat){.kind = LAGraph_PAT_RIGHT,
									   .group = LAGraph_PAT_GROUP_XOR,
									   .masks = (uint32_t)num_states,
									   .states = R,
									   .idx_count = (uint32_t)n_par,
									   .role = {G_role, S_role, plain}});

		const LAGraph_rule_Pat *d_alpha =
			pat_new((LAGraph_rule_Pat){.kind = LAGraph_PAT_CONV,
									   .group = LAGraph_PAT_GROUP_XOR,
									   .masks = (uint32_t)num_states,
									   .states = R,
									   .idx_count = (uint32_t)n_par,
									   .role = {S_role, E_role, S_role},
									   .element = LAGraph_PAT_ALL_ELEMENTS});

		// E_{im,i} -> A_i G_{im,i}
		rules[r++] = (LAGraph_rule_EWCNF){E_base,
										  A_0,
										  G_base,
										  (uint32_t)n_par,
										  LAGraph_EWNCF_INDEX_NONTERM |
											  LAGraph_EWNCF_INDEX_PROD_A |
											  LAGraph_EWNCF_INDEX_PROD_B,
										  d_left};

		// G_{im,i} -> S_im B_i
		rules[r++] = (LAGraph_rule_EWCNF){G_base,
										  S_base,
										  B_0,
										  (uint32_t)n_par,
										  LAGraph_EWNCF_INDEX_NONTERM |
											  LAGraph_EWNCF_INDEX_PROD_B,
										  d_right};

		// S_m -> E_{im,i} S_{m^im}
		rules[r++] = (LAGraph_rule_EWCNF){S_base, E_base, S_base, 0, 0, d_alpha};
	}

	return gr;
}
