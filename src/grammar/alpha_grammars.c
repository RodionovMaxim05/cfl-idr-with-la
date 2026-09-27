#include "grammar.h"

/**
 * S -> `(_i` S `)_i` S | `[_i` S | `]_i` S | `normal` S | eps in WCNF :
 *
 * S -> eps
 * S -> N S
 * N -> `normal`
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

	int64_t rules_count = 1;
	if (has_normal) {
		rules_count += 2;
	}
	if (n_par > 0) {
		rules_count += 5;
	}
	if (n_bra > 0) {
		rules_count += 4;
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

	int32_t offset = 1;

	// Blocks for parenthesis
	int32_t A_0 = (n_par > 0) ? offset : -1;
	offset += (int32_t)n_par;
	int32_t B_0 = (n_par > 0) ? offset : -1;
	offset += (int32_t)n_par;
	int32_t E_0 = (n_par > 0) ? offset : -1;
	offset += (int32_t)n_par;
	int32_t G_0 = (n_par > 0) ? offset : -1;
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

		// S -> A_i E_i
		rules[r++] = (LAGraph_rule_EWCNF){NT_START, A_0, E_0, count,
										  LAGraph_EWNCF_INDEX_PROD_A |
											  LAGraph_EWNCF_INDEX_PROD_B};
		// E_i -> S G_i
		rules[r++] = (LAGraph_rule_EWCNF){E_0, NT_START, G_0, count,
										  LAGraph_EWNCF_INDEX_NONTERM |
											  LAGraph_EWNCF_INDEX_PROD_B};
		// G_i -> B_i S
		rules[r++] = (LAGraph_rule_EWCNF){G_0, B_0, NT_START, count,
										  LAGraph_EWNCF_INDEX_NONTERM |
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

		// S -> C_i S
		rules[r++] = (LAGraph_rule_EWCNF){NT_START, C_0, NT_START, count,
										  LAGraph_EWNCF_INDEX_PROD_A};
		// S -> D_i S
		rules[r++] = (LAGraph_rule_EWCNF){NT_START, D_0, NT_START, count,
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

MRGrammar_t dyck_alpha_grammar_k_parity(int64_t n_par, int64_t n_bra,
										bool has_normal, int64_t k) {
	int64_t num_states = (int64_t)1 << k; // 2^k

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
	// S_mask:          num_states
	// A_i, B_i:        2 * n_par
	// C_i, D_i:        2 * n_bra
	// E_{im,i}, G_{im,i}: for each par i, for each innerMask im,
	//   we need E and G.
	//   Total: 2 * n_par * num_states (E and G for each im)
	// N (if has_normal): 1
	int64_t nonterms_count = num_states + 2 * n_par + 2 * n_bra +
							 2 * (n_par * num_states) + (has_normal ? 1 : 0);

	// Terminal counts:
	// par_open_i, par_close_i: 2 * n_par
	// bra_open_i, bra_close_i: 2 * n_bra
	// `normal` (if has_normal): 1
	int64_t terms_count = 2 * n_par + 2 * n_bra + (has_normal ? 1 : 0);

	// Rules count:
	// S_0 -> eps:                           1
	// S_m -> N S_m (for each m):            num_states (if has_normal)
	// N -> `normal`:                        1          (if has_normal)
	// Terminal rules for A_i, B_i:          2          (if n_par>0)
	// Terminal rules for C, D by groups:    2 * active_groups
	// Beta transitions by groups and masks: 2 * active_groups (XOR OPT)
	// Alpha transitions (if n_par>0):
	//   E_{im,i} -> A_i G_{im,i}:           num_states
	//   G_{im,i} -> S_im B_i:               num_states
	//   S_m -> E_{im,i} S_next:             1 (XOR OPT)
	//   Total: 1 + 2 * num_states
	int64_t rules_to_alloc = 1 + (has_normal ? num_states + 1 : 0) +
							 (n_par > 0 ? 2 : 0) + 2 * active_groups +
							 2 * active_groups +
							 (n_par > 0 ? (1 + 2 * num_states) : 0);

	MRGrammarConfig config = {.kind = MR_GRAMMAR_ALPHA, .exclude_index = -1};
	MRGrammar_t gr =
		make_grammar(rules_to_alloc, terms_count, nonterms_count, config, k);
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

	// Label nonterminals C and D with continuous groups
	for (int64_t b = 0; b < k; b++) {
		C_group_base[b] = offset;
		offset += bra_group_size[b];
	}
	for (int64_t b = 0; b < k; b++) {
		D_group_base[b] = offset;
		offset += bra_group_size[b];
	}

	int32_t E_base = offset;
	offset += (int32_t)(n_par * num_states);
	int32_t G_base = offset;
	offset += (int32_t)(n_par * num_states);

#define E_ID(im) (E_base + (int32_t)((im) * n_par))
#define G_ID(im) (G_base + (int32_t)((im) * n_par))

	int32_t NT_N = (has_normal) ? (int32_t)(nonterms_count - 1) : -1;

	// --- Terminal IDs ---

	int32_t t_idx = 0;
	int32_t t_par_open_0 = 0;
	t_idx += (int32_t)n_par;
	int32_t t_par_close_0 = t_idx;
	t_idx += (int32_t)n_par;

	// Label terminals for opening and closing brackets by groups
	for (int64_t b = 0; b < k; b++) {
		t_bra_open_group_base[b] = t_idx;
		t_idx += bra_group_size[b];
	}
	for (int64_t b = 0; b < k; b++) {
		t_bra_close_group_base[b] = t_idx;
		t_idx += bra_group_size[b];
	}

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
	if (n_par > 0) {
		// A_i -> `(_i`
		rules[r++] = (LAGraph_rule_EWCNF){
			A_0, TERM(t_par_open_0), -1, (uint32_t)n_par,
			LAGraph_EWNCF_INDEX_NONTERM | LAGraph_EWNCF_INDEX_PROD_A};
		// B_i -> `)_i`
		rules[r++] = (LAGraph_rule_EWCNF){
			B_0, TERM(t_par_close_0), -1, (uint32_t)n_par,
			LAGraph_EWNCF_INDEX_NONTERM | LAGraph_EWNCF_INDEX_PROD_A};
	}

	// Terminal rules for brackets
	for (int64_t b = 0; b < k; b++) {
		if (bra_group_size[b] == 0) {
			continue;
		}

		// C_{b, j} -> `[_{b, j}`
		rules[r++] = (LAGraph_rule_EWCNF){
			C_group_base[b], TERM(t_bra_open_group_base[b]), -1,
			(uint32_t)bra_group_size[b],
			LAGraph_EWNCF_INDEX_NONTERM | LAGraph_EWNCF_INDEX_PROD_A};
		// D_{b, j} -> `]_{b, j}`
		rules[r++] = (LAGraph_rule_EWCNF){
			D_group_base[b], TERM(t_bra_close_group_base[b]), -1,
			(uint32_t)bra_group_size[b],
			LAGraph_EWNCF_INDEX_NONTERM | LAGraph_EWNCF_INDEX_PROD_A};
	}

	// --- Beta transition rules ---
	for (int64_t b = 0; b < k; b++) {
		if (bra_group_size[b] == 0) {
			continue;
		}
		int64_t bit = (int64_t)1 << b;
		int64_t cnt = bra_group_size[b];

		// S_m -> C_{b, j} S_next
		XorFamilyRouting *rt_open =
			(XorFamilyRouting *)calloc(1, sizeof(XorFamilyRouting));
		rt_open->S_operand_base = NT_START;
		rt_open->S_target_base = NT_START;
		rt_open->S_stride = 1;
		rt_open->N = num_states;
		rt_open->active_masks_count = cnt;
		rt_open->active_masks = (int64_t *)malloc(cnt * sizeof(int64_t));
		rt_open->p_ids = (int32_t *)malloc(cnt * sizeof(int32_t));
		for (int64_t j = 0; j < cnt; j++) {
			rt_open->active_masks[j] = bit;
			rt_open->p_ids[j] = C_group_base[b] + (int32_t)j;
		}
		LAGraph_rule_EWCNF rule_open = {0};
		rule_open.nonterm = NT_START;
		rule_open.prod_A = C_group_base[b];
		rule_open.prod_B = NT_START;
		rule_open.xor_routing = rt_open;
		rules[r++] = rule_open;

		// S_m -> D_{b, j} S_next
		XorFamilyRouting *rt_close =
			(XorFamilyRouting *)calloc(1, sizeof(XorFamilyRouting));
		rt_close->S_operand_base = NT_START;
		rt_close->S_target_base = NT_START;
		rt_close->S_stride = 1;
		rt_close->N = num_states;
		rt_close->active_masks_count = cnt;
		rt_close->active_masks = (int64_t *)malloc(cnt * sizeof(int64_t));
		rt_close->p_ids = (int32_t *)malloc(cnt * sizeof(int32_t));
		for (int64_t j = 0; j < cnt; j++) {
			rt_close->active_masks[j] = bit;
			rt_close->p_ids[j] = D_group_base[b] + (int32_t)j;
		}

		LAGraph_rule_EWCNF rule_close = {0};
		rule_close.nonterm = NT_START;
		rule_close.prod_A = D_group_base[b];
		rule_close.prod_B = NT_START;
		rule_close.xor_routing = rt_close;
		rules[r++] = rule_close;
	}

	// --- Alpha transition rules ---
	if (n_par > 0) {
		// E_{im,i} -> A_i G_{im,i}
		for (int64_t im = 0; im < num_states; im++) {
			rules[r++] = (LAGraph_rule_EWCNF){
				E_ID(im), A_0, G_ID(im), (uint32_t)n_par,
				LAGraph_EWNCF_INDEX_NONTERM | LAGraph_EWNCF_INDEX_PROD_A |
					LAGraph_EWNCF_INDEX_PROD_B};
		}

		// G_{im,i} -> S_im B_i
		for (int64_t im = 0; im < num_states; im++) {
			rules[r++] = (LAGraph_rule_EWCNF){
				G_ID(im), NT_START + (int32_t)im, B_0, (uint32_t)n_par,
				LAGraph_EWNCF_INDEX_NONTERM | LAGraph_EWNCF_INDEX_PROD_B};
		}

		int64_t total_count = num_states * n_par;

		// S_m -> E_{im,i} S_next
		XorFamilyRouting *rt_alpha =
			(XorFamilyRouting *)calloc(1, sizeof(XorFamilyRouting));
		rt_alpha->S_operand_base = NT_START;
		rt_alpha->S_target_base = NT_START;
		rt_alpha->S_stride = 1;
		rt_alpha->N = num_states;
		rt_alpha->active_masks_count = total_count;
		rt_alpha->active_masks = (int64_t *)malloc(total_count * sizeof(int64_t));
		rt_alpha->p_ids = (int32_t *)malloc(total_count * sizeof(int32_t));
		for (int64_t im = 0; im < num_states; im++) {
			for (int64_t i = 0; i < n_par; i++) {
				int64_t idx = im * n_par + i;
				rt_alpha->active_masks[idx] = im;
				rt_alpha->p_ids[idx] = E_ID(im) + (int32_t)i;
			}
		}

		LAGraph_rule_EWCNF rule_alpha = {0};
		rule_alpha.nonterm = NT_START;
		rule_alpha.prod_A = E_base;
		rule_alpha.prod_B = NT_START;
		rule_alpha.xor_routing = rt_alpha;
		rules[r++] = rule_alpha;
	}

// Cleanup macros
#undef E_ID
#undef G_ID

	return gr;
}

MRGrammar_t dyck_alpha_grammar_k_parity_se(int64_t n_par, int64_t n_bra,
										   bool has_normal, int64_t k) {
	int64_t num_parity_states = (int64_t)1 << k;
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
	// S (start):                1
	// Eps:                      1
	// S_{mask, qStart, qEnd}:   s_count
	// A_i, B_i (par aux):       2 * n_par
	// C_i, D_i (bra aux):       2 * n_bra
	// E_{im,qs,qmid,i}, G_{im,qs,qmid,i}: for each par i, for each (im,qs,qmid)
	//   triple, we need E and G.
	//   Total: 2 * n_par * num_parity_states * RSTATE_COUNT^2 (im, qs, qmid)
	// N (if has_normal):        1
	int64_t eg_count = n_par * num_parity_states * RSTATE_COUNT * RSTATE_COUNT;

	int64_t nonterms_count =
		2 + s_count + 2 * n_par + 2 * n_bra + 2 * eg_count + (has_normal ? 1 : 0);

	int64_t terms_count = 2 * n_par + 2 * n_bra + (has_normal ? 1 : 0);

	// Rules count:
	// Eps -> eps:                                        1
	// S -> S_{0,QE,QE} Eps, S -> S_{0,QE,QC} Eps:        2
	// S_{0,q,q} -> eps (for each q):                     RSTATE_COUNT
	// N -> `normal`, S_{m,qs,qe} -> N S_{m,qs,qe}:       1 + s_count (if has_normal)
	// C_i -> `[_i`, D_i -> `]_i`:                        2 * n_bra
	// S_{m,qs,qe} -> C_i S_{next,QO,qe}:                 active_groups *
	//                                                       RSTATE_COUNT *
	//                                                       RSTATE_COUNT (XOR OPT)
	// S_{m,qs,qe} -> D_i S_{next,QC,qe} (qs!=QE):        active_groups *
	//                                                       RSTATE_COUNT *
	//                                                       (RSTATE_COUNT - 1)
	//                                                       (XOR OPT)
	// A_i -> `(_i`, B_i -> `)_i`:                        2
	// E_{im,qs,qmid,i} -> A_i G_{...}:                   eg_base_count
	// G_{im,qs,qmid,i} -> S_{im,qs,qmid} B_i:            eg_base_count
	// S_{m,qs,qe} -> E_{im,qs,qmid,i} S_{next,qmid,qe}:  RSTATE_COUNT *
	//                                                       RSTATE_COUNT *
	//                                                       RSTATE_COUNT
	//                                                       (XOR OPT)
	int64_t bra_close_count =
		active_groups * num_parity_states * RSTATE_COUNT * (RSTATE_COUNT - 1);
	int64_t par_eg_count = s_count * num_parity_states * RSTATE_COUNT;
	int64_t eg_base_count = num_parity_states * RSTATE_COUNT * RSTATE_COUNT;
	int64_t rules_to_alloc =
		3 + RSTATE_COUNT + (has_normal ? s_count + 1 : 0) + 2 * active_groups +
		(n_par > 0 ? 2 : 0) + active_groups * RSTATE_COUNT * RSTATE_COUNT +
		active_groups * RSTATE_COUNT * (RSTATE_COUNT - 1) +
		(n_par > 0 ? (2 * eg_base_count + RSTATE_COUNT * RSTATE_COUNT * RSTATE_COUNT)
				   : 0);

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
	int32_t NT_N = (has_normal) ? (int32_t)(nonterms_count - 1) : -1;

	// --- Nonterminal IDs helpers ---

#define NT_S(mask, qs, qe)                                                          \
	((int32_t)(S_base + (mask) * RSTATE_COUNT * RSTATE_COUNT +                      \
			   (qs) * RSTATE_COUNT + (qe)))

#define EG_IDX(im, qs, qmid)                                                        \
	((int64_t)(im) * RSTATE_COUNT * RSTATE_COUNT * n_par +                          \
	 (int64_t)(qs) * RSTATE_COUNT * n_par + (int64_t)(qmid) * n_par)

#define E_ID(im, qs, qmid) (E_base + (int32_t)EG_IDX(im, qs, qmid))
#define G_ID(im, qs, qmid) (G_base + (int32_t)EG_IDX(im, qs, qmid))

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
													  NT_S(m, qs, qe), 0, 0};
				}
			}
		}
	}

	// --- Terminal rules ---
	if (n_par > 0) {
		// A_i -> `(_i`
		rules[r++] = (LAGraph_rule_EWCNF){
			A_0, TERM(t_par_open_0), -1, (uint32_t)n_par,
			LAGraph_EWNCF_INDEX_NONTERM | LAGraph_EWNCF_INDEX_PROD_A};
		// B_i -> `)_i`
		rules[r++] = (LAGraph_rule_EWCNF){
			B_0, TERM(t_par_close_0), -1, (uint32_t)n_par,
			LAGraph_EWNCF_INDEX_NONTERM | LAGraph_EWNCF_INDEX_PROD_A};
	}

	// Terminal rules for brackets by groups
	for (int64_t b = 0; b < k; b++) {
		if (bra_group_size[b] == 0) {
			continue;
		}

		// C_{b, j} -> `[_{b, j}`
		rules[r++] = (LAGraph_rule_EWCNF){
			C_group_base[b], TERM(t_bra_open_group_base[b]), -1,
			(uint32_t)bra_group_size[b],
			LAGraph_EWNCF_INDEX_NONTERM | LAGraph_EWNCF_INDEX_PROD_A};
		// D_{b, j} -> `]_{b, j}`
		rules[r++] = (LAGraph_rule_EWCNF){
			D_group_base[b], TERM(t_bra_close_group_base[b]), -1,
			(uint32_t)bra_group_size[b],
			LAGraph_EWNCF_INDEX_NONTERM | LAGraph_EWNCF_INDEX_PROD_A};
	}

	// --- Beta transitions rules ---
	for (int64_t b = 0; b < k; b++) {
		if (bra_group_size[b] == 0) {
			continue;
		}
		int64_t bit = (int64_t)1 << b;
		int64_t cnt = bra_group_size[b];

		for (int qs = 0; qs < RSTATE_COUNT; qs++) {
			for (int qe = 0; qe < RSTATE_COUNT; qe++) {
				// S_{m,qs,qe} -> C_{b,j} S_{next,QO,qe}  (open bra, always allowed)
				XorFamilyRouting *rt_open = calloc(1, sizeof(XorFamilyRouting));
				rt_open->S_operand_base = NT_S(0, RSTATE_QO, qe);
				rt_open->S_target_base = NT_S(0, qs, qe);
				rt_open->S_stride = RSTATE_COUNT * RSTATE_COUNT;
				rt_open->N = num_parity_states;
				rt_open->active_masks_count = cnt;
				rt_open->active_masks = malloc(cnt * sizeof(int64_t));
				rt_open->p_ids = malloc(cnt * sizeof(int32_t));
				for (int64_t j = 0; j < cnt; j++) {
					rt_open->active_masks[j] = bit;
					rt_open->p_ids[j] = C_group_base[b] + (int32_t)j;
				}

				LAGraph_rule_EWCNF rule_open = {0};
				rule_open.nonterm = NT_START;
				rule_open.prod_A = C_group_base[b];
				rule_open.prod_B = NT_START;
				rule_open.xor_routing = rt_open;
				rules[r++] = rule_open;

				// S_{m,qs,qe} -> D_{b,j} S_{next,QC,qe}  (close bra, forbidden from
				//   QE)
				if (qs != RSTATE_QE) {
					XorFamilyRouting *rt_close = calloc(1, sizeof(XorFamilyRouting));
					rt_close->S_operand_base = NT_S(0, RSTATE_QC, qe);
					rt_close->S_target_base = NT_S(0, qs, qe);
					rt_close->S_stride = RSTATE_COUNT * RSTATE_COUNT;
					rt_close->N = num_parity_states;
					rt_close->active_masks_count = cnt;
					rt_close->active_masks = malloc(cnt * sizeof(int64_t));
					rt_close->p_ids = malloc(cnt * sizeof(int32_t));
					for (int64_t j = 0; j < cnt; j++) {
						rt_close->active_masks[j] = bit;
						rt_close->p_ids[j] = D_group_base[b] + (int32_t)j;
					}

					LAGraph_rule_EWCNF rule_close = {0};
					rule_close.nonterm = NT_START;
					rule_close.prod_A = D_group_base[b];
					rule_close.prod_B = NT_START;
					rule_close.xor_routing = rt_close;
					rules[r++] = rule_close;
				}
			}
		}
	}

	// --- Alpha transitions rules ---
	if (n_par > 0) {
		// E_{im,qs,qmid,i} -> A_i G_{im,qs,qmid,i}
		for (int64_t im = 0; im < num_parity_states; im++) {
			for (int qs = 0; qs < RSTATE_COUNT; qs++) {
				for (int qmid = 0; qmid < RSTATE_COUNT; qmid++) {
					rules[r++] = (LAGraph_rule_EWCNF){
						E_ID(im, qs, qmid), A_0, G_ID(im, qs, qmid), (uint32_t)n_par,
						LAGraph_EWNCF_INDEX_NONTERM | LAGraph_EWNCF_INDEX_PROD_A |
							LAGraph_EWNCF_INDEX_PROD_B};
				}
			}
		}

		// G_{im,qs,qmid,i} -> S_{im,qs,qmid} B_i
		for (int64_t im = 0; im < num_parity_states; im++) {
			for (int qs = 0; qs < RSTATE_COUNT; qs++) {
				for (int qmid = 0; qmid < RSTATE_COUNT; qmid++) {
					rules[r++] = (LAGraph_rule_EWCNF){
						G_ID(im, qs, qmid), NT_S(im, qs, qmid), B_0, (uint32_t)n_par,
						LAGraph_EWNCF_INDEX_NONTERM | LAGraph_EWNCF_INDEX_PROD_B};
				}
			}
		}

		// S_{m,qs,qe} -> E_{im,qs,qmid,i} S_{next,qmid,qe}
		for (int qs = 0; qs < RSTATE_COUNT; qs++) {
			for (int qmid = 0; qmid < RSTATE_COUNT; qmid++) {
				for (int qe = 0; qe < RSTATE_COUNT; qe++) {
					int64_t total_count = num_parity_states * n_par;

					XorFamilyRouting *rt = calloc(1, sizeof(XorFamilyRouting));
					rt->S_operand_base = NT_S(0, qmid, qe);
					rt->S_target_base = NT_S(0, qs, qe);
					rt->S_stride = RSTATE_COUNT * RSTATE_COUNT;
					rt->N = num_parity_states;
					rt->active_masks_count = total_count;
					rt->active_masks = malloc(total_count * sizeof(int64_t));
					rt->p_ids = malloc(total_count * sizeof(int32_t));
					for (int64_t im = 0; im < num_parity_states; im++) {
						for (int64_t i = 0; i < n_par; i++) {
							int64_t idx = im * n_par + i;
							rt->active_masks[idx] = im;
							rt->p_ids[idx] = E_ID(im, qs, qmid) + (int32_t)i;
						}
					}

					LAGraph_rule_EWCNF rule = {0};
					rule.nonterm = NT_START;
					rule.prod_A = E_base;
					rule.prod_B = NT_START;
					rule.xor_routing = rt;
					rules[r++] = rule;
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

MRGrammar_t dyck_alpha_grammar_k_parity_exclude(int64_t n_par, int64_t n_bra,
												bool has_normal, int64_t k,
												int64_t ex_bra) {
	int64_t num_states = (int64_t)1 << k;

	// Arrays for storing the sizes of brackets and their base IDs.
	// Size is 65 to safely accommodate k normal groups + 1 exception group
	int32_t bra_group_size[65] = {0};
	int32_t C_group_base[65] = {0};
	int32_t D_group_base[65] = {0};
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
	int64_t nonterms_count = num_states + 2 * n_par + 2 * n_bra +
							 2 * (n_par * num_states) + (has_normal ? 1 : 0);

	int64_t terms_count = 2 * n_par + 2 * n_bra + (has_normal ? 1 : 0);

	// Rules count (same layout/reasoning as dyck_alpha_grammar_k_parity):
	// S_0 -> eps:                         1
	// S_m -> N S_m (for each m):          num_states (if has_normal)
	// N -> `normal`:                      1          (if has_normal)
	// Terminal rules for A_i, B_i:        2          (if n_par>0)
	// Terminal rules for C, D by groups:  2 * active_groups
	// Beta transitions by groups:         2 * active_groups (XOR OPT)
	// Alpha transitions (if n_par>0):
	//   E_{im,i} -> A_i G_{im,i}:         num_states
	//   G_{im,i} -> S_im B_i:             num_states
	//   S_m -> E_{im,i} S_next:           1 (XOR OPT)
	//   Total: 1 + 2 * num_states
	int64_t rules_to_alloc = 1 + (has_normal ? num_states + 1 : 0) +
							 (n_par > 0 ? 2 : 0) + 2 * active_groups +
							 2 * active_groups +
							 (n_par > 0 ? (1 + 2 * num_states) : 0);

	MRGrammarConfig config = {.kind = MR_GRAMMAR_ALPHA, .exclude_index = ex_bra};
	MRGrammar_t gr =
		make_grammar(rules_to_alloc, terms_count, nonterms_count, config, k);
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

	// Label nonterminals C and D with continuous groups (including exception group)
	for (int b = 0; b <= k; b++) {
		C_group_base[b] = (bra_group_size[b] > 0) ? offset : -1;
		offset += bra_group_size[b];
	}
	for (int b = 0; b <= k; b++) {
		D_group_base[b] = (bra_group_size[b] > 0) ? offset : -1;
		offset += bra_group_size[b];
	}

	int32_t E_base = offset;
	offset += (int32_t)(n_par * num_states);
	int32_t G_base = offset;
	offset += (int32_t)(n_par * num_states);

#define E_ID(im) (E_base + (int32_t)((im) * n_par))
#define G_ID(im) (G_base + (int32_t)((im) * n_par))

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

	// --- Rules ---

	// S -> eps
	rules[r++] = (LAGraph_rule_EWCNF){NT_START, -1, -1, 0, 0};

	if (has_normal) {
		// N -> `normal`
		rules[r++] = (LAGraph_rule_EWCNF){NT_N, TERM(t_normal_idx), -1, 0, 0};
		// S_m -> N S_m
		for (int64_t m = 0; m < num_states; m++) {
			rules[r++] = (LAGraph_rule_EWCNF){(int32_t)m, NT_N, (int32_t)m, 0, 0};
		}
	}

	// --- Terminal rules ---
	if (n_par > 0) {
		rules[r++] = (LAGraph_rule_EWCNF){
			A_0, TERM(t_par_open_0), -1, (uint32_t)n_par,
			LAGraph_EWNCF_INDEX_NONTERM | LAGraph_EWNCF_INDEX_PROD_A};
		rules[r++] = (LAGraph_rule_EWCNF){
			B_0, TERM(t_par_close_0), -1, (uint32_t)n_par,
			LAGraph_EWNCF_INDEX_NONTERM | LAGraph_EWNCF_INDEX_PROD_A};
	}

	// Terminal rules for brackets by groups
	for (int b = 0; b <= k; b++) {
		if (bra_group_size[b] == 0) {
			continue;
		}

		// C_{b, j} -> `[_{b, j}`
		rules[r++] = (LAGraph_rule_EWCNF){
			C_group_base[b], TERM(t_bra_open_group_base[b]), -1,
			(uint32_t)bra_group_size[b],
			LAGraph_EWNCF_INDEX_NONTERM | LAGraph_EWNCF_INDEX_PROD_A};
		// D_{b, j} -> `]_{b, j}`
		rules[r++] = (LAGraph_rule_EWCNF){
			D_group_base[b], TERM(t_bra_close_group_base[b]), -1,
			(uint32_t)bra_group_size[b],
			LAGraph_EWNCF_INDEX_NONTERM | LAGraph_EWNCF_INDEX_PROD_A};
	}

	// --- Beta transition rules ---
	for (int b = 0; b <= k; b++) {
		if (bra_group_size[b] == 0) {
			continue;
		}

		// If it's the exception group (b == k), bit change is 0. Otherwise 1 << b.
		int64_t bit = (b == k) ? 0 : (int64_t)1 << b;
		int64_t cnt = bra_group_size[b];

		// S_m -> C_{b, j} S_next
		XorFamilyRouting *rt_open =
			(XorFamilyRouting *)calloc(1, sizeof(XorFamilyRouting));
		rt_open->S_operand_base = NT_START;
		rt_open->S_target_base = NT_START;
		rt_open->S_stride = 1;
		rt_open->N = num_states;
		rt_open->active_masks_count = cnt;
		rt_open->active_masks = (int64_t *)malloc(cnt * sizeof(int64_t));
		rt_open->p_ids = (int32_t *)malloc(cnt * sizeof(int32_t));
		for (int64_t j = 0; j < cnt; j++) {
			rt_open->active_masks[j] = bit;
			rt_open->p_ids[j] = C_group_base[b] + (int32_t)j;
		}

		LAGraph_rule_EWCNF rule_open = {0};
		rule_open.nonterm = NT_START;
		rule_open.prod_A = C_group_base[b];
		rule_open.prod_B = NT_START;
		rule_open.xor_routing = rt_open;
		rules[r++] = rule_open;

		// S_m -> D_{b, j} S_next
		XorFamilyRouting *rt_close =
			(XorFamilyRouting *)calloc(1, sizeof(XorFamilyRouting));
		rt_close->S_operand_base = NT_START;
		rt_close->S_target_base = NT_START;
		rt_close->S_stride = 1;
		rt_close->N = num_states;
		rt_close->active_masks_count = cnt;
		rt_close->active_masks = (int64_t *)malloc(cnt * sizeof(int64_t));
		rt_close->p_ids = (int32_t *)malloc(cnt * sizeof(int32_t));
		for (int64_t j = 0; j < cnt; j++) {
			rt_close->active_masks[j] = bit;
			rt_close->p_ids[j] = D_group_base[b] + (int32_t)j;
		}

		LAGraph_rule_EWCNF rule_close = {0};
		rule_close.nonterm = NT_START;
		rule_close.prod_A = D_group_base[b];
		rule_close.prod_B = NT_START;
		rule_close.xor_routing = rt_close;
		rules[r++] = rule_close;
	}

	// --- Alpha transition rules ---
	if (n_par > 0) {
		// E_{im,i} -> A_i G_{im,i}
		for (int64_t im = 0; im < num_states; im++) {
			rules[r++] = (LAGraph_rule_EWCNF){
				E_ID(im), A_0, G_ID(im), (uint32_t)n_par,
				LAGraph_EWNCF_INDEX_NONTERM | LAGraph_EWNCF_INDEX_PROD_A |
					LAGraph_EWNCF_INDEX_PROD_B};
		}

		// G_{im,i} -> S_im B_i
		for (int64_t im = 0; im < num_states; im++) {
			rules[r++] = (LAGraph_rule_EWCNF){
				G_ID(im), NT_START + (int32_t)im, B_0, (uint32_t)n_par,
				LAGraph_EWNCF_INDEX_NONTERM | LAGraph_EWNCF_INDEX_PROD_B};
		}

		int64_t total_count = num_states * n_par;

		// S_m -> E_{im,i} S_next
		XorFamilyRouting *rt_alpha =
			(XorFamilyRouting *)calloc(1, sizeof(XorFamilyRouting));
		rt_alpha->S_operand_base = NT_START;
		rt_alpha->S_target_base = NT_START;
		rt_alpha->S_stride = 1;
		rt_alpha->N = num_states;
		rt_alpha->active_masks_count = total_count;
		rt_alpha->active_masks = (int64_t *)malloc(total_count * sizeof(int64_t));
		rt_alpha->p_ids = (int32_t *)malloc(total_count * sizeof(int32_t));
		for (int64_t im = 0; im < num_states; im++) {
			for (int64_t i = 0; i < n_par; i++) {
				int64_t idx = im * n_par + i;
				rt_alpha->active_masks[idx] = im;
				rt_alpha->p_ids[idx] = E_ID(im) + (int32_t)i;
			}
		}

		LAGraph_rule_EWCNF rule_alpha = {0};
		rule_alpha.nonterm = NT_START;
		rule_alpha.prod_A = E_base;
		rule_alpha.prod_B = NT_START;
		rule_alpha.xor_routing = rt_alpha;
		rules[r++] = rule_alpha;
	}

#undef E_ID
#undef G_ID

	return gr;
}
