#include "grammar.h"

/**
 * S -> `(_i` S `)_i` S | `[_i` S `]_i` S | `normal` S | eps in WCNF :
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
 * S -> C_i F_i
 * F_i -> S H_i
 * H_i -> D_i S
 * C_i -> `[_i`
 * D_i -> `]_i`
 */
MRGrammar_t dyck_grammar(int64_t n_par, int64_t n_bra, bool has_normal) {
	int64_t nonterms_count = 1 + 4 * n_par + 4 * n_bra + (has_normal ? 1 : 0);
	int64_t terms_count = 2 * n_par + 2 * n_bra + (has_normal ? 1 : 0);

	int64_t rules_count = 1; // S -> eps
	if (has_normal) {
		rules_count += 2;
	}
	if (n_par > 0) {
		rules_count += 5;
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

	// Allocation of Nonterminal IDs

	int32_t offset = 1;

	int32_t A_0 = (n_par > 0) ? offset : -1;
	offset += (int32_t)n_par;
	int32_t B_0 = (n_par > 0) ? offset : -1;
	offset += (int32_t)n_par;
	int32_t E_0 = (n_par > 0) ? offset : -1;
	offset += (int32_t)n_par;
	int32_t G_0 = (n_par > 0) ? offset : -1;
	offset += (int32_t)n_par;

	int32_t C_0 = (n_bra > 0) ? offset : -1;
	offset += (int32_t)n_bra;
	int32_t D_0 = (n_bra > 0) ? offset : -1;
	offset += (int32_t)n_bra;
	int32_t F_0 = (n_bra > 0) ? offset : -1;
	offset += (int32_t)n_bra;
	int32_t H_0 = (n_bra > 0) ? offset : -1;
	offset += (int32_t)n_bra;

	int32_t NT_N = has_normal ? (int32_t)(nonterms_count - 1) : -1;

	// Terminal IDs

	int32_t t_idx = 0;
	int32_t t_par_open_0 = t_idx;
	t_idx += (int32_t)n_par;
	int32_t t_par_close_0 = t_idx;
	t_idx += (int32_t)n_par;
	int32_t t_bra_open_0 = t_idx;
	t_idx += (int32_t)n_bra;
	int32_t t_bra_close_0 = t_idx;
	t_idx += (int32_t)n_bra;

	int32_t t_normal_idx = has_normal ? t_idx : -1;

	// Rules

	// S -> eps
	rules[r++] = (LAGraph_rule_EWCNF){NT_START, -1, -1, 0, 0};

	if (has_normal) {
		// S -> N S
		rules[r++] = (LAGraph_rule_EWCNF){NT_START, NT_N, NT_START, 0, 0};
		// N -> normal
		rules[r++] = (LAGraph_rule_EWCNF){NT_N, TERM(t_normal_idx), -1, 0, 0};
	}

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
		// A_i -> (_i
		rules[r++] = (LAGraph_rule_EWCNF){A_0, TERM(t_par_open_0), -1, count,
										  LAGraph_EWNCF_INDEX_NONTERM |
											  LAGraph_EWNCF_INDEX_PROD_A};
		// B_i -> )_i
		rules[r++] = (LAGraph_rule_EWCNF){B_0, TERM(t_par_close_0), -1, count,
										  LAGraph_EWNCF_INDEX_NONTERM |
											  LAGraph_EWNCF_INDEX_PROD_A};
	}

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
		// C_i -> [_i
		rules[r++] = (LAGraph_rule_EWCNF){C_0, TERM(t_bra_open_0), -1, count,
										  LAGraph_EWNCF_INDEX_NONTERM |
											  LAGraph_EWNCF_INDEX_PROD_A};
		// D_i -> ]_i
		rules[r++] = (LAGraph_rule_EWCNF){D_0, TERM(t_bra_close_0), -1, count,
										  LAGraph_EWNCF_INDEX_NONTERM |
											  LAGraph_EWNCF_INDEX_PROD_A};
	}

	return gr;
}

MRGrammar_t dyck_project_grammar(int64_t n_par, int64_t n_bra, bool has_normal) {
	// Nonterminals: S, E, G, AnyOpen, AnyClose, [N]
	int64_t nonterms_count = 5 + (has_normal ? 1 : 0);
	int64_t terms_count = 2 * n_par + 2 * n_bra + (has_normal ? 1 : 0);

	// Rules:
	// S -> eps:                         1
	// S -> N S, N -> `normal`:          2  (if has_normal)
	// S -> AnyOpen E:                   1
	// E -> S G:                         1
	// G -> AnyClose S:                  1
	// AnyOpen  -> open_i  (par + bra):  2
	// AnyClose -> close_j (par + bra):  2
	int64_t rules_count = 4 + (has_normal ? 2 : 0);
	if (n_par > 0) {
		rules_count += 2;
	}
	if (n_bra > 0) {
		rules_count += 2;
	}

	MRGrammar_t gr = make_grammar(rules_count, terms_count, nonterms_count);
	if (!gr.rules) {
		return gr;
	}

	LAGraph_rule_EWCNF *rules = gr.rules;
	int64_t r = 0;

	// Nonterminal IDs

	int32_t NT_E = 1;
	int32_t NT_G = 2;
	int32_t NT_ANY_OPEN = 3;
	int32_t NT_ANY_CLOSE = 4;
	int32_t NT_N = 5; // if has_normal

	// Terminal offset

	int32_t t_idx = 0;
	int32_t t_par_open_0 = t_idx;
	t_idx += (int32_t)n_par;
	int32_t t_par_close_0 = t_idx;
	t_idx += (int32_t)n_par;
	int32_t t_bra_open_0 = t_idx;
	t_idx += (int32_t)n_bra;
	int32_t t_bra_close_0 = t_idx;
	t_idx += (int32_t)n_bra;
	int32_t t_normal_idx = has_normal ? t_idx : -1;

	// S -> eps
	rules[r++] = (LAGraph_rule_EWCNF){NT_START, -1, -1, 0, 0};

	if (has_normal) {
		// S -> N S
		rules[r++] = (LAGraph_rule_EWCNF){NT_START, NT_N, NT_START, 0, 0};
		// N -> `normal`
		rules[r++] = (LAGraph_rule_EWCNF){NT_N, TERM(t_normal_idx), -1, 0, 0};
	}

	// S -> AnyOpen E
	rules[r++] = (LAGraph_rule_EWCNF){NT_START, NT_ANY_OPEN, NT_E, 0, 0};
	// E -> S G
	rules[r++] = (LAGraph_rule_EWCNF){NT_E, NT_START, NT_G, 0, 0};
	// G -> AnyClose S
	rules[r++] = (LAGraph_rule_EWCNF){NT_G, NT_ANY_CLOSE, NT_START, 0, 0};

	// Rules for terminals

	if (n_par > 0) {
		uint32_t count = (uint32_t)n_par;
		// AnyOpen -> (_i
		rules[r++] = (LAGraph_rule_EWCNF){NT_ANY_OPEN, TERM(t_par_open_0), -1, count,
										  LAGraph_EWNCF_INDEX_PROD_A};
		// AnyClose -> )_i
		rules[r++] = (LAGraph_rule_EWCNF){NT_ANY_CLOSE, TERM(t_par_close_0), -1,
										  count, LAGraph_EWNCF_INDEX_PROD_A};
	}

	if (n_bra > 0) {
		uint32_t count = (uint32_t)n_bra;
		// AnyOpen -> [_i
		rules[r++] = (LAGraph_rule_EWCNF){NT_ANY_OPEN, TERM(t_bra_open_0), -1, count,
										  LAGraph_EWNCF_INDEX_PROD_A};
		// AnyClose -> ]_i
		rules[r++] = (LAGraph_rule_EWCNF){NT_ANY_CLOSE, TERM(t_bra_close_0), -1,
										  count, LAGraph_EWNCF_INDEX_PROD_A};
	}

	return gr;
}
