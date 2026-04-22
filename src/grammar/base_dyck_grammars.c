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
	int64_t rules_count = 1 + 5 * n_par + 5 * n_bra + (has_normal ? 2 : 0);

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

	// Adding rules for each pair of square brackets
	for (int64_t i = 0; i < n_bra; i++) {
		int32_t C_i = (int32_t)(1 + 4 * n_par + 4 * i);
		int32_t D_i = (int32_t)(1 + 4 * n_par + 4 * i + 1);
		int32_t F_i = (int32_t)(1 + 4 * n_par + 4 * i + 2);
		int32_t H_i = (int32_t)(1 + 4 * n_par + 4 * i + 3);

		int32_t bra_open = (int32_t)(2 * n_par + 2 * i);
		int32_t bra_close = (int32_t)(2 * n_par + 2 * i + 1);

		// S -> C_i F_i
		rules[r++] = (LAGraph_rule_WCNF){NT_START, C_i, F_i, 0};
		// F_i -> S H_i
		rules[r++] = (LAGraph_rule_WCNF){F_i, NT_START, H_i, 0};
		// H_i -> D_i S
		rules[r++] = (LAGraph_rule_WCNF){H_i, D_i, NT_START, 0};
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
	// AnyOpen  -> open_i  (par + bra):  n_par + n_bra
	// AnyClose -> close_j (par + bra):  n_par + n_bra
	int64_t rules_count = 4 + (has_normal ? 2 : 0) + 2 * (n_par + n_bra);

	LAGraph_rule_WCNF *rules = malloc(rules_count * sizeof(LAGraph_rule_WCNF));
	int64_t r = 0;

	int32_t NT_E = 1;
	int32_t NT_G = 2;
	int32_t NT_ANY_OPEN = 3;
	int32_t NT_ANY_CLOSE = 4;
	int32_t NT_N = 5; // only if has_normal

#define TERM_PAR_OPEN(i) ((int32_t)(2 * (i)))
#define TERM_PAR_CLOSE(i) ((int32_t)(2 * (i) + 1))
#define TERM_BRA_OPEN(i) ((int32_t)(2 * n_par + 2 * (i)))
#define TERM_BRA_CLOSE(i) ((int32_t)(2 * n_par + 2 * (i) + 1))

	// S -> eps
	rules[r++] = (LAGraph_rule_WCNF){NT_START, -1, -1, 0};

	if (has_normal) {
		// N -> `normal`
		rules[r++] = (LAGraph_rule_WCNF){NT_N, (int32_t)(terms_count - 1), -1, 0};
		// S -> N S
		rules[r++] = (LAGraph_rule_WCNF){NT_START, NT_N, NT_START, 0};
	}

	// S -> AnyOpen E
	rules[r++] = (LAGraph_rule_WCNF){NT_START, NT_ANY_OPEN, NT_E, 0};
	// E -> S G
	rules[r++] = (LAGraph_rule_WCNF){NT_E, NT_START, NT_G, 0};
	// G -> AnyClose S
	rules[r++] = (LAGraph_rule_WCNF){NT_G, NT_ANY_CLOSE, NT_START, 0};

	for (int64_t i = 0; i < n_par; i++) {
		rules[r++] = (LAGraph_rule_WCNF){NT_ANY_OPEN, TERM_PAR_OPEN(i), -1, 0};
		rules[r++] = (LAGraph_rule_WCNF){NT_ANY_CLOSE, TERM_PAR_CLOSE(i), -1, 0};
	}
	for (int64_t i = 0; i < n_bra; i++) {
		rules[r++] = (LAGraph_rule_WCNF){NT_ANY_OPEN, TERM_BRA_OPEN(i), -1, 0};
		rules[r++] = (LAGraph_rule_WCNF){NT_ANY_CLOSE, TERM_BRA_CLOSE(i), -1, 0};
	}

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

void grammar_free(MRGrammar_t *gr) {
	free(gr->rules);
	gr->rules = NULL;
	gr->rules_count = 0;
	gr->terms_count = 0;
	gr->nonterms_count = 0;
}
