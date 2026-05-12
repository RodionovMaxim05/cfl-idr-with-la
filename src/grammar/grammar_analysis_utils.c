#include "grammar_analysis_utils.h"

#include "grammar.h"

MRGrammar_t get_dyck_grammar(int64_t n_par, int64_t n_bra, bool has_normal) {
	return dyck_grammar(n_par, n_bra, has_normal);
}

MRGrammar_t get_project_grammar(int64_t n_par, int64_t n_bra, bool has_normal) {
	return dyck_project_grammar(n_par, n_bra, has_normal);
}

MRGrammar_t get_exclude_grammar(int64_t n_par, int64_t n_bra, bool has_normal,
								int64_t ex_bra) {
	return dyck_alpha_grammar_k_parity_exclude(n_par, n_bra, has_normal,
											   /*k=*/2, ex_bra);
}

MRGrammar_t get_alpha_grammar(IdrGrammarType grammar_type, int64_t n_par,
							  int64_t n_bra, bool has_normal) {
	switch (grammar_type) {
		case IDR_PARITY:
		case IDR_PROJECT:
		case IDR_EXCLUDE:
			return dyck_alpha_grammar_k_parity(n_par, n_bra, has_normal, 1);
		case IDR_PARITY2:
			return dyck_alpha_grammar_k_parity(n_par, n_bra, has_normal, 2);
		case IDR_SE:
		case IDR_ALL:
			return dyck_alpha_grammar_k_parity_se(n_par, n_bra, has_normal, 2);
		default:
			return dyck_alpha_grammar(n_par, n_bra, has_normal);
	}
}

MRGrammar_t get_beta_grammar(IdrGrammarType grammar_type, int64_t n_par,
							 int64_t n_bra, bool has_normal) {
	switch (grammar_type) {
		case IDR_PARITY:
		case IDR_PROJECT:
		case IDR_EXCLUDE:
			return dyck_beta_grammar_k_parity(n_par, n_bra, has_normal, 1);
		case IDR_PARITY2:
			return dyck_beta_grammar_k_parity(n_par, n_bra, has_normal, 2);
		case IDR_SE:
		case IDR_ALL:
			return dyck_beta_grammar_k_parity_se(n_par, n_bra, has_normal, 2);
		default:
			return dyck_beta_grammar(n_par, n_bra, has_normal);
	}
}

void grammar_free(MRGrammar_t *gr) {
	if (!gr) {
		return;
	}
	free(gr->rules);
	gr->rules = NULL;
	gr->rules_count = 0;
	gr->terms_count = 0;
	gr->nonterms_count = 0;
}
