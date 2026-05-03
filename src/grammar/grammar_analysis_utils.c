#include "grammar_analysis_utils.h"

MRGrammar_t get_alpha_grammar(IdrGrammarType curGrammar, int64_t n_par,
							  int64_t n_bra, bool has_normal) {
	switch (curGrammar) {
		case IDR_PARITY:
			return dyck_alpha_grammar_k_parity(n_par, n_bra, has_normal, 1);
		case IDR_PARITY2:
			return dyck_alpha_grammar_k_parity(n_par, n_bra, has_normal, 2);
		case IDR_SE:
			return dyck_alpha_grammar_k_parity_se(n_par, n_bra, has_normal, 2);
		case IDR_PROJECT:
			return dyck_alpha_grammar_k_parity(n_par, n_bra, has_normal, 1);
		case IDR_EXCLUDE:
			return dyck_alpha_grammar_k_parity(n_par, n_bra, has_normal, 1);
		case IDR_ALL:
			return dyck_alpha_grammar_k_parity_se(n_par, n_bra, has_normal, 2);
		default:
			return dyck_alpha_grammar(n_par, n_bra, has_normal);
	}
}

MRGrammar_t get_beta_grammar(IdrGrammarType curGrammar, int64_t n_par, int64_t n_bra,
							 bool has_normal) {
	switch (curGrammar) {
		case IDR_PARITY:
			return dyck_beta_grammar_k_parity(n_par, n_bra, has_normal, 1);
		case IDR_PARITY2:
			return dyck_beta_grammar_k_parity(n_par, n_bra, has_normal, 2);
		case IDR_SE:
			return dyck_beta_grammar_k_parity_se(n_par, n_bra, has_normal, 2);
		case IDR_PROJECT:
			return dyck_beta_grammar_k_parity(n_par, n_bra, has_normal, 1);
		case IDR_EXCLUDE:
			return dyck_beta_grammar_k_parity(n_par, n_bra, has_normal, 1);
		case IDR_ALL:
			return dyck_beta_grammar_k_parity_se(n_par, n_bra, has_normal, 2);
		default:
			return dyck_beta_grammar(n_par, n_bra, has_normal);
	}
}
