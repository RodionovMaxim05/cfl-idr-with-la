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

int64_t MR_grammar_get_group_size(const MRGrammarConfig *config, int64_t b,
								  int64_t n, int64_t k) {
	bool has_exclude = (config->exclude_index >= 0);
	if (b == k) {
		return has_exclude ? 1 : 0;
	}

	int64_t size = (b < n) ? (n - b + k - 1) / k : 0;

	// If this group suffered from the removal of a parenthesis, reduce its size
	if (has_exclude && b == (config->exclude_index % k)) {
		size--;
	}
	return size;
}

void MR_grammar_get_bracket_layout(const MRGrammarConfig *config, int64_t i,
								   int64_t n, int64_t k, int64_t *out_group,
								   int64_t *out_sub_idx) {
	// If this is the same excluded parenthesis, send it to the special group 'k'
	if (config->exclude_index >= 0 && i == config->exclude_index) {
		*out_group = k;
		*out_sub_idx = 0;
		return;
	}

	int64_t g = i % k;
	int64_t j = i / k;

	// Collapse the hole if the current element comes after the excluded one in the
	// same group
	if (config->exclude_index >= 0 && i > config->exclude_index &&
		(i % k == config->exclude_index % k)) {
		j--;
	}

	*out_group = g;
	*out_sub_idx = j;
}
