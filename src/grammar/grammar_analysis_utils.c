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

static void pat_free(LAGraph_rule_Pat *pat) {
	if (!pat) {
		return;
	}

	if (pat->state) {
		free((void *)pat->state->entries);
		free((void *)pat->state);
	}
	free(pat);
}

void grammar_free(MRGrammar_t *gr) {
	if (!gr || !gr->rules) {
		return;
	}

	for (int64_t i = 0; i < gr->rules_count; i++) {
		LAGraph_rule_Pat *pat = (LAGraph_rule_Pat *)gr->rules[i].pat;
		if (pat != NULL) {
			for (int64_t j = i + 1; j < gr->rules_count; j++) {
				if (gr->rules[j].pat == pat) {
					gr->rules[j].pat = NULL;
				}
			}
			pat_free(pat);
			gr->rules[i].pat = NULL;
		}
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

/**
 * @brief Fills the interleaved mapping table for grammar terminals to standard
 * terminal order.
 *
 * @param[out] map         Output mapping array.
 * @param[in]  n           Number of symbols.
 * @param[in]  k           Interleaving factor/k-parameter.
 * @param[in]  start_idx   Offset index in output map.
 * @param[in]  out_base    Base output index.
 * @param[in]  config      Grammar configuration.
 */
static void fill_interleaved_map(int32_t *map, int64_t n, int64_t k,
								 int32_t start_idx, int32_t out_base,
								 const MRGrammarConfig *config) {
	int32_t open_group_base[65] = {0};
	int32_t close_group_base[65] = {0};
	int32_t t_idx = start_idx;

	for (int64_t b = 0; b <= k; b++) {
		open_group_base[b] = t_idx;
		t_idx += (int32_t)MR_grammar_get_group_size(config, b, n, k);
	}
	for (int64_t b = 0; b <= k; b++) {
		close_group_base[b] = t_idx;
		t_idx += (int32_t)MR_grammar_get_group_size(config, b, n, k);
	}

	for (int64_t i = 0; i < n; i++) {
		int64_t b, j;
		MR_grammar_get_bracket_layout(config, i, n, k, &b, &j);

		map[open_group_base[b] + j] = out_base + (int32_t)i;
		map[close_group_base[b] + j] = out_base + (int32_t)n + (int32_t)i;
	}
}

int32_t *build_grammar_to_straight_map(const MRGrammarConfig *config,
									   int64_t terms_count, int64_t n_par,
									   int64_t n_bra, int64_t k) {
	int32_t *map = (int32_t *)malloc((size_t)terms_count * sizeof(int32_t));
	if (!map) {
		return NULL;
	}

	switch (config->kind) {
		case MR_GRAMMAR_ALPHA: {
			for (int64_t i = 0; i < n_par; i++) {
				map[i] = (int32_t)i;
				map[n_par + i] = (int32_t)(n_par + i);
			}
			fill_interleaved_map(map, n_bra, k, 2 * (int32_t)n_par,
								 2 * (int32_t)n_par, config);
			break;
		}
		case MR_GRAMMAR_BETA: {
			fill_interleaved_map(map, n_par, k, 0, 0, config);

			int32_t bra_base = 2 * (int32_t)n_par;
			for (int64_t i = 0; i < n_bra; i++) {
				map[bra_base + i] = bra_base + (int32_t)i;
				map[bra_base + (int32_t)n_bra + i] =
					bra_base + (int32_t)n_bra + (int32_t)i;
			}
			break;
		}
	}
	if (terms_count > 2 * n_par + 2 * n_bra) {
		map[2 * n_par + 2 * n_bra] = (int32_t)(2 * n_par + 2 * n_bra);
	}
	return map;
}
