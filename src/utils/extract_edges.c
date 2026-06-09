#include "extract_edges.h"

#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "internal/grb_utils.h"
#include "uthash.h"

typedef struct {
	GrB_Index row_idx;
	GrB_Index col_idx;
	int32_t nonterm;
} StackElem;

typedef struct {
	StackElem *elems;
	int64_t top;
	int64_t capacity;
} Stack;

static Stack stack_new(void) {
	return (Stack){.elems = NULL, .top = 0, .capacity = 0};
}

static void stack_free(Stack *s) { free(s->elems); }

static bool stack_push(Stack *s, GrB_Index i, GrB_Index j, int32_t A) {
	if (s->top == s->capacity) {
		int64_t new_cap = s->capacity == 0 ? 64 : s->capacity * 2;
		StackElem *tmp = realloc(s->elems, (size_t)new_cap * sizeof(StackElem));
		if (!tmp) {
			return false;
		}
		s->elems = tmp;
		s->capacity = new_cap;
	}
	s->elems[s->top++] = (StackElem){i, j, A};
	return true;
}

static StackElem stack_pop(Stack *s) { return s->elems[--s->top]; }

static bool stack_empty(Stack *s) { return s->top == 0; }

typedef struct {
	GrB_Index i;
	GrB_Index j;
	int32_t A;
	int32_t _pad;
} VisitedKey;

typedef struct {
	VisitedKey key;
	UT_hash_handle hh;
} VisitedEntry;

typedef struct {
	int32_t term; // terminal index (prod_A of the rule)
} TermRule;

typedef struct {
	int32_t B; // left nonterminal (prod_A of the rule)
	int32_t C; // right nonterminal (prod_B of the rule)
} BinaryRule;

typedef struct {
	TermRule *term_rules;
	int32_t term_count;
	BinaryRule *binary_rules;
	int32_t binary_count;
} NontermRules;

static void fill_interleaved_map(int32_t *map, int64_t n_elems, int64_t k,
								 int32_t start_idx, int32_t out_base) {
	int32_t open_group_base[64] = {0};
	int32_t close_group_base[64] = {0};
	int32_t t_idx = start_idx;

	for (int64_t b = 0; b < k; b++) {
		open_group_base[b] = t_idx;
		t_idx += (b < n_elems) ? (int32_t)((n_elems - b + k - 1) / k) : 0;
	}
	for (int64_t b = 0; b < k; b++) {
		close_group_base[b] = t_idx;
		t_idx += (b < n_elems) ? (int32_t)((n_elems - b + k - 1) / k) : 0;
	}
	for (int64_t i = 0; i < n_elems; i++) {
		int64_t b = i % k;
		int64_t j = i / k;
		map[open_group_base[b] + j] = out_base + (int32_t)i;
		map[close_group_base[b] + j] = out_base + (int32_t)n_elems + (int32_t)i;
	}
}

static NontermRules *build_grammar_index(const MRGrammar_t *g) {
	NontermRules *idx = calloc((size_t)g->nonterms_count, sizeof(NontermRules));
	if (!idx) {
		return NULL;
	}

	// Count pass
	for (int64_t r = 0; r < g->rules_count; r++) {
		if (g->rules[r].prod_A == -1) {
			continue;
		}

		int32_t A = g->rules[r].nonterm;
		uint32_t count = g->rules[r].indexed_count;
		uint32_t flags = g->rules[r].indexed;
		int32_t n = (count > 0) ? (int32_t)count : 1;

		for (int32_t k = 0; k < n; k++) {
			int32_t cur_A = A + ((flags & LAGraph_EWNCF_INDEX_NONTERM) ? k : 0);
			if (g->rules[r].prod_B == -1) {
				// Terminal rule
				idx[cur_A].term_count++;
			} else {
				// Binary rule
				idx[cur_A].binary_count++;
			}
		}
	}

	// Alloc pass
	for (int32_t a = 0; a < g->nonterms_count; a++) {
		if (idx[a].term_count) {
			idx[a].term_rules = malloc((size_t)idx[a].term_count * sizeof(TermRule));
			if (!idx[a].term_rules) {
				goto oom;
			}
		}
		if (idx[a].binary_count) {
			idx[a].binary_rules =
				malloc((size_t)idx[a].binary_count * sizeof(BinaryRule));
			if (!idx[a].binary_rules) {
				goto oom;
			}
		}
		// Reset - reuse as write cursor
		idx[a].term_count = 0;
		idx[a].binary_count = 0;
	}

	// Fill pass
	for (int64_t r = 0; r < g->rules_count; r++) {
		if (g->rules[r].prod_A == -1) {
			continue;
		}

		int32_t A = g->rules[r].nonterm;
		int32_t prod_A = g->rules[r].prod_A;
		int32_t prod_B = g->rules[r].prod_B;
		uint32_t count = g->rules[r].indexed_count;
		uint32_t flags = g->rules[r].indexed;
		int32_t n = (count > 0) ? (int32_t)count : 1;

		for (int32_t k = 0; k < n; k++) {
			int32_t cur_A = A + ((flags & LAGraph_EWNCF_INDEX_NONTERM) ? k : 0);
			int32_t cur_pA = prod_A + ((flags & LAGraph_EWNCF_INDEX_PROD_A) ? k : 0);
			int32_t cur_pB = (prod_B != -1 && (flags & LAGraph_EWNCF_INDEX_PROD_B))
								 ? prod_B + k
								 : prod_B;

			if (cur_pB == -1) {
				idx[cur_A].term_rules[idx[cur_A].term_count++].term = cur_pA;
			} else {
				BinaryRule *br = &idx[cur_A].binary_rules[idx[cur_A].binary_count++];
				br->B = cur_pA;
				br->C = cur_pB;
			}
		}
	}

	return idx;

oom:
	for (int32_t a = 0; a < g->nonterms_count; a++) {
		free(idx[a].term_rules);
		free(idx[a].binary_rules);
	}
	free(idx);
	return NULL;
}

static void free_grammar_index(NontermRules *idx, int64_t nonterms_count) {
	if (!idx) {
		return;
	}
	for (int32_t a = 0; a < nonterms_count; a++) {
		free(idx[a].term_rules);
		free(idx[a].binary_rules);
	}
	free(idx);
}

GrB_Info extract_edges_from_outputs(GrB_Matrix *out, GrB_Matrix *paths,
									GrB_Matrix *adj_matrices, MRGrammar_t grammar,
									GrB_Index n, int64_t n_par, int64_t n_bra,
									bool is_beta_parity_group,
									const TargetPath *target_path) {
	GrB_Info info = GrB_SUCCESS;

	NontermRules *grammar_idx = NULL;
	VisitedEntry *visited_ht = NULL; // visited - hash table
	Stack stack = stack_new();
	GrB_Index *rows = NULL;
	GrB_Index *cols = NULL;
	int32_t *grammar_to_straight = NULL;

	// Initialize out (result matrices)
	for (int64_t t = 0; t < grammar.terms_count; t++) {
		GRB_TRY(GrB_Matrix_new(&out[t], GrB_BOOL, n, n));
	}

	grammar_to_straight = malloc((size_t)grammar.terms_count * sizeof(int32_t));
	if (!grammar_to_straight) {
		return GrB_OUT_OF_MEMORY;
	}

	if (!is_beta_parity_group) {
		for (int64_t i = 0; i < n_par; i++) {
			grammar_to_straight[i] = (int32_t)i;
			grammar_to_straight[n_par + i] = (int32_t)(n_par + i);
		}
		fill_interleaved_map(grammar_to_straight, n_bra, grammar.k,
							 2 * (int32_t)n_par, 2 * (int32_t)n_par);
	} else {
		fill_interleaved_map(grammar_to_straight, n_par, grammar.k, 0, 0);
		for (int64_t i = 0; i < n_bra; i++) {
			grammar_to_straight[2 * n_par + i] = (int32_t)(2 * n_par + i);
			grammar_to_straight[2 * n_par + n_bra + i] =
				(int32_t)(2 * n_par + n_bra + i);
		}
	}
	if (grammar.terms_count > 2 * n_par + 2 * n_bra) {
		grammar_to_straight[2 * n_par + 2 * n_bra] =
			(int32_t)(2 * n_par + 2 * n_bra);
	}

	grammar_idx = build_grammar_index(&grammar);
	if (!grammar_idx) {
		info = GrB_OUT_OF_MEMORY;
		goto cleanup;
	}

	// Stack initialization
	if (target_path != NULL) {
		// On-demand analysis

		if (!stack_push(&stack, target_path->src, target_path->tgt, NT_START)) {
			info = GrB_OUT_OF_MEMORY;
			goto cleanup;
		}
	} else {
		// Default: all pairs from paths[NT_START]

		GrB_Index nnz = 0;
		GrB_Matrix_nvals(&nnz, paths[NT_START]);
		if (nnz == 0) {
			goto cleanup;
		}

		rows = malloc(nnz * sizeof(GrB_Index));
		cols = malloc(nnz * sizeof(GrB_Index));
		if (!rows || !cols) {
			info = GrB_OUT_OF_MEMORY;
			goto cleanup;
		}

		GRB_TRY(GrB_Matrix_extractTuples(rows, cols, NULL, &nnz, paths[NT_START]));

		for (GrB_Index k = 0; k < nnz; k++) {
			if (rows[k] != cols[k]) {
				if (!stack_push(&stack, rows[k], cols[k], NT_START)) {
					info = GrB_OUT_OF_MEMORY;
					goto cleanup;
				}
			}
		}

		free((void *)rows);
		free((void *)cols);
		rows = NULL;
		cols = NULL;
	}

	// Main DFS loop
	while (!stack_empty(&stack)) {
		StackElem cur = stack_pop(&stack);
		GrB_Index i = cur.row_idx;
		GrB_Index j = cur.col_idx;
		int32_t A = cur.nonterm;

		// Hash table search
		VisitedKey lk = {.i = i, .j = j, .A = A, ._pad = 0};
		VisitedEntry *found;
		HASH_FIND(hh, visited_ht, &lk, sizeof(VisitedKey), found);
		if (found) {
			continue;
		}

		// Add to the hash table
		VisitedEntry *new_entry = malloc(sizeof(VisitedEntry));
		if (!new_entry) {
			info = GrB_OUT_OF_MEMORY;
			goto cleanup;
		}
		new_entry->key = lk;
		HASH_ADD(hh, visited_ht, key, sizeof(VisitedKey), new_entry);

		AllPathsElem elem;
		if (GrB_Matrix_extractElement_UDT(&elem, paths[A], i, j) == GrB_NO_VALUE) {
			continue;
		}

		const NontermRules *rules = &grammar_idx[A];

		// Iterate over intermediate vertices
		for (GrB_Index ki = 0; ki < elem.n; ki++) {
			GrB_Index mid =
				(elem.n == 1) ? elem.data.single_elem : elem.data.middle[ki];

			if (mid == GrB_INDEX_MAX) {
				// Terminal edge: A -> term
				for (int32_t ri = 0; ri < rules->term_count; ri++) {
					int32_t term = rules->term_rules[ri].term;

					bool has_edge = false;
					GrB_Matrix_extractElement_BOOL(&has_edge, adj_matrices[term], i,
												   j);
					if (has_edge) {
						int32_t straight_idx =
							grammar_to_straight[term -
												(int32_t)grammar.nonterms_count];
						GRB_TRY(GrB_Matrix_setElement_BOOL(out[straight_idx], true,
														   i, j));
					}
				}
			} else {
				// Binary split: A -> B C via mid
				for (int32_t ri = 0; ri < rules->binary_count; ri++) {
					int32_t B = rules->binary_rules[ri].B;
					int32_t C = rules->binary_rules[ri].C;

					// Check that both subpaths exist in paths
					AllPathsElem dummy;
					if (GrB_Matrix_extractElement_UDT(&dummy, paths[B], i, mid) ==
						GrB_NO_VALUE) {
						continue;
					}
					if (GrB_Matrix_extractElement_UDT(&dummy, paths[C], mid, j) ==
						GrB_NO_VALUE) {
						continue;
					}

					// Only push if not already visited - avoids stacking the
					// same pair multiple times
					VisitedKey lkB = {i, mid, B, ._pad = 0};
					VisitedKey lkC = {mid, j, C, ._pad = 0};
					VisitedEntry *fB, *fC;
					HASH_FIND(hh, visited_ht, &lkB, sizeof(VisitedKey), fB);
					HASH_FIND(hh, visited_ht, &lkC, sizeof(VisitedKey), fC);

					if (!fB) {
						if (!stack_push(&stack, i, mid, B)) {
							info = GrB_OUT_OF_MEMORY;
							goto cleanup;
						}
					}
					if (!fC) {
						if (!stack_push(&stack, mid, j, C)) {
							info = GrB_OUT_OF_MEMORY;
							goto cleanup;
						}
					}
				}
			}
		}
	}

cleanup:
	free(grammar_to_straight);
	free(rows);
	free(cols);
	VisitedEntry *curr_entry, *tmp_entry;
	HASH_ITER(hh, visited_ht, curr_entry, tmp_entry) {
		HASH_DEL(visited_ht, curr_entry);
		free(curr_entry);
	}
	free_grammar_index(grammar_idx, grammar.nonterms_count);
	stack_free(&stack);
	return info;
}
