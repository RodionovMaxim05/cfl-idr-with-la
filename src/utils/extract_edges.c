#include "extract_edges.h"

#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "internal/grb_utils.h"
#include "uthash.h"

typedef struct {
	GrB_Index col;
	AllPathsElem val;
} FastEntry;

typedef struct {
	uint32_t count;
	FastEntry *entries;
} FastRow;

typedef struct {
	FastRow *rows;
	GrB_Index n;
} FastPathsMatrix;

// The threshold for the number of elements in a row to switch from linear to binary
// search.
static const uint32_t FAST_PATHS_LINEAR_THRESHOLD = 16;

/**
 * @brief Comparator for `qsort()` to sort `FastEntry` items by column index.
 */
static int compare_fast_entries(const void *a, const void *b) {
	GrB_Index ca = ((const FastEntry *)a)->col;
	GrB_Index cb = ((const FastEntry *)b)->col;
	return (ca > cb) - (ca < cb);
}

static void fast_paths_free(FastPathsMatrix *fmat) {
	if (!fmat) {
		return;
	}
	if (fmat->rows) {
		for (GrB_Index r = 0; r < fmat->n; r++) {
			free(fmat->rows[r].entries);
		}
		free(fmat->rows);
	}
	free(fmat);
}

/**
 * @brief Builds a fast row-oriented lookup structure for a sparse path matrix.
 *
 * Extracts tuples from `mat`, groups entries by row, and sorts each row by column
 * index to allow fast binary search lookups.
 *
 * @param[in] mat Source GraphBLAS matrix containing path elements.
 * @param[in] n   Matrix dimension (number of rows/cols).
 *
 * @return Pointer to newly allocated `FastPathsMatrix`, or `NULL` on allocation
 * failure.
 */
static FastPathsMatrix *fast_paths_build(GrB_Matrix mat, GrB_Index n) {
	FastPathsMatrix *fmat = malloc(sizeof(FastPathsMatrix));
	if (!fmat) {
		return NULL;
	}
	fmat->n = n;
	fmat->rows = calloc(n, sizeof(FastRow));
	if (!fmat->rows && n > 0) {
		free(fmat);
		return NULL;
	}

	GrB_Index nnz = 0;
	GrB_Matrix_nvals(&nnz, mat);
	if (nnz == 0) {
		return fmat;
	}

	GrB_Index *row_indices = malloc(nnz * sizeof(GrB_Index));
	GrB_Index *col_indices = malloc(nnz * sizeof(GrB_Index));
	AllPathsElem *values = malloc(nnz * sizeof(AllPathsElem));
	if (!row_indices || !col_indices || !values) {
		goto fail;
	}

	GrB_Matrix_extractTuples_UDT(row_indices, col_indices, values, &nnz, mat);

	// Count the elements in each row
	for (GrB_Index k = 0; k < nnz; k++) {
		fmat->rows[row_indices[k]].count++;
	}

	// Allocate storage for each row's entries
	for (GrB_Index r = 0; r < n; r++) {
		if (fmat->rows[r].count > 0) {
			fmat->rows[r].entries = malloc(fmat->rows[r].count * sizeof(FastEntry));
			if (!fmat->rows[r].entries) {
				goto fail;
			}
			fmat->rows[r].count = 0; // Reset to use as a write cursor below
		}
	}

	// Fill each row's entries from the extracted tuples
	for (GrB_Index k = 0; k < nnz; k++) {
		GrB_Index r = row_indices[k];
		uint32_t idx = fmat->rows[r].count++;
		fmat->rows[r].entries[idx] =
			(FastEntry){.col = col_indices[k], .val = values[k]};
	}

	// Sorting columns for binary search
	for (GrB_Index r = 0; r < n; r++) {
		if (fmat->rows[r].count > 1) {
			qsort(fmat->rows[r].entries, fmat->rows[r].count, sizeof(FastEntry),
				  compare_fast_entries);
		}
	}

	free(row_indices);
	free(col_indices);
	free(values);
	return fmat;

fail:
	free(row_indices);
	free(col_indices);
	free(values);
	fast_paths_free(fmat);
	return NULL;
}

/**
 * @brief Looks up value at coordinate (i, j) in a `FastPathsMatrix`.
 *
 * @param[in]  fmat Input lookup matrix.
 * @param[in]  i    Row index.
 * @param[in]  j    Column index.
 * @param[out] out  Pointer to store the retrieved `AllPathsElem`.
 *
 * @return `true` if element exists and was written to `out`, `false` otherwise.
 */
static bool fast_paths_get(const FastPathsMatrix *fmat, GrB_Index i, GrB_Index j,
						   AllPathsElem *out) {
	if (!fmat || i >= fmat->n) {
		return false;
	}
	const FastRow *row = &fmat->rows[i];
	uint32_t count = row->count;

	if (count == 0) {
		return false;
	}

	if (count <= FAST_PATHS_LINEAR_THRESHOLD) {
		for (uint32_t k = 0; k < count; k++) {
			if (row->entries[k].col == j) {
				*out = row->entries[k].val;
				return true;
			}
		}
		return false;
	}

	// Binary search
	int low = 0, high = (int)count - 1;
	while (low <= high) {
		int mid = low + (high - low) / 2;
		if (row->entries[mid].col == j) {
			*out = row->entries[mid].val;
			return true;
		}
		if (row->entries[mid].col < j) {
			low = mid + 1;
		} else {
			high = mid - 1;
		}
	}
	return false;
}

typedef struct {
	uint32_t count;
	GrB_Index *cols;
} FastBoolRow;

typedef struct {
	FastBoolRow *rows;
	GrB_Index n;
} FastBoolMatrix;

static void fast_bool_free(FastBoolMatrix *fmat) {
	if (!fmat) {
		return;
	}
	if (fmat->rows) {
		for (GrB_Index r = 0; r < fmat->n; r++) {
			free(fmat->rows[r].cols);
		}
		free(fmat->rows);
	}
	free(fmat);
}

/**
 * @brief Builds a fast row-oriented lookup structure for a boolean adjacency matrix.
 *
 * @param[in] mat Source GraphBLAS boolean matrix.
 * @param[in] n   Matrix dimension.
 *
 * @return Pointer to `FastBoolMatrix`, or `NULL` on allocation failure.
 */
static FastBoolMatrix *fast_bool_build(GrB_Matrix mat, GrB_Index n) {
	FastBoolMatrix *fmat = malloc(sizeof(FastBoolMatrix));
	if (!fmat) {
		return NULL;
	}
	fmat->n = n;
	fmat->rows = calloc(n, sizeof(FastBoolRow));
	if (!fmat->rows && n > 0) {
		free(fmat);
		return NULL;
	}

	GrB_Index nnz = 0;
	GrB_Matrix_nvals(&nnz, mat);
	if (nnz == 0) {
		return fmat;
	}

	GrB_Index *row_indices = malloc(nnz * sizeof(GrB_Index));
	GrB_Index *col_indices = malloc(nnz * sizeof(GrB_Index));
	if (!row_indices || !col_indices) {
		goto fail;
	}

	GrB_Matrix_extractTuples_BOOL(row_indices, col_indices, NULL, &nnz, mat);

	for (GrB_Index k = 0; k < nnz; k++) {
		fmat->rows[row_indices[k]].count++;
	}
	for (GrB_Index r = 0; r < n; r++) {
		if (fmat->rows[r].count > 0) {
			fmat->rows[r].cols = malloc(fmat->rows[r].count * sizeof(GrB_Index));
			if (!fmat->rows[r].cols) {
				goto fail;
			}
			fmat->rows[r].count = 0;
		}
	}
	for (GrB_Index k = 0; k < nnz; k++) {
		GrB_Index r = row_indices[k];
		fmat->rows[r].cols[fmat->rows[r].count++] = col_indices[k];
	}
	free(row_indices);
	free(col_indices);
	return fmat;

fail:
	free(row_indices);
	free(col_indices);
	fast_bool_free(fmat);
	return NULL;
}

/**
 * @brief Checks if edge (i, j) exists in a `FastBoolMatrix`.
 *
 * @param[in] fmat Input boolean lookup matrix.
 * @param[in] i    Row index.
 * @param[in] j    Column index.
 *
 * @return `true` if edge exists, `false` otherwise.
 */
static bool fast_bool_has(const FastBoolMatrix *fmat, GrB_Index i, GrB_Index j) {
	if (!fmat || i >= fmat->n) {
		return false;
	}
	const FastBoolRow *row = &fmat->rows[i];
	uint32_t count = row->count;
	if (count == 0) {
		return false;
	}
	for (uint32_t k = 0; k < count; k++) {
		if (row->cols[k] == j) {
			return true;
		}
	}
	return false;
}

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

/**
 * @brief Builds indexed grammar rules lookup array for nonterminals.
 *
 * @param[in] g Input grammar pointer.
 *
 * @return Pointer to array of `NontermRules`, or `NULL` on memory allocation
 * failure.
 */
static NontermRules *build_grammar_index(const MRGrammar_t *g) {
	NontermRules *idx = calloc((size_t)g->nonterms_count, sizeof(NontermRules));
	if (!idx) {
		return NULL;
	}

	// Count pass
	for (int64_t r = 0; r < g->rules_count; r++) {
		if (g->rules[r].prod_A == -1) {
			// epsilon - ignore
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

/**
 * @brief Lazily retrieves or builds the `FastPathsMatrix` index for `paths[a]`.
 *
 * @param[out]    out_fmat   Pointer to receive the target `FastPathsMatrix`.
 * @param[in,out] fast_paths Array of cached `FastPathsMatrix` structures.
 * @param[in]     paths      Array of raw path matrices.
 * @param[in]     a          Nonterminal symbol index.
 * @param[in]     n          Matrix dimension.
 *
 * @return `GrB_SUCCESS` on success, or `GrB_OUT_OF_MEMORY` on allocation failure.
 */
static inline GrB_Info get_fast_paths(FastPathsMatrix **out_fmat,
									  FastPathsMatrix **fast_paths,
									  GrB_Matrix *paths, int32_t a, GrB_Index n) {
	if (fast_paths[a] == NULL) {
		fast_paths[a] = fast_paths_build(paths[a], n);
		if (!fast_paths[a]) {
			*out_fmat = NULL;
			return GrB_OUT_OF_MEMORY;
		}
	}
	*out_fmat = fast_paths[a];
	return GrB_SUCCESS;
}

/**
 * @brief Lazily retrieves or builds the `FastBoolMatrix` index for terminal
 * adjacency.
 *
 * @param[out]    out_fadj       Pointer to receive the target `FastBoolMatrix`.
 * @param[in,out] fast_adj       Array of cached `FastBoolMatrix` structures.
 * @param[in]     adj_matrices   Array of terminal adjacency matrices.
 * @param[in]     nonterms_count Number of nonterminal symbols.
 * @param[in]     term_idx       Terminal symbol index.
 * @param[in]     n              Matrix dimension.
 *
 * @return `GrB_SUCCESS` on success, or `GrB_OUT_OF_MEMORY` on allocation failure.
 */
static inline GrB_Info get_fast_adj(FastBoolMatrix **out_fadj,
									FastBoolMatrix **fast_adj,
									GrB_Matrix *adj_matrices, int32_t nonterms_count,
									int32_t term_idx, GrB_Index n) {
	if (fast_adj[term_idx] == NULL) {
		fast_adj[term_idx] =
			fast_bool_build(adj_matrices[nonterms_count + term_idx], n);
		if (!fast_adj[term_idx]) {
			*out_fadj = NULL;
			return GrB_OUT_OF_MEMORY;
		}
	}
	*out_fadj = fast_adj[term_idx];
	return GrB_SUCCESS;
}

GrB_Info extract_edges_from_outputs(GrB_Matrix *out, GrB_Matrix *paths,
									GrB_Matrix *adj_matrices, MRGrammar_t grammar,
									GrB_Index n, int64_t n_par, int64_t n_bra,
									const MRGrammarConfig *config,
									const TargetPath *target_path) {
	GrB_Info info = GrB_SUCCESS;

	NontermRules *grammar_idx = NULL;
	VisitedEntry *visited_ht = NULL;
	Stack stack = stack_new();
	GrB_Index *rows = NULL;
	GrB_Index *cols = NULL;
	int32_t *grammar_to_straight = NULL;
	FastPathsMatrix **fast_paths = NULL;
	FastBoolMatrix **fast_adj = NULL;

	// Initialize out (result matrices)
	for (int64_t t = 0; t < grammar.terms_count; t++) {
		GRB_TRY(GrB_Matrix_new(&out[t], GrB_BOOL, n, n));
	}

	grammar_to_straight = malloc((size_t)grammar.terms_count * sizeof(int32_t));
	if (!grammar_to_straight) {
		return GrB_OUT_OF_MEMORY;
	}

	switch (config->kind) {
		case MR_GRAMMAR_ALPHA: {
			for (int64_t i = 0; i < n_par; i++) {
				grammar_to_straight[i] = (int32_t)i;
				grammar_to_straight[n_par + i] = (int32_t)(n_par + i);
			}

			fill_interleaved_map(grammar_to_straight, n_bra, grammar.k,
								 2 * (int32_t)n_par, 2 * (int32_t)n_par, config);
			break;
		}
		case MR_GRAMMAR_BETA: {
			fill_interleaved_map(grammar_to_straight, n_par, grammar.k, 0, 0,
								 config);

			int32_t bra_base = 2 * (int32_t)n_par;
			for (int64_t i = 0; i < n_bra; i++) {
				grammar_to_straight[bra_base + i] = bra_base + (int32_t)i;
				grammar_to_straight[bra_base + (int32_t)n_bra + i] =
					bra_base + (int32_t)n_bra + (int32_t)i;
			}
			break;
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

	fast_paths = calloc(grammar.nonterms_count, sizeof(FastPathsMatrix *));
	fast_adj = calloc(grammar.terms_count, sizeof(FastBoolMatrix *));
	if (!fast_paths || !fast_adj) {
		info = GrB_OUT_OF_MEMORY;
		goto cleanup;
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
		FastPathsMatrix *cur_paths = NULL;
		GRB_TRY(get_fast_paths(&cur_paths, fast_paths, paths, A, n));
		if (!fast_paths_get(cur_paths, i, j, &elem)) {
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
					int32_t term_idx = term - (int32_t)grammar.nonterms_count;

					FastBoolMatrix *adj = NULL;
					GRB_TRY(get_fast_adj(&adj, fast_adj, adj_matrices,
										 grammar.nonterms_count, term_idx, n));
					if (fast_bool_has(adj, i, j)) {
						int32_t straight_idx = grammar_to_straight[term_idx];
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
					FastPathsMatrix *B_paths = NULL;
					GRB_TRY(get_fast_paths(&B_paths, fast_paths, paths, B, n));
					if (!fast_paths_get(B_paths, i, mid, &dummy)) {
						continue;
					}
					FastPathsMatrix *C_paths = NULL;
					GRB_TRY(get_fast_paths(&C_paths, fast_paths, paths, C, n));
					if (!fast_paths_get(C_paths, mid, j, &dummy)) {
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
	if (fast_paths) {
		for (int32_t a = 0; a < grammar.nonterms_count; a++) {
			fast_paths_free(fast_paths[a]);
		}
		free(fast_paths);
	}
	if (fast_adj) {
		for (int64_t t = 0; t < grammar.terms_count; t++) {
			fast_bool_free(fast_adj[t]);
		}
		free(fast_adj);
	}
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
