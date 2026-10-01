#include "extract_edges.h"

#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "internal/grb_utils.h"

typedef struct {
	GrB_Index col;
	AllPathsElem val;
} FastEntry;

typedef struct {
	uint32_t count;
	GrB_Index row_offset;
	FastEntry *entries;
} FastRow;

typedef struct {
	FastRow *rows;
	GrB_Index n;
	GrB_Index nnz;
	uint64_t *visited;
} FastPathsMatrix;

// The threshold for the number of elements in a row to switch from linear to binary
// search.
static const uint32_t FAST_PATHS_LINEAR_THRESHOLD = 16;

/**
 * @brief Tests if a specific bit is set in a bitset array.
 *
 * @param[in] b Pointer to the bitset array.
 * @param[in] p Bit index to test.
 *
 * @return `true` if the bit is set, `false` otherwise or if `b` is `NULL`.
 */
static inline bool test_bit(const uint64_t *b, GrB_Index p) {
	if (!b) {
		return false;
	}
	return (b[p >> 6] & (1ULL << (p & 63))) != 0;
}

/**
 * @brief Tests and sets a specific bit in a bitset array.
 *
 * @param[in,out] b Pointer to the bitset array.
 * @param[in]     p Bit index to set.
 *
 * @return Previous state of the bit before setting (`true` if set, `false`
 * otherwise).
 */
static inline bool test_and_set(uint64_t *b, GrB_Index p) {
	uint64_t m = 1ULL << (p & 63);
	bool old = (b[p >> 6] & m) != 0;
	b[p >> 6] |= m;
	return old;
}

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
	free(fmat->visited);
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
	fmat->visited = NULL;
	fmat->rows = calloc(n, sizeof(FastRow));
	if (!fmat->rows && n > 0) {
		free(fmat);
		return NULL;
	}

	GrB_Index nnz = 0;
	GrB_Matrix_nvals(&nnz, mat);
	fmat->nnz = nnz;
	if (nnz == 0) {
		return fmat;
	}

	fmat->visited = calloc((size_t)((nnz + 63) / 64), sizeof(uint64_t));
	if (!fmat->visited) {
		free(fmat->rows);
		free(fmat);
		return NULL;
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

	// Calculating row_offset and allocating memory for records
	GrB_Index current_offset = 0;
	for (GrB_Index r = 0; r < n; r++) {
		fmat->rows[r].row_offset = current_offset;
		if (fmat->rows[r].count > 0) {
			fmat->rows[r].entries = malloc(fmat->rows[r].count * sizeof(FastEntry));
			if (!fmat->rows[r].entries) {
				goto fail;
			}
			current_offset += fmat->rows[r].count;
			fmat->rows[r].count = 0;
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
 * @brief Looks up value at coordinate (i, j) in a `FastPathsMatrix` and returns its
 * linear index.
 *
 * @param[in]  fmat     Input lookup matrix.
 * @param[in]  i        Row index.
 * @param[in]  j        Column index.
 * @param[out] out_idx  Pointer to index of retrieved `AllPathsElem`.
 * @param[out] out_val  Pointer to value of retrieved `AllPathsElem`.
 *
 * @return `true` if element exists, `false` otherwise.
 */
static bool fast_paths_get_with_index(const FastPathsMatrix *fmat, GrB_Index i,
									  GrB_Index j, GrB_Index *out_idx,
									  AllPathsElem *out_val) {
	if (!fmat || i >= fmat->n) {
		return false;
	}
	const FastRow *row = &fmat->rows[i];
	uint32_t count = row->count;

	if (count == 0) {
		return false;
	}

	int found_k = -1;

	if (count <= FAST_PATHS_LINEAR_THRESHOLD) {
		for (uint32_t k = 0; k < count; k++) {
			if (row->entries[k].col == j) {
				found_k = (int)k;
				break;
			}
		}
	} else {
		int low = 0, high = (int)count - 1;
		while (low <= high) {
			int mid = low + (high - low) / 2;
			if (row->entries[mid].col == j) {
				found_k = mid;
				break;
			}
			if (row->entries[mid].col < j) {
				low = mid + 1;
			} else {
				high = mid - 1;
			}
		}
	}

	if (found_k < 0) {
		return false;
	}

	if (out_idx) {
		*out_idx = row->row_offset + (GrB_Index)found_k;
	}
	if (out_val) {
		*out_val = row->entries[found_k].val;
	}
	return true;
}

/**
 * @brief Checks if element at coordinate (i, j) has already been visited.
 *
 * @param[in] fmat Input lookup matrix.
 * @param[in] i    Row index.
 * @param[in] j    Column index.
 *
 * @return `true` if element exists and is marked visited, `false` otherwise.
 */
static inline bool is_visited(const FastPathsMatrix *fmat, GrB_Index i,
							  GrB_Index j) {
	if (!fmat) {
		return false;
	}
	GrB_Index idx;
	if (fast_paths_get_with_index(fmat, i, j, &idx, NULL)) {
		return test_bit(fmat->visited, idx);
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
	int32_t term; // terminal index (prod_A of the rule)
} TermRule;

typedef struct {
	TermRule *term_rules;
	int32_t term_count;
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
		LAGraph_rule_EWCNF rule = g->rules[r];

		if (rule.xor_routing == NULL && rule.prod_A != -1) {
			int32_t A = rule.nonterm;
			uint32_t count = rule.indexed_count;
			uint32_t flags = rule.indexed;
			int32_t n = (count > 0) ? (int32_t)count : 1;

			if (rule.prod_B == -1) {
				for (int32_t k = 0; k < n; k++) {
					int32_t cur_A =
						A + ((flags & LAGraph_EWNCF_INDEX_NONTERM) ? k : 0);
					idx[cur_A].term_count++;
				}
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
		// Reset - reuse as write cursor
		idx[a].term_count = 0;
	}

	// Fill pass
	for (int64_t r = 0; r < g->rules_count; r++) {
		LAGraph_rule_EWCNF rule = g->rules[r];

		if (rule.xor_routing == NULL && rule.prod_A != -1 && rule.prod_B == -1) {
			int32_t A = rule.nonterm;
			int32_t prod_A = rule.prod_A;
			uint32_t count = rule.indexed_count;
			uint32_t flags = rule.indexed;
			int32_t n = (count > 0) ? (int32_t)count : 1;

			for (int32_t k = 0; k < n; k++) {
				int32_t cur_A = A + ((flags & LAGraph_EWNCF_INDEX_NONTERM) ? k : 0);
				int32_t cur_pA =
					prod_A + ((flags & LAGraph_EWNCF_INDEX_PROD_A) ? k : 0);

				idx[cur_A].term_rules[idx[cur_A].term_count++].term = cur_pA;
			}
		}
	}

	return idx;

oom:
	for (int32_t a = 0; a < g->nonterms_count; a++) {
		free(idx[a].term_rules);
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

static inline bool mid_entry_uses_heap(const MidEntry *e) {
	return e->rule_count > MID_ENTRY_INLINE_CAP + 1;
}

static inline int32_t mid_entry_rest_id(const MidEntry *e, uint32_t k) {
	if (mid_entry_uses_heap(e)) {
		return e->rest.rule_ids_rest[k];
	}
	return e->rest.inline_ids[k];
}

GrB_Info extract_edges_from_outputs(GrB_Matrix *out, GrB_Matrix *paths,
									const BinaryRuleInfo *rule_table,
									GrB_Matrix *adj_matrices, MRGrammar_t grammar,
									GrB_Index n, int64_t n_par, int64_t n_bra,
									const MRGrammarConfig *config,
									const TargetPath *target_path) {
	GrB_Info info = GrB_SUCCESS;

	NontermRules *grammar_idx = NULL;
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

	fast_paths = (FastPathsMatrix **)calloc(grammar.nonterms_count,
											sizeof(FastPathsMatrix *));
	fast_adj =
		(FastPathsMatrix **)calloc(grammar.terms_count, sizeof(FastBoolMatrix *));
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

		FastPathsMatrix *cur_paths = NULL;
		GRB_TRY(get_fast_paths(&cur_paths, fast_paths, paths, A, n));

		GrB_Index elem_idx;
		AllPathsElem elem;
		if (!fast_paths_get_with_index(cur_paths, i, j, &elem_idx, &elem)) {
			continue;
		}

		// Check and set the 'visited' bit
		if (test_and_set(cur_paths->visited, elem_idx)) {
			continue;
		}

		const NontermRules *rules = &grammar_idx[A];

		// Iterate over intermediate vertices
		for (GrB_Index ki = 0; ki < elem.n; ki++) {
			MidEntry entry =
				(elem.n == 1) ? elem.data.single_elem : elem.data.middle[ki];
			GrB_Index mid = entry.mid;

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
				continue;
			}

			// Binary split: A -> B C via mid

			uint32_t rule_count = entry.rule_count;

			for (uint32_t part = 0; part < rule_count; part++) {
				int32_t rid = (part == 0) ? entry.rule_id0
										  : mid_entry_rest_id(&entry, part - 1);
				BinaryRuleInfo ri = rule_table[rid];

				// Only push if not already visited - avoids stacking the
				// same pair multiple times
				if (!is_visited(fast_paths[ri.B], i, mid)) {
					if (!stack_push(&stack, i, mid, ri.B)) {
						info = GrB_OUT_OF_MEMORY;
						goto cleanup;
					}
				}
				if (!is_visited(fast_paths[ri.C], mid, j)) {
					if (!stack_push(&stack, mid, j, ri.C)) {
						info = GrB_OUT_OF_MEMORY;
						goto cleanup;
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
		free((void *)fast_paths);
	}
	if (fast_adj) {
		for (int64_t t = 0; t < grammar.terms_count; t++) {
			fast_bool_free(fast_adj[t]);
		}
		free((void *)fast_adj);
	}
	free(grammar_to_straight);
	free(rows);
	free(cols);
	free_grammar_index(grammar_idx, grammar.nonterms_count);
	stack_free(&stack);
	return info;
}
