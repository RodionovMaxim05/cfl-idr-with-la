#include "extract_edges.h"

#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "internal/grb_utils.h"

/**
 * @brief Memory budget divisor for bucket index allocation.
 *
 * Caps bucket index memory footprint to 1 / BUCKET_MEM_BUDGET_DIVISOR (25%)
 * of the underlying HyperCSR matrix structure size.
 */
#define BUCKET_MEM_BUDGET_DIVISOR 4

/**
 * @brief Maximum range size threshold for switching from binary search to linear
 * scan.
 *
 * Ranges with a length less than or equal to BINARY_SEARCH_THRESHOLD are scanned
 * linearly to minimize branch mispredictions and binary search overheads.
 */
#define BINARY_SEARCH_THRESHOLD 8

/**
 * @brief Lightweight zero-copy HyperCSR view of a nonterminal path matrix.
 *
 * Wraps unpacked GraphBLAS HyperCSR array pointers for hyper-sparse element lookup
 * with optional bucket indexing, and tracks visited matrix entries using a bitset to
 * prevent redundant traversals during DFS.
 */
typedef struct {
	GrB_Matrix mat;	   // Handle to the underlying GraphBLAS matrix
	GrB_Index *Ap;	   // HyperCSR column start pointers array
	GrB_Index *Ah;	   // HyperCSR active row indices array
	GrB_Index *Aj;	   // HyperCSR column indices array
	AllPathsElem *Ax;  // Array of path reconstruction elements
	GrB_Index Ap_size; // Allocated size of column pointers array
	GrB_Index Ah_size; // Allocated size of active row indices array
	GrB_Index Aj_size; // Allocated size of column indices array
	GrB_Index Ax_size; // Allocated size of element values array
	GrB_Index nvec;	   // Number of non-empty rows in the hyper-sparse matrix
	GrB_Index nnz;	   // Total number of non-zero entries in the matrix
	bool iso;		   // True if all non-zero entries share a single iso-value
	bool open;		   // True if the view currently holds unpacked pointers
	uint64_t *visited; // Bitset tracking visited non-zero entries (1 bit per entry)
	uint32_t
		*bstart;	// Array of bucket start indices in Ah for accelerated row lookup
	uint32_t shift; // Bitwise shift amount for bucket index calculation
	GrB_Index nbuckets; // Total number of allocated lookup buckets
} PathsView;

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
 * @brief Packs the HyperCSR matrix back to GraphBLAS ownership and frees internal
 * view resources.
 *
 * @param[in,out] v Pointer to the `PathsView` to close.
 *
 * @return `GrB_SUCCESS` on success, or GraphBLAS error code on packing failure.
 */
static GrB_Info paths_view_close(PathsView *v) {
	if (!v || !v->open) {
		return GrB_SUCCESS;
	}
	GrB_Info info = GxB_Matrix_pack_HyperCSR(
		v->mat, &v->Ap, &v->Ah, &v->Aj, (void **)&v->Ax, v->Ap_size, v->Ah_size,
		v->Aj_size, v->Ax_size, v->iso, v->nvec, false, NULL);
	free(v->visited);
	free(v->bstart);
	memset(v, 0, sizeof(*v));
	return info;
}

/**
 * @brief Builds a bucket index for fast estimation of row boundaries in the HyperCSR
 *
 * @param[in,out] v Pointer to the open `PathsView`.
 * @param[in]     n Total number of vertices in the graph.
 *
 * @return `GrB_SUCCESS` on success, or `GrB_OUT_OF_MEMORY` on allocation failure.
 */
static GrB_Info paths_view_build_buckets(PathsView *v, GrB_Index n) {
	if (v->nvec < 32 || v->nvec > UINT32_MAX) {
		return GrB_SUCCESS;
	}

	// Memory budget: limit size relative to HyperCSR base structure (Ap + Ah + Aj)
	size_t base = ((size_t)v->nvec * 2 + v->nnz) * sizeof(GrB_Index);
	GrB_Index max_nb =
		(GrB_Index)(base / BUCKET_MEM_BUDGET_DIVISOR / sizeof(uint32_t));
	if (max_nb < v->nvec) {
		max_nb = v->nvec; // Fallback to baseline nvec density
	}

	uint32_t shift = 0;
	while (((n >> shift) + 1) > max_nb) {
		shift++;
	}

	GrB_Index nb = (n >> shift) + 1;
	v->bstart = malloc((nb + 1) * sizeof(uint32_t));
	if (!v->bstart) {
		return GrB_OUT_OF_MEMORY;
	}

	v->shift = shift;
	v->nbuckets = nb;

	GrB_Index m = 0;
	for (GrB_Index b = 0; b <= nb; b++) {
		while (m < v->nvec && (v->Ah[m] >> shift) < b) {
			m++;
		}
		v->bstart[b] = (uint32_t)m;
	}

	return GrB_SUCCESS;
}

/**
 * @brief Unpacks a GraphBLAS HyperCSR structure into a zero-copy view without
 * converting to standard CSR.
 *
 * @param[out] v   Pointer to the target `PathsView` structure.
 * @param[in]  mat GraphBLAS matrix to unpack.
 * @param[in]  n   Total number of vertices in the graph.
 *
 * @return `GrB_SUCCESS` on success, or an error code on failure.
 */
static GrB_Info paths_view_open(PathsView *v, GrB_Matrix mat, GrB_Index n) {
	GrB_Info info = GrB_SUCCESS;
	memset(v, 0, sizeof(*v));
	v->mat = mat;

	GRB_TRY(GrB_Matrix_wait(mat, GrB_MATERIALIZE));

	GRB_TRY(GxB_Matrix_unpack_HyperCSR(mat, &v->Ap, &v->Ah, &v->Aj, (void **)&v->Ax,
									   &v->Ap_size, &v->Ah_size, &v->Aj_size,
									   &v->Ax_size, &v->iso, &v->nvec, NULL, NULL));

	v->open = true;
	v->nnz = (v->Ap && v->nvec) ? v->Ap[v->nvec] : 0;

	if (v->nnz > 0) {
		v->visited = calloc((size_t)((v->nnz + 63) / 64), sizeof(uint64_t));
		if (!v->visited) {
			info = GrB_OUT_OF_MEMORY;
			goto cleanup;
		}
	}

	GRB_TRY(paths_view_build_buckets(v, n));

	return GrB_SUCCESS;

cleanup:
	paths_view_close(v);
	return info;
}

/**
 * @brief Looks up a column index `col_idx` in row `row_idx` using bucket indexing,
 * direct row resolution when shift==0, combined binary search, and linear scan.
 *
 * Performs bucket-assisted search on `Ah` to locate the active row, followed by a
 * search on `Aj` to locate the column.
 *
 * @param[in]  v       Pointer to the open `PathsView`.
 * @param[in]  row_idx Index of the row (source vertex).
 * @param[in]  col_idx Index of the column (target vertex).
 * @param[out] out_idx Optional pointer to store the linear element index if found.
 *
 * @return `true` if the entry exists in the matrix view, `false` otherwise.
 */
static inline bool paths_view_find_elem(const PathsView *v, GrB_Index row_idx,
										GrB_Index col_idx, GrB_Index *out_idx) {
	if (!v->open || v->nvec == 0) {
		return false;
	}

	// Step 1: Bucket-assisted search for row_idx in Ah array
	GrB_Index lo, hi;
	if (v->bstart) {
		GrB_Index b = row_idx >> v->shift;
		if (b >= v->nbuckets) {
			return false;
		}
		lo = v->bstart[b];
		hi = v->bstart[b + 1];
		if (v->shift == 0) {
			if (lo == hi) {
				return false; // Row row_idx is not present
			}
			goto have_row; // Ah[lo] == row_idx guaranteed; skip Ah lookup
		}
	} else {
		lo = 0;
		hi = v->nvec;
	}

	// Binary search while range > BINARY_SEARCH_THRESHOLD
	while (hi - lo > BINARY_SEARCH_THRESHOLD) {
		GrB_Index m = lo + (hi - lo) / 2;
		if (v->Ah[m] < row_idx) {
			lo = m + 1;
		} else {
			hi = m;
		}
	}
	// Linear scan for remaining range <= BINARY_SEARCH_THRESHOLD
	while (lo < hi && v->Ah[lo] < row_idx) {
		lo++;
	}

	if (lo >= v->nvec || v->Ah[lo] != row_idx) {
		return false;
	}

have_row:
	// Step 2: Search for col_idx in Aj array
	GrB_Index l = v->Ap[lo], h = v->Ap[lo + 1], e = h;
	// Binary search while range > BINARY_SEARCH_THRESHOLD
	while (h - l > 8) {
		GrB_Index m = l + (h - l) / 2;
		if (v->Aj[m] < col_idx) {
			l = m + 1;
		} else {
			h = m;
		}
	}
	// Linear scan for remaining range <= BINARY_SEARCH_THRESHOLD
	while (l < h && v->Aj[l] < col_idx) {
		l++;
	}

	if (l < e && v->Aj[l] == col_idx) {
		if (out_idx) {
			*out_idx = l;
		}
		return true;
	}
	return false;
}

/**
 * @brief Element stored in the explicit DFS extraction stack.
 *
 * Captures a single step in the CFG path reconstruction traversal, representing an
 * entry $(i, j)$ in the matrix of a nonterminal symbol $A$.
 */
typedef struct {
	GrB_Index row_idx;	// Source vertex index (row i)
	GrB_Index col_idx;	// Target vertex index (column j)
	int32_t nonterm;	// Nonterminal symbol ID
	GrB_Index elem_idx; // Linear index of the element in CSR/HyperCSR arrays
} StackElem;

/**
 * @brief Dynamic stack container for explicit DFS path extraction.
 */
typedef struct {
	StackElem *elems; // Pointer to heap-allocated elements buffer
	int64_t top;	  // Current top index (number of pushed elements)
	int64_t cap;	  // Allocated buffer capacity
} Stack;

static Stack stack_new(void) { return (Stack){.elems = NULL, .top = 0, .cap = 0}; }

static void stack_free(Stack *s) { free(s->elems); }

static bool stack_push(Stack *s, GrB_Index row_idx, GrB_Index col_idx,
					   int32_t nonterm, GrB_Index idx) {
	if (s->top == s->cap) {
		int64_t new_cap = s->cap == 0 ? 64 : s->cap * 2;
		StackElem *tmp = realloc(s->elems, (size_t)new_cap * sizeof(StackElem));
		if (!tmp) {
			return false;
		}
		s->elems = tmp;
		s->cap = new_cap;
	}
	s->elems[s->top++] = (StackElem){row_idx, col_idx, nonterm, idx};
	return true;
}

static StackElem stack_pop(Stack *s) { return s->elems[--s->top]; }

static bool stack_empty(Stack *s) { return s->top == 0; }

/**
 * @brief Represents a terminal production rule (A -> a).
 */
typedef struct {
	int32_t term; // Target terminal index (prod_A of the rule)
} TermRule;

/**
 * @brief Lookup entry containing all terminal production rules for a specific
 * nonterminal.
 */
typedef struct {
	TermRule *term_rules; // Array of terminal production rules
	int32_t term_count;	  // Number of terminal rules for this nonterminal
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

		if (rule.pat == NULL && rule.prod_A != -1 && rule.prod_B == -1) {
			int32_t A = rule.nonterm;
			uint32_t count = rule.indexed_count;
			uint32_t flags = rule.indexed;
			int32_t n = (count > 0) ? (int32_t)count : 1;

			for (int32_t k = 0; k < n; k++) {
				int32_t cur_A = A + ((flags & LAGraph_EWNCF_INDEX_NONTERM) ? k : 0);
				idx[cur_A].term_count++;
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

		if (rule.pat == NULL && rule.prod_A != -1 && rule.prod_B == -1) {
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
 * @brief Lazily retrieves or opens the zero-copy HyperCSR `PathsView` for
 * nonterminal `paths[nonterm_idx]`.
 *
 * @param[out] out_view    Output pointer to the opened `PathsView`.
 * @param[in]  fast_paths  Array of views per nonterminal symbol.
 * @param[in]  paths       Array of original GraphBLAS nonterminal matrices.
 * @param[in]  nonterm_idx Index of the nonterminal symbol.
 * @param[in]  n           Total number of vertices in the graph.
 *
 * @return `GrB_SUCCESS` on success, or an error code on failure.
 */
static inline GrB_Info get_or_open_paths_view(PathsView **out_view,
											  PathsView *fast_paths,
											  GrB_Matrix *paths, int32_t nonterm_idx,
											  GrB_Index n) {
	if (!fast_paths[nonterm_idx].open) {
		GrB_Info info =
			paths_view_open(&fast_paths[nonterm_idx], paths[nonterm_idx], n);
		if (info != GrB_SUCCESS) {
			*out_view = NULL;
			return info;
		}
	}
	*out_view = &fast_paths[nonterm_idx];
	return GrB_SUCCESS;
}

/**
 * @brief Binary rule A -> B C without the left-hand side (stored grouped by A).
 */
typedef struct {
	int32_t B; // First symbol on the right-hand side
	int32_t C; // Second symbol on the right-hand side
} BinRule;

/**
 * @brief CSR-like index of binary rules grouped by the left-hand nonterminal.
 *
 * Rules of nonterminal A are `rules[offset[A] .. offset[A + 1])`.
 */
typedef struct {
	BinRule *rules;	 // Binary rules sorted by left-hand nonterminal
	size_t *offset;	 // nonterms_count + 1 offsets into `rules`
	size_t nonterms; // Number of nonterminals
} BinaryIndex;

static void free_binary_index(BinaryIndex *bi) {
	if (!bi) {
		return;
	}
	free(bi->rules);
	free(bi->offset);
	memset(bi, 0, sizeof(*bi));
}

/**
 * @brief Role-to-nonterminal-id mapping used by grammar patterns.
 *
 * Must stay identical to `pat_role_id` in LAGraph_CFL_AllPaths_adv.
 */
static inline int32_t pat_role_id(const LAGraph_Pat_role *role, int64_t mask,
								  int64_t qr, int64_t qc, int64_t i) {
	if (role->kind == LAGraph_PAT_SYM_PLAIN) {
		return role->base + (int32_t)i;
	}
	return role->base + (int32_t)(mask * role->mask_stride + qr * role->qrow_stride +
								  qc * role->qcol_stride + i);
}

/**
 * @brief State of the two-pass (count / fill) enumeration of concrete binary rules.
 */
typedef struct {
	BinaryIndex *bi;	  // Index being built
	const GrB_Index *nnz; // nvals of paths[a] for every nonterminal a
	size_t nonterms;	  // Number of nonterminals
	size_t *cursor;		  // Write cursors (fill pass only)
	bool fill;			  // false: count pass, true: fill pass
} EmitCtx;

/**
 * @brief Records a concrete rule A -> B C, skipping out-of-range symbols and rules
 * with an empty matrix (such rules could never have produced a middle vertex).
 */
static inline void emit_rule(EmitCtx *c, int32_t A, int32_t B, int32_t C) {
	size_t n = c->nonterms;
	if (A < 0 || B < 0 || C < 0 || (size_t)A >= n || (size_t)B >= n ||
		(size_t)C >= n) {
		return;
	}
	if (c->nnz[A] == 0 || c->nnz[B] == 0 || c->nnz[C] == 0) {
		return;
	}
	if (!c->fill) {
		c->bi->offset[A + 1]++;
	} else {
		c->bi->rules[c->cursor[A]++] = (BinRule){.B = B, .C = C};
	}
}

/**
 * @brief Enumerates all concrete binary rules A -> B C of the grammar.
 *
 * Mirrors the expansion done by LAGraph_CFL_AllPaths_adv: indexed rules are
 * exploded first, then CONV / SCALE / LEFT / RIGHT patterns are unrolled.
 */
static void enumerate_binary_rules(const MRGrammar_t *g, EmitCtx *ctx) {
	for (int64_t r = 0; r < g->rules_count; r++) {
		const LAGraph_rule_EWCNF base = g->rules[r];
		if (base.prod_A == -1 || base.prod_B == -1) {
			continue; // eps / terminal rule
		}

		int64_t reps = (base.indexed_count == 0) ? 1 : (int64_t)base.indexed_count;
		for (int64_t ri = 0; ri < reps; ri++) {
			LAGraph_rule_EWCNF bin_rule = base;
			if (base.indexed_count != 0) {
				bin_rule.indexed_count = 0;
				bin_rule.indexed = 0;
				if (base.nonterm != -1 &&
					(base.indexed & LAGraph_EWNCF_INDEX_NONTERM)) {
					bin_rule.nonterm = base.nonterm + (int32_t)ri;
				}
				if (base.prod_A != -1 &&
					(base.indexed & LAGraph_EWNCF_INDEX_PROD_A)) {
					bin_rule.prod_A = base.prod_A + (int32_t)ri;
				}
				if (base.prod_B != -1 &&
					(base.indexed & LAGraph_EWNCF_INDEX_PROD_B)) {
					bin_rule.prod_B = base.prod_B + (int32_t)ri;
				}
			}

			const LAGraph_rule_Pat *pat = bin_rule.pat;

			if (pat != NULL && pat->kind == LAGraph_PAT_CONV) {
				int64_t N = pat->masks, R = pat->states, K = pat->idx_count;
				if (pat->element == LAGraph_PAT_ALL_ELEMENTS) {
					for (int64_t m = 0; m < N; m++) {
						for (int64_t im = 0; im < N; im++) {
							for (int64_t qs = 0; qs < R; qs++) {
								for (int64_t qm = 0; qm < R; qm++) {
									for (int64_t qe = 0; qe < R; qe++) {
										for (int64_t k = 0; k < K; k++) {
											emit_rule(ctx,
													  pat_role_id(&pat->role[0], m,
																  qs, qe, 0),
													  pat_role_id(&pat->role[1], im,
																  qs, qm, k),
													  pat_role_id(&pat->role[2],
																  m ^ im, qm, qe,
																  0));
										}
									}
								}
							}
						}
					}
				} else {
					int64_t el = pat->element;
					uint32_t st_count = pat->state ? pat->state->nnz : (uint32_t)R;
					for (int64_t m = 0; m < N; m++) {
						for (uint32_t st = 0; st < st_count; st++) {
							int64_t qs =
								pat->state ? pat->state->entries[st].row : st;
							int64_t qmid =
								pat->state ? pat->state->entries[st].col : st;
							for (int64_t qe = 0; qe < R; qe++) {
								emit_rule(
									ctx, pat_role_id(&pat->role[0], m, qs, qe, 0),
									bin_rule.prod_A,
									pat_role_id(&pat->role[2], m ^ el, qmid, qe, 0));
							}
						}
					}
				}
				continue;
			}

			if (pat != NULL && pat->kind == LAGraph_PAT_SCALE) {
				int64_t N = pat->masks, R = pat->states;
				for (int64_t m = 0; m < N; m++) {
					for (int64_t qs = 0; qs < R; qs++) {
						for (int64_t qe = 0; qe < R; qe++) {
							emit_rule(ctx, pat_role_id(&pat->role[0], m, qs, qe, 0),
									  bin_rule.prod_A,
									  pat_role_id(&pat->role[2], m, qs, qe, 0));
						}
					}
				}
				continue;
			}

			if (pat != NULL &&
				(pat->kind == LAGraph_PAT_LEFT || pat->kind == LAGraph_PAT_RIGHT)) {
				int64_t N = pat->masks, R = pat->states;
				int64_t off = bin_rule.nonterm - pat->role[0].base;
				if (off < 0 || off >= (int64_t)pat->idx_count) {
					continue;
				}
				for (int64_t im = 0; im < N; im++) {
					for (int64_t qs = 0; qs < R; qs++) {
						for (int64_t qm = 0; qm < R; qm++) {
							int32_t cA = pat_role_id(&pat->role[0], im, qs, qm, off);
							if (pat->kind == LAGraph_PAT_LEFT) {
								emit_rule(
									ctx, cA, bin_rule.prod_A,
									pat_role_id(&pat->role[2], im, qs, qm, off));
							} else {
								emit_rule(ctx, cA,
										  pat_role_id(&pat->role[1], im, qs, qm, 0),
										  bin_rule.prod_B);
							}
						}
					}
				}
				continue;
			}

			// Plain binary rule
			emit_rule(ctx, bin_rule.nonterm, bin_rule.prod_A, bin_rule.prod_B);
		}
	}
}

/**
 * @brief Builds the per-nonterminal index of concrete binary rules directly from
 * the grammar (the rule table returned by the AllPaths algorithm is not used).
 *
 * @param[out] bi   Resulting index (free with free_binary_index).
 * @param[in]  g    Grammar.
 * @param[in]  paths AllPaths matrices, used only to skip rules with empty
 * matrices.
 */
static GrB_Info build_binary_index(BinaryIndex *bi, const MRGrammar_t *g,
								   GrB_Matrix *paths) {
	size_t nonterms = (size_t)g->nonterms_count;
	memset(bi, 0, sizeof(*bi));
	bi->nonterms = nonterms;

	GrB_Index *nnz = calloc(nonterms ? nonterms : 1, sizeof(GrB_Index));
	size_t *cursor = calloc(nonterms ? nonterms : 1, sizeof(size_t));
	bi->offset = calloc(nonterms + 1, sizeof(size_t));
	if (!nnz || !cursor || !bi->offset) {
		free(nnz);
		free(cursor);
		free_binary_index(bi);
		return GrB_OUT_OF_MEMORY;
	}
	for (size_t a = 0; a < nonterms; a++) {
		if (paths[a] != NULL) {
			GrB_Matrix_nvals(&nnz[a], paths[a]);
		}
	}

	EmitCtx ctx = {
		.bi = bi, .nnz = nnz, .nonterms = nonterms, .cursor = cursor, .fill = false};
	enumerate_binary_rules(g, &ctx); // count pass

	for (size_t a = 0; a < nonterms; a++) {
		bi->offset[a + 1] += bi->offset[a];
	}
	size_t total = bi->offset[nonterms];
	bi->rules = malloc((total ? total : 1) * sizeof(BinRule));
	if (!bi->rules) {
		free(nnz);
		free(cursor);
		free_binary_index(bi);
		return GrB_OUT_OF_MEMORY;
	}
	memcpy(cursor, bi->offset, nonterms * sizeof(size_t));

	ctx.fill = true;
	enumerate_binary_rules(g, &ctx); // fill pass

	free(nnz);
	free(cursor);
	return GrB_SUCCESS;
}

/**
 * @brief Attempts to push a path derivation step $(i, j)$ for nonterminal $A$
 * onto the DFS stack.
 *
 * Performs element presence lookup in the corresponding `PathsView` and checks
 * if the entry was already visited using the bitset.
 *
 * @param[in,out] st         Pointer to the active stack.
 * @param[in,out] fast_paths Array of cached `PathsView` objects.
 * @param[in]     paths      Array of raw GraphBLAS matrices.
 * @param[in]     row_idx    Source vertex index (row i).
 * @param[in]     col_idx    Target vertex index (column j).
 * @param[in]     nonterm    Nonterminal symbol ID.
 * @param[in]     n          Total number of vertices in the graph.
 *
 * @return `GrB_SUCCESS` on success, or an error code on failure.
 */
static GrB_Info push_unvisited_path_step(Stack *st, PathsView *fast_paths,
										 GrB_Matrix *paths, GrB_Index row_idx,
										 GrB_Index col_idx, int32_t nonterm,
										 GrB_Index n) {
	GrB_Info info = GrB_SUCCESS;
	PathsView *v = NULL;

	GRB_TRY(get_or_open_paths_view(&v, fast_paths, paths, nonterm, n));

	GrB_Index idx;
	if (!paths_view_find_elem(v, row_idx, col_idx, &idx)) {
		return GrB_SUCCESS;
	}
	if (test_and_set(v->visited, idx)) {
		return GrB_SUCCESS;
	}

	if (!stack_push(st, row_idx, col_idx, nonterm, idx)) {
		info = GrB_OUT_OF_MEMORY;
		goto cleanup;
	}

	return GrB_SUCCESS;

cleanup:
	return info;
}

/**
 * @brief Marks an already located element as visited and pushes it onto the
 * stack.
 *
 * Unlike push_unvisited_path_step, the element must have been found beforehand
 * (so that both halves of a binary split can be validated before pushing
 * either).
 *
 * @return `GrB_SUCCESS` or `GrB_OUT_OF_MEMORY`.
 */
static inline GrB_Info claim_and_push(Stack *st, PathsView *v, GrB_Index row_idx,
									  GrB_Index col_idx, int32_t nonterm,
									  GrB_Index idx) {
	if (test_and_set(v->visited, idx)) {
		return GrB_SUCCESS;
	}
	return stack_push(st, row_idx, col_idx, nonterm, idx) ? GrB_SUCCESS
														  : GrB_OUT_OF_MEMORY;
}

/**
 * @brief Bundles runtime context and resources required for edge extraction.
 *
 * Passed across helper routines to minimize parameter passing during DFS stack
 * draining.
 */
typedef struct {
	Stack *stack;		   // Active DFS stack
	PathsView *fast_paths; // Array of HyperCSR views per nonterminal symbol
	GrB_Matrix *paths;	   // Raw GraphBLAS matrices for nonterminals
	GrB_Matrix *out;	   // Output GraphBLAS matrices for extracted terminal edges
	GrB_Matrix *adj_matrices;		 // Input graph adjacency matrices for terminals
	const BinaryIndex *bin_index;	 // Binary rules (A -> B C) grouped by A
	const NontermRules *grammar_idx; // Index of terminal grammar rules (A -> a)
	const int32_t *grammar_to_straight; // Interleaved to linear index mapping
										// table for terminals
	MRGrammar_t grammar;				// Context-free grammar definition
	GrB_Index n;						// Number of vertices in the graph
} ExtractionContext;

/**
 * @brief Drains current stack completely in DFS order to reconstruct derivation
 * trees.
 *
 * Recursively pops derivation elements, splits binary rules $A \to B C$, and
 * records terminal edges into output matrices when terminal derivations $A \to
 * a$ are reached.
 *
 * @param[in,out] ctx Pointer to the active extraction context.
 *
 * @return `GrB_SUCCESS` on success, or an error code on failure.
 */
static GrB_Info drain_extraction_stack(ExtractionContext *ctx) {
	GrB_Info info = GrB_SUCCESS;

	while (!stack_empty(ctx->stack)) {
		StackElem cur = stack_pop(ctx->stack);
		GrB_Index i = cur.row_idx;
		GrB_Index j = cur.col_idx;
		int32_t A = cur.nonterm;
		GrB_Index elem_idx = cur.elem_idx;

		PathsView *cur_view = &ctx->fast_paths[A];
		AllPathsElem elem = cur_view->iso ? cur_view->Ax[0] : cur_view->Ax[elem_idx];

		const NontermRules *rules = &ctx->grammar_idx[A];

		const BinaryIndex *bi = ctx->bin_index;
		size_t r_begin = bi->offset[A];
		size_t r_end = bi->offset[A + 1];

		// Iterate over intermediate vertices
		for (size_t ki = 0; ki < elem.n; ki++) {
			GrB_Index mid =
				(elem.n == 1) ? elem.data.single_elem : elem.data.middle[ki];

			if (mid == GrB_INDEX_MAX) {
				// Terminal edge: A -> term
				for (int32_t ri = 0; ri < rules->term_count; ri++) {
					int32_t term = rules->term_rules[ri].term;
					int32_t term_idx = term - (int32_t)ctx->grammar.nonterms_count;

					bool has_edge = false;
					GrB_Info err = GrB_Matrix_extractElement_BOOL(
						&has_edge,
						ctx->adj_matrices[ctx->grammar.nonterms_count + term_idx], i,
						j);
					if (err == GrB_SUCCESS && has_edge) {
						int32_t straight_idx = ctx->grammar_to_straight[term_idx];
						GRB_TRY(GrB_Matrix_setElement_BOOL(ctx->out[straight_idx],
														   true, i, j));
					}
				}
				continue;
			}

			// Binary split A -> B C via mid. The rule id is not stored in the
			// element, so every rule of A is tested: it is a valid split iff
			// both B(i, mid) and C(mid, j) are present.
			for (size_t r = r_begin; r < r_end; r++) {
				BinRule br = bi->rules[r];
				PathsView *vB = NULL, *vC = NULL;
				GrB_Index idxB, idxC;

				GRB_TRY(get_or_open_paths_view(&vB, ctx->fast_paths, ctx->paths,
											   br.B, ctx->n));
				if (!paths_view_find_elem(vB, i, mid, &idxB)) {
					continue;
				}
				GRB_TRY(get_or_open_paths_view(&vC, ctx->fast_paths, ctx->paths,
											   br.C, ctx->n));
				if (!paths_view_find_elem(vC, mid, j, &idxC)) {
					continue;
				}

				GRB_TRY(claim_and_push(ctx->stack, vB, i, mid, br.B, idxB));
				GRB_TRY(claim_and_push(ctx->stack, vC, mid, j, br.C, idxC));
			}
		}
	}
	return GrB_SUCCESS;

cleanup:
	return info;
}

GrB_Info extract_edges_from_outputs(GrB_Matrix *out, GrB_Matrix *paths,
									GrB_Matrix *adj_matrices, MRGrammar_t grammar,
									GrB_Index n, int64_t n_par, int64_t n_bra,
									const MRGrammarConfig *config,
									const TargetPath *target_path) {
	GrB_Info info = GrB_SUCCESS;

	NontermRules *grammar_idx = NULL;
	BinaryIndex bin_index = {0};
	Stack stack = stack_new();
	int32_t *grammar_to_straight = NULL;
	PathsView *fast_paths = NULL;

	// Pre-initialize out pointers to NULL for safe cleanup
	for (int64_t t = 0; t < grammar.terms_count; t++) {
		out[t] = NULL;
	}

	// Initialize out & force materialize adj_matrices
	for (int64_t t = 0; t < grammar.terms_count; t++) {
		GRB_TRY(GrB_Matrix_new(&out[t], GrB_BOOL, n, n));
		GRB_TRY(GrB_Matrix_wait(adj_matrices[grammar.nonterms_count + t],
								GrB_MATERIALIZE));
	}

	grammar_to_straight = malloc((size_t)grammar.terms_count * sizeof(int32_t));
	if (!grammar_to_straight) {
		info = GrB_OUT_OF_MEMORY;
		goto cleanup;
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

	GRB_TRY(build_binary_index(&bin_index, &grammar, paths));

	fast_paths = calloc((size_t)grammar.nonterms_count, sizeof(PathsView));
	if (!fast_paths) {
		info = GrB_OUT_OF_MEMORY;
		goto cleanup;
	}

	ExtractionContext ctx = {.stack = &stack,
							 .fast_paths = fast_paths,
							 .paths = paths,
							 .out = out,
							 .adj_matrices = adj_matrices,
							 .bin_index = &bin_index,
							 .grammar_idx = grammar_idx,
							 .grammar_to_straight = grammar_to_straight,
							 .grammar = grammar,
							 .n = n};

	// Stack initialization & execution
	if (target_path != NULL) {
		// On-demand analysis
		GRB_TRY(push_unvisited_path_step(&stack, fast_paths, paths, target_path->src,
										 target_path->tgt, NT_START, n));
		GRB_TRY(drain_extraction_stack(&ctx));
	} else {
		// Default: all pairs from paths[NT_START]
		PathsView *start_view = NULL;
		GRB_TRY(get_or_open_paths_view(&start_view, fast_paths, paths, NT_START, n));

		for (GrB_Index rr = 0; rr < start_view->nvec; rr++) {
			GrB_Index r = start_view->Ah[rr];
			for (GrB_Index e = start_view->Ap[rr]; e < start_view->Ap[rr + 1]; e++) {
				GrB_Index c = start_view->Aj[e];
				if (r != c) {
					GRB_TRY(push_unvisited_path_step(&stack, fast_paths, paths, r, c,
													 NT_START, n));
					GRB_TRY(drain_extraction_stack(&ctx));
				}
			}
		}
	}

	for (int64_t t = 0; t < grammar.terms_count; t++) {
		if (out[t]) {
			GrB_Matrix_wait(out[t], GrB_MATERIALIZE);
		}
	}

cleanup:
	if (info != GrB_SUCCESS) {
		for (int64_t t = 0; t < grammar.terms_count; t++) {
			if (out[t]) {
				GrB_Matrix_free(&out[t]);
			}
		}
	}
	if (fast_paths) {
		for (int32_t a = 0; a < grammar.nonterms_count; a++) {
			GrB_Info close_info = paths_view_close(&fast_paths[a]);
			if (info == GrB_SUCCESS && close_info != GrB_SUCCESS) {
				info = close_info;
			}
		}
		free(fast_paths);
	}
	free(grammar_to_straight);
	free_grammar_index(grammar_idx, grammar.nonterms_count);
	free_binary_index(&bin_index);
	stack_free(&stack);
	return info;
}
