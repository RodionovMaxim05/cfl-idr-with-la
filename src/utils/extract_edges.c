#include "extract_edges.h"

#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

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
		if (tmp == NULL) {
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
	int32_t B; // left nonterminal (prod_A of the rule)
	int32_t C; // right nonterminal (prod_B of the rule)
} BinaryRule;

typedef struct {
	TermRule *term_rules;
	int32_t term_count;
	BinaryRule *binary_rules;
	int32_t binary_count;
} NontermRules;

static NontermRules *build_grammar_index(const MRGrammar_t *g) {
	NontermRules *idx = calloc((size_t)g->nonterms_count, sizeof(NontermRules));
	if (!idx) {
		return NULL;
	}

	// Count pass
	for (int64_t r = 0; r < g->rules_count; r++) {
		int32_t A = g->rules[r].nonterm;
		int32_t prod_A = g->rules[r].prod_A;
		int32_t prod_B = g->rules[r].prod_B;

		if (prod_A == -1) {
			// epsilon - ignore
			continue;
		}
		if (prod_B == -1) {
			// Terminal rule
			idx[A].term_count++;
		} else {
			// Binary rule
			idx[A].binary_count++;
		}
	}

	// Alloc pass
	for (int32_t a = 0; a < g->nonterms_count; a++) {
		if (idx[a].term_count) {
			idx[a].term_rules = malloc((size_t)idx[a].term_count * sizeof(TermRule));
		}
		if (idx[a].binary_count) {
			idx[a].binary_rules =
				malloc((size_t)idx[a].binary_count * sizeof(BinaryRule));
		}
		// Reset - reuse as write cursor
		idx[a].term_count = 0;
		idx[a].binary_count = 0;
	}

	// Fill pass
	for (int64_t r = 0; r < g->rules_count; r++) {
		int32_t A = g->rules[r].nonterm;
		int32_t prod_A = g->rules[r].prod_A;
		int32_t prod_B = g->rules[r].prod_B;

		if (prod_A == -1) {
			continue;
		}
		if (prod_B == -1) {
			idx[A].term_rules[idx[A].term_count++].term = prod_A;
		} else {
			BinaryRule *br = &idx[A].binary_rules[idx[A].binary_count++];
			br->B = prod_A;
			br->C = prod_B;
		}
	}

	return idx;
}

static void free_grammar_index(NontermRules *idx, int32_t nonterms_count) {
	if (!idx) {
		return;
	}
	for (int32_t a = 0; a < nonterms_count; a++) {
		free(idx[a].term_rules);
		free(idx[a].binary_rules);
	}
	free(idx);
}

#define VISITED_IDX(A, i, j, n) ((size_t)(A) * (n) * (n) + (size_t)(i) * (n) + (j))
#define VISITED_GET(v, A, i, j, n) ((v)[VISITED_IDX(A, i, j, n)])
#define VISITED_SET(v, A, i, j, n) ((v)[VISITED_IDX(A, i, j, n)] = 1)

GrB_Info extractEdgesFromOutputs(GrB_Matrix *paths, GrB_Matrix *adj_matrices,
								 MRGrammar_t grammar, GrB_Index n,
								 TargetPath *target_path,
								 GrB_Matrix *result_matrices, char *msg) {
	GrB_Info info = GrB_SUCCESS;

	// Initialize result_matrices
	for (int64_t t = 0; t < grammar.terms_count; t++) {
		GrB_Matrix_new(&result_matrices[t], GrB_BOOL, n, n);
	}

	// Build grammar index
	NontermRules *grammar_idx = build_grammar_index(&grammar);
	if (!grammar_idx) {
		return GrB_OUT_OF_MEMORY;
	}

	// visited - flat array
	size_t visited_size = (size_t)grammar.nonterms_count * n * n;
	uint8_t *visited = calloc(visited_size, 1);
	if (!visited) {
		free_grammar_index(grammar_idx, grammar.nonterms_count);
		return GrB_OUT_OF_MEMORY;
	}

	Stack stack = stack_new();

	// Stack initialization
	if (target_path != NULL) {
		// On-demand analysis

		stack_push(&stack, target_path->src, target_path->tgt, NT_START);
	} else {
		// Default: all pairs from paths[NT_START]

		GrB_Index nnz = 0;
		GrB_Matrix_nvals(&nnz, paths[NT_START]);
		if (nnz == 0) {
			goto cleanup;
		}

		GrB_Index *rows = NULL;
		GrB_Index *cols = NULL;
		void *val_void = NULL;
		LAGraph_Malloc((void **)&rows, nnz, sizeof(GrB_Index), msg);
		LAGraph_Malloc((void **)&cols, nnz, sizeof(GrB_Index), msg);
		LAGraph_Malloc(&val_void, nnz, sizeof(AllPathsElem), msg);

		GrB_Matrix_extractTuples(rows, cols, val_void, &nnz, paths[NT_START]);

		for (GrB_Index k = 0; k < nnz; k++) {
			if (rows[k] != cols[k]) {
				stack_push(&stack, rows[k], cols[k], NT_START);
			}
		}

		LAGraph_Free((void **)&rows, msg);
		LAGraph_Free((void **)&cols, msg);
		LAGraph_Free(&val_void, msg);
	}

	// Main DFS loop
	while (!stack_empty(&stack)) {
		StackElem cur = stack_pop(&stack);
		GrB_Index i = cur.row_idx;
		GrB_Index j = cur.col_idx;
		int32_t A = cur.nonterm;

		if (VISITED_GET(visited, A, i, j, n)) {
			continue;
		}
		VISITED_SET(visited, A, i, j, n);

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
						GrB_Matrix_setElement_BOOL(result_matrices[term], true, i,
												   j);
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

					// Only push if not already visited - avoids stacking the same
					// pair multiple times
					if (!VISITED_GET(visited, B, i, mid, n)) {
						stack_push(&stack, i, mid, B);
					}
					if (!VISITED_GET(visited, C, mid, j, n)) {
						stack_push(&stack, mid, j, C);
					}
				}
			}
		}
	}

cleanup:
	free(visited);
	free_grammar_index(grammar_idx, grammar.nonterms_count);
	stack_free(&stack);
	return GrB_SUCCESS;
}
