#include "extract_edges.h"
#include "grammar/grammar.h"
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
	Stack s = {.elems = NULL, .top = 0, .capacity = 0};
	return s;
}

static void stack_free(Stack *s) { free(s->elems); }

static void stack_push(Stack *s, GrB_Index i, GrB_Index j, int32_t A) {
	if (s->top == s->capacity) {
		s->capacity = s->capacity == 0 ? 16 : s->capacity * 2;
		StackElem *tmp = realloc(s->elems, s->capacity * sizeof(StackElem));
		if (tmp == NULL) {
			free(s->elems);
			s->elems = NULL;
		}
		s->elems = tmp;
	}
	s->elems[s->top++] = (StackElem){i, j, A};
}

static StackElem stack_pop(Stack *s) { return s->elems[--s->top]; }

static bool stack_empty(Stack *s) { return s->top == 0; }

GrB_Info extractEdgesFromOutputs(GrB_Matrix *paths, GrB_Matrix *adj_matrices,
								 MRGrammar_t grammar, GrB_Index n, GrB_Index start_i,
								 GrB_Index start_j, GrB_Matrix *result_matrices,
								 char *msg) {
	GrB_Info info = GrB_SUCCESS;

	// Initialize result_matrices
	for (int64_t t = 0; t < grammar.terms_count; t++) {
		GrB_Matrix_new(&result_matrices[t], GrB_BOOL, n, n);
	}

	// visited - one Boolean matrix for each nonterminal
	GrB_Matrix *visited =
		(GrB_Matrix *)malloc(grammar.nonterms_count * sizeof(GrB_Matrix));
	for (int32_t a = 0; a < grammar.nonterms_count; a++) {
		GrB_Matrix_new(&visited[a], GrB_BOOL, n, n);
	}

	Stack stack = stack_new();

	// Stack initialization
	if (start_i != GrB_INDEX_MAX && start_j != GrB_INDEX_MAX) {
		// On-demand analysis

		stack_push(&stack, start_i, start_j, NT_S);
	} else {
		// Default: all pairs from paths[NT_S]

		GrB_Index nnz = 0;
		GrB_Matrix_nvals(&nnz, paths[NT_S]);
		if (nnz == 0) {
			goto cleanup;
		}

		GrB_Index *rows = NULL;
		GrB_Index *cols = NULL;
		void *val_void = NULL;
		LAGraph_Malloc((void **)&rows, nnz, sizeof(GrB_Index), msg);
		LAGraph_Malloc((void **)&cols, nnz, sizeof(GrB_Index), msg);
		LAGraph_Malloc(&val_void, nnz, sizeof(AllPathsElem), msg);

		GrB_Matrix_extractTuples(rows, cols, val_void, &nnz, paths[NT_S]);

		for (GrB_Index k = 0; k < nnz; k++) {
			if (rows[k] != cols[k]) {
				stack_push(&stack, rows[k], cols[k], NT_S);
			}
		}

		LAGraph_Free((void **)&rows, msg);
		LAGraph_Free((void **)&cols, msg);
		LAGraph_Free(&val_void, msg);
	}

	while (!stack_empty(&stack)) {
		StackElem cur_node = stack_pop(&stack);
		GrB_Index i = cur_node.row_idx;
		GrB_Index j = cur_node.col_idx;
		int32_t A = cur_node.nonterm;

		bool seen = false;
		info = GrB_Matrix_extractElement_BOOL(&seen, visited[A], i, j);
		if (info == GrB_SUCCESS && seen) {
			continue;
		}

		// Mark as visited
		GrB_Matrix_setElement_BOOL(visited[A], true, i, j);

		AllPathsElem elem;
		GrB_Info elem_info = GrB_Matrix_extractElement_UDT(&elem, paths[A], i, j);

		if (elem_info == GrB_NO_VALUE) {
			continue;
		}

		// Iterate over intermediate vertices
		for (GrB_Index k = 0; k < elem.n; k++) {
			GrB_Index m = elem.middle[k];

			if (m == GrB_INDEX_MAX) {
				// Add single-edge or empty paths

				for (int64_t r = 0; r < grammar.rules_count; r++) {
					// Looking for rules A -> terminal
					if (grammar.rules[r].nonterm != A) {
						continue;
					}
					if (grammar.rules[r].prod_B != -1) {
						continue;
					}
					if (grammar.rules[r].prod_A == -1 &&
						grammar.rules[r].prod_B == -1) {
						continue;
					}

					int32_t term = grammar.rules[r].prod_A;

					bool has_edge = false;
					GrB_Matrix_extractElement_BOOL(&has_edge, adj_matrices[term], i,
												   j);

					if (has_edge) {
						GrB_Matrix_setElement_BOOL(result_matrices[term], true, i,
												   j);
					}
				}
			} else {
				// Add to result the concatenated paths from i to m and from m to j

				for (int64_t r = 0; r < grammar.rules_count; r++) {
					// Looking for rules A -> B C
					if (grammar.rules[r].nonterm != A) {
						continue;
					}
					if (grammar.rules[r].prod_B == -1) {
						continue;
					}

					int32_t B = grammar.rules[r].prod_A;
					int32_t C = grammar.rules[r].prod_B;

					// Check that both subpaths exist in paths
					AllPathsElem dummy;
					GrB_Info b_info =
						GrB_Matrix_extractElement_UDT(&dummy, paths[B], i, m);
					GrB_Info c_info =
						GrB_Matrix_extractElement_UDT(&dummy, paths[C], m, j);

					if (b_info == GrB_NO_VALUE || c_info == GrB_NO_VALUE) {
						continue;
					}

					stack_push(&stack, i, m, B);
					stack_push(&stack, m, j, C);
				}
			}
		}
	}

cleanup:
	for (int32_t a = 0; a < grammar.nonterms_count; a++) {
		GrB_Matrix_free(&visited[a]);
	}
	free((void *)visited);

	stack_free(&stack);
	return GrB_SUCCESS;
}
