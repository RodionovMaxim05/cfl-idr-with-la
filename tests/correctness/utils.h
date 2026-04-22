#include "LAGraph.h"
#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct {
	GrB_Index row;
	GrB_Index col;
} Pair;

static int cmp_pair(const void *a, const void *b) {
	const Pair *pa = (const Pair *)a;
	const Pair *pb = (const Pair *)b;
	if (pa->row != pb->row)
		return (pa->row > pb->row) - (pa->row < pb->row);
	return (pa->col > pb->col) - (pa->col < pb->col);
}

static size_t extract_pairs(GrB_Matrix m, Pair **out) {
	GrB_Matrix_wait(m, GrB_MATERIALIZE);
	GrB_Index nvals = 0;
	GrB_Matrix_nvals(&nvals, m);

	if (nvals == 0) {
		*out = NULL;
		return 0;
	}

	GrB_Index *rows = malloc(nvals * sizeof(GrB_Index));
	GrB_Index *cols = malloc(nvals * sizeof(GrB_Index));
	void *vals = malloc(nvals * sizeof(double));
	GrB_Matrix_extractTuples_BOOL(rows, cols, (bool *)vals, &nvals, m);

	Pair *pairs = malloc(nvals * sizeof(Pair));
	for (GrB_Index i = 0; i < nvals; i++) {
		pairs[i] = (Pair){rows[i], cols[i]};
	}

	free(rows);
	free(cols);
	free(vals);

	qsort(pairs, nvals, sizeof(Pair), cmp_pair);
	*out = pairs;
	return (size_t)nvals;
}

static size_t load_expected(const char *path, Pair **out) {
	FILE *f = fopen(path, "r");
	if (!f) {
		fprintf(stderr, "File not found: %s\n", path);
		assert(f != NULL);
	}

	size_t capacity = 128;
	size_t count = 0;
	Pair *pairs = malloc(capacity * sizeof(Pair));

	GrB_Index row, col;
	while (fscanf(f, "%lu %lu", &row, &col) == 2) {
		if (count == capacity) {
			capacity *= 2;
			pairs = realloc(pairs, capacity * sizeof(Pair));
		}
		pairs[count++] = (Pair){row, col};
	}
	fclose(f);

	qsort(pairs, count, sizeof(Pair), cmp_pair);
	*out = pairs;
	return count;
}

static void assert_pairs_equal(Pair *actual, size_t actual_count, Pair *expected,
							   size_t expected_count) {
	if (actual_count != expected_count) {
		fprintf(stderr, "Size mismatch: expected %zu, got %zu\n", expected_count,
				actual_count);
	}
	assert(actual_count == expected_count);
	for (size_t i = 0; i < actual_count; i++) {
		assert(actual[i].row == expected[i].row);
		assert(actual[i].col == expected[i].col);
	}
}
