#include "LAGraph.h"
#include "approximation/approximation.h"
#include "parse_utils.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static char msg[LAGRAPH_MSG_LEN];

#define RUN_FIG_TEST(fig_name, ext_expected)                                        \
	do {                                                                            \
		char g[1024], e[1024];                                                      \
		snprintf(g, 1024, "%s/" #fig_name "/graph.g", TEST_DATA_DIR);               \
		snprintf(e, 1024, "%s/" #fig_name "/" ext_expected, TEST_DATA_DIR);         \
		run_test(g, e);                                                             \
	} while (0)

// Utils

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

	GrB_Index *rows = malloc(nvals * sizeof(GrB_Index));
	GrB_Index *cols = malloc(nvals * sizeof(GrB_Index));
	bool *vals = malloc(nvals * sizeof(bool));
	GrB_Matrix_extractTuples_BOOL(rows, cols, vals, &nvals, m);

	Pair *pairs = malloc(nvals * sizeof(Pair));
	for (GrB_Index i = 0; i < nvals; i++)
		pairs[i] = (Pair){rows[i], cols[i]};

	free(rows);
	free(cols);
	free(vals);

	qsort(pairs, nvals, sizeof(Pair), cmp_pair);
	*out = pairs;
	return nvals;
}

static size_t load_expected(const char *path, Pair **out) {
	FILE *f = fopen(path, "r");
	assert(f != NULL);

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
	assert(actual_count == expected_count);
	for (size_t i = 0; i < actual_count; i++) {
		assert(actual[i].row == expected[i].row);
		assert(actual[i].col == expected[i].col);
	}
}

static void run_test(const char *graph_path, const char *expected_path) {
	MRGraph graph = {0};
	GrB_Info info = parse_graph(graph_path, &graph);
	assert(info == GrB_SUCCESS);

	GrB_Matrix result = NULL;
	info = get_under_approx(&graph, &result);
	assert(info == GrB_SUCCESS);

	Pair *actual = NULL;
	size_t actual_count = extract_pairs(result, &actual);
	Pair *expected = NULL;
	size_t expected_count = load_expected(expected_path, &expected);

	assert_pairs_equal(actual, actual_count, expected, expected_count);

	free(actual);
	free(expected);
	GrB_Matrix_free(&result);
	mr_graph_free(&graph);
}

int main(void) {
	LAGraph_Init(msg);

	RUN_FIG_TEST(figure5, "under_approx.txt");
	RUN_FIG_TEST(figure9, "under_approx.txt");
	RUN_FIG_TEST(figure10, "under_approx.txt");
	RUN_FIG_TEST(loozfon, "under_approx.txt");
	RUN_FIG_TEST(faketaobao, "under_approx.txt");
	RUN_FIG_TEST(jollyserv, "under_approx.txt");
	RUN_FIG_TEST(zertsecurity, "under_approx.txt");
	RUN_FIG_TEST(fakebanker, "under_approx.txt");
	RUN_FIG_TEST(uranai, "under_approx.txt");

	LAGraph_Finalize(msg);
	return 0;
}
