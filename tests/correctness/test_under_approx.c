#include "utils.h"

static char msg[LAGRAPH_MSG_LEN];

#define RUN_UNDER_APPROX_TEST(folder, valueflow, ext_expected)                      \
	do {                                                                            \
		char g[1024], e[1024];                                                      \
		snprintf(g, 1024, "%s/" #folder "/graph.g", TEST_DATA_DIR);                 \
		snprintf(e, 1024, "%s/" #folder "/" ext_expected, TEST_DATA_DIR);           \
		printf("Running test: under approx on " #folder "%s...\n",                  \
			   valueflow ? " (valueflow)" : "");                                    \
		fflush(stdout);                                                             \
		run_test_logic(g, valueflow, e);                                            \
	} while (0)

static void run_test_logic(const char *graph_path, bool valueflow,
						   const char *expected_path) {
	IdrGraph graph = {0};
	GrB_Info info = parse_graph(graph_path, &graph);
	assert(info == GrB_SUCCESS);

	IdrGraph new_graph = {0};
	if (valueflow) {
		info = idr_remove_valueflow_unreachable(&new_graph, &graph);
		assert(info == GrB_SUCCESS);
		idr_graph_free(&graph);
	} else {
		new_graph = graph;
	}

	GrB_Matrix result = NULL;
	info = idr_get_under_approx(&result, &new_graph, valueflow);
	assert(info == GrB_SUCCESS);

	Pair *actual = NULL;
	size_t actual_count = extract_pairs(result, &actual);
	Pair *expected = NULL;
	size_t expected_count = load_expected(expected_path, &expected);

	assert_pairs_equal(actual, actual_count, expected, expected_count);

	free(actual);
	free(expected);
	GrB_Matrix_free(&result);
	idr_graph_free(&new_graph);
}

int main(void) {
	LAGraph_Init(msg);

	RUN_UNDER_APPROX_TEST(figure5, false, "under_approx.txt");
	RUN_UNDER_APPROX_TEST(figure9, false, "under_approx.txt");
	RUN_UNDER_APPROX_TEST(figure10, false, "under_approx.txt");
	RUN_UNDER_APPROX_TEST(loozfon, false, "under_approx.txt");
	RUN_UNDER_APPROX_TEST(faketaobao, false, "under_approx.txt");
	RUN_UNDER_APPROX_TEST(jollyserv, false, "under_approx.txt");
	RUN_UNDER_APPROX_TEST(zertsecurity, false, "under_approx.txt");
	RUN_UNDER_APPROX_TEST(fakebanker, false, "under_approx.txt");
	RUN_UNDER_APPROX_TEST(uranai, false, "under_approx.txt");

	// Value-flow analysis
	RUN_UNDER_APPROX_TEST(xz, /*valueflow=*/true, "paths.txt");
	RUN_UNDER_APPROX_TEST(nab, /*valueflow=*/true, "paths.txt");
	RUN_UNDER_APPROX_TEST(leela, /*valueflow=*/true, "paths.txt");

	LAGraph_Finalize(msg);
	return 0;
}
