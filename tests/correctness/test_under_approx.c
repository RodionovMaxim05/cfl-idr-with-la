#include "utils.h"

#include "approximation/approximation.h"
#include "graph/parse_utils.h"

static char msg[LAGRAPH_MSG_LEN];

#define RUN_UNDER_APPROX_TEST(folder, ext_expected)                                 \
	do {                                                                            \
		char g[1024], e[1024];                                                      \
		snprintf(g, 1024, "%s/" #folder "/graph.g", TEST_DATA_DIR);                 \
		snprintf(e, 1024, "%s/" #folder "/" ext_expected, TEST_DATA_DIR);           \
		printf("Running test: under approx on " #folder "...\n");                   \
		fflush(stdout);                                                             \
		run_test_logic(g, e);                                                       \
	} while (0)

static void run_test_logic(const char *graph_path, const char *expected_path) {
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

	RUN_UNDER_APPROX_TEST(figure5, "under_approx.txt");
	RUN_UNDER_APPROX_TEST(figure9, "under_approx.txt");
	RUN_UNDER_APPROX_TEST(figure10, "under_approx.txt");
	RUN_UNDER_APPROX_TEST(loozfon, "under_approx.txt");
	RUN_UNDER_APPROX_TEST(faketaobao, "under_approx.txt");
	RUN_UNDER_APPROX_TEST(jollyserv, "under_approx.txt");
	RUN_UNDER_APPROX_TEST(zertsecurity, "under_approx.txt");
	RUN_UNDER_APPROX_TEST(fakebanker, "under_approx.txt");
	RUN_UNDER_APPROX_TEST(uranai, "under_approx.txt");

	LAGraph_Finalize(msg);
	return 0;
}
