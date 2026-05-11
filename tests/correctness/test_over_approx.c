#include "utils.h"

static char msg[LAGRAPH_MSG_LEN];

#define RUN_OVER_APPROX_TEST(folder, grammar_enum, valueflow, expected_count,       \
							 expected_file)                                         \
	do {                                                                            \
		char g_path[1024], e_path[1024];                                            \
		snprintf(g_path, 1024, "%s/" #folder "/graph.g", TEST_DATA_DIR);            \
		snprintf(e_path, 1024, "%s/" #folder "/" expected_file, TEST_DATA_DIR);     \
		printf("Running test: " #grammar_enum " on " #folder "%s...\n",             \
			   valueflow ? " (valueflow)" : "");                                    \
		fflush(stdout);                                                             \
		run_test_logic(g_path, grammar_enum, valueflow, expected_count, e_path);    \
	} while (0)

static void run_test_logic(const char *graph_path, IdrGrammarType grammar_type,
						   bool valueflow, size_t expected_size,
						   const char *expected_path) {
	IdrGraph graph = {0};
	GrB_Info info = parse_graph(&graph, graph_path);
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
	info = idr_get_over_approx(&result, &new_graph, grammar_type, NULL, valueflow,
							   true);
	assert(info == GrB_SUCCESS);

	Pair *actual = NULL;
	size_t actual_count = extract_pairs(result, &actual);

	if (actual_count != expected_size) {
		fprintf(stderr, "Count mismatch: expected %zu, found %zu\n", expected_size,
				actual_count);
	}
	assert(actual_count == expected_size);

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

	// figure 5
	RUN_OVER_APPROX_TEST(figure5, IDR_PARITY2, false, 2, "parity2_paths.txt");
	RUN_OVER_APPROX_TEST(figure5, IDR_SE, false, 2, "se_paths.txt");

	// figure 9
	RUN_OVER_APPROX_TEST(figure9, IDR_PARITY2, false, 7, "parity2_paths.txt");
	RUN_OVER_APPROX_TEST(figure9, IDR_SE, false, 6, "se_paths.txt");

	// figure 10
	RUN_OVER_APPROX_TEST(figure10, IDR_PARITY2, false, 9, "parity2_paths.txt");
	RUN_OVER_APPROX_TEST(figure10, IDR_SE, false, 4, "se_paths.txt");
	RUN_OVER_APPROX_TEST(figure10, IDR_PROJECT, false, 2, "project_paths.txt");

	// figure 11
	RUN_OVER_APPROX_TEST(figure11, IDR_PARITY2, false, 9, "parity2_paths.txt");
	RUN_OVER_APPROX_TEST(figure11, IDR_SE, false, 7, "se_paths.txt");
	RUN_OVER_APPROX_TEST(figure11, IDR_PROJECT, false, 8, "project_paths.txt");
	RUN_OVER_APPROX_TEST(figure11, IDR_EXCLUDE, false, 9, "exclude_paths.txt");

	// loozfon
	RUN_OVER_APPROX_TEST(loozfon, IDR_PARITY, false, 212, "parity1_paths.txt");
	RUN_OVER_APPROX_TEST(loozfon, IDR_PARITY2, false, 211, "parity2_paths.txt");
	RUN_OVER_APPROX_TEST(loozfon, IDR_SE, false, 93, "se_paths.txt");
	RUN_OVER_APPROX_TEST(loozfon, IDR_PROJECT, false, 153, "project_paths.txt");
	RUN_OVER_APPROX_TEST(loozfon, IDR_EXCLUDE, false, 212, "exclude_paths.txt");
	RUN_OVER_APPROX_TEST(loozfon, IDR_ALL, false, 93, "all_paths.txt");

	// faketaobao
	RUN_OVER_APPROX_TEST(faketaobao, IDR_PARITY, false, 151, "parity1_paths.txt");
	RUN_OVER_APPROX_TEST(faketaobao, IDR_PARITY2, false, 151, "parity2_paths.txt");
	RUN_OVER_APPROX_TEST(faketaobao, IDR_SE, false, 64, "se_paths.txt");
	RUN_OVER_APPROX_TEST(faketaobao, IDR_PROJECT, false, 59, "project_paths.txt");
	RUN_OVER_APPROX_TEST(faketaobao, IDR_EXCLUDE, false, 151, "exclude_paths.txt");
	RUN_OVER_APPROX_TEST(faketaobao, IDR_ALL, false, 59, "all_paths.txt");

	// jollyserv
	RUN_OVER_APPROX_TEST(jollyserv, IDR_PARITY, false, 176, "parity1_paths.txt");
	RUN_OVER_APPROX_TEST(jollyserv, IDR_PARITY2, false, 175, "parity2_paths.txt");
	RUN_OVER_APPROX_TEST(jollyserv, IDR_SE, false, 164, "se_paths.txt");
	RUN_OVER_APPROX_TEST(jollyserv, IDR_PROJECT, false, 174, "project_paths.txt");
	RUN_OVER_APPROX_TEST(jollyserv, IDR_EXCLUDE, false, 176, "exclude_paths.txt");
	RUN_OVER_APPROX_TEST(jollyserv, IDR_ALL, false, 164, "all_paths.txt");

	// zertsecurity
	RUN_OVER_APPROX_TEST(zertsecurity, IDR_PARITY, false, 1081, "parity1_paths.txt");
	RUN_OVER_APPROX_TEST(zertsecurity, IDR_PARITY2, false, 1045,
						 "parity2_paths.txt");
	RUN_OVER_APPROX_TEST(zertsecurity, IDR_SE, false, 883, "se_paths.txt");
	RUN_OVER_APPROX_TEST(zertsecurity, IDR_PROJECT, false, 1081,
						 "project_paths.txt");
	RUN_OVER_APPROX_TEST(zertsecurity, IDR_EXCLUDE, false, 1080,
						 "exclude_paths.txt");
	RUN_OVER_APPROX_TEST(zertsecurity, IDR_ALL, false, 883, "all_paths.txt");

	// fakebanker
	RUN_OVER_APPROX_TEST(fakebanker, IDR_PARITY, false, 590, "parity1_paths.txt");
	RUN_OVER_APPROX_TEST(fakebanker, IDR_PARITY2, false, 590, "parity2_paths.txt");
	RUN_OVER_APPROX_TEST(fakebanker, IDR_SE, false, 279, "se_paths.txt");
	RUN_OVER_APPROX_TEST(fakebanker, IDR_PROJECT, false, 354, "project_paths.txt");
	RUN_OVER_APPROX_TEST(fakebanker, IDR_EXCLUDE, false, 590, "exclude_paths.txt");
	RUN_OVER_APPROX_TEST(fakebanker, IDR_ALL, false, 272, "all_paths.txt");

	// uranai
	RUN_OVER_APPROX_TEST(uranai, IDR_PARITY, false, 143, "parity1_paths.txt");
	RUN_OVER_APPROX_TEST(uranai, IDR_PARITY2, false, 143, "parity2_paths.txt");
	RUN_OVER_APPROX_TEST(uranai, IDR_SE, false, 143, "se_paths.txt");
	RUN_OVER_APPROX_TEST(uranai, IDR_PROJECT, false, 143, "project_paths.txt");
	RUN_OVER_APPROX_TEST(uranai, IDR_EXCLUDE, false, 143, "exclude_paths.txt");
	RUN_OVER_APPROX_TEST(uranai, IDR_ALL, false, 143, "all_paths.txt");

	// Value-flow analysis
	RUN_OVER_APPROX_TEST(xz, IDR_PARITY2, /*valueflow=*/true, 211, "paths.txt");
	RUN_OVER_APPROX_TEST(nab, IDR_PARITY2, /*valueflow=*/true, 1788, "paths.txt");
	RUN_OVER_APPROX_TEST(leela, IDR_PARITY2, /*valueflow=*/true, 392, "paths.txt");

	LAGraph_Finalize(msg);
	return 0;
}
