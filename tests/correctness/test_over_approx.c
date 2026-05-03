#include "utils.h"

#include "approximation/approximation.h"
#include "graph/parse_utils.h"

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

static void run_test_logic(const char *graph_path, MRGrammarType grammar_type,
						   bool valueflow, size_t expected_size,
						   const char *expected_path) {
	MRGraph graph = {0};
	GrB_Info info = parse_graph(graph_path, &graph);
	assert(info == GrB_SUCCESS);

	GrB_Matrix result = NULL;
	info = get_over_approx(&graph, grammar_type, NULL, &result, valueflow, true);
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
	mr_graph_free(&graph);
}

int main(void) {
	LAGraph_Init(msg);

	// figure 5
	RUN_OVER_APPROX_TEST(figure5, PARITY2, false, 2, "parity2_paths.txt");
	RUN_OVER_APPROX_TEST(figure5, SE, false, 2, "se_paths.txt");

	// figure 9
	RUN_OVER_APPROX_TEST(figure9, PARITY2, false, 7, "parity2_paths.txt");
	RUN_OVER_APPROX_TEST(figure9, SE, false, 6, "se_paths.txt");

	// figure 10
	RUN_OVER_APPROX_TEST(figure10, PARITY2, false, 9, "parity2_paths.txt");
	RUN_OVER_APPROX_TEST(figure10, SE, false, 4, "se_paths.txt");
	RUN_OVER_APPROX_TEST(figure10, PROJECT, false, 2, "project_paths.txt");

	// figure 11
	RUN_OVER_APPROX_TEST(figure11, PARITY2, false, 9, "parity2_paths.txt");
	RUN_OVER_APPROX_TEST(figure11, SE, false, 7, "se_paths.txt");
	RUN_OVER_APPROX_TEST(figure11, PROJECT, false, 8, "project_paths.txt");
	RUN_OVER_APPROX_TEST(figure11, EXCLUDE, false, 9, "exclude_paths.txt");

	// loozfon
	RUN_OVER_APPROX_TEST(loozfon, PARITY, false, 212, "parity1_paths.txt");
	RUN_OVER_APPROX_TEST(loozfon, PARITY2, false, 211, "parity2_paths.txt");
	RUN_OVER_APPROX_TEST(loozfon, SE, false, 93, "se_paths.txt");
	RUN_OVER_APPROX_TEST(loozfon, PROJECT, false, 153, "project_paths.txt");
	RUN_OVER_APPROX_TEST(loozfon, EXCLUDE, false, 212, "exclude_paths.txt");
	RUN_OVER_APPROX_TEST(loozfon, ALL, false, 93, "all_paths.txt");

	// faketaobao
	RUN_OVER_APPROX_TEST(faketaobao, PARITY, false, 151, "parity1_paths.txt");
	RUN_OVER_APPROX_TEST(faketaobao, PARITY2, false, 151, "parity2_paths.txt");
	RUN_OVER_APPROX_TEST(faketaobao, SE, false, 64, "se_paths.txt");
	RUN_OVER_APPROX_TEST(faketaobao, PROJECT, false, 59, "project_paths.txt");
	RUN_OVER_APPROX_TEST(faketaobao, EXCLUDE, false, 151, "exclude_paths.txt");
	RUN_OVER_APPROX_TEST(faketaobao, ALL, false, 59, "all_paths.txt");

	// jollyserv
	RUN_OVER_APPROX_TEST(jollyserv, PARITY, false, 176, "parity1_paths.txt");
	RUN_OVER_APPROX_TEST(jollyserv, PARITY2, false, 175, "parity2_paths.txt");
	RUN_OVER_APPROX_TEST(jollyserv, SE, false, 164, "se_paths.txt");
	RUN_OVER_APPROX_TEST(jollyserv, PROJECT, false, 174, "project_paths.txt");
	RUN_OVER_APPROX_TEST(jollyserv, EXCLUDE, false, 176, "exclude_paths.txt");
	RUN_OVER_APPROX_TEST(jollyserv, ALL, false, 164, "all_paths.txt");

	// zertsecurity
	RUN_OVER_APPROX_TEST(zertsecurity, PARITY, false, 1081, "parity1_paths.txt");
	RUN_OVER_APPROX_TEST(zertsecurity, PARITY2, false, 1045, "parity2_paths.txt");
	RUN_OVER_APPROX_TEST(zertsecurity, SE, false, 883, "se_paths.txt");
	RUN_OVER_APPROX_TEST(zertsecurity, PROJECT, false, 1081, "project_paths.txt");
	RUN_OVER_APPROX_TEST(zertsecurity, EXCLUDE, false, 1080, "exclude_paths.txt");
	RUN_OVER_APPROX_TEST(zertsecurity, ALL, false, 883, "all_paths.txt");

	// fakebanker
	RUN_OVER_APPROX_TEST(fakebanker, PARITY, false, 590, "parity1_paths.txt");
	RUN_OVER_APPROX_TEST(fakebanker, PARITY2, false, 590, "parity2_paths.txt");
	RUN_OVER_APPROX_TEST(fakebanker, SE, false, 279, "se_paths.txt");
	RUN_OVER_APPROX_TEST(fakebanker, PROJECT, false, 354, "project_paths.txt");
	RUN_OVER_APPROX_TEST(fakebanker, EXCLUDE, false, 590, "exclude_paths.txt");
	RUN_OVER_APPROX_TEST(fakebanker, ALL, false, 272, "all_paths.txt");

	// uranai
	RUN_OVER_APPROX_TEST(uranai, PARITY, false, 143, "parity1_paths.txt");
	RUN_OVER_APPROX_TEST(uranai, PARITY2, false, 143, "parity2_paths.txt");
	RUN_OVER_APPROX_TEST(uranai, SE, false, 143, "se_paths.txt");
	RUN_OVER_APPROX_TEST(uranai, PROJECT, false, 143, "project_paths.txt");
	RUN_OVER_APPROX_TEST(uranai, EXCLUDE, false, 143, "exclude_paths.txt");
	RUN_OVER_APPROX_TEST(uranai, ALL, false, 143, "all_paths.txt");

	// Value-flow analysis
	RUN_OVER_APPROX_TEST(xz, PARITY2, /*valueflow=*/true, 211, "paths.txt");
	RUN_OVER_APPROX_TEST(nab, PARITY2, /*valueflow=*/true, 1788, "paths.txt");
	RUN_OVER_APPROX_TEST(leela, PARITY2, /*valueflow=*/true, 392, "paths.txt");

	LAGraph_Finalize(msg);
	return 0;
}
