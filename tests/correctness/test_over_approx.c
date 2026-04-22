#include "utils.h"

#include "approximation/approximation.h"
#include "parse_utils.h"

static char msg[LAGRAPH_MSG_LEN];

#define RUN_OVER_APPROX_TEST(folder, grammar_enum, expected_count, expected_file)   \
	do {                                                                            \
		char g_path[1024], e_path[1024];                                            \
		snprintf(g_path, 1024, "%s/" #folder "/graph.g", TEST_DATA_DIR);            \
		snprintf(e_path, 1024, "%s/" #folder "/" expected_file, TEST_DATA_DIR);     \
		printf("Running test: " #grammar_enum " on " #folder "...\n");              \
		run_test_logic(g_path, grammar_enum, expected_count, e_path);               \
	} while (0)

static void run_test_logic(const char *graph_path, MRGrammarType grammar_type,
						   size_t expected_size, const char *expected_path) {
	MRGraph graph = {0};
	GrB_Info info = parse_graph(graph_path, &graph);
	assert(info == GrB_SUCCESS);

	GrB_Matrix result = NULL;
	info = get_over_approx(&graph, grammar_type, NULL, &result, true);
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
	RUN_OVER_APPROX_TEST(figure5, PARITY2, 2, "parity2_paths.txt");
	RUN_OVER_APPROX_TEST(figure5, SE, 2, "se_paths.txt");

	// figure 9
	RUN_OVER_APPROX_TEST(figure9, PARITY2, 7, "parity2_paths.txt");
	RUN_OVER_APPROX_TEST(figure9, SE, 6, "se_paths.txt");

	// figure 10
	RUN_OVER_APPROX_TEST(figure10, PARITY2, 9, "parity2_paths.txt");
	RUN_OVER_APPROX_TEST(figure10, SE, 4, "se_paths.txt");
	// RUN_OVER_APPROX_TEST(figure10, PROJECT, 2, "project_paths.txt");

	// figure 11
	RUN_OVER_APPROX_TEST(figure11, PARITY2, 9, "parity2_paths.txt");
	RUN_OVER_APPROX_TEST(figure11, SE, 7, "se_paths.txt");
	// RUN_OVER_APPROX_TEST(figure11, PROJECT, 8, "project_paths.txt");
	// RUN_OVER_APPROX_TEST(figure11, EXCLUDE, 9, "exclude_paths.txt");

	// loozfon
	RUN_OVER_APPROX_TEST(loozfon, PARITY, 212, "parity1_paths.txt");
	RUN_OVER_APPROX_TEST(loozfon, PARITY2, 211, "parity2_paths.txt");
	RUN_OVER_APPROX_TEST(loozfon, SE, 93, "se_paths.txt");

	// faketaobao
	RUN_OVER_APPROX_TEST(faketaobao, PARITY, 151, "parity1_paths.txt");
	RUN_OVER_APPROX_TEST(faketaobao, PARITY2, 151, "parity2_paths.txt");
	RUN_OVER_APPROX_TEST(faketaobao, SE, 64, "se_paths.txt");

	// jollyserv
	RUN_OVER_APPROX_TEST(jollyserv, PARITY, 176, "parity1_paths.txt");
	RUN_OVER_APPROX_TEST(jollyserv, PARITY2, 175, "parity2_paths.txt");
	RUN_OVER_APPROX_TEST(jollyserv, SE, 164, "se_paths.txt");

	// zertsecurity
	RUN_OVER_APPROX_TEST(zertsecurity, PARITY, 1081, "parity1_paths.txt");
	RUN_OVER_APPROX_TEST(zertsecurity, PARITY2, 1045, "parity2_paths.txt");
	RUN_OVER_APPROX_TEST(zertsecurity, SE, 883, "se_paths.txt");

	// fakebanker
	RUN_OVER_APPROX_TEST(fakebanker, PARITY, 590, "parity1_paths.txt");
	RUN_OVER_APPROX_TEST(fakebanker, PARITY2, 590, "parity2_paths.txt");
	RUN_OVER_APPROX_TEST(fakebanker, SE, 279, "se_paths.txt");

	// uranai
	RUN_OVER_APPROX_TEST(uranai, PARITY, 143, "parity1_paths.txt");
	RUN_OVER_APPROX_TEST(uranai, PARITY2, 143, "parity2_paths.txt");
	RUN_OVER_APPROX_TEST(uranai, SE, 143, "se_paths.txt");

	LAGraph_Finalize(msg);
	return 0;
}
