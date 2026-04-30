#include "utils.h"

#include "approximation/approximation.h"
#include "approximation/on_demand.h"
#include "parse_utils.h"

static char msg[LAGRAPH_MSG_LEN];

#define RUN_ON_DEMAND_TEST(folder, parity_d, expected_count, expected_file)         \
	do {                                                                            \
		char g_path[1024], e_path[1024];                                            \
		snprintf(g_path, 1024, "%s/" #folder "/graph.g", TEST_DATA_DIR);            \
		snprintf(e_path, 1024, "%s/" #folder "/" expected_file, TEST_DATA_DIR);     \
		printf("Running test: on_demand parity_d=%d on " #folder "...\n",           \
			   parity_d);                                                           \
		run_test_logic(g_path, parity_d, expected_count, e_path);                   \
	} while (0)

static void run_test_logic(const char *graph_path, bool parity_d,
						   size_t expected_size, const char *expected_path) {
	MRGraph graph = {0};
	GrB_Info info = parse_graph(graph_path, &graph);
	assert(info == GrB_SUCCESS);

	GrB_Matrix under_approx = NULL;
	info = get_under_approx(&graph, &under_approx);
	assert(info == GrB_SUCCESS);

	MRGrammarType grammar_type;
	if (parity_d) {
		grammar_type = PARITY;
	} else {
		grammar_type = ALL;
	}

	GrB_Matrix over_approx = NULL;
	info = get_over_approx(&graph, grammar_type, NULL, &over_approx, true);
	assert(info == GrB_SUCCESS);

	GrB_Matrix result = NULL;
	info = get_on_demand(&graph, under_approx, over_approx, parity_d, &result, true,
						 msg);
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
	GrB_Matrix_free(&under_approx);
	GrB_Matrix_free(&over_approx);
	mr_graph_free(&graph);
}

int main(void) {
	LAGraph_Init(msg);

	// figure 5
	RUN_ON_DEMAND_TEST(figure5, false, 2, "onDemand_paths.txt");

	// figure 9
	RUN_ON_DEMAND_TEST(figure9, false, 6, "onDemand_paths.txt");

	// figure 10
	RUN_ON_DEMAND_TEST(figure10, false, 2, "onDemand_paths.txt");

	// figure 11
	RUN_ON_DEMAND_TEST(figure11, false, 6, "onDemand_paths.txt");

	// loozfon
	RUN_ON_DEMAND_TEST(loozfon, true, 93, "parityD_paths.txt");
	RUN_ON_DEMAND_TEST(loozfon, false, 93, "onDemand_paths.txt");

	// faketaobao
	RUN_ON_DEMAND_TEST(faketaobao, true, 61, "parityD_paths.txt");
	RUN_ON_DEMAND_TEST(faketaobao, false, 59, "onDemand_paths.txt");

	// jollyserv
	RUN_ON_DEMAND_TEST(jollyserv, true, 164, "parityD_paths.txt");
	RUN_ON_DEMAND_TEST(jollyserv, false, 164, "onDemand_paths.txt");

	// zertsecurity
	RUN_ON_DEMAND_TEST(zertsecurity, true, 808, "parityD_paths.txt");
	RUN_ON_DEMAND_TEST(zertsecurity, false, 794, "onDemand_paths.txt");

	// fakebanker
	RUN_ON_DEMAND_TEST(fakebanker, true, 254, "parityD_paths.txt");
	RUN_ON_DEMAND_TEST(fakebanker, false, 251, "onDemand_paths.txt");

	// uranai
	RUN_ON_DEMAND_TEST(uranai, true, 143, "parityD_paths.txt");
	RUN_ON_DEMAND_TEST(uranai, false, 143, "onDemand_paths.txt");

	LAGraph_Finalize(msg);
	return 0;
}
