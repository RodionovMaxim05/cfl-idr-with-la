#include "utils.h"

static char msg[LAGRAPH_MSG_LEN];

#define RUN_ON_DEMAND_TEST(folder, parity_d, valueflow, expected_count,             \
						   expected_file)                                           \
	do {                                                                            \
		char g_path[1024], e_path[1024];                                            \
		snprintf(g_path, 1024, "%s/" #folder "/graph.g", TEST_DATA_DIR);            \
		snprintf(e_path, 1024, "%s/" #folder "/" expected_file, TEST_DATA_DIR);     \
		char *grammar = parity_d ? "PARITYD" : "ON-DEMAND";                         \
		printf("Running test: %s on " #folder "%s...\n", grammar,                   \
			   valueflow ? " (valueflow)" : "");                                    \
		fflush(stdout);                                                             \
		run_test_logic(g_path, parity_d, valueflow, expected_count, e_path);        \
	} while (0)

static void run_test_logic(const char *graph_path, bool parity_d, bool valueflow,
						   size_t expected_size, const char *expected_path) {
	IdrGraph graph = {0};
	GrB_Info info = parse_graph(graph_path, &graph);
	assert(info == GrB_SUCCESS);

	IdrGraph new_graph = {0};
	if (valueflow) {
		info = idr_remove_valueflow_unreachable(&graph, &new_graph);
		assert(info == GrB_SUCCESS);
		idr_graph_free(&graph);
	} else {
		new_graph = graph;
	}

	GrB_Matrix under_approx = NULL;
	info = idr_get_under_approx(&new_graph, valueflow, &under_approx);
	assert(info == GrB_SUCCESS);

	IdrGrammarType grammar_type;
	if (parity_d) {
		grammar_type = IDR_PARITY;
	} else {
		grammar_type = IDR_ALL;
	}

	GrB_Matrix over_approx = NULL;
	info = idr_get_over_approx(&new_graph, grammar_type, NULL, &over_approx,
							   valueflow, true);
	assert(info == GrB_SUCCESS);

	GrB_Matrix result = NULL;
	info = idr_get_on_demand(&new_graph, under_approx, over_approx, parity_d,
							 &result, valueflow, true);
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
	idr_graph_free(&new_graph);
}

int main(void) {
	LAGraph_Init(msg);

	// figure 5
	RUN_ON_DEMAND_TEST(figure5, /*parityD=*/false, /*valueflow=*/false, 2,
					   "onDemand_paths.txt");

	// figure 9
	RUN_ON_DEMAND_TEST(figure9, /*parityD=*/false, /*valueflow=*/false, 6,
					   "onDemand_paths.txt");

	// figure 10
	RUN_ON_DEMAND_TEST(figure10, /*parityD=*/false, /*valueflow=*/false, 2,
					   "onDemand_paths.txt");

	// figure 11
	RUN_ON_DEMAND_TEST(figure11, /*parityD=*/false, /*valueflow=*/false, 6,
					   "onDemand_paths.txt");

	// loozfon
	RUN_ON_DEMAND_TEST(loozfon, /*parityD=*/true, /*valueflow=*/false, 93,
					   "parityD_paths.txt");
	RUN_ON_DEMAND_TEST(loozfon, /*parityD=*/false, /*valueflow=*/false, 93,
					   "onDemand_paths.txt");

	// faketaobao
	RUN_ON_DEMAND_TEST(faketaobao, /*parityD=*/true, /*valueflow=*/false, 61,
					   "parityD_paths.txt");
	RUN_ON_DEMAND_TEST(faketaobao, /*parityD=*/false, /*valueflow=*/false, 59,
					   "onDemand_paths.txt");

	// jollyserv
	RUN_ON_DEMAND_TEST(jollyserv, /*parityD=*/true, /*valueflow=*/false, 164,
					   "parityD_paths.txt");
	RUN_ON_DEMAND_TEST(jollyserv, /*parityD=*/false, /*valueflow=*/false, 164,
					   "onDemand_paths.txt");

	// zertsecurity
	RUN_ON_DEMAND_TEST(zertsecurity, /*parityD=*/true, /*valueflow=*/false, 808,
					   "parityD_paths.txt");
	RUN_ON_DEMAND_TEST(zertsecurity, /*parityD=*/false, /*valueflow=*/false, 794,
					   "onDemand_paths.txt");

	// fakebanker
	RUN_ON_DEMAND_TEST(fakebanker, /*parityD=*/true, /*valueflow=*/false, 254,
					   "parityD_paths.txt");
	RUN_ON_DEMAND_TEST(fakebanker, /*parityD=*/false, /*valueflow=*/false, 251,
					   "onDemand_paths.txt");

	// uranai
	RUN_ON_DEMAND_TEST(uranai, /*parityD=*/true, /*valueflow=*/false, 143,
					   "parityD_paths.txt");
	RUN_ON_DEMAND_TEST(uranai, /*parityD=*/false, /*valueflow=*/false, 143,
					   "onDemand_paths.txt");

	// Value-flow analysis
	RUN_ON_DEMAND_TEST(xz, /*parityD=*/false, /*valueflow=*/true, 211, "paths.txt");
	RUN_ON_DEMAND_TEST(nab, /*parityD=*/false, /*valueflow=*/true, 1788,
					   "paths.txt");
	RUN_ON_DEMAND_TEST(leela, /*parityD=*/false, /*valueflow=*/true, 392,
					   "paths.txt");

	LAGraph_Finalize(msg);
	return 0;
}
