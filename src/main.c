#include "approximation/approximation.h"
#include "approximation/on_demand.h"
#include "graph/parse_utils.h"
#include "graph/valueflow_extensions.h"
#include "io/cli.h"
#include "io/output.h"

int main(int argc, char *argv[]) {
	char msg[LAGRAPH_MSG_LEN];

	Args args = {0};
	if (!parse_args(argc, argv, &args)) {
		print_usage(argv[0]);
		return EXIT_FAILURE;
	}

	char output_file[512];
	resolve_output_path(&args, output_file, sizeof(output_file));

	LAGraph_Init(msg);

	MRGraph refined_graph = {0};
	GrB_Info info = parse_graph(args.graph_file_path, &refined_graph);
	if (info != GrB_SUCCESS) {
		fprintf(stderr, "Error: Failed to parse graph '%s': %d\n",
				args.graph_file_path, info);
		LAGraph_Finalize(msg);
		return EXIT_FAILURE;
	}

	MRGraph filtered_graph = {0};
	if (args.valueflow) {
		remove_valueflow_unreachable(&refined_graph, &filtered_graph, msg);
		mr_graph_free(&refined_graph);
	} else {
		filtered_graph = refined_graph;
	}

	// Under approximation
	GrB_Matrix under_result = NULL;
	info = get_under_approx(&filtered_graph, args.valueflow, &under_result);
	if (info != GrB_SUCCESS) {
		fprintf(stderr, "Error in under-approximation: %d\n", info);
		mr_graph_free(&filtered_graph);
		LAGraph_Finalize(msg);
		return EXIT_FAILURE;
	}
	GrB_Matrix_wait(under_result, GrB_MATERIALIZE);

	GrB_Index under_nvals;
	GrB_Matrix_nvals(&under_nvals, under_result);

	// Over approximation
	GrB_Matrix over_result = NULL;
	info = get_over_approx(&filtered_graph, args.grammar_type, under_result,
						   &over_result, args.valueflow, true);
	if (info != GrB_SUCCESS) {
		fprintf(stderr, "Error in over-approximation: %d\n", info);
		GrB_Matrix_free(&under_result);
		mr_graph_free(&filtered_graph);
		LAGraph_Finalize(msg);
		return EXIT_FAILURE;
	}
	GrB_Matrix_wait(over_result, GrB_MATERIALIZE);

	GrB_Index over_nvals;
	GrB_Matrix_nvals(&over_nvals, over_result);

	// Write results
	FILE *out = fopen(output_file, "w");
	if (out == NULL) {
		fprintf(stderr, "Error: Cannot open output file '%s'\n", output_file);
		GrB_Matrix_free(&under_result);
		GrB_Matrix_free(&over_result);
		mr_graph_free(&filtered_graph);
		LAGraph_Finalize(msg);
		return EXIT_FAILURE;
	}

	fprintf(out, "Under approximation paths: %lu\n", (unsigned long)under_nvals);
	fprintf(out, "\nOver approximation paths (%s): %lu\n", args.grammar_str,
			(unsigned long)over_nvals);

	if (!args.on_demand && !args.parity_d) {
		if (!args.quiet) {
			write_matrix_pairs(out, over_result, over_nvals);
		}
		fclose(out);
		GrB_Matrix_free(&under_result);
		GrB_Matrix_free(&over_result);
		mr_graph_free(&filtered_graph);
		LAGraph_Finalize(msg);
		printf("Analysis completed. Results written to %s\n", output_file);
		return EXIT_SUCCESS;
	}

	// On-demand refinement
	GrB_Matrix on_demand_result = NULL;
	info = get_on_demand(&filtered_graph, under_result, over_result, args.parity_d,
						 &on_demand_result, args.valueflow, true, msg);
	if (info != GrB_SUCCESS) {
		fprintf(stderr, "Error in on-demand refinement: %d\n", info);
		fclose(out);
		GrB_Matrix_free(&under_result);
		GrB_Matrix_free(&over_result);
		mr_graph_free(&filtered_graph);
		LAGraph_Finalize(msg);
		return EXIT_FAILURE;
	}
	GrB_Matrix_wait(on_demand_result, GrB_MATERIALIZE);

	GrB_Index od_nvals;
	GrB_Matrix_nvals(&od_nvals, on_demand_result);

	fprintf(out, "\nOn-Demand paths (%s): %lu\n",
			args.parity_d ? "parityD" : "on-demand", (unsigned long)od_nvals);
	if (!args.quiet) {
		write_matrix_pairs(out, on_demand_result, od_nvals);
	}

	fclose(out);
	GrB_Matrix_free(&under_result);
	GrB_Matrix_free(&over_result);
	GrB_Matrix_free(&on_demand_result);
	mr_graph_free(&filtered_graph);
	LAGraph_Finalize(msg);

	printf("On-demand refinement completed. Results written to %s\n", output_file);
	return EXIT_SUCCESS;
}
