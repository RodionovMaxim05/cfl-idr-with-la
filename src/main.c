#include "cfl_idr.h"

#include "graph/parse_utils.h"
#include "io/cli.h"
#include "io/output.h"

int main(int argc, char *argv[]) {
	char msg[LAGRAPH_MSG_LEN];
	int exit_code = EXIT_SUCCESS;

	GrB_Matrix under_result = NULL;
	GrB_Matrix over_result = NULL;
	GrB_Matrix on_demand_result = NULL;
	IdrGraph parsed_graph = {0};
	IdrGraph working_graph = {0};
	bool working_graph_owned = false;
	FILE *out = NULL;

	Args args = {0};
	if (!parse_args(argc, argv, &args)) {
		print_usage(argv[0]);
		return EXIT_FAILURE;
	}

	char output_file[512];
	resolve_output_path(&args, output_file, sizeof(output_file));

	LAGraph_Init(msg);

	GrB_Info info = parse_graph(args.graph_file_path, &parsed_graph);
	if (info != GrB_SUCCESS) {
		fprintf(stderr, "Error: Failed to parse graph '%s': %d\n",
				args.graph_file_path, info);
		exit_code = EXIT_FAILURE;
		goto cleanup;
	}

	if (args.valueflow) {
		info = idr_remove_valueflow_unreachable(&working_graph, &parsed_graph);
		if (info != GrB_SUCCESS) {
			fprintf(stderr, "Error: valueflow filtering failed: %d\n", info);
			exit_code = EXIT_FAILURE;
			goto cleanup;
		}
		idr_graph_free(&parsed_graph);
		working_graph_owned = true;
	} else {
		working_graph = parsed_graph;
	}

	// Under approximation
	info = idr_get_under_approx(&under_result, &working_graph, args.valueflow);
	if (info != GrB_SUCCESS) {
		fprintf(stderr, "Error: under-approximation failed: %d\n", info);
		exit_code = EXIT_FAILURE;
		goto cleanup;
	}
	GrB_Matrix_wait(under_result, GrB_MATERIALIZE);

	GrB_Index under_nvals = 0;
	GrB_Matrix_nvals(&under_nvals, under_result);

	// Over approximation
	info = idr_get_over_approx(&over_result, &working_graph, args.grammar_type,
							   under_result, args.valueflow, true);
	if (info != GrB_SUCCESS) {
		fprintf(stderr, "Error: over-approximation failed: %d\n", info);
		exit_code = EXIT_FAILURE;
		goto cleanup;
	}
	GrB_Matrix_wait(over_result, GrB_MATERIALIZE);

	GrB_Index over_nvals = 0;
	GrB_Matrix_nvals(&over_nvals, over_result);

	// Write results
	out = fopen(output_file, "w");
	if (!out) {
		fprintf(stderr, "Error: cannot open output file '%s'\n", output_file);
		exit_code = EXIT_FAILURE;
		goto cleanup;
	}

	fprintf(out, "Under approximation paths: %lu\n", (unsigned long)under_nvals);
	fprintf(out, "\nOver approximation paths (%s): %lu\n", args.grammar_str,
			(unsigned long)over_nvals);

	if (!args.on_demand && !args.parity_d) {
		if (!args.quiet) {
			write_matrix_pairs(out, over_result, over_nvals);
		}
		printf("Analysis completed. Results written to %s\n", output_file);
		goto cleanup;
	}

	// On-demand refinement
	info = idr_get_on_demand(&on_demand_result, &working_graph, under_result,
							 over_result, args.parity_d, args.valueflow, true);
	if (info != GrB_SUCCESS) {
		fprintf(stderr, "Error: on-demand refinement failed: %d\n", info);
		exit_code = EXIT_FAILURE;
		goto cleanup;
	}
	GrB_Matrix_wait(on_demand_result, GrB_MATERIALIZE);

	GrB_Index od_nvals = 0;
	GrB_Matrix_nvals(&od_nvals, on_demand_result);

	fprintf(out, "\nOn-Demand paths (%s): %lu\n",
			args.parity_d ? "parityD" : "on-demand", (unsigned long)od_nvals);
	if (!args.quiet) {
		write_matrix_pairs(out, on_demand_result, od_nvals);
	}

	printf("On-demand refinement completed. Results written to %s\n", output_file);

cleanup:
	if (out) {
		fclose(out);
	}
	GrB_Matrix_free(&under_result);
	GrB_Matrix_free(&over_result);
	GrB_Matrix_free(&on_demand_result);
	if (working_graph_owned) {
		idr_graph_free(&working_graph);
	} else {
		idr_graph_free(&parsed_graph);
	}
	LAGraph_Finalize(msg);
	return exit_code;
}
