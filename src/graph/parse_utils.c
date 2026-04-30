#include "parse_utils.h"

#include "parser.h"
#include "symbol_list.h"
#include "terminal/convert_graph.h"

GrB_Info parse_graph(const char *filename, MRGraph *out) {
	FILE *graph_file = fopen(filename, "r");
	if (graph_file == NULL) {
		fprintf(stderr, "\x1B[31m[ERROR]\033[0m Could not open graph file: %s\n",
				filename);
		return GrB_INVALID_VALUE;
	}

	SymbolList symbol_list = symbol_list_create();
	Graph graph = process_graph(graph_file, &symbol_list);
	fclose(graph_file);

	GrB_Matrix *matrices = get_grb_matrices_from_graph(graph, &symbol_list);
	free(graph.edges);

	GrB_Info info = build_mr_graph(matrices, &symbol_list, graph.node_count,
								   &DefaultTerminalFormat, out);

	free((void *)matrices);
	symbol_list_free(&symbol_list);

	return info;
}
