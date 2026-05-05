#include "parse_utils.h"

#include "convert_graph.h"
#include "parser.h"
#include "symbol_list.h"

GrB_Info parse_graph(const char *filename, IdrGraph *out) {
	FILE *graph_file = fopen(filename, "r");
	if (graph_file == NULL) {
		fprintf(stderr, "\x1B[31m[ERROR]\033[0m Could not open graph file: %s\n",
				filename);
		return GrB_INVALID_VALUE;
	}

	SymbolList symbol_list = symbol_list_create();
	Graph graph = process_graph(graph_file, &symbol_list);
	fclose(graph_file);

	GraphMatrices gm = get_grb_matrices_from_graph(graph, &symbol_list);
	free(graph.edges);

	GrB_Info info = get_idr_graph(&gm, &symbol_list, graph.node_count,
								  &DefaultTerminalFormat, out);

	free(gm.matrix_symbols);
	symbol_list_free(&symbol_list);

	return info;
}
