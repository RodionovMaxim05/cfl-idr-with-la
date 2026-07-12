#include "parse_utils.h"

#include "cfl_idr_graph_builder.h"
#include "parser.h"
#include "symbol_list.h"

GrB_Info parse_graph(IdrGraph *out, const char *filename) {
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

	GrB_Info info = get_idr_graph(out, &gm, &symbol_list, graph.node_count,
								  &DefaultTerminalFormat);

	if (gm.matrices != NULL) {
		for (size_t i = 0; i < gm.count; i++) {
			if (gm.matrices[i] != NULL) {
				GrB_Matrix_free(&gm.matrices[i]);
			}
		}
		free((void *)gm.matrices);
	}
	free(gm.matrix_symbols);
	symbol_list_free(&symbol_list);

	return info;
}
