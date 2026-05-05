#include "output.h"

#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include "internal/grb_utils.h"

GrB_Info write_matrix_pairs(FILE *out, GrB_Matrix matrix, GrB_Index nvals) {
	GrB_Info info = GrB_SUCCESS;

	GrB_Index *rows = NULL;
	GrB_Index *cols = NULL;
	bool *vals = NULL;

	if (nvals == 0) {
		return GrB_SUCCESS;
	}

	rows = malloc(nvals * sizeof(GrB_Index));
	cols = malloc(nvals * sizeof(GrB_Index));
	vals = malloc(nvals * sizeof(bool));
	if (!rows || !cols || !vals) {
		info = GrB_OUT_OF_MEMORY;
		goto cleanup;
	}

	GRB_TRY(GrB_Matrix_extractTuples_BOOL(rows, cols, vals, &nvals, matrix));
	for (GrB_Index k = 0; k < nvals; k++) {
		fprintf(out, "\t%lu %lu\n", (unsigned long)rows[k], (unsigned long)cols[k]);
	}

cleanup:
	free(rows);
	free(cols);
	free(vals);
	return info;
}

void resolve_output_path(const Args *args, char *output_file, size_t size) {
	const char *base = strrchr(args->graph_file_path, '/');
	base = base ? base + 1 : args->graph_file_path;

	char benchmark_name[256];
	strncpy(benchmark_name, base, sizeof(benchmark_name) - 1);
	benchmark_name[sizeof(benchmark_name) - 1] = '\0';
	char *dot = strrchr(benchmark_name, '.');
	if (dot) {
		*dot = '\0';
	}

	const char *dir = args->output_path ? args->output_path : DEFAULT_OUTPUT_DIR;
	mkdir(dir, 0755);
	snprintf(output_file, size, "%s/%s.out", dir, benchmark_name);
}
