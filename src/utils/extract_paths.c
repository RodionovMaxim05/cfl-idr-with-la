#include "extract_paths.h"

#include <LAGraphX.h>
#include <stdlib.h>

#include "internal/grb_utils.h"

GrB_Info extract_non_trivial_paths(GrB_Matrix *out, GrB_Matrix paths) {
	GrB_Info info = GrB_SUCCESS;

	GrB_Index *rows = NULL;
	GrB_Index *cols = NULL;
	AllPathsElem *vals = NULL;

	GrB_Index nvals = 0;
	GrB_Matrix_nvals(&nvals, paths);
	if (nvals == 0) {
		goto cleanup;
	}

	rows = malloc(nvals * sizeof(GrB_Index));
	cols = malloc(nvals * sizeof(GrB_Index));
	vals = malloc(nvals * sizeof(AllPathsElem));
	if (!rows || !cols || !vals) {
		info = GrB_OUT_OF_MEMORY;
		goto cleanup;
	}

	GRB_TRY(GrB_Matrix_extractTuples_UDT(rows, cols, vals, &nvals, paths));

	for (GrB_Index k = 0; k < nvals; k++) {
		if (rows[k] != cols[k]) {
			GRB_TRY(GrB_Matrix_setElement_BOOL(*out, true, rows[k], cols[k]));
		}
	}

cleanup:
	free(rows);
	free(cols);
	free(vals);
	return info;
}
