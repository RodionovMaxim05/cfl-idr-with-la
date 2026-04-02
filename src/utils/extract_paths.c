#include "extract_paths.h"

#include "LAGraphX.h"

void extractNonTrivialPaths(GrB_Matrix paths, const MRGraph *graph,
							GrB_Matrix *result) {
	GrB_Index nvals = 0;
	GrB_Matrix_nvals(&nvals, paths);

	GrB_Index *rows = malloc(nvals * sizeof(GrB_Index));
	GrB_Index *cols = malloc(nvals * sizeof(GrB_Index));
	AllPathsElem *vals = malloc(nvals * sizeof(AllPathsElem));
	GrB_Matrix_extractTuples_UDT(rows, cols, vals, &nvals, paths);

	for (GrB_Index k = 0; k < nvals; k++) {
		if (rows[k] != cols[k]) {
			GrB_Matrix_setElement_BOOL(*result, true, rows[k], cols[k]);
		}
	}

	free(rows);
	free(cols);
	free(vals);
}
