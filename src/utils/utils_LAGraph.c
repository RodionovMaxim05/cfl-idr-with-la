#include "utils_LAGraph.h"

#include "LAGraphX.h"

// Cleaning of internal elements before free matrix
void free_AllPaths_matrix(GrB_Matrix *ptr_output) {
	GrB_Matrix output = *ptr_output;
	GxB_Iterator iterator;
	GxB_Iterator_new(&iterator);
	GrB_Info info = GxB_Matrix_Iterator_attach(iterator, output, NULL);
	info = GxB_Matrix_Iterator_seek(iterator, 0);
	AllPathsElem val;
	while (info != GxB_EXHAUSTED) {
		GxB_Iterator_get_UDT(iterator, (void *)&val);
		if (val.middle) {
			free(val.middle);
		}
		info = GxB_Matrix_Iterator_next(iterator);
	}

	GrB_free(&iterator);
	GrB_free(ptr_output);
}

void setup(char *msg) { LAGraph_Init(msg); }

void teardown(char *msg) { LAGraph_Finalize(msg); }
