#include "extract_paths.h"

#include "internal/grb_utils.h"

GrB_Info extract_non_trivial_paths(GrB_Matrix *out, GrB_Matrix paths) {
	GrB_Info info = GrB_SUCCESS;
	GrB_Scalar thunk_zero = NULL;

	GRB_TRY(GrB_Scalar_new(&thunk_zero, GrB_INT64));
	GRB_TRY(GrB_Scalar_setElement_INT64(thunk_zero, 0));
	GRB_TRY(GrB_select(*out, NULL, NULL, GrB_OFFDIAG, paths, thunk_zero, NULL));

cleanup:
	GrB_free(&thunk_zero);
	return info;
}
