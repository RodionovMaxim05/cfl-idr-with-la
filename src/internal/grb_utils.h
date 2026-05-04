#pragma once

#include <GraphBLAS.h>
#include <stdio.h>

#ifndef IDR_GRB_CATCH
#define IDR_GRB_CATCH(_info)                                                        \
	{                                                                               \
		info = (_info);                                                             \
		goto cleanup;                                                               \
	}
#endif

#define GRB_TRY(method)                                                             \
	do {                                                                            \
		GrB_Info _grb_info = (method);                                              \
		if (_grb_info < GrB_SUCCESS) {                                              \
			fprintf(stderr, "GraphBLAS error (%s:%d): %d\n", __FILE__, __LINE__,    \
					(int)_grb_info);                                                \
			IDR_GRB_CATCH(_grb_info);                                               \
		}                                                                           \
	} while (0)
