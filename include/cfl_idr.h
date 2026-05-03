#pragma once

#include <GraphBLAS.h>
#include <stdbool.h>
#include <stdint.h>

typedef struct {
	GrB_Index n;
	int64_t n_par;
	GrB_Matrix *open_par;
	GrB_Matrix *close_par;
	int64_t n_bra;
	GrB_Matrix *open_bra;
	GrB_Matrix *close_bra;
	GrB_Matrix normal;
} IdrGraph;

void idr_graph_free(const IdrGraph *g);

void idr_graph_print(const IdrGraph *g);

typedef enum {
	IDR_UNKNOWN = -1,
	IDR_DEFAULT,
	IDR_PARITY,
	IDR_PARITY2,
	IDR_SE,
	IDR_PROJECT,
	IDR_EXCLUDE,
	IDR_ALL
} IdrGrammarType;

GrB_Info idr_get_under_approx(const IdrGraph *graph, bool valueflow,
							  GrB_Matrix *result);

GrB_Info idr_get_over_approx(const IdrGraph *graph, IdrGrammarType grammar_type,
							 GrB_Matrix under_approx, // NULL if not
							 GrB_Matrix *result, bool valueflow,
							 bool filter_empty); // whether to filter out empty paths

GrB_Info idr_remove_valueflow_unreachable(const IdrGraph *graph, IdrGraph *out,
										  char *msg);

GrB_Info idr_get_on_demand(const IdrGraph *graph, GrB_Matrix under_approx,
						   GrB_Matrix over_approx, bool parityD, GrB_Matrix *result,
						   bool valueflow, bool filter_empty, char *msg);
