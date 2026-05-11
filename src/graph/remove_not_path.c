#include "remove_not_path.h"

#include <LAGraph.h>
#include <LAGraphX.h>
#include <stdlib.h>
#include <string.h>

#include "internal/grb_utils.h"

GrB_Info compute_sccs(SccResult *out, const IdrGraph *graph, GrB_Matrix adj) {
	GrB_Info info = GrB_SUCCESS;
	char msg[LAGRAPH_MSG_LEN];
	GrB_Vector scc_vec = NULL;
	GrB_Matrix reach = NULL;
	GrB_Matrix next = NULL;
	GrB_Matrix S = NULL;
	GrB_Matrix tmp_scc = NULL;
	GrB_Scalar s_true = NULL;
	GrB_Index *scc_ids = NULL;
	GrB_Index *remap = NULL;
	GrB_Index *raw_inds = NULL;
	GrB_Index *raw_vals = NULL;
	GrB_Index *v_ids = NULL;

	GrB_Index n = graph->n;

	if (n == 0) {
		out->n_scc = 0;
		out->scc_ids = NULL;
		GRB_TRY(GrB_Matrix_new(&out->scc_reach, GrB_BOOL, 0, 0));
		return GrB_SUCCESS;
	}

	// LAGraph_scc can't handle n=1
	if (n == 1) {
		out->scc_ids = malloc(sizeof(GrB_Index));
		if (!out->scc_ids) {
			info = GrB_OUT_OF_MEMORY;
			goto cleanup;
		}
		out->scc_ids[0] = 0;
		out->n_scc = 1;
		GRB_TRY(GrB_Matrix_new(&out->scc_reach, GrB_BOOL, 1, 1));
		GRB_TRY(GrB_Matrix_setElement_BOOL(out->scc_reach, true, 0, 0));
		return GrB_SUCCESS;
	}

	GRB_TRY(LAGraph_scc(&scc_vec, adj, msg));

	scc_ids = malloc(n * sizeof(GrB_Index));
	if (!scc_ids) {
		info = GrB_OUT_OF_MEMORY;
		goto cleanup;
	}

	GrB_Index nvals = 0;
	GRB_TRY(GrB_Vector_nvals(&nvals, scc_vec));

	raw_inds = malloc(nvals * sizeof(GrB_Index));
	raw_vals = malloc(nvals * sizeof(GrB_Index));
	if (!raw_inds || !raw_vals) {
		info = GrB_OUT_OF_MEMORY;
		goto cleanup;
	}

	GRB_TRY(GrB_Vector_extractTuples_UINT64(raw_inds, raw_vals, &nvals, scc_vec));

	// Remapping

	remap = malloc(n * sizeof(GrB_Index));
	if (!remap) {
		info = GrB_OUT_OF_MEMORY;
		goto cleanup;
	}
	memset(remap, 0xFF, n * sizeof(GrB_Index)); // 0xFF... == (GrB_Index)-1

	GrB_Index n_scc = 0;
	for (GrB_Index i = 0; i < nvals; i++) {
		GrB_Index raw = raw_vals[i];
		if (remap[raw] == (GrB_Index)-1) {
			remap[raw] = n_scc++;
		}
		scc_ids[raw_inds[i]] = remap[raw];
	}
	free(raw_inds);
	free(raw_vals);
	free(remap);
	raw_inds = NULL;
	raw_vals = NULL;
	remap = NULL;

	// Construction of the condensation matrix

	GRB_TRY(GrB_Matrix_new(&S, GrB_BOOL, n, n_scc));

	v_ids = malloc(n * sizeof(GrB_Index));
	if (!v_ids) {
		info = GrB_OUT_OF_MEMORY;
		goto cleanup;
	}
	for (GrB_Index i = 0; i < n; i++) {
		v_ids[i] = i;
	}

	GRB_TRY(GrB_Scalar_new(&s_true, GrB_BOOL));
	GRB_TRY(GrB_Scalar_setElement_BOOL(s_true, true));
	GRB_TRY(GxB_Matrix_build_Scalar(S, v_ids, scc_ids, s_true, n));
	free(v_ids);
	v_ids = NULL;
	GrB_Scalar_free(&s_true);

	// reach = S^T * adj * S
	GRB_TRY(GrB_Matrix_new(&tmp_scc, GrB_BOOL, n_scc, n));
	GRB_TRY(GrB_Matrix_new(&reach, GrB_BOOL, n_scc, n_scc));

	// tmp = S^T * adj
	GRB_TRY(GrB_mxm(tmp_scc, NULL, NULL, GrB_LOR_LAND_SEMIRING_BOOL, S, adj,
					GrB_DESC_T0));
	// reach = tmp * S
	GRB_TRY(
		GrB_mxm(reach, NULL, NULL, GrB_LOR_LAND_SEMIRING_BOOL, tmp_scc, S, NULL));

	GrB_Matrix_free(&S);
	GrB_Matrix_free(&tmp_scc);

	// Transitive closure

	// Reflexive pairs
	for (GrB_Index i = 0; i < n_scc; i++) {
		GrB_Matrix_setElement_BOOL(reach, true, i, i);
	}

	GrB_Index nnz_old = 0, nnz_new = 0;
	GRB_TRY(GrB_Matrix_nvals(&nnz_new, reach));
	GRB_TRY(GrB_Matrix_new(&next, GrB_BOOL, n_scc, n_scc));

	while (nnz_new > nnz_old) {
		nnz_old = nnz_new;

		GRB_TRY(GrB_mxm(next, NULL, NULL, GrB_LOR_LAND_SEMIRING_BOOL, reach, reach,
						NULL));

		GRB_TRY(GrB_Matrix_nvals(&nnz_new, next));

		GrB_Matrix temp_ptr = reach;
		reach = next;
		next = temp_ptr;
	}

	out->n_scc = n_scc;
	out->scc_ids = scc_ids;
	out->scc_reach = reach;
	scc_ids = NULL;
	reach = NULL;

cleanup:
	GrB_Vector_free(&scc_vec);
	GrB_Matrix_free(&reach);
	GrB_Matrix_free(&next);
	GrB_Matrix_free(&S);
	GrB_Matrix_free(&tmp_scc);
	GrB_Scalar_free(&s_true);
	free(scc_ids);
	free(remap);
	free(raw_inds);
	free(raw_vals);
	free(v_ids);
	return info;
}

void scc_result_free(SccResult *r) {
	if (!r) {
		return;
	}
	GrB_Matrix_free(&r->scc_reach);
	free(r->scc_ids);
	memset(r, 0, sizeof(*r));
}

bool is_all_pairs(GrB_Matrix over_approx, GrB_Index n) {
	GrB_Index nnz_approx = 0;
	GrB_Matrix_nvals(&nnz_approx, over_approx);
	return nnz_approx >= (n * n);
}

GrB_Info remove_not_path(IdrGraph *out, const IdrGraph *graph,
						 GrB_Matrix over_approx) {
	GrB_Info info = GrB_SUCCESS;

	SccResult sr = {0};
	GrB_Matrix S = NULL;
	GrB_Matrix allowed_scc = NULL;
	GrB_Matrix indirect = NULL;
	GrB_Matrix M = NULL;
	GrB_Matrix adj = NULL;
	GrB_Matrix keep_mask = NULL;
	GrB_Matrix tmp_scc = NULL;
	GrB_Matrix tmp_n = NULL;
	GrB_Scalar s_true = NULL;
	GrB_Index *v_ids = NULL;

	GrB_Index n = graph->n;

	if (n == 0) {
		GRB_TRY(GrB_Matrix_new(&keep_mask, GrB_BOOL, n, n));
		goto build_output;
	}

	GRB_TRY(idr_graph_to_adjacency(&adj, graph));

	// Calculating SCC
	GRB_TRY(compute_sccs(&sr, graph, adj));
	GrB_Index n_scc = sr.n_scc;

	// Build selection matrix S (n × n_scc): S[i][sr.scc_ids[i]] = true
	GRB_TRY(GrB_Matrix_new(&S, GrB_BOOL, n, n_scc));
	v_ids = malloc(n * sizeof(GrB_Index));
	if (!v_ids) {
		info = GrB_OUT_OF_MEMORY;
		goto cleanup;
	}
	for (GrB_Index i = 0; i < n; i++) {
		v_ids[i] = i;
	}

	GRB_TRY(GrB_Scalar_new(&s_true, GrB_BOOL));
	GRB_TRY(GrB_Scalar_setElement_BOOL(s_true, true));
	GRB_TRY(GxB_Matrix_build_Scalar(S, v_ids, sr.scc_ids, s_true, n));
	free(v_ids);
	v_ids = NULL;

	// Move over_approx to SCC level:
	// allowed_scc = S^T * over_approx * S
	GRB_TRY(GrB_Matrix_new(&allowed_scc, GrB_BOOL, n_scc, n_scc));
	GRB_TRY(GrB_Matrix_new(&tmp_scc, GrB_BOOL, n_scc, n));

	// tmp = S^T * over_approx
	GRB_TRY(GrB_mxm(tmp_scc, NULL, NULL, GrB_LOR_LAND_SEMIRING_BOOL, S, over_approx,
					GrB_DESC_T0));
	// allowed_scc = tmp * S
	GRB_TRY(GrB_mxm(allowed_scc, NULL, NULL, GrB_LOR_LAND_SEMIRING_BOOL, tmp_scc, S,
					NULL));
	GrB_Matrix_free(&tmp_scc);

	// Calculate indirect reachability:
	// indirect = sr.scc_reach^T * allowed_scc * sr.scc_reach^T
	GRB_TRY(GrB_Matrix_new(&indirect, GrB_BOOL, n_scc, n_scc));
	GRB_TRY(GrB_Matrix_new(&tmp_scc, GrB_BOOL, n_scc, n_scc));

	GRB_TRY(GrB_mxm(tmp_scc, NULL, NULL, GrB_LOR_LAND_SEMIRING_BOOL, sr.scc_reach,
					allowed_scc, GrB_DESC_T0));
	GRB_TRY(GrB_mxm(indirect, NULL, NULL, GrB_LOR_LAND_SEMIRING_BOOL, tmp_scc,
					sr.scc_reach, GrB_DESC_T1));
	GrB_Matrix_free(&tmp_scc);

	// Raise result back to n × n level:
	// M = S * indirect * S^T
	GRB_TRY(GrB_Matrix_new(&M, GrB_BOOL, n, n));
	GRB_TRY(GrB_Matrix_new(&tmp_n, GrB_BOOL, n, n_scc));

	GRB_TRY(
		GrB_mxm(tmp_n, NULL, NULL, GrB_LOR_LAND_SEMIRING_BOOL, S, indirect, NULL));
	GRB_TRY(
		GrB_mxm(M, NULL, NULL, GrB_LOR_LAND_SEMIRING_BOOL, tmp_n, S, GrB_DESC_T1));
	GrB_Matrix_free(&tmp_n);

	// Final mask: leaving only those edges that are in the original graph
	GRB_TRY(GrB_Matrix_new(&keep_mask, GrB_BOOL, n, n));
	GRB_TRY(GrB_eWiseMult(keep_mask, NULL, NULL, GrB_LAND, adj, M, NULL));

build_output:
	out->n = n;
	out->n_par = graph->n_par;
	out->n_bra = graph->n_bra;
	out->open_par = graph->n_par
						? (GrB_Matrix *)malloc(graph->n_par * sizeof(GrB_Matrix))
						: NULL;
	out->close_par = graph->n_par
						 ? (GrB_Matrix *)malloc(graph->n_par * sizeof(GrB_Matrix))
						 : NULL;
	out->open_bra = graph->n_par
						? (GrB_Matrix *)malloc(graph->n_bra * sizeof(GrB_Matrix))
						: NULL;
	out->close_bra = graph->n_par
						 ? (GrB_Matrix *)malloc(graph->n_bra * sizeof(GrB_Matrix))
						 : NULL;

	// Filtering all matrices of a graph

#define FILTER(field)                                                               \
	do {                                                                            \
		if (graph->field) {                                                         \
			GRB_TRY(GrB_Matrix_new(&out->field, GrB_BOOL, n, n));                   \
			GRB_TRY(GrB_eWiseMult(out->field, NULL, NULL, GrB_LAND, graph->field,   \
								  keep_mask, NULL));                                \
		} else {                                                                    \
			out->field = NULL;                                                      \
		}                                                                           \
	} while (0)

	for (int64_t i = 0; i < graph->n_par; i++) {
		FILTER(open_par[i]);
		FILTER(close_par[i]);
	}
	for (int64_t i = 0; i < graph->n_bra; i++) {
		FILTER(open_bra[i]);
		FILTER(close_bra[i]);
	}
	if (graph->normal) {
		GRB_TRY(GrB_Matrix_new(&out->normal, GrB_BOOL, n, n));
		GRB_TRY(GrB_eWiseMult(out->normal, NULL, NULL, GrB_LAND, graph->normal,
							  keep_mask, NULL));
	} else {
		out->normal = NULL;
	}

#undef FILTER

cleanup:
	scc_result_free(&sr);
	GrB_Matrix_free(&S);
	GrB_Matrix_free(&allowed_scc);
	GrB_Matrix_free(&indirect);
	GrB_Matrix_free(&M);
	GrB_Matrix_free(&adj);
	GrB_Matrix_free(&keep_mask);
	GrB_Matrix_free(&tmp_scc);
	GrB_Matrix_free(&tmp_n);
	GrB_Scalar_free(&s_true);
	free(v_ids);
	return info;
}
