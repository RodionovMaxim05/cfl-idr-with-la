#include "remove_not_path.h"

#include <stdlib.h>
#include <string.h>

#include "GraphBLAS.h"
#include "LAGraph.h"
#include "LAGraphX.h"

#define OK(f)                                                                       \
	do {                                                                            \
		info = (f);                                                                 \
		if (info != GrB_SUCCESS)                                                    \
			goto cleanup;                                                           \
	} while (0)

static GrB_Info build_adjacency(const MRGraph *graph, GrB_Matrix *adj_out) {
	GrB_Info info;
	GrB_Matrix adj = NULL;

	OK(GrB_Matrix_new(&adj, GrB_BOOL, graph->n, graph->n));

	for (int64_t i = 0; i < graph->n_par; i++) {
		OK(GrB_eWiseAdd(adj, NULL, GrB_LOR, GrB_LOR, adj, graph->open_par[i], NULL));
		OK(GrB_eWiseAdd(adj, NULL, GrB_LOR, GrB_LOR, adj, graph->close_par[i],
						NULL));
	}

	for (int64_t i = 0; i < graph->n_bra; i++) {
		OK(GrB_eWiseAdd(adj, NULL, GrB_LOR, GrB_LOR, adj, graph->open_bra[i], NULL));

		OK(GrB_eWiseAdd(adj, NULL, GrB_LOR, GrB_LOR, adj, graph->close_bra[i],
						NULL));
	}

	if (graph->normal != NULL) {
		OK(GrB_eWiseAdd(adj, NULL, GrB_LOR, GrB_LOR, adj, graph->normal, NULL));
	}

	*adj_out = adj;
	return GrB_SUCCESS;

cleanup:
	GrB_Matrix_free(&adj);
	return info;
}

GrB_Info compute_sccs(const MRGraph *graph, SccResult *out, char *msg) {
	GrB_Info info = GrB_SUCCESS;
	GrB_Matrix adj = NULL;
	GrB_Vector scc_vec = NULL;
	GrB_Matrix reach = NULL;
	GrB_Matrix next = NULL;

	GrB_Index *scc_ids = NULL;
	GrB_Index *rows = NULL, *cols = NULL;
	GrB_Index *remap = NULL;

	GrB_Index n = graph->n;

	OK(build_adjacency(graph, &adj));

	// LAGraph_scc can't handle n=1
	if (n == 1) {
		out->n_scc = 1;
		out->scc_ids = malloc(sizeof(GrB_Index));
		out->scc_ids[0] = 0;
		GrB_Matrix_new(&out->scc_reach, GrB_BOOL, 1, 1);
		GrB_Matrix_setElement_BOOL(out->scc_reach, true, 0, 0);
		GrB_Matrix_free(&adj);
		return GrB_SUCCESS;
	}

	OK(LAGraph_scc(&scc_vec, adj, msg));

	scc_ids = malloc(n * sizeof(GrB_Index));
	if (!scc_ids) {
		info = GrB_OUT_OF_MEMORY;
		goto cleanup;
	}

	GrB_Index nvals = 0;
	OK(GrB_Vector_nvals(&nvals, scc_vec));
	GrB_Index *raw_inds = malloc(nvals * sizeof(GrB_Index));
	GrB_Index *raw_vals = malloc(nvals * sizeof(GrB_Index));
	OK(GrB_Vector_extractTuples_UINT64(raw_inds, raw_vals, &nvals, scc_vec));

	// Remapping

	remap = malloc(n * sizeof(GrB_Index));
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

	OK(GrB_Matrix_nvals(&nvals, adj));
	rows = malloc(nvals * sizeof(GrB_Index));
	cols = malloc(nvals * sizeof(GrB_Index));
	bool *vals = malloc(nvals * sizeof(bool));

	OK(GrB_Matrix_extractTuples_BOOL(rows, cols, vals, &nvals, adj));

	GrB_Index k_out = 0;
	for (GrB_Index k = 0; k < nvals; k++) {
		GrB_Index cu = scc_ids[rows[k]];
		GrB_Index cv = scc_ids[cols[k]];
		if (cu != cv) {
			rows[k_out] = cu;
			cols[k_out] = cv;
			vals[k_out] = true;
			k_out++;
		}
	}

	OK(GrB_Matrix_new(&reach, GrB_BOOL, n_scc, n_scc));
	if (k_out > 0) {
		OK(GrB_Matrix_build_BOOL(reach, rows, cols, vals, k_out, GrB_LOR));
	}
	free(rows);
	free(cols);
	free(vals);
	rows = NULL;
	cols = NULL;
	vals = NULL;

	// Transitive closure

	// Reflexive pairs
	for (GrB_Index i = 0; i < n_scc; i++) {
		GrB_Matrix_setElement_BOOL(reach, true, i, i);
	}

	GrB_Index nnz_old = 0, nnz_new = 0;
	OK(GrB_Matrix_nvals(&nnz_new, reach));
	OK(GrB_Matrix_new(&next, GrB_BOOL, n_scc, n_scc));

	while (nnz_new > nnz_old) {
		nnz_old = nnz_new;

		OK(GrB_mxm(next, NULL, NULL, GrB_LOR_LAND_SEMIRING_BOOL, reach, reach,
				   NULL));

		OK(GrB_Matrix_nvals(&nnz_new, next));

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
	GrB_Matrix_free(&adj);
	GrB_Matrix_free(&reach);
	GrB_Vector_free(&scc_vec);
	GrB_Matrix_free(&next);
	free(scc_ids);
	free(remap);
	free(rows);
	free(cols);
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

GrB_Info remove_not_path(const MRGraph *graph, GrB_Matrix over_approx, MRGraph *out,
						 char *msg) {
	GrB_Info info = GrB_SUCCESS;
	SccResult sr = {0};
	GrB_Matrix S = NULL;
	GrB_Matrix allowed_scc = NULL;
	GrB_Matrix indirect = NULL;
	GrB_Matrix M = NULL;
	GrB_Matrix adj = NULL;
	GrB_Matrix keep_mask = NULL;

	GrB_Index n = graph->n;

	// Calculating SCC
	OK(compute_sccs(graph, &sr, msg));
	GrB_Index n_scc = sr.n_scc;

	// Build selection matrix S (n × n_scc): S[i][sr.scc_ids[i]] = true
	OK(GrB_Matrix_new(&S, GrB_BOOL, n, n_scc));
	GrB_Index *v_ids = malloc(n * sizeof(GrB_Index));
	for (GrB_Index i = 0; i < n; i++) {
		v_ids[i] = i;
	}

	GrB_Scalar s_true = NULL;
	OK(GrB_Scalar_new(&s_true, GrB_BOOL));
	OK(GrB_Scalar_setElement_BOOL(s_true, true));

	OK(GxB_Matrix_build_Scalar(S, v_ids, sr.scc_ids, s_true, n));
	free(v_ids);

	// Move over_approx to SCC level:
	// allowed_scc = S^T * over_approx * S
	OK(GrB_Matrix_new(&allowed_scc, GrB_BOOL, n_scc, n_scc));
	GrB_Matrix tmp_scc = NULL;
	OK(GrB_Matrix_new(&tmp_scc, GrB_BOOL, n_scc, n));

	// tmp = S^T * over_approx
	OK(GrB_mxm(tmp_scc, NULL, NULL, GrB_LOR_LAND_SEMIRING_BOOL, S, over_approx,
			   GrB_DESC_T0));
	// allowed_scc = tmp * S
	OK(GrB_mxm(allowed_scc, NULL, NULL, GrB_LOR_LAND_SEMIRING_BOOL, tmp_scc, S,
			   NULL));
	GrB_Matrix_free(&tmp_scc);

	// Calculate indirect reachability:
	// indirect = sr.scc_reach^T * allowed_scc * sr.scc_reach^T
	OK(GrB_Matrix_new(&indirect, GrB_BOOL, n_scc, n_scc));
	OK(GrB_Matrix_new(&tmp_scc, GrB_BOOL, n_scc, n_scc));

	OK(GrB_mxm(tmp_scc, NULL, NULL, GrB_LOR_LAND_SEMIRING_BOOL, sr.scc_reach,
			   allowed_scc, GrB_DESC_T0));
	OK(GrB_mxm(indirect, NULL, NULL, GrB_LOR_LAND_SEMIRING_BOOL, tmp_scc,
			   sr.scc_reach, GrB_DESC_T1));
	GrB_Matrix_free(&tmp_scc);

	// Raise result back to n × n level:
	// M = S * indirect * S^T
	OK(GrB_Matrix_new(&M, GrB_BOOL, n, n));
	GrB_Matrix tmp_n = NULL;
	OK(GrB_Matrix_new(&tmp_n, GrB_BOOL, n, n_scc));

	OK(GrB_mxm(tmp_n, NULL, NULL, GrB_LOR_LAND_SEMIRING_BOOL, S, indirect, NULL));
	OK(GrB_mxm(M, NULL, NULL, GrB_LOR_LAND_SEMIRING_BOOL, tmp_n, S, GrB_DESC_T1));
	GrB_Matrix_free(&tmp_n);

	// Final mask: leaving only those edges that are in the original graph
	OK(build_adjacency(graph, &adj));
	OK(GrB_Matrix_new(&keep_mask, GrB_BOOL, n, n));
	OK(GrB_eWiseMult(keep_mask, NULL, NULL, GrB_LAND, adj, M, NULL));

	out->n = n;
	out->n_par = graph->n_par;
	out->n_bra = graph->n_bra;
	out->open_par = (GrB_Matrix *)malloc(graph->n_par * sizeof(GrB_Matrix));
	out->close_par = (GrB_Matrix *)malloc(graph->n_par * sizeof(GrB_Matrix));
	out->open_bra = (GrB_Matrix *)malloc(graph->n_bra * sizeof(GrB_Matrix));
	out->close_bra = (GrB_Matrix *)malloc(graph->n_bra * sizeof(GrB_Matrix));

	// Filtering all matrices of a graph

#define FILTER(field)                                                               \
	do {                                                                            \
		if (graph->field) {                                                         \
			OK(GrB_Matrix_new(&out->field, GrB_BOOL, n, n));                        \
			OK(GrB_eWiseMult(out->field, NULL, NULL, GrB_LAND, graph->field,        \
							 keep_mask, NULL));                                     \
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
		OK(GrB_Matrix_new(&out->normal, GrB_BOOL, n, n));
		OK(GrB_eWiseMult(out->normal, NULL, NULL, GrB_LAND, graph->normal, keep_mask,
						 NULL));
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
	return info;
}
