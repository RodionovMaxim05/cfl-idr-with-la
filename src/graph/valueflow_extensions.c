#include "valueflow_extensions.h"

#include <LAGraph.h>

#define OK(f)                                                                       \
	do {                                                                            \
		info = (f);                                                                 \
		if (info != GrB_SUCCESS)                                                    \
			goto cleanup;                                                           \
	} while (0)

static GrB_Info get_scc_selector(GrB_Index n, GrB_Index n_scc,
								 const GrB_Index *scc_ids, GrB_Matrix *S_out) {
	GrB_Matrix S;
	GrB_Info info;
	GrB_Matrix_new(&S, GrB_BOOL, n, n_scc);
	for (GrB_Index i = 0; i < n; i++) {
		GrB_Matrix_setElement_BOOL(S, true, i, scc_ids[i]);
	}
	*S_out = S;
	return GrB_SUCCESS;
}

GrB_Info idr_remove_valueflow_unreachable(const IdrGraph *graph, IdrGraph *out,
										  char *msg) {
	GrB_Info info = GrB_SUCCESS;

	GrB_Vector v_has_open = NULL, v_has_close = NULL;
	GrB_Vector scc_sources = NULL, scc_sinks = NULL;
	GrB_Vector can_reach_sink = NULL, reached_from_source = NULL;
	GrB_Vector v_keep = NULL;
	GrB_Matrix S = NULL;
	GrB_Vector v_tmp = NULL;
	GrB_Matrix M_keep = NULL;

	GrB_Matrix adj = NULL;
	SccResult scc = {0};

	OK(build_adjacency(graph, &adj));
	OK(compute_sccs(graph, adj, &scc, msg));
	GrB_Index n_scc = scc.n_scc;

	GrB_Index n = graph->n;

	GrB_Vector_new(&v_has_open, GrB_BOOL, n);
	GrB_Vector_new(&v_has_close, GrB_BOOL, n);

	for (int64_t i = 0; i < graph->n_bra; i++) {
		GrB_Matrix_wait(graph->open_bra[i], GrB_MATERIALIZE);
		GrB_Matrix_wait(graph->close_bra[i], GrB_MATERIALIZE);

		// v_has_open[u] = true if exists edge u -(`[_i`)-> v
		GrB_reduce(v_has_open, NULL, GrB_LOR, GrB_LOR_MONOID_BOOL,
				   graph->open_bra[i], NULL);
		// v_has_close[v] = true if there is an edge u -(`]_i`)-> v
		// Use transpose to find the incoming elements of v
		GrB_reduce(v_has_close, NULL, GrB_LOR, GrB_LOR_MONOID_BOOL,
				   graph->close_bra[i], GrB_DESC_T0);
	}

	// Mapping nodes to SCC space
	OK(get_scc_selector(n, n_scc, scc.scc_ids, &S));
	GrB_Vector_new(&scc_sources, GrB_BOOL, n_scc);
	GrB_Vector_new(&scc_sinks, GrB_BOOL, n_scc);

	// scc_sources = v_has_open * S
	GrB_vxm(scc_sources, NULL, NULL, GrB_LOR_LAND_SEMIRING_BOOL, v_has_open, S,
			NULL);
	// scc_sinks = v_has_close * S
	GrB_vxm(scc_sinks, NULL, NULL, GrB_LOR_LAND_SEMIRING_BOOL, v_has_close, S, NULL);

	// Extending reachability in SCC space
	GrB_Vector_new(&reached_from_source, GrB_BOOL, n_scc);
	GrB_Vector_new(&can_reach_sink, GrB_BOOL, n_scc);

	// reached_from_source = scc_sources * SCC_reach
	GrB_vxm(reached_from_source, NULL, NULL, GrB_LOR_LAND_SEMIRING_BOOL, scc_sources,
			scc.scc_reach, NULL);
	// can_reach_sink = SCC_reach * scc_sinks
	GrB_mxv(can_reach_sink, NULL, NULL, GrB_LOR_LAND_SEMIRING_BOOL, scc.scc_reach,
			scc_sinks, NULL);

	// Final node mask: (v ∈ SCC_i) AND (Source -> SCC_i) AND (SCC_i -> Sink)
	GrB_Vector_new(&v_keep, GrB_BOOL, n);
	// Transferring from SCC back to nodes
	// v_keep = (S * reached_from_source) & (S * can_reach_sink)
	GrB_Vector_new(&v_tmp, GrB_BOOL, n);
	GrB_mxv(v_tmp, NULL, NULL, GrB_LOR_LAND_SEMIRING_BOOL, S, reached_from_source,
			NULL);
	GrB_mxv(v_keep, v_tmp, GrB_LAND, GrB_LOR_LAND_SEMIRING_BOOL, S, can_reach_sink,
			NULL);

	GrB_Matrix_new(&M_keep, GrB_BOOL, n, n);
	GrB_Matrix_diag(&M_keep, v_keep, 0);

	out->n = n;
	out->n_par = graph->n_par;
	out->n_bra = graph->n_bra;
	out->open_par = graph->n_par
						? (GrB_Matrix *)malloc(graph->n_par * sizeof(GrB_Matrix))
						: NULL;
	out->close_par = graph->n_par
						 ? (GrB_Matrix *)malloc(graph->n_par * sizeof(GrB_Matrix))
						 : NULL;
	out->open_bra = graph->n_bra
						? (GrB_Matrix *)malloc(graph->n_bra * sizeof(GrB_Matrix))
						: NULL;
	out->close_bra = graph->n_bra
						 ? (GrB_Matrix *)malloc(graph->n_bra * sizeof(GrB_Matrix))
						 : NULL;

	// Filtering all matrices of a graph

#define FILTER(field, out_field)                                                    \
	do {                                                                            \
		if (field) {                                                                \
			GrB_Matrix tmp_m = NULL;                                                \
			OK(GrB_Matrix_new(&(out_field), GrB_BOOL, n, n));                       \
			OK(GrB_Matrix_new(&tmp_m, GrB_BOOL, n, n));                             \
			OK(GrB_mxm(tmp_m, NULL, NULL, GrB_LOR_LAND_SEMIRING_BOOL, M_keep,       \
					   field, NULL));                                               \
			OK(GrB_mxm(out_field, NULL, NULL, GrB_LOR_LAND_SEMIRING_BOOL, tmp_m,    \
					   M_keep, NULL));                                              \
			GrB_Matrix_free(&tmp_m);                                                \
		} else {                                                                    \
			(out_field) = NULL;                                                     \
		}                                                                           \
	} while (0)

	for (int64_t i = 0; i < graph->n_par; i++) {
		FILTER(graph->open_par[i], out->open_par[i]);
		FILTER(graph->close_par[i], out->close_par[i]);
	}
	for (int64_t i = 0; i < graph->n_bra; i++) {
		FILTER(graph->open_bra[i], out->open_bra[i]);
		FILTER(graph->close_bra[i], out->close_bra[i]);
	}
	FILTER(graph->normal, out->normal);

#undef FILTER

cleanup:
	GrB_free(&v_has_open);
	GrB_free(&v_has_close);
	GrB_free(&scc_sources);
	GrB_free(&scc_sinks);
	GrB_free(&reached_from_source);
	GrB_free(&can_reach_sink);
	GrB_free(&S);
	GrB_free(&v_keep);
	GrB_free(&v_tmp);
	GrB_Matrix_free(&M_keep);

	return info;
}

GrB_Info filter_bracket_paths(const IdrGraph *graph, GrB_Matrix paths,
							  GrB_Matrix *filtered_out, char *msg) {
	GrB_Info info;
	GrB_Matrix S = NULL;
	GrB_Matrix R_node = NULL;
	GrB_Matrix Tmp = NULL;
	GrB_Matrix Combined = NULL;
	GrB_Index n = graph->n;

	GrB_Matrix adj = NULL;
	SccResult scc = {0};

	OK(build_adjacency(graph, &adj));
	OK(compute_sccs(graph, adj, &scc, msg));
	GrB_Index n_scc = scc.n_scc;

	// Expanding SCC reachability to all nodes
	OK(get_scc_selector(n, n_scc, scc.scc_ids, &S));
	GrB_Matrix_new(&Tmp, GrB_BOOL, n, n_scc);
	GrB_Matrix_new(&R_node, GrB_BOOL, n, n);

	// R_node = S * R_scc * S^T
	OK(GrB_mxm(Tmp, NULL, NULL, GrB_LOR_LAND_SEMIRING_BOOL, S, scc.scc_reach, NULL));
	OK(GrB_mxm(R_node, NULL, NULL, GrB_LOR_LAND_SEMIRING_BOOL, Tmp, S, GrB_DESC_T1));

	GrB_Matrix_new(&Combined, GrB_BOOL, n, n);

	// For each bracket ID, we check the condition: Open_i * R_node * Close_i
	for (int64_t i = 0; i < graph->n_bra; i++) {
		GrB_Matrix step1 = NULL, step2 = NULL;
		GrB_Matrix_new(&step1, GrB_BOOL, n, n);
		GrB_Matrix_new(&step2, GrB_BOOL, n, n);

		// step1 = open_bra[i] * R_node
		GrB_mxm(step1, NULL, NULL, GrB_LOR_LAND_SEMIRING_BOOL, graph->open_bra[i],
				R_node, NULL);
		// step2 = step1 * close_bra[i]
		GrB_mxm(step2, NULL, NULL, GrB_LOR_LAND_SEMIRING_BOOL, step1,
				graph->close_bra[i], NULL);

		// Accumulating result: Combined |= step2
		GrB_eWiseAdd(Combined, NULL, NULL, GrB_LOR, Combined, step2, NULL);

		GrB_free(&step1);
		GrB_free(&step2);
	}

	// Leave only those paths that were in original set
	// filtered = Combined & paths
	GrB_Matrix_new(filtered_out, GrB_BOOL, n, n);
	GrB_eWiseMult(*filtered_out, NULL, NULL, GrB_LAND, Combined, paths, NULL);

cleanup:
	GrB_free(&S);
	GrB_free(&R_node);
	GrB_free(&Tmp);
	GrB_free(&Combined);
	return info;
}
