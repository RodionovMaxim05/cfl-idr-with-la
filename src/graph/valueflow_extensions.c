#include "valueflow_extensions.h"

#include <LAGraph.h>
#include <stdlib.h>

#include "internal/grb_utils.h"
#include "remove_not_path.h"

static GrB_Info get_scc_selector(GrB_Matrix *S_out, GrB_Index n, GrB_Index n_scc,
								 const GrB_Index *scc_ids) {
	GrB_Info info = GrB_SUCCESS;
	GrB_Matrix S = NULL;

	GRB_TRY(GrB_Matrix_new(&S, GrB_BOOL, n, n_scc));
	for (GrB_Index i = 0; i < n; i++) {
		GRB_TRY(GrB_Matrix_setElement_BOOL(S, true, i, scc_ids[i]));
	}

	*S_out = S;
	S = NULL;

cleanup:
	GrB_Matrix_free(&S);
	return info;
}

GrB_Info idr_remove_valueflow_unreachable(IdrGraph *out, const IdrGraph *graph) {
	GrB_Info info = GrB_SUCCESS;

	GrB_Vector v_has_open = NULL, v_has_close = NULL;
	GrB_Vector scc_sources = NULL, scc_sinks = NULL;
	GrB_Vector can_reach_sink = NULL, reached_from_source = NULL;
	GrB_Vector v_keep = NULL, v_tmp = NULL;
	GrB_Matrix scc_selector = NULL;
	GrB_Matrix keep_mask = NULL;
	GrB_Matrix adj = NULL;
	GrB_Matrix tmp_m = NULL;
	SccResult scc = {0};

	GRB_TRY(idr_graph_to_adjacency(&adj, graph));
	GRB_TRY(compute_sccs(&scc, graph, adj));
	GrB_Index n_scc = scc.n_scc;
	GrB_Index n = graph->n;

	GRB_TRY(GrB_Vector_new(&v_has_open, GrB_BOOL, n));
	GRB_TRY(GrB_Vector_new(&v_has_close, GrB_BOOL, n));

	for (int64_t i = 0; i < graph->n_bra; i++) {
		GRB_TRY(GrB_Matrix_wait(graph->open_bra[i], GrB_MATERIALIZE));
		GRB_TRY(GrB_Matrix_wait(graph->close_bra[i], GrB_MATERIALIZE));

		// v_has_open[u] = true if exists edge u -(`[_i`)-> v
		GRB_TRY(GrB_reduce(v_has_open, NULL, GrB_LOR, GrB_LOR_MONOID_BOOL,
						   graph->open_bra[i], NULL));
		// v_has_close[v] = true if there is an edge u -(`]_i`)-> v
		// Use transpose to find the incoming elements of v
		GRB_TRY(GrB_reduce(v_has_close, NULL, GrB_LOR, GrB_LOR_MONOID_BOOL,
						   graph->close_bra[i], GrB_DESC_T0));
	}

	// Mapping nodes to SCC space
	GRB_TRY(get_scc_selector(&scc_selector, n, n_scc, scc.scc_ids));

	GRB_TRY(GrB_Vector_new(&scc_sources, GrB_BOOL, n_scc));
	GRB_TRY(GrB_Vector_new(&scc_sinks, GrB_BOOL, n_scc));

	// scc_sources = v_has_open * S
	GRB_TRY(GrB_vxm(scc_sources, NULL, NULL, GrB_LOR_LAND_SEMIRING_BOOL, v_has_open,
					scc_selector, NULL));
	// scc_sinks = v_has_close * S
	GRB_TRY(GrB_vxm(scc_sinks, NULL, NULL, GrB_LOR_LAND_SEMIRING_BOOL, v_has_close,
					scc_selector, NULL));

	// Extending reachability in SCC space
	GRB_TRY(GrB_Vector_new(&reached_from_source, GrB_BOOL, n_scc));
	GRB_TRY(GrB_Vector_new(&can_reach_sink, GrB_BOOL, n_scc));

	// reached_from_source = scc_sources * SCC_reach
	GRB_TRY(GrB_vxm(reached_from_source, NULL, NULL, GrB_LOR_LAND_SEMIRING_BOOL,
					scc_sources, scc.scc_reach, NULL));
	// can_reach_sink = SCC_reach * scc_sinks
	GRB_TRY(GrB_mxv(can_reach_sink, NULL, NULL, GrB_LOR_LAND_SEMIRING_BOOL,
					scc.scc_reach, scc_sinks, NULL));

	// Final node mask: (v ∈ SCC_i) AND (Source -> SCC_i) AND (SCC_i -> Sink)
	GRB_TRY(GrB_Vector_new(&v_keep, GrB_BOOL, n));
	// Transferring from SCC back to nodes
	// v_keep = (S * reached_from_source) & (S * can_reach_sink)
	GRB_TRY(GrB_Vector_new(&v_tmp, GrB_BOOL, n));
	GRB_TRY(GrB_mxv(v_tmp, NULL, NULL, GrB_LOR_LAND_SEMIRING_BOOL, scc_selector,
					reached_from_source, NULL));
	GRB_TRY(GrB_mxv(v_keep, v_tmp, GrB_LAND, GrB_LOR_LAND_SEMIRING_BOOL,
					scc_selector, can_reach_sink, NULL));

	GRB_TRY(GrB_Matrix_diag(&keep_mask, v_keep, 0));

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
			GrB_Matrix_free(&tmp_m);                                                \
			GRB_TRY(GrB_Matrix_new(&(out_field), GrB_BOOL, n, n));                  \
			GRB_TRY(GrB_Matrix_new(&tmp_m, GrB_BOOL, n, n));                        \
			GRB_TRY(GrB_mxm(tmp_m, NULL, NULL, GrB_LOR_LAND_SEMIRING_BOOL,          \
							keep_mask, field, NULL));                               \
			GRB_TRY(GrB_mxm(out_field, NULL, NULL, GrB_LOR_LAND_SEMIRING_BOOL,      \
							tmp_m, keep_mask, NULL));                               \
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
	GrB_Vector_free(&v_has_open);
	GrB_Vector_free(&v_has_close);
	GrB_Vector_free(&scc_sources);
	GrB_Vector_free(&scc_sinks);
	GrB_Vector_free(&reached_from_source);
	GrB_Vector_free(&can_reach_sink);
	GrB_Vector_free(&v_keep);
	GrB_Vector_free(&v_tmp);
	GrB_Matrix_free(&scc_selector);
	GrB_Matrix_free(&keep_mask);
	GrB_Matrix_free(&adj);
	scc_result_free(&scc);
	GrB_Matrix_free(&tmp_m);

	return info;
}

GrB_Info filter_bracket_paths(GrB_Matrix *out, const IdrGraph *graph,
							  GrB_Matrix paths) {
	GrB_Info info = GrB_SUCCESS;

	GrB_Matrix scc_selector = NULL;
	GrB_Matrix node_reach = NULL;
	GrB_Matrix scc_tmp = NULL;
	GrB_Matrix bracket_reachable = NULL;
	GrB_Matrix open_reach = NULL, open_close_reach = NULL;
	GrB_Matrix adj = NULL;
	SccResult scc = {0};
	GrB_Index n = graph->n;

	GRB_TRY(idr_graph_to_adjacency(&adj, graph));
	GRB_TRY(compute_sccs(&scc, graph, adj));
	GrB_Index n_scc = scc.n_scc;

	// Expanding SCC reachability to all nodes
	GRB_TRY(get_scc_selector(&scc_selector, n, n_scc, scc.scc_ids));

	GRB_TRY(GrB_Matrix_new(&scc_tmp, GrB_BOOL, n, n_scc));
	GRB_TRY(GrB_Matrix_new(&node_reach, GrB_BOOL, n, n));
	// R_node = S * R_scc * S^T
	GRB_TRY(GrB_mxm(scc_tmp, NULL, NULL, GrB_LOR_LAND_SEMIRING_BOOL, scc_selector,
					scc.scc_reach, NULL));
	GRB_TRY(GrB_mxm(node_reach, NULL, NULL, GrB_LOR_LAND_SEMIRING_BOOL, scc_tmp,
					scc_selector, GrB_DESC_T1));

	GRB_TRY(GrB_Matrix_new(&bracket_reachable, GrB_BOOL, n, n));

	// For each bracket ID, we check the condition: Open_i * R_node * Close_i
	for (int64_t i = 0; i < graph->n_bra; i++) {
		open_reach = NULL;
		open_close_reach = NULL;

		GRB_TRY(GrB_Matrix_new(&open_reach, GrB_BOOL, n, n));
		GRB_TRY(GrB_Matrix_new(&open_close_reach, GrB_BOOL, n, n));

		// step1 = open_bra[i] * R_node
		GRB_TRY(GrB_mxm(open_reach, NULL, NULL, GrB_LOR_LAND_SEMIRING_BOOL,
						graph->open_bra[i], node_reach, NULL));
		// step2 = step1 * close_bra[i]
		GRB_TRY(GrB_mxm(open_close_reach, NULL, NULL, GrB_LOR_LAND_SEMIRING_BOOL,
						open_reach, graph->close_bra[i], NULL));

		// Accumulating result: Combined |= step2
		GRB_TRY(GrB_eWiseAdd(bracket_reachable, NULL, NULL, GrB_LOR,
							 bracket_reachable, open_close_reach, NULL));

		GrB_Matrix_free(&open_reach);
		GrB_Matrix_free(&open_close_reach);
	}

	// Leave only those paths that were in original set
	// out = combined & paths
	GRB_TRY(GrB_Matrix_new(out, GrB_BOOL, n, n));
	GRB_TRY(
		GrB_eWiseMult(*out, NULL, NULL, GrB_LAND, bracket_reachable, paths, NULL));

cleanup:
	GrB_Matrix_free(&scc_selector);
	GrB_Matrix_free(&node_reach);
	GrB_Matrix_free(&scc_tmp);
	GrB_Matrix_free(&bracket_reachable);
	GrB_Matrix_free(&open_reach);
	GrB_Matrix_free(&open_close_reach);
	GrB_Matrix_free(&adj);
	scc_result_free(&scc);

	return info;
}
