#include "split_into_components.h"

#include <LAGraph.h>
#include <stdlib.h>

#include "internal/grb_utils.h"

typedef struct {
	GrB_Index comp_id;
	GrB_Index original_idx;
} VertexComp;

static int compare_vertex_comp(const void *a, const void *b) {
	GrB_Index comp1 = ((VertexComp *)a)->comp_id;
	GrB_Index comp2 = ((VertexComp *)b)->comp_id;
	return (comp1 > comp2) - (comp1 < comp2);
}

static GrB_Info IdrGraph_to_adjacency(GrB_Matrix *A_out, const IdrGraph *graph) {
	GrB_Info info = GrB_SUCCESS;
	GrB_Matrix A = NULL;

	GRB_TRY(GrB_Matrix_new(&A, GrB_BOOL, graph->n, graph->n));

	for (int64_t i = 0; i < graph->n_par; i++) {
		GRB_TRY(GrB_eWiseAdd(A, NULL, NULL, GrB_LOR, A, graph->open_par[i], NULL));
		GRB_TRY(GrB_eWiseAdd(A, NULL, NULL, GrB_LOR, A, graph->close_par[i], NULL));
	}

	for (int64_t i = 0; i < graph->n_bra; i++) {
		GRB_TRY(GrB_eWiseAdd(A, NULL, NULL, GrB_LOR, A, graph->open_bra[i], NULL));
		GRB_TRY(GrB_eWiseAdd(A, NULL, NULL, GrB_LOR, A, graph->close_bra[i], NULL));
	}

	if (graph->normal != NULL) {
		GRB_TRY(GrB_eWiseAdd(A, NULL, NULL, GrB_LOR, A, graph->normal, NULL));
	}

	*A_out = A;
	A = NULL;

cleanup:
	GrB_Matrix_free(&A);
	return info;
}

static GrB_Info extract_component_IdrGraph(IdrGraph *out, const IdrGraph *graph,
										   GrB_Index *verts, GrB_Index verts_count) {
	GrB_Info info = GrB_SUCCESS;

	out->n = verts_count;
	out->n_par = 0;
	out->n_bra = 0;
	out->normal = NULL;
	out->open_par = NULL;
	out->close_par = NULL;
	out->open_bra = NULL;
	out->close_bra = NULL;

	GrB_Matrix *open_par = (GrB_Matrix *)malloc(graph->n_par * sizeof(GrB_Matrix));
	GrB_Matrix *close_par = (GrB_Matrix *)malloc(graph->n_par * sizeof(GrB_Matrix));
	GrB_Matrix *open_bra = (GrB_Matrix *)malloc(graph->n_bra * sizeof(GrB_Matrix));
	GrB_Matrix *close_bra = (GrB_Matrix *)malloc(graph->n_bra * sizeof(GrB_Matrix));

	if ((graph->n_par > 0 && (!open_par || !close_par)) ||
		(graph->n_bra > 0 && (!open_bra || !close_bra))) {
		info = GrB_OUT_OF_MEMORY;
		goto cleanup;
	}

	for (int64_t i = 0; i < graph->n_par; i++) {
		GrB_Matrix op = NULL, cp = NULL;
		GRB_TRY(GrB_Matrix_new(&op, GrB_BOOL, verts_count, verts_count));
		GRB_TRY(GrB_Matrix_new(&cp, GrB_BOOL, verts_count, verts_count));

		GRB_TRY(GrB_Matrix_extract(op, NULL, NULL, graph->open_par[i], verts,
								   verts_count, verts, verts_count, NULL));
		GRB_TRY(GrB_Matrix_extract(cp, NULL, NULL, graph->close_par[i], verts,
								   verts_count, verts, verts_count, NULL));

		GrB_Index op_nnz = 0, cp_nnz = 0;
		GRB_TRY(GrB_Matrix_nvals(&op_nnz, op));
		GRB_TRY(GrB_Matrix_nvals(&cp_nnz, cp));

		if (op_nnz > 0 || cp_nnz > 0) {
			open_par[out->n_par] = op;
			close_par[out->n_par] = cp;
			out->n_par++;
		} else {
			GrB_Matrix_free(&op);
			GrB_Matrix_free(&cp);
		}
	}

	for (int64_t i = 0; i < graph->n_bra; i++) {
		GrB_Matrix ob = NULL, cb = NULL;
		GRB_TRY(GrB_Matrix_new(&ob, GrB_BOOL, verts_count, verts_count));
		GRB_TRY(GrB_Matrix_new(&cb, GrB_BOOL, verts_count, verts_count));

		GRB_TRY(GrB_Matrix_extract(ob, NULL, NULL, graph->open_bra[i], verts,
								   verts_count, verts, verts_count, NULL));
		GRB_TRY(GrB_Matrix_extract(cb, NULL, NULL, graph->close_bra[i], verts,
								   verts_count, verts, verts_count, NULL));

		GrB_Index ob_nnz = 0, cb_nnz = 0;
		GRB_TRY(GrB_Matrix_nvals(&ob_nnz, ob));
		GRB_TRY(GrB_Matrix_nvals(&cb_nnz, cb));

		if (ob_nnz > 0 || cb_nnz > 0) {
			open_bra[out->n_bra] = ob;
			close_bra[out->n_bra] = cb;
			out->n_bra++;
		} else {
			GrB_Matrix_free(&ob);
			GrB_Matrix_free(&cb);
		}
	}

	if (graph->normal != NULL) {
		GrB_Matrix nm = NULL;
		GRB_TRY(GrB_Matrix_new(&nm, GrB_BOOL, verts_count, verts_count));
		GRB_TRY(GrB_Matrix_extract(nm, NULL, NULL, graph->normal, verts, verts_count,
								   verts, verts_count, NULL));

		GrB_Index nm_nnz = 0;
		GRB_TRY(GrB_Matrix_nvals(&nm_nnz, nm));
		if (nm_nnz > 0) {
			out->normal = nm;
		} else {
			GrB_Matrix_free(&nm);
		}
	}

	out->open_par = open_par;
	out->close_par = close_par;
	out->open_bra = open_bra;
	out->close_bra = close_bra;
	return info;

cleanup:
	for (int64_t i = 0; i < out->n_par; i++) {
		GrB_Matrix_free(&open_par[i]);
		GrB_Matrix_free(&close_par[i]);
	}
	for (int64_t i = 0; i < out->n_bra; i++) {
		GrB_Matrix_free(&open_bra[i]);
		GrB_Matrix_free(&close_bra[i]);
	}
	free((void *)open_par);
	free((void *)close_par);
	free((void *)open_bra);
	free((void *)close_bra);
	return info;
}

GrB_Info split_IdrGraph_into_components(const IdrGraph *graph,
										IdrGraph **out_components,
										GrB_Index ***out_vertex_maps,
										GrB_Index *out_count) {
	GrB_Info info = GrB_SUCCESS;
	char msg[LAGRAPH_MSG_LEN];
	GrB_Index n = graph->n;

	GrB_Matrix A = NULL;
	GrB_Matrix A_T = NULL;
	LAGraph_Graph G = NULL;
	GrB_Vector component = NULL;
	GrB_Index *comp_ids = NULL;
	VertexComp *v_list = NULL;
	IdrGraph *components = NULL;
	GrB_Index **vertex_maps = NULL;

	// LAGr_ConnectedComponents can't handle n=1
	if (n <= 1) {
		*out_components = NULL;
		*out_vertex_maps = NULL;
		*out_count = 0;
		return GrB_SUCCESS;
	}

	GRB_TRY(IdrGraph_to_adjacency(&A, graph));

	GRB_TRY(GrB_Matrix_new(&A_T, GrB_BOOL, n, n));
	GRB_TRY(GrB_transpose(A_T, NULL, GrB_LOR, A, NULL));
	GRB_TRY(GrB_eWiseAdd(A, NULL, NULL, GrB_LOR, A, A_T, NULL));
	GrB_Matrix_free(&A_T);

	GRB_TRY(LAGraph_New(&G, &A, LAGraph_ADJACENCY_UNDIRECTED, msg));
	GRB_TRY(LAGr_ConnectedComponents(&component, G, msg));
	GRB_TRY(LAGraph_Delete(&G, msg));

	comp_ids = malloc(n * sizeof(GrB_Index));
	if (!comp_ids) {
		info = GrB_OUT_OF_MEMORY;
		goto cleanup;
	}

	GrB_Index actual_n = n;
	GRB_TRY(GrB_Vector_extractTuples_UINT64(NULL, comp_ids, &actual_n, component));
	GrB_Vector_free(&component);

	v_list = malloc(n * sizeof(VertexComp));
	if (!v_list) {
		info = GrB_OUT_OF_MEMORY;
		goto cleanup;
	}

	for (GrB_Index i = 0; i < n; i++) {
		v_list[i].comp_id = comp_ids[i];
		v_list[i].original_idx = i;
	}
	qsort(v_list, n, sizeof(VertexComp), compare_vertex_comp);

	GrB_Index unique_count = (n > 0) ? 1 : 0;
	for (GrB_Index i = 1; i < n; i++) {
		if (v_list[i].comp_id != v_list[i - 1].comp_id) {
			unique_count++;
		}
	}

	components = (IdrGraph *)malloc(unique_count * sizeof(IdrGraph));
	vertex_maps = (GrB_Index **)malloc(unique_count * sizeof(GrB_Index *));
	if (!components || !vertex_maps) {
		info = GrB_OUT_OF_MEMORY;
		goto cleanup;
	}

	GrB_Index valid_count = 0;
	GrB_Index start = 0;

	for (GrB_Index c = 0; c < unique_count; c++) {
		GrB_Index end = start;
		while (end < n && v_list[end].comp_id == v_list[start].comp_id) {
			end++;
		}

		GrB_Index count = end - start;
		// Ignore isolated vertices
		if (count > 1) {
			GrB_Index *verts = malloc(count * sizeof(GrB_Index));
			if (!verts) {
				info = GrB_OUT_OF_MEMORY;
				goto cleanup;
			}

			for (GrB_Index k = 0; k < count; k++) {
				verts[k] = v_list[start + k].original_idx;
			}

			GrB_Info extract_info = extract_component_IdrGraph(
				&components[valid_count], graph, verts, count);
			if (extract_info < GrB_SUCCESS) {
				free(verts);
				info = extract_info;
				goto cleanup;
			}

			vertex_maps[valid_count] = verts;
			valid_count++;
		}
		start = end;
	}

	if (valid_count == 0) {
		free(components);
		free((void *)vertex_maps);
		*out_components = NULL;
		*out_vertex_maps = NULL;
		*out_count = 0;
		components = NULL;
		vertex_maps = NULL;
		goto cleanup;
	}

	if (valid_count != unique_count) {
		IdrGraph *tmp_comp = realloc(components, valid_count * sizeof(IdrGraph));
		if (tmp_comp == NULL) {
			info = GrB_OUT_OF_MEMORY;
			goto cleanup;
		}
		components = tmp_comp;

		GrB_Index **tmp_vert_maps = (GrB_Index **)realloc(
			(void *)vertex_maps, valid_count * sizeof(GrB_Index *));
		if (tmp_vert_maps == NULL) {
			info = GrB_OUT_OF_MEMORY;
			goto cleanup;
		}
		vertex_maps = tmp_vert_maps;
	}

	*out_components = components;
	*out_vertex_maps = vertex_maps;
	*out_count = valid_count;
	components = NULL;
	vertex_maps = NULL;

cleanup:
	GrB_Matrix_free(&A);
	GrB_Matrix_free(&A_T);
	if (G) {
		LAGraph_Delete(&G, msg);
	}
	GrB_Vector_free(&component);
	free(comp_ids);
	free(v_list);
	if (components) {
		for (GrB_Index i = 0; i < valid_count; i++) {
			idr_graph_free(&components[i]);
		}
		free(components);
	}
	if (vertex_maps) {
		for (GrB_Index i = 0; i < valid_count; i++) {
			free(vertex_maps[i]);
		}
		free((void *)vertex_maps);
	}
	return info;
}
