#include "split_into_components.h"

#include <LAGraph.h>

GrB_Info MRGraph_to_adjacency(MRGraph *graph, GrB_Matrix *A_out) {
	GrB_Matrix A = NULL;
	GrB_Matrix_new(&A, GrB_BOOL, graph->n, graph->n);

	for (int64_t i = 0; i < graph->n_par; i++) {
		GrB_eWiseAdd(A, NULL, NULL, GrB_LOR, A, graph->open_par[i], NULL);
		GrB_eWiseAdd(A, NULL, NULL, GrB_LOR, A, graph->close_par[i], NULL);
	}

	for (int64_t i = 0; i < graph->n_bra; i++) {
		GrB_eWiseAdd(A, NULL, NULL, GrB_LOR, A, graph->open_bra[i], NULL);
		GrB_eWiseAdd(A, NULL, NULL, GrB_LOR, A, graph->close_bra[i], NULL);
	}

	if (graph->normal != NULL) {
		GrB_eWiseAdd(A, NULL, NULL, GrB_LOR, A, graph->normal, NULL);
	}

	*A_out = A;
	return GrB_SUCCESS;
}

GrB_Info extract_component_MRgraph(MRGraph *graph, GrB_Index *verts,
								   GrB_Index verts_count, MRGraph *out, char *msg) {
	out->n = verts_count;
	out->n_par = 0;
	out->n_bra = 0;
	out->normal = NULL;

	GrB_Matrix *open_par = (GrB_Matrix *)malloc(graph->n_par * sizeof(GrB_Matrix));
	GrB_Matrix *close_par = (GrB_Matrix *)malloc(graph->n_par * sizeof(GrB_Matrix));
	GrB_Matrix *open_bra = (GrB_Matrix *)malloc(graph->n_bra * sizeof(GrB_Matrix));
	GrB_Matrix *close_bra = (GrB_Matrix *)malloc(graph->n_bra * sizeof(GrB_Matrix));

	for (int64_t i = 0; i < graph->n_par; i++) {
		GrB_Matrix op = NULL, cp = NULL;
		GrB_Matrix_new(&op, GrB_BOOL, verts_count, verts_count);
		GrB_Matrix_new(&cp, GrB_BOOL, verts_count, verts_count);

		GrB_Matrix_extract(op, NULL, NULL, graph->open_par[i], verts, verts_count,
						   verts, verts_count, NULL);
		GrB_Matrix_extract(cp, NULL, NULL, graph->close_par[i], verts, verts_count,
						   verts, verts_count, NULL);

		GrB_Index op_nnz = 0, cp_nnz = 0;
		GrB_Matrix_nvals(&op_nnz, op);
		GrB_Matrix_nvals(&cp_nnz, cp);
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
		GrB_Matrix_new(&ob, GrB_BOOL, verts_count, verts_count);
		GrB_Matrix_new(&cb, GrB_BOOL, verts_count, verts_count);

		GrB_Matrix_extract(ob, NULL, NULL, graph->open_bra[i], verts, verts_count,
						   verts, verts_count, NULL);
		GrB_Matrix_extract(cb, NULL, NULL, graph->close_bra[i], verts, verts_count,
						   verts, verts_count, NULL);

		GrB_Index ob_nnz = 0, cb_nnz = 0;
		GrB_Matrix_nvals(&ob_nnz, ob);
		GrB_Matrix_nvals(&cb_nnz, cb);
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
		GrB_Matrix_new(&nm, GrB_BOOL, verts_count, verts_count);
		GrB_Matrix_extract(nm, NULL, NULL, graph->normal, verts, verts_count, verts,
						   verts_count, NULL);

		GrB_Index nm_nnz = 0;
		GrB_Matrix_nvals(&nm_nnz, nm);
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

	return GrB_SUCCESS;
}

GrB_Info split_MRGraph_into_components(MRGraph *graph, MRGraph **out_components,
									   GrB_Index ***out_vertex_maps,
									   GrB_Index *out_count, char *msg) {
	GrB_Index n = graph->n;

	GrB_Matrix A = NULL;
	MRGraph_to_adjacency(graph, &A);

	GrB_Matrix A_T = NULL, A_sym = NULL;
	GrB_Matrix_new(&A_T, GrB_BOOL, n, n);
	GrB_Matrix_new(&A_sym, GrB_BOOL, n, n);
	GrB_transpose(A_T, NULL, NULL, A, NULL);
	GrB_eWiseAdd(A_sym, NULL, NULL, GrB_LOR, A, A_T, NULL);
	GrB_Matrix_free(&A);
	GrB_Matrix_free(&A_T);

	LAGraph_Graph G = NULL;
	LAGraph_New(&G, &A_sym, LAGraph_ADJACENCY_UNDIRECTED, msg);
	LAGraph_Cached_AT(G, msg);

	GrB_Vector component = NULL;
	LAGr_ConnectedComponents(&component, G, msg);
	LAGraph_Delete(&G, msg);

	GrB_Index *comp_ids = malloc(n * sizeof(GrB_Index));
	for (GrB_Index i = 0; i < n; i++) {
		GrB_Vector_extractElement_UINT64(&comp_ids[i], component, i);
	}
	GrB_Vector_free(&component);

	GrB_Index *unique_ids = malloc(n * sizeof(GrB_Index));
	GrB_Index unique_count = 0;
	for (GrB_Index i = 0; i < n; i++) {
		bool found = false;
		for (GrB_Index u = 0; u < unique_count; u++) {
			if (unique_ids[u] == comp_ids[i]) {
				found = true;
				break;
			}
		}
		if (!found) {
			unique_ids[unique_count++] = comp_ids[i];
		}
	}

	MRGraph *components = (MRGraph *)malloc(unique_count * sizeof(MRGraph));
	GrB_Index **vertex_maps =
		(GrB_Index **)malloc(unique_count * sizeof(GrB_Index *));
	GrB_Index valid_count = 0;

	for (GrB_Index c = 0; c < unique_count; c++) {
		GrB_Index comp_id = unique_ids[c];

		GrB_Index *verts = malloc(n * sizeof(GrB_Index));
		GrB_Index verts_count = 0;
		for (GrB_Index i = 0; i < n; i++) {
			if (comp_ids[i] == comp_id) {
				verts[verts_count++] = i;
			}
		}

		// Remove isolated vertices
		if (verts_count <= 1) {
			free(verts);
			continue;
		}

		extract_component_MRgraph(graph, verts, verts_count,
								  &components[valid_count], msg);
		vertex_maps[valid_count] = verts;
		valid_count++;
	}

	if (valid_count == 0) {
		free(components);
		free((void *)vertex_maps);
		*out_components = NULL;
		*out_vertex_maps = NULL;
		*out_count = 0;

		goto cleanup;
	}

	if (valid_count != unique_count) {
		MRGraph *tmp_comp = realloc(components, valid_count * sizeof(MRGraph));
		if (tmp_comp == NULL) {
			free(components);
		}
		components = tmp_comp;

		GrB_Index **tmp_vert_maps = (GrB_Index **)realloc(
			(void *)vertex_maps, valid_count * sizeof(GrB_Index *));
		if (tmp_vert_maps == NULL) {
			free((void *)vertex_maps);
		}
		vertex_maps = tmp_vert_maps;
	}

	*out_components = components;
	*out_vertex_maps = vertex_maps;
	*out_count = valid_count;

cleanup:
	free(comp_ids);
	free(unique_ids);

	return GrB_SUCCESS;
}
