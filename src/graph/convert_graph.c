#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "cfl_idr_graph_builder.h"
#include "internal/grb_utils.h"
#include "uthash.h"

#define ID_MAX_LEN 256

typedef struct {
	char id[ID_MAX_LEN];
	bool has_open;
	bool has_close;
	GrB_Matrix open_mat;
	GrB_Matrix close_mat;
	int idx;
	UT_hash_handle hh;
} BracketEntry;

static void hashmap_free(BracketEntry **map) {
	BracketEntry *entry, *tmp;
	HASH_ITER(hh, *map, entry, tmp) {
		HASH_DEL(*map, entry);
		if (!entry->has_open || !entry->has_close) {
			if (entry->has_open) {
				GrB_Matrix_free(&entry->open_mat);
			}
			if (entry->has_close) {
				GrB_Matrix_free(&entry->close_mat);
			}
		}
		free(entry);
	}
}

GrB_Info get_idr_graph(IdrGraph *out, GraphMatrices *gm,
					   const SymbolList *symbol_list, GrB_Index n,
					   const TerminalFormat *fmt) {
	if (!gm || !symbol_list || !fmt || !out) {
		return GrB_INVALID_VALUE;
	}

	GrB_Info info = GrB_SUCCESS;
	BracketEntry *par_map = NULL;
	BracketEntry *bra_map = NULL;
	bool has_normal = false;
	GrB_Matrix normal_mat = NULL;

	for (size_t i = 0; i < gm->count; i++) {
		MatrixSymbolInfo msi = gm->matrix_symbols[i];
		const char *label = symbol_list_get_str(symbol_list, msi.symbol_index);
		BracketType type = fmt->get_type(label);

		if (type == BRACKET_TYPE_UNKNOWN) {
			if (strcmp(label, "normal") == 0) {
				has_normal = true;
				normal_mat = gm->matrices[i];
				gm->matrices[i] = NULL;
			}
			continue;
		}

		GrB_Index nvals = 0;
		GRB_TRY(GrB_Matrix_nvals(&nvals, gm->matrices[i]));
		if (nvals == 0) {
			continue;
		}

		char id[ID_MAX_LEN];
		snprintf(id, sizeof(id), "%zu", msi.block_index);

		BracketEntry **map =
			(type == BRACKET_TYPE_PARENTHESES) ? &par_map : &bra_map;
		BracketEntry *entry = NULL;
		HASH_FIND_STR(*map, id, entry);
		if (!entry) {
			entry = calloc(1, sizeof(BracketEntry));
			strncpy(entry->id, id, sizeof(entry->id) - 1);
			HASH_ADD_STR(*map, id, entry);
		}

		if (fmt->is_open(label)) {
			entry->has_open = true;
			entry->open_mat = gm->matrices[i];
			gm->matrices[i] = NULL;
		} else {
			entry->has_close = true;
			entry->close_mat = gm->matrices[i];
			gm->matrices[i] = NULL;
		}
	}

	int64_t valid_par = 0, valid_bra = 0;
	BracketEntry *e, *tmp;
	HASH_ITER(hh, par_map, e, tmp) {
		if (e->has_open && e->has_close) {
			e->idx = valid_par++;
		}
	}
	HASH_ITER(hh, bra_map, e, tmp) {
		if (e->has_open && e->has_close) {
			e->idx = valid_bra++;
		}
	}

	out->n = n;
	out->n_par = valid_par;
	out->n_bra = valid_bra;
	out->normal = NULL;

	out->open_par = valid_par ? malloc(valid_par * sizeof(GrB_Matrix)) : NULL;
	out->close_par = valid_par ? malloc(valid_par * sizeof(GrB_Matrix)) : NULL;
	out->open_bra = valid_bra ? malloc(valid_bra * sizeof(GrB_Matrix)) : NULL;
	out->close_bra = valid_bra ? malloc(valid_bra * sizeof(GrB_Matrix)) : NULL;

	if ((valid_par && (!out->open_par || !out->close_par)) ||
		(valid_bra && (!out->open_bra || !out->close_bra))) {
		info = GrB_OUT_OF_MEMORY;
		goto cleanup;
	}

	HASH_ITER(hh, par_map, e, tmp) {
		if (e->has_open && e->has_close) {
			out->open_par[e->idx] = e->open_mat;
			out->close_par[e->idx] = e->close_mat;
		}
	}
	HASH_ITER(hh, bra_map, e, tmp) {
		if (e->has_open && e->has_close) {
			out->open_bra[e->idx] = e->open_mat;
			out->close_bra[e->idx] = e->close_mat;
		}
	}

	if (has_normal) {
		out->normal = normal_mat;
	}

cleanup:
	hashmap_free(&par_map);
	hashmap_free(&bra_map);
	if (info < GrB_SUCCESS) {
		free(out->open_par);
		free(out->close_par);
		free(out->open_bra);
		free(out->close_bra);
		out->open_par = out->close_par = NULL;
		out->open_bra = out->close_bra = NULL;
	}
	return info;
}
