#include "convert_graph.h"

#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

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
		free(entry);
	}
}

GrB_Info build_mr_graph(GrB_Matrix *matrices, const SymbolList *symbol_list,
						GrB_Index n, const TerminalFormat *fmt, MRGraph *out) {
	if (!matrices || !symbol_list || !fmt || !out) {
		return GrB_INVALID_VALUE;
	}

	BracketEntry *par_map = NULL;
	BracketEntry *bra_map = NULL;
	bool has_normal = false;

	for (size_t i = 0; i < symbol_list->count; i++) {
		const char *label = symbol_list->symbols[i].label;
		BracketType type = fmt->get_type(label);

		if (type == BRACKET_TYPE_UNKNOWN) {
			if (strcmp(label, "normal") == 0) {
				has_normal = true;
			}
			continue;
		}

		GrB_Index nvals = 0;
		GrB_Matrix_nvals(&nvals, matrices[i]);
		if (nvals == 0) {
			continue;
		}

		char id[ID_MAX_LEN];
		fmt->extract_id(label, id, sizeof(id));
		if (id[0] == '\0') {
			continue;
		}

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
			entry->open_mat = matrices[i];
		} else {
			entry->has_close = true;
			entry->close_mat = matrices[i];
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

	out->open_par = valid_par ? malloc(valid_par * sizeof(GrB_Matrix)) : NULL;
	out->close_par = valid_par ? malloc(valid_par * sizeof(GrB_Matrix)) : NULL;
	out->open_bra = valid_bra ? malloc(valid_bra * sizeof(GrB_Matrix)) : NULL;
	out->close_bra = valid_bra ? malloc(valid_bra * sizeof(GrB_Matrix)) : NULL;
	out->normal = NULL;

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
		for (size_t i = 0; i < symbol_list->count; i++) {
			if (strcmp(symbol_list->symbols[i].label, "normal") == 0) {
				out->normal = matrices[i];
				break;
			}
		}
	}

	hashmap_free(&par_map);
	hashmap_free(&bra_map);

	return GrB_SUCCESS;
}
