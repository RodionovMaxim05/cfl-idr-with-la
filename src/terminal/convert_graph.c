#include "convert_graph.h"
#include <stdlib.h>
#include <string.h>

#include "uthash.h"

#define ID_MAX_LEN 256

typedef struct {
	char id[ID_MAX_LEN];
	int idx;
	UT_hash_handle hh;
} IdEntry;

static void hashmap_free(IdEntry **map) {
	IdEntry *entry, *tmp;
	HASH_ITER(hh, *map, entry, tmp) {
		HASH_DEL(*map, entry);
		free(entry);
	}
}

static IdEntry *get_or_create(IdEntry **map, const char *id, size_t *count) {
	IdEntry *entry = NULL;
	HASH_FIND_STR(*map, id, entry);
	if (!entry) {
		entry = malloc(sizeof(IdEntry));
		strncpy(entry->id, id, sizeof(entry->id) - 1);
		entry->idx = (int)(*count)++;
		HASH_ADD_STR(*map, id, entry);
	}
	return entry;
}

GrB_Info build_mr_graph(GrB_Matrix *matrices, const SymbolList *symbol_list,
						GrB_Index n, const TerminalFormat *fmt, MRGraph *out) {
	if (!matrices || !symbol_list || !fmt || !out) {
		return GrB_INVALID_VALUE;
	}

	IdEntry *par_map = NULL;
	IdEntry *bra_map = NULL;
	size_t par_count = 0;
	size_t bra_count = 0;
	bool has_normal = false;

	for (size_t i = 0; i < symbol_list->count; i++) {
		const char *label = symbol_list->symbols[i].label;
		BracketType type = fmt->get_type(label);
		int is_open = fmt->is_open(label);

		if (type == BRACKET_TYPE_UNKNOWN) {
			if (strcmp(label, "normal") == 0) {
				has_normal = true;
			}
			continue;
		}
		if (is_open != 1) {
			continue;
		}

		GrB_Index nvals = 0;
		GrB_Matrix_nvals(&nvals, matrices[i]);
		if (nvals == 0) {
			continue;
		}

		char id[256];
		fmt->extract_id(label, id, sizeof(id));
		if (id[0] == '\0') {
			continue;
		}

		if (type == BRACKET_TYPE_PARENTHESES) {
			get_or_create(&par_map, id, &par_count);
		} else if (type == BRACKET_TYPE_BRACKETS) {
			get_or_create(&bra_map, id, &bra_count);
		}
	}

	out->n = n;
	out->n_par = (int64_t)par_count;
	out->n_bra = (int64_t)bra_count;
	out->open_par =
		par_count ? (GrB_Matrix *)malloc(par_count * sizeof(GrB_Matrix)) : NULL;
	out->close_par =
		par_count ? (GrB_Matrix *)malloc(par_count * sizeof(GrB_Matrix)) : NULL;
	out->open_bra =
		bra_count ? (GrB_Matrix *)malloc(bra_count * sizeof(GrB_Matrix)) : NULL;
	out->close_bra =
		bra_count ? (GrB_Matrix *)malloc(bra_count * sizeof(GrB_Matrix)) : NULL;
	out->normal = NULL;

	for (size_t i = 0; i < par_count; i++) {
		out->open_par[i] = out->close_par[i] = NULL;
	}
	for (size_t i = 0; i < bra_count; i++) {
		out->open_bra[i] = out->close_bra[i] = NULL;
	}

	for (size_t i = 0; i < symbol_list->count; i++) {
		const char *label = symbol_list->symbols[i].label;
		BracketType type = fmt->get_type(label);
		int is_open = fmt->is_open(label);

		if (type == BRACKET_TYPE_UNKNOWN) {
			if (has_normal && strcmp(label, "normal") == 0) {
				out->normal = matrices[i];
			}
			continue;
		}
		if (is_open == -1) {
			continue;
		}

		char id[256];
		fmt->extract_id(label, id, sizeof(id));
		if (id[0] == '\0') {
			continue;
		}

		if (type == BRACKET_TYPE_PARENTHESES) {
			IdEntry *entry = NULL;
			HASH_FIND_STR(par_map, id, entry);
			if (!entry) {
				continue;
			}
			if (is_open) {
				out->open_par[entry->idx] = matrices[i];
			} else {
				out->close_par[entry->idx] = matrices[i];
			}
		} else if (type == BRACKET_TYPE_BRACKETS) {
			IdEntry *entry = NULL;
			HASH_FIND_STR(bra_map, id, entry);
			if (!entry) {
				continue;
			}
			if (is_open) {
				out->open_bra[entry->idx] = matrices[i];
			} else {
				out->close_bra[entry->idx] = matrices[i];
			}
		}
	}

	hashmap_free(&par_map);
	hashmap_free(&bra_map);

	return GrB_SUCCESS;
}
