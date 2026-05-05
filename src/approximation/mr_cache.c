#include "mr_cache.h"

#include <LAGraphX.h>
#include <stdlib.h>
#include <string.h>

#define FNV_OFFSET_BASIS 0xCBF29CE484222325ULL
#define FNV_PRIME 0x00000100000001B3ULL

static inline uint64_t fnv1a_update(uint64_t hash, const void *data, size_t len) {
	const uint8_t *p = (const uint8_t *)data;
	for (size_t i = 0; i < len; i++) {
		hash ^= (uint64_t)p[i];
		hash *= FNV_PRIME;
	}
	return hash;
}

static int cmp_u64(const void *a, const void *b) {
	uint64_t x = *(const uint64_t *)a;
	uint64_t y = *(const uint64_t *)b;
	return (x > y) - (x < y);
}

static uint64_t fnv1a_feed_matrix(uint64_t hash, GrB_Matrix m, const char *tag,
								  int64_t idx) {
	hash = fnv1a_update(hash, tag, strlen(tag));
	hash = fnv1a_update(hash, &idx, sizeof(idx));

	if (m == NULL) {
		return hash;
	}

	GrB_Index nvals = 0;
	GrB_Matrix_nvals(&nvals, m);
	if (nvals == 0) {
		return hash;
	}

	GrB_Index *rows = malloc(nvals * sizeof(GrB_Index));
	GrB_Index *cols = malloc(nvals * sizeof(GrB_Index));
	uint64_t *pairs = malloc(nvals * sizeof(uint64_t));
	if (!rows || !cols || !pairs) {
		goto cleanup;
	}

	GrB_Matrix_extractTuples_BOOL(rows, cols, NULL, &nvals, m);
	for (GrB_Index k = 0; k < nvals; k++) {
		pairs[k] = ((uint64_t)rows[k] << 32) | (uint32_t)cols[k];
	}
	qsort(pairs, nvals, sizeof(uint64_t), cmp_u64);
	hash = fnv1a_update(hash, pairs, nvals * sizeof(uint64_t));

cleanup:
	free(rows);
	free(cols);
	free(pairs);
	return hash;
}

uint64_t get_graph_cache_hash(const IdrGraph *graph) {
	uint64_t h = FNV_OFFSET_BASIS;
	h = fnv1a_update(h, &graph->n, sizeof(graph->n));
	for (int64_t i = 0; i < graph->n_par; i++) {
		h = fnv1a_feed_matrix(h, graph->open_par[i], "op", i);
		h = fnv1a_feed_matrix(h, graph->close_par[i], "cp", i);
	}
	for (int64_t i = 0; i < graph->n_bra; i++) {
		h = fnv1a_feed_matrix(h, graph->open_bra[i], "ob", i);
		h = fnv1a_feed_matrix(h, graph->close_bra[i], "cb", i);
	}
	if (graph->normal != NULL) {
		h = fnv1a_feed_matrix(h, graph->normal, "nm", 0);
	}

	return h;
}

void mr_cache_init(MRCache *c) { memset(c, 0, sizeof(*c)); }

static void free_step_result(MRStepResult *r) {
	if (!r || !r->matrices) {
		return;
	}
	LAGraph_CFL_AllPaths_free_outputs(r->matrices, r->count, &r->all_paths_t);
	GrB_free(&r->all_paths_t);
	memset(r, 0, sizeof(*r));
}

void mr_cache_free(MRCache *c) {
	if (!c) {
		return;
	}
	for (size_t i = 0; i < c->count; i++) {
		free_step_result(&c->entries[i].value);
	}
	free(c->entries);
	memset(c, 0, sizeof(*c));
}

const MRStepResult *mr_cache_lookup(const MRCache *c, uint64_t graph_key,
									uint32_t grammar_tag) {
	if (!c) {
		return NULL;
	}
	for (size_t i = 0; i < c->count; i++) {
		const MRCacheEntry *e = &c->entries[i];
		if (e->graph_key == graph_key && e->grammar_tag == grammar_tag) {
			return &e->value;
		}
	}
	return NULL;
}

GrB_Info mr_cache_insert(MRCache *c, uint64_t graph_key, uint32_t grammar_tag,
						 GrB_Matrix *paths_matrices, int64_t nonterms_count,
						 GrB_Type all_paths_t) {
	if (!c) {
		return GrB_SUCCESS;
	}

	if (c->count == c->capacity) {
		size_t new_cap = c->capacity ? c->capacity * 2 : 8;
		MRCacheEntry *tmp_buf = realloc(c->entries, new_cap * sizeof(MRCacheEntry));
		if (!tmp_buf) {
			return GrB_OUT_OF_MEMORY;
		}
		c->entries = tmp_buf;
		c->capacity = new_cap;
	}

	MRCacheEntry *e = &c->entries[c->count++];
	e->graph_key = graph_key;
	e->grammar_tag = grammar_tag;
	e->value.matrices = paths_matrices;
	e->value.count = nonterms_count;
	e->value.all_paths_t = all_paths_t;

	return GrB_SUCCESS;
}
