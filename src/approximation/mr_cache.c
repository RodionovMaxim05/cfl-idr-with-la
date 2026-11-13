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
	if (!rows || !cols) {
		goto cleanup;
	}

	GrB_Matrix_extractTuples_BOOL(rows, cols, NULL, &nvals, m);
	for (GrB_Index k = 0; k < nvals; k++) {
		rows[k] = ((uint64_t)rows[k] << 32) | (uint32_t)cols[k];
	}
	qsort(rows, nvals, sizeof(uint64_t), cmp_u64);
	hash = fnv1a_update(hash, rows, nvals * sizeof(uint64_t));

cleanup:
	free(rows);
	free(cols);
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

static void free_reachability_result(MRReachabilityResult *r) {
	if (!r || !r->matrices) {
		return;
	}
	for (int64_t i = 0; i < r->count; i++) {
		GrB_free(&r->matrices[i]);
	}
	LAGraph_Free((void **)&r->matrices, NULL);
	memset(r, 0, sizeof(*r));
}

static void free_extracted_result(MRExtractedResult *r) {
	if (!r) {
		return;
	}
	if (r->edges) {
		for (int64_t i = 0; i < r->edges_count; i++) {
			if (r->edges[i] != NULL) {
				GrB_free(&r->edges[i]);
			}
		}
		free((void *)r->edges);
	}
	if (r->reachability != NULL) {
		GrB_free(&r->reachability);
	}
	memset(r, 0, sizeof(*r));
}

void mr_cache_free(MRCache *c) {
	if (!c) {
		return;
	}

	for (size_t i = 0; i < c->reachability_count; i++) {
		free_reachability_result(&c->reachability[i].value);
	}
	free(c->reachability);

	for (size_t i = 0; i < c->extracted_count; i++) {
		free_extracted_result(&c->extracted[i].value);
	}
	free(c->extracted);

	memset(c, 0, sizeof(*c));
}

const MRReachabilityResult *mr_cache_lookup_reachability(const MRCache *c,
														 uint64_t graph_key,
														 uint32_t grammar_tag) {
	if (!c) {
		return NULL;
	}

	for (size_t i = 0; i < c->reachability_count; i++) {
		const MRReachabilityEntry *e = &c->reachability[i];
		if (e->graph_key == graph_key && e->grammar_tag == grammar_tag) {
			return &e->value;
		}
	}

	return NULL;
}

const MRExtractedResult *mr_cache_lookup_extracted(const MRCache *c,
												   uint64_t graph_key,
												   uint32_t grammar_tag) {
	if (!c) {
		return NULL;
	}

	for (size_t i = 0; i < c->extracted_count; i++) {
		const MRExtractedEntry *e = &c->extracted[i];
		if (e->graph_key == graph_key && e->grammar_tag == grammar_tag) {
			return &e->value;
		}
	}

	return NULL;
}

GrB_Info mr_cache_insert_reachability(MRCache *c, uint64_t graph_key,
									  uint32_t grammar_tag,
									  GrB_Matrix *reachability_matrices,
									  int64_t symbols_amount) {
	if (!c) {
		return GrB_SUCCESS;
	}

	if (c->reachability_count == c->reachability_cap) {
		size_t new_cap = c->reachability_cap ? c->reachability_cap * 2 : 8;
		MRReachabilityEntry *tmp_buf =
			realloc(c->reachability, new_cap * sizeof(MRReachabilityEntry));
		if (!tmp_buf) {
			return GrB_OUT_OF_MEMORY;
		}
		c->reachability = tmp_buf;
		c->reachability_cap = new_cap;
	}

	MRReachabilityEntry *e = &c->reachability[c->reachability_count++];
	e->graph_key = graph_key;
	e->grammar_tag = grammar_tag;
	e->value.matrices = reachability_matrices;
	e->value.count = symbols_amount;

	return GrB_SUCCESS;
}

GrB_Info mr_cache_insert_extracted(MRCache *c, uint64_t graph_key,
								   uint32_t grammar_tag, GrB_Matrix *edges,
								   int64_t edges_count, GrB_Matrix reachability) {
	if (!c) {
		return GrB_SUCCESS;
	}

	if (c->extracted_count == c->extracted_cap) {
		size_t new_cap = c->extracted_cap ? c->extracted_cap * 2 : 8;
		MRExtractedEntry *tmp_buf =
			realloc(c->extracted, new_cap * sizeof(MRExtractedEntry));
		if (!tmp_buf) {
			return GrB_OUT_OF_MEMORY;
		}
		c->extracted = tmp_buf;
		c->extracted_cap = new_cap;
	}

	MRExtractedEntry *e = &c->extracted[c->extracted_count++];
	e->graph_key = graph_key;
	e->grammar_tag = grammar_tag;
	e->value.edges = edges;
	e->value.edges_count = edges_count;
	e->value.reachability = reachability;

	return GrB_SUCCESS;
}
