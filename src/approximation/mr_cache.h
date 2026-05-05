#pragma once

#include <GraphBLAS.h>
#include <stdint.h>

#include "cfl_idr.h"

typedef enum {
	GRAMMAR_TAG_ALPHA = 0,
	GRAMMAR_TAG_BETA = 1,
	GRAMMAR_TAG_PROJECT = 2,
	GRAMMAR_TAG_EXCLUDE_BASE = 100,
} GrammarTag;

typedef struct {
	GrB_Matrix *matrices;
	int64_t count;
	GrB_Type all_paths_t;
} MRStepResult;

typedef struct {
	uint64_t graph_key;
	uint32_t grammar_tag;
	MRStepResult value;
} MRCacheEntry;

typedef struct {
	MRCacheEntry *entries;
	size_t count;
	size_t capacity;
} MRCache;

uint64_t get_graph_cache_hash(const IdrGraph *graph);

void mr_cache_init(MRCache *c);

void mr_cache_free(MRCache *c);

const MRStepResult *mr_cache_lookup(const MRCache *c, uint64_t graph_key,
									uint32_t grammar_tag);

GrB_Info mr_cache_insert(MRCache *c, uint64_t graph_key, uint32_t grammar_tag,
						 GrB_Matrix *paths_matrices, int64_t nonterms_count,
						 GrB_Type all_paths_t);
