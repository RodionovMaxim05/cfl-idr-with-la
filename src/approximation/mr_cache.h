#pragma once

#include <GraphBLAS.h>
#include <LAGraphX.h>
#include <stdint.h>

#include "cfl_idr.h"

/**
 * @brief Tags identifying the grammar variant used in mutual-refinement caching.
 *
 * Each tag corresponds to a specific Dyck-CFL grammar configuration. The cache
 * uses `(graph_key, grammar_tag)` as a composite key to avoid recomputing
 * expensive CFL-reachability and edge extraction results for identical inputs.
 */
typedef enum {
	GRAMMAR_TAG_ALPHA = 0,			// Alpha grammar tag
	GRAMMAR_TAG_BETA = 1,			// Beta grammar tag
	GRAMMAR_TAG_PROJECT = 2,		// Project grammar tag
	GRAMMAR_TAG_EXCLUDE_BASE = 100, // Base tag for exclusion grammars
} GrammarTag;

/**
 * @brief Refined edge matrices and reachability matrix resulting from a CFL step.
 *
 * Stores extracted terminal edge matrices and precomputed boolean reachability
 * matrix.
 */
typedef struct {
	GrB_Matrix *edges;		 // Array of `edges_count` boolean edge matrices
	int64_t edges_count;	 // Number of edge matrices
	GrB_Matrix reachability; // Precomputed reachability matrix
} MRRefinedResult;

/**
 * @brief Cache entry storing refined terminal edges and reachability.
 */
typedef struct {
	uint64_t graph_key;	   // Hash identifying the input graph structure
	uint32_t grammar_tag;  // Grammar variant tag (see `GrammarTag`)
	MRRefinedResult value; // Cached refined edges result
} MREntry;

/**
 * @brief Cache storing refined edge results per (graph_key, grammar_tag).
 */
typedef struct {
	MREntry *entries; // Array of cached entries
	size_t count;	  // Number of active entries
	size_t cap;		  // Allocated capacity for entries
} MRCache;

/**
 * @brief Computes a deterministic hash key for an `IdrGraph` to enable cache lookup.
 *
 * Uses the FNV-1a 64-bit hash algorithm over:
 * - Graph size (`n`)
 * - Structural fingerprints of each adjacency matrix (parenthesis/bracket pairs
 *   and optional `normal` matrix), computed from sorted edge tuples `(row, col)`.
 *
 * @param[in] graph Input graph to hash. Must not be `NULL`.
 *
 * @return 64-bit hash value suitable for use as `graph_key` in cache operations.
 */
uint64_t get_graph_cache_hash(const IdrGraph *graph);

/**
 * @brief Initializes an `MRCache` structure to an empty state.
 *
 * Sets all fields to zero; no memory is allocated at this stage.
 *
 * @param[out] c Cache structure to initialize. Must not be `NULL`.
 */
void mr_cache_init(MRCache *c);

/**
 * @brief Releases all resources owned by an `MRCache`.
 *
 * For each cached entry, frees all boolean edge matrices and the reachability
 * matrix via `GrB_free`, frees internal arrays, and resets the cache to zero.
 *
 * @param[in,out] c Cache to free. May be `NULL` (no-op).
 */
void mr_cache_free(MRCache *c);

/**
 * @brief Looks up a cached refined edges and reachability result.
 *
 * Performs a linear scan over entries matching both `graph_key` and `grammar_tag`.
 *
 * @param[in] c           Cache to search. May be `NULL`.
 * @param[in] graph_key   Hash of the input graph (from `get_graph_cache_hash`).
 * @param[in] grammar_tag Grammar variant tag (see `GrammarTag`).
 *
 * @return Pointer to `const MRRefinedResult`, or `NULL` if not found.
 *
 * @note The returned pointer remains valid until the cache is modified or freed.
 */
const MRRefinedResult *mr_cache_lookup_refined(const MRCache *c, uint64_t graph_key,
											   uint32_t grammar_tag);

/**
 * @brief Inserts a refined edges result and reachability matrix into the cache.
 *
 * Takes ownership of `edges` array (and contained matrices) and `reachability`.
 * Dynamically resizes the internal entry array when capacity is exceeded.
 *
 * @param[in,out] c            Cache instance. May be `NULL` (no-op).
 * @param[in]     graph_key    Hash of the input graph.
 * @param[in]     grammar_tag  Grammar variant tag.
 * @param[in]     edges        Array of extracted terminal edge matrices.
 * @param[in]     edges_count  Number of edge matrices in `edges`.
 * @param[in]     reachability Reachability matrix.
 *
 * @return `GrB_SUCCESS` on success, or `GrB_OUT_OF_MEMORY` on allocation failure.
 */
GrB_Info mr_cache_insert_refined(MRCache *c, uint64_t graph_key,
								 uint32_t grammar_tag, GrB_Matrix *edges,
								 int64_t edges_count, GrB_Matrix reachability);
