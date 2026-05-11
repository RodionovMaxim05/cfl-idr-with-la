#pragma once

#include <GraphBLAS.h>
#include <stdint.h>

#include "cfl_idr.h"

/**
 * @brief Tags identifying the grammar variant used in mutual-refinement caching.
 *
 * Each tag corresponds to a specific Dyck-CFL grammar configuration. The cache
 * uses `(graph_key, grammar_tag)` as a composite key to avoid recomputing
 * expensive CFL-reachability results for identical inputs.
 */
typedef enum {
	GRAMMAR_TAG_ALPHA = 0,			// Alpha grammar tag
	GRAMMAR_TAG_BETA = 1,			// Beta grammar tag
	GRAMMAR_TAG_PROJECT = 2,		// Project grammar tag
	GRAMMAR_TAG_EXCLUDE_BASE = 100, // Base tag for exclusion grammars
} GrammarTag;

/**
 * @brief Result of a single mutual-refinement CFL-reachability step.
 *
 * Holds the output matrices for all nonterminals produced by `LAGraph_CFL_AllPaths`,
 * along with metadata required for proper cleanup and reuse.
 */
typedef struct {
	GrB_Matrix *matrices; // Array of `count` matrices, one per grammar nonterminal
	int64_t count;		  // Number of nonterminals (length of `matrices`)
	GrB_Type all_paths_t; // GraphBLAS type used for path values
} MRStepResult;

/**
 * @brief Single entry in the mutual-refinement result cache.
 *
 * Maps a `(graph_key, grammar_tag)` pair to a precomputed `MRStepResult`.
 * Used to avoid redundant CFL-reachability computations during iterative
 * mutual-refinement passes.
 */
typedef struct {
	uint64_t graph_key;	  // Hash identifying the input graph structure
	uint32_t grammar_tag; // Grammar variant tag (see `GrammarTag`)
	MRStepResult value;	  // Cached CFL-reachability result
} MRCacheEntry;

/**
 * @brief Simple linear-probe cache for mutual-refinement intermediate results.
 *
 * Stores `MRCacheEntry` records in a dynamically-resized array.
 */
typedef struct {
	MRCacheEntry *entries; // Array of cached entries
	size_t count;		   // Number of valid entries currently stored
	size_t capacity;	   // Allocated capacity of `entries` array
} MRCache;

/**
 * @brief Computes a deterministic hash key for an `IdrGraph` to enable cache lookup.
 *
 * Uses the FNV-1a 64-bit hash algorithm over:
 * - Graph size (`n`)
 * - Structural fingerprints of each adjacency matrix (parenthesis/bracket pairs
 *   and optional `normal` matrix), computed from sorted edge tuples `(row, col)`.
 *
 * @param[in] graph  Input graph to hash. Must not be `NULL`.
 *
 * @return 64-bit hash value suitable for use as `graph_key` in cache operations.
 */
uint64_t get_graph_cache_hash(const IdrGraph *graph);

/**
 * @brief Initializes an `MRCache` structure to an empty state.
 *
 * Sets all fields to zero; no memory is allocated at this stage.
 *
 * @param[out] c  Cache structure to initialize. Must not be `NULL`.
 *
 * @note This function does not allocate heap memory. Allocation occurs on
 *       first `mr_cache_insert` call.
 */
void mr_cache_init(MRCache *c);

/**
 * @brief Releases all resources owned by an `MRCache`.
 *
 * For each cached entry:
 * - Frees CFL-reachability output matrices via `LAGraph_CFL_AllPaths_free_outputs`
 * - Frees the `all_paths_t` type object via `GrB_free`
 *
 * Finally, frees the internal `entries` array and resets the cache to zero.
 *
 * @param[in] c  Cache to free. May be `NULL` (no-op).
 */
void mr_cache_free(MRCache *c);

/**
 * @brief Looks up a cached result by graph hash and grammar tag.
 *
 * Performs a linear scan of the cache entries to find an exact match for both
 * `graph_key` and `grammar_tag`. Returns a pointer to the stored `MRStepResult`
 * if found; otherwise returns `NULL`.
 *
 * @param[in] c           Cache to search. Must not be `NULL`.
 * @param[in] graph_key   Hash of the input graph (from `get_graph_cache_hash`).
 * @param[in] grammar_tag Grammar variant tag (see `GrammarTag`).
 *
 * @return Pointer to the cached `MRStepResult`, or `NULL` if not found.
 *
 * @note The returned pointer remains valid until the cache is modified or freed.
 *       Do not modify the contents of the returned result.
 */
const MRStepResult *mr_cache_lookup(const MRCache *c, uint64_t graph_key,
									uint32_t grammar_tag);

/**
 * @brief Inserts a new CFL-reachability result into the cache.
 *
 * Takes ownership of the provided `paths_matrices` array and `all_paths_t` type.
 * The cache dynamically resizes its internal storage (doubling strategy) when
 * capacity is exceeded.
 *
 * @param[in,out] c              Cache to insert into. Must not be `NULL`.
 * @param[in]     graph_key      Hash of the input graph.
 * @param[in]     grammar_tag    Grammar variant tag.
 * @param[in]     paths_matrices Array of `nonterms_count` result matrices.
 * @param[in]     nonterms_count Number of nonterminals (length of `paths_matrices`).
 * @param[in]     all_paths_t    GraphBLAS type for path values.
 *
 * @return `GrB_SUCCESS` on success, or `GrB_OUT_OF_MEMORY` if reallocation fails.
 *
 * @note The caller must not free `paths_matrices` or `all_paths_t` after insertion;
 *       they will be freed automatically when the cache entry is evicted or
 *       `mr_cache_free` is called.
 */
GrB_Info mr_cache_insert(MRCache *c, uint64_t graph_key, uint32_t grammar_tag,
						 GrB_Matrix *paths_matrices, int64_t nonterms_count,
						 GrB_Type all_paths_t);
