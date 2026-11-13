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
 * @brief Per-symbol CFL-reachability matrices produced by
 * `LAGraph_CFL_reachability_adv`.
 *
 * Holds the reachability matrices for all grammar symbols (nonterminals and
 * terminals), allowing repeated edge extraction queries for different target
 * vertex pairs without recomputing reachability.
 */
typedef struct {
	GrB_Matrix
		*matrices; // Array of `count` reachability matrices, one per grammar symbol
	int64_t count; // Number of symbols (length of `matrices`)
} MRReachabilityResult;

/**
 * @brief Extracted terminal edge matrices and start-symbol reachability matrix.
 *
 * Stores the final extracted terminal edge matrices (mapped to straight indices)
 * and the non-trivial boolean reachability matrix of the start nonterminal for
 * full-graph extraction.
 */
typedef struct {
	GrB_Matrix *edges;	 // Array of `edges_count` boolean edge matrices
	int64_t edges_count; // Number of edge matrices
	GrB_Matrix
		reachability; // Non-trivial reachability matrix of the start nonterminal
} MRExtractedResult;

/**
 * @brief Cache entry for per-symbol CFL-reachability results.
 */
typedef struct {
	uint64_t graph_key;			// Hash identifying the input graph structure
	uint32_t grammar_tag;		// Grammar variant tag (see `GrammarTag`)
	MRReachabilityResult value; // Cached per-symbol reachability matrices
} MRReachabilityEntry;

/**
 * @brief Cache entry for full-graph extracted terminal edges and reachability.
 */
typedef struct {
	uint64_t graph_key;		 // Hash identifying the input graph structure
	uint32_t grammar_tag;	 // Grammar variant tag (see `GrammarTag`)
	MRExtractedResult value; // Cached extracted edges and reachability result
} MRExtractedEntry;

/**
 * @brief Two-level cache storing per-symbol CFL reachability and full extracted
 * edges.
 */
typedef struct {
	MRReachabilityEntry *reachability; // Array of per-symbol reachability entries
	size_t reachability_count;		   // Number of active reachability entries
	size_t reachability_cap;		   // Allocated capacity for reachability entries

	MRExtractedEntry *extracted; // Array of full extracted edge entries
	size_t extracted_count;		 // Number of active extracted entries
	size_t extracted_cap;		 // Allocated capacity for extracted entries
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
 * For each cached reachability entry:
 * - Frees all symbol reachability matrices via `GrB_free` and releases the array.
 *
 * For each cached extracted entry:
 * - Frees all terminal edge matrices and the start-symbol reachability matrix
 *   via `GrB_free`, and releases the `edges` array.
 *
 * Finally, frees internal entry arrays and resets the cache to zero.
 *
 * @param[in,out] c Cache to free. May be `NULL` (no-op).
 */
void mr_cache_free(MRCache *c);

/**
 * @brief Looks up cached per-symbol CFL-reachability matrices.
 *
 * Performs a linear scan over reachability entries matching both `graph_key`
 * and `grammar_tag`.
 *
 * @param[in] c           Cache to search. May be `NULL`.
 * @param[in] graph_key   Hash of the input graph (from `get_graph_cache_hash`).
 * @param[in] grammar_tag Grammar variant tag (see `GrammarTag`).
 *
 * @return Pointer to `const MRReachabilityResult`, or `NULL` if not found.
 *
 * @note The returned pointer remains valid until the cache is modified or freed.
 */
const MRReachabilityResult *mr_cache_lookup_reachability(const MRCache *c,
														 uint64_t graph_key,
														 uint32_t grammar_tag);

/**
 * @brief Looks up a cached full-graph extracted edges and reachability result.
 *
 * Performs a linear scan over extracted entries matching both `graph_key`
 * and `grammar_tag`.
 *
 * @param[in] c           Cache to search. May be `NULL`.
 * @param[in] graph_key   Hash of the input graph (from `get_graph_cache_hash`).
 * @param[in] grammar_tag Grammar variant tag (see `GrammarTag`).
 *
 * @return Pointer to `const MRExtractedResult`, or `NULL` if not found.
 *
 * @note The returned pointer remains valid until the cache is modified or freed.
 */
const MRExtractedResult *mr_cache_lookup_extracted(const MRCache *c,
												   uint64_t graph_key,
												   uint32_t grammar_tag);

/**
 * @brief Inserts per-symbol CFL-reachability matrices into the cache.
 *
 * Takes ownership of `reachability_matrices` and all matrices contained within it.
 * Dynamically resizes the `reachability` entry array when capacity is exceeded.
 *
 * @param[in,out] c                     Cache to insert into. May be `NULL` (no-op).
 * @param[in]     graph_key             Hash of the input graph.
 * @param[in]     grammar_tag           Grammar variant tag.
 * @param[in]     reachability_matrices Array of per-symbol reachability matrices.
 * @param[in]     symbols_amount        Number of symbols in `reachability_matrices`.
 *
 * @return `GrB_SUCCESS` on success, or `GrB_OUT_OF_MEMORY` if reallocation fails.
 *
 * @note On `GrB_SUCCESS`, the caller must not free `reachability_matrices`; it will
 *       be freed automatically when `mr_cache_free` is called.
 */
GrB_Info mr_cache_insert_reachability(MRCache *c, uint64_t graph_key,
									  uint32_t grammar_tag,
									  GrB_Matrix *reachability_matrices,
									  int64_t symbols_amount);

/**
 * @brief Inserts extracted terminal edge matrices and reachability matrix into the
 * cache.
 *
 * Takes ownership of the `edges` array (and contained matrices) and `reachability`.
 * Dynamically resizes the `extracted` entry array when capacity is exceeded.
 *
 * @param[in,out] c            Cache instance. May be `NULL` (no-op).
 * @param[in]     graph_key    Hash of the input graph.
 * @param[in]     grammar_tag  Grammar variant tag.
 * @param[in]     edges        Array of extracted terminal edge matrices.
 * @param[in]     edges_count  Number of edge matrices in `edges`.
 * @param[in]     reachability Start-symbol reachability matrix.
 *
 * @return `GrB_SUCCESS` on success, or `GrB_OUT_OF_MEMORY` on allocation failure.
 */
GrB_Info mr_cache_insert_extracted(MRCache *c, uint64_t graph_key,
								   uint32_t grammar_tag, GrB_Matrix *edges,
								   int64_t edges_count, GrB_Matrix reachability);
