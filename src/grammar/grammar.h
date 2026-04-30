#pragma once

#include "GraphBLAS.h"
#include "LAGraph.h"
#include "LAGraphX.h"
#include <stdbool.h>

#define NT_START 0

typedef enum {
	UNKNOWN = -1,
	DEFAULT,
	PARITY,
	PARITY2,
	SE,
	PROJECT,
	EXCLUDE,
	ALL
} MRGrammarType;

typedef struct {
	int64_t nonterms_count;
	int64_t terms_count;
	int64_t rules_count;
	LAGraph_rule_WCNF *rules;
} MRGrammar_t;

#define RSTATE_QE 0
#define RSTATE_QO 1
#define RSTATE_QC 2
#define RSTATE_COUNT (int64_t)3

MRGrammar_t dyck_grammar(int64_t n_par, int64_t n_bra, bool has_normal);

MRGrammar_t dyck_project_grammar(int64_t n_par, int64_t n_bra, bool has_normal);

MRGrammar_t dyck_alpha_grammar(int64_t n_par, int64_t n_bra, bool has_normal);

MRGrammar_t dyck_alpha_grammar_k_parity(int64_t n_par, int64_t n_bra,
										bool has_normal, int64_t k);

MRGrammar_t dyck_alpha_grammar_k_parity_exclude(int64_t n_par, int64_t n_bra,
												bool has_normal, int64_t k,
												int64_t ex_bra);

MRGrammar_t dyck_beta_grammar(int64_t n_par, int64_t n_bra, bool has_normal);

MRGrammar_t dyck_beta_grammar_k_parity(int64_t n_par, int64_t n_bra, bool has_normal,
									   int64_t k);

MRGrammar_t dyck_alpha_grammar_k_parity_se(int64_t n_par, int64_t n_bra,
										   bool has_normal, int64_t k);

MRGrammar_t dyck_beta_grammar_k_parity_se(int64_t n_par, int64_t n_bra,
										  bool has_normal, int64_t k);

void grammar_free(MRGrammar_t *gr);
