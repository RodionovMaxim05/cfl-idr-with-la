#pragma once

#include "GraphBLAS.h"
#include "LAGraph.h"
#include "LAGraphX.h"
#include <stdbool.h>

#define NT_S 0

typedef enum {
	DEFAULT,
	PARITY,
	PARITY2,
	SE,
	PROJECT,
	EXCLUDE,
	ALL,
	PARITYD,
	ON_DEMAND
} MRGrammarType;

typedef struct {
	int64_t nonterms_count;
	int64_t terms_count;
	int64_t rules_count;
	LAGraph_rule_WCNF *rules;
} MRGrammar_t;

MRGrammar_t dyck_grammar(int64_t n_par, int64_t n_bra, bool has_normal);

MRGrammar_t dyck_alpha_grammar(int64_t n_par, int64_t n_bra, bool has_normal);

MRGrammar_t dyck_beta_grammar(int64_t n_par, int64_t n_bra, bool has_normal);

void grammar_free(MRGrammar_t *gr);
