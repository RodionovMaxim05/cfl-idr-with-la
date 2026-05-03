#pragma once

#include "grammar.h"

MRGrammar_t get_alpha_grammar(IdrGrammarType curGrammar, int64_t n_par,
							  int64_t n_bra, bool has_normal);

MRGrammar_t get_beta_grammar(IdrGrammarType curGrammar, int64_t n_par, int64_t n_bra,
							 bool has_normal);
