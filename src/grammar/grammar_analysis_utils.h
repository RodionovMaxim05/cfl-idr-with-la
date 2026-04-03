#pragma once

#include "grammar.h"

MRGrammar_t get_alpha_grammar(MRGrammarType curGrammar, int n_par, int n_bra,
							  bool has_normal);

MRGrammar_t get_beta_grammar(MRGrammarType curGrammar, int n_par, int n_bra,
							 bool has_normal);
