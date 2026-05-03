#pragma once

#include "grammar/grammar.h"

#define DEFAULT_OUTPUT_DIR "output"

typedef struct {
	const char *graph_file_path;
	const char *grammar_str;
	const char *output_path;
	IdrGrammarType grammar_type;
	int quiet;
	int valueflow;
	int on_demand;
	int parity_d;
} Args;

void print_usage(const char *prog);

IdrGrammarType grammar_name_to_type(const char *grammar);

int parse_args(int argc, char *argv[], Args *out);
