#include "cli.h"

#include <stdio.h>

void print_usage(const char *prog) {
	fprintf(stderr,
			"Usage: %s [options] <input.g> <grammar>\n"
			"\n"
			"Arguments:\n"
			"  input.g      Path to the input graph file\n"
			"  grammar      Grammar type (see below)\n"
			"\n"
			"Grammars:\n"
			"  parity       Parity grammar (k=1)\n"
			"  parity2      Parity grammar (k=2)\n"
			"  se           Structured equality grammar\n"
			"  project      Projection grammar\n"
			"  exclude      Exclusion grammar\n"
			"  all          Comprehensive grammar\n"
			"  parityD      On-demand parity grammar\n"
			"  on-demand    On-demand combined method\n"
			"\n"
			"Options:\n"
			"  -o <path>    Output directory (default: " DEFAULT_OUTPUT_DIR ")\n"
			"  -q           Quiet mode — suppress pair output\n"
			"  -valueflow   Enable value-flow specific optimizations and "
			"constraints (`{s | s = [i ∗ ]i}`)\n"
			"  -h           Show this help message\n",
			prog);
}

MRGrammarType grammar_name_to_type(const char *grammar) {
	if (strcmp(grammar, "default") == 0) {
		return DEFAULT;
	}
	if (strcmp(grammar, "parity") == 0) {
		return PARITY;
	}
	if (strcmp(grammar, "parity2") == 0) {
		return PARITY2;
	}
	if (strcmp(grammar, "se") == 0) {
		return SE;
	}
	if (strcmp(grammar, "project") == 0) {
		return PROJECT;
	}
	if (strcmp(grammar, "exclude") == 0) {
		return EXCLUDE;
	}
	if (strcmp(grammar, "all") == 0) {
		return ALL;
	}
	if (strcmp(grammar, "parityD") == 0) {
		return PARITY;
	}
	if (strcmp(grammar, "on-demand") == 0) {
		return ALL;
	}
	return UNKNOWN;
}

int parse_args(int argc, char *argv[], Args *out) {
	const char *grammar_str = NULL;
	const char *graph_file_path = NULL;
	const char *output_path = NULL;
	int quiet = 0;
	int valueflow = 0;

	for (int i = 1; i < argc; i++) {
		if (strcmp(argv[i], "-h") == 0) {
			print_usage(argv[0]);
			exit(EXIT_SUCCESS);
		} else if (strcmp(argv[i], "-o") == 0) {
			if (i + 1 >= argc) {
				fprintf(stderr, "Error: Missing value for -o\n");
				return 0;
			}
			output_path = argv[++i];
		} else if (strcmp(argv[i], "-q") == 0) {
			quiet = 1;
		} else if (strcmp(argv[i], "-valueflow") == 0) {
			valueflow = 1;
		} else if (graph_file_path == NULL) {
			graph_file_path = argv[i];
		} else if (grammar_str == NULL) {
			grammar_str = argv[i];
		} else {
			fprintf(stderr, "Error: Unexpected extra argument: %s\n", argv[i]);
			return 0;
		}
	}

	if (graph_file_path == NULL) {
		fprintf(stderr, "Error: Missing graph file path\n");
		return 0;
	}
	if (grammar_str == NULL) {
		fprintf(stderr, "Error: Missing grammar\n");
		return 0;
	}

	int on_demand = 0;
	int parity_d = 0;
	if (strcmp(grammar_str, "parityD") == 0) {
		grammar_str = "parity";
		parity_d = 1;
	} else if (strcmp(grammar_str, "on-demand") == 0) {
		grammar_str = "all";
		on_demand = 1;
	}

	MRGrammarType grammar_type = grammar_name_to_type(grammar_str);
	if (grammar_type == UNKNOWN) {
		fprintf(stderr, "Error: Invalid grammar '%s'\n", grammar_str);
		fprintf(stderr, "Run with -h to see supported grammars\n");
		return 0;
	}

	out->graph_file_path = graph_file_path;
	out->grammar_str = grammar_str;
	out->output_path = output_path;
	out->grammar_type = grammar_type;
	out->quiet = quiet;
	out->valueflow = valueflow;
	out->on_demand = on_demand;
	out->parity_d = parity_d;
	return 1;
}
