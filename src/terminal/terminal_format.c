#include "terminal_format.h"
#include <stdio.h>
#include <string.h>

static BracketType default_get_type(const char *label) {
	if (label == NULL || label[1] == '\0') {
		return BRACKET_TYPE_UNKNOWN;
	}
	if (label[1] == 'p') {
		return BRACKET_TYPE_PARENTHESES;
	}
	if (label[1] == 'b') {
		return BRACKET_TYPE_BRACKETS;
	}
	return BRACKET_TYPE_UNKNOWN;
}

static int default_is_open(const char *label) {
	if (label == NULL || label[0] == '\0') {
		return -1;
	}
	if (label[0] == 'o') {
		return 1;
	}
	if (label[0] == 'c') {
		return 0;
	}
	return -1;
}

static void default_extract_id(const char *label, char *out, size_t out_size) {
	const char *sep = strstr(label, "_");
	if (sep == NULL) {
		out[0] = '\0';
		return;
	}
	strncpy(out, sep + 1, out_size - 1);
	out[out_size - 1] = '\0';
}

const TerminalFormat DefaultTerminalFormat = {
	.get_type = default_get_type,
	.is_open = default_is_open,
	.extract_id = default_extract_id,
};
