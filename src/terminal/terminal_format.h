#pragma once

#include <stdbool.h>
#include <stddef.h>

typedef enum {
	BRACKET_TYPE_PARENTHESES,
	BRACKET_TYPE_BRACKETS,
	BRACKET_TYPE_UNKNOWN
} BracketType;

typedef struct {
	BracketType (*get_type)(const char *label);

	int (*is_open)(const char *label);

	void (*extract_id)(const char *label, char *out, size_t out_size);
} TerminalFormat;

extern const TerminalFormat DefaultTerminalFormat;
