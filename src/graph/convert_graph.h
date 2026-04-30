#pragma once

#include <GraphBLAS.h>

#include "approximation/mr_graph.h"
#include "parser.h"
#include "symbol_list.h"
#include "terminal/terminal_format.h"

GrB_Info build_mr_graph(const GraphMatrices *gm, const SymbolList *symbol_list,
						GrB_Index n, const TerminalFormat *fmt, MRGraph *out);
