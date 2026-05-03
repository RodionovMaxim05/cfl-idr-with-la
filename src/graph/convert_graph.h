#pragma once

#include <GraphBLAS.h>

#include "cfl_idr.h"
#include "parser.h"
#include "symbol_list.h"
#include "terminal/terminal_format.h"

GrB_Info build_idr_graph(const GraphMatrices *gm, const SymbolList *symbol_list,
						 GrB_Index n, const TerminalFormat *fmt, IdrGraph *out);
