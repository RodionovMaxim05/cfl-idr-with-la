#pragma once

#include "GraphBLAS.h"
#include "approximations/mr_graph.h"
#include "symbol_list.h"
#include "terminal_format.h"

GrB_Info build_mr_graph(GrB_Matrix *matrices, const SymbolList *symbol_list,
						GrB_Index n, const TerminalFormat *fmt, MRGraph *out);
