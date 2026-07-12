# cfl-idr-with-la

A C implementation of CFL-based approximation methods for interleaved Dyck reachability using GraphBLAS and LAGraph.

## Problem Overview

Interleaved Dyck Reachability is a fundamental problem in static program analysis where path constraints are expressed as the shuffle of two Dyck languages — typically representing context sensitivity (parentheses) and field sensitivity (brackets). While this formulation enables high-precision analyses, the problem is known to be undecidable in general, and thus existing approaches resort to clever overapproximations.

This project implements approximation techniques based on the research paper ["CFL-based methods for approximating interleaved Dyck reachability" (Conrado & Pavlogiannis, 2025)](https://www.researchgate.net/publication/390804794_CFL-based_methods_for_approximating_interleaved_Dyck_reachability), which transforms the undecidable problem into a series of tractable CFL reachability problems solvable via the [GraphBLAS-based algorithm from the LAGraph library](https://github.com/SparseLinearAlgebra/LAGraph/blob/homka122/all_algorithms_benchmark/include/LAGraphX.h#L1117).

### Supported Grammars

| Grammar Type | Description                                                                                                                                                                       | Article designation |
|--------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------------------|
| `parity`     | Mutual refinement with k=1 parity condition on parentheses/brackets: requires the total count of parentheses/bracket symbols to be even                                           | `PAR`               |
| `parity2`    | Mutual refinement with k=2 parity condition: parentheses/brackets are split into 2 groups, count in each group must be even                                                       | `PAR2`              |
| `se`         | Extends `PAR2` with a "valid endpoints" (structured equality) condition: parentheses/bracket substrings must start with `(`/`[` and end with `)`/`]`                              | `PAR2E`             |
| `project`    | Extends `PAR` with unlabeled projection: additionally runs reachability on the graph with all labels stripped, catching unreachable pairs where total nesting depth is unbalanced | `PARUnl`            |
| `exclude`    | Extends `PAR` with label erasure: iterates over each bracket label, erases it from the graph, and reruns `PAR` - exposing unreachable pairs that were hidden by that label        | `PARErase`          |
| `all`        | Combines `PAR2E` grammar with both the unlabeled projection phase and the label-erasure phase in the mutual refinement loop                                                       | `COM`               |
| `parityD`    | Same grammar as `PAR`, but run in on-demand mode: mutual refinement is performed independently for each queried node pair, filtering out irrelevant graph parts                   | `PARD`              |
| `on-demand`  | Same grammar as `COM`, but run in on-demand mode: combines all stronger grammars and graph simplifications with per-pair refinement                                               | `COMD`              |

## Architecture

The project is organized as a modular C library with a command-line interface:

### Programmatic API

All functionality is accessible via [include/cfl_idr.h](./include/cfl_idr.h), which provides:

- `IdrGraph`: Sparse-matrix representation of a Dyck-reachability graph with parenthesis/bracket-labeled edges
- `IdrGrammarType`: Enum controlling the precision/cost trade-off in over-approximation
- Core analysis functions:
  - `idr_get_under_approx()` — conservative under-approximation of reachable paths
  - `idr_get_over_approx()` — refined over-approximation via mutual refinement
  - `idr_get_on_demand()` — efficient per-pair refinement for incremental queries
- Utility functions:
  - `idr_remove_valueflow_unreachable()` — prune vertices irrelevant to value-flow analysis
  - `idr_graph_free()` — release all resources owned by an IdrGraph

### Documentation

Full API reference generated with Doxygen from the public headers in [`include/`](./include) is available **[here](https://rodionovmaxim05.github.io/cfl-idr-with-la/)**.

### Command-Line Interface

```bash
./build/cfl-idr-with-la [options] <input.g> <grammar>
```

## Value-Flow Analysis Mode

The library provides a specialized value-flow analysis mode (`-valueflow` flag) optimized for memory value-flow analysis: a valid path must start with `[` (similar to the store operation) and end with `]` (similar to the load operation) with the same id

`Recommendation`: For manual preprocessing before analysis:

```C
// Remove vertices unreachable from store/load operations
IdrGraph filtered_graph;
GrB_Info info = idr_remove_valueflow_unreachable(&filtered_graph, &original_graph);
```

## Quick Start

### Requirements

- C11-compatible compiler (developed and tested with GCC 13.3)
- [GraphBLAS library](https://github.com/DrTimothyAldenDavis/GraphBLAS) installed
- [LAGraph library](https://github.com/SparseLinearAlgebra/LAGraph/tree/homka122/all_algorithms_benchmark) is cloned and LAGRAPH_DIR variable points to the library's root directory (**note that `all_algorithms_benchmark` branch is required**)
- CMake 3.20+
- Optional: valgrind for memory checking, clang-format, clang-tidy for code quality

### Building from Source

Before building, you need to specify the path to your LAGraph directory. You can do this in two ways:

1) Set it directly in the Makefile:  
    Open the [Makefile](./Makefile) and define the variable at the top:

    ```makefile
    LAGRAPH_DIR=/path/to/lagraph
    ```

    After setting it, you can simply run `make`, `make release`, etc.

2) Pass it via command line:

    ````bash
    make LAGRAPH_DIR=/path/to/lagraph
    ````

Once configured, proceed with the build:

```bash
# Build debug version
make

# Or build release version
make release

# Build only the static library (no executable)
make lib
```

### Running Tests

```bash
# Run unit tests (requires CTest)
make test

# Run programm with Valgrind memory checking
make test-memcheck ARGS="path/to/test.g parity"

# Run tests with Valgrind memory checking (takes a very long time)
make test-memcheck
```

### Code Quality

```bash
# Check code formatting (clang-format)
make format

# Run static analysis (clang-tidy)
make lint
```

## Usage Examples

### CLI

```bash
# Basic analysis with parity grammar
./build/cfl-idr-with-la graph.g parity

# Analysis with structured equality grammar and custom output
./build/cfl-idr-with-la -o results graph.g se

# Quiet mode (suppress pair output)
./build/cfl-idr-with-la -q graph.g parity2 -o output

# Value-flow analysis with exclusion grammar
./build/cfl-idr-with-la -valueflow graph.g exclude

# On-demand parity analysis (PARD mode)
./build/cfl-idr-with-la graph.g parityD

# Comprehensive on-demand analysis (COMD mode)
./build/cfl-idr-with-la graph.g on-demand
```

### API Integration

```C
#include "cfl_idr.h"
#include "graph/parse_utils.h"

// Load graph from file
IdrGraph graph;
GrB_Info info = parse_graph(&graph, "graph.g");
if (info != GrB_SUCCESS) { /* handle error */ }

// Compute under-approximation
GrB_Matrix under;
info = idr_get_under_approx(&under, &graph, /*valueflow=*/false);

// Compute over-approximation with mutual refinement
GrB_Matrix over;
info = idr_get_over_approx(&over, &graph, IDR_PARITY2, under, /*valueflow=*/false, /*filter_empty=*/true);

// On-demand refinement for specific query
GrB_Matrix refined;
info = idr_get_on_demand(&refined, &graph, under, over, /*parityD=*/false, /*valueflow=*/false, /*filter_empty=*/true);

// Cleanup
GrB_Matrix_free(&under);
GrB_Matrix_free(&over);
GrB_Matrix_free(&refined);
idr_graph_free(&graph);
```

## Benchmarks

Benchmarks for this library are maintained in a separate repository: [CFG_bench](https://github.com/RodionovMaxim05/CFG_bench/tree/cfl-idr).

Benchmark results, along with a comparison against other CFL-reachability implementations, are available **[here](https://github.com/RodionovMaxim05/cfl-idr-with-gll#performance-analysis)**.

## License

Distributed under the [MIT License](https://choosealicense.com/licenses/mit/). See [`LICENSE`](LICENSE) for more information.
