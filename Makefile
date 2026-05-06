BUILD_DIR = build
TARGET = cfl-idr-with-la

ifndef LAGRAPH_DIR
$(error LAGRAPH_DIR is not set. Usage: make LAGRAPH_DIR=/path/to/lagraph)
endif

.PHONY: all release lib run clean rebuild test memcheck test-memcheck format lint

all: $(BUILD_DIR)/CMakeCache.txt
	cmake --build $(BUILD_DIR) -j$(nproc)

$(BUILD_DIR)/CMakeCache.txt:
	cmake -B $(BUILD_DIR) -DCMAKE_BUILD_TYPE=Debug -DLAGRAPH_DIR=$(LAGRAPH_DIR)

release:
	cmake -B $(BUILD_DIR) -DCMAKE_BUILD_TYPE=Release -DLAGRAPH_DIR=$(LAGRAPH_DIR)
	cmake --build $(BUILD_DIR) -j$(nproc)

lib:
	cmake -B $(BUILD_DIR) -DCMAKE_BUILD_TYPE=Release -DLAGRAPH_DIR=$(LAGRAPH_DIR) -DBUILD_EXECUTABLE=OFF
	cmake --build $(BUILD_DIR) --target cfl_lib -j$(nproc)

run: all
	./$(BUILD_DIR)/$(TARGET) $(ARGS)

clean:
	rm -rf $(BUILD_DIR)

rebuild: clean all

test: all
	cd $(BUILD_DIR) && ctest --verbose --output-on-failure -L cfl-idr

memcheck: all
	valgrind \
		--leak-check=full \
		--show-leak-kinds=all \
		--track-origins=yes \
		--error-exitcode=1 \
		./$(BUILD_DIR)/$(TARGET) $(ARGS)

test-memcheck: clean all
	cd $(BUILD_DIR) && ctest -T memcheck --verbose --output-on-failure -L cfl-idr

format:
	find . -type f \( -name "*.c" -o -name "*.h" \) \
		-not \( -path "./build/*" -o -path "./vendor/*" \) | \
	xargs -n1 sh -c \
		'clang-format --dry-run --Werror "$$1" || { echo "Format error: $$1"; exit 1; }' _

lint:
	find . -type f -name "*.c" \
		-not \( -path "./tests/*" -o -path "./build/*" -o -path "./vendor/*" \
		        -o -name "convert_graph.c" -o -name "extract_edges.c" \) \
	| xargs clang-tidy -p ./build
