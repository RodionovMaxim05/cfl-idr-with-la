BUILD_DIR = build
TARGET = cfl-idr-with-la

.PHONY: all release run clean rebuild test memcheck test-memcheck format lint

all:
	cmake -B $(BUILD_DIR) -DCMAKE_BUILD_TYPE=Debug
	cmake --build $(BUILD_DIR) -j$(nproc)

release:
	cmake -B $(BUILD_DIR) -DCMAKE_BUILD_TYPE=Release
	cmake --build $(BUILD_DIR) -j$(nproc)

run: all
	./$(BUILD_DIR)/$(TARGET)

clean:
	rm -rf $(BUILD_DIR)

rebuild: clean all

test:
	cmake -B $(BUILD_DIR) -DCMAKE_BUILD_TYPE=Debug
	cmake --build $(BUILD_DIR) -j$(nproc)
	cd $(BUILD_DIR) && ctest --verbose --output-on-failure

memcheck:
	cmake -B $(BUILD_DIR) -DCMAKE_BUILD_TYPE=Debug
	cmake --build $(BUILD_DIR) -j$(nproc)
	valgrind \
		--leak-check=full \
		--show-leak-kinds=all \
		--track-origins=yes \
		--error-exitcode=1 \
		./$(BUILD_DIR)/$(TARGET)

test-memcheck: clean
	cmake -B $(BUILD_DIR) -DCMAKE_BUILD_TYPE=Debug
	cmake --build $(BUILD_DIR) -j$(nproc)
	cd $(BUILD_DIR) && ctest -T memcheck --verbose --output-on-failure

format:
	find . -type f \( -name "*.c" -o -name "*.h" \) \
		-not \( -path "./build/*" -o -path "./vendor/*" \) | \
	xargs -n1 sh -c \
		'clang-format --dry-run --Werror "$$1" || { echo "Format error: $$1"; exit 1; }' _

lint:
	find . -type f -name "*.c" \
		-not \( -path "./tests/*" -o -path "./build/*" -o -path "./vendor/*" \
		        -o -name "convert_graph.c" \) \
	| xargs clang-tidy -p ./build
