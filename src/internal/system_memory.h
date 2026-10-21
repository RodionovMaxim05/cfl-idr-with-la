#pragma once

/**
 * @brief Global threshold for system memory critical check (in percent).
 *
 * TODO: Consider refactoring this to pass it as a parameter.
 */
#define RAM_BUDGET_PERCENT 90

/**
 * @brief Checks whether the memory in use exceeds the configured budget.
 *
 * "In use" is measured as the complement of the memory that the OS can still hand
 * out without swapping (free memory plus reclaimable caches). On Linux the limit of
 * the enclosing cgroup (Docker, systemd) is honoured when it is lower than the
 * physical RAM.
 *
 * Supported platforms: Linux, macOS and Windows. On any other platform the function
 * always returns -1.
 *
 * @param[in] ram_budget_percent Maximum percent of RAM the program is allowed to
 *                               have in use (1..100).
 *
 * @return 1 if memory use is above the budget, 0 if it is within the budget,
 *         -1 if the platform is unsupported, the OS memory information cannot
 *         be read, or `ram_budget_percent` is out of range.
 */
int is_system_memory_critical(unsigned int ram_budget_percent);

/**
 * @brief Asks the allocator to return unused heap memory to the OS.
 *
 * Freeing memory with `free()` does not necessarily lower the system's memory usage.
 * Currently implemented for glibc (`malloc_trim`); on other platforms it does
 * nothing.
 */
void release_unused_memory_to_os(void);
