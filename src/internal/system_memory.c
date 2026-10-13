#include "system_memory.h"

#include <ctype.h>
#include <errno.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(__linux__) && defined(__GLIBC__)
#include <malloc.h>
#endif

#if defined(_WIN32)
#include <windows.h>
#elif defined(__APPLE__) && defined(__MACH__)
#include <mach/mach.h>
#include <mach/mach_host.h>
#include <mach/vm_statistics.h>
#include <sys/sysctl.h>
#endif

#if defined(__linux__)

/**
 * @brief Reads a single unsigned integer from the first line of a file.
 *
 * Files whose content is not a number (for example, cgroup v2 `memory.max` contains
 * the literal string "max" when no limit is set) are treated as unreadable.
 *
 * @param[in]  path  Path of the file to read.
 * @param[out] out   Receives the parsed value on success. Left untouched on failure.
 *
 * @return 0 on success, -1 if the file cannot be opened or doesn`t start with a
 *         number.
 */
static int read_ull_file(const char *path, unsigned long long *out) {
	FILE *f = fopen(path, "r");
	if (f == NULL) {
		return -1;
	}

	char buf[64];
	int ok = -1;
	if (fgets(buf, sizeof(buf), f) != NULL && isdigit((unsigned char)buf[0])) {
		char *end = NULL;
		errno = 0;
		unsigned long long v = strtoull(buf, &end, 10);
		if (errno == 0 && end != buf) {
			*out = v;
			ok = 0;
		}
	}
	fclose(f);
	return ok;
}

/**
 * @brief Looks up `key <value>` in a "key value" style file (e.g. `memory.stat`).
 *
 * @param[in]  path  Path of the file to read.
 * @param[in]  key   Exact key to look for (without trailing space or colon).
 * @param[out] out   Receives the value of the first matching line. Left untouched on
 *                   failure.
 *
 * @return 0 on success, -1 if the file cannot be opened, the key is missing, or its
 *         value cannot be parsed.
 */
static int read_stat_value(const char *path, const char *key,
						   unsigned long long *out) {
	FILE *f = fopen(path, "r");
	if (f == NULL) {
		return -1;
	}

	size_t key_len = strlen(key);
	char line[256];
	int ok = -1;
	while (fgets(line, sizeof(line), f) != NULL) {
		if (strncmp(line, key, key_len) == 0 && line[key_len] == ' ') {
			if (sscanf(line + key_len, " %llu", out) == 1) {
				ok = 0;
			}
			break;
		}
	}
	fclose(f);
	return ok;
}

/**
 * @brief Reads host-wide memory numbers from `/proc/meminfo`.
 *
 * Requires the `MemAvailable` field (Linux 3.14 and newer).
 *
 * @param[out] total_bytes      Total physical RAM, in bytes.
 * @param[out] available_bytes  RAM available without swapping, in bytes. Never
 *                              exceeds `total_bytes`.
 *
 * @return 0 on success, -1 if `/proc/meminfo` cannot be read or lacks
 *         `MemTotal` or `MemAvailable`.
 */
static int read_linux_host_memory(unsigned long long *total_bytes,
								  unsigned long long *available_bytes) {
	FILE *f = fopen("/proc/meminfo", "r");
	if (f == NULL) {
		return -1;
	}

	unsigned long long total_kb = 0, available_kb = 0;
	bool have_total = false, have_available = false;
	char line[256];

	while (fgets(line, sizeof(line), f) != NULL) {
		if (strncmp(line, "MemTotal:", 9) == 0) {
			have_total = sscanf(line + 9, " %llu", &total_kb) == 1;
		} else if (strncmp(line, "MemAvailable:", 13) == 0) {
			have_available = sscanf(line + 13, " %llu", &available_kb) == 1;
			break;
		}
	}
	fclose(f);

	if (!have_total || !have_available || total_kb == 0) {
		return -1;
	}
	if (available_kb > total_kb) {
		available_kb = total_kb;
	}

	*total_bytes = total_kb * 1024ULL;
	*available_bytes = available_kb * 1024ULL;
	return 0;
}

/**
 * @brief Builds the sysfs directory of the current process's cgroup v2 group.
 *
 * Inside a container with its own cgroup namespace the group path is "/" and the
 * directory is `/sys/fs/cgroup`; otherwise the path reported by `/proc/self/cgroup`
 * is appended to it.
 *
 * @param[out] buf   Buffer receiving the NUL-terminated directory path.
 * @param[in]  size  Size of `buf`, in bytes.
 *
 * @return 0 on success, -1 if the process is not in a cgroup v2 hierarchy or the
 *         path does not fit into `buf`.
 */
static int cgroup2_dir(char *buf, size_t size) {
	FILE *f = fopen("/proc/self/cgroup", "r");
	if (f == NULL) {
		return -1;
	}

	char line[512];
	int ok = -1;
	while (fgets(line, sizeof(line), f) != NULL) {
		if (strncmp(line, "0::", 3) == 0) {
			line[strcspn(line, "\r\n")] = '\0';
			const char *path = line + 3;
			int written = (strcmp(path, "/") == 0)
							  ? snprintf(buf, size, "/sys/fs/cgroup")
							  : snprintf(buf, size, "/sys/fs/cgroup%s", path);
			ok = (written > 0 && (size_t)written < size) ? 0 : -1;
			break;
		}
	}
	fclose(f);
	return ok;
}

/**
 * @brief Reads the memory limit, current usage and reclaimable page cache of
 * the current cgroup.
 *
 * Tries cgroup v2 first and falls back to cgroup v1.
 *
 * @param[out] limit          Memory limit, in bytes. In cgroup v1 an unlimited
 *                            group reports a very large number instead of an
 *                            error.
 * @param[out] current        Current memory usage, in bytes.
 * @param[out] inactive_file  Inactive file-backed pages, in bytes. These can be
 *                            dropped by the kernel and are not counted as used.
 *                            Set to 0 if the statistic is unavailable.
 *
 * @return 0 on success, -1 if no cgroup memory controller with a limit is
 *         found.
 */
static int read_cgroup_memory(unsigned long long *limit, unsigned long long *current,
							  unsigned long long *inactive_file) {
	char dir[300];
	char path[400];
	*inactive_file = 0;

	// cgroup v2
	if (cgroup2_dir(dir, sizeof(dir)) == 0) {
		snprintf(path, sizeof(path), "%s/memory.max", dir);
		if (read_ull_file(path, limit) == 0) {
			snprintf(path, sizeof(path), "%s/memory.current", dir);
			if (read_ull_file(path, current) == 0) {
				snprintf(path, sizeof(path), "%s/memory.stat", dir);
				(void)read_stat_value(path, "inactive_file", inactive_file);
				return 0;
			}
		}
	}

	// cgroup v1
	if (read_ull_file("/sys/fs/cgroup/memory/memory.limit_in_bytes", limit) == 0 &&
		read_ull_file("/sys/fs/cgroup/memory/memory.usage_in_bytes", current) == 0) {
		if (read_stat_value("/sys/fs/cgroup/memory/memory.stat",
							"total_inactive_file", inactive_file) != 0) {
			(void)read_stat_value("/sys/fs/cgroup/memory/memory.stat",
								  "inactive_file", inactive_file);
		}
		return 0;
	}

	return -1;
}

/**
 * @brief Computes the percent of available memory on Linux, taking cgroup limits
 * (Docker, systemd) into account.
 *
 * If the cgroup limit is lower than the host's physical RAM, the limit is used as
 * the total and the available amount is the smaller of the host's and the cgroup's
 * own estimate.
 *
 * @param[out] available_percent  Receives the available memory, as a percent
 *                                (0..100) of the total.
 *
 * @return 0 on success, -1 if the host memory information cannot be read.
 */
static int linux_available_percent(double *available_percent) {
	unsigned long long total = 0, available = 0;
	if (read_linux_host_memory(&total, &available) != 0) {
		return -1;
	}

	unsigned long long limit = 0, current = 0, inactive_file = 0;
	// An "unlimited" cgroup v1 group reports a huge limit, so `limit < total` covers
	// that case as well
	if (read_cgroup_memory(&limit, &current, &inactive_file) == 0 && limit > 0 &&
		limit < total) {
		unsigned long long used =
			current > inactive_file ? current - inactive_file : 0;
		unsigned long long cg_available = limit > used ? limit - used : 0;
		total = limit;
		if (cg_available < available) {
			available = cg_available;
		}
	}

	*available_percent = 100.0 * (double)available / (double)total;
	return 0;
}

#endif /* __linux__ */

int is_system_memory_critical(unsigned int ram_budget_percent) {
	if (ram_budget_percent < 1 || ram_budget_percent > 100) {
		return -1;
	}

	double available_percent = -1.0;

#if defined(_WIN32)
	MEMORYSTATUSEX mem = {0};
	mem.dwLength = sizeof(mem);
	if (!GlobalMemoryStatusEx(&mem) || mem.ullTotalPhys == 0) {
		return -1;
	}
	available_percent = 100.0 * (double)mem.ullAvailPhys / (double)mem.ullTotalPhys;

#elif defined(__linux__)
	if (linux_available_percent(&available_percent) != 0) {
		return -1;
	}

#elif defined(__APPLE__) && defined(__MACH__)
	// Preferred: the kernel's own "percent of memory available"
	int level = 0;
	size_t level_size = sizeof(level);
	if (sysctlbyname("kern.memorystatus_level", &level, &level_size, NULL, 0) == 0 &&
		level >= 0 && level <= 100) {
		available_percent = (double)level;
	} else {
		// Fallback: count free pages, file-backed (droppable) pages and purgeable
		// pages

		uint64_t total_memory = 0;
		size_t total_memory_size = sizeof(total_memory);
		if (sysctlbyname("hw.memsize", &total_memory, &total_memory_size, NULL, 0) !=
				0 ||
			total_memory == 0) {
			return -1;
		}

		// mach_host_self adds a send right on every call; release it.
		mach_port_t host = mach_host_self();
		vm_statistics64_data_t vm_stats = {0};
		mach_msg_type_number_t count = HOST_VM_INFO64_COUNT;
		vm_size_t page_size = 0;
		kern_return_t stats_kr = host_statistics64(host, HOST_VM_INFO64,
												   (host_info64_t)&vm_stats, &count);
		kern_return_t page_kr = host_page_size(host, &page_size);
		mach_port_deallocate(mach_task_self(), host);

		if (stats_kr != KERN_SUCCESS || page_kr != KERN_SUCCESS || page_size == 0) {
			return -1;
		}

		unsigned long long free_pages = vm_stats.free_count;
		unsigned long long speculative = vm_stats.speculative_count;
		if (speculative > free_pages) {
			speculative = free_pages;
		}
		// free_count includes speculative pages, which are file-backed and are
		// already part of external_page_count
		unsigned long long available_pages = (free_pages - speculative) +
											 vm_stats.external_page_count +
											 vm_stats.purgeable_count;
		unsigned long long available_bytes =
			available_pages * (unsigned long long)page_size;
		if (available_bytes > total_memory) {
			available_bytes = total_memory;
		}
		available_percent = 100.0 * (double)available_bytes / (double)total_memory;
	}
#else
	return -1;
#endif

	if (available_percent < 0.0) {
		available_percent = 0.0;
	} else if (available_percent > 100.0) {
		available_percent = 100.0;
	}

	double used_percent = 100.0 - available_percent;
	return used_percent > (double)ram_budget_percent;
}

void release_unused_memory_to_os(void) {
#if defined(__linux__) && defined(__GLIBC__)
	malloc_trim(0);
#endif
}
