#ifndef MEMORY_UTILS_HPP
#define MEMORY_UTILS_HPP

// Platform-specific includes for memory tracking
#ifdef _WIN32
  #include "gflow/win_compat.hpp"  // << replaces <windows.h> and <psapi.h>
#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
  #include <unistd.h>
  #include <sys/resource.h>
  #if defined(__APPLE__) && defined(__MACH__)
    #include <mach/mach.h>
  #endif
#endif

#include <cstddef>   // size_t
#include <cinttypes> // PRIu64, PRId64 (if you switch to uint64_t/int64_t)
#include <R.h>       // Rprintf, R_FlushConsole

using std::size_t;

size_t get_current_rss();  // defined in .cpp

struct memory_tracker_t {
    size_t initial_usage;
    const char* context;

    explicit memory_tracker_t(const char* ctx) : initial_usage(get_current_rss()), context(ctx) {}

    void report() {
        size_t current = get_current_rss();
        // compute a signed diff to avoid size_t underflow on transient drops
        ptrdiff_t diff = (ptrdiff_t)current - (ptrdiff_t)initial_usage;
        Rprintf("%s Memory: Current=%zu MB, Diff=%+td MB\n",
                context, current / 1024 / 1024, diff / (ptrdiff_t)1024 / 1024);
        R_FlushConsole();
    }
};

#endif // MEMORY_UTILS_HPP
