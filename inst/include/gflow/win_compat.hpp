#pragma once
#if defined(_WIN32)

// keep headers small + avoid min/max macros
#ifndef WIN32_LEAN_AND_MEAN
#  define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#  define NOMINMAX
#endif

// Temporarily neutralize R's macros that collide with COM
#ifdef Realloc
#  define GFLOW__R_REALLOC_SAVED Realloc
#  undef Realloc
#endif
#ifdef Free
#  define GFLOW__R_FREE_SAVED Free
#  undef Free
#endif

#  include <windows.h>   // pulls base Win32 stuff
#  include <psapi.h>     // GetProcessMemoryInfo

// Restore R macros
#ifdef GFLOW__R_REALLOC_SAVED
#  define Realloc GFLOW__R_REALLOC_SAVED
#  undef GFLOW__R_REALLOC_SAVED
#endif
#ifdef GFLOW__R_FREE_SAVED
#  define Free GFLOW__R_FREE_SAVED
#  undef GFLOW__R_FREE_SAVED
#endif

#endif // _WIN32
