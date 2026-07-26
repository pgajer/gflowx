#ifndef EIGEN_WARNINGS_DISABLED
#define EIGEN_WARNINGS_DISABLED

/* gflow local patch:
 * Keep Eigen include semantics while avoiding compiler-specific diagnostic
 * suppression pragmas that trigger CRAN NOTE checks.
 */

#else
// warnings already disabled:
# ifndef EIGEN_WARNINGS_DISABLED_2
#  define EIGEN_WARNINGS_DISABLED_2
# elif defined(EIGEN_INTERNAL_DEBUGGING)
#  error "Do not include \"DisableStupidWarnings.h\" recursively more than twice!"
# endif

#endif // not EIGEN_WARNINGS_DISABLED
