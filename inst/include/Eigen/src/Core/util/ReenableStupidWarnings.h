#ifdef EIGEN_WARNINGS_DISABLED_2
// "DisableStupidWarnings.h" was included twice recursively: Do not reenable warnings yet!
#  undef EIGEN_WARNINGS_DISABLED_2

#elif defined(EIGEN_WARNINGS_DISABLED)
#undef EIGEN_WARNINGS_DISABLED

/* gflow local patch:
 * DisableStupidWarnings.h no longer pushes diagnostic state.
 * Keep this include as a no-op for compatibility.
 */

#endif // EIGEN_WARNINGS_DISABLED
