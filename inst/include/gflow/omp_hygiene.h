// This to be drop in right before pragma-heavy code - ideally the last include in any file using OpenMP
#pragma once
#ifdef match
#  undef match
#endif
#ifdef check
#  undef check
#endif
