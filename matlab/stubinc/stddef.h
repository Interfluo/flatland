/*
 * matlab/stubinc/stddef.h — minimal <stddef.h> shadow for MATLAB's loadlibrary.
 *
 * WHY THIS EXISTS
 * ---------------
 * flatland.h includes <stddef.h> and <stdint.h>. On Linux, GCC's real
 * <stddef.h> defines max_align_t as
 *
 *     typedef struct {
 *       long long   __max_align_ll __attribute__((__aligned__(...)));
 *       long double __max_align_ld __attribute__((__aligned__(...)));
 *     } max_align_t;
 *
 * which contains THREE constructs MATLAB's restricted C parser is documented
 * not to handle: long double, __attribute__, and an unnamed struct in a
 * typedef. glibc's <stdint.h> additionally drags in a large pile of
 * feature-test machinery for no benefit here.
 *
 * flatland_load.m passes this directory to loadlibrary via 'includepath'.
 * Preprocessors search -I directories before the standard system directories
 * for <...> includes, so this file shadows the real one for the duration of
 * the header parse only. Nothing compiled or linked uses it -- the .so on disk
 * was built against the real system headers, and these typedefs are chosen to
 * be ABI-identical to them.
 *
 * FlatLand needs exactly one name from <stddef.h>: size_t.
 */
#ifndef FLATLAND_MATLAB_STDDEF_H
#define FLATLAND_MATLAB_STDDEF_H

#ifndef NULL
#define NULL 0
#endif

/* size_t must have the same width and signedness as the platform's real one,
 * or every count argument would marshal wrong. */
#if defined(_WIN64)
typedef unsigned long long size_t;
typedef long long          ptrdiff_t;
#elif defined(_WIN32)
typedef unsigned int size_t;
typedef int          ptrdiff_t;
#elif defined(__LP64__) || defined(_LP64) || defined(__x86_64__) || defined(__aarch64__)
typedef unsigned long size_t;    /* LP64: 8 bytes, matches glibc/musl/Darwin */
typedef long          ptrdiff_t;
#else
typedef unsigned int size_t;
typedef int          ptrdiff_t;
#endif

#endif /* FLATLAND_MATLAB_STDDEF_H */
