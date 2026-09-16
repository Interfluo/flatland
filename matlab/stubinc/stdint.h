/*
 * matlab/stubinc/stdint.h — minimal <stdint.h> shadow for MATLAB's loadlibrary.
 * See stddef.h in this directory for why the shadow exists at all.
 *
 * FlatLand's public header uses exactly four names from <stdint.h>:
 *     int32_t  int64_t  uint8_t  uint32_t
 * The rest are provided so the file is a usable <stdint.h> if anything else
 * in the parse reaches for them. Widths and signedness match the real header
 * on every platform FlatLand builds for; only the parser ever sees this file.
 */
#ifndef FLATLAND_MATLAB_STDINT_H
#define FLATLAND_MATLAB_STDINT_H

typedef signed char        int8_t;
typedef unsigned char      uint8_t;
typedef short              int16_t;
typedef unsigned short     uint16_t;
typedef int                int32_t;
typedef unsigned int       uint32_t;

/* int64_t is "long" on LP64 (Linux/macOS) and "long long" on Windows.
 * Both are 8 bytes, so the ABI is identical either way; this just keeps the
 * spelling honest. */
#if defined(_WIN32) || defined(_WIN64)
typedef long long          int64_t;
typedef unsigned long long uint64_t;
#elif defined(__LP64__) || defined(_LP64) || defined(__x86_64__) || defined(__aarch64__)
typedef long               int64_t;
typedef unsigned long      uint64_t;
#else
typedef long long          int64_t;
typedef unsigned long long uint64_t;
#endif

#endif /* FLATLAND_MATLAB_STDINT_H */
