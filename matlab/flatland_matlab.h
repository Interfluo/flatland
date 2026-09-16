/*
 * flatland_matlab.h — the header MATLAB's loadlibrary parses.
 *
 * This file adds NOTHING to the ABI. It is a two-line shim whose only job is to
 * present include/flatland.h to loadlibrary's restricted C parser in a form that
 * parser accepts. The declarations, types and calling conventions all come from
 * the real header, so this shim cannot drift away from the contract.
 *
 * THE PROBLEM IT SOLVES
 * ---------------------
 * Every function in flatland.h is declared FL_API. On Linux/macOS with GCC or
 * Clang that macro expands to
 *
 *     __attribute__((visibility("default")))
 *
 * and on Windows to __declspec(dllimport). loadlibrary preprocesses the header
 * with a real C compiler and then parses the RESULT with its own restricted C
 * parser, so what that parser sees is a declaration beginning with a GCC
 * attribute. Verified with `gcc -E` against this repo:
 *
 *     __attribute__((visibility("default"))) fl_status fl_mesh_create(...)
 *
 * Attributes are not part of the C subset loadlibrary documents, so the shim
 * makes FL_API expand to nothing before the header is read.
 *
 * Doing that takes two steps, because flatland.h's linkage block is
 *
 *     #if   defined(_WIN32) || defined(__CYGWIN__)
 *     #  if defined(FLATLAND_STATIC)   -> FL_API empty
 *     ...
 *     #elif defined(__GNUC__) && __GNUC__ >= 4
 *                                      -> FL_API __attribute__((visibility(...)))
 *     #else                            -> FL_API empty
 *     #endif
 *
 * FLATLAND_STATIC alone only helps on Windows; on Linux the __GNUC__ arm wins
 * regardless. So:
 *
 *   1. define FLATLAND_STATIC   — neutralises FL_API if MATLAB preprocesses as
 *                                 Windows (also correct on Cygwin/MinGW), and
 *   2. #undef __GNUC__          — neutralises it if MATLAB preprocesses with
 *                                 GCC or Clang.
 *
 * With neither macro set the header's #else arm already leaves FL_API empty, so
 * all three preprocessor personalities land on an empty FL_API. This is a
 * header-parsing concern only: it does not affect the symbols the .so exports,
 * which are chosen at build time by the real FL_API and are all present.
 *
 * WHY THE SYSTEM HEADERS COME FIRST
 * ---------------------------------
 * flatland.h includes <stddef.h> and <stdint.h>. Undefining __GNUC__ before
 * those are read would break them, since GCC's own headers branch on it. They
 * are therefore pulled in here FIRST, while __GNUC__ is still intact; the
 * copies in this directory's stubinc/ shadow the system ones (flatland_load.m
 * passes stubinc via loadlibrary's 'includepath') and exist because GCC's real
 * <stddef.h> defines max_align_t using long double, __attribute__ and an
 * unnamed struct, all of which are outside loadlibrary's subset. See
 * stubinc/stddef.h for the details. flatland.h's own includes then no-op on
 * their include guards.
 */
#ifndef FLATLAND_MATLAB_H
#define FLATLAND_MATLAB_H

#include <stddef.h>
#include <stdint.h>

#undef  __GNUC__
#define FLATLAND_STATIC 1

/* Resolved relative to THIS file, so no include path is required to find the
 * real header in a source checkout. To parse an installed copy instead, define
 * FLATLAND_MATLAB_USE_INSTALLED and put its directory on the include path. */
#if defined(FLATLAND_MATLAB_USE_INSTALLED)
#  include <flatland.h>
#else
#  include "../include/flatland.h"
#endif

#endif /* FLATLAND_MATLAB_H */
