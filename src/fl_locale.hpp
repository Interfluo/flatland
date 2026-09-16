// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Interfluo
//
// FlatLand is dual-licensed: GNU AGPL v3 (see LICENSE) or a commercial licence
// for closed-source or hosted use (see COMMERCIAL-LICENSE.md).

#pragma once

/*
 * Locale-independent numeric parsing.
 *
 * strtod() is specified to use the decimal separator of the current LC_NUMERIC
 * locale. Every file format FlatLand reads — OBJ, ASCII STL, field matrices,
 * batch files — writes numbers with a '.', so parsing them inside a process
 * whose locale uses ',' reads "0.5" as 0 and stops at the separator.
 *
 * The CLI has always guarded against this by calling setlocale(LC_ALL, "C") in
 * main(). The LIBRARY cannot: a shared library that changes its host process's
 * global locale would be corrupting state it does not own, and would do it to
 * every other thread at the same time. That matters here specifically because
 * the bindings put FlatLand inside exactly the kind of host that sets a locale —
 * MATLAB and many desktop Python applications call setlocale(LC_ALL, "") at
 * startup, and comma-decimal locales are the default across most of Europe.
 *
 * So parsing is instead wrapped in a scope that switches only THIS THREAD to
 * the C locale and restores it on the way out, leaving the rest of the process
 * untouched. POSIX 2008 does this with uselocale(); Windows with
 * _configthreadlocale(). Where neither is available the scope compiles to
 * nothing and behaviour is as it was before — no worse, and not silently
 * pretending to be fixed.
 */

#if defined(_WIN32)
#  include <clocale>
#  include <cstdlib>
#  include <cstring>
#  define FL_HAVE_THREAD_LOCALE 1
#elif defined(__unix__) || defined(__APPLE__) || defined(__linux__)
#  include <locale.h>
#  define FL_HAVE_THREAD_LOCALE 1
#else
#  define FL_HAVE_THREAD_LOCALE 0
#endif

namespace flatland {

// RAII: for as long as one of these is alive, numeric conversions on the
// calling thread use the C locale. Nothing outside this thread is affected.
class CNumericScope {
public:
    CNumericScope() {
#if defined(_WIN32)
        old_mode_ = _configthreadlocale(_ENABLE_PER_THREAD_LOCALE);
        const char* cur = std::setlocale(LC_NUMERIC, nullptr);
        if (cur && std::strcmp(cur, "C") != 0) {
            saved_ = _strdup(cur);
            std::setlocale(LC_NUMERIC, "C");
        }
#elif FL_HAVE_THREAD_LOCALE
        c_locale_ = newlocale(LC_NUMERIC_MASK, "C", (locale_t)0);
        if (c_locale_ != (locale_t)0) previous_ = uselocale(c_locale_);
#endif
    }

    ~CNumericScope() {
#if defined(_WIN32)
        if (saved_) {
            std::setlocale(LC_NUMERIC, saved_);
            std::free(saved_);
        }
        if (old_mode_ != -1) _configthreadlocale(old_mode_);
#elif FL_HAVE_THREAD_LOCALE
        if (c_locale_ != (locale_t)0) {
            if (previous_ != (locale_t)0) uselocale(previous_);
            freelocale(c_locale_);
        }
#endif
    }

    CNumericScope(const CNumericScope&) = delete;
    CNumericScope& operator=(const CNumericScope&) = delete;

private:
#if defined(_WIN32)
    int   old_mode_ = -1;
    char* saved_    = nullptr;
#elif FL_HAVE_THREAD_LOCALE
    locale_t c_locale_ = (locale_t)0;
    locale_t previous_ = (locale_t)0;
#endif
};

} // namespace flatland
