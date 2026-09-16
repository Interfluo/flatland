#!/usr/bin/env bash
#
# Builds and runs the C ABI test the way a real consumer would: compiled with
# the C compiler against include/flatland.h, then linked against BOTH the static
# and the shared library. Those two can fail independently — a missing export
# only shows up in the shared case — so both are exercised.
#
# Sourced by run_tests.sh with ROOT and TMP already set.

sect "C ABI"

CC_BIN="${CC:-cc}"
A="$TMP/capi"; mkdir -p "$A"

if make -C "$ROOT" lib >"$A/build.log" 2>&1; then
    ok "libflatland builds"
else
    bad "libflatland builds ($(tail -3 "$A/build.log" | tr '\n' ' '))"
    return 0 2>/dev/null || true
fi

CUBE="$ROOT/examples/cube_area/cube.obj"
# A single-column field with one row per cube vertex, for the file-backed checks.
awk '/^v /{print NR % 7}' "$CUBE" > "$A/field.txt"

SHARED=""
for cand in "$ROOT/libflatland.so" "$ROOT/libflatland.dylib"; do
    [ -f "$cand" ] && SHARED="$cand" && break
done

# Replay the C test's own PASS/FAIL lines into this suite's counters, so a
# failure inside the C program is reported here rather than as one opaque
# non-zero exit.
replay() {
    while IFS= read -r line; do
        case "$line" in
            "  PASS "*) pass=$((pass+1)) ;;
            "  FAIL "*) bad "${line#  FAIL }" ;;
        esac
    done < "$1"
}

# --- static link ------------------------------------------------------------
# -lstdc++ because the engine underneath is C++; that is an implementation
# detail of static linking and does not leak into the shared library's users.
if "$CC_BIN" -std=c99 -Wall -Wextra -I"$ROOT/include" "$ROOT/tests/test_capi.c" \
        "$ROOT/libflatland.a" -lstdc++ -lm -pthread -o "$A/capi_static" 2>"$A/cc_static.log"; then
    ok "the C test compiles as C99 against the static library"
    if [ -s "$A/cc_static.log" ]; then
        bad "the C test compiles without warnings ($(head -1 "$A/cc_static.log"))"
    else
        ok "the C test compiles without warnings"
    fi
    ( cd "$A" && ./capi_static "$CUBE" "$A/field.txt" > "$A/static.out" 2>&1 )
    rc=$?
    replay "$A/static.out"
    [ "$rc" -eq 0 ] && ok "static-linked C test exits 0" \
                    || bad "static-linked C test exits $rc"
else
    bad "the C test compiles as C99 against the static library ($(tail -2 "$A/cc_static.log" | tr '\n' ' '))"
fi

# --- shared link ------------------------------------------------------------
if [ -z "$SHARED" ]; then
    bad "a shared library was built"
elif "$CC_BIN" -std=c99 -Wall -Wextra -I"$ROOT/include" "$ROOT/tests/test_capi.c" \
        -L"$ROOT" -lflatland -lm -pthread -o "$A/capi_shared" 2>"$A/cc_shared.log"; then
    ok "the C test links against the shared library"
    ( cd "$A" && LD_LIBRARY_PATH="$ROOT:${LD_LIBRARY_PATH:-}" \
                 DYLD_LIBRARY_PATH="$ROOT:${DYLD_LIBRARY_PATH:-}" \
        ./capi_shared "$CUBE" "$A/field.txt" > "$A/shared.out" 2>&1 )
    rc=$?
    [ "$rc" -eq 0 ] && ok "shared-linked C test exits 0" \
                    || { bad "shared-linked C test exits $rc"; tail -5 "$A/shared.out" | sed 's/^/      /'; }
else
    bad "the C test links against the shared library ($(tail -2 "$A/cc_shared.log" 2>/dev/null | tr '\n' ' '))"
fi

# --- ABI surface ------------------------------------------------------------
# The shared library must export the C ABI and nothing else. A leaked C++ symbol
# is an ABI-stability hazard: it ties every consumer to this exact compiler and
# standard-library version.
if command -v nm >/dev/null 2>&1 && [ -n "$SHARED" ]; then
    n_c=$(nm -D --defined-only "$SHARED" 2>/dev/null | grep -c ' T fl_')
    n_cxx=$(nm -D --defined-only "$SHARED" 2>/dev/null | grep ' T ' | grep -c '_ZN')
    [ "${n_c:-0}" -ge 25 ] && ok "the shared library exports the C ABI ($n_c symbols)" \
                           || bad "the shared library exports only ${n_c:-0} fl_ symbols"
    [ "${n_cxx:-0}" -eq 0 ] && ok "no C++ symbols leak from the shared library" \
                            || bad "${n_cxx} C++ symbols leak from the shared library"
fi

# --- header hygiene ---------------------------------------------------------
# The header is the contract; it must stay clean C under strict ISO rules and
# must also be includable from C++ (the extern "C" guard).
printf '#include "flatland.h"\nint main(void){fl_options o; fl_options_init(&o); return o.cull==1?0:1;}\n' > "$A/hdr.c"
if "$CC_BIN" -std=c99 -pedantic -Wall -Wextra -Werror -I"$ROOT/include" -fsyntax-only "$A/hdr.c" 2>"$A/hdr.log"; then
    ok "flatland.h is clean under -std=c99 -pedantic -Werror"
else
    bad "flatland.h under -std=c99 -pedantic -Werror ($(head -2 "$A/hdr.log" | tr '\n' ' '))"
fi
printf '#include "flatland.h"\nint main(){fl_options o; fl_options_init(&o); return o.cull==1?0:1;}\n' > "$A/hdr.cpp"
if "${CXX:-c++}" -std=c++17 -Wall -Wextra -Werror -I"$ROOT/include" -fsyntax-only "$A/hdr.cpp" 2>"$A/hdrxx.log"; then
    ok "flatland.h is clean when included from C++"
else
    bad "flatland.h from C++ ($(head -2 "$A/hdrxx.log" | tr '\n' ' '))"
fi

# --- locale independence ----------------------------------------------------
# Every format FlatLand reads writes numbers with a '.', but strtod uses the
# decimal separator of the process's LC_NUMERIC locale. The CLI has always
# guarded this with setlocale(LC_ALL,"C") in main(); a library cannot do that
# without corrupting state it does not own, so it scopes the C locale to the
# parsing thread instead. This matters concretely for the bindings: MATLAB and
# many desktop Python applications call setlocale(LC_ALL,"") at startup, and
# comma-decimal locales are the default across much of Europe.
COMMA_LOCALE=""
for cand in de_DE.utf8 de_DE.UTF-8 fr_FR.utf8 fr_FR.UTF-8 es_ES.utf8 ru_RU.utf8 nl_NL.utf8; do
    if locale -a 2>/dev/null | grep -qxF "$cand"; then COMMA_LOCALE="$cand"; break; fi
done

if [ -z "$COMMA_LOCALE" ]; then
    ok "locale independence (skipped: no comma-decimal locale installed)"
else
    cat > "$A/loc.c" <<'EOF'
#include "flatland.h"
#include <locale.h>
#include <stdio.h>
#include <string.h>
int main(int argc, char** argv) {
    if (argc < 3) return 2;
    if (!setlocale(LC_ALL, argv[1])) { printf("SKIP\n"); return 0; }
    char before[128];
    snprintf(before, sizeof before, "%s", setlocale(LC_NUMERIC, NULL));

    fl_mesh* m = NULL;
    if (fl_mesh_load(argv[2], &m) != FL_OK) { printf("LOAD_FAILED %s\n", fl_last_error()); return 1; }
    fl_options o; fl_options_init(&o);
    o.resolution = 0.01; o.precision = FL_PRECISION_DOUBLE;
    double v[3] = {1,0,0}; fl_result r;
    if (fl_project(m, v, NULL, 0, &o, &r) != FL_OK) { printf("PROJECT_FAILED\n"); return 1; }
    fl_mesh_destroy(m);

    /* setlocale returns a pointer into its own static buffer, which the next
       setlocale call overwrites — so copy it out before switching to C. */
    char after[128];
    snprintf(after, sizeof after, "%s", setlocale(LC_NUMERIC, NULL));

    /* Print in the C locale so this program's own output stays parseable. */
    setlocale(LC_ALL, "C");
    printf("AREA %.6f VERTS_OK %d LOCALE_PRESERVED %d\n",
           r.area, 1, strcmp(before, after) == 0 ? 1 : 0);
    return 0;
}
EOF
    if "$CC_BIN" -std=c99 -I"$ROOT/include" "$A/loc.c" "$ROOT/libflatland.a" \
            -lstdc++ -lm -pthread -o "$A/loc" 2>"$A/loc_cc.log"; then
        OUT=$("$A/loc" "$COMMA_LOCALE" "$CUBE" 2>&1)
        case "$OUT" in
            SKIP*)
                ok "locale independence (skipped: $COMMA_LOCALE not usable)" ;;
            LOAD_FAILED*)
                bad "a mesh loads under $COMMA_LOCALE (${OUT#LOAD_FAILED })" ;;
            AREA*)
                AREA=$(echo "$OUT" | awk '{print $2}')
                PRES=$(echo "$OUT" | awk '{print $6}')
                near "$AREA" 1.0 0.01 "a mesh parses correctly under $COMMA_LOCALE"
                equal "$PRES" "1" "...and the host process's locale is left unchanged" ;;
            *)
                bad "locale test produced unexpected output ($OUT)" ;;
        esac
    else
        bad "the locale test compiles ($(tail -1 "$A/loc_cc.log"))"
    fi
fi
