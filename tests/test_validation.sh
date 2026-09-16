#!/usr/bin/env bash
#
# Validation against closed-form results.
#
# Distinct in kind from the other suites: those pin BEHAVIOUR, this one checks
# CORRECTNESS. Every expectation is derived analytically in validation/cases.py —
# projected areas of convex polyhedra, Cauchy's mean-projection identity, area
# integrals of linear and Lambertian fields — so nothing here can be satisfied by
# blessing whatever FlatLand happened to print.
#
# Runs the quick subset. The full study, including the convergence sweeps behind
# docs/VALIDATION.md, is `make validate`.
#
# Sourced by run_tests.sh with BIN, ROOT and TMP already set.

sect "Validation against closed-form results"

PY_BIN="${PYTHON:-python3}"
V="$TMP/validation"; mkdir -p "$V"

if ! command -v "$PY_BIN" >/dev/null 2>&1; then
    ok "validation studies (skipped: $PY_BIN not available)"
    return 0 2>/dev/null || true
fi

# The analytic formulas are themselves checked first. If cases.py were wrong,
# every study built on it would agree with the wrong answer and report success.
if "$PY_BIN" "$ROOT/validation/self_check.py" > "$V/self.out" 2>&1; then
    n=$(grep -c '^  ok ' "$V/self.out" || echo 0)
    ok "the analytic formulas self-check ($n identities)"
else
    bad "the analytic formulas self-check"
    grep -E '^  BAD' "$V/self.out" | head -5 | sed 's/^/      /'
fi

if "$PY_BIN" "$ROOT/validation/run_validation.py" --flatland "$BIN" --quick \
        --json "$V/results.json" > "$V/run.out" 2>&1; then
    RC=0
else
    RC=$?
fi

# Report each study group as one assertion, naming the worst case inside it, so a
# failure points at the physics rather than at a line number.
if [ -f "$V/results.json" ]; then
    "$PY_BIN" - "$V/results.json" <<'PY' > "$V/summary.txt"
import json, sys
rows = json.load(open(sys.argv[1]))
groups = {}
for r in rows:
    groups.setdefault(r["study"], []).append(r)
for study, rs in groups.items():
    bad = [r for r in rs if not r["ok"]]
    worst = max(rs, key=lambda r: r["rel_error"])
    if bad:
        w = max(bad, key=lambda r: r["rel_error"])
        print("FAIL\t%s\t%d/%d outside tolerance, worst: %s (%.2e > %.2e)"
              % (study, len(bad), len(rs), w["case"], w["rel_error"], w["tol"]))
    else:
        print("PASS\t%s\t%d checks, worst relative error %.2e (%s)"
              % (study, len(rs), worst["rel_error"], worst["case"]))
PY
    while IFS=$'\t' read -r verdict study detail; do
        [ -z "$verdict" ] && continue
        if [ "$verdict" = "PASS" ]; then
            ok "$study — $detail"
        else
            bad "$study — $detail"
        fi
    done < "$V/summary.txt"
else
    bad "validation produced no results (exit $RC)"
    tail -12 "$V/run.out" | sed 's/^/      /'
fi

[ "$RC" -eq 0 ] && ok "all validation studies within tolerance" \
                || bad "validation exited $RC"
