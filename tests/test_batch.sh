#!/usr/bin/env bash
# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Interfluo
#
# Batch files and the time-series workflow.
#
# Sourced by run_tests.sh with BIN, ROOT and TMP already set.

CUBE="$ROOT/examples/cube_area/cube.obj"
B="$TMP/batch"; mkdir -p "$B/case"
NV=$(grep -c '^v ' "$CUBE")

sect "Batch basics"

seq 1 "$NV" | awk '{print 2.0}' > "$B/case/f.txt"
printf '1 0 0 0.01 f.txt\n0 1 0 0.01 f.txt\n' > "$B/case/ts.txt"
"$BIN" "$CUBE" -b "$B/case/ts.txt" -j 2>/dev/null > "$B/b.json"
equal "$(grep -c '"idx"' "$B/b.json")" "2" "batch produced 2 results"
near "$(jget "$B/b.json" average)" 2.0 0.001 "relative data path resolves against the batch dir"

# A data token with no resolution must still attach the field.
printf '1 0 0 f.txt\n' > "$B/case/nores.txt"
"$BIN" "$CUBE" -b "$B/case/nores.txt" -j 2>/dev/null > "$B/nr.json"
near "$(jget "$B/nr.json" average)" 2.0 0.001 "data file without a resolution still applies"

# Comments, blank lines and CRLF.
printf '# leading comment\n\n1 0 0 0.01 f.txt\r\n\n# trailing\n' > "$B/case/messy.txt"
"$BIN" "$CUBE" -b "$B/case/messy.txt" -j 2>/dev/null > "$B/messy.json"
equal "$(grep -c '"idx"' "$B/messy.json")" "1" "comments, blanks and CRLF are tolerated"

# Angle lines.
printf 'a 90 0 0.01\n' > "$B/case/ang.txt"
"$BIN" "$CUBE" -b "$B/case/ang.txt" -j 2>/dev/null > "$B/ab.json"
near "$(jget "$B/ab.json" area)" 1.0 0.01 "batch angle line area == 1"

# -v/-a and -b are additive.
"$BIN" "$CUBE" -v 1 0 0 -b "$B/case/ang.txt" -r 0.01 -j 2>/dev/null > "$B/add.json"
equal "$(grep -c '"idx"' "$B/add.json")" "2" "-v and -b views are additive"

sect "Batch token handling"

# The help text and README both promise that [resolution] and [data] are
# order-independent. They were not: a data token in first position left the
# resolution silently unconsumed and the view fell back to the default.
printf '1 0 0 0.01 f.txt\n' > "$B/case/rd.txt"
printf '1 0 0 f.txt 0.01\n' > "$B/case/dr.txt"
"$BIN" "$CUBE" -b "$B/case/rd.txt" -j 2>/dev/null > "$B/rd.json"
"$BIN" "$CUBE" -b "$B/case/dr.txt" -j 2>/dev/null > "$B/dr.json"
equal "$(jget "$B/dr.json" pixels)" "$(jget "$B/rd.json" pixels)" "batch tokens are order-independent (resolution applied)"
near "$(jget "$B/dr.json" average)" "$(jget "$B/rd.json" average)" 0.0001 "batch tokens are order-independent (field)"

# Unrecognised extra tokens are a typo, not something to ignore.
printf '1 0 0 0.01 f.txt surplus\n' > "$B/case/extra.txt"
"$BIN" "$CUBE" -b "$B/case/extra.txt" >/dev/null 2>&1; expect_fail $? "surplus batch token rejected"

# Time-series data files are routinely named for their timestep. Sniffing "a
# pure number means resolution" silently ate the filename and dropped the field.
cp "$B/case/f.txt" "$B/case/0100"
printf '1 0 0 0100\n' > "$B/case/numname.txt"
"$BIN" "$CUBE" -b "$B/case/numname.txt" -r 0.02 -j 2>/dev/null > "$B/nn.json"
near "$(jget "$B/nn.json" average)" 2.0 0.001 "a numerically-named data file is treated as data"

# A '#' inside a path is a path character, not the start of a comment.
mkdir -p "$B/case/h#d" && cp "$B/case/f.txt" "$B/case/h#d/f.txt"
printf '1 0 0 0.01 h#d/f.txt\n' > "$B/case/hash.txt"
"$BIN" "$CUBE" -b "$B/case/hash.txt" -j 2>/dev/null > "$B/hash.json"
near "$(jget "$B/hash.json" average)" 2.0 0.001 "a '#' inside a path does not start a comment"

# Malformed batch lines.
printf '1 0\n' > "$B/case/short.txt"
"$BIN" "$CUBE" -b "$B/case/short.txt" >/dev/null 2>&1; expect_fail $? "short batch line rejected"
printf '0 0 0\n' > "$B/case/zero.txt"
"$BIN" "$CUBE" -b "$B/case/zero.txt" >/dev/null 2>&1; expect_fail $? "zero view vector in a batch rejected"
printf '1 0 0 nan\n' > "$B/case/nanres.txt"
"$BIN" "$CUBE" -b "$B/case/nanres.txt" >/dev/null 2>&1; expect_fail $? "non-finite batch resolution rejected"
"$BIN" "$CUBE" -b "$B/case/missing.txt" >/dev/null 2>&1; expect_fail $? "missing batch file rejected"

sect "Time series over one field matrix"

# The advertised workflow: one matrix, one column per timestep, one batch line
# per timestep. Column c holds the constant c+1, so view c must average c+1.
awk -v n="$NV" 'BEGIN{for(i=0;i<n;i++){for(c=0;c<6;c++) printf "%d ", c+1; print ""}}' > "$B/case/matrix.txt"
: > "$B/case/series.txt"
for c in 0 1 2 3 4 5; do printf '1 0 0 0.02 matrix.txt@%d\n' "$c" >> "$B/case/series.txt"; done
"$BIN" "$CUBE" -b "$B/case/series.txt" -j 2>/dev/null > "$B/series.json"
equal "$(grep -c '"idx"' "$B/series.json")" "6" "6-step time series produced 6 results"
for c in 1 2 3 4 5 6; do
    near "$(jget_n "$B/series.json" average "$c")" "$c" 0.001 "timestep $((c-1)) reads column $((c-1))"
done

# Results must be identical however many workers ran them, and must stay in
# view order regardless of completion order.
"$BIN" "$CUBE" -b "$B/case/series.txt" -t 1 -j 2>/dev/null | grep -v '"time"' > "$B/t1.json"
"$BIN" "$CUBE" -b "$B/case/series.txt" -t 8 -j 2>/dev/null | grep -v '"time"' > "$B/t8.json"
if diff -q "$B/t1.json" "$B/t8.json" >/dev/null; then
    ok "batch output is independent of thread count"
else
    bad "batch output differs between -t 1 and -t 8"
fi

# Many distinct field files in one batch: exercises the matrix cache under
# eviction pressure. Each file is constant k, so view k must average k.
mkdir -p "$B/many"
for k in $(seq 0 19); do
    seq 1 "$NV" | awk -v k="$k" '{print k}' > "$B/many/f$k.txt"
    printf '1 0 0 0.05 f%d.txt\n' "$k" >> "$B/many/batch.txt"
done
"$BIN" "$CUBE" -b "$B/many/batch.txt" -t 4 -j 2>/dev/null > "$B/many.json"
equal "$(grep -c '"idx"' "$B/many.json")" "20" "20 distinct field files all processed"
near "$(jget_n "$B/many.json" average 1)"  0 0.001 "first of 20 distinct fields is correct"
near "$(jget_n "$B/many.json" average 20)" 19 0.001 "last of 20 distinct fields is correct"

sect "Batch failure handling"

# The reported failure must be the lowest-indexed one, or the same input yields
# a different message run to run and CI output becomes unreproducible.
: > "$B/case/errs.txt"
for k in $(seq 0 7); do printf '1 0 0 0.02 no_such_file_%d.txt\n' "$k" >> "$B/case/errs.txt"; done
FIRST=""
SAME=1
for _ in $(seq 1 10); do
    MSG=$("$BIN" "$CUBE" -b "$B/case/errs.txt" -t 8 2>&1 >/dev/null | grep -oE 'timestep [0-9]+' | head -1)
    [ -z "$FIRST" ] && FIRST="$MSG"
    [ "$MSG" = "$FIRST" ] || SAME=0
done
equal "$FIRST" "timestep 0" "a failing batch reports the lowest-indexed failure"
equal "$SAME" "1" "...and reports it deterministically across 10 runs"

# An aborted batch must not leave half-written images behind.
rm -rf "$B/abortout"; mkdir -p "$B/abortout"
: > "$B/case/partial.txt"
for k in $(seq 0 5); do printf '1 0 0 0.05 f.txt\n' >> "$B/case/partial.txt"; done
printf '1 0 0 0.05 no_such_file.txt\n' >> "$B/case/partial.txt"
"$BIN" "$CUBE" -b "$B/case/partial.txt" -t 2 -o "$B/abortout/p" >/dev/null 2>&1
expect_fail $? "batch with a bad data file fails"
LEFT=$(find "$B/abortout" -name '*.png' | wc -l | tr -d ' ')
equal "$LEFT" "0" "...leaving no partial images behind"

sect "Committed example case"

"$BIN" "$ROOT/examples/sphere_areas/sphere_coarse.obj" \
       -b "$ROOT/examples/with_fields/timeseries_demo/timeseries.txt" -j 2>/dev/null > "$B/demo.json"
expect_ok $? "committed time-series demo runs"
equal "$(grep -c '"idx"' "$B/demo.json")" "24" "committed demo produced 24 timesteps"
