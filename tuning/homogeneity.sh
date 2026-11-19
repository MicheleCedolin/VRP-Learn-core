#!/usr/bin/env bash
# Does the ranking of parameter configurations depend on the instance size?  (homogeneity check in the sense of Schneider & Hoos 2012)
# N configurations (the default one + N-1 drawn uniformly from the ranges of the parameter table; log scale where the table says so)
# are run once on every instance of a list with two size groups; per instance the configurations are ranked by cost, and the rank
# vectors of the instances are correlated (Spearman).  If configurations rank alike on small and on large instances, one tuning for both
# is plausible; if the between-group correlation is far below the within-group correlations, size classes are worth tuning separately.
# usage: tuning/homogeneity.sh METHOD [-N 20] [-F 6] [-P 4] [-S 1] [-R 1] [-I tuning/split/homog.txt] [-O results/homog]
#   list file: path TAB group (two groups; the first one in the file is reported first).   -R: seed of the random configurations
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."; ROOT=$PWD
METHOD=${1:?usage: tuning/homogeneity.sh METHOD [-N n] [-F s] [-P n] [-S seed] [-R seed] [-I list] [-O dir]}; shift
N=20; F=6; P=4; SEED=1; RSEED=1; LIST=tuning/split/homog.txt; OUT=results/homog
while [ $# -gt 0 ]; do
  case "$1" in -N) N=$2; shift 2;; -F) F=$2; shift 2;; -P) P=$2; shift 2;; -S) SEED=$2; shift 2;; -R) RSEED=$2; shift 2;; -I) LIST=$2; shift 2;; -O) OUT=$2; shift 2;;
    *) echo "unknown option $1" >&2; exit 2;; esac
done
OUT=$(realpath -m "$OUT"); LIST=$(realpath "$LIST")
[ -x cpp/build/vrp_experiments ] && [ -x cpp/build/vrp ] || { echo "build first (see README)" >&2; exit 1; }
[ "${SKIP_DATA_GUARD:-0}" = 1 ] || tuning/check_data.sh
mkdir -p "$OUT/params"; : > "$OUT/trlog.txt"
cpp/build/vrp_experiments params --methods "$METHOD" --out "$OUT/params" > "$OUT/params/table.txt"
PT="$OUT/params/$METHOD.txt"; PD="$OUT/params/${METHOD}_default.txt"
[ -s "$PT" ] || { echo "$METHOD has no tunable parameters" >&2; exit 2; }

# configurations: id 1 = default, ids 2..N random;  configs.tsv: id TAB name=value ...
awk -v n="$N" -v seed="$RSEED" '
  BEGIN { srand(seed) }
  FILENAME == ARGV[1] { if (NF == 0) next; k++; name[k] = $1; line = $0; sub(/^[^"]*"[^"]*"[ \t]*/, "", line);
                         split(line, a, /[ \t(),]+/); type[k] = a[1]; lo[k] = a[2] + 0; hi[k] = a[3] + 0; next }
  FNR == 1 { next }
  { printf "1"; for (i = 1; i <= k; i++) printf "\t%s=%s", name[i], $i; printf "\n";
    for (c = 2; c <= n; c++) {
      printf "%d", c
      for (i = 1; i <= k; i++) {
        u = rand(); islog = (type[i] ~ /log/); isint = (type[i] ~ /^i/)
        if (islog) x = exp(log(lo[i]) + u * (log(hi[i]) - log(lo[i]))); else x = lo[i] + u * (hi[i] - lo[i])
        if (isint) { x = int(x + 0.5); if (x < lo[i]) x = lo[i]; if (x > hi[i]) x = hi[i]; printf "\t%s=%d", name[i], x }
        else printf "\t%s=%.4g", name[i], x
      }
      printf "\n"
    } }' "$PT" "$PD" > "$OUT/configs.tsv"
echo "$(wc -l < "$OUT/configs.tsv") configurations (1 = default) in $OUT/configs.tsv"

# jobs: largest instances first, so that the cores are busy until the end
while IFS=$'\t' read -r path group; do
  [ -s "$path" ] || { echo "instance missing: $path" >&2; exit 1; }
  printf '%s\t%s\t%s\n' "$(cpp/build/vrp info "$path" | head -1 | cut -f3)" "$path" "$group"
done < "$LIST" | sort -rn > "$OUT/instances_by_size.tsv"
id=0; : > "$OUT/jobs.txt"
while IFS=$'\t' read -r nn path group; do
  id=$((id + 1))
  while IFS=$'\t' read -r cid rest; do
    args=$(echo "$rest" | awk -F'\t' '{for (i = 1; i <= NF; i++) { split($i, kv, "="); printf " --%s %s", kv[1], kv[2] }}')
    echo "$cid $id $SEED $path$args" >> "$OUT/jobs.txt"
  done < "$OUT/configs.tsv"
done < "$OUT/instances_by_size.tsv"
export METHOD TIME_PER_100=$F VRP_ROOT=$ROOT TRLOG=$OUT/trlog.txt ALARMS=$OUT/alarms.txt FAILDIR=$OUT/failed
echo "$(date -Is) $METHOD: $(wc -l < "$OUT/jobs.txt") runs ($N configurations x $(wc -l < "$OUT/instances_by_size.tsv") instances), F=$F, parallel $P, host $(hostname)"
xargs -a "$OUT/jobs.txt" -L1 -P "$P" tuning/target-runner > /dev/null || echo "some runs failed (see trlog)"
echo "$(date -Is) done; runs by status:"; awk '{c[$7]++} END{for (s in c) print "  " s, c[s]}' "$OUT/trlog.txt"
echo "=== raw results (conf instance-id instance status cost) ==="; awk '{print $3, $4, $6, $7, $8}' "$OUT/trlog.txt" | sort -k2,2n -k1,1n

# analysis: ranks per instance, Spearman correlations between instances
awk -v nconf="$N" '
  function rankvec(inst,   i, j, c, lt, eq) {                 # rank 1 = lowest cost, ties get the mean rank
    for (i = 1; i <= nc; i++) { lt = 0; eq = 0
      for (j = 1; j <= nc; j++) { if (cost[cf[j], inst] < cost[cf[i], inst]) lt++; else if (cost[cf[j], inst] == cost[cf[i], inst]) eq++ }
      rk[inst, i] = lt + (eq + 1) / 2 } }
  function corr(a, b,   i, mean_a, mean_b, cab, caa, cbb) {
    for (i = 1; i <= nc; i++) { mean_a += rk[a, i]; mean_b += rk[b, i] } mean_a /= nc; mean_b /= nc
    for (i = 1; i <= nc; i++) { cab += (rk[a, i] - mean_a) * (rk[b, i] - mean_b); caa += (rk[a, i] - mean_a) ^ 2; cbb += (rk[b, i] - mean_b) ^ 2 }
    return (caa > 0 && cbb > 0) ? cab / sqrt(caa * cbb) : 0 }
  FILENAME == ARGV[1] { ni++; n[ni] = $1; path[ni] = $2; grp[ni] = $3; next }
  FILENAME == ARGV[2] { if ($7 == "ok") { ok[$3, $4] = 1; cost[$3, $4] = $8 + 0 } next }
  END {
    for (c = 1; c <= nconf; c++) { good = 1; for (i = 1; i <= ni; i++) if (!ok[c, i]) good = 0; if (good) cf[++nc] = c }
    printf "configurations with all %d runs ok: %d of %d\n", ni, nc, nconf
    if (nc < 5) { print "too few complete configurations for a correlation"; exit }
    for (i = 1; i <= ni; i++) rankvec(i)
    for (i = 1; i <= ni; i++) {                               # an instance on which all configurations tie carries no ranking information
      mn = 1e300; mx = -1e300
      for (k = 1; k <= nc; k++) { v = cost[cf[k], i]; if (v < mn) mn = v; if (v > mx) mx = v }
      flat[i] = (mx == mn); if (flat[i]) nflat++ }
    print "=== Spearman correlation of the rankings of the configurations, instance by instance (- : all configurations tie, left out) ==="
    printf "%-6s", ""; for (j = 1; j <= ni; j++) printf "%6s", "i" j; printf "\n"
    for (i = 1; i <= ni; i++) {
      printf "i%-2d %-3s", i, substr(grp[i], 1, 3)
      for (j = 1; j <= ni; j++) { if (flat[i] || flat[j]) { r[i, j] = 0; printf "%6s", "-" } else { r[i, j] = corr(i, j); printf "%6.2f", r[i, j] } }
      printf "  n=%d %s\n", n[i], path[i] }
    for (i = 1; i <= ni; i++) g[grp[i]] = 1
    for (i = 1; i <= ni; i++) for (j = i + 1; j <= ni; j++) {
      if (flat[i] || flat[j]) continue
      key = (grp[i] == grp[j]) ? "within " grp[i] : "between groups"; s[key] += r[i, j]; m[key]++ }
    print "=== mean pairwise correlation ==="
    for (key in s) printf "  %-22s %6.3f   (%d pairs)\n", key, s[key] / m[key], m[key]
    # mean rank per group and configuration; correlation of the group profiles
    printf "  instances on which all configurations tie: %d\n", nflat + 0
    for (i = 1; i <= ni; i++) { if (flat[i]) continue; gc[grp[i]]++; for (k = 1; k <= nc; k++) gs[grp[i], k] += rk[i, k] }
    ng = 0; for (x in g) { ng++; gname[ng] = x }
    if (ng == 2 && gc[gname[1]] > 0 && gc[gname[2]] > 0) {
      for (k = 1; k <= nc; k++) { ma[k] = gs[gname[1], k] / gc[gname[1]]; mb[k] = gs[gname[2], k] / gc[gname[2]] }
      xa = 0; xb = 0; for (k = 1; k <= nc; k++) { xa += ma[k]; xb += mb[k] } xa /= nc; xb /= nc
      for (k = 1; k <= nc; k++) { sab += (ma[k] - xa) * (mb[k] - xb); saa += (ma[k] - xa) ^ 2; sbb += (mb[k] - xb) ^ 2 }
      printf "correlation of the mean ranks (%s vs %s): %.3f (on average ranks, noisier than single instances)\n", gname[1], gname[2], sab / sqrt(saa * sbb)
      ba = 1; bb = 1; for (k = 2; k <= nc; k++) { if (ma[k] < ma[ba]) ba = k; if (mb[k] < mb[bb]) bb = k }
      printf "best configuration on %s: %d (mean rank %.1f on %s, of %d);  best on %s: %d (mean rank %.1f on %s)\n", gname[1], cf[ba], mb[ba], gname[2], nc, gname[2], cf[bb], ma[bb], gname[1]
    }
  }' "$OUT/instances_by_size.tsv" "$OUT/trlog.txt" | tee "$OUT/report.txt"
