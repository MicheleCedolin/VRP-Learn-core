#!/usr/bin/env bash
# ONE irace run = the tuning of one method, on this machine (workstation, CI runner, or inside a cluster job).
# usage: tuning/run_irace.sh METHOD [-E max_experiments] [-F seconds_per_100_customers] [-P parallel] [-T train_list] [-S seed]
#                                   [-O outdir] [--check]
#   defaults: -E 2000  -F 24  -P 4  -T tuning/split/train.txt  -S 1  -O results/irace/METHOD
#   --check : only check the scenario (the target runner is executed twice), no tuning
# Needs: cpp/build (see README), the instance files (tuning/fetch_data.sh), R with the irace package.
# Writes OUTDIR/{params/, scenario.txt, irace.log, irace.Rdata, trlog.txt, alarms.txt, failed/, final_elites.tsv, best_config.txt}.
# Never submit this on a cluster login node (irace --check and the tuning execute the solver): use tuning/run_irace.sge.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."; ROOT=$PWD
METHOD=${1:?usage: tuning/run_irace.sh METHOD [-E n] [-F s] [-P n] [-T list] [-S seed] [-O dir] [--check]}; shift
E=2000; F=24; P=4; TRAIN=tuning/split/train.txt; SEED=1; OUT=""; CHECK=0
while [ $# -gt 0 ]; do
  case "$1" in
    -E) E=$2; shift 2;; -F) F=$2; shift 2;; -P) P=$2; shift 2;; -T) TRAIN=$2; shift 2;; -S) SEED=$2; shift 2;; -O) OUT=$2; shift 2;;
    --check) CHECK=1; shift;; *) echo "unknown option $1" >&2; exit 2;;
  esac
done
OUT=$(realpath -m "${OUT:-results/irace/$METHOD}"); TRAIN=$(realpath "$TRAIN")
[ -x cpp/build/vrp_experiments ] && [ -x cpp/build/vrp ] || { echo "build first: cmake -S cpp -B cpp/build -DCMAKE_BUILD_TYPE=Release && cmake --build cpp/build -j" >&2; exit 1; }
command -v Rscript >/dev/null || { echo "Rscript not found (R with the irace package is needed)" >&2; exit 1; }
Rscript -e 'if (!requireNamespace("irace", quietly = TRUE)) quit(status = 1)' || { echo "R package irace is not installed" >&2; exit 1; }
[ "${SKIP_DATA_GUARD:-0}" = 1 ] || tuning/check_data.sh
while read -r f; do [ -s "$f" ] || { echo "training instance missing: $f" >&2; exit 1; }; done < "$TRAIN"
mkdir -p "$OUT/params"
cpp/build/vrp_experiments params --methods "$METHOD" --out "$OUT/params" > "$OUT/params/table.txt"
[ -s "$OUT/params/$METHOD.txt" ] || { echo "$METHOD has no tunable parameters (see $OUT/params/table.txt)" >&2; exit 2; }
cat > "$OUT/scenario.txt" <<SCEN
targetRunner = "$ROOT/tuning/target-runner"
parameterFile = "$OUT/params/$METHOD.txt"
configurationsFile = "$OUT/params/${METHOD}_default.txt"
trainInstancesDir = ""
trainInstancesFile = "$TRAIN"
execDir = "$ROOT"
logFile = "$OUT/irace.Rdata"
maxExperiments = $E
parallel = $P
seed = $SEED
deterministic = FALSE
SCEN
export METHOD TIME_PER_100=$F VRP_ROOT=$ROOT TRLOG=$OUT/trlog.txt ALARMS=$OUT/alarms.txt FAILDIR=$OUT/failed
echo "$(date -Is) $METHOD: maxExperiments $E, $F s per 100 customers, parallel $P, seed $SEED, train list $(wc -l < "$TRAIN") instances, host $(hostname)"
if [ $CHECK = 1 ]; then Rscript tuning/irace_run.R "$OUT/scenario.txt" --check; exit $?; fi
rc=0; Rscript tuning/irace_run.R "$OUT/scenario.txt" 2>&1 | tee "$OUT/irace.log" || rc=$?
echo "$(date -Is) irace finished with status $rc"
if [ -s "$OUT/trlog.txt" ]; then
  echo "runs of the target runner by status (trlog.txt):"
  awk '{c[$7]++; n++} END{for (s in c) printf "  %-10s %d\n", s, c[s]; printf "  total      %d\n", n}' "$OUT/trlog.txt"
fi
[ -s "$OUT/alarms.txt" ] && { echo "ALARMS ($(wc -l < "$OUT/alarms.txt")), see $OUT/alarms.txt:"; head -3 "$OUT/alarms.txt"; }
if [ -s "$OUT/final_elites.tsv" ]; then                       # best configuration in the form of the --param options of vrp_experiments
  awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) h[i]=$i; next} NR==2{s=""; for(i=1;i<=NF;i++) if (h[i] !~ /^\./) s=s " --param " h[i] "=" $i; print substr(s,2)}' \
    "$OUT/final_elites.tsv" > "$OUT/best_config.txt"
  echo "best configuration: $(cat "$OUT/best_config.txt")"
fi
exit $rc
