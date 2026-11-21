#!/usr/bin/env bash
# ONE irace run = the tuning of one method, on this machine (workstation, CI runner, or inside a cluster job).
# usage: tuning/run_irace.sh METHOD [-E max_experiments] [-F seconds_per_100_customers] [-P parallel] [-T train_list] [-S seed]
#                                   [-O outdir] [--check] [--recover] [--window MINUTES]
#   defaults: -E 2000  -F 24  -P 4  -T tuning/split/train.txt  -S 1  -O results/irace/METHOD
#   --check   : only check the scenario (the target runner is executed twice), no tuning
#   --recover : continue a run that was stopped (uses the last saved state OUTDIR/irace_good.Rdata or irace.Rdata)
#   --window  : stop irace after MINUTES (at a saved state if one is due, else at the limit) and exit with status 75 so that the
#               run can be continued with --recover in a later job (a hosted runner job is capped at 6 h).
# Needs: cpp/build (see README), the instance files (tuning/fetch_data.sh), R with the irace package.
# Writes OUTDIR/{params/, scenario.txt, env.txt, irace.log, irace.Rdata, trlog.txt, alarms.txt, failed/, final_elites.tsv, best_config.txt}.
# Never submit this on a cluster login node (irace --check and the tuning execute the solver): use tuning/run_irace.sge.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."; ROOT=$PWD
METHOD=${1:?usage: tuning/run_irace.sh METHOD [-E n] [-F s] [-P n] [-T list] [-S seed] [-O dir] [--check] [--recover] [--window min]}; shift
E=2000; F=24; P=4; TRAIN=tuning/split/train.txt; SEED=1; OUT=""; CHECK=0; RECOVER=0; WINDOW=0
while [ $# -gt 0 ]; do
  case "$1" in
    -E) E=$2; shift 2;; -F) F=$2; shift 2;; -P) P=$2; shift 2;; -T) TRAIN=$2; shift 2;; -S) SEED=$2; shift 2;; -O) OUT=$2; shift 2;;
    --check) CHECK=1; shift;; --recover) RECOVER=1; shift;; --window) WINDOW=$2; shift 2;;
    *) echo "unknown option $1" >&2; exit 2;;
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
# Wider ranges than the first proposals of the C++ table (tuning/ranges.tsv: method parameter type lower upper): the parameter file is
# rewritten for these rows only; the default configuration is checked to lie inside the new range.  Applied again in every stage.
if [ -s tuning/ranges.tsv ]; then
  python3 - "$METHOD" tuning/ranges.tsv "$OUT/params" <<'PY' | tee -a "$OUT/env.txt"
import sys
m, rf, pdir = sys.argv[1:4]
rows = {}
for ln in open(rf):
    if ln.startswith('#') or not ln.strip():
        continue
    a = ln.rstrip('\n').split('\t')
    if a[0] == m:
        rows[a[1]] = (a[2], a[3], a[4])
if rows:
    pf = f"{pdir}/{m}.txt"
    names, vals = [l.split() for l in open(f"{pdir}/{m}_default.txt")][:2]
    default = dict(zip(names, map(float, vals)))
    out = []
    for ln in open(pf):
        t = ln.split()                                # name "--name " type (lo, hi)
        if t and t[0] in rows:
            typ, slo, shi = rows.pop(t[0])
            lo, hi = float(slo), float(shi)
            assert lo <= default[t[0]] <= hi, f"{m}.{t[0]}: default {default[t[0]]} outside ({slo}, {shi})"
            assert 'log' not in typ or lo > 0, f"{m}.{t[0]}: log scale needs a positive lower bound"
            assert typ.split(',')[0] == t[3].split(',')[0], f"{m}.{t[0]}: type {typ} differs from the table type {t[3]}"
            print(f"[ranges] {m}.{t[0]}: {t[3]} {' '.join(t[4:])} -> {typ} ({slo}, {shi}) (default {default[t[0]]:g})")
            ln = f'{t[0]:<12} "--{t[0]} "  {typ}  ({slo}, {shi})\n'
        out.append(ln)
    assert not rows, f"tuning/ranges.tsv names unknown parameters of {m}: {sorted(rows)}"
    open(pf, 'w').writelines(out)
PY
fi
RECFILE=""
if [ $RECOVER = 1 ]; then                                     # continue from the last saved state of an earlier job
  for c in "$OUT/irace_good.Rdata" "$OUT/irace.Rdata"; do [ -s "$c" ] && { RECFILE=$c; break; }; done
  [ -n "$RECFILE" ] || { echo "--recover: no saved state in $OUT" >&2; exit 1; }
  cp "$RECFILE" "$OUT/irace_resume.Rdata"; RECFILE="$OUT/irace_resume.Rdata"
fi
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
[ -n "$RECFILE" ] && echo "recoveryFile = \"$RECFILE\"" >> "$OUT/scenario.txt"
export METHOD TIME_PER_100=$F VRP_ROOT=$ROOT TRLOG=$OUT/trlog.txt ALARMS=$OUT/alarms.txt FAILDIR=$OUT/failed
{ echo "== $(date -Is) host $(hostname), stage start (recover=$RECOVER window=$WINDOW)"
  lscpu 2>/dev/null | grep -E '^(Model name|Thread|Core|Socket|CPU\(s\)|CPU MHz)' || true
  echo "$(nproc) usable cores; $(uname -sr); $(free -g 2>/dev/null | awk '/Mem:/{print $2" GB RAM"}')"; } >> "$OUT/env.txt"
echo "$(date -Is) $METHOD: maxExperiments $E, $F s per 100 customers, parallel $P, seed $SEED, train list $(wc -l < "$TRAIN") instances, host $(hostname), recover $RECOVER, window $WINDOW min"
if [ $CHECK = 1 ]; then Rscript tuning/irace_run.R "$OUT/scenario.txt" --check; exit $?; fi

rc=0
if [ "$WINDOW" -le 0 ]; then
  [ $RECOVER = 1 ] || : > "$OUT/irace.log"
  Rscript tuning/irace_run.R "$OUT/scenario.txt" 2>&1 | tee -a "$OUT/irace.log" || rc=$?
else
  [ $RECOVER = 1 ] || : > "$OUT/irace.log"
  hard=$((WINDOW * 60)); t0=$SECONDS; t_prev=$SECONDS; iter_len=0; ncomp=0; armed=0
  [ $RECOVER = 1 ] || armed=1
  off=$(wc -l < "$OUT/irace.log")                              # this stage only reads the lines irace appends from here on
  # irace saves its state often but can only resume at the end of an iteration, so the stop is placed there.
  # An iteration ends with the line "Elite configurations"; in a resumed stage the replay of the old iterations is skipped (up to "Recovery completed").
  count_comp() { awk -v off="$off" -v rec="$RECOVER" 'BEGIN{n=0; seen=(rec==0)} NR>off { if ($0 ~ /Recovery completed/) {seen=1; n=0} else if ($0 ~ /UTC: Elite configurations/) n++ } END{print (seen ? n : -1)}' "$OUT/irace.log"; }
  setsid Rscript tuning/irace_run.R "$OUT/scenario.txt" >> "$OUT/irace.log" 2>&1 & ipid=$!
  tail -n 0 -F "$OUT/irace.log" 2>/dev/null & tpid=$!
  stopped=0
  while kill -0 "$ipid" 2>/dev/null; do
    sleep 15
    el=$((SECONDS - t0)); c=$(count_comp)
    if [ "$armed" = 0 ] && [ "$c" -ge 0 ]; then armed=1; t_prev=$SECONDS; fi
    if [ "$armed" = 1 ] && [ "$c" -gt "$ncomp" ]; then          # an iteration has just ended and irace has saved its state
      iter_len=$((SECONDS - t_prev)); t_prev=$SECONDS; ncomp=$c
      sleep 10; [ -s "$OUT/irace.Rdata" ] && { cp "$OUT/irace.Rdata" "$OUT/irace_good.Rdata.tmp" && mv "$OUT/irace_good.Rdata.tmp" "$OUT/irace_good.Rdata"; }
      echo "$(date -Is) [window] iteration done in this stage: #$c after $((el / 60)) min (this one took $((iter_len / 60)) min), state saved" | tee -a "$OUT/env.txt"
      # time the next iteration needs: 1.15 x the last one, scaled by (budget of the next iteration / experiments the last one used) when that is above 1
      # (the last iteration of a run often gets the left-over budget and is longer than the one before; a stop in the middle of an iteration loses its work)
      scale=$(awk '/^# experimentsUsed:/{a=b; b=$3} /^# currentBudget:/{cb=$3} END{ if (a != "" && b != "" && cb != "" && b > a) { r=cb/(b-a); if (r < 1) r=1; print r } else print 1 }' "$OUT/irace.log")
      need=$(awk -v n="$iter_len" -v r="$scale" 'BEGIN{printf "%d", n * 1.15 * r}')
      echo "[window] next iteration needs about $((need / 60)) min (scale $scale), $(( (hard - el) / 60 )) min left" | tee -a "$OUT/env.txt"
      if [ $((hard - el)) -lt "$need" ]; then stopped=1; echo "[window] no time left for another iteration: stopping at the saved state"; break; fi
    fi
    if [ $el -ge $hard ]; then stopped=1; echo "[window] limit of $WINDOW min reached: stopping (the unfinished iteration is lost)"; break; fi
  done
  if [ $stopped = 1 ] && kill -0 "$ipid" 2>/dev/null; then       # irace is still running: stop it, the next job continues
    kill -TERM -- "-$ipid" 2>/dev/null || true; sleep 5; kill -KILL -- "-$ipid" 2>/dev/null || true
    wait "$ipid" 2>/dev/null || true; rc=75
  else
    wait "$ipid" || rc=$?                                        # irace has ended by itself: its own status
  fi
  kill "$tpid" 2>/dev/null || true
fi
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
