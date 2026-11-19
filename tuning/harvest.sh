#!/usr/bin/env bash
# Turns the result bundle that the last job of a finished irace-tune run printed to its job log back into files.
# usage: tuning/harvest.sh RUN_ID METHOD [REPO]          (needs gh with read access to the repository)
# Writes tuning/results/METHOD/{final_elites.tsv, best_config.txt, scenario.txt, env.txt, irace.log, irace.Rdata, trlog.txt, alarms.txt, params/}.
# Why: GitHub deletes logs and artifacts after 90 days at the latest (public repositories), so the results are committed to the repository.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
RUN=${1:?usage: tuning/harvest.sh RUN_ID METHOD [REPO]}; METHOD=${2:?method}; REPO=${3:-MicheleCedolin/VRP-Learn-core}
JID=$(gh api "repos/$REPO/actions/runs/$RUN/jobs" --paginate --jq '.jobs[] | select(.name=="finish") | .id' | head -1)
[ -n "$JID" ] || { echo "no finish job in run $RUN (is the run completed?)" >&2; exit 1; }
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
gh api "repos/$REPO/actions/jobs/$JID/logs" > "$TMP/log.txt"
sed 's/^[0-9T:.-]*Z //' "$TMP/log.txt" | awk -v m="$METHOD" '
  $0=="-----BEGIN IRACE BUNDLE " m "-----" {on=1; next}
  $0=="-----END IRACE BUNDLE " m "-----"   {on=0}
  on' > "$TMP/b64.txt"
[ -s "$TMP/b64.txt" ] || { echo "no bundle for $METHOD in the log of job $JID" >&2; exit 1; }
base64 -d "$TMP/b64.txt" > "$TMP/bundle.tgz"
WANT=$(sed 's/^[0-9T:.-]*Z //' "$TMP/log.txt" | awk '/bundle.tgz$/ && NF==2 && length($1)==64 {print $1; exit}')
GOT=$(sha256sum "$TMP/bundle.tgz" | cut -d' ' -f1)
[ "$WANT" = "$GOT" ] || { echo "sha256 mismatch: log says $WANT, decoded $GOT" >&2; exit 1; }
mkdir -p tuning/results
rm -rf "tuning/results/$METHOD"
tar -C tuning/results -xzf "$TMP/bundle.tgz"
echo "wrote tuning/results/$METHOD ($(du -sh "tuning/results/$METHOD" | cut -f1)), sha256 ok"
ls -l "tuning/results/$METHOD"
