#!/usr/bin/env bash
# Data guard: every instance file of the study is present and byte-identical to the files the tuning was prepared with.
# Exit status 0 = all 167 files match tuning/data.md5.  Run it on the cluster before submitting (a different file there would
# mean that the tuning and the main study use different instances).
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
n=$(wc -l < tuning/data.md5)
bad=$(md5sum -c tuning/data.md5 2>&1 | grep -v ': OK$' || true)
if [ -n "$bad" ]; then echo "DATA GUARD FAILED:"; echo "$bad" | head -20; exit 1; fi
echo "data guard: $n instance files match tuning/data.md5"
