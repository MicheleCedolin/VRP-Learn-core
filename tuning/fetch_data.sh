#!/usr/bin/env bash
# The 167 instance files of the study into data/external/ (X, CMT, Golden: vidalt/HGS-CVRP; Cordeau: STOM-Group/MDVRP-Instances), taken at
# the commits recorded in tuning/data_sources.txt, then the check against tuning/data.md5.  Needs git and internet access (login node or CI runner).
# The Cordeau files are GPL-2.0 and stay out of this repository: they are downloaded, never committed (data/external/ is git-ignored).
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
HGS_SHA=1a927955cd2861a29d978f0d359d6e647db9319c
MDVRP_SHA=4a7738ccc11d872520e3820a9639c23b4026dda9
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
mkdir -p data/external/cvrp data/external/mdvrp
fetch() {   # fetch REPO_URL SHA DIR : shallow checkout of exactly that commit
  git init -q "$3" && git -C "$3" remote add origin "$1" && git -C "$3" fetch -q --depth 1 origin "$2" && git -C "$3" checkout -q FETCH_HEAD
}
fetch https://github.com/vidalt/HGS-CVRP.git "$HGS_SHA" "$TMP/hgs"
fetch https://github.com/STOM-Group/MDVRP-Instances.git "$MDVRP_SHA" "$TMP/mdvrp"
cp "$TMP"/hgs/Instances/CVRP/X-n*.vrp "$TMP"/hgs/Instances/CVRP/CMT*.vrp "$TMP"/hgs/Instances/CVRP/Golden_*.vrp data/external/cvrp/
for n in p01 p02 p03 p04 p05 p06 p07 p08 p09 p10 p11 p12 p13 p14 p15 p16 p17 p18 p19 p20 p21 p22 p23 pr01 pr02 pr03 pr04 pr05 pr06 pr07 pr08 pr09 pr10; do
  cp "$TMP/mdvrp/dat/$n" data/external/mdvrp/
done
exec tuning/check_data.sh
