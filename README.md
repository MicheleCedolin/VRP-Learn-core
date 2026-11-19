# VRP-Learn core (irace tuning pipeline)

C++ core of the VRP-Learn comparison study (CVRP and MDVRP, 20 heuristics, 167 benchmark instances) together with the
[irace](https://mlopez-ibanez.github.io/irace/) parameter-tuning set-up for its 18 tunable heuristics (68 parameters).

This repository is a **reduced extract** of the full project: the solver library and the two command-line programs, the unit tests
that do not need frozen reference tables, the best-known-solution table and the tuning scripts.  The browser (JS) code, the frozen
reference tables and the remaining unit tests are not part of it.

```
cpp/                C++17 library (vrplearn), `vrp` and `vrp_experiments`, unit tests
data/bks/           best-known solutions (X, CMT, Golden; proven_optimal flag)
data/catalogue.tsv  instance catalogue read by the programs
data/external/      instance files, NOT in git (tuning/fetch_data.sh downloads them; Cordeau files are GPL-2.0)
scripts/hpc/        methods_main.txt, the 20 methods of the main study
scripts/irace/      the original irace target runner
tuning/             guarded target runner, run scripts, train/test split, data guard  -> tuning/README.md
.github/workflows/  irace-smoke: pipeline smoke test (and ga-f pilot) on a hosted runner
```

Build and test:

```
cmake -S cpp -B cpp/build -DCMAKE_BUILD_TYPE=Release && cmake --build cpp/build -j
ctest --test-dir cpp/build --output-on-failure
tuning/fetch_data.sh
```

Licence: MIT (code).  Benchmark instances keep the licences of their sources (see `tuning/data_sources.txt`).
