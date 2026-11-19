# irace tuning of the VRP-Learn heuristics

## Design decisions
| item | decision |
|---|---|
| methods | the 18 tunable heuristics of `scripts/hpc/methods_main.txt` (68 parameters); no `-f` twins for ta, rtr, grasp, lns, aco, pso |
| split | hold-out: 28 training / 139 test instances (CVRP and MDVRP together, systematic sampling), `tuning/split/{train,test}.txt`, `split.tsv` |
| training budget | `maxExperiments = 2000` per method, `F = 24` s per 100 customers (the time budget of the main study) |
| cost | raw tour cost; irace's default Friedman test is rank-based, so different cost scales across instances are tolerable |
| comparison | tuned vs default parameters, reported **separately for training and test instances** |

## Pieces
| file | role |
|---|---|
| `target-runner` | one run of one method with one configuration on one instance; prints the cost. Hard time guard (`timeout -s KILL`, 3 x budget + 30 s), +1e9 penalty for infeasible / killed / crashed runs, best-known-value alarm (below a *proven* optimum: penalty), one line per run in `trlog.txt` (cost, solver seconds, iterations, wall time, limit) |
| `irace_run.R` | irace through its R API (`readScenario`, `checkIraceScenario`, `irace`) |
| `run_irace.sh` | one tuning run on the current machine: writes parameter file, default configuration, scenario; runs the data guard; summarises statuses and alarms; `--check` only checks the scenario |
| `run_irace.sge` | Grid Engine job around `run_irace.sh` (one job per method) |
| `fetch_data.sh`, `check_data.sh`, `data.md5` | instance download at pinned commits and the md5 guard (167 files) |
| `split/smoke.txt` | three small instances for the pipeline test |

## Run
```
cpp/build/vrp_experiments params --methods all --out /tmp/params      # parameter tables, 68 parameters
tuning/run_irace.sh ga-f --check -T tuning/split/smoke.txt -F 1        # irace executes the runner twice
tuning/run_irace.sh ga-f -E 500 -F 6 -P 4                              # pilot
qsub -N irace_ga-f -v METHOD=ga-f,E=2000,F=24 tuning/run_irace.sge     # cluster, one job per method
```

Output in `results/irace/METHOD/` (git-ignored): `irace.Rdata`, `final_elites.tsv`, `best_config.txt` (the `--param name=value` options of
`vrp_experiments`), `trlog.txt`, `alarms.txt`, `failed/`.

## Notes for interpretation
* Elite configurations are descriptive; the evidence for tuned vs default is the test-set comparison, not the training result.
* The training runs use the same wall-time budget as the main study, so tuning cost is `maxExperiments x F x mean(n)/100` core-seconds.
* A cost below an unproven best-known value only writes an alarm line; below a proven optimum it is treated as a bug.
