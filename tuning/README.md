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
`vrp_experiments`), `trlog.txt`, `alarms.txt`, `env.txt` (CPU model of every stage), `failed/`.

## Tuning on GitHub-hosted runners (staged)
A hosted-runner job is capped at 6 h; one method needs about 5-9 h at F=24, E=2000 on 4 vCPU. The workflow `irace-tune` therefore runs a method
as a chain of stages (`irace-stage.yml`, up to four jobs): each stage stops itself at the end of an iteration of irace (the only point irace can
resume from; `run_irace.sh --window`), hands the saved state to the next stage as an artifact, and the next stage continues it (`--recover`).
Start: Actions > irace-tune > Run workflow (method, optional `next` list of methods to run afterwards, E, F, window). The last job prints the
results and a complete bundle (tar.gz, base64) to its log. Logs and artifacts are deleted after 90 days at the latest (public repositories), so
the results are brought into the repository: `tuning/harvest.sh RUN_ID METHOD` writes `tuning/results/METHOD/` (committed).

Ranges: the parameter files come from the C++ table (`vrp_experiments params`); `tuning/ranges.tsv` (method, parameter, type, lower, upper) replaces
the range of single parameters whose elites sat at a range edge in the first pass (2026-11-21). `run_irace.sh` rewrites those lines in every
stage, checks that the default lies inside, and prints the changes to `env.txt`. The first-pass results of ils-f, sa-f, ta and tabu-f are kept in
`tuning/results/METHOD-narrow/`; the wide-range runs write `tuning/results/METHOD/`. Several methods can run at the same time: start one
`irace-tune` run per chain with disjoint `next` lists (hosted runners are separate machines; the Free plan allows 20 concurrent jobs).

Hosted runners are not uniform: stages of one method can run on different CPUs (seen: Intel Xeon 6975P-C and AMD EPYC 7763), and a 4 vCPU runner
has 2 physical cores with 2 threads each, so four parallel solver runs share two cores. `env.txt` records the CPU of every stage and `trlog.txt`
the iterations reached per run. The target runner converts the time budget into an iteration budget per run (two short calibration runs), so the
budget is equal in wall time, not in iterations.

## Notes for interpretation
* Elite configurations are descriptive; the evidence for tuned vs default is the test-set comparison, not the training result.
* The training runs use the same wall-time budget as the main study, so tuning cost is `maxExperiments x F x mean(n)/100` core-seconds.
* A cost below an unproven best-known value only writes an alarm line; below a proven optimum it is treated as a bug.
