// Reproducible experiments.  Every run is seeded; results go to results/*.csv.
//
//   vrp_experiments construct
//   vrp_experiments ops
//   vrp_experiments depots [--iters 8000 --seeds 3]
//   vrp_experiments meta   [--iters 8000 --seeds 3]
//   vrp_experiments methods [--time 2 --seeds 3]    all twenty methods with an equal time budget (iterations calibrated per method)
//   vrp_experiments bench --instances ... [--methods all] [--time 10 --seeds 10]    all twenty methods at equal time, for cluster runs
//   options: --out DIR (default results)   --only NAME[,NAME...] (restrict instances)
//
// Instances must be present in data/external (scripts/fetch_instances.sh).  Timings are CPU seconds of
// whatever machine runs it.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <cctype>
#include <iostream>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <variant>
#include <vector>

#include "vrplearn/benchmarks.hpp"
#include "vrplearn/construct.hpp"
#include "vrplearn/localsearch.hpp"
#include "vrplearn/methods.hpp"
#include "vrplearn/metaheuristics.hpp"

using namespace vrp;

namespace {

const std::vector<std::string> CVRP = {"X-n101-k25", "X-n153-k22", "X-n261-k13"};
const std::vector<std::string> MD = {"p01", "p02", "p03", "p04", "p05", "p06", "p07"};

std::set<std::string> only;                  // empty = all
std::string outDir = "results";

bool wanted(const std::string& name) { return only.empty() || only.count(name); }

double now() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }

// ---- CSV values formatted like Python's csv module / float repr -----------------------------------
double roundTo(double x, int nd) {
    char buf[64];
    std::snprintf(buf, sizeof buf, "%.*f", nd, x);
    return std::strtod(buf, nullptr);
}

std::string pyFloat(double x) {                                    // repr(float) for ordinary magnitudes
    char buf[64];                                                   // shortest text that reads back as x (no <charconv>: old compilers)
    int p = 1;
    for (; p <= 17; p++) { std::snprintf(buf, sizeof buf, "%.*e", p - 1, x); if (std::strtod(buf, nullptr) == x) break; }
    const double ax = std::fabs(x);
    if (ax == 0 || (ax >= 1e-4 && ax < 1e16)) {                      // fixed notation with exactly the significant digits found
        const int e10 = ax == 0 ? 0 : (int)std::floor(std::log10(ax));
        for (int d = std::max(0, p - 1 - e10);; d++) { std::snprintf(buf, sizeof buf, "%.*f", d, x); if (std::strtod(buf, nullptr) == x || d > 30) break; }
    }
    std::string s(buf);
    if (s.find_first_of(".eEn") == std::string::npos) s += ".0";
    return s;
}

using Cell = std::variant<std::string, long, double, bool>;

std::string fmt(const Cell& c) {
    if (auto s = std::get_if<std::string>(&c)) return *s;
    if (auto l = std::get_if<long>(&c)) return std::to_string(*l);
    if (auto b = std::get_if<bool>(&c)) return *b ? "True" : "False";
    return pyFloat(std::get<double>(c));
}

struct Table {
    std::vector<std::string> cols;
    std::vector<std::vector<Cell>> rows;
    void add(std::vector<Cell> r) { rows.push_back(std::move(r)); }
    void write(const std::string& name) const {
        mkdir(outDir.c_str(), 0755);
        std::string path = outDir + "/" + name;
        std::ofstream f(path);
        for (size_t i = 0; i < cols.size(); i++) f << (i ? "," : "") << cols[i];
        f << "\r\n";                                                // Python's csv module uses \r\n
        for (const auto& r : rows) {
            for (size_t i = 0; i < r.size(); i++) f << (i ? "," : "") << fmt(r[i]);
            f << "\r\n";
        }
        std::printf("wrote %s\n", path.c_str());
    }
};

Cell gapCell(double c, const std::string& name) {
    auto g = gapPct(c, name);
    if (!g) return std::string();
    return roundTo(*g, 3);
}

double checkedCost(const Instance& inst, const Solution& sol) {
    auto v = sol.violations(inst);
    if (!v.empty()) { std::fprintf(stderr, "infeasible solution: %s\n", v[0].c_str()); std::exit(1); }
    return sol.cost(inst);
}

// ---------------------------------------------------------------------------------------------------
void expConstruct() {
    Table T{{"instance", "kind", "method", "cost", "routes", "gap_bks_pct", "seconds"}, {}};
    std::vector<std::string> names = CVRP;
    names.insert(names.end(), MD.begin(), MD.end());
    for (const auto& name : names) {
        if (!wanted(name)) continue;
        Instance inst = loadBenchmark(name);
        struct M { const char* label; Solution (*fn)(const Instance&); };
        const M methods[] = {
            {"clarke_wright", [](const Instance& i) { return savings(i); }},
            {"sweep", [](const Instance& i) { return sweep(i); }},
            {"nearest_neighbor", [](const Instance& i) { return nearestNeighbor(i); }},
            {"cheapest_insertion", [](const Instance& i) { return cheapestInsertion(i); }},
        };
        for (const M& m : methods) {
            double t0 = now();
            Solution sol = m.fn(inst);
            double dt = now() - t0;
            double c = checkedCost(inst, sol);
            T.add({name, findEntry(name)->kind, std::string(m.label), roundTo(c, 2), (long)sol.nRoutes(), gapCell(c, name), roundTo(dt, 3)});
            std::printf("%s %s %.2f\n", name.c_str(), m.label, c);
        }
    }
    T.write("construct.csv");
}

// ---------------------------------------------------------------------------------------------------
const std::vector<std::string> FULL = {"relocate", "swap", "two_opt", "two_opt_star", "or_opt", "swap_star", "route_depot", "depot_swap"};

void expOps() {
    Table T{{"instance", "start", "config", "start_cost", "cost", "improvement_pct", "routes", "gap_bks_pct", "seconds", "applied"}, {}};
    std::vector<std::string> names = {CVRP[0], CVRP[1]};
    names.insert(names.end(), MD.begin(), MD.end());
    for (const auto& name : names) {
        if (!wanted(name)) continue;
        Instance inst = loadBenchmark(name);
        LSContext ctx(inst, 20);
        struct Start { const char* label; Solution sol; };
        std::vector<Start> starts = {{"savings", savingsMd(inst)}, {"nearest_neighbor", nearestNeighbor(inst)}};
        for (const Start& st : starts) {
            double c0 = checkedCost(inst, st.sol);
            struct Cfg { std::string name; std::string preset; std::vector<std::string> ops; bool presetOps; };
            std::vector<Cfg> configs;
            for (const char* o : {"relocate", "swap", "two_opt", "or_opt", "two_opt_star", "swap_star"})
                configs.push_back({o, "full", {o}, false});
            for (const char* p : {"basic", "hgs", "full"}) configs.push_back({p, p, {}, true});
            if (inst.t > 1) configs.push_back({"route_depot+depot_swap", "full", {"route_depot", "depot_swap"}, false});
            for (const auto& op : FULL) {
                if (inst.t == 1 && (op == "route_depot" || op == "depot_swap")) continue;
                std::vector<std::string> rest;
                for (const auto& o : FULL) if (o != op) rest.push_back(o);
                configs.push_back({"full minus " + op, "full", rest, false});
            }
            for (const Cfg& cf : configs) {
                double t0 = now();
                LSResult r = localSearch(inst, st.sol, cf.preset, cf.presetOps ? nullptr : &cf.ops, 20, 0, &ctx);
                double dt = now() - t0;
                double c = checkedCost(inst, r.sol);
                long applied = 0;
                for (long s : r.stats) applied += s;
                T.add({name, std::string(st.label), cf.name, roundTo(c0, 2), roundTo(c, 2), roundTo(100 * (c0 - c) / c0, 3),
                       (long)r.sol.nRoutes(), gapCell(c, name), roundTo(dt, 3), applied});
            }
        }
        std::printf("%s done\n", name.c_str());
    }
    T.write("ops.csv");
}

// ---------------------------------------------------------------------------------------------------
Instance subinstance(const Instance& inst, const std::vector<int>& cust, int d) {
    Instance s;
    s.name = inst.name + "-d" + std::to_string(d);
    s.n = (int)cust.size();
    s.t = 1;
    std::vector<int> idx = cust;
    idx.push_back(inst.n + d);
    for (int i : idx) { s.x.push_back(inst.x[i]); s.y.push_back(inst.y[i]); }
    for (int c : cust) { s.demand.push_back(inst.demand[c]); s.service.push_back(inst.service[c]); }
    s.capacity = {inst.capacity[d]};
    s.maxDuration = {inst.maxDuration[d]};
    s.fleet = {inst.fleet[d]};
    s.rounded = inst.rounded;
    const int N = (int)idx.size();
    s.dist.resize(static_cast<size_t>(N) * N);
    for (int a = 0; a < N; a++) for (int b = 0; b < N; b++) s.dist[static_cast<size_t>(a) * N + b] = inst.d(idx[a], idx[b]);
    return s;
}

// Cluster first / route second: every customer goes to its nearest depot (fixed), then t independent CVRPs.
Solution assignThenRoute(const Instance& inst, int iters, unsigned long long seed) {
    std::vector<std::vector<int>> groups(inst.t);
    for (int c = 0; c < inst.n; c++) {
        int d = 0;
        for (int x = 1; x < inst.t; x++) if (inst.d(c, inst.n + x) < inst.d(c, inst.n + d)) d = x;
        groups[d].push_back(c);
    }
    Solution out;
    for (int d = 0; d < inst.t; d++) {
        const auto& cust = groups[d];
        if (cust.empty()) continue;
        Instance sub = subinstance(inst, cust, d);
        SisrOptions o;
        o.iters = std::max(200, (int)((double)(iters * (long)cust.size()) / inst.n));
        o.seed = seed;
        Result r = sisr(sub, o);
        for (const auto& rt : r.best.routes) {
            std::vector<int> mapped;
            for (int v : rt) mapped.push_back(cust[v]);
            out.routes.push_back(mapped);
            out.depots.push_back(d);
        }
    }
    return out;
}

void expDepots(int iters, int seeds) {
    Table T{{"instance", "method", "seed", "cost", "routes", "feasible", "seconds", "gap_bks_pct"}, {}};
    for (const auto& name : MD) {
        if (!wanted(name)) continue;
        Instance inst = loadBenchmark(name);
        for (int seed = 0; seed < seeds; seed++) {
            double t0 = now();
            Solution sol = assignThenRoute(inst, iters, seed);
            double dt = now() - t0;
            bool ok = sol.violations(inst).empty();
            double c = sol.cost(inst);
            T.add({name, std::string("assign-nearest-depot then SISR per depot"), (long)seed, roundTo(c, 2), (long)sol.nRoutes(), ok,
                   roundTo(dt, 2), gapCell(c, name)});
            t0 = now();
            SisrOptions o;
            o.iters = iters;
            o.seed = seed;
            Result r = sisr(inst, o);
            double c2 = checkedCost(inst, r.best);
            T.add({name, std::string("joint SISR (customers may switch depot)"), (long)seed, roundTo(c2, 2), (long)r.best.nRoutes(), true,
                   roundTo(now() - t0, 2), gapCell(c2, name)});
            std::printf("%.2f %.2f %s %d\n", c, c2, name.c_str(), seed);
        }
    }
    T.write("depots.csv");
}

// ---------------------------------------------------------------------------------------------------
void expMeta(int iters, int seeds) {
    Table T{{"instance", "method", "seed", "iters", "cost", "gap_bks_pct", "seconds"}, {}};
    std::vector<std::string> names = {CVRP[0], CVRP[1]};
    names.insert(names.end(), MD.begin(), MD.end());
    for (const auto& name : names) {
        if (!wanted(name)) continue;
        Instance inst = loadBenchmark(name);
        for (int seed = 0; seed < seeds; seed++) {
            double t0 = now();
            SisrOptions so;
            so.iters = iters;
            so.seed = seed;
            Result r = sisr(inst, so);
            double c = checkedCost(inst, r.best);
            T.add({name, std::string("SISR-style"), (long)seed, (long)iters, roundTo(c, 2), gapCell(c, name), roundTo(now() - t0, 2)});
            t0 = now();
            IlsOptions io;
            io.iters = std::max(30, iters / 40);
            io.seed = seed;
            Result r2 = ils(inst, io);
            double c2 = checkedCost(inst, r2.best);
            T.add({name, std::string("ILS (R&R + local search)"), (long)seed, (long)io.iters, roundTo(c2, 2), gapCell(c2, name),
                   roundTo(now() - t0, 2)});
            std::printf("%s %d %.2f %.2f\n", name.c_str(), seed, c, c2);
        }
    }
    T.write("meta.csv");
}

// Equal-time comparison: a short timing probe sets the iteration budget of each method (the schedules of SA, ALNS and SISR depend on
// the number of iterations), then every method runs that many iterations with the same seeds.  The iteration budget is written to the CSV.
int calibrated(const std::string& id, const Instance& inst, double seconds) {
    const double probe = std::min(0.4, std::max(0.1, seconds / 5));
    double it = 0, sec = 0;
    if (id == "sisr") {
        SisrOptions o; o.iters = 1000000; o.timeLimit = probe; Result r = sisr(inst, o); it = r.iters; sec = r.seconds;
    } else if (id == "ils") {
        IlsOptions o; o.iters = 1000000; o.timeLimit = probe; o.preset = inst.t > 1 ? "full" : "hgs"; Result r = ils(inst, o); it = r.iters; sec = r.seconds;
    } else {
        MethodOptions o; o.iters = 100000000; o.timeLimit = probe; Result r = runMethod(id, inst, o); it = r.iters; sec = r.seconds;
    }
    return std::max(10, (int)(it / std::max(sec, 1e-3) * seconds));
}

void expMethods(double seconds, int seeds) {
    Table T{{"instance", "method", "seed", "budget_s", "iters", "cost", "gap_bks_pct", "seconds"}, {}};
    std::vector<std::string> names = {CVRP[0], CVRP[1]};
    names.insert(names.end(), MD.begin(), MD.end());
    std::vector<std::string> ids = {"sisr", "ils"};
    for (const auto& m : methodIds()) ids.push_back(m);
    for (const auto& name : names) {
        if (!wanted(name)) continue;
        Instance inst = loadBenchmark(name);
        for (const auto& id : ids) {
            const int iters = calibrated(id, inst, seconds);
            for (int seed = 0; seed < seeds; seed++) {
                const double t0 = now();
                Result r;
                if (id == "sisr") { SisrOptions o; o.iters = iters; o.seed = seed; r = sisr(inst, o); }
                else if (id == "ils") { IlsOptions o; o.iters = iters; o.seed = seed; o.preset = inst.t > 1 ? "full" : "hgs"; r = ils(inst, o); }
                else { MethodOptions o; o.iters = iters; o.seed = seed; r = runMethod(id, inst, o); }
                const double c = checkedCost(inst, r.best);
                T.add({name, id, (long)seed, roundTo(seconds, 2), (long)iters, roundTo(c, 2), gapCell(c, name), roundTo(now() - t0, 2)});
                std::printf("%s %s %d %.2f\n", name.c_str(), id.c_str(), seed, c);
            }
        }
    }
    T.write("methods.csv");
}

// ---------------------------------------------------------------------------------------------------
// bench: equal-time comparison of all twenty methods, written for cluster runs (see scripts/hpc/).
//   The time budget is turned into an iteration budget per (instance, method) by two short timing runs, because the
//   cooling / adaptation schedules of several methods depend on the number of iterations.  Budgets can be computed once
//   (--calibrate-only) and reused (--budgets FILE) so that every seed and every node runs the same amount of work.
//   Rows are appended to <out>/bench_<tag>.csv as soon as a run finishes, so an interrupted job keeps its results.
const std::vector<std::string>& allMethods() {
    static std::vector<std::string> ids;
    if (ids.empty()) { ids = {"sisr", "ils"}; for (const auto& m : methodIds()) ids.push_back(m); for (const auto& m : fullMethodIds()) ids.push_back(m); }
    return ids;
}

// target > 0: stop as soon as a feasible solution reaches it (pass proven optima only).  init: common starting solution.
std::map<std::string, double> tuned;                               // --param name=value: tunable parameters of the method (methodParams)
OpSet excludedOps = 0;                                             // --without: operators switched off in every descent (ablation)
std::string excludedNames;

Result runAny(const std::string& id, const Instance& inst, int iters, unsigned long long seed, double timeLimit, double target = 0,
              const Solution* init = nullptr) {
    if (id == "ils" && excludedOps) {                                // the compact ILS takes its operators as a list of names
        IlsOptions o; o.iters = iters; o.seed = seed; o.timeLimit = timeLimit; o.target = target; o.init = init; o.preset = inst.t > 1 ? "full" : "hgs";
        std::vector<std::string> names;
        for (int k = 0; k < OP_NEW_ROUTE; k++) if ((preset(o.preset).ops & ~excludedOps) >> k & 1) names.push_back(opName(k));
        o.ops = &names;
        return ils(inst, o);
    }
    if (id == "sisr") {
        SisrOptions o; o.iters = iters; o.seed = seed; o.timeLimit = timeLimit; o.target = target; o.init = init;
        for (const auto& kv : tuned) {
            if (kv.first == "cbar") o.cbar = kv.second; else if (kv.first == "lmax") o.lmax = (int)std::lround(kv.second); else if (kv.first == "alpha") o.alpha = kv.second;
            else if (kv.first == "blink") o.blink = kv.second; else if (kv.first == "t0") o.t0Factor = kv.second; else if (kv.first == "tf") o.tfFactor = kv.second;
            else throw std::runtime_error("method sisr has no parameter " + kv.first);
        }
        return sisr(inst, o);
    }
    if (id == "ils" && !tuned.empty()) throw std::runtime_error("method ils has no tunable parameters in this driver");
    if (id == "ils") { IlsOptions o; o.iters = iters; o.seed = seed; o.timeLimit = timeLimit; o.target = target; o.init = init; o.preset = inst.t > 1 ? "full" : "hgs"; return ils(inst, o); }
    MethodOptions o; o.iters = iters; o.seed = seed; o.timeLimit = timeLimit; o.target = target; o.init = init; o.excludeOps = excludedOps; o.params = tuned;
    return runMethod(id, inst, o);
}

// per-operator rows of a run: moves applied and total cost decrease in the descents
void writeOps(std::ostream& f, const std::string& prefix, const Result& r) {
    char b[160];
    for (int k = 0; k < OP_COUNT; k++)
        if (r.opStats[k]) { std::snprintf(b, sizeof b, "%s,%s,%ld,%.4f\n", prefix.c_str(), opName(k), r.opStats[k], r.opGain[k]); f << b; }
    for (const NamedStat& e : r.extra) { std::snprintf(b, sizeof b, "%s,%s,%ld,%.4f\n", prefix.c_str(), e.name.c_str(), e.count, e.value); f << b; }
}

// --resume: drop an incomplete last line (a job that was killed while writing) so that appended rows start on a new line
void trimPartialLine(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) return;
    std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    in.close();
    if (text.empty() || text.back() == '\n') return;
    const size_t cut = text.find_last_of('\n');
    std::ofstream out(path, std::ios::binary | std::ios::trunc);
    out << (cut == std::string::npos ? std::string() : text.substr(0, cut + 1));
}

// Common starting solution of every method.  The savings solution can use more vehicles than a depot owns on tight-fleet
// instances; methods that accept moves by cost only would then finish fleet-infeasible.  In that case the start is repaired
// once, deterministically, by a short penalised ruin & recreate run (fixed seed), and all methods start from the repaired
// solution.  Returns the savings solution unchanged when it is feasible (the usual case).
Solution commonStart(const Instance& inst, bool repair, bool* repaired) {
    Solution s = savingsMd(inst);
    *repaired = false;
    if (!repair || s.violations(inst).empty()) return s;
    for (int iters : {20000, 80000, 320000}) {
        SisrOptions o; o.iters = iters; o.seed = 987654321ULL; o.init = &s;
        Result r = sisr(inst, o);
        if (r.best.violations(inst).empty()) { *repaired = true; return r.best; }
    }
    return s;
}

// iterations that take about `seconds`: the rate is measured between a short and a longer run of the same trajectory, so the
// start-up of population methods (initial population, first descents) is not mistaken for the steady-state cost of an iteration
int calibrateIters(const std::string& id, const Instance& inst, double seconds) {
    const double p = std::min(0.5, std::max(0.1, seconds / 8));
    const int huge = 1000000000;
    Result a = runAny(id, inst, huge, 12345, p), b = runAny(id, inst, huge, 12345, 3 * p);
    double rate;                                                   // iterations per second
    if (b.iters > a.iters && b.seconds > a.seconds + 1e-6) rate = (b.iters - a.iters) / (b.seconds - a.seconds);
    else rate = b.iters / std::max(b.seconds, 1e-6);
    const double it = a.seconds < seconds ? a.iters + rate * (seconds - a.seconds) : rate * seconds;
    return (int)std::max(10.0, std::min(1.9e9, it));
}

std::string baseName(const std::string& path) {
    std::string s = path.substr(path.find_last_of("/\\") + 1);
    const size_t dot = s.rfind('.');
    if (dot != std::string::npos && dot > 0) s = s.substr(0, dot);
    return s;
}
bool fileExists(const std::string& p) { struct stat st; return stat(p.c_str(), &st) == 0; }
std::vector<std::string> splitList(const std::string& arg) {       // a,b,c  or  @file (one entry per line, # comments)
    std::vector<std::string> out;
    if (!arg.empty() && arg[0] == '@') {
        std::ifstream f(arg.substr(1));
        if (!f) throw std::runtime_error("cannot read " + arg.substr(1));
        std::string line;
        while (std::getline(f, line)) {
            const size_t a = line.find_first_not_of(" \t\r"), b = line.find_last_not_of(" \t\r");
            if (a == std::string::npos || line[a] == '#') continue;
            out.push_back(line.substr(a, b - a + 1));
        }
    } else {
        std::stringstream ss(arg);
        std::string t;
        while (std::getline(ss, t, ',')) if (!t.empty()) out.push_back(t);
    }
    return out;
}

struct BenchArgs {
    std::string instances, methods = "all", budgets, bks, tag = "run";
    double seconds = 10, per100 = 0;                               // per100 > 0: budget = customers x per100 / 100 seconds
    int seeds = 10, seed0 = 0;
    bool calibrateOnly = false;
    bool earlyStop = true;                                         // stop a run when it reaches a proven optimum (single depot)
    bool feasibleStart = true;                                     // repair a fleet-infeasible savings start once, for all methods
    bool traces = true;                                            // write every improvement with its timestamp to trace_<tag>.csv
    bool resume = false;                                           // skip runs that are already in bench_<tag>.csv
    bool routes = true;                                            // write the routes of every final solution to routes_<tag>.csv
    bool ops = true;                                               // write the local-search operator statistics to ops_<tag>.csv
};

int expBench(const BenchArgs& A) {
    if (A.instances.empty()) throw std::runtime_error("bench needs --instances NAME[,NAME...] | PATH[,PATH...] | @FILE");
    // all: the twenty compact methods; full: the upgraded variants; everything: both
    std::vector<std::string> ids = A.methods == "everything" ? allMethods() : A.methods == "full" ? fullMethodIds() : splitList(A.methods);
    if (A.methods == "all") { ids = {"sisr", "ils"}; for (const auto& m : methodIds()) ids.push_back(m); }
    for (const auto& id : ids)
        if (std::find(allMethods().begin(), allMethods().end(), id) == allMethods().end()) throw std::runtime_error("unknown method " + id);

    std::map<std::string, int> budget;                             // "instance,method" -> iterations
    if (!A.budgets.empty()) {
        std::ifstream f(A.budgets);
        if (!f) throw std::runtime_error("cannot read " + A.budgets);
        std::string line;
        std::getline(f, line);                                     // header
        while (std::getline(f, line)) {
            std::stringstream ss(line);
            std::string inst, m, sec, it;
            if (std::getline(ss, inst, ',') && std::getline(ss, m, ',') && std::getline(ss, sec, ',') && std::getline(ss, it, ',')) budget[inst + "," + m] = std::stoi(it);
        }
    }
    // Reference values: a CSV or TSV table with a header row.  Columns are found by name: instance | name;
    // bks | ub | upper_bound | best_known | value;  optional proven_optimal | optimal | proven | opt (yes / true / 1).
    // Without a recognisable header the first two columns are instance and value.
    std::map<std::string, double> bks;
    std::map<std::string, bool> proven;
    if (!A.bks.empty()) {
        std::ifstream f(A.bks);
        if (!f) throw std::runtime_error("cannot read " + A.bks);
        std::string line;
        int cName = 0, cVal = 1, cOpt = -1;
        bool first = true;
        while (std::getline(f, line)) {
            if (!line.empty() && line.back() == '\r') line.pop_back();
            if (line.empty() || line[0] == '#') continue;
            const char sep = line.find('\t') != std::string::npos ? '\t' : ',';
            std::vector<std::string> cell;
            std::stringstream ss(line);
            std::string t;
            while (std::getline(ss, t, sep)) cell.push_back(t);
            auto lower = [](std::string x) { for (char& ch : x) ch = (char)std::tolower((unsigned char)ch); return x; };
            if (first) {
                first = false;
                bool header = false;
                for (int i = 0; i < (int)cell.size(); i++) {
                    const std::string h = lower(cell[i]);
                    if (h == "instance" || h == "name") { cName = i; header = true; }
                    else if (h == "bks" || h == "ub" || h == "upper_bound" || h == "best_known" || h == "value") { cVal = i; header = true; }
                    else if (h == "proven_optimal" || h == "optimal" || h == "proven" || h == "opt") { cOpt = i; header = true; }
                }
                if (header) continue;
            }
            if ((int)cell.size() <= std::max(cName, cVal)) continue;
            try { bks[cell[cName]] = std::stod(cell[cVal]); } catch (const std::exception&) { continue; }
            if (cOpt >= 0 && cOpt < (int)cell.size()) { const std::string v = lower(cell[cOpt]); proven[cell[cName]] = v == "yes" || v == "true" || v == "1" || v == "y"; }
        }
    }

    mkdir(outDir.c_str(), 0755);
    const std::string path = outDir + "/" + (A.calibrateOnly ? "budgets_" : "bench_") + A.tag + ".csv";
    std::set<std::string> done;                                    // instance,method,seed,without of the finished runs
    if (A.resume && !A.calibrateOnly) {
        for (const char* pre : {"bench_", "trace_", "routes_", "ops_"}) trimPartialLine(outDir + "/" + pre + A.tag + ".csv");
        std::ifstream f(path);
        std::string line;
        std::getline(f, line);
        while (std::getline(f, line)) {
            std::vector<std::string> c;
            std::stringstream ss(line);
            std::string t;
            while (std::getline(ss, t, ',')) c.push_back(t);
            if (c.size() >= 28) done.insert(c[0] + "," + c[3] + "," + c[4] + "," + (c.size() > 28 ? c[28] : std::string()));
        }
    }
    const bool fresh = !fileExists(path);
    std::ofstream out(path, std::ios::app);
    if (!out) throw std::runtime_error("cannot write " + path);
    if (fresh) out << (A.calibrateOnly ? "instance,method,budget_s,iters\n"
                                       : "instance,n,depots,method,seed,budget_s,iters,cost,best_known,gap_pct,routes,seconds,seconds_to_best,iters_to_best,feasible,"
                                         "c01,c02,c05,c10,c15,c20,c30,c50,c75,c100,proven_optimal,reached_proven_optimum,start_repaired,without\n");
    // timestamped improvement log: one row per new best solution of a run (iteration, seconds since the start of the run, cost)
    std::ofstream trace;
    if (A.traces && !A.calibrateOnly) {
        const std::string tp = outDir + "/trace_" + A.tag + ".csv";
        const bool tfresh = !fileExists(tp);
        trace.open(tp, std::ios::app);
        if (!trace) throw std::runtime_error("cannot write " + tp);
        if (tfresh) trace << "instance,method,seed,iter,seconds,cost\n";
    }
    std::ofstream routesOut, opsOut;
    if (A.routes && !A.calibrateOnly) {
        const std::string rp = outDir + "/routes_" + A.tag + ".csv";
        const bool rfresh = !fileExists(rp);
        routesOut.open(rp, std::ios::app);
        if (!routesOut) throw std::runtime_error("cannot write " + rp);
        if (rfresh) routesOut << "instance,method,seed,without,cost,route,depot,load,customers\n";
    }
    if (A.ops && !A.calibrateOnly) {
        const std::string op = outDir + "/ops_" + A.tag + ".csv";
        const bool ofresh = !fileExists(op);
        opsOut.open(op, std::ios::app);
        if (!opsOut) throw std::runtime_error("cannot write " + op);
        if (ofresh) opsOut << "instance,method,seed,without,operator,moves,gain\n";
    }
    char buf[512];
    int bad = 0;
    for (const std::string& spec : splitList(A.instances)) {
        const bool isFile = fileExists(spec);
        const std::string name = isFile ? baseName(spec) : spec;
        Instance inst = isFile ? readInstance(spec) : loadBenchmark(spec);
        std::optional<double> ref;
        if (bks.count(name)) ref = bks[name];
        else if (const BenchEntry* e = findEntry(name)) { if (e->bks && e->bksStatus.rfind("verified", 0) == 0) ref = e->bks; }
        const double secs = A.per100 > 0 ? inst.n * A.per100 / 100.0 : A.seconds;
        // early stop only at a PROVEN optimum of a single-depot instance; never at an unproven best-known value
        const bool isProven = ref && proven.count(name) && proven[name];
        // reference values of real-valued sets are published with two decimals: allow half a unit of the last digit
        const double target = (A.earlyStop && isProven && inst.t == 1) ? *ref + (inst.rounded ? 0.0 : 0.005) : 0.0;
        bool startRepaired = false;
        const Solution start = commonStart(inst, A.feasibleStart, &startRepaired);
        if (startRepaired) std::printf("%s: the savings start exceeded the fleet; all methods start from a repaired solution (%.2f)\n", name.c_str(), start.cost(inst));
        for (const std::string& id : ids) {
            int iters;
            auto f = budget.find(name + "," + id);
            if (f != budget.end()) iters = f->second;
            else if (!A.budgets.empty()) throw std::runtime_error("no budget for " + name + "," + id + " in " + A.budgets);
            else iters = calibrateIters(id, inst, secs);
            if (A.calibrateOnly) {
                std::snprintf(buf, sizeof buf, "%s,%s,%g,%d\n", name.c_str(), id.c_str(), secs, iters);
                out << buf << std::flush;
                std::printf("%s %s %d iterations for %g s\n", name.c_str(), id.c_str(), iters, secs);
                continue;
            }
            for (int s = A.seed0; s < A.seed0 + A.seeds; s++) {
                if (done.count(name + "," + id + "," + std::to_string(s) + "," + excludedNames)) continue;
                const double t0 = now();
                Result r = runAny(id, inst, iters, (unsigned long long)s, 0, target, &start);
                const double dt = now() - t0;
                const bool feas = r.best.violations(inst).empty();
                if (!feas) bad++;
                const double c = r.best.cost(inst);
                const HistoryEntry& h = r.history.back();
                std::string refs, gaps;
                if (ref) { std::snprintf(buf, sizeof buf, "%.6g", *ref); refs = buf; std::snprintf(buf, sizeof buf, "%.4f", 100.0 * (c / *ref - 1.0)); gaps = buf; }
                std::snprintf(buf, sizeof buf, "%s,%d,%d,%s,%d,%g,%d,%.4f,%s,%s,%d,%.3f,%.3f,%d,%d", name.c_str(), inst.n, inst.t, id.c_str(), s, secs, iters, c,
                              refs.c_str(), gaps.c_str(), r.best.nRoutes(), dt, h.seconds, h.iter, feas ? 1 : 0);
                std::string row = buf;                               // written last: the row marks the run as finished
                // convergence profile: best cost after 1, 2, 5, ... 100 % of the run time of this run (the fractions of Vidal 2022)
                for (double f : {0.01, 0.02, 0.05, 0.10, 0.15, 0.20, 0.30, 0.50, 0.75, 1.0}) {
                    double cf = r.history.front().cost;
                    // a run that stopped early at the proven optimum is measured against its planned budget
                    const double span = r.reachedTarget ? std::max(secs, r.seconds) : r.seconds;
                    for (const HistoryEntry& e : r.history) { if (e.seconds <= f * span + 1e-12) cf = e.cost; else break; }
                    std::snprintf(buf, sizeof buf, ",%.4f", cf);
                    row += buf;
                }
                row += std::string(",") + (isProven ? "1" : "0") + "," + (r.reachedTarget ? "1" : "0") + "," + (startRepaired ? "1" : "0") + "," + excludedNames + "\n";
                const std::string runKey = name + "," + id + "," + std::to_string(s) + "," + excludedNames;
                if (routesOut.is_open()) {                           // one row per route, customers numbered from 1 as in the instance file order
                    int k = 0;
                    for (size_t q = 0; q < r.best.routes.size(); q++) {
                        if (r.best.routes[q].empty()) continue;
                        double load = 0;
                        for (int v : r.best.routes[q]) load += inst.demand[v];
                        std::snprintf(buf, sizeof buf, ",%.4f,%d,%d,%g,", c, ++k, r.best.depots[q] + 1, load);
                        routesOut << runKey << buf;
                        for (size_t x = 0; x < r.best.routes[q].size(); x++) routesOut << (x ? " " : "") << r.best.routes[q][x] + 1;
                        routesOut << "\n";
                    }
                    routesOut << std::flush;
                }
                if (opsOut.is_open()) { writeOps(opsOut, runKey, r); opsOut << std::flush; }
                if (trace.is_open()) {
                    for (const HistoryEntry& e : r.history) {
                        std::snprintf(buf, sizeof buf, "%s,%s,%d,%d,%.6f,%.4f\n", name.c_str(), id.c_str(), s, e.iter, e.seconds, e.cost);
                        trace << buf;
                    }
                    trace << std::flush;
                }
                out << row << std::flush;
                std::printf("%s %s seed %d: %.2f in %.2f s%s\n", name.c_str(), id.c_str(), s, c, dt, feas ? "" : "  INFEASIBLE");
                std::fflush(stdout);
            }
        }
    }
    std::printf("wrote %s\n", path.c_str());
    return bad ? 1 : 0;
}


// ---------------------------------------------------------------------------------------------------
// run: ONE run of one method on one instance, for external experiment managers (one task = one process).
//   vrp_experiments run --method ID --instance FILE|NAME --seed K --time SEC [--iters N] [--round auto|yes|no]
//                       [--stop-at VALUE [--tol EPS]] [--no-feasible-start] --out PREFIX
// Writes PREFIX.sol (CVRPLIB style "Route #k: ids" + "Cost v" for one depot, "Route #k (depot d): ids" for several),
// PREFIX.trace.csv (iter,seconds,cost: every improvement, time-stamped from the start of the run) and PREFIX.json (summary).
// --time is converted into an iteration budget by two short timing runs unless --iters gives the budget (use --iters with the
// value of a budgets file to make every seed and node do the same work).  --stop-at: stop as soon as a feasible solution costs
// at most VALUE + EPS (pass proven optima only).  Exit status 0: feasible solution written; 1: infeasible; 2: usage error.
struct RunArgs {
    std::string method, instance, out, round = "auto";
    unsigned long long seed = 0;
    double seconds = 0, stopAt = 0, tol = 0, per100 = 0;
    int iters = -1;
    bool feasibleStart = true;
};

int expRun(const RunArgs& A) {
    if (A.method.empty() || A.instance.empty() || A.out.empty()) throw std::runtime_error("run needs --method, --instance and --out");
    if (std::find(allMethods().begin(), allMethods().end(), A.method) == allMethods().end()) throw std::runtime_error("unknown method " + A.method);
    if (A.iters < 1 && A.seconds <= 0 && A.per100 <= 0) throw std::runtime_error("run needs --time SEC, --time-per-100 SEC or --iters N");
    Instance inst;
    if (!fileExists(A.instance)) inst = loadBenchmark(A.instance);
    else if (A.round == "auto") inst = readInstance(A.instance);
    else if (A.round == "yes" || A.round == "1" || A.round == "no" || A.round == "0") {
        std::ifstream probe(A.instance);
        std::string text((std::istreambuf_iterator<char>(probe)), std::istreambuf_iterator<char>());
        if (text.find("NODE_COORD_SECTION") == std::string::npos) inst = readInstance(A.instance);      // Cordeau files have no rounding option
        else inst = readVrplib(A.instance, A.round == "yes" || A.round == "1");
    } else throw std::runtime_error("--round must be auto, yes or no");
    bool repaired = false;
    const Solution start = commonStart(inst, A.feasibleStart, &repaired);
    const double budget = A.seconds > 0 ? A.seconds : inst.n * A.per100 / 100.0;
    const int iters = A.iters >= 1 ? A.iters : calibrateIters(A.method, inst, budget);
    const double target = A.stopAt > 0 ? A.stopAt + A.tol : 0.0;
    Result r = runAny(A.method, inst, iters, A.seed, 0, target, &start);
    const bool feas = r.best.violations(inst).empty();
    const double cost = r.best.cost(inst);
    char buf[512];
    {
        std::ofstream f(A.out + ".sol");
        if (!f) throw std::runtime_error("cannot write " + A.out + ".sol");
        if (inst.t == 1) {
            int k = 0;
            for (size_t i = 0; i < r.best.routes.size(); i++) {
                if (r.best.routes[i].empty()) continue;
                f << "Route #" << ++k << ":";
                for (int c : r.best.routes[i]) f << ' ' << c + 1;
                f << '\n';
            }
            std::snprintf(buf, sizeof buf, inst.rounded ? "Cost %.0f\n" : "Cost %.6f\n", cost);
            f << buf;
        } else f << r.best.toText(inst);
    }
    {
        std::ofstream g(A.out + ".ops.csv");
        g << "operator,moves,gain\n";
        char ob[160];
        for (int k = 0; k < OP_COUNT; k++)
            if (r.opStats[k]) { std::snprintf(ob, sizeof ob, "%s,%ld,%.4f\n", opName(k), r.opStats[k], r.opGain[k]); g << ob; }
        std::ofstream f(A.out + ".trace.csv");
        f << "iter,seconds,cost\n";
        for (const HistoryEntry& e : r.history) { std::snprintf(buf, sizeof buf, "%d,%.6f,%.6f\n", e.iter, e.seconds, e.cost); f << buf; }
    }
    {
        std::ofstream f(A.out + ".json");
        std::snprintf(buf, sizeof buf,
                      "{\"instance\":\"%s\",\"n\":%d,\"depots\":%d,\"method\":\"%s\",\"seed\":%llu,\"budget_s\":%g,\"iters_budget\":%d,\"iters\":%d,"
                      "\"cost\":%.6f,\"routes\":%d,\"seconds\":%.3f,\"seconds_to_best\":%.3f,\"feasible\":%s,\"stop_at\":%g,\"reached_target\":%s,"
                      "\"start_repaired\":%s,\"rounded\":%s}\n",
                      inst.name.c_str(), inst.n, inst.t, A.method.c_str(), A.seed, A.seconds, iters, r.iters, cost, r.best.nRoutes(), r.seconds,
                      r.history.back().seconds, feas ? "true" : "false", A.stopAt, r.reachedTarget ? "true" : "false", repaired ? "true" : "false",
                      inst.rounded ? "true" : "false");
        f << buf;
    }
    std::printf("%s %s seed %llu: cost %.4f, %d routes, %.2f s, %s%s\n", inst.name.c_str(), A.method.c_str(), A.seed, cost, r.best.nRoutes(), r.seconds,
                feas ? "FEASIBLE" : "INFEASIBLE", r.reachedTarget ? ", stopped at the target" : "");
    return feas ? 0 : 1;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc < 2) { std::puts("usage: vrp_experiments construct|ops|depots|meta|methods [--iters N] [--time SEC] [--seeds S] [--out DIR] [--only a,b]\n"
                                "       vrp_experiments bench --instances NAMES|PATHS|@FILE [--methods all|a,b] [--time SEC | --time-per-100 SEC] [--seeds S] [--seed0 K]\n"
                                "                             [--budgets FILE | --calibrate-only] [--bks FILE] [--tag NAME] [--out DIR] [--no-traces]\n"
                                "                             [--no-early-stop] [--no-feasible-start] [--no-routes] [--no-ops] [--without OP[,OP]] [--resume]\n"
                                "       vrp_experiments run --method ID --instance FILE --seed K --time SEC [--iters N] [--round auto|yes|no]\n"
                                "                           [--stop-at VALUE [--tol EPS]] --out PREFIX   (writes PREFIX.sol, .trace.csv, .json)"); return 2; }
        std::string which = argv[1];
        int iters = 8000, seeds = 3;
        double seconds = 2;
        BenchArgs B;
        RunArgs R;
        bool timeSet = false, seedsSet = false, outSet = false;
        for (int a = 2; a < argc; a++) {
            std::string s = argv[a];
            auto next = [&]() { if (a + 1 >= argc) throw std::runtime_error("missing value after " + s); return std::string(argv[++a]); };
            if (s == "--iters") { iters = std::stoi(next()); R.iters = iters; }
            else if (s == "--method") R.method = next();
            else if (s == "--instance") R.instance = next();
            else if (s == "--seed") R.seed = std::stoull(next());
            else if (s == "--round") R.round = next();
            else if (s == "--stop-at") R.stopAt = std::stod(next());
            else if (s == "--tol") R.tol = std::stod(next());
            else if (s == "--seeds") { seeds = std::stoi(next()); seedsSet = true; }
            else if (s == "--time") { seconds = std::stod(next()); timeSet = true; }
            else if (s == "--instances") B.instances = next();
            else if (s == "--methods") B.methods = next();
            else if (s == "--budgets") B.budgets = next();
            else if (s == "--bks") B.bks = next();
            else if (s == "--tag") B.tag = next();
            else if (s == "--seed0") B.seed0 = std::stoi(next());
            else if (s == "--time-per-100") B.per100 = std::stod(next());
            else if (s == "--calibrate-only") B.calibrateOnly = true;
            else if (s == "--no-traces") B.traces = false;
            else if (s == "--no-routes") B.routes = false;
            else if (s == "--resume") B.resume = true;
            else if (s == "--param") {                              // name=value, may be given several times
                const std::string kv = next();
                const size_t eq = kv.find('=');
                if (eq == std::string::npos) throw std::runtime_error("--param needs name=value");
                tuned[kv.substr(0, eq)] = std::stod(kv.substr(eq + 1));
            }
            else if (s == "--no-ops") B.ops = false;
            else if (s == "--without") { excludedNames = next(); excludedOps = opsFromNames(splitList(excludedNames)); std::replace(excludedNames.begin(), excludedNames.end(), ',', '+'); }
            else if (s == "--no-early-stop") B.earlyStop = false;
            else if (s == "--no-feasible-start") B.feasibleStart = false;
            else if (s == "--out") { outDir = next(); outSet = true; }
            else if (s == "--only") { std::stringstream ss(next()); std::string t; while (std::getline(ss, t, ',')) only.insert(t); }
            else throw std::runtime_error("unknown option " + s);
        }
        if (which == "construct") expConstruct();
        else if (which == "ops") expOps();
        else if (which == "depots") expDepots(iters, seeds);
        else if (which == "meta") expMeta(iters, seeds);
        else if (which == "methods") expMethods(seconds, seeds);
        else if (which == "params") {                               // tunable parameters; with --out DIR also irace parameter files
            const std::vector<std::string> ids = B.methods == "all" ? splitList("@scripts/hpc/methods_main.txt") : splitList(B.methods);
            if (outSet) mkdir(outDir.c_str(), 0755);
            for (const std::string& id : ids) {
                const std::vector<ParamSpec>& ps = methodParams(id);
                std::printf("%s: %zu parameter(s)\n", id.c_str(), ps.size());
                std::ofstream f, d;
                if (outSet && !ps.empty()) { f.open(outDir + "/" + id + ".txt"); d.open(outDir + "/" + id + "_default.txt"); }
                std::string names, vals;
                for (const ParamSpec& p : ps) {
                    std::printf("  %-12s %c  default %-8g range [%g, %g]%s  %s\n", p.name.c_str(), p.type, p.def, p.lo, p.hi, p.logScale ? " log" : "", p.doc.c_str());
                    char line[256];
                    std::snprintf(line, sizeof line, "%-12s \"--%s \"  %c%s  (%g, %g)\n", p.name.c_str(), p.name.c_str(), p.type, p.logScale && p.lo > 0 ? ",log" : "", p.lo, p.hi);
                    if (f.is_open()) f << line;
                    char v[32];
                    std::snprintf(v, sizeof v, "%g", p.def);
                    names += (names.empty() ? "" : " ") + p.name;
                    vals += (vals.empty() ? "" : " ") + std::string(v);
                }
                if (d.is_open()) d << names << "\n" << vals << "\n";
            }
        }
        else if (which == "run") { R.per100 = B.per100; if (timeSet) R.seconds = seconds; R.out = outSet ? outDir : std::string(); R.feasibleStart = B.feasibleStart; return expRun(R); }
        else if (which == "bench") { if (timeSet) B.seconds = seconds; if (seedsSet) B.seeds = seeds; return expBench(B); }
        else throw std::runtime_error("unknown experiment " + which);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 2;
    }
    return 0;
}
