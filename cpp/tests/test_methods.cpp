// Tests of the C++ twins of the browser methods (sa, tabu, vns, alns, ga) and of the optimal split.
#include <algorithm>
#include <cmath>
#include <numeric>
#include <set>
#include <sys/stat.h>
#include <fstream>
#include <string>

#include "testing.hpp"
#include "vrplearn/benchmarks.hpp"
#include "vrplearn/construct.hpp"
#include "vrplearn/exact.hpp"
#include "vrplearn/methods.hpp"

using namespace vrp;

namespace {
GenParams gen(int n, int t, unsigned long long seed, double avgRoute = 6, double maxDur = INF) {
    GenParams p;
    p.n = n; p.t = t; p.seed = seed; p.avgRouteSize = avgRoute; p.maxDuration = maxDur;
    return p;
}
// small budgets so that the whole file runs in a few seconds
int smallBudget(const std::string& id) {
    if (id == "sa" || id == "ta" || id == "rtr") return 60000;
    if (id == "tabu") return 3000;
    if (id == "alns" || id == "lns" || id == "ig") return 1500;
    if (id == "brkga" || id == "pso" || id == "de" || id == "rkrs") return 2500;
    if (id == "vns") return 40;
    return 60;                                     // one descent per iteration: grasp, msls, ga, aco, grasppr, scatter
}
}  // namespace

TEST(split_is_optimal_and_feasible) {
    // Prins split vs. enumeration of all 2^(m-1) ways to cut the giant tour (capacity and duration respected).
    for (double maxDur : {INF, 900.0}) {
        Instance inst = generate(gen(9, 2, 4, 3, maxDur));
        PyRandom rng(2);
        for (int rep = 0; rep < 12; rep++) {
            std::vector<int> tour(inst.n);
            std::iota(tour.begin(), tour.end(), 0);
            rng.shuffle(tour);
            tour.resize(5 + rep % 5);
            for (int d = 0; d < inst.t; d++) {
                auto routes = splitTour(inst, d, tour);
                std::vector<int> flat;
                Solution s;
                for (auto& r : routes) { flat.insert(flat.end(), r.begin(), r.end()); s.routes.push_back(r); s.depots.push_back(d); }
                CHECK(flat == tour);                                                  // order preserved, every customer once
                double best = INF;
                const int m = (int)tour.size();
                for (int cuts = 0; cuts < (1 << (m - 1)); cuts++) {
                    Solution c;
                    c.routes.assign(1, {});
                    c.depots.assign(1, d);
                    for (int k = 0; k < m; k++) {
                        c.routes.back().push_back(tour[k]);
                        if (k < m - 1 && (cuts >> k & 1)) { c.routes.emplace_back(); c.depots.push_back(d); }
                    }
                    bool ok = true;
                    for (size_t r = 0; r < c.routes.size() && ok; r++) {
                        double q = 0, sv = 0;
                        for (int v : c.routes[r]) { q += inst.demand[v]; sv += inst.service[v]; }
                        ok = q <= inst.capacity[d] + 1e-9 && c.routeCost(inst, r) + sv <= inst.maxDuration[d] + 1e-9;
                    }
                    if (ok) best = std::min(best, c.cost(inst));
                }
                if (best < INF) {
                    CHECK_NEAR(s.cost(inst), best, 1e-6);
                    // a partial tour: the checker reports the unvisited customers (and maybe the fleet), but never load or duration
                    const std::vector<std::string> viol = s.violations(inst);
                    CHECK(std::none_of(viol.begin(), viol.end(), [](const std::string& v) {
                        return v.find("load") != std::string::npos || v.find("duration") != std::string::npos; }));
                }
            }
        }
    }
}

TEST(catalogue_has_eighteen_methods_with_defaults) {
    CHECK(methodIds().size() == 18);
    std::set<std::string> seen(methodIds().begin(), methodIds().end());
    CHECK(seen.size() == 18);                      // no duplicates
    for (const std::string& id : methodIds()) CHECK(defaultIters(id) > 0);
    bool thrown = false;
    try { defaultIters("nope"); } catch (const std::exception&) { thrown = true; }
    CHECK(thrown);
}

TEST(methods_are_feasible_monotone_and_deterministic) {
    for (const std::string& id : methodIds()) {
        for (int t : {1, 3}) {
            Instance inst = generate(gen(45, t, 11));
            Solution init = savingsMd(inst);
            MethodOptions o;
            o.iters = smallBudget(id); o.seed = 5; o.init = &init;
            Result r1 = runMethod(id, inst, o), r2 = runMethod(id, inst, o);
            CHECK_MSG(r1.best.violations(inst).empty(), id << " t=" << t);
            CHECK_MSG(r1.cost <= init.cost(inst) + 1e-9, id << " t=" << t);
            CHECK_NEAR(r1.cost, r1.best.cost(inst), 1e-9);
            CHECK_NEAR(r1.cost, r2.cost, 1e-9);                                            // same seed, same trajectory
            CHECK(r1.iters >= 1 && r1.iters <= o.iters);
            CHECK(!r1.history.empty());
            for (size_t k = 0; k + 1 < r1.history.size(); k++) {
                CHECK_MSG(r1.history[k].cost > r1.history[k + 1].cost - 1e-12, id);        // best-so-far never worsens
                CHECK_MSG(r1.history[k].iter < r1.history[k + 1].iter, id);                // one point per iteration
            }
            CHECK_NEAR(r1.history.back().cost, r1.cost, 1e-9);
        }
    }
}

TEST(methods_depend_on_the_seed) {
    Instance inst = generate(gen(60, 2, 3));
    for (const std::string& id : methodIds()) {
        if (id == "rkrs") continue;                // the control baseline rarely improves on the savings start, whatever the seed
        bool differ = false;
        std::vector<double> costs;
        for (unsigned long long s = 1; s <= 4; s++) {
            MethodOptions o;
            o.iters = smallBudget(id) / 3 + 5; o.seed = s;
            Result r = runMethod(id, inst, o);
            costs.push_back(r.cost);
            if (std::fabs(r.cost - costs[0]) > 1e-9) differ = true;
        }
        CHECK_MSG(differ, id << " gives identical results for four seeds");
    }
}

TEST(methods_respect_the_time_limit) {
    Instance inst = generate(gen(150, 2, 1));
    for (const std::string& id : methodIds()) {
        MethodOptions o;
        o.iters = 100000000; o.timeLimit = 0.3;
        Result r = runMethod(id, inst, o);
        CHECK_MSG(r.seconds < 3.0 && r.iters < 100000000, id << " ran " << r.seconds << " s");
        CHECK(r.best.violations(inst).empty());
    }
}

TEST(methods_never_beat_the_exact_optimum) {
    // On tiny instances nothing can be below the DP optimum; ALNS and GA with a modest budget reach it.
    for (auto cfg : {std::make_tuple(8, 1, 21ULL), std::make_tuple(8, 2, 22ULL), std::make_tuple(9, 3, 23ULL)}) {
        auto [n, t, seed] = cfg;
        Instance inst = generate(gen(n, t, seed, 3));
        ExactResult ex = solveExact(inst);
        CHECK(ex.feasible);
        for (const std::string& id : methodIds()) {
            MethodOptions o;
            o.iters = smallBudget(id) * 3;
            o.seed = 9;
            Result r = runMethod(id, inst, o);
            CHECK_MSG(r.cost >= ex.cost - 1e-6, id << " n=" << n << " t=" << t << ": " << r.cost << " < " << ex.cost);
            if (id == "alns" || id == "ga") CHECK_MSG(r.cost <= ex.cost + 1e-6, id << " n=" << n << " t=" << t << ": " << r.cost << " > " << ex.cost);
        }
    }
}

TEST(target_stops_a_run_at_the_optimum_and_only_then) {
    Instance inst = generate(gen(9, 1, 31ULL, 3));
    const double opt = solveExact(inst).cost;
    // every method that reaches the optimum stops there, with the flag set and iterations to spare
    int stopped = 0;
    std::vector<std::string> ids = methodIds();
    ids.push_back("sisr");
    ids.push_back("ils");
    for (const std::string& id : ids) {
        Result r;
        const int budget = id == "sisr" ? 200000 : id == "ils" ? 2000 : defaultIters(id);
        if (id == "sisr") { SisrOptions o; o.iters = budget; o.seed = 4; o.target = opt + 1e-7; r = sisr(inst, o); }
        else if (id == "ils") { IlsOptions o; o.iters = budget; o.seed = 4; o.target = opt + 1e-7; r = ils(inst, o); }
        else { MethodOptions o; o.iters = std::min(budget, 20000); o.seed = 4; o.target = opt + 1e-7; r = runMethod(id, inst, o); }
        CHECK(r.best.violations(inst).empty());
        if (r.reachedTarget) {
            stopped++;
            CHECK_MSG(std::fabs(r.cost - opt) < 1e-6, id << " stopped at " << r.cost << " instead of " << opt);
            CHECK_MSG(r.iters < budget, id << " reached the target but used the whole budget");
        } else CHECK_MSG(r.cost > opt + 1e-7, id << " reached the optimum without setting the flag");
    }
    CHECK_MSG(stopped >= 10, "only " << stopped << " methods reached the optimum of a 9-customer instance");
    // an unreachable target changes nothing: same result as without a target
    for (const std::string& id : {std::string("sa"), std::string("ga"), std::string("brkga")}) {
        MethodOptions a, b;
        a.iters = b.iters = smallBudget(id);
        a.seed = b.seed = 2;
        b.target = opt * 0.5;
        Result ra = runMethod(id, inst, a), rb = runMethod(id, inst, b);
        CHECK(!rb.reachedTarget);
        CHECK_MSG(ra.cost == rb.cost && ra.iters == rb.iters, id << ": an unreachable target changed the run");
    }
}

TEST(vrplib_reader_reads_duration_limit_service_time_and_real_distances) {
    const std::string dir = "/tmp/vrplearn_test_reader";
    mkdir(dir.c_str(), 0755);
    auto write = [&](const std::string& name, const std::string& extra) {
        const std::string p = dir + "/" + name + ".vrp";
        std::ofstream f(p);
        f << "NAME : " << name << "\nTYPE : CVRP\nDIMENSION : 4\nEDGE_WEIGHT_TYPE : EUC_2D\nCAPACITY : 10\n" << extra
          << "NODE_COORD_SECTION\n1 0 0\n2 1 1\n3 3 0\n4 0 2.5\nDEMAND_SECTION\n1 0\n2 3\n3 4\n4 5\nDEPOT_SECTION\n1\n-1\nEOF\n";
        return p;
    };
    Instance a = readInstance(write("T-n4-k2", ""));                       // ordinary set: rounded, no duration limit
    CHECK(a.rounded && a.maxDuration[0] == INF && a.service[0] == 0.0);
    CHECK(a.d(0, a.n) == 1.0);                                             // sqrt(2) rounds to 1
    Instance b = readInstance(write("CMT99", "DISTANCE : 20.5\nSERVICE_TIME : 2\n"));
    CHECK(!b.rounded && b.maxDuration[0] == 20.5 && b.service[2] == 2.0);
    CHECK(std::fabs(b.d(0, b.n) - std::sqrt(2.0)) < 1e-12);
    // a route whose travel + service time exceeds the limit is reported by the checker
    Solution s;
    s.routes = {{0, 1, 2}};
    s.depots = {0};
    Instance c = readInstance(write("Golden_99", "DISTANCE : 9\nSERVICE_TIME : 1\n"));
    bool dur = false;
    for (const std::string& v : s.violations(c)) if (v.find("duration") != std::string::npos) dur = true;
    CHECK(dur);
    CHECK(readVrplib(dir + "/CMT99.vrp", true).rounded);                   // the caller can override the convention
}

TEST(methods_never_below_published_best_known) {
    for (const char* name : {"p01", "X-n101-k25"}) {
        if (!available(name)) SKIP("run scripts/fetch_instances.sh");
        Instance inst = loadBenchmark(name);
        Solution init = savingsMd(inst);
        for (const std::string& id : methodIds()) {
            MethodOptions o;
            o.iters = smallBudget(id); o.seed = 1; o.init = &init;
            Result r = runMethod(id, inst, o);
            CHECK_MSG(r.best.violations(inst).empty(), id << " " << name);
            CHECK_MSG(r.cost >= *findEntry(name)->bks - 0.006, id << " " << name << ": " << r.cost);
            CHECK(r.cost <= init.cost(inst) + 1e-9);
        }
    }
}

// ---- upgraded ("full") variants: the same properties as the compact methods ------------------------------------------------
namespace {
int fullBudget(const std::string& id) {
    if (id == "tabu-f") return 400;
    if (id == "alns-f" || id == "lns-f") return 600;
    if (id == "sa-f") return 1500000;                  // the hyperbolic cooling of Osman needs a long run to pay off
    if (id == "ta-f" || id == "rtr-f") return 60000;
    if (id == "ig-f") return 150;
    return 60;
}
}  // namespace

TEST(full_methods_are_listed_and_distinct_from_the_compact_ones) {
    CHECK(!fullMethodIds().empty());
    for (const std::string& id : fullMethodIds()) {
        CHECK(std::find(methodIds().begin(), methodIds().end(), id) == methodIds().end());
        CHECK(defaultIters(id) > 0);
    }
}

TEST(full_methods_are_feasible_monotone_and_deterministic) {
    for (const std::string& id : fullMethodIds())
        for (int t : {1, 3})
            for (double maxDur : {INF, 1500.0}) {
                Instance inst = generate(gen(45, t, 11, 6, maxDur));
                Solution init = savingsMd(inst);
                MethodOptions o;
                o.iters = fullBudget(id); o.seed = 5; o.init = &init;
                Result r1 = runMethod(id, inst, o), r2 = runMethod(id, inst, o);
                CHECK_MSG(r1.best.violations(inst).empty(), id << " t=" << t << " D=" << maxDur);
                CHECK_MSG(r1.cost <= init.cost(inst) + 1e-9, id << " t=" << t);
                CHECK_NEAR(r1.cost, r1.best.cost(inst), 1e-9);
                CHECK_NEAR(r1.cost, r2.cost, 1e-9);
                CHECK(r1.best.routes == r2.best.routes);
                CHECK(r1.iters >= 1 && r1.iters <= o.iters);
                for (size_t k = 0; k + 1 < r1.history.size(); k++) {
                    CHECK_MSG(r1.history[k].cost > r1.history[k + 1].cost - 1e-12, id);
                    CHECK_MSG(r1.history[k].iter < r1.history[k + 1].iter, id);
                }
                CHECK_NEAR(r1.history.back().cost, r1.cost, 1e-9);
            }
}

TEST(full_methods_improve_depend_on_the_seed_and_respect_the_time_limit) {
    Instance inst = generate(gen(80, 2, 3));
    Solution init = savingsMd(inst);
    for (const std::string& id : fullMethodIds()) {
        std::set<long long> costs;
        for (unsigned long long s = 1; s <= 4; s++) {
            MethodOptions o;
            o.iters = fullBudget(id); o.seed = s; o.init = &init;
            Result r = runMethod(id, inst, o);
            CHECK_MSG(r.cost < init.cost(inst) - 1e-6, id << " did not improve the savings start");
            costs.insert(std::llround(r.cost * 1e6));
        }
        CHECK_MSG(costs.size() > 1, id << " gives identical results for four seeds");
        Instance big = generate(gen(200, 2, 1));
        MethodOptions o;
        o.iters = 100000000; o.timeLimit = 0.3;
        Result r = runMethod(id, big, o);
        CHECK_MSG(r.seconds < 3.0 && r.iters < 100000000, id << " ran " << r.seconds << " s");
        CHECK(r.best.violations(big).empty());
    }
}

TEST(full_methods_reach_but_never_beat_the_exact_optimum) {
    for (auto cfg : {std::make_tuple(8, 1, 21ULL), std::make_tuple(8, 2, 22ULL), std::make_tuple(9, 3, 23ULL)}) {
        auto [n, t, seed] = cfg;
        Instance inst = generate(gen(n, t, seed, 3));
        ExactResult ex = solveExact(inst);
        CHECK(ex.feasible);
        for (const std::string& id : fullMethodIds()) {
            MethodOptions o;
            o.iters = fullBudget(id) * 3; o.seed = 9;
            Result r = runMethod(id, inst, o);
            CHECK_MSG(r.cost >= ex.cost - 1e-6, id << " n=" << n << " t=" << t << ": " << r.cost << " < " << ex.cost);
            CHECK_MSG(r.cost <= ex.cost * 1.02 + 1e-6, id << " n=" << n << " t=" << t << ": " << r.cost << " vs optimum " << ex.cost);
        }
    }
}

TEST(operator_statistics_and_ablation_switch) {
    Instance inst = generate(gen(60, 2, 4));
    for (const char* id : {"ga", "vns", "ga-f", "ils-f", "ig-f"}) {
        MethodOptions o;
        o.iters = 40; o.seed = 2;
        Result all = runMethod(id, inst, o);
        long moves = 0;
        for (long m : all.opStats) moves += m;
        CHECK_MSG(moves > 0, id << " reports no local-search moves");
        o.excludeOps = opBit(OP_RELOCATE) | opBit(OP_TWO_OPT_STAR);
        Result part = runMethod(id, inst, o);
        CHECK_MSG(part.opStats[OP_RELOCATE] == 0 && part.opStats[OP_TWO_OPT_STAR] == 0 && part.opStats[OP_NEW_ROUTE] == 0, id << " used a switched-off operator");
        CHECK(part.best.violations(inst).empty());
    }
    MethodOptions o;                                  // methods without descents report nothing
    o.iters = 2000;
    Result r = runMethod("sa", inst, o);
    for (long m : r.opStats) CHECK(m == 0);
}

TEST(tunable_parameters_have_valid_tables_and_take_effect) {
    Instance inst = generate(gen(60, 2, 4));
    for (const char* id : {"sisr", "ils-f", "sa-f", "ta", "rtr", "tabu-f", "vns-f", "grasp", "msls", "alns-f", "lns", "ig-f", "ga-f", "aco", "grasppr-f",
                           "scatter-f", "brkga-f", "pso", "de-f", "rkrs"}) {
        const std::string sid = id;
        const std::vector<ParamSpec>& specs = methodParams(sid);
        std::set<std::string> names;
        for (const ParamSpec& p : specs) {
            CHECK_MSG(names.insert(p.name).second, id << ": duplicate parameter " << p.name);
            CHECK_MSG(p.lo <= p.def && p.def <= p.hi && (p.type == 'r' || p.type == 'i'), id << "." << p.name);
            CHECK_MSG(!p.logScale || p.lo > 0, id << "." << p.name << ": log scale needs a positive lower bound");
            if (p.type == 'i') CHECK(p.def == std::floor(p.def) && p.lo == std::floor(p.lo) && p.hi == std::floor(p.hi));
        }
        if (std::string(id) == "sisr") continue;                    // its parameters are fields of SisrOptions
        // both ends of every range give a feasible run; an explicit default reproduces the run without parameters
        MethodOptions o;
        o.iters = fullBudget(id) / 4 + 20; o.seed = 3;
        const Result base = runMethod(id, inst, o);
        MethodOptions d = o;
        for (const ParamSpec& p : specs) d.params[p.name] = p.def;
        CHECK_MSG(std::fabs(runMethod(id, inst, d).cost - base.cost) < 1e-9, id << ": explicit defaults change the result");
        for (int end = 0; end < 2; end++) {
            MethodOptions e = o;
            for (const ParamSpec& p : specs) e.params[p.name] = end ? p.hi : p.lo;
            CHECK_MSG(runMethod(id, inst, e).best.violations(inst).empty(), id << " at the " << (end ? "upper" : "lower") << " ends of its ranges");
        }
    }
    MethodOptions bad;
    bad.iters = 10;
    bad.params["nope"] = 1;
    bool thrown = false;
    try { runMethod("ga-f", inst, bad); } catch (const std::exception&) { thrown = true; }
    CHECK(thrown);
}
