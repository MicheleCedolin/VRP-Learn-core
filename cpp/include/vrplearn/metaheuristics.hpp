// Ruin-and-recreate machinery and two metaheuristics built on it.
//   sisr : Slack-Induction-by-String-Removals-style ruin & recreate + simulated annealing
//          (Christiaens & Vanden Berghe 2020), generalised to several depots: customers may be
//          re-inserted into any route of any depot, new routes open at the cheapest depot.
//   ils  : iterated local search: ruin & recreate perturbation + the granular local search of
//          localsearch.hpp (relocate / swap / 2-opt* / SWAP* / depot moves), SA-style acceptance.
// Teaching-grade reference implementations: not tuned to reproduce published results.
#pragma once
#include <string>
#include <tuple>
#include <vector>

#include "vrplearn/localsearch.hpp"

namespace vrp {

// Mutable solution used by ruin & recreate (lists of lists + cached route data).
class RRState {
public:
    RRState(const LSContext& ctx, const Solution& sol);
    int fleetExcess() const;
    Solution toSolution() const;
    // Remove adjacent strings (Christiaens & Vanden Berghe); returns the removed customers.
    std::vector<int> ruin(PyRandom& rng, double cbar = 10.0, int lmax = 10, double alpha = 0.01);
    // Greedy insertion with random blinks; `removed` is reordered in place.
    void recreate(std::vector<int>& removed, PyRandom& rng, double blink = 0.01);

    std::vector<std::vector<int>> routes;
    std::vector<int> depots;
    std::vector<double> load, costR, dur;
    std::vector<int> rid;                       // route index of every customer
    std::vector<int> cnt;                       // routes per depot
    double cost = 0;

private:
    const LSContext* ctx_;
    int n_, t_;
    void recompute(size_t r);
};

struct HistoryEntry { int iter; double seconds; double cost; };

struct NamedStat { std::string name; long count; double value; };

struct Result {
    Solution best;
    double cost = 0;
    int iters = 0;
    double seconds = 0;
    std::vector<HistoryEntry> history;
    bool reachedTarget = false;                 // the run stopped early because a feasible solution reached Options::target
    std::array<long, OP_COUNT> opStats{};       // ILS only: moves applied per local-search operator
    std::array<double, OP_COUNT> opGain{};
    std::vector<NamedStat> extra;               // method-specific counters (ALNS: selections and final weight of every operator)
};

struct SisrOptions {
    int iters = 20000;
    unsigned long long seed = 0;
    double timeLimit = 0;                       // seconds, 0 = none
    double target = 0;                          // > 0: stop as soon as a feasible solution costs <= target (use proven optima only)
    double cbar = 10.0;
    int lmax = 10;
    double alpha = 0.01;
    double blink = 0.01;
    double t0Factor = 0.45;                     // temperatures are relative to the mean edge length
    double tfFactor = 0.0045;
    double lamFleet = 1e6;
    const Solution* init = nullptr;             // default: savings
    LSContext* ctx = nullptr;                   // default: built internally (k = 20)
};
Result sisr(const Instance& inst, const SisrOptions& opt = SisrOptions());

struct IlsOptions {
    int iters = 500;
    unsigned long long seed = 0;
    double timeLimit = 0;
    double target = 0;                          // as in SisrOptions
    std::string preset = "hgs";
    const std::vector<std::string>* ops = nullptr;
    int k = 20;
    double cbar = 6.0;
    int lmax = 6;
    double blink = 0.01;
    double t0Factor = 0.05;
    double tfFactor = 0.001;
    double lamFleet = 1e6;
    const Solution* init = nullptr;
    LSContext* ctx = nullptr;
};
Result ils(const Instance& inst, const IlsOptions& opt = IlsOptions());

}  // namespace vrp
