#include "vrplearn/construct.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>

#include "vrplearn/pyrandom.hpp"

namespace vrp {

std::vector<int> assignNearestDepot(const Instance& inst) {
    std::vector<int> out(inst.n);
    for (int i = 0; i < inst.n; i++) {
        int best = 0;
        for (int d = 1; d < inst.t; d++)
            if (inst.d(i, inst.n + d) < inst.d(i, inst.n + best)) best = d;
        out[i] = best;
    }
    return out;
}

namespace {

constexpr double TOL = 1e-9;

// ---- Clarke & Wright for one depot and a customer subset ---------------------------------------
std::vector<std::vector<int>> cwSingle(const Instance& inst, const std::vector<int>& customers, int dep, double lam) {
    const int dn = inst.depotNode(dep);
    const double Q = inst.capacity[dep], dmax = inst.maxDuration[dep];
    const int n = inst.n;
    // route id = id of the customer that started it; insertion order of the ids = `customers` order
    std::vector<std::vector<int>> routes(n);
    std::vector<int> where(n, -1);
    std::vector<double> load(n, 0.0), dur(n, 0.0);
    std::vector<char> alive(n, 0);
    for (int c : customers) {
        routes[c] = {c};
        where[c] = c;
        load[c] = inst.demand[c];
        dur[c] = 2 * inst.d(dn, c) + inst.service[c];
        alive[c] = 1;
    }
    struct Sav { double s; int i, j; };
    std::vector<Sav> sav;
    for (size_t a = 0; a < customers.size(); a++)
        for (size_t b = a + 1; b < customers.size(); b++) {
            int i = customers[a], j = customers[b];
            sav.push_back({inst.d(dn, i) + inst.d(dn, j) - lam * inst.d(i, j), i, j});
        }
    std::sort(sav.begin(), sav.end(), [](const Sav& p, const Sav& q) {          // descending (s, i, j)
        if (p.s != q.s) return p.s > q.s;
        if (p.i != q.i) return p.i > q.i;
        return p.j > q.j;
    });
    for (const Sav& e : sav) {
        if (e.s <= 0) break;
        int i = e.i, j = e.j;
        int ri = where[i], rj = where[j];
        if (ri == rj) continue;
        std::vector<int>& a = routes[ri];
        std::vector<int>& b = routes[rj];
        if (a.front() != i && a.back() != i) continue;                          // i, j must be route endpoints
        if (b.front() != j && b.back() != j) continue;
        if (load[ri] + load[rj] > Q + TOL) continue;
        if (a.back() != i) std::reverse(a.begin(), a.end());
        if (b.front() != j) std::reverse(b.begin(), b.end());
        double newDur = dur[ri] + dur[rj] - inst.d(dn, i) - inst.d(dn, j) + inst.d(i, j);
        if (newDur > dmax + TOL) continue;
        a.insert(a.end(), b.begin(), b.end());
        load[ri] += load[rj];
        dur[ri] = newDur;
        for (int v : b) where[v] = ri;
        alive[rj] = 0;
        b.clear();
    }
    std::vector<std::vector<int>> out;
    for (int c : customers)
        if (alive[c]) out.push_back(routes[c]);
    return out;
}

double routeLoad(const Instance& inst, const std::vector<int>& nodes) {
    double s = 0;
    for (int v : nodes) s += inst.demand[v];
    return s;
}

double routeDuration(const Instance& inst, int dn, const std::vector<int>& nodes) {
    double s = 0;
    for (int v : nodes) s += inst.service[v];
    return tourCost(inst, dn, nodes) + s;
}

// ---- fleet repair -----------------------------------------------------------------------------
// Move every customer of route r0 into other routes (cheapest feasible insertion); all-or-nothing.
bool tryDissolve(const Instance& inst, Solution& sol, size_t r0) {
    std::vector<std::vector<int>> routes = sol.routes;
    std::vector<int> todo = routes[r0];
    for (int c : todo) {
        int bestR = -1, bestP = -1;
        double bd = INF;
        for (size_t r = 0; r < routes.size(); r++) {
            if (r == r0 || routes[r].empty()) continue;
            int d = sol.depots[r];
            if (routeLoad(inst, routes[r]) + inst.demand[c] > inst.capacity[d] + TOL) continue;
            int dn = inst.depotNode(d);
            const auto& nodes = routes[r];
            for (size_t p = 0; p <= nodes.size(); p++) {
                int a = p == 0 ? dn : nodes[p - 1], b = p == nodes.size() ? dn : nodes[p];
                double delta = inst.d(a, c) + inst.d(c, b) - inst.d(a, b);
                if (delta < bd) {
                    std::vector<int> trial(nodes);
                    trial.insert(trial.begin() + p, c);
                    if (routeDuration(inst, dn, trial) <= inst.maxDuration[d] + TOL) { bestR = (int)r; bestP = (int)p; bd = delta; }
                }
            }
        }
        if (bestR < 0) return false;
        routes[r0].erase(std::find(routes[r0].begin(), routes[r0].end(), c));
        routes[bestR].insert(routes[bestR].begin() + bestP, c);
    }
    sol.routes = routes;
    return true;
}

Solution repairFleet(const Instance& inst, Solution sol) {
    for (int it = 0; it < 4 * inst.n; it++) {
        std::vector<int> cnt = sol.routesPerDepot(inst);
        int d0 = -1;
        for (int d = 0; d < inst.t; d++)
            if (cnt[d] > inst.fleet[d]) { d0 = d; break; }
        if (d0 < 0) break;
        std::vector<size_t> cand;
        for (size_t r = 0; r < sol.routes.size(); r++)
            if (sol.depots[r] == d0 && !sol.routes[r].empty()) cand.push_back(r);
        std::stable_sort(cand.begin(), cand.end(), [&](size_t a, size_t b) {
            return routeLoad(inst, sol.routes[a]) < routeLoad(inst, sol.routes[b]);
        });
        bool ok = false;
        for (size_t r : cand)
            if (tryDissolve(inst, sol, r)) { ok = true; break; }
        if (!ok) break;
        sol = sol.nonEmpty();
    }
    return sol.nonEmpty();
}

// ---- sweep ------------------------------------------------------------------------------------
std::vector<std::vector<int>> sweepOne(const Instance& inst, const std::vector<int>& customers, int dep, double start) {
    const int dn = inst.depotNode(dep);
    const double cx = inst.x[dn], cy = inst.y[dn];
    const double twoPi = 2 * M_PI;
    std::vector<std::pair<double, int>> keyed;
    for (int c : customers) keyed.push_back({pymod(std::atan2(inst.y[c] - cy, inst.x[c] - cx) - start, twoPi), c});
    std::stable_sort(keyed.begin(), keyed.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
    std::vector<std::vector<int>> routes;
    std::vector<int> cur;
    double load = 0;
    const double Q = inst.capacity[dep];
    for (auto& kc : keyed) {
        int c = kc.second;
        double extra = inst.demand[c];
        std::vector<int> trial = cur;
        trial.push_back(c);
        double pairs = 0, svc = 0;
        for (size_t k = 0; k + 1 < trial.size(); k++) pairs += inst.d(trial[k], trial[k + 1]);
        for (int v : trial) svc += inst.service[v];
        double dur = inst.d(dn, trial.front()) + inst.d(trial.back(), dn) + pairs + svc;
        if (!cur.empty() && (load + extra > Q + TOL || dur > inst.maxDuration[dep] + TOL)) {
            routes.push_back(cur);
            cur.clear();
            load = 0;
        }
        cur.push_back(c);
        load += extra;
    }
    if (!cur.empty()) routes.push_back(cur);
    return routes;
}

}  // namespace

Solution savingsMd(const Instance& inst, double lam, const std::vector<int>* assignment) {
    std::vector<int> asg = assignment ? *assignment : assignNearestDepot(inst);
    Solution sol;
    for (int d = 0; d < inst.t; d++) {
        std::vector<int> cs;
        for (int i = 0; i < inst.n; i++) if (asg[i] == d) cs.push_back(i);
        for (auto& r : cwSingle(inst, cs, d, lam)) { sol.routes.push_back(r); sol.depots.push_back(d); }
    }
    return repairFleet(inst, sol);
}

Solution savings(const Instance& inst, double lam) { return savingsMd(inst, lam); }

Solution sweep(const Instance& inst, int nStarts) {
    std::vector<int> asg = assignNearestDepot(inst);
    Solution best;
    double bestCost = INF;
    bool have = false;
    for (int s = 0; s < nStarts; s++) {
        double start = 2 * M_PI * s / nStarts;
        Solution sol;
        for (int d = 0; d < inst.t; d++) {
            std::vector<int> cs;
            for (int i = 0; i < inst.n; i++) if (asg[i] == d) cs.push_back(i);
            for (auto& r : sweepOne(inst, cs, d, start)) { sol.routes.push_back(r); sol.depots.push_back(d); }
        }
        double c = sol.cost(inst);
        if (!have || c < bestCost) { best = sol; bestCost = c; have = true; }
    }
    return repairFleet(inst, best);
}

Solution nearestNeighbor(const Instance& inst) {
    std::vector<int> asg = assignNearestDepot(inst);
    Solution sol;
    for (int d = 0; d < inst.t; d++) {
        PySetInt todo;                                  // same iteration order as the reference implementation
        for (int i = 0; i < inst.n; i++) if (asg[i] == d) todo.add(i);
        const int dn = inst.depotNode(d);
        while (!todo.empty()) {
            std::vector<int> cur;
            double load = 0, dur = 0;
            int pos = dn;
            while (true) {
                int best = -1;
                double bd = INF;
                for (long long cc : todo.items()) {
                    int c = static_cast<int>(cc);
                    if (load + inst.demand[c] > inst.capacity[d] + TOL) continue;
                    double nd = dur + inst.d(pos, c) + inst.service[c] + inst.d(c, dn);
                    if (nd > inst.maxDuration[d] + TOL) continue;
                    if (inst.d(pos, c) < bd) { best = c; bd = inst.d(pos, c); }
                }
                if (best < 0) break;
                cur.push_back(best);
                todo.discard(best);
                load += inst.demand[best];
                dur += inst.d(pos, best) + inst.service[best];
                pos = best;
            }
            if (cur.empty()) cur = {static_cast<int>(todo.pop())};     // a single customer violates D/Q: route alone
            sol.routes.push_back(cur);
            sol.depots.push_back(d);
        }
    }
    return repairFleet(inst, sol);
}

Solution cheapestInsertion(const Instance& inst, std::optional<unsigned long long> shuffleSeed) {
    const int n = inst.n, t = inst.t;
    std::vector<double> nearest(n);
    for (int i = 0; i < n; i++) {
        double m = INF;
        for (int d = 0; d < t; d++) m = std::min(m, inst.d(i, n + d));
        nearest[i] = m;
    }
    std::vector<int> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(), [&](int a, int b) { return -nearest[a] < -nearest[b]; });
    if (shuffleSeed) { PyRandom rng(*shuffleSeed); rng.shuffle(order); }

    std::vector<std::vector<int>> routes;
    std::vector<int> depots;
    std::vector<double> load;
    std::vector<int> cnt(t, 0);
    for (int c : order) {
        int bestR = -1, bestP = -1;
        double bestDelta = INF;
        for (size_t r = 0; r < routes.size(); r++) {
            int d = depots[r];
            if (load[r] + inst.demand[c] > inst.capacity[d] + TOL) continue;
            int dn = inst.depotNode(d);
            const auto& nodes = routes[r];
            for (size_t p = 0; p <= nodes.size(); p++) {
                int a = p == 0 ? dn : nodes[p - 1], b = p == nodes.size() ? dn : nodes[p];
                double delta = inst.d(a, c) + inst.d(c, b) - inst.d(a, b);
                if (delta < bestDelta) {
                    std::vector<int> trial(nodes);
                    trial.insert(trial.begin() + p, c);
                    if (routeDuration(inst, dn, trial) <= inst.maxDuration[d] + TOL) { bestR = (int)r; bestP = (int)p; bestDelta = delta; }
                }
            }
        }
        double openCost = INF;
        int openD = -1;
        for (int d = 0; d < t; d++) {
            double oc = 2 * inst.d(c, inst.depotNode(d));
            if (oc < openCost && (cnt[d] < inst.fleet[d] || bestR < 0)) { openCost = oc; openD = d; }
        }
        if (bestR >= 0 && (bestDelta <= openCost || openD < 0)) {
            routes[bestR].insert(routes[bestR].begin() + bestP, c);
            load[bestR] += inst.demand[c];
        } else {
            if (openD < 0) {
                openD = 0;
                for (int d = 1; d < t; d++) if (inst.d(c, n + d) < inst.d(c, n + openD)) openD = d;
            }
            routes.push_back({c});
            depots.push_back(openD);
            load.push_back(inst.demand[c]);
            cnt[openD]++;
        }
    }
    Solution sol;
    sol.routes = routes;
    sol.depots = depots;
    return repairFleet(inst, sol);
}

Solution repairFleetLimit(const Instance& inst, Solution sol) { return repairFleet(inst, std::move(sol)); }

}  // namespace vrp
