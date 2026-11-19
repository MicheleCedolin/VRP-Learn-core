#include "vrplearn/methods.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <functional>
#include <limits>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <tuple>
#include <unordered_map>

#include "vrplearn/construct.hpp"

namespace vrp {

namespace {
constexpr double EPS = 1e-9;

double now() {
    using namespace std::chrono;
    return duration<double>(steady_clock::now().time_since_epoch()).count();
}

int excessOf(const Instance& I, const std::vector<int>& cnt) {
    int e = 0;
    for (int d = 0; d < I.t; d++) e += std::max(0, cnt[d] - I.fleet[d]);
    return e;
}
double objective(const Instance& I, const Solution& s, double lam) { return s.cost(I) + lam * excessOf(I, s.routesPerDepot(I)); }

// Best-so-far bookkeeping shared by all methods: one history point per iteration, as in the JavaScript engine.
struct Tracker {
    const Instance& I;
    double lam, t0, bestO;
    Solution best;
    Result res;
    int it = 0;
    double target;                                 // > 0: stop once a feasible solution reaches it (bestO includes the fleet penalty)
    bool reached = false;
    Tracker(const Instance& inst, const Solution& init, double lamFleet, double tgt = 0) : I(inst), lam(lamFleet), t0(now()), target(tgt) {
        best = init.nonEmpty();
        bestO = objective(I, best, lam);
        if (target > 0 && bestO <= target * (1 + 1e-9)) reached = true;
        res.history.push_back({0, 0.0, best.cost(I)});
    }
    bool offer(const Solution& s, double o) {
        if (o < bestO - EPS) {
            best = s.nonEmpty();
            bestO = o;
            if (target > 0 && bestO <= target * (1 + 1e-9)) reached = true;
            HistoryEntry h{it, now() - t0, best.cost(I)};
            if (res.history.size() > 1 && res.history.back().iter == it) res.history.back() = h;
            else res.history.push_back(h);
            return true;
        }
        return false;
    }
    bool offer(const Solution& s) { return offer(s, objective(I, s, lam)); }
    // The clock is read about every 5 ms of work: every iteration for methods whose iteration is a whole descent, at most
    // every 64 iterations for cheap random moves.  Runs limited by iterations only are not affected.
    mutable int nextCheck = 1, lastIt = 0;
    mutable double lastT = 0;
    bool timeUp(double limit) const {
        if (reached) return true;                  // target reached: every method leaves its loop at the end of the iteration
        if (limit <= 0 || it < nextCheck) return false;
        const double t = now() - t0, per = (t - lastT) / std::max(1, it - lastIt);
        const int stride = per > 0 ? std::max(1, std::min(64, (int)(0.005 / per))) : 64;
        lastT = t;
        lastIt = it;
        nextCheck = it + stride;
        return t > limit;
    }
    Result finish(int done) {
        res.best = best;
        res.cost = best.cost(I);
        res.iters = done;
        res.seconds = now() - t0;
        res.reachedTarget = reached;
        return res;
    }
};

// ===========================================================================================
// Random-move state (sa, tabu, vns): routes as lists with cached cost / load / service, six kinds of random moves
// between a random customer u and a random granular neighbour v; infeasible or degenerate moves are rejected.
// ===========================================================================================
struct WR {
    int d = 0;
    std::vector<int> nodes;
    double cost = 0, load = 0, svc = 0;
    bool live = false;
};
enum Kind { K_RELOCATE, K_OROPT, K_SWAP, K_TWOOPT, K_TWOOPTSTAR, K_NEWROUTE };
const std::vector<Kind> KINDS = {K_RELOCATE, K_OROPT, K_SWAP, K_TWOOPT, K_TWOOPTSTAR};
const std::vector<Kind> KINDS_NEW = {K_RELOCATE, K_OROPT, K_SWAP, K_TWOOPT, K_TWOOPTSTAR, K_NEWROUTE};

struct Cand {
    Kind kind = K_RELOCATE;
    double delta = 0;
    int ra = -1, rb = -1, d = -1;                          // rb = -1: the move touches one route only; d: depot of the new route
    std::vector<int> na, nb;
    std::vector<std::pair<int, int>> puts, leaves;         // (customer, route) attributes for the tabu memory (-1: a new route)
};

class Walk {
public:
    Walk(const LSContext& ctx, const Solution& sol) : I(ctx.inst), ctx_(ctx) {
        where.assign(I.n, -1);
        pos.assign(I.n, 0);
        cnt.assign(I.t, 0);
        for (size_t r = 0; r < sol.routes.size(); r++) {
            if (sol.routes[r].empty()) continue;
            WR w;
            w.d = sol.depots[r];
            w.nodes = sol.routes[r];
            w.live = true;
            refresh(w);
            R.push_back(std::move(w));
            cnt[R.back().d]++;
            index((int)R.size() - 1);
        }
        cost = recomputeCost();
    }

    const Instance& I;
    std::vector<WR> R;
    std::vector<int> where, pos, cnt, freeSlots;
    double cost = 0;

    int excess() const { return excessOf(I, cnt); }
    double recomputeCost() const {
        double c = 0;
        for (const WR& w : R) if (w.live) c += w.cost;
        return c;
    }
    Solution toSolution() const {
        Solution s;
        for (const WR& w : R) if (w.live) { s.routes.push_back(w.nodes); s.depots.push_back(w.d); }
        return s;
    }

    bool randomMove(PyRandom& rng, const std::vector<Kind>& kinds, Cand& m) {
        m = Cand();
        const int u = rng.randrange(I.n), ra = where[u], i = pos[u];
        const std::vector<int>& an = R[ra].nodes;
        const int LA = (int)an.size();
        const Kind kind = kinds[rng.randrange((int)kinds.size())];
        m.kind = kind;
        m.ra = ra;
        if (kind == K_NEWROUTE) {                                       // open a new route at one of the two nearest depots
            if (LA < 2) return false;
            std::vector<int> order(I.t);
            std::iota(order.begin(), order.end(), 0);
            std::stable_sort(order.begin(), order.end(), [&](int a, int b) { return I.d(u, I.n + a) < I.d(u, I.n + b); });
            const int d = order[rng.randrange(std::min(2, I.t))];
            if (cnt[d] >= I.fleet[d]) return false;
            m.na = an;
            m.na.erase(m.na.begin() + i);
            m.nb = {u};
            double ca, cb, q, s;
            if (!eval(R[ra].d, m.na, ca, q, s) || !eval(d, m.nb, cb, q, s)) return false;
            m.d = d;
            m.rb = -1;
            m.delta = ca + cb - R[ra].cost;
            m.puts = {{u, -1}};
            m.leaves = {{u, ra}};
            return true;
        }
        const std::vector<int>& neigh = ctx_.nbrs[u];
        if (neigh.empty()) return false;
        const int v = neigh[rng.randrange((int)neigh.size())], rb = where[v], j = pos[v];
        const std::vector<int>& bn = R[rb].nodes;
        m.rb = rb;
        if (kind == K_RELOCATE || kind == K_OROPT) {
            const int l = kind == K_RELOCATE ? 1 : 2 + rng.randrange(2);
            if (i + l > LA) return false;
            std::vector<int> seg(an.begin() + i, an.begin() + i + l), s = seg;
            if (l > 1 && rng.random() < 0.5) std::reverse(s.begin(), s.end());
            const int after = rng.random() < 0.5 ? 1 : 0;
            if (rb == ra) {
                if (j >= i && j < i + l) return false;
                std::vector<int> rest(an.begin(), an.begin() + i);
                rest.insert(rest.end(), an.begin() + i + l, an.end());
                const int jj = j > i ? j - l : j, p = jj + after;
                m.na.assign(rest.begin(), rest.begin() + p);
                m.na.insert(m.na.end(), s.begin(), s.end());
                m.na.insert(m.na.end(), rest.begin() + p, rest.end());
                m.rb = -1;
            } else {
                m.na.assign(an.begin(), an.begin() + i);
                m.na.insert(m.na.end(), an.begin() + i + l, an.end());
                const int p = j + after;
                m.nb.assign(bn.begin(), bn.begin() + p);
                m.nb.insert(m.nb.end(), s.begin(), s.end());
                m.nb.insert(m.nb.end(), bn.begin() + p, bn.end());
                for (int c : seg) { m.puts.push_back({c, rb}); m.leaves.push_back({c, ra}); }
            }
        } else if (kind == K_SWAP) {
            if (rb == ra) {
                if (i == j) return false;
                m.na = an;
                m.na[i] = v;
                m.na[j] = u;
                m.rb = -1;
            } else {
                m.na = an;
                m.nb = bn;
                m.na[i] = v;
                m.nb[j] = u;
                m.puts = {{u, rb}, {v, ra}};
                m.leaves = {{u, ra}, {v, rb}};
            }
        } else if (kind == K_TWOOPT) {
            if (rb != ra || i == j) return false;
            const int lo = std::min(i, j), hi = std::max(i, j);
            if (hi - lo < 2) return false;
            m.na.assign(an.begin(), an.begin() + lo + 1);
            m.na.insert(m.na.end(), an.rend() - hi - 1, an.rend() - lo - 1);       // reversed an[lo+1 .. hi]
            m.na.insert(m.na.end(), an.begin() + hi + 1, an.end());
            m.rb = -1;
        } else {                                                                    // K_TWOOPTSTAR: exchange the tails of two routes
            if (rb == ra) return false;
            if (R[ra].d == R[rb].d && rng.random() < 0.5) {                         // reversed variant only for a common depot
                m.na.assign(an.begin(), an.begin() + i + 1);
                m.na.insert(m.na.end(), bn.rend() - j - 1, bn.rend());              // reversed bn[0 .. j]
                m.nb.assign(an.rbegin(), an.rend() - i - 1);                        // reversed an[i+1 ..]
                m.nb.insert(m.nb.end(), bn.begin() + j + 1, bn.end());
            } else {
                m.na.assign(an.begin(), an.begin() + i + 1);
                m.na.insert(m.na.end(), bn.begin() + j + 1, bn.end());
                m.nb.assign(bn.begin(), bn.begin() + j + 1);
                m.nb.insert(m.nb.end(), an.begin() + i + 1, an.end());
            }
        }
        double ca, cb = 0, q, s;
        if (!eval(R[ra].d, m.na, ca, q, s)) return false;
        if (m.rb >= 0 && !eval(R[m.rb].d, m.nb, cb, q, s)) return false;
        m.delta = ca + cb - R[ra].cost - (m.rb >= 0 ? R[m.rb].cost : 0.0);
        return true;
    }

    void apply(const Cand& m) {
        if (m.kind == K_NEWROUTE) {
            setRoute(m.ra, m.na);
            int slot;
            if (!freeSlots.empty()) { slot = freeSlots.back(); freeSlots.pop_back(); }
            else { slot = (int)R.size(); R.emplace_back(); }
            WR& w = R[slot];
            w.d = m.d;
            w.nodes = m.nb;
            w.live = true;
            refresh(w);
            cnt[m.d]++;
            index(slot);
        } else {
            setRoute(m.ra, m.na);
            if (m.rb >= 0) setRoute(m.rb, m.nb);
        }
        cost += m.delta;
    }

private:
    const LSContext& ctx_;

    void refresh(WR& w) {
        w.cost = w.nodes.empty() ? 0.0 : tourCost(I, I.n + w.d, w.nodes);
        w.load = w.svc = 0;
        for (int c : w.nodes) { w.load += I.demand[c]; w.svc += I.service[c]; }
    }
    // cost of the route and feasibility (capacity, duration = travel + service) at depot d
    bool eval(int d, const std::vector<int>& nodes, double& c, double& q, double& s) const {
        c = nodes.empty() ? 0.0 : tourCost(I, I.n + d, nodes);
        q = s = 0;
        for (int v : nodes) { q += I.demand[v]; s += I.service[v]; }
        return q <= I.capacity[d] + EPS && c + s <= I.maxDuration[d] + EPS;
    }
    void index(int r) {
        const std::vector<int>& nodes = R[r].nodes;
        for (int k = 0; k < (int)nodes.size(); k++) { where[nodes[k]] = r; pos[nodes[k]] = k; }
    }
    void setRoute(int r, const std::vector<int>& nodes) {
        WR& w = R[r];
        w.nodes = nodes;
        refresh(w);
        if (nodes.empty()) { w.live = false; cnt[w.d]--; freeSlots.push_back(r); }
        else index(r);
    }
};

// apply s random feasible moves regardless of their cost (shaking); returns the customers of the touched routes, in order of first touch
std::vector<int> shake(Walk& w, int s, PyRandom& rng, const std::vector<Kind>& kinds) {
    std::vector<char> seen(w.I.n, 0);
    std::vector<int> touched;
    Cand m;
    int done = 0, tries = 0;
    while (done < s && tries < 30 * s) {
        tries++;
        if (!w.randomMove(rng, kinds, m)) continue;
        w.apply(m);
        done++;
        for (const auto* lst : {&m.na, &m.nb})
            for (int c : *lst) if (!seen[c]) { seen[c] = 1; touched.push_back(c); }
    }
    return touched;
}

// ---------------------------------------------------------------------------------------------
Result runSA(const Instance& I, const MethodOptions& o, const LSContext& ctx, const Solution& init) {
    PyRandom rng(o.seed);
    Walk w(ctx, init);
    Tracker T(I, init, o.lamFleet, o.target);
    const int iters = o.iters, cycles = 3, len = (iters + cycles - 1) / cycles;
    const double t0 = 0.1, tf = 0.002, reheat = 0.5, avg = w.cost / (I.n + init.nRoutes());
    double cur = w.cost;
    Cand m;
    int it;
    for (it = 1; it <= iters; it++) {
        T.it = it;
        if (w.randomMove(rng, KINDS_NEW, m)) {
            bool take = m.delta <= 0;
            if (!take) {                                                // geometric cooling in `cycles` cycles, each reheated to reheat x the previous start
                const int c = (it - 1) / len, k = (it - 1) - c * len;
                const double T0 = t0 * avg * std::pow(reheat, c), Tf = std::min(tf * avg, T0 * 0.05);
                take = rng.random() < std::exp(-m.delta / (T0 * std::pow(Tf / T0, (double)k / len)));
            }
            if (take) {
                w.apply(m);
                cur += m.delta;
                if (cur + o.lamFleet * w.excess() < T.bestO - 1e-7) {
                    Solution s = w.toSolution();
                    T.offer(s);
                    w.cost = w.recomputeCost();                          // re-synchronise the running cost
                    cur = w.cost;
                }
            }
        }
        if (T.timeUp(o.timeLimit)) break;
    }
    return T.finish(std::min(it, iters));
}

Result runTabu(const Instance& I, const MethodOptions& o, const LSContext& ctx, const Solution& init) {
    PyRandom rng(o.seed);
    Walk w(ctx, init);
    Tracker T(I, init, o.lamFleet, o.target);
    const int cs = 30, tenure = 8, spread = 8;
    std::unordered_map<long long, int> tabu;
    auto key = [](int c, int r) { return (long long)c * 1000003LL + r; };
    double cur = w.cost;
    Cand m, best;
    int it;
    for (it = 1; it <= o.iters; it++) {
        T.it = it;
        bool have = false;
        for (int q = 0; q < cs; q++) {                                  // randomised candidate list: the best admissible of cs random moves
            if (!w.randomMove(rng, KINDS_NEW, m)) continue;
            bool isTabu = false;
            for (const auto& pr : m.puts) {
                auto f = tabu.find(key(pr.first, pr.second));
                if (f != tabu.end() && f->second > it) { isTabu = true; break; }
            }
            if (isTabu && !(cur + m.delta < T.bestO - EPS)) continue;   // aspiration: a tabu move is allowed if it beats the best solution
            if (!have || m.delta < best.delta) { best = m; have = true; }
        }
        if (have) {
            w.apply(best);
            cur += best.delta;
            for (const auto& lv : best.leaves) tabu[key(lv.first, lv.second)] = it + tenure + rng.randrange(spread + 1);
            if (tabu.size() > 60000)
                for (auto f = tabu.begin(); f != tabu.end();) f = f->second <= it ? tabu.erase(f) : std::next(f);
            if (cur + o.lamFleet * w.excess() < T.bestO - 1e-7) {
                Solution s = w.toSolution();
                T.offer(s);
                w.cost = w.recomputeCost();
                cur = w.cost;
            }
        }
        if (T.timeUp(o.timeLimit)) break;
    }
    return T.finish(std::min(it, o.iters));
}

std::string presetOf(const Instance& I, const MethodOptions& o) { return o.preset.empty() ? (I.t > 1 ? "full" : "hgs") : o.preset; }

// granular local search descent from `s` (all customers, or only `dirty`)
// tunable parameters: value given in MethodOptions::params, else the default
double par(const MethodOptions& o, const char* name, double def) { auto f = o.params.find(name); return f == o.params.end() ? def : f->second; }
int ipar(const MethodOptions& o, const char* name, int def) { return (int)std::lround(par(o, name, def)); }

// Per-run accounting of the local-search operators (set by runMethod; the library is single-threaded per process).
std::array<long, OP_COUNT>* g_opStats = nullptr;
std::array<double, OP_COUNT>* g_opGain = nullptr;
void account(const LSState& ls) {
    if (!g_opStats) return;
    for (int k = 0; k < OP_COUNT; k++) { (*g_opStats)[k] += ls.stats[k]; (*g_opGain)[k] += ls.gain[k]; }
}

Solution descend(const LSContext& ctx, const Solution& s, const MethodOptions& o, PyRandom& rng, const std::vector<int>* dirty = nullptr) {
    LSState ls(ctx, s, o.lamFleet);
    RunOptions ro;
    const Preset p = preset(presetOf(ctx.inst, o));
    ro.ops = p.ops & ~o.excludeOps;
    if (ro.ops == 0) return s.nonEmpty();                           // every operator switched off
    ro.relocLen = p.relocLen;
    ro.swapLen = p.swapLen;
    ro.dirty = dirty;
    ls.run(ro, rng);
    account(ls);
    return ls.toSolution();
}

Result runVNS(const Instance& I, const MethodOptions& o, const LSContext& ctx, const Solution& init) {
    PyRandom rng(o.seed);
    Tracker T(I, init, o.lamFleet, o.target);
    const int kmax = 5;
    Solution cur = init.nonEmpty();
    double curO = objective(I, cur, o.lamFleet);
    int k = 1, it;
    for (it = 1; it <= o.iters; it++) {
        T.it = it;
        Walk w(ctx, cur);
        std::vector<int> touched = shake(w, 2 * k, rng, KINDS);          // shaking with 2k random moves, then repair by local search
        Solution cand = touched.empty() ? cur : descend(ctx, w.toSolution(), o, rng, &touched);
        const double oc = objective(I, cand, o.lamFleet);
        if (oc < curO - EPS) { cur = cand; curO = oc; k = 1; T.offer(cand, oc); }
        else k = k % kmax + 1;
        if (T.timeUp(o.timeLimit)) break;
    }
    return T.finish(std::min(it, o.iters));
}

// ===========================================================================================
// Destroy and repair on route lists (ALNS)
// ===========================================================================================
struct PR {
    int d = 0;
    std::vector<int> nodes;
    double cost = 0, load = 0, svc = 0;
};
struct Part {
    const Instance& I;
    std::vector<PR> routes;
    std::vector<int> cnt, where;
    explicit Part(const Instance& inst, const Solution& s) : I(inst) {
        cnt.assign(I.t, 0);
        where.assign(I.n, -1);
        for (size_t r = 0; r < s.routes.size(); r++) {
            if (s.routes[r].empty()) continue;
            PR p;
            p.d = s.depots[r];
            p.nodes = s.routes[r];
            fix(p);
            routes.push_back(std::move(p));
        }
        reindex();
    }
    void fix(PR& p) const {
        p.cost = p.nodes.empty() ? 0.0 : tourCost(I, I.n + p.d, p.nodes);
        p.load = p.svc = 0;
        for (int c : p.nodes) { p.load += I.demand[c]; p.svc += I.service[c]; }
    }
    void reindex() {                                                    // drop emptied routes, rebuild the counters
        routes.erase(std::remove_if(routes.begin(), routes.end(), [](const PR& p) { return p.nodes.empty(); }), routes.end());
        std::fill(cnt.begin(), cnt.end(), 0);
        for (size_t k = 0; k < routes.size(); k++) {
            cnt[routes[k].d]++;
            for (int c : routes[k].nodes) where[c] = (int)k;
        }
    }
    void remove(int c) {
        PR& p = routes[where[c]];
        p.nodes.erase(std::find(p.nodes.begin(), p.nodes.end(), c));
        fix(p);
        where[c] = -1;
    }
    Solution toSolution() const {
        Solution s;
        for (const PR& p : routes) if (!p.nodes.empty()) { s.routes.push_back(p.nodes); s.depots.push_back(p.d); }
        return s;
    }
};

double removalGain(const Instance& I, const PR& R, int i) {
    const int dn = I.n + R.d, c = R.nodes[i], pv = i ? R.nodes[i - 1] : dn, nx = i + 1 < (int)R.nodes.size() ? R.nodes[i + 1] : dn;
    return I.d(pv, c) + I.d(c, nx) - I.d(pv, nx);
}

std::vector<int> destroyRandom(Part& P, int q, PyRandom& rng) {
    std::vector<int> all(P.I.n);
    std::iota(all.begin(), all.end(), 0);
    rng.shuffle(all);
    std::vector<int> rm(all.begin(), all.begin() + std::min<int>(q, (int)all.size()));
    for (int c : rm) P.remove(c);
    return rm;
}
std::vector<int> destroyWorst(Part& P, int q, PyRandom& rng) {         // largest removal gain first, randomised with p = 3 (Ropke & Pisinger)
    std::vector<int> rm;
    for (int k = 0; k < q; k++) {
        std::vector<std::pair<double, int>> cand;
        for (const PR& R : P.routes)
            for (int i = 0; i < (int)R.nodes.size(); i++) cand.push_back({removalGain(P.I, R, i), R.nodes[i]});
        if (cand.empty()) break;
        std::stable_sort(cand.begin(), cand.end(), [](const auto& a, const auto& b) { return a.first > b.first; });
        const int idx = std::min((int)cand.size() - 1, (int)std::floor(std::pow(rng.random(), 3.0) * cand.size()));
        const int c = cand[idx].second;
        P.remove(c);
        rm.push_back(c);
    }
    return rm;
}
std::vector<int> destroyRelated(Part& P, int q, PyRandom& rng, const LSContext& ctx) {   // Shaw removal by distance, randomised with p = 3
    std::vector<int> rm = {rng.randrange(P.I.n)};
    std::vector<char> gone(P.I.n, 0);
    gone[rm[0]] = 1;
    while ((int)rm.size() < q) {
        const int r = rm[rng.randrange((int)rm.size())];
        std::vector<int> cand;
        for (int c : ctx.adj[r]) {
            if (gone[c]) continue;
            cand.push_back(c);
            if (cand.size() >= 24) break;
        }
        if (cand.empty()) break;
        const int c = cand[std::min((int)cand.size() - 1, (int)std::floor(std::pow(rng.random(), 3.0) * cand.size()))];
        gone[c] = 1;
        rm.push_back(c);
    }
    for (int c : rm) P.remove(c);
    return rm;
}

// cheapest feasible insertion of c into route R: position and delta, or false
bool bestPos(const Instance& I, const PR& R, int c, int& pos, double& delta) {
    if (R.load + I.demand[c] > I.capacity[R.d] + EPS) return false;
    const int dn = I.n + R.d, L = (int)R.nodes.size();
    int prev = dn;
    double bd = INF;
    int bp = -1;
    for (int p = 0; p <= L; p++) {
        const int nx = p < L ? R.nodes[p] : dn;
        const double dl = I.d(prev, c) + I.d(c, nx) - I.d(prev, nx);
        if (dl < bd && R.cost + dl + R.svc + I.service[c] <= I.maxDuration[R.d] + EPS) { bd = dl; bp = p; }
        prev = nx;
    }
    if (bp < 0) return false;
    pos = bp;
    delta = bd;
    return true;
}
void openRoute(Part& P, int c) {                                        // new route at the nearest depot that still has a vehicle
    const Instance& I = P.I;
    std::vector<int> o(I.t);
    std::iota(o.begin(), o.end(), 0);
    std::stable_sort(o.begin(), o.end(), [&](int a, int b) { return I.d(c, I.n + a) < I.d(c, I.n + b); });
    int d = o[0];
    for (int x : o) if (P.cnt[x] < I.fleet[x]) { d = x; break; }
    PR R;
    R.d = d;
    R.nodes = {c};
    P.fix(R);
    P.routes.push_back(std::move(R));
    P.cnt[d]++;
    P.where[c] = (int)P.routes.size() - 1;
}
void insertAt(Part& P, int r, int p, int c) {
    PR& R = P.routes[r];
    R.nodes.insert(R.nodes.begin() + p, c);
    P.fix(R);
    P.where[c] = r;
}
void repairGreedy(Part& P, std::vector<int> rm, PyRandom& rng) {         // random order, cheapest feasible insertion over all routes of all depots
    rng.shuffle(rm);
    for (int c : rm) {
        int br = -1, bp = -1;
        double bd = INF;
        for (int r = 0; r < (int)P.routes.size(); r++) {
            int p;
            double dl;
            if (bestPos(P.I, P.routes[r], c, p, dl) && dl < bd) { bd = dl; br = r; bp = p; }
        }
        if (br >= 0) insertAt(P, br, bp, c); else openRoute(P, c);
    }
}
void repairRegret2(Part& P, std::vector<int> left) {                    // insert the customer whose best and second-best routes differ most
    while (!left.empty()) {
        int pc = -1, pr = -1, pp = -1;
        double bestReg = -INF;
        int firstStuck = -1;
        for (int k = 0; k < (int)left.size(); k++) {
            const int c = left[k];
            double b1 = INF, b2 = INF;
            int r1 = -1, p1 = -1;
            for (int r = 0; r < (int)P.routes.size(); r++) {
                int p;
                double dl;
                if (!bestPos(P.I, P.routes[r], c, p, dl)) continue;
                if (dl < b1) { b2 = b1; b1 = dl; r1 = r; p1 = p; } else if (dl < b2) b2 = dl;
            }
            if (r1 < 0) { if (firstStuck < 0) firstStuck = k; continue; }
            const double reg = b2 == INF ? 1e9 - b1 : b2 - b1;
            if (reg > bestReg) { bestReg = reg; pc = k; pr = r1; pp = p1; }
        }
        if (pc < 0) pc = firstStuck;                                     // nothing fits anywhere: open a route for the first such customer
        const int c = left[pc];
        left.erase(left.begin() + pc);
        if (pr >= 0 && bestReg > -INF) insertAt(P, pr, pp, c); else openRoute(P, c);
    }
}

Result runALNS(const Instance& I, const MethodOptions& o, LSContext& ctx, const Solution& init) {
    ctx.ensureAdjacency();
    PyRandom rng(o.seed);
    Tracker T(I, init, o.lamFleet, o.target);
    const double t0f = 0.1, tff = 0.001, rho = 0.1, qlo = 0.05, qhi = 0.2;
    const int seg = 100, s1 = 33, s2 = 9, s3 = 3;
    const int qa = std::max(2, (int)std::ceil(qlo * I.n)), qb = std::max(3, std::min(60, (int)std::ceil(qhi * I.n)));
    Solution cur = init.nonEmpty();
    double curO = objective(I, cur, o.lamFleet);
    const double avg = cur.cost(I) / (I.n + cur.nRoutes()), Tstart = t0f * avg, Tend = tff * avg;
    const double cool = std::pow(Tend / Tstart, 1.0 / std::max(1, o.iters));
    double temp = Tstart;
    std::vector<double> wd(3, 1.0), wr(2, 1.0), scd(3, 0.0), scr(2, 0.0), usd(3, 0.0), usr(2, 0.0);
    std::vector<long> totD(3, 0), totR(2, 0);
    auto roulette = [&](const std::vector<double>& ws) {
        double sum = 0;
        for (double x : ws) sum += x;
        double x = rng.random() * sum;
        for (int i = 0; i < (int)ws.size(); i++) { x -= ws[i]; if (x <= 0) return i; }
        return (int)ws.size() - 1;
    };
    int it;
    for (it = 1; it <= o.iters; it++) {
        T.it = it;
        const int di = roulette(wd), ri = roulette(wr), q = qa + rng.randrange(qb - qa + 1);
        Part P(I, cur);
        std::vector<int> rm = di == 0 ? destroyRandom(P, q, rng) : di == 1 ? destroyWorst(P, q, rng) : destroyRelated(P, q, rng, ctx);
        P.reindex();
        if (ri == 0) repairGreedy(P, rm, rng); else repairRegret2(P, rm);
        P.reindex();
        Solution cand = P.toSolution();
        const double oc = objective(I, cand, o.lamFleet);
        double score = 0;
        if (oc < curO - temp * std::log(1.0 - rng.random())) {           // simulated-annealing acceptance
            if (oc < T.bestO - EPS) { score = s1; T.offer(cand, oc); } else if (oc < curO - EPS) score = s2; else score = s3;
            cur = cand;
            curO = oc;
        }
        scd[di] += score; scr[ri] += score; usd[di]++; usr[ri]++; totD[di]++; totR[ri]++;
        if (it % seg == 0) {                                             // adapt the roulette weights from the scores of the last segment
            for (int i = 0; i < 3; i++) { if (usd[i] > 0) wd[i] = std::max(0.05, (1 - rho) * wd[i] + rho * scd[i] / usd[i]); scd[i] = usd[i] = 0; }
            for (int i = 0; i < 2; i++) { if (usr[i] > 0) wr[i] = std::max(0.05, (1 - rho) * wr[i] + rho * scr[i] / usr[i]); scr[i] = usr[i] = 0; }
        }
        temp *= cool;
        if (T.timeUp(o.timeLimit)) break;
    }
    Result res = T.finish(std::min(it, o.iters));
    const char* dn[] = {"removal:random", "removal:worst", "removal:related"};
    const char* rn[] = {"insertion:greedy", "insertion:regret2"};
    for (int i = 0; i < 3; i++) res.extra.push_back({dn[i], totD[i], wd[i]});
    for (int i = 0; i < 2; i++) res.extra.push_back({rn[i], totR[i], wr[i]});
    return res;
}

}  // namespace

// ===========================================================================================
// Giant tours: optimal split (Prins 2004), order crossover, memetic GA
// ===========================================================================================
std::vector<std::vector<int>> splitTour(const Instance& I, int d, const std::vector<int>& tour) {
    const int m = (int)tour.size(), dn = I.n + d;
    const double Q = I.capacity[d], Dm = I.maxDuration[d];
    std::vector<double> f(m + 1, INF);
    std::vector<int> pred(m + 1, -1);
    f[0] = 0;
    for (int i = 0; i < m; i++) {
        if (f[i] == INF) continue;
        double q = 0, sv = 0, c = 0;
        for (int j = i; j < m; j++) {
            const int v = tour[j];
            q += I.demand[v];
            sv += I.service[v];
            c += j == i ? I.d(dn, v) : I.d(tour[j - 1], v);
            if (q > Q + EPS) break;
            const double full = c + I.d(v, dn);
            if (full + sv > Dm + EPS && j > i) break;
            if (f[i] + full < f[j + 1]) { f[j + 1] = f[i] + full; pred[j + 1] = i; }
        }
    }
    std::vector<std::vector<int>> routes;
    if (f[m] == INF) {                                                   // a single customer is infeasible on its own: one route each
        for (int c : tour) routes.push_back({c});
        return routes;
    }
    for (int j = m; j > 0;) {
        const int i = pred[j];
        routes.emplace_back(tour.begin() + i, tour.begin() + j);
        j = i;
    }
    std::reverse(routes.begin(), routes.end());
    return routes;
}

namespace {
using DepotOrder = std::vector<std::vector<int>>;                        // for every customer the depots by distance

DepotOrder buildDepotOrder(const Instance& I) {
    DepotOrder o(I.n, std::vector<int>(I.t));
    for (int c = 0; c < I.n; c++) {
        std::iota(o[c].begin(), o[c].end(), 0);
        std::stable_sort(o[c].begin(), o[c].end(), [&](int a, int b) { return I.d(c, I.n + a) < I.d(c, I.n + b); });
    }
    return o;
}
// depot genes: 0 = nearest depot, 1 = second nearest
std::vector<char> depGenes(const Instance& I, const DepotOrder& dord, const Solution& s) {
    std::vector<char> g(I.n, 0);
    if (I.t > 1)
        for (size_t r = 0; r < s.routes.size(); r++)
            for (int c : s.routes[r]) g[c] = dord[c].size() > 1 && dord[c][1] == s.depots[r] ? 1 : 0;
    return g;
}
// giant tour of a solution: routes ordered by depot, then by the polar angle of their centroid around the depot
std::vector<int> giantTour(const Instance& I, const Solution& s) {
    struct K { int d; double ang; size_t r; };
    std::vector<K> ks;
    for (size_t r = 0; r < s.routes.size(); r++) {
        if (s.routes[r].empty()) continue;
        const int dn = I.n + s.depots[r];
        double sx = 0, sy = 0;
        for (int c : s.routes[r]) { sx += I.x[c]; sy += I.y[c]; }
        const double L = (double)s.routes[r].size();
        ks.push_back({s.depots[r], std::atan2(sy / L - I.y[dn], sx / L - I.x[dn]), r});
    }
    std::stable_sort(ks.begin(), ks.end(), [](const K& a, const K& b) { return a.d != b.d ? a.d < b.d : a.ang < b.ang; });
    std::vector<int> t;
    for (const K& k : ks) t.insert(t.end(), s.routes[k.r].begin(), s.routes[k.r].end());
    return t;
}
Solution decodeChrom(const Instance& I, const DepotOrder& dord, const std::vector<int>& tour, const std::vector<char>& dep) {
    std::vector<std::vector<int>> per(I.t);
    for (int c : tour) {
        const auto& o = dord[c];
        per[I.t > 1 && o.size() > 1 && dep[c] ? o[1] : o[0]].push_back(c);
    }
    Solution s;
    for (int d = 0; d < I.t; d++) {
        if (per[d].empty()) continue;
        for (auto& r : splitTour(I, d, per[d])) { s.routes.push_back(std::move(r)); s.depots.push_back(d); }
    }
    return repairFleetLimit(I, s);
}
std::vector<int> orderCrossover(const std::vector<int>& a, const std::vector<int>& b, PyRandom& rng) {
    const int n = (int)a.size();
    int i = rng.randrange(n), j = rng.randrange(n);
    if (i > j) std::swap(i, j);
    std::vector<int> child(n, -1);
    std::vector<char> used(n, 0);
    for (int k = i; k <= j; k++) { child[k] = a[k]; used[a[k]] = 1; }
    int p = (j + 1) % n;
    for (int q = 0; q < n; q++) {
        const int c = b[(j + 1 + q) % n];
        if (!used[c]) { child[p] = c; p = (p + 1) % n; }
    }
    return child;
}

struct Ind {
    Solution sol;
    double o = 0;
    std::vector<int> tour;
    std::vector<char> dep;
};

Result runGA(const Instance& I, const MethodOptions& o, const LSContext& ctx, const Solution& init) {
    PyRandom rng(o.seed);
    Tracker T(I, init, o.lamFleet, o.target);
    const DepotOrder dord = buildDepotOrder(I);
    const int popSize = 20, tourn = 2;
    const double pmut = 0.1;
    auto mk = [&](const Solution& s) { Ind x; x.sol = s; x.o = objective(I, s, o.lamFleet); x.tour = giantTour(I, s); x.dep = depGenes(I, dord, s); return x; };
    std::vector<Ind> pop = {mk(init.nonEmpty())};
    auto pick = [&]() -> const Ind* {
        const Ind* b = &pop[rng.randrange((int)pop.size())];
        for (int k = 1; k < tourn; k++) { const Ind* c = &pop[rng.randrange((int)pop.size())]; if (c->o < b->o) b = c; }
        return b;
    };
    int it;
    for (it = 1; it <= o.iters; it++) {
        T.it = it;
        if ((int)pop.size() < popSize) {                                 // initial population: random giant tours, random depot genes, descent
            std::vector<int> tour(I.n);
            std::iota(tour.begin(), tour.end(), 0);
            rng.shuffle(tour);
            std::vector<char> dep(I.n, 0);
            if (I.t > 1) for (int c = 0; c < I.n; c++) dep[c] = rng.random() >= 0.75 ? 1 : 0;
            Solution s = descend(ctx, decodeChrom(I, dord, tour, dep), o, rng);
            T.offer(s);
            pop.push_back(mk(s));
        } else {
            const Ind* p1 = pick();
            const Ind* p2 = pick();
            for (int t = 0; t < 5 && p2 == p1; t++) p2 = pick();
            std::vector<int> tour = orderCrossover(p1->tour, p2->tour, rng);
            std::vector<char> dep(I.n);
            for (int c = 0; c < I.n; c++) dep[c] = rng.random() < 0.5 ? p1->dep[c] : p2->dep[c];
            if (rng.random() < pmut) {
                const int i = rng.randrange(I.n), j = rng.randrange(I.n);
                std::swap(tour[i], tour[j]);
                if (I.t > 1) dep[tour[i]] ^= 1;
            }
            Solution s = descend(ctx, decodeChrom(I, dord, tour, dep), o, rng);
            T.offer(s);
            Ind ind = mk(s);
            int w = 0;
            for (int i = 1; i < (int)pop.size(); i++) if (pop[i].o > pop[w].o) w = i;
            bool dup = false;
            for (const Ind& x : pop) if (std::fabs(x.o - ind.o) < 1e-7) { dup = true; break; }
            if (ind.o < pop[w].o - EPS && !dup) pop[w] = std::move(ind);
        }
        if (T.timeUp(o.timeLimit)) break;
    }
    return T.finish(std::min(it, o.iters));
}

// ===========================================================================================
// The remaining thirteen methods of the browser Solver (docs/assets/methods.js): ta, rtr, grasp, msls, lns, ig,
// aco, grasppr, scatter, brkga, pso, de, rkrs.  Same schemes and default parameters as the JavaScript versions.
// ===========================================================================================

// ---- threshold accepting (Dueck & Scheuer 1990) and record-to-record travel (Dueck 1993): the walk of runSA with another acceptance rule
Result runWalkRule(const Instance& I, const MethodOptions& o, const LSContext& ctx, const Solution& init, bool recordToRecord) {
    PyRandom rng(o.seed);
    Walk w(ctx, init);
    Tracker T(I, init, o.lamFleet, o.target);
    const double th0 = par(o, "th0", 0.05), dev0 = par(o, "dev0", 0.02), avg = w.cost / (I.n + init.nRoutes());
    double cur = w.cost;
    Cand m;
    int it;
    for (it = 1; it <= o.iters; it++) {
        T.it = it;
        if (w.randomMove(rng, KINDS_NEW, m)) {
            const double left = 1.0 - (double)(it - 1) / o.iters;                    // the threshold / the band shrinks linearly to zero
            const bool take = recordToRecord ? cur + m.delta <= T.bestO * (1 + dev0 * left) + 1e-9 : m.delta <= th0 * avg * left;
            if (take) {
                w.apply(m);
                cur += m.delta;
                if (cur + o.lamFleet * w.excess() < T.bestO - 1e-7) {
                    Solution s = w.toSolution();
                    T.offer(s);
                    w.cost = w.recomputeCost();
                    cur = w.cost;
                }
            }
        }
        if (T.timeUp(o.timeLimit)) break;
    }
    return T.finish(std::min(it, o.iters));
}

// ---- random keys: n customer priorities, then (several depots) n depot keys; a customer goes to its nearest depot unless its
//      depot key is >= 0.75 (second nearest); the customers of a depot sorted by key form a giant tour that is split optimally
using Keys = std::vector<double>;
int keyLen(const Instance& I) { return I.t > 1 ? 2 * I.n : I.n; }
Keys randomKeys(const Instance& I, PyRandom& rng) {
    Keys k(keyLen(I));
    for (double& x : k) x = rng.random();
    return k;
}
Solution decodeKeys(const Instance& I, const DepotOrder& dord, const Keys& keys) {
    std::vector<std::vector<int>> per(I.t);
    for (int c = 0; c < I.n; c++) {
        const auto& o = dord[c];
        per[I.t > 1 && o.size() > 1 && keys[I.n + c] >= 0.75 ? o[1] : o[0]].push_back(c);
    }
    Solution s;
    for (int d = 0; d < I.t; d++) {
        if (per[d].empty()) continue;
        std::stable_sort(per[d].begin(), per[d].end(), [&](int a, int b) { return keys[a] < keys[b]; });
        for (auto& r : splitTour(I, d, per[d])) { s.routes.push_back(std::move(r)); s.depots.push_back(d); }
    }
    return repairFleetLimit(I, s);
}
// a permutation of all customers -> solution (nearest depot, optimal split)
Solution decodeTour(const Instance& I, const DepotOrder& dord, const std::vector<int>& tour) {
    return decodeChrom(I, dord, tour, std::vector<char>(I.n, 0));
}
// encode a solution as keys: position in the giant tour, depot keys by nearest / second nearest
Keys seedKeys(const Instance& I, const DepotOrder& dord, const Solution& s) {
    Keys k(keyLen(I), 0.0);
    const std::vector<int> t = giantTour(I, s);
    const std::vector<char> g = depGenes(I, dord, s);
    for (int p = 0; p < (int)t.size(); p++) k[t[p]] = (p + 0.5) / I.n;
    if (I.t > 1) for (int c = 0; c < I.n; c++) k[I.n + c] = g[c] ? 0.9 : 0.1;
    return k;
}
Solution randomStart(const Instance& I, const LSContext& ctx, const DepotOrder& dord, const MethodOptions& o, PyRandom& rng) {
    return descend(ctx, decodeKeys(I, dord, randomKeys(I, rng)), o, rng);
}

// ---- randomised nearest-neighbour construction (GRASP): the next customer is drawn from the best alpha-fraction of the feasible ones
Solution graspConstruct(const Instance& I, const DepotOrder& dord, PyRandom& rng, double alpha) {
    std::vector<std::vector<int>> todo(I.t);
    for (int c = 0; c < I.n; c++) todo[dord[c][0]].push_back(c);
    Solution s;
    for (int d = 0; d < I.t; d++) {
        const int dn = I.n + d;
        std::vector<int>& left = todo[d];
        while (!left.empty()) {
            std::vector<int> cur;
            double q = 0, dur = 0;
            int pos = dn;
            for (;;) {
                std::vector<std::pair<double, int>> cand;
                for (int k = 0; k < (int)left.size(); k++) {
                    const int c = left[k];
                    if (q + I.demand[c] > I.capacity[d] + EPS || dur + I.d(pos, c) + I.service[c] + I.d(c, dn) > I.maxDuration[d] + EPS) continue;
                    cand.push_back({I.d(pos, c), k});
                }
                if (cand.empty()) break;
                std::stable_sort(cand.begin(), cand.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
                const int rcl = std::max(1, (int)std::ceil(alpha * cand.size()));
                const int k = cand[rng.randrange(rcl)].second, c = left[k];
                left.erase(left.begin() + k);
                cur.push_back(c);
                q += I.demand[c];
                dur += I.d(pos, c) + I.service[c];
                pos = c;
            }
            if (cur.empty()) { cur.push_back(left.back()); left.pop_back(); }       // a customer that fits no route on its own
            s.routes.push_back(std::move(cur));
            s.depots.push_back(d);
        }
    }
    return repairFleetLimit(I, s);
}

// ---- multi-start: grasp (Feo & Resende 1995) and msls (random decodes + local search)
Result runMultiStart(const Instance& I, const MethodOptions& o, const LSContext& ctx, const Solution& init, bool grasp) {
    PyRandom rng(o.seed);
    Tracker T(I, init, o.lamFleet, o.target);
    const DepotOrder dord = buildDepotOrder(I);
    int it;
    for (it = 1; it <= o.iters; it++) {
        T.it = it;
        T.offer(grasp ? descend(ctx, graspConstruct(I, dord, rng, par(o, "alpha", 0.2)), o, rng) : randomStart(I, ctx, dord, o, rng));
        if (T.timeUp(o.timeLimit)) break;
    }
    return T.finish(std::min(it, o.iters));
}

// ---- destroy & repair with the operators of the ALNS: di = 0 random, 1 worst, 2 related; ri = 0 greedy, 1 regret-2
Solution destroyRepair(const Instance& I, const LSContext& ctx, PyRandom& rng, const Solution& cur, int di, int ri, int q) {
    Part P(I, cur);
    std::vector<int> rm = di == 0 ? destroyRandom(P, q, rng) : di == 1 ? destroyWorst(P, q, rng) : destroyRelated(P, q, rng, ctx);
    P.reindex();
    if (ri == 0) repairGreedy(P, rm, rng); else repairRegret2(P, rm);
    P.reindex();
    return P.toSolution();
}
// lns (Shaw 1998): operators chosen uniformly, only non-worsening solutions are accepted
// ig  (Ruiz & Stuetzle 2007, here for routing): remove a fixed fraction at random, reinsert greedily, accept within a band
Result runLnsIg(const Instance& I, const MethodOptions& o, LSContext& ctx, const Solution& init, bool ig) {
    ctx.ensureAdjacency();
    PyRandom rng(o.seed);
    Tracker T(I, init, o.lamFleet, o.target);
    const int qa = std::max(2, (int)std::ceil(par(o, "qlo", 0.05) * I.n)), qb = std::max(qa + 1, std::min(60, (int)std::ceil(par(o, "qhi", 0.2) * I.n)));
    const int qig = std::max(2, (int)std::floor(0.1 * I.n + 0.5));
    const double band = 0.002;
    Solution cur = init.nonEmpty();
    double curO = objective(I, cur, o.lamFleet);
    int it;
    for (it = 1; it <= o.iters; it++) {
        T.it = it;
        Solution cand;
        if (ig) cand = destroyRepair(I, ctx, rng, cur, 0, 0, qig);
        else {
            const int di = rng.randrange(3), ri = rng.randrange(2), q = qa + rng.randrange(qb - qa + 1);
            cand = destroyRepair(I, ctx, rng, cur, di, ri, q);
        }
        const double oc = objective(I, cand, o.lamFleet);
        if (oc <= (ig ? curO * (1 + band) : curO) + 1e-9) { cur = cand; curO = oc; T.offer(cand, oc); }
        if (T.timeUp(o.timeLimit)) break;
    }
    return T.finish(std::min(it, o.iters));
}

// ---- giant-tour tools of the elite-pool methods
// broken-pairs distance: adjacent pairs of a that are not adjacent in b
int brokenPairs(const std::vector<int>& a, const std::vector<int>& b) {
    const int n = (int)b.size();
    std::vector<int> l(n, -1), r(n, -1);
    for (int i = 0; i < n; i++) { if (i) l[b[i]] = b[i - 1]; if (i + 1 < n) r[b[i]] = b[i + 1]; }
    int d = 0;
    for (int i = 0; i + 1 < (int)a.size(); i++) { const int x = a[i], y = a[i + 1]; if (l[x] != y && r[x] != y) d++; }
    return d;
}
// path relinking on permutations: move a towards b by position swaps; the tours at the given fractions of the path
std::vector<std::vector<int>> relink(const std::vector<int>& a, const std::vector<int>& b, const std::vector<double>& fracs) {
    const int n = (int)a.size();
    std::vector<int> cur = a, pos(n);
    for (int p = 0; p < n; p++) pos[cur[p]] = p;
    std::vector<std::pair<int, int>> sw;
    for (int i = 0; i < n; i++)
        if (cur[i] != b[i]) {
            const int j = pos[b[i]];
            std::swap(cur[i], cur[j]);
            pos[cur[i]] = i;
            pos[cur[j]] = j;
            sw.push_back({i, j});
        }
    std::vector<std::vector<int>> out;
    const int S = (int)sw.size();
    if (S < 2) return out;
    for (double f : fracs) {
        const int k = std::max(1, std::min(S - 1, (int)std::floor(f * S + 0.5)));
        std::vector<int> t = a;
        for (int s = 0; s < k; s++) std::swap(t[sw[s].first], t[sw[s].second]);
        out.push_back(std::move(t));
    }
    return out;
}
Ind makeInd(const Instance& I, const DepotOrder& dord, const Solution& s, double lam) {
    Ind x;
    x.sol = s;
    x.o = objective(I, s, lam);
    x.tour = giantTour(I, s);
    x.dep = depGenes(I, dord, s);
    return x;
}

// ---- aco: MAX-MIN ant system (Stuetzle & Hoos 2000) with arc pheromones for routing (Bullnheimer, Hartl & Strauss 1999).
// One iteration = one ant: probabilistic nearest-neighbour construction (weights tau^alpha (1/d)^beta over the unrouted feasible
// customers among the k nearest), then descent.  Every `ants` ants the pheromone evaporates and the colony best (every 5th cycle
// the global best) deposits 1/cost; tau is kept in [tmin, tmax].
Result runACO(const Instance& I, const MethodOptions& o, LSContext& ctx, const Solution& init) {
    ctx.ensureAdjacency();
    PyRandom rng(o.seed);
    Tracker T(I, init, o.lamFleet, o.target);
    const DepotOrder dord = buildDepotOrder(I);
    const double alpha = par(o, "alpha", 1), beta = par(o, "beta", 3), rho = par(o, "rho", 0.1);
    const int ants = ipar(o, "ants", 10), kn = 20, N = I.n + I.t;
    Solution globS = init.nonEmpty();
    double globO = objective(I, globS, o.lamFleet);
    double tmax = 1.0 / (rho * globO), tmin = tmax / (2.0 * I.n);
    std::vector<float> tau((size_t)N * N, (float)tmax);
    std::vector<int> grp(I.n);
    std::vector<std::vector<int>> groups(I.t);
    for (int c = 0; c < I.n; c++) { grp[c] = dord[c][0]; groups[grp[c]].push_back(c); }
    auto deposit = [&](const Solution& s, double amt) {
        for (size_t r = 0; r < s.routes.size(); r++) {
            if (s.routes[r].empty()) continue;
            const int e = I.n + s.depots[r];
            int a = e;
            for (int c : s.routes[r]) { tau[(size_t)a * N + c] += (float)amt; tau[(size_t)c * N + a] += (float)amt; a = c; }
            tau[(size_t)a * N + e] += (float)amt;
            tau[(size_t)e * N + a] += (float)amt;
        }
    };
    std::vector<char> un(I.n);
    std::vector<int> cand;
    std::vector<double> wgt;
    auto construct = [&]() {
        std::fill(un.begin(), un.end(), 1);
        Solution s;
        for (int d = 0; d < I.t; d++) {
            const int dn = I.n + d;
            const std::vector<int>& G = groups[d];
            int remaining = (int)G.size();
            while (remaining > 0) {
                std::vector<int> cur;
                int pos = dn;
                double q = 0, dur = 0;
                for (;;) {
                    cand.clear();
                    wgt.clear();
                    double tot = 0;
                    auto scan = [&](int c) {
                        if (!un[c] || grp[c] != d) return;
                        if (q + I.demand[c] > I.capacity[d] + EPS || dur + I.d(pos, c) + I.service[c] + I.d(c, dn) > I.maxDuration[d] + EPS) return;
                        const double x = std::pow((double)tau[(size_t)pos * N + c], alpha) * std::pow(1.0 / (I.d(pos, c) + 1e-9), beta);
                        cand.push_back(c);
                        wgt.push_back(x);
                        tot += x;
                    };
                    if (pos != dn) {                                              // the k nearest customers first (adj lists the customer itself first)
                        const std::vector<int>& a = ctx.adj[pos];
                        const int m = std::min(kn + 1, (int)a.size());
                        for (int i = 0; i < m; i++) if (a[i] != pos) scan(a[i]);
                    }
                    if (cand.empty()) for (int c : G) scan(c);
                    if (cand.empty()) break;
                    double x = rng.random() * tot;
                    int pi = (int)cand.size() - 1;
                    for (int i = 0; i < (int)cand.size(); i++) { x -= wgt[i]; if (x <= 0) { pi = i; break; } }
                    const int c = cand[pi];
                    un[c] = 0;
                    remaining--;
                    cur.push_back(c);
                    q += I.demand[c];
                    dur += I.d(pos, c) + I.service[c];
                    pos = c;
                }
                if (cur.empty()) {
                    for (int c : G) if (un[c]) { un[c] = 0; remaining--; cur.push_back(c); break; }
                }
                s.routes.push_back(std::move(cur));
                s.depots.push_back(d);
            }
        }
        return repairFleetLimit(I, s);
    };
    std::vector<std::pair<Solution, double>> colony;
    int cycle = 0, it;
    for (it = 1; it <= o.iters; it++) {
        T.it = it;
        Solution s = descend(ctx, construct(), o, rng);
        const double os = objective(I, s, o.lamFleet);
        T.offer(s, os);
        colony.push_back({std::move(s), os});
        if ((int)colony.size() >= ants) {
            size_t b = 0;
            for (size_t i = 1; i < colony.size(); i++) if (colony[i].second < colony[b].second) b = i;
            if (colony[b].second < globO) { globS = colony[b].first; globO = colony[b].second; }
            cycle++;
            for (float& x : tau) x *= (float)(1 - rho);
            if (cycle % 5 == 0) deposit(globS, 1.0 / globO); else deposit(colony[b].first, 1.0 / colony[b].second);
            tmax = 1.0 / (rho * globO);
            tmin = tmax / (2.0 * I.n);
            for (float& x : tau) { if (x > tmax) x = (float)tmax; else if (x < tmin) x = (float)tmin; }
            colony.clear();
        }
        if (T.timeUp(o.timeLimit)) break;
    }
    return T.finish(std::min(it, o.iters));
}

// ---- grasppr: GRASP with path relinking and an elite pool (Resende & Ribeiro 2005; Prins, Prodhon & Wolfler Calvo 2006).
// One iteration = one randomised construction + descent, then a relinking of the local optimum towards a random elite solution on
// giant tours; the tours at 25 / 50 / 75 % of the path are split optimally, the best is improved by descent and offered to the pool.
Result runGraspPR(const Instance& I, const MethodOptions& o, const LSContext& ctx, const Solution& init) {
    PyRandom rng(o.seed);
    Tracker T(I, init, o.lamFleet, o.target);
    const DepotOrder dord = buildDepotOrder(I);
    const int eliteSize = 8;
    std::vector<Ind> elite = {makeInd(I, dord, init.nonEmpty(), o.lamFleet)};
    auto add = [&](Ind ind) {
        for (const Ind& x : elite) if (std::fabs(x.o - ind.o) < 1e-7) return;
        if ((int)elite.size() < eliteSize) { elite.push_back(std::move(ind)); return; }
        int w = 0;
        for (int i = 1; i < (int)elite.size(); i++) if (elite[i].o > elite[w].o) w = i;
        if (ind.o < elite[w].o - EPS) elite[w] = std::move(ind);
    };
    int it;
    for (it = 1; it <= o.iters; it++) {
        T.it = it;
        Solution x = descend(ctx, graspConstruct(I, dord, rng, 0.2), o, rng);
        T.offer(x);
        Ind X = makeInd(I, dord, x, o.lamFleet);
        const std::vector<int> ytour = elite[rng.randrange((int)elite.size())].tour;
        Solution bs;
        double bo = INF;
        for (const auto& t : relink(X.tour, ytour, {0.25, 0.5, 0.75})) {
            Solution s = decodeTour(I, dord, t);
            const double os = objective(I, s, o.lamFleet);
            if (os < bo) { bo = os; bs = std::move(s); }
        }
        add(std::move(X));
        if (bo < INF) {
            Solution s = descend(ctx, bs, o, rng);
            T.offer(s);
            add(makeInd(I, dord, s, o.lamFleet));
        }
        if (T.timeUp(o.timeLimit)) break;
    }
    return T.finish(std::min(it, o.iters));
}

// ---- scatter: scatter search (Glover 1998; Marti, Laguna & Glover 2006).  Diversified pool of random-decode local optima ->
// reference set of b1 best + b2 most distant (broken-pairs distance) -> every pair is combined by path relinking (midpoint of the
// path in a random direction) -> descent -> the reference set is rebuilt from the old set plus the children; when it no longer
// changes, the b2 diverse members are regenerated.  One iteration = one combination, one new pool member or one set update.
Result runScatter(const Instance& I, const MethodOptions& o, const LSContext& ctx, const Solution& init) {
    PyRandom rng(o.seed);
    Tracker T(I, init, o.lamFleet, o.target);
    const DepotOrder dord = buildDepotOrder(I);
    const int b1 = 6, b2 = 4, poolSize = 16;
    struct SInd { Ind ind; bool isNew = false; };
    using SP = std::shared_ptr<SInd>;
    auto mk = [&](const Solution& s) { auto p = std::make_shared<SInd>(); p->ind = makeInd(I, dord, s, o.lamFleet); return p; };
    auto contains = [](const std::vector<SP>& v, const SP& x) { return std::find(v.begin(), v.end(), x) != v.end(); };
    auto select = [&](const std::vector<SP>& set) {
        std::vector<SP> sorted = set, out, rest;
        std::stable_sort(sorted.begin(), sorted.end(), [](const SP& a, const SP& b) { return a->ind.o < b->ind.o; });
        for (const SP& x : sorted) {
            if ((int)out.size() >= b1) break;
            bool dup = false;
            for (const SP& y : out) if (std::fabs(y->ind.o - x->ind.o) < 1e-7) { dup = true; break; }
            if (!dup) out.push_back(x);
        }
        for (const SP& x : sorted) if (!contains(out, x)) rest.push_back(x);
        while ((int)out.size() < b1 + b2 && !rest.empty()) {
            int bi = 0, bd = -1;
            for (int i = 0; i < (int)rest.size(); i++) {
                int m = std::numeric_limits<int>::max();
                for (const SP& y : out) m = std::min(m, brokenPairs(rest[i]->ind.tour, y->ind.tour));
                if (m > bd) { bd = m; bi = i; }
            }
            out.push_back(rest[bi]);
            rest.erase(rest.begin() + bi);
        }
        return out;
    };
    std::vector<SP> pool = {mk(init.nonEmpty())}, ref, cands;
    std::vector<std::pair<SP, SP>> queue;
    auto makeQueue = [&]() {
        queue.clear();
        for (size_t i = 0; i < ref.size(); i++)
            for (size_t j = i + 1; j < ref.size(); j++)
                if (ref[i]->isNew || ref[j]->isNew) queue.push_back({ref[i], ref[j]});
        for (int i = (int)queue.size() - 1; i > 0; i--) std::swap(queue[i], queue[rng.randrange(i + 1)]);
    };
    enum Phase { INIT, COMB, REGEN } phase = INIT;
    int regen = 0, it;
    for (it = 1; it <= o.iters; it++) {
        T.it = it;
        if (phase == INIT) {
            Solution s = randomStart(I, ctx, dord, o, rng);
            T.offer(s);
            pool.push_back(mk(s));
            if ((int)pool.size() >= poolSize) {
                ref = select(pool);
                for (SP& x : ref) x->isNew = true;
                makeQueue();
                phase = COMB;
            }
        } else if (phase == REGEN) {
            Solution s = randomStart(I, ctx, dord, o, rng);
            T.offer(s);
            SP ind = mk(s);
            ind->isNew = true;
            ref.push_back(ind);
            if (++regen >= b2) { ref = select(ref); makeQueue(); phase = COMB; }
        } else if (queue.empty()) {                                                // reference-set update
            const std::vector<SP> old = ref;
            std::vector<SP> merged = ref;
            merged.insert(merged.end(), cands.begin(), cands.end());
            cands.clear();
            std::vector<SP> next = select(merged);
            bool changed = false;
            for (SP& x : next) { if (!contains(old, x)) { x->isNew = true; changed = true; } else x->isNew = false; }
            ref = next;
            if (changed) makeQueue();
            else {
                if ((int)ref.size() > b1) ref.resize(b1);
                for (SP& x : ref) x->isNew = false;
                regen = 0;
                phase = REGEN;
            }
            if (queue.empty() && phase == COMB) { for (SP& x : ref) x->isNew = true; makeQueue(); }
        } else {
            const std::pair<SP, SP> pr = queue.back();
            queue.pop_back();
            const bool fwd = rng.random() < 0.5;
            const Ind& A = fwd ? pr.first->ind : pr.second->ind;
            const Ind& B = fwd ? pr.second->ind : pr.first->ind;
            const auto mid = relink(A.tour, B.tour, {0.5});
            const std::vector<int> tour = !mid.empty() ? mid[0] : orderCrossover(A.tour, B.tour, rng);
            Solution s = descend(ctx, decodeTour(I, dord, tour), o, rng);
            T.offer(s);
            cands.push_back(mk(s));
        }
        if (T.timeUp(o.timeLimit)) break;
    }
    return T.finish(std::min(it, o.iters));
}

// ---- random-key methods.  All share decodeKeys; the fitness is the decoded objective, the optional descent (probability pls)
// only improves the reported solution (the keys are not changed).
struct KeyEval {
    const Instance& I;
    const LSContext& ctx;
    const DepotOrder& dord;
    const MethodOptions& o;
    PyRandom& rng;
    Tracker& T;
    double pls;
    double operator()(const Keys& k) {
        Solution s = decodeKeys(I, dord, k);
        const double os = objective(I, s, o.lamFleet);
        T.offer(s, os);
        if (pls > 0 && rng.random() < pls) T.offer(descend(ctx, s, o, rng));
        return os;
    }
};

// brkga (Goncalves & Resende 2011): population 30, elite 20 %, mutants 15 %, a child gene comes from the elite parent with
// probability 0.7; one iteration = one decoded individual (the elite is copied into the next generation for free)
// pso (Kennedy & Eberhart 1995, random-key coding): inertia 0.9 -> 0.4, c1 = c2 = 1.5, velocity clamp 0.5, positions modulo 1
// de  (Storn & Price 1997): DE/rand/1/bin, F = 0.5, CR = 0.9, mutant keys modulo 1; one iteration = one trial vector
// rkrs: a fresh uniform key vector every iteration; the control baseline of the family (no evolution, no descent)
Result runKeys(const Instance& I, const MethodOptions& o, const LSContext& ctx, const Solution& init, const std::string& id) {
    PyRandom rng(o.seed);
    Tracker T(I, init, o.lamFleet, o.target);
    const DepotOrder dord = buildDepotOrder(I);
    const int dim = keyLen(I);
    KeyEval ev{I, ctx, dord, o, rng, T, id == "rkrs" ? 0.0 : 0.05};
    struct KI { Keys keys; double o; };
    auto byO = [](const KI& a, const KI& b) { return a.o < b.o; };
    int it = 1;
    if (id == "rkrs") {
        for (; it <= o.iters; it++) { T.it = it; ev(randomKeys(I, rng)); if (T.timeUp(o.timeLimit)) break; }
    } else if (id == "brkga") {
        const int pop = 30, ne = std::max(1, (int)std::floor(0.2 * pop + 0.5)), nm = (int)std::floor(0.15 * pop + 0.5);
        const double rho = 0.7;
        std::vector<KI> cur, nxt;
        bool evolving = false;
        for (; it <= o.iters; it++) {
            T.it = it;
            if (!evolving) {
                Keys k = cur.empty() ? seedKeys(I, dord, init) : randomKeys(I, rng);
                const double f = ev(k);
                cur.push_back({std::move(k), f});
                if ((int)cur.size() >= pop) {
                    std::stable_sort(cur.begin(), cur.end(), byO);
                    nxt.assign(cur.begin(), cur.begin() + ne);
                    evolving = true;
                }
            } else {
                const int idx = (int)nxt.size() - ne, nOff = pop - ne;
                Keys k;
                if (idx >= nOff - nm) k = randomKeys(I, rng);
                else {
                    const KI& e = cur[rng.randrange(ne)];
                    const KI& x = cur[ne + rng.randrange(pop - ne)];
                    k.resize(dim);
                    for (int i = 0; i < dim; i++) k[i] = rng.random() < rho ? e.keys[i] : x.keys[i];
                }
                const double f = ev(k);
                nxt.push_back({std::move(k), f});
                if ((int)nxt.size() >= pop) {
                    cur = std::move(nxt);
                    std::stable_sort(cur.begin(), cur.end(), byO);
                    nxt.assign(cur.begin(), cur.begin() + ne);
                }
            }
            if (T.timeUp(o.timeLimit)) break;
        }
    } else if (id == "pso") {
        const int swarm = ipar(o, "swarm", 20);
        const double w0 = par(o, "w0", 0.9), w1 = par(o, "w1", 0.4), c1 = par(o, "c1", 1.5), c2 = par(o, "c2", 1.5), vmax = par(o, "vmax", 0.5);
        struct Pt { Keys x, v, pb; double po; };
        std::vector<Pt> sw;
        Keys gb;
        double gbo = INF;
        for (; it <= o.iters; it++) {
            T.it = it;
            if ((int)sw.size() < swarm) {
                Pt p;
                p.x = sw.empty() ? seedKeys(I, dord, init) : randomKeys(I, rng);
                p.v.resize(dim);
                for (int i = 0; i < dim; i++) p.v[i] = (rng.random() - 0.5) * vmax;
                p.po = ev(p.x);
                p.pb = p.x;
                if (p.po < gbo) { gbo = p.po; gb = p.x; }
                sw.push_back(std::move(p));
            } else {
                Pt& p = sw[(it - 1) % swarm];
                const double w = w0 + (w1 - w0) * std::min(1.0, (double)it / o.iters);
                for (int i = 0; i < dim; i++) {
                    double v = w * p.v[i];
                    v += c1 * rng.random() * (p.pb[i] - p.x[i]);
                    v += c2 * rng.random() * (gb[i] - p.x[i]);
                    v = std::max(-vmax, std::min(vmax, v));
                    p.v[i] = v;
                    const double x = p.x[i] + v;
                    p.x[i] = x - std::floor(x);
                }
                const double f = ev(p.x);
                if (f < p.po) { p.po = f; p.pb = p.x; }
                if (f < gbo) { gbo = f; gb = p.x; }
            }
            if (T.timeUp(o.timeLimit)) break;
        }
    } else {                                                                       // de
        const int np = 20;
        const double F = 0.5, CR = 0.9;
        std::vector<KI> pop;
        int tgt = 0;
        for (; it <= o.iters; it++) {
            T.it = it;
            if ((int)pop.size() < np) {
                Keys x = pop.empty() ? seedKeys(I, dord, init) : randomKeys(I, rng);
                const double f = ev(x);
                pop.push_back({std::move(x), f});
            } else {
                const int i = tgt;
                tgt = (tgt + 1) % np;
                int a, b, c;
                do { a = rng.randrange(np); } while (a == i);
                do { b = rng.randrange(np); } while (b == i || b == a);
                do { c = rng.randrange(np); } while (c == i || c == a || c == b);
                const int jr = rng.randrange(dim);
                Keys t(dim);
                for (int j = 0; j < dim; j++) {
                    if (rng.random() < CR || j == jr) {
                        const double v = pop[a].keys[j] + F * (pop[b].keys[j] - pop[c].keys[j]);
                        t[j] = v - std::floor(v);
                    } else t[j] = pop[i].keys[j];
                }
                const double f = ev(t);
                if (f <= pop[i].o) pop[i] = {std::move(t), f};
            }
            if (T.timeUp(o.timeLimit)) break;
        }
    }
    return T.finish(std::min(it, o.iters));
}
#include "methods_full.inc"                      // the upgraded ("full") variants: tabu-f, ...
}  // namespace

// ---------------------------------------------------------------------------------------------
const std::vector<ParamSpec>& methodParams(const std::string& id) {
    static const std::map<std::string, std::vector<ParamSpec>> table = {
        {"sisr", {{"cbar", 'r', 10, 5, 20, false, "average number of removed customers"}, {"lmax", 'i', 10, 5, 20, false, "maximum string length"},
                  {"alpha", 'r', 0.01, 0.001, 0.5, true, "split-string preservation probability"}, {"blink", 'r', 0.01, 0, 0.1, false, "blink rate of the insertion"},
                  {"t0", 'r', 0.45, 0.05, 2, true, "start temperature / mean arc"}, {"tf", 'r', 0.0045, 0.0005, 0.05, true, "final temperature / mean arc"}}},
        {"ils-f", {{"pmin", 'i', 2, 1, 4, false, "minimum number of perturbation moves"}, {"prange", 'i', 3, 1, 6, false, "number of perturbation moves = pmin + U{0..prange-1}"},
                   {"failfactor", 'r', 5, 1, 20, false, "restart after n + failfactor x routes failures"}, {"alpha", 'r', 0.3, 0.05, 0.6, false, "restricted candidate list of the restart construction"}}},
        {"sa-f", {{"resetfactor", 'r', 20, 2, 100, true, "temperature reset after resetfactor x n rejections"}, {"tfloor", 'r', 0.002, 0.0002, 0.02, true, "lower bound of the final temperature / mean arc"}}},
        {"ta", {{"th0", 'r', 0.05, 0.005, 0.3, true, "start threshold / mean arc"}}},
        {"rtr", {{"dev0", 'r', 0.02, 0.002, 0.1, true, "start band relative to the record"}}},
        {"tabu-f", {{"beta0", 'r', 1.5, 1, 2.5, false, "sparsification factor"}, {"betadiv", 'r', 3.0, 2, 5, false, "sparsification factor during diversification"},
                    {"tmin", 'i', 5, 2, 15, false, "minimum tabu tenure"}, {"tspan", 'i', 5, 0, 15, false, "tenure = tmin + U{0..tspan}"},
                    {"h", 'i', 10, 5, 30, false, "penalty update period"}, {"kmax", 'i', 40, 10, 60, false, "maximum neighbours per customer"},
                    {"divfactor", 'r', 2, 0.5, 10, true, "diversification after divfactor x n iterations without improvement"}}},
        {"vns-f", {{"kmax", 'i', 6, 2, 10, false, "largest shaking neighbourhood"}, {"t0", 'r', 0.5, 0, 2, false, "start temperature / mean arc"}}},
        {"grasp", {{"alpha", 'r', 0.2, 0.05, 0.6, false, "restricted candidate list (fraction)"}}},
        {"msls", {}},
        {"alns-f", {{"w", 'r', 0.05, 0.005, 0.3, true, "start temperature control (worsening accepted with probability 0.5)"}, {"cool", 'r', 0.00193, 0.0001, 0.05, true, "final / start temperature"},
                    {"react", 'r', 0.1, 0.01, 0.5, false, "reaction factor"}, {"eta", 'r', 0.025, 0, 0.1, false, "noise amplitude"}, {"xi", 'r', 0.4, 0.05, 0.5, false, "largest removed fraction"},
                    {"qmax", 'i', 100, 20, 150, false, "largest number of removed customers"}, {"s1", 'r', 33, 1, 50, false, "score: new best"}, {"s2", 'r', 9, 1, 50, false, "score: better than current"},
                    {"s3", 'r', 13, 1, 50, false, "score: accepted worse"}, {"seg", 'i', 100, 50, 300, false, "segment length"}}},
        {"lns", {{"qlo", 'r', 0.05, 0.01, 0.2, false, "smallest removed fraction"}, {"qhi", 'r', 0.2, 0.05, 0.5, false, "largest removed fraction"}}},
        {"ig-f", {{"frac", 'r', 0.1, 0.02, 0.3, false, "removed fraction"}, {"temp", 'r', 0.05, 0.001, 0.5, true, "constant temperature / mean arc"}}},
        {"ga-f", {{"mu", 'i', 25, 10, 60, false, "minimum sub-population size"}, {"lambda", 'i', 40, 10, 100, false, "generation size"}, {"nelite", 'i', 4, 1, 10, false, "elite individuals"},
                  {"nclose", 'i', 5, 2, 10, false, "closest individuals in the diversity measure"}, {"target", 'r', 0.2, 0.05, 0.5, false, "target share of feasible offspring"},
                  {"repair", 'r', 0.5, 0, 1, false, "repair probability"}, {"divafter", 'i', 4000, 500, 20000, true, "diversification after this many offspring without improvement"}}},
        {"aco", {{"alpha", 'r', 1, 0.5, 5, false, "pheromone exponent"}, {"beta", 'r', 3, 1, 8, false, "visibility exponent"}, {"rho", 'r', 0.1, 0.01, 0.5, true, "evaporation rate"}, {"ants", 'i', 10, 5, 30, false, "ants per cycle"}}},
        {"grasppr-f", {{"alpha", 'r', 0.2, 0.05, 0.5, false, "skip probability of the randomised savings"}, {"elite", 'i', 10, 4, 20, false, "elite pool size"}, {"pevo", 'r', 0.1, 0, 0.5, false, "probability of an elite-elite relinking"}}},
        {"scatter-f", {{"b", 'i', 10, 6, 20, false, "reference set size"}, {"bdiv", 'i', 4, 2, 8, false, "members replaced at regeneration"}}},
        {"brkga-f", {{"pop", 'i', 30, 15, 80, false, "population size"}, {"pe", 'r', 0.2, 0.1, 0.4, false, "elite fraction"}, {"pm", 'r', 0.15, 0.05, 0.3, false, "mutant fraction"}, {"restart", 'i', 30, 10, 100, false, "restart after this many generations without improvement"}}},
        {"pso", {{"w0", 'r', 0.9, 0.5, 1.2, false, "start inertia"}, {"w1", 'r', 0.4, 0.1, 0.8, false, "final inertia"}, {"c1", 'r', 1.5, 0.5, 2.5, false, "cognitive coefficient"},
                 {"c2", 'r', 1.5, 0.5, 2.5, false, "social coefficient"}, {"vmax", 'r', 0.5, 0.1, 1, false, "velocity clamp"}, {"swarm", 'i', 20, 10, 50, false, "swarm size"}}},
        {"de-f", {{"np", 'i', 20, 10, 50, false, "population size"}, {"cr", 'r', 0.9, 0.5, 1, false, "crossover rate"}, {"fmin", 'r', 0.4, 0.1, 0.8, false, "smallest scale factor"}, {"fspan", 'r', 0.5, 0, 0.8, false, "scale factor = fmin + U(0, fspan)"}}},
        {"rkrs", {}},
    };
    static const std::vector<ParamSpec> none;
    auto f = table.find(id);
    return f == table.end() ? none : f->second;
}

const std::vector<std::string>& fullMethodIds() {
    static const std::vector<std::string> ids = {"ils-f", "sa-f", "ta-f", "rtr-f", "tabu-f", "vns-f", "grasp-f", "alns-f", "lns-f", "ig-f", "ga-f", "aco-f", "grasppr-f", "scatter-f", "brkga-f", "pso-f", "de-f"};
    return ids;
}

const std::vector<std::string>& methodIds() {
    static const std::vector<std::string> ids = {"sa", "ta", "rtr", "tabu", "vns", "grasp", "msls", "alns", "lns", "ig",
                                                 "ga", "aco", "grasppr", "scatter", "brkga", "pso", "de", "rkrs"};
    return ids;
}

int defaultIters(const std::string& id) {
    static const std::unordered_map<std::string, int> d = {
        {"sa", 2000000}, {"ta", 2000000}, {"rtr", 2000000}, {"tabu", 50000}, {"vns", 300}, {"grasp", 200}, {"msls", 200},
        {"alns", 10000}, {"lns", 10000}, {"ig", 10000}, {"ga", 400}, {"aco", 300}, {"grasppr", 150}, {"scatter", 300},
        {"brkga", 20000}, {"pso", 20000}, {"de", 20000}, {"rkrs", 20000},
        {"tabu-f", 20000}, {"alns-f", 25000}, {"ga-f", 2000}, {"ils-f", 500}, {"sa-f", 2000000}, {"ta-f", 2000000}, {"rtr-f", 2000000}, {"vns-f", 300}, {"grasp-f", 200},
        {"lns-f", 10000}, {"ig-f", 2000}, {"aco-f", 300}, {"grasppr-f", 150}, {"scatter-f", 300}, {"brkga-f", 600}, {"pso-f", 600}, {"de-f", 600}};
    auto f = d.find(id);
    if (f == d.end()) throw std::runtime_error("unknown method " + id);
    return f->second;
}

static Result runMethodImpl(const std::string& id, const Instance& inst, const MethodOptions& opt);
Result runMethod(const std::string& id, const Instance& inst, const MethodOptions& opt) {
    std::array<long, OP_COUNT> st{};
    std::array<double, OP_COUNT> gn{};
    struct Guard { ~Guard() { g_opStats = nullptr; g_opGain = nullptr; } } guard;
    for (const auto& kv : opt.params) {                              // unknown parameter names are errors, not silently ignored
        bool known = false;
        for (const ParamSpec& p : methodParams(id)) known = known || p.name == kv.first;
        if (!known) throw std::runtime_error("method " + id + " has no parameter " + kv.first);
    }
    g_opStats = &st;
    g_opGain = &gn;
    Result r = runMethodImpl(id, inst, opt);
    r.opStats = st;
    r.opGain = gn;
    return r;
}

static Result runMethodImpl(const std::string& id, const Instance& inst, const MethodOptions& opt) {
    MethodOptions o = opt;
    if (o.iters < 0) o.iters = defaultIters(id);
    if (o.iters < 1) throw std::runtime_error("the number of iterations must be positive");
    std::unique_ptr<LSContext> own;
    LSContext* ctx = o.ctx;
    if (!ctx) { own = std::make_unique<LSContext>(inst, o.k); ctx = own.get(); }
    const Solution init = o.init ? o.init->nonEmpty() : savingsMd(inst);
    if (id == "tabu-f") return runTabuFull(inst, o, *ctx, init);
    if (id == "alns-f") return runAlnsFull(inst, o, *ctx, init);
    if (id == "ga-f") return runGaFull(inst, o, *ctx, init);
    if (id == "ils-f") return runIlsFull(inst, o, *ctx, init);
    if (id == "sa-f") return runSaFull(inst, o, *ctx, init);
    if (id == "ta-f") return runTaFull(inst, o, *ctx, init);
    if (id == "rtr-f") return runRtrFull(inst, o, *ctx, init);
    if (id == "vns-f") return runVnsFull(inst, o, *ctx, init);
    if (id == "grasp-f") return runGraspFull(inst, o, *ctx, init);
    if (id == "lns-f") return runLnsFull(inst, o, *ctx, init);
    if (id == "ig-f") return runIgFull(inst, o, *ctx, init);
    if (id == "aco-f") return runAcoFull(inst, o, *ctx, init);
    if (id == "grasppr-f") return runGraspPrFull(inst, o, *ctx, init);
    if (id == "scatter-f") return runScatterFull(inst, o, *ctx, init);
    if (id == "brkga-f" || id == "pso-f" || id == "de-f") return runKeysFull(inst, o, *ctx, init, id);
    if (id == "sa") return runSA(inst, o, *ctx, init);
    if (id == "ta") return runWalkRule(inst, o, *ctx, init, false);
    if (id == "rtr") return runWalkRule(inst, o, *ctx, init, true);
    if (id == "tabu") return runTabu(inst, o, *ctx, init);
    if (id == "vns") return runVNS(inst, o, *ctx, init);
    if (id == "grasp") return runMultiStart(inst, o, *ctx, init, true);
    if (id == "msls") return runMultiStart(inst, o, *ctx, init, false);
    if (id == "alns") return runALNS(inst, o, *ctx, init);
    if (id == "lns") return runLnsIg(inst, o, *ctx, init, false);
    if (id == "ig") return runLnsIg(inst, o, *ctx, init, true);
    if (id == "ga") return runGA(inst, o, *ctx, init);
    if (id == "aco") return runACO(inst, o, *ctx, init);
    if (id == "grasppr") return runGraspPR(inst, o, *ctx, init);
    if (id == "scatter") return runScatter(inst, o, *ctx, init);
    if (id == "brkga" || id == "pso" || id == "de" || id == "rkrs") return runKeys(inst, o, *ctx, init, id);
    throw std::runtime_error("unknown method " + id);
}

}  // namespace vrp
