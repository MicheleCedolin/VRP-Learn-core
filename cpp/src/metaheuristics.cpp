#include "vrplearn/metaheuristics.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <memory>

#include "vrplearn/construct.hpp"

namespace vrp {

namespace {
constexpr double EPS = 1e-9;

double now() {
    using namespace std::chrono;
    return duration<double>(steady_clock::now().time_since_epoch()).count();
}
}  // namespace

// ---------------------------------------------------------------------------------------------
RRState::RRState(const LSContext& ctx, const Solution& sol) : ctx_(&ctx), n_(ctx.inst.n), t_(ctx.inst.t) {
    for (size_t r = 0; r < sol.routes.size(); r++)
        if (!sol.routes[r].empty()) { routes.push_back(sol.routes[r]); depots.push_back(sol.depots[r]); }
    load.assign(routes.size(), 0.0);
    costR.assign(routes.size(), 0.0);
    dur.assign(routes.size(), 0.0);
    rid.assign(n_, -1);
    for (size_t r = 0; r < routes.size(); r++) recompute(r);
    cost = pysum(costR);
    cnt.assign(t_, 0);
    for (size_t r = 0; r < routes.size(); r++)
        if (!routes[r].empty()) cnt[depots[r]]++;
}

void RRState::recompute(size_t r) {
    const Instance& I = ctx_->inst;
    const std::vector<int>& nodes = routes[r];
    const int dn = n_ + depots[r];
    if (!nodes.empty()) {
        double c = I.d(dn, nodes.front()) + I.d(nodes.back(), dn);
        for (size_t k = 0; k + 1 < nodes.size(); k++) c += I.d(nodes[k], nodes[k + 1]);
        costR[r] = c;
        std::vector<double> q, s;
        for (int v : nodes) { q.push_back(I.demand[v]); s.push_back(I.service[v]); }
        load[r] = pysum(q);
        dur[r] = c + pysum(s);
    } else {
        costR[r] = load[r] = dur[r] = 0.0;
    }
    for (int v : nodes) rid[v] = (int)r;
}

int RRState::fleetExcess() const {
    int e = 0;
    for (int d = 0; d < t_; d++) e += std::max(0, cnt[d] - ctx_->inst.fleet[d]);
    return e;
}

Solution RRState::toSolution() const {
    Solution s;
    for (size_t r = 0; r < routes.size(); r++)
        if (!routes[r].empty()) { s.routes.push_back(routes[r]); s.depots.push_back(depots[r]); }
    return s;
}

std::vector<int> RRState::ruin(PyRandom& rng, double cbar, int lmax, double alpha) {
    const LSContext& ctx = *ctx_;
    long liveRoutes = 0, totalLen = 0;
    for (const auto& r : routes) if (!r.empty()) { liveRoutes++; totalLen += (long)r.size(); }
    double avgCard = (double)totalLen / std::max(1L, liveRoutes);
    double lsMax = std::min((double)lmax, avgCard);
    double ksMax = 4.0 * cbar / (1.0 + lsMax) - 1.0;
    int ks = static_cast<int>(rng.uniform(1.0, ksMax + 1.0));
    int seed = rng.randrange(n_);
    std::vector<int> removed;
    std::vector<char> gone(n_, 0);
    std::vector<int> ruinedList;
    std::vector<char> ruined(routes.size(), 0);
    size_t nRuined = 0;
    for (int c : ctx.adj[seed]) {
        if ((int)nRuined >= ks) break;
        if (gone[c]) continue;
        int r = rid[c];
        if (r < 0 || ruined[r]) continue;
        const std::vector<int> route = routes[r];                   // copy: routes[r] is replaced below
        const int L = (int)route.size();
        double ltMax = std::min((double)L, lsMax);
        int lt = std::max(1, static_cast<int>(rng.uniform(1.0, ltMax + 1.0)));
        int pos = (int)(std::find(route.begin(), route.end(), c) - route.begin());
        std::vector<int> take;
        if (lt == L || rng.random() < 0.5) {                        // plain string
            int lo = std::max(0, pos - lt + 1), hi = std::min(pos, L - lt);
            int s = rng.randint(lo, hi);
            take.assign(route.begin() + s, route.begin() + s + lt);
            std::vector<int> nr(route.begin(), route.begin() + s);
            nr.insert(nr.end(), route.begin() + s + lt, route.end());
            routes[r] = nr;
        } else {                                                    // split string
            int m = 1;
            while (m < L - lt && rng.random() > alpha) m++;
            int w = lt + m;
            int lo = std::max(0, pos - w + 1), hi = std::min(pos, L - w);
            int s = rng.randint(lo, hi);
            std::vector<int> window(route.begin() + s, route.begin() + s + w);
            int off = rng.randint(0, lt);
            take.assign(window.begin(), window.begin() + off);
            take.insert(take.end(), window.begin() + off + m, window.end());
            std::vector<int> nr(route.begin(), route.begin() + s);
            nr.insert(nr.end(), window.begin() + off, window.begin() + off + m);
            nr.insert(nr.end(), route.begin() + s + w, route.end());
            routes[r] = nr;
        }
        ruined[r] = 1;
        ruinedList.push_back(r);
        nRuined++;
        for (int v : take) { removed.push_back(v); gone[v] = 1; }
    }
    for (int r : ruinedList) recompute(r);
    bool anyEmpty = false;
    for (const auto& rt : routes) if (rt.empty()) anyEmpty = true;
    if (anyEmpty) {                                                  // compact: emptied routes must not accumulate
        std::vector<size_t> keep;
        for (size_t r = 0; r < routes.size(); r++) if (!routes[r].empty()) keep.push_back(r);
        std::vector<std::vector<int>> nroutes;
        std::vector<int> ndepots;
        std::vector<double> nload, ncost, ndur;
        for (size_t r : keep) {
            nroutes.push_back(routes[r]); ndepots.push_back(depots[r]);
            nload.push_back(load[r]); ncost.push_back(costR[r]); ndur.push_back(dur[r]);
        }
        routes.swap(nroutes); depots.swap(ndepots); load.swap(nload); costR.swap(ncost); dur.swap(ndur);
        for (size_t r = 0; r < routes.size(); r++) for (int v : routes[r]) rid[v] = (int)r;
    }
    cost = pysum(costR);
    cnt.assign(t_, 0);
    for (int d : depots) cnt[d]++;
    return removed;
}

void RRState::recreate(std::vector<int>& removed, PyRandom& rng, double blink) {
    const Instance& I = ctx_->inst;
    const int n = n_;
    int mode = rng.choice({4, 4, 2, 1});                             // random, demand, far, close
    auto nd = [&](int c) {
        double m = INF;
        for (int d = 0; d < t_; d++) m = std::min(m, I.d(c, n + d));
        return m;
    };
    if (mode == 0) rng.shuffle(removed);
    else if (mode == 1) std::stable_sort(removed.begin(), removed.end(), [&](int a, int b) { return -I.demand[a] < -I.demand[b]; });
    else if (mode == 2) std::stable_sort(removed.begin(), removed.end(), [&](int a, int b) { return -nd(a) < -nd(b); });
    else std::stable_sort(removed.begin(), removed.end(), [&](int a, int b) { return nd(a) < nd(b); });

    for (int c : removed) {
        int bestR = -1, bestP = -1;
        double bd = INF;
        for (size_t r = 0; r < routes.size(); r++) {
            const std::vector<int>& nodes = routes[r];
            if (nodes.empty()) continue;
            int d = depots[r];
            if (load[r] + I.demand[c] > I.capacity[d] + EPS) continue;
            const int dn = n + d;
            int prev = dn;
            const int L = (int)nodes.size();
            for (int p = 0; p <= L; p++) {
                int nxt = p < L ? nodes[p] : dn;
                if (rng.random() >= blink) {
                    double delta = I.d(prev, c) + I.d(c, nxt) - I.d(prev, nxt);
                    if (delta < bd && dur[r] + delta + I.service[c] <= I.maxDuration[d] + EPS) { bestR = (int)r; bestP = p; bd = delta; }
                }
                prev = nxt;
            }
        }
        if (bestR >= 0) {
            routes[bestR].insert(routes[bestR].begin() + bestP, c);
            costR[bestR] += bd;
            load[bestR] += I.demand[c];
            dur[bestR] += bd + I.service[c];
            rid[c] = bestR;
            cost += bd;
        } else {                                                     // open a new route
            std::vector<int> cand(t_);
            for (int d = 0; d < t_; d++) cand[d] = d;
            std::stable_sort(cand.begin(), cand.end(), [&](int a, int b) { return 2 * I.d(c, n + a) < 2 * I.d(c, n + b); });
            int d = cand[0];
            for (int x : cand) if (cnt[x] < I.fleet[x]) { d = x; break; }
            routes.push_back({c});
            depots.push_back(d);
            double oc = 2 * I.d(c, n + d);
            costR.push_back(oc);
            load.push_back(I.demand[c]);
            dur.push_back(oc + I.service[c]);
            rid[c] = (int)routes.size() - 1;
            cnt[d]++;
            cost += oc;
        }
    }
}

// ---------------------------------------------------------------------------------------------
namespace {
double avgEdge(const Instance& inst, const Solution& sol) {
    return sol.cost(inst) / std::max(1, inst.n + sol.nRoutes());
}
}  // namespace

Result sisr(const Instance& inst, const SisrOptions& o) {
    std::unique_ptr<LSContext> own;
    LSContext* ctx = o.ctx;
    if (!ctx) { own = std::make_unique<LSContext>(inst, 20); ctx = own.get(); }
    ctx->ensureAdjacency();
    PyRandom rng(o.seed);
    Solution initSol = o.init ? *o.init : savingsMd(inst);
    RRState cur(*ctx, initSol);
    RRState best = cur;
    const double avg = avgEdge(inst, initSol);
    const double T0 = o.t0Factor * avg, Tf = o.tfFactor * avg;
    const double cool = std::pow(Tf / T0, 1.0 / std::max(1, o.iters));
    double T = T0;
    const double tStart = now();
    Result res;
    res.history.push_back({0, 0.0, best.cost});
    auto pen = [&](const RRState& s) { return o.lamFleet * s.fleetExcess(); };
    double curObj = cur.cost + pen(cur), bestObj = best.cost + pen(best);
    if (o.target > 0 && pen(best) == 0 && best.cost <= o.target * (1 + 1e-9)) res.reachedTarget = true;   // the start is already there
    int it = 0;
    for (it = 1; it <= o.iters && !res.reachedTarget; it++) {
        RRState cand = cur;
        std::vector<int> removed = cand.ruin(rng, o.cbar, o.lmax, o.alpha);
        cand.recreate(removed, rng, o.blink);
        double obj = cand.cost + pen(cand);
        if (obj < curObj - T * std::log(1.0 - rng.random())) {
            cur = std::move(cand);
            curObj = obj;
            if (obj < bestObj - EPS) {
                best = cur;
                bestObj = obj;
                res.history.push_back({it, now() - tStart, best.cost});
                if (o.target > 0 && pen(best) == 0 && best.cost <= o.target * (1 + 1e-9)) { res.reachedTarget = true; break; }
            }
        }
        T *= cool;
        if (o.timeLimit > 0 && (it & 63) == 0 && now() - tStart > o.timeLimit) break;
    }
    if (it > o.iters) it = o.iters;
    res.best = best.toSolution();
    res.cost = res.best.cost(inst);
    res.iters = it;
    res.seconds = now() - tStart;
    return res;
}

Result ils(const Instance& inst, const IlsOptions& o) {
    std::unique_ptr<LSContext> own;
    LSContext* ctx = o.ctx;
    if (!ctx) { own = std::make_unique<LSContext>(inst, o.k); ctx = own.get(); }
    ctx->ensureAdjacency();
    Preset ps = preset(o.preset);
    if (o.ops) ps.ops = opsFromNames(*o.ops);
    PyRandom rng(o.seed);
    Solution initSol = o.init ? *o.init : savingsMd(inst);

    RunOptions ro;
    ro.ops = ps.ops;
    ro.relocLen = ps.relocLen;
    ro.swapLen = ps.swapLen;

    LSState st0(*ctx, initSol, o.lamFleet);
    st0.run(ro, rng);
    Solution cur = st0.toSolution();
    double curObj = st0.objective();
    Solution best = cur;
    double bestObj = curObj;
    Result res;
    res.opStats = st0.stats;
    res.opGain = st0.gain;
    const double avg = avgEdge(inst, cur);
    const double T0 = o.t0Factor * avg, Tf = o.tfFactor * avg;
    const double cool = std::pow(Tf / T0, 1.0 / std::max(1, o.iters));
    double T = T0;
    const double tStart = now();
    res.history.push_back({0, 0.0, best.cost(inst)});
    if (o.target > 0 && bestObj - best.cost(inst) < 1e-6 && best.cost(inst) <= o.target * (1 + 1e-9)) res.reachedTarget = true;
    int it = 0;
    for (it = 1; it <= o.iters && !res.reachedTarget; it++) {
        RRState rr(*ctx, cur);
        std::vector<int> removed = rr.ruin(rng, o.cbar, o.lmax);
        rr.recreate(removed, rng, o.blink);
        Solution cand = rr.toSolution();
        std::vector<int> dirty = removed;
        PySetInt routeIds;                                           // iteration order of Python's set of route ids
        for (int c : removed) routeIds.add(rr.rid[c]);
        for (long long r : routeIds.items()) dirty.insert(dirty.end(), rr.routes[r].begin(), rr.routes[r].end());
        PySetInt uniq;
        for (int c : dirty) uniq.add(c);
        std::vector<int> dirtyList;
        for (long long c : uniq.items()) dirtyList.push_back((int)c);
        LSState st(*ctx, cand, o.lamFleet);
        ro.dirty = &dirtyList;
        st.run(ro, rng);
        for (int k = 0; k < OP_COUNT; k++) { res.opStats[k] += st.stats[k]; res.opGain[k] += st.gain[k]; }
        double obj = st.objective();
        if (obj < curObj - T * std::log(1.0 - rng.random())) {
            cur = st.toSolution();
            curObj = obj;
            if (obj < bestObj - EPS) {
                best = cur;
                bestObj = obj;
                const double bc = best.cost(inst);
                res.history.push_back({it, now() - tStart, bc});
                if (o.target > 0 && bestObj - bc < 1e-6 && bc <= o.target * (1 + 1e-9)) { res.reachedTarget = true; break; }
            }
        }
        T *= cool;
        if (o.timeLimit > 0 && now() - tStart > o.timeLimit) break;
    }
    if (it > o.iters) it = o.iters;
    res.best = best;
    res.cost = best.cost(inst);
    res.iters = it;
    res.seconds = now() - tStart;
    return res;
}

}  // namespace vrp
