#include "vrplearn/localsearch.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace vrp {

namespace {
constexpr double EPS = 1e-9;
const char* const OP_NAMES[OP_COUNT] = {"relocate", "swap", "two_opt", "two_opt_star", "or_opt",
                                        "swap_star", "route_depot", "depot_swap", "new_route"};
}  // namespace

const char* opName(int op) { return OP_NAMES[op]; }

std::optional<Op> opFromName(const std::string& name) {
    for (int o = 0; o < OP_COUNT; o++)
        if (name == OP_NAMES[o]) return static_cast<Op>(o);
    return std::nullopt;
}

Preset preset(const std::string& name) {
    const OpSet basic = opBit(OP_RELOCATE) | opBit(OP_SWAP) | opBit(OP_TWO_OPT) | opBit(OP_TWO_OPT_STAR) | opBit(OP_OR_OPT);
    if (name == "basic") return {basic, 1, 1};
    if (name == "hgs") return {basic | opBit(OP_SWAP_STAR), 2, 2};
    if (name == "full") return {basic | opBit(OP_SWAP_STAR) | opBit(OP_ROUTE_DEPOT) | opBit(OP_DEPOT_SWAP), 3, 2};
    throw std::invalid_argument("unknown preset " + name);
}

OpSet opsFromNames(const std::vector<std::string>& names) {
    OpSet s = 0;
    for (const auto& nm : names) {
        auto o = opFromName(nm);
        if (!o || *o == OP_NEW_ROUTE) throw std::invalid_argument("unknown operator " + nm);
        s |= opBit(*o);
    }
    return s;
}

// ---------------------------------------------------------------------------------------------
LSContext::LSContext(const Instance& I, int kk) : inst(I) {
    const int n = I.n;
    k = std::min(kk, n - 1);
    nbrs.resize(n);
    std::vector<int> idx(n);
    for (int i = 0; i < n; i++) {
        std::iota(idx.begin(), idx.end(), 0);
        std::stable_sort(idx.begin(), idx.end(), [&](int a, int b) {
            double da = a == i ? 1e18 : I.d(i, a), db = b == i ? 1e18 : I.d(i, b);
            return da < db;
        });
        nbrs[i].assign(idx.begin(), idx.begin() + k);
    }
}

void LSContext::ensureAdjacency() {
    if (!adj.empty()) return;
    const int n = inst.n;
    adj.resize(n);
    std::vector<int> idx(n);
    for (int i = 0; i < n; i++) {
        std::iota(idx.begin(), idx.end(), 0);
        std::stable_sort(idx.begin(), idx.end(), [&](int a, int b) { return inst.d(i, a) < inst.d(i, b); });
        adj[i] = idx;
    }
}

// ---------------------------------------------------------------------------------------------
long RouteData::nextUid = 0;

RouteData::RouteData(int d_, int dn_, std::vector<int> nodes_) : d(d_), dn(dn_), nodes(std::move(nodes_)), uid(nextUid++) {}

void RouteData::refresh(const Instance& I, long clock) {
    const size_t L = nodes.size();
    cq.assign(L + 1, 0.0);
    cs.assign(L + 1, 0.0);
    ce.assign(L + 1, 0.0);
    for (size_t k = 0; k < L; k++) {
        int v = nodes[k];
        cq[k + 1] = cq[k] + I.demand[v];
        cs[k + 1] = cs[k] + I.service[v];
        if (k >= 1) ce[k + 1] = ce[k] + I.d(nodes[k - 1], v);
    }
    load = cq[L];
    svc = cs[L];
    cost = L ? (I.d(dn, nodes[0]) + ce[L] + I.d(nodes[L - 1], dn)) : 0.0;
    stamp = clock;
}

// ---------------------------------------------------------------------------------------------
LSState::LSState(const LSContext& ctx, const Solution& sol, double lamFleet)
    : ctx_(ctx), inst_(ctx.inst), lam_(lamFleet), where_(ctx.inst.n, nullptr), pos_(ctx.inst.n, 0), count_(ctx.inst.t, 0) {
    const int n = inst_.n;
    for (size_t r = 0; r < sol.routes.size(); r++)
        if (!sol.routes[r].empty()) {
            int d = sol.depots[r];
            addRoute(std::make_unique<RouteData>(d, n + d, sol.routes[r]));
        }
}

void LSState::addRoute(std::unique_ptr<RouteData> R) {
    R->refresh(inst_, clock_);
    count_[R->d]++;
    RouteData* p = R.get();
    routes_.push_back(std::move(R));
    for (size_t k = 0; k < p->nodes.size(); k++) { where_[p->nodes[k]] = p; pos_[p->nodes[k]] = (int)k; }
}

void LSState::reindex(RouteData* R) {
    for (size_t k = 0; k < R->nodes.size(); k++) { where_[R->nodes[k]] = R; pos_[R->nodes[k]] = (int)k; }
}

Solution LSState::toSolution() const {
    Solution s;
    for (const auto& R : routes_)
        if (!R->nodes.empty() && !R->dead) { s.routes.push_back(R->nodes); s.depots.push_back(R->d); }
    return s;
}

double LSState::cost() const {
    std::vector<double> c;
    for (const auto& R : routes_)
        if (!R->nodes.empty() && !R->dead) c.push_back(R->cost);
    return pysum(c);
}

int LSState::excess() const {
    int e = 0;
    for (int d = 0; d < inst_.t; d++) e += std::max(0, count_[d] - inst_.fleet[d]);
    return e;
}

double LSState::softPen(int d, double q, double dur) const {
    double p = penQ_ * std::max(0.0, q - inst_.capacity[d]);
    if (std::isfinite(inst_.maxDuration[d])) p += penD_ * std::max(0.0, dur - inst_.maxDuration[d]);
    return p;
}

double LSState::excessLoad() const {
    double e = 0;
    for (const auto& R : routes_)
        if (!R->nodes.empty() && !R->dead) e += std::max(0.0, R->load - inst_.capacity[R->d]);
    return e;
}

double LSState::excessDuration() const {
    double e = 0;
    for (const auto& R : routes_)
        if (!R->nodes.empty() && !R->dead && std::isfinite(inst_.maxDuration[R->d])) e += std::max(0.0, R->cost + R->svc - inst_.maxDuration[R->d]);
    return e;
}

// ---- piece evaluation ------------------------------------------------------------------------
bool LSState::eval(int dn, const Pieces& ps, double& cost, double& q, double& s, int& cnt) const {
    const double* D = inst_.dist.data();
    const size_t N = inst_.size();
    int prev = dn;
    cost = 0; q = 0; s = 0; cnt = 0;
    for (int k = 0; k < ps.k; k++) {
        const Piece& pc = ps.p[k];
        if (pc.j <= pc.i) continue;
        const std::vector<int>& nd = pc.R->nodes;
        int first, last;
        if (pc.rev) { first = nd[pc.j - 1]; last = nd[pc.i]; }
        else { first = nd[pc.i]; last = nd[pc.j - 1]; }
        cost += D[prev * N + first] + pc.R->ce[pc.j] - pc.R->ce[pc.i + 1];
        prev = last;
        q += pc.R->cq[pc.j] - pc.R->cq[pc.i];
        s += pc.R->cs[pc.j] - pc.R->cs[pc.i];
        cnt += pc.j - pc.i;
    }
    cost += D[prev * N + dn];
    return true;
}

double LSState::fleetPen(const RouteData* A, int na, const RouteData* B, int nb) const {
    const std::vector<int>& m = inst_.fleet;
    int da = (na == 0 && !A->nodes.empty()) ? -1 : 0;
    int db = (nb > 0 && B->nodes.empty()) ? 1 : 0;
    int dA = A->d, dB = B->d;
    if (dA == dB)
        return lam_ * (std::max(0, count_[dA] + da + db - m[dA]) - std::max(0, count_[dA] - m[dA]));
    int pen = 0;
    if (da) pen += std::max(0, count_[dA] + da - m[dA]) - std::max(0, count_[dA] - m[dA]);
    if (db) pen += std::max(0, count_[dB] + db - m[dB]) - std::max(0, count_[dB] - m[dB]);
    return lam_ * pen;
}

bool LSState::delta2(RouteData* A, const Pieces& pa, RouteData* B, const Pieces& pb, double& out) const {
    double ca, qa, sa, cb, qb, sb;
    int na, nb;
    eval(A->dn, pa, ca, qa, sa, na);
    if (soft_) {                                                    // violations are priced, not rejected
        eval(B->dn, pb, cb, qb, sb, nb);
        double delta = ca + cb - A->cost - B->cost + softPen(A->d, qa, ca + sa) + softPen(B->d, qb, cb + sb) -
                       softPen(A->d, A->load, A->cost + A->svc) - softPen(B->d, B->load, B->cost + B->svc);
        if ((na == 0 && !A->nodes.empty()) || (nb > 0 && B->nodes.empty())) delta += fleetPen(A, na, B, nb);
        out = delta;
        return true;
    }
    if (qa > inst_.capacity[A->d] + EPS || ca + sa > inst_.maxDuration[A->d] + EPS) return false;
    eval(B->dn, pb, cb, qb, sb, nb);
    if (qb > inst_.capacity[B->d] + EPS || cb + sb > inst_.maxDuration[B->d] + EPS) return false;
    double delta = ca + cb - A->cost - B->cost;
    if ((na == 0 && !A->nodes.empty()) || (nb > 0 && B->nodes.empty())) delta += fleetPen(A, na, B, nb);
    out = delta;
    return true;
}

bool LSState::delta1(RouteData* A, const Pieces& pa, double& out) const {
    double ca, qa, sa;
    int na;
    eval(A->dn, pa, ca, qa, sa, na);
    if (soft_) { out = ca - A->cost + softPen(A->d, qa, ca + sa) - softPen(A->d, A->load, A->cost + A->svc); return true; }
    if (qa > inst_.capacity[A->d] + EPS || ca + sa > inst_.maxDuration[A->d] + EPS) return false;
    out = ca - A->cost;
    return true;
}

std::vector<int> LSState::build(const Pieces& ps) {
    std::vector<int> out;
    for (int k = 0; k < ps.k; k++) {
        const Piece& pc = ps.p[k];
        if (pc.j > pc.i) {
            if (pc.rev) for (int x = pc.j - 1; x >= pc.i; x--) out.push_back(pc.R->nodes[x]);
            else for (int x = pc.i; x < pc.j; x++) out.push_back(pc.R->nodes[x]);
        }
    }
    return out;
}

std::vector<RouteData*> LSState::apply(const Move& m) {
    ++clock_;
    RouteData *A = m.A, *B = m.B;
    const bool inter = B != A;
    double before = 0;
    int excessBefore = 0;
    if (verify) { before = cost(); excessBefore = excess(); }
    std::vector<int> newA = build(m.pa);
    std::vector<int> newB;
    if (inter) newB = build(m.pb);
    const bool hadA = !A->nodes.empty();
    const bool hadB = inter ? !B->nodes.empty() : true;
    A->nodes = newA;
    A->refresh(inst_, clock_);
    std::vector<RouteData*> touched{A};
    if (hadA && newA.empty()) { A->dead = true; count_[A->d]--; }
    if (inter) {
        B->nodes = newB;
        B->refresh(inst_, clock_);
        touched.push_back(B);
        if (!hadB && !newB.empty()) {
            // B is a fresh route object created by newRouteMove: take ownership
            for (auto& up : tmp_)
                if (up.get() == B) { routes_.push_back(std::move(up)); break; }
            count_[B->d]++;
        }
        if (hadB && newB.empty()) { B->dead = true; count_[B->d]--; }
    }
    for (RouteData* R : touched) if (!R->nodes.empty()) reindex(R);
    stats[m.op]++;
    gain[m.op] += -m.delta;
    if (verify) {
        double expected = m.delta - lam_ * (excess() - excessBefore);
        maxVerifyError = std::max(maxVerifyError, std::fabs((cost() - before) - expected));
    }
    return touched;
}

// ---- candidate generation for one customer pair (u, v) -----------------------------------------
void LSState::bestPair(int u, int v, OpSet ops, int relocLen, int swapLen, Move& best) {
    RouteData* A = where_[u];
    RouteData* B = where_[v];
    const int i = pos_[u], j = pos_[v];
    const int LA = (int)A->nodes.size(), LB = (int)B->nodes.size();

    auto consider = [&](Op op, double delta, RouteData* R1, const Pieces& p1, RouteData* R2, const Pieces& p2) {
        if (delta < best.delta) {
            if (admissible) {
                Move c;
                c.op = op; c.A = R1; c.B = R2; c.pa = p1; c.pb = p2; c.delta = delta; c.valid = true;
                if (!admissible(c)) return;
            }
            best.op = op; best.A = R1; best.B = R2; best.pa = p1; best.pb = p2; best.delta = delta; best.valid = true;
        }
    };
    double dl;

    if (A == B) {
        if ((ops & opBit(OP_TWO_OPT)) && i != j) {
            int lo = std::min(i, j), hi = std::max(i, j);
            if (hi - lo >= 2) {
                Pieces pa;
                pa.add(A, 0, lo + 1, false); pa.add(A, lo + 1, hi + 1, true); pa.add(A, hi + 1, LA, false);
                if (delta1(A, pa, dl)) consider(OP_TWO_OPT, dl, A, pa, A, Pieces());
            }
        }
        if (ops & opBit(OP_OR_OPT)) {
            for (int l = 1; l <= 3; l++) {
                if (i + l > LA) break;
                if (i <= j && j < i + l) continue;
                for (int rv = 0; rv < (l > 1 ? 2 : 1); rv++) {
                    bool rev = rv == 1;
                    for (int p : {j, j + 1}) {                         // before v / after v
                        Pieces pa;
                        if (j > i) {                                    // v after the segment
                            pa.add(A, 0, i, false); pa.add(A, i + l, p, false); pa.add(A, i, i + l, rev); pa.add(A, p, LA, false);
                        } else {
                            pa.add(A, 0, p, false); pa.add(A, i, i + l, rev); pa.add(A, p, i, false); pa.add(A, i + l, LA, false);
                        }
                        if (delta1(A, pa, dl)) consider(OP_OR_OPT, dl, A, pa, A, Pieces());
                    }
                }
            }
        }
        return;
    }

    // ---- inter-route ----------------------------------------------------------------------------
    if (ops & opBit(OP_RELOCATE)) {
        for (int l = 1; l <= relocLen; l++) {
            if (i + l > LA) break;
            for (int rv = 0; rv < (l > 1 ? 2 : 1); rv++) {
                bool rev = rv == 1;
                for (int p : {j, j + 1}) {
                    Pieces pa, pb;
                    pa.add(A, 0, i, false); pa.add(A, i + l, LA, false);
                    pb.add(B, 0, p, false); pb.add(A, i, i + l, rev); pb.add(B, p, LB, false);
                    if (delta2(A, pa, B, pb, dl)) consider(OP_RELOCATE, dl, A, pa, B, pb);
                }
            }
        }
    }
    if (ops & opBit(OP_SWAP)) {
        for (int l1 = 1; l1 <= swapLen; l1++) {
            if (i + l1 > LA) break;
            for (int l2 = 1; l2 <= swapLen; l2++) {
                if (j + l2 > LB) break;
                for (int rv = 0; rv < ((l1 > 1 || l2 > 1) ? 2 : 1); rv++) {
                    bool rev = rv == 1;
                    Pieces pa, pb;
                    pa.add(A, 0, i, false); pa.add(B, j, j + l2, rev); pa.add(A, i + l1, LA, false);
                    pb.add(B, 0, j, false); pb.add(A, i, i + l1, rev); pb.add(B, j + l2, LB, false);
                    if (delta2(A, pa, B, pb, dl)) consider(OP_SWAP, dl, A, pa, B, pb);
                }
            }
        }
    }
    if (ops & opBit(OP_TWO_OPT_STAR)) {
        Pieces pa, pb;
        pa.add(A, 0, i + 1, false); pa.add(B, j + 1, LB, false);       // exchange tails
        pb.add(B, 0, j + 1, false); pb.add(A, i + 1, LA, false);
        if (delta2(A, pa, B, pb, dl)) consider(OP_TWO_OPT_STAR, dl, A, pa, B, pb);
        if (A->d == B->d) {                                             // reversed heads: same depot only
            Pieces qa, qb;
            qa.add(A, 0, i + 1, false); qa.add(B, 0, j + 1, true);
            qb.add(A, i + 1, LA, true); qb.add(B, j + 1, LB, false);
            if (delta2(A, qa, B, qb, dl)) consider(OP_TWO_OPT_STAR, dl, A, qa, B, qb);
        }
    }
}

void LSState::newRouteMove(int u, int relocLen, Move& best) {
    RouteData* A = where_[u];
    const int i = pos_[u], LA = (int)A->nodes.size(), n = inst_.n;
    std::vector<int> ds(inst_.t);
    std::iota(ds.begin(), ds.end(), 0);
    std::stable_sort(ds.begin(), ds.end(), [&](int a, int b) { return inst_.d(u, n + a) < inst_.d(u, n + b); });
    if (ds.size() > 2) ds.resize(2);
    tmp_.clear();
    for (int d : ds) {
        tmp_.push_back(std::make_unique<RouteData>(d, n + d, std::vector<int>()));
        RouteData* B = tmp_.back().get();
        B->refresh(inst_, clock_);
        for (int l = 1; l <= relocLen; l++) {
            if (i + l > LA) break;
            Pieces pa, pb;
            pa.add(A, 0, i, false); pa.add(A, i + l, LA, false);
            pb.add(A, i, i + l, false);
            double dl;
            if (delta2(A, pa, B, pb, dl) && dl < best.delta) {
                best.op = OP_NEW_ROUTE; best.A = A; best.B = B; best.pa = pa; best.pb = pb; best.delta = dl; best.valid = true;
            }
        }
    }
}

// ---- SWAP* (Vidal 2022) -------------------------------------------------------------------------
namespace {
using Cand = std::pair<double, int>;
std::vector<Cand> top3(const Instance& I, int c, const RouteData* R) {
    std::vector<Cand> out;
    int prev = R->dn;
    const int L = (int)R->nodes.size();
    for (int p = 0; p <= L; p++) {
        int nxt = p < L ? R->nodes[p] : R->dn;
        out.push_back({I.d(prev, c) + I.d(c, nxt) - I.d(prev, nxt), p});
        prev = nxt;
    }
    std::sort(out.begin(), out.end());
    if (out.size() > 3) out.resize(3);
    return out;
}
}  // namespace

void LSState::swapStar(RouteData* A, RouteData* B, Move& best) const {
    const Instance& I = inst_;
    const std::vector<int>& an = A->nodes;
    const std::vector<int>& bn = B->nodes;
    const int LA = (int)an.size(), LB = (int)bn.size();
    std::vector<std::vector<Cand>> topB, topA;
    for (int u : an) topB.push_back(top3(I, u, B));          // insertion of A's customers into B
    for (int v : bn) topA.push_back(top3(I, v, A));          // insertion of B's customers into A
    best = Move();
    best.delta = -EPS;
    for (int i = 0; i < LA; i++) {
        int u = an[i];
        int pu = i ? an[i - 1] : A->dn;
        int nu = i + 1 < LA ? an[i + 1] : A->dn;
        double remA = I.d(pu, nu) - I.d(pu, u) - I.d(u, nu);
        for (int j = 0; j < LB; j++) {
            int v = bn[j];
            if (A->load - I.demand[u] + I.demand[v] > I.capacity[A->d] + EPS) continue;
            if (B->load - I.demand[v] + I.demand[u] > I.capacity[B->d] + EPS) continue;
            int pv = j ? bn[j - 1] : B->dn;
            int nv = j + 1 < LB ? bn[j + 1] : B->dn;
            double remB = I.d(pv, nv) - I.d(pv, v) - I.d(v, nv);
            std::vector<Cand> cands{{I.d(pu, v) + I.d(v, nu) - I.d(pu, nu), i}};
            for (const Cand& c : topA[j]) if (c.second != i && c.second != i + 1) cands.push_back(c);
            std::sort(cands.begin(), cands.end());
            bool okA = false;
            Cand insA{0, 0};
            for (const Cand& c : cands) {
                double cA = A->cost + remA + c.first;
                if (cA + A->svc - I.service[u] + I.service[v] <= I.maxDuration[A->d] + EPS) { insA = c; okA = true; break; }
            }
            if (!okA) continue;
            cands.assign(1, Cand{I.d(pv, u) + I.d(u, nv) - I.d(pv, nv), j});
            for (const Cand& c : topB[i]) if (c.second != j && c.second != j + 1) cands.push_back(c);
            std::sort(cands.begin(), cands.end());
            bool okB = false;
            Cand insB{0, 0};
            for (const Cand& c : cands) {
                double cB = B->cost + remB + c.first;
                if (cB + B->svc - I.service[v] + I.service[u] <= I.maxDuration[B->d] + EPS) { insB = c; okB = true; break; }
            }
            if (!okB) continue;
            double delta = remA + insA.first + remB + insB.first;
            if (delta < best.delta) {
                best.op = OP_SWAP_STAR; best.A = A; best.B = B; best.delta = delta; best.valid = true;
                best.ss[0] = i; best.ss[1] = j; best.ss[2] = insA.second; best.ss[3] = insB.second;
            }
        }
    }
}

std::vector<RouteData*> LSState::applySwapStar(const Move& m) {
    RouteData *A = m.A, *B = m.B;
    const int i = m.ss[0], j = m.ss[1], pA = m.ss[2], pB = m.ss[3];
    ++clock_;
    int u = A->nodes[i], v = B->nodes[j];
    double before = verify ? cost() : 0.0;
    std::vector<int> remA(A->nodes), remB(B->nodes);
    remA.erase(remA.begin() + i);
    remB.erase(remB.begin() + j);
    int pA2 = pA <= i ? pA : pA - 1;
    int pB2 = pB <= j ? pB : pB - 1;
    remA.insert(remA.begin() + pA2, v);
    remB.insert(remB.begin() + pB2, u);
    A->nodes = remA;
    B->nodes = remB;
    A->refresh(inst_, clock_);
    B->refresh(inst_, clock_);
    reindex(A);
    reindex(B);
    stats[OP_SWAP_STAR]++;
    gain[OP_SWAP_STAR] += -m.delta;
    if (verify) maxVerifyError = std::max(maxVerifyError, std::fabs((cost() - before) - m.delta));
    return {A, B};
}

// ---- multi-depot route-level moves -----------------------------------------------------------------
bool LSState::routeDepotPass(OpSet ops) {
    bool improved = false;
    const int t = inst_.t, n = inst_.n;
    if (t < 2) return false;
    auto liveRoutes = [&]() {
        std::vector<RouteData*> live;
        for (auto& R : routes_) if (!R->nodes.empty() && !R->dead) live.push_back(R.get());
        return live;
    };
    if (ops & opBit(OP_ROUTE_DEPOT)) {
        for (RouteData* R : liveRoutes()) {
            int first = R->nodes.front(), last = R->nodes.back();
            const int L = (int)R->nodes.size();
            double inner = R->ce[L];
            double bestDl = -EPS;
            int bestD = -1;
            for (int d2 = 0; d2 < t; d2++) {
                if (d2 == R->d || R->load > inst_.capacity[d2] + EPS) continue;
                int dn2 = n + d2;
                double c2 = inst_.d(dn2, first) + inner + inst_.d(last, dn2);
                if (c2 + R->svc > inst_.maxDuration[d2] + EPS) continue;
                const std::vector<int>& m = inst_.fleet;
                int penI = (std::max(0, count_[R->d] - 1 - m[R->d]) - std::max(0, count_[R->d] - m[R->d])) +
                           std::max(0, count_[d2] + 1 - m[d2]) - std::max(0, count_[d2] - m[d2]);
                double dl = c2 - R->cost + lam_ * penI;
                if (dl < bestDl) { bestDl = dl; bestD = d2; }
            }
            if (bestD >= 0) {
                ++clock_;
                count_[R->d]--;
                count_[bestD]++;
                R->d = bestD;
                R->dn = n + bestD;
                R->refresh(inst_, clock_);
                stats[OP_ROUTE_DEPOT]++;
                gain[OP_ROUTE_DEPOT] += -bestDl;
                improved = true;
            }
        }
    }
    if (ops & opBit(OP_DEPOT_SWAP)) {
        std::vector<RouteData*> live = liveRoutes();
        for (size_t a = 0; a < live.size(); a++)
            for (size_t b = a + 1; b < live.size(); b++) {
                RouteData *R1 = live[a], *R2 = live[b];
                if (R1->d == R2->d || R1->nodes.empty() || R2->nodes.empty()) continue;
                if (R1->load > inst_.capacity[R2->d] + EPS || R2->load > inst_.capacity[R1->d] + EPS) continue;
                double c1 = inst_.d(R2->dn, R1->nodes.front()) + R1->ce[R1->nodes.size()] + inst_.d(R1->nodes.back(), R2->dn);
                double c2 = inst_.d(R1->dn, R2->nodes.front()) + R2->ce[R2->nodes.size()] + inst_.d(R2->nodes.back(), R1->dn);
                if (c1 + R1->svc > inst_.maxDuration[R2->d] + EPS || c2 + R2->svc > inst_.maxDuration[R1->d] + EPS) continue;
                double dl = c1 + c2 - R1->cost - R2->cost;
                if (dl < -EPS) {
                    ++clock_;
                    std::swap(R1->d, R2->d);
                    std::swap(R1->dn, R2->dn);
                    R1->refresh(inst_, clock_);
                    R2->refresh(inst_, clock_);
                    stats[OP_DEPOT_SWAP]++;
                    gain[OP_DEPOT_SWAP] += -dl;
                    improved = true;
                }
            }
    }
    return improved;
}

// ---- main loop -----------------------------------------------------------------------------------------
void LSState::run(const RunOptions& opt, PyRandom& rng) {
    const OpSet ops = opt.ops ? opt.ops : preset("hgs").ops;
    const int n = inst_.n;
    std::vector<int> queue;
    if (opt.dirty) queue = *opt.dirty;
    else { queue.resize(n); std::iota(queue.begin(), queue.end(), 0); }
    rng.shuffle(queue);
    std::vector<char> inq(n, 0);
    for (int c : queue) inq[c] = 1;
    const auto& nbrs = ctx_.nbrs;
    std::unordered_map<unsigned long long, std::pair<long, long>> seen;       // SWAP* pair -> stamps at last scan
    long steps = 0;
    const long maxSteps = opt.maxPasses * n;
    auto enqueueRoute = [&](const RouteData* R) {
        for (int w : R->nodes) if (!inq[w]) { inq[w] = 1; queue.push_back(w); }
    };

    while (true) {
        while (!queue.empty() && steps < maxSteps) {
            ++steps;
            int u = queue.back();
            queue.pop_back();
            inq[u] = 0;
            RouteData* A = where_[u];
            if (!A || A->dead) continue;
            Move best;
            best.delta = -EPS;
            for (int v : nbrs[u]) bestPair(u, v, ops, opt.relocLen, opt.swapLen, best);
            if (opt.allowNewRoute && (ops & opBit(OP_RELOCATE))) newRouteMove(u, opt.relocLen, best);
            if (best.valid) {
                std::vector<RouteData*> touched = apply(best);
                for (RouteData* R : touched) enqueueRoute(R);
                if (!inq[u]) { inq[u] = 1; queue.push_back(u); }
            }
        }
        bool moved = false;
        if (ops & opBit(OP_SWAP_STAR)) {
            std::vector<std::pair<RouteData*, RouteData*>> pairs;
            std::unordered_set<unsigned long long> keys;
            const size_t take = std::max(5, ctx_.k / 2);
            auto key = [](long a, long b) { return (static_cast<unsigned long long>(a) << 32) | static_cast<unsigned long long>(b); };
            for (int u = 0; u < n; u++) {
                RouteData* Au = where_[u];
                if (!Au || Au->dead) continue;
                for (size_t q = 0; q < std::min(take, nbrs[u].size()); q++) {
                    RouteData* Bv = where_[nbrs[u][q]];
                    if (Bv == Au || !Bv || Bv->dead) continue;
                    RouteData *lo = Au->uid < Bv->uid ? Au : Bv, *hi = Au->uid < Bv->uid ? Bv : Au;
                    if (keys.insert(key(lo->uid, hi->uid)).second) pairs.push_back({lo, hi});
                }
            }
            for (auto& pr : pairs) {
                RouteData *A = pr.first, *B = pr.second;
                if (A->dead || B->dead || A->nodes.empty() || B->nodes.empty()) continue;
                unsigned long long k = key(A->uid, B->uid);
                auto it = seen.find(k);
                if (it != seen.end() && it->second == std::make_pair(A->stamp, B->stamp)) continue;
                Move mv;
                swapStar(A, B, mv);
                if (mv.valid) {
                    std::vector<RouteData*> touched = applySwapStar(mv);
                    for (RouteData* R : touched) enqueueRoute(R);
                    moved = true;
                }
                seen[k] = {A->stamp, B->stamp};
            }
        }
        if (routeDepotPass(ops)) {
            moved = true;
            for (auto& R : routes_)
                if (!R->dead) enqueueRoute(R.get());
        }
        if (!moved && queue.empty()) break;
        if (steps >= maxSteps) break;
    }
}

LSResult localSearch(const Instance& inst, const Solution& sol, const std::string& presetName,
                     const std::vector<std::string>* ops, int k, unsigned long long seed, const LSContext* ctx) {
    Preset ps = preset(presetName);
    if (ops) ps.ops = opsFromNames(*ops);
    std::unique_ptr<LSContext> own;
    if (!ctx) { own = std::make_unique<LSContext>(inst, k); ctx = own.get(); }
    LSState st(*ctx, sol);
    PyRandom rng(seed);
    RunOptions o;
    o.ops = ps.ops;
    o.relocLen = ps.relocLen;
    o.swapLen = ps.swapLen;
    st.run(o, rng);
    return {st.toSolution(), st.stats, st.gain};
}

}  // namespace vrp
