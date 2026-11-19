// Upgraded ("full") variants of the search methods.  Included by methods.cpp inside its anonymous namespace, so the helpers of
// the compact methods (Tracker, objective, descend, destroy / repair, split, ...) are shared.  The compact methods are not
// changed by this file; each variant here has its own id with the suffix "-f".
//
// Parameters are first values in the range of the cited papers; they are meant to be tuned (the compact and the full variants
// are tuned with the same budget before they are compared).

// ===========================================================================================
// tabu-f: granular tabu search for the CVRP / MDVRP.
//   Toth & Vigo (2003): only moves that create at least one "short" arc are examined; an arc is short if its cost is at most
//     theta = beta * z / (n + K)  (z, K: cost and number of routes of the start solution; beta: sparsification factor).
//     The whole granular neighbourhood is evaluated in every iteration and the best admissible move is executed, even if it
//     worsens the solution.  When the best solution has not improved for a while, beta is raised for some iterations
//     (diversification) and then reset.
//   Gendreau, Hertz & Laporte (1994): solutions that violate capacity or duration are visited; the excess is priced by two
//     penalties that are doubled when the last h solutions were all infeasible and halved when they were all feasible.
//   Tabu rule: an arc removed by a move may not be re-created for a random number of iterations in [tmin, tmax]; a tabu move is
//     accepted if it leads below the best cost found (aspiration).
//   Moves (all evaluated in O(1) from prefix sums): relocate of 1-3 consecutive customers (also reversed), exchange of segments
//     of 1-2 customers, 2-opt, 2-opt* (tail exchange; reversed heads for a common depot), Or-opt, opening a new route.
//   Implementation: the best move of every customer is cached and recomputed only when one of the routes it involves has
//     changed; all caches are refreshed every h iterations (tabu expiry, new penalties).  Moves that leave the solution
//     unchanged (a segment re-inserted where it was) are never selected.
//   One iteration = one executed move.
// ===========================================================================================
Result runTabuFull(const Instance& I, const MethodOptions& o, LSContext& ctx, const Solution& init) {
    ctx.ensureAdjacency();
    PyRandom rng(o.seed);
    Tracker T(I, init, o.lamFleet, o.target);
    const int n = I.n, N = I.n + I.t;
    // ---- parameters
    const double beta0 = par(o, "beta0", 1.5), betaDiv = par(o, "betadiv", 3.0);   // sparsification factor: normal / during diversification
    const int kMin = 3, kMax = std::max(kMin, ipar(o, "kmax", 40));              // neighbours per customer are clipped to [kMin, kMax]
    const int tMin = ipar(o, "tmin", 5), tMax = tMin + std::max(0, ipar(o, "tspan", 5));   // tabu tenure
    const int h = std::max(1, ipar(o, "h", 10));                   // penalties are updated every h iterations
    const int divAfter = std::max(100, (int)(par(o, "divfactor", 2) * n)), divLen = std::max(25, n / 2);
    const int relocLen = 3, swapLen = 2;
    const OpSet ops = (opBit(OP_RELOCATE) | opBit(OP_SWAP) | opBit(OP_TWO_OPT) | opBit(OP_TWO_OPT_STAR) | opBit(OP_OR_OPT)) & ~o.excludeOps;

    LSState st(ctx, init, o.lamFleet);
    const double z0 = std::max(init.cost(I), 1e-9), arc = z0 / (n + std::max(1, init.nRoutes()));
    const double penQ0 = arc / std::max(I.totalDemand() / n, 1e-9);  // one average demand of excess costs one average arc
    double penQ = penQ0, penD = 1.0;
    st.setSoft(penQ, penD);

    // granular neighbour lists for a sparsification factor
    std::vector<std::vector<int>> gran(n), revGran(n);            // revGran[v]: customers that have v as a neighbour
    std::vector<char> dirty(n, 1);
    auto buildGran = [&](double beta) {
        const double theta = beta * arc;
        for (int u = 0; u < n; u++) revGran[u].clear();
        std::fill(dirty.begin(), dirty.end(), 1);
        for (int u = 0; u < n; u++) {
            gran[u].clear();
            for (int c : ctx.adj[u]) {
                if (c == u) continue;
                if ((int)gran[u].size() >= kMax || ((int)gran[u].size() >= kMin && I.d(u, c) > theta)) break;
                gran[u].push_back(c);
                revGran[c].push_back(u);
            }
        }
    };
    buildGran(beta0);

    // arc-based tabu memory (undirected arcs between nodes; depots are nodes n .. n+t-1)
    const bool dense = (long long)N * N <= 30000000LL;
    std::vector<int> tabuDense(dense ? (size_t)N * N : 0, 0);
    std::unordered_map<long long, int> tabuMap;
    auto key = [&](int a, int b) { return a < b ? (long long)a * N + b : (long long)b * N + a; };
    auto tabuUntil = [&](int a, int b) -> int {
        const long long k = key(a, b);
        if (dense) return tabuDense[(size_t)k];
        auto f = tabuMap.find(k);
        return f == tabuMap.end() ? 0 : f->second;
    };
    auto setTabu = [&](long long k, int until) { if (dense) tabuDense[(size_t)k] = until; else tabuMap[k] = until; };
    auto adjacent = [&](int a, int b) {                             // is the arc {a, b} part of the current solution?
        if (a >= n && b >= n) return true;                          // depot - depot: an empty route, never tabu
        if (a >= n) std::swap(a, b);                                // a is a customer
        const RouteData* R = st.routeOf(a);
        if (b >= n) return R->dn == b && (st.posOf(a) == 0 || st.posOf(a) == (int)R->nodes.size() - 1);
        return st.routeOf(b) == R && std::abs(st.posOf(a) - st.posOf(b)) == 1;
    };
    int it = 0;
    double curObj = 0;
    // arcs created by a route given as pieces: is any of them new, is any of them tabu?
    auto scanArcs = [&](int dn, const Pieces& ps, bool& anyNew, bool& anyTabu) {
        int prev = dn;
        auto look = [&](int a, int b) {
            if (adjacent(a, b)) return;
            anyNew = true;
            if (tabuUntil(a, b) > it) anyTabu = true;
        };
        for (int k = 0; k < ps.k; k++) {
            const Piece& pc = ps.p[k];
            if (pc.j <= pc.i) continue;
            const std::vector<int>& nd = pc.R->nodes;
            look(prev, pc.rev ? nd[pc.j - 1] : nd[pc.i]);
            prev = pc.rev ? nd[pc.i] : nd[pc.j - 1];
        }
        if (prev != dn) look(prev, dn);
    };
    auto admissible = [&](const Move& m) {
        bool anyNew = false, anyTabu = false;
        scanArcs(m.A->dn, m.pa, anyNew, anyTabu);
        if (m.B != m.A) scanArcs(m.B->dn, m.pb, anyNew, anyTabu);
        if (!anyNew) return false;                                  // the move would leave the solution unchanged
        return !anyTabu || curObj + m.delta < T.bestO - EPS;        // aspiration
    };
    st.admissible = admissible;
    auto arcsOf = [&](const RouteData* R, std::vector<long long>& out) {
        if (R->nodes.empty()) return;
        out.push_back(key(R->dn, R->nodes.front()));
        for (size_t k = 0; k + 1 < R->nodes.size(); k++) out.push_back(key(R->nodes[k], R->nodes[k + 1]));
        out.push_back(key(R->nodes.back(), R->dn));
    };
    auto evaluate = [&](double& exQ, double& exD) {
        exQ = st.excessLoad();
        exD = st.excessDuration();
        return st.cost() + o.lamFleet * st.excess();
    };
    double exQ, exD;
    curObj = evaluate(exQ, exD);                                    // the start is feasible for capacity and duration

    int infQ = 0, infD = 0, nonImp = 0, divLeft = 0;
    std::vector<long long> before, after;
    std::vector<Move> cache(n);
    std::vector<double> nrCache(n, INF);
    for (it = 1; it <= o.iters; it++) {
        T.it = it;
        int bu = -1, nrU = -1;
        for (int u = 0; u < n; u++) {
            if (dirty[u]) {
                dirty[u] = 0;
                Move& b = cache[u];
                b = Move();
                b.delta = INF;
                for (int v : gran[u]) st.scanPair(u, v, ops, relocLen, swapLen, b);
                Move m;
                m.delta = INF;
                st.scanNewRoute(u, relocLen, m);
                nrCache[u] = m.valid && admissible(m) ? m.delta : INF;
            }
            if (cache[u].valid && (bu < 0 || cache[u].delta < cache[bu].delta)) bu = u;
            if (nrCache[u] < INF && (nrU < 0 || nrCache[u] < nrCache[nrU])) nrU = u;
        }
        Move nr;
        nr.delta = INF;
        if (nrU >= 0 && (bu < 0 || nrCache[nrU] < cache[bu].delta)) st.scanNewRoute(nrU, relocLen, nr);   // regenerate: its route object is temporary
        const Move mv = nr.valid ? nr : bu >= 0 ? cache[bu] : Move();
        if (mv.valid) {
            before.clear();
            after.clear();
            arcsOf(mv.A, before);
            if (mv.B != mv.A) arcsOf(mv.B, before);
            const std::vector<RouteData*> touched = st.applyMove(mv);
            bool fleetChanged = mv.op == OP_NEW_ROUTE;
            for (const RouteData* R : touched) {
                if (R->dead) { fleetChanged = true; continue; }
                arcsOf(R, after);
                for (int c : R->nodes) { dirty[c] = 1; for (int w : revGran[c]) dirty[w] = 1; }
            }
            if (fleetChanged) std::fill(dirty.begin(), dirty.end(), 1);          // the fleet penalty of other moves changed
            std::sort(after.begin(), after.end());
            for (long long k : before)
                if (!std::binary_search(after.begin(), after.end(), k)) setTabu(k, it + tMin + rng.randrange(tMax - tMin + 1));
        }
        const double feasObj = evaluate(exQ, exD);
        const bool okQ = exQ < 1e-9, okD = exD < 1e-9;
        if (okQ && okD && feasObj < T.bestO - EPS) { T.offer(st.toSolution(), feasObj); nonImp = 0; }
        else nonImp++;
        infQ += !okQ;
        infD += !okD;
        if (it % h == 0) {                                          // self-adjusting penalties
            if (infQ == h) penQ = std::min(penQ * 2, penQ0 * 1e4); else if (infQ == 0) penQ = std::max(penQ / 2, penQ0 * 1e-3);
            if (infD == h) penD = std::min(penD * 2, 1e4); else if (infD == 0) penD = std::max(penD / 2, 1e-3);
            infQ = infD = 0;
            st.setSoft(penQ, penD);
            std::fill(dirty.begin(), dirty.end(), 1);                // periodic refresh: new penalties, expired tabu status
        }
        curObj = feasObj + penQ * exQ + penD * exD;
        if (divLeft > 0) { if (--divLeft == 0) buildGran(beta0); }
        else if (nonImp >= divAfter) { buildGran(betaDiv); divLeft = divLen; nonImp = 0; }
        if (!dense && tabuMap.size() > 200000)
            for (auto f = tabuMap.begin(); f != tabuMap.end();) f = f->second <= it ? tabuMap.erase(f) : std::next(f);
        if (T.timeUp(o.timeLimit)) break;
    }
    account(st);                                                    // moves executed per operator (the gain includes worsening moves)
    return T.finish(std::min(it, o.iters));
}

// ===========================================================================================
// alns-f: adaptive large neighbourhood search with the operator set of Pisinger & Ropke (2007) / Ropke & Pisinger (2006).
//   Removal:   random; worst (largest removal gain, randomised with p = 3); Shaw (relatedness = 9 x distance + 2 x demand
//              difference, both normalised, + 3 if the two customers are on different routes; p = 6); cluster (a route is split
//              into two clusters by Kruskal's algorithm and one of them is removed, then the search moves to a nearby route);
//              historical node-pair (every arc remembers the best objective of a solution it appeared in; customers whose two
//              arcs have the worst record are removed, p = 3).
//   Insertion: greedy (cheapest insertion over all customers), regret-2, regret-3, regret-m (m = number of routes); each with or
//              without noise on the insertion cost (uniform in +- 0.025 x the largest distance), chosen adaptively as well.
//   The number of removed customers is uniform in [4, min(100, 0.4 n)].
//   Acceptance: simulated annealing; the start temperature accepts a solution 5 % worse than the start with probability 0.5
//              and falls geometrically to 0.00193 times that value (= 0.99975^25000, the schedule of the paper) at the last iteration.
//   Adaptation: roulette weights per removal operator, insertion operator and noise setting, updated every 100 iterations with
//              reaction 0.1 from the scores 33 (new best), 9 (new solution better than the current one), 13 (new solution, worse,
//              accepted); only solutions that were not accepted before score.
//   Insertion deltas are cached per (customer, route) and only the column of the changed route is recomputed.
//   One iteration = one removal + insertion.
// ===========================================================================================
// insertion with regret over the K best routes (K = 1: greedy); noise > 0 perturbs the insertion costs
void repairRegretFull(Part& P, const std::vector<int>& left, int K, double noise, PyRandom& rng) {
    const Instance& I = P.I;
    const int q = (int)left.size();
    K = std::max(1, std::min(K, (int)P.routes.size()));
    auto calc = [&](int c, int r) -> std::pair<double, int> {
        int p;
        double dl;
        if (!bestPos(I, P.routes[r], c, p, dl)) return {INF, -1};
        if (noise > 0) dl = std::max(0.0, dl + noise * (2 * rng.random() - 1));
        return {dl, p};
    };
    std::vector<std::vector<std::pair<double, int>>> ins(q);
    for (int k = 0; k < q; k++) {
        ins[k].resize(P.routes.size());
        for (int r = 0; r < (int)P.routes.size(); r++) ins[k][r] = calc(left[k], r);
    }
    std::vector<char> done(q, 0);
    std::vector<double> b(K);
    for (int remaining = q; remaining > 0; remaining--) {
        int pk = -1, pr = -1, stuck = -1;
        double pv = -INF;
        for (int k = 0; k < q; k++) {
            if (done[k]) continue;
            std::fill(b.begin(), b.end(), INF);
            int r1 = -1;
            for (int r = 0; r < (int)ins[k].size(); r++) {
                const double d = ins[k][r].first;
                if (d >= b[K - 1]) continue;
                int j = K - 1;
                while (j > 0 && b[j - 1] > d) { b[j] = b[j - 1]; j--; }
                b[j] = d;
                if (j == 0) r1 = r;
            }
            if (r1 < 0) { if (stuck < 0) stuck = k; continue; }
            double val = K == 1 ? -b[0] : -1e-6 * b[0];              // ties of the regret: cheapest insertion first
            for (int j = 1; j < K; j++) val += (b[j] < INF ? b[j] : 1e9) - b[0];
            if (val > pv) { pv = val; pk = k; pr = r1; }
        }
        if (pk < 0) {                                                // no customer fits anywhere: open a route
            pk = stuck;
            openRoute(P, left[pk]);
            done[pk] = 1;
            const int r = (int)P.routes.size() - 1;
            for (int k = 0; k < q; k++) ins[k].push_back(done[k] ? std::make_pair(INF, -1) : calc(left[k], r));
            continue;
        }
        insertAt(P, pr, ins[pk][pr].second, left[pk]);
        done[pk] = 1;
        for (int k = 0; k < q; k++) if (!done[k]) ins[k][pr] = calc(left[k], pr);
    }
}

std::vector<int> destroyWorstOnce(Part& P, int q, PyRandom& rng, double p) {
    std::vector<std::pair<double, int>> cand;
    for (const PR& R : P.routes)
        for (int i = 0; i < (int)R.nodes.size(); i++) cand.push_back({removalGain(P.I, R, i), R.nodes[i]});
    std::stable_sort(cand.begin(), cand.end(), [](const auto& a, const auto& b) { return a.first > b.first; });
    std::vector<int> rm;
    while ((int)rm.size() < q && !cand.empty()) {
        const int idx = std::min((int)cand.size() - 1, (int)std::floor(std::pow(rng.random(), p) * cand.size()));
        rm.push_back(cand[idx].second);
        cand.erase(cand.begin() + idx);
    }
    for (int c : rm) P.remove(c);
    return rm;
}

std::vector<int> destroyShawFull(Part& P, int q, PyRandom& rng, const LSContext& ctx, double maxD, double maxQ) {
    const Instance& I = P.I;
    std::vector<int> rm = {rng.randrange(I.n)};
    std::vector<char> gone(I.n, 0);
    gone[rm[0]] = 1;
    std::vector<std::pair<double, int>> cand;
    while ((int)rm.size() < q) {
        const int r = rm[rng.randrange((int)rm.size())];
        cand.clear();
        for (int c : ctx.adj[r]) {
            if (gone[c]) continue;
            cand.push_back({9 * I.d(r, c) / maxD + 2 * std::fabs(I.demand[r] - I.demand[c]) / maxQ + (P.where[c] == P.where[r] ? 0.0 : 3.0), c});
            if (cand.size() >= 100) break;
        }
        if (cand.empty()) break;
        std::stable_sort(cand.begin(), cand.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
        const int c = cand[std::min((int)cand.size() - 1, (int)std::floor(std::pow(rng.random(), 6.0) * cand.size()))].second;
        gone[c] = 1;
        rm.push_back(c);
    }
    for (int c : rm) P.remove(c);
    return rm;
}

std::vector<int> destroyCluster(Part& P, int q, PyRandom& rng, const LSContext& ctx) {
    const Instance& I = P.I;
    std::vector<int> rm;
    std::vector<char> used(P.routes.size(), 0), gone(I.n, 0);
    int seed = rng.randrange(I.n);
    while ((int)rm.size() < q) {
        int r = used[P.where[seed]] ? -1 : P.where[seed];
        if (r < 0)                                                   // the nearest customer on a route that was not split yet
            for (int c : ctx.adj[seed]) if (!gone[c] && !used[P.where[c]]) { r = P.where[c]; break; }
        if (r < 0) break;
        used[r] = 1;
        const std::vector<int>& nodes = P.routes[r].nodes;
        const int L = (int)nodes.size();
        std::vector<int> take;
        if (L <= 2) take = nodes;
        else {                                                       // Kruskal until two components remain
            struct E { double d; int a, b; };
            std::vector<E> es;
            for (int a = 0; a < L; a++) for (int b = a + 1; b < L; b++) es.push_back({I.d(nodes[a], nodes[b]), a, b});
            std::sort(es.begin(), es.end(), [](const E& x, const E& y) { return x.d < y.d; });
            std::vector<int> comp(L);
            std::iota(comp.begin(), comp.end(), 0);
            std::function<int(int)> find = [&](int x) { return comp[x] == x ? x : comp[x] = find(comp[x]); };
            int comps = L;
            for (const E& e : es) {
                if (comps == 2) break;
                const int a = find(e.a), b = find(e.b);
                if (a != b) { comp[a] = b; comps--; }
            }
            const int pick = find(rng.randrange(L));
            for (int a = 0; a < L; a++) if (find(a) == pick) take.push_back(nodes[a]);
        }
        for (int c : take) { gone[c] = 1; rm.push_back(c); }
        seed = take[rng.randrange((int)take.size())];
    }
    for (int c : rm) P.remove(c);
    return rm;
}

Result runAlnsFull(const Instance& I, const MethodOptions& o, LSContext& ctx, const Solution& init) {
    ctx.ensureAdjacency();
    PyRandom rng(o.seed);
    Tracker T(I, init, o.lamFleet, o.target);
    const int n = I.n, N = I.n + I.t;
    // ---- parameters (Ropke & Pisinger 2006, Pisinger & Ropke 2007)
    const double w = par(o, "w", 0.05), coolTotal = par(o, "cool", 0.00193), react = par(o, "react", 0.1), eta = par(o, "eta", 0.025), xi = par(o, "xi", 0.4);
    const int seg = std::max(1, ipar(o, "seg", 100)), qMinAbs = 4, qMaxAbs = ipar(o, "qmax", 100);
    const double s1 = par(o, "s1", 33), s2 = par(o, "s2", 9), s3 = par(o, "s3", 13);
    double maxD = 0, maxQ = 1e-9;
    for (double d : I.dist) maxD = std::max(maxD, d);
    for (double q : I.demand) maxQ = std::max(maxQ, q);
    const int qa = std::max(1, std::min(qMinAbs, n - 1)), qb = std::max(qa, std::min(qMaxAbs, (int)(xi * n)));

    Solution cur = init.nonEmpty();
    double curO = objective(I, cur, o.lamFleet);
    const double Tstart = -w * std::max(cur.cost(I), 1e-9) / std::log(0.5), cool = std::pow(coolTotal, 1.0 / std::max(1, o.iters));
    double temp = Tstart;

    const bool useHist = (long long)N * N <= 30000000LL;             // historical node-pair removal needs an N x N table
    std::vector<float> H(useHist ? (size_t)N * N : 0, std::numeric_limits<float>::max());
    auto record = [&](const Solution& s, double obj) {
        if (!useHist) return;
        for (size_t r = 0; r < s.routes.size(); r++) {
            int a = I.n + s.depots[r];
            for (int c : s.routes[r]) { float& h = H[(size_t)a * N + c]; if (obj < h) h = (float)obj; a = c; }
            float& h = H[(size_t)a * N + I.n + s.depots[r]];
            if (obj < h) h = (float)obj;
        }
    };
    auto destroyHistory = [&](Part& P, int q) {
        std::vector<std::pair<double, int>> cand;
        for (const PR& R : P.routes) {
            const int dn = I.n + R.d, L = (int)R.nodes.size();
            for (int i = 0; i < L; i++) {
                const int c = R.nodes[i], pv = i ? R.nodes[i - 1] : dn, nx = i + 1 < L ? R.nodes[i + 1] : dn;
                cand.push_back({(double)H[(size_t)pv * N + c] + (double)H[(size_t)c * N + nx], c});
            }
        }
        std::stable_sort(cand.begin(), cand.end(), [](const auto& a, const auto& b) { return a.first > b.first; });
        std::vector<int> rm;
        while ((int)rm.size() < q && !cand.empty()) {
            const int idx = std::min((int)cand.size() - 1, (int)std::floor(std::pow(rng.random(), 3.0) * cand.size()));
            rm.push_back(cand[idx].second);
            cand.erase(cand.begin() + idx);
        }
        for (int c : rm) P.remove(c);
        return rm;
    };
    auto hashOf = [&](const Solution& s) {
        unsigned long long hsh = 0;
        auto mix = [](unsigned long long x) { x += 0x9e3779b97f4a7c15ULL; x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL; x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL; return x ^ (x >> 31); };
        for (size_t r = 0; r < s.routes.size(); r++) {
            int a = I.n + s.depots[r];
            for (int c : s.routes[r]) { hsh += mix((unsigned long long)a * N + c); a = c; }
            hsh += mix((unsigned long long)a * N + I.n + s.depots[r]);
        }
        return hsh;
    };
    std::unordered_map<unsigned long long, char> visited;
    visited[hashOf(cur)] = 1;
    record(cur, curO);

    const int ND = useHist ? 5 : 4, NR = 4, NN = 2;
    std::vector<double> wd(ND, 1.0), wr(NR, 1.0), wn(NN, 1.0), scd(ND, 0.0), scr(NR, 0.0), scn(NN, 0.0), usd(ND, 0.0), usr(NR, 0.0), usn(NN, 0.0);
    std::vector<long> totD(ND, 0), totR(NR, 0), totN(NN, 0);
    auto roulette = [&](const std::vector<double>& ws) {
        double sum = 0;
        for (double x : ws) sum += x;
        double x = rng.random() * sum;
        for (int i = 0; i < (int)ws.size(); i++) { x -= ws[i]; if (x <= 0) return i; }
        return (int)ws.size() - 1;
    };
    auto adapt = [&](std::vector<double>& ws, std::vector<double>& sc, std::vector<double>& us) {
        for (size_t i = 0; i < ws.size(); i++) { if (us[i] > 0) ws[i] = std::max(0.01, (1 - react) * ws[i] + react * sc[i] / us[i]); sc[i] = us[i] = 0; }
    };
    int it;
    for (it = 1; it <= o.iters; it++) {
        T.it = it;
        const int di = roulette(wd), ri = roulette(wr), ni = roulette(wn), q = qa + rng.randrange(qb - qa + 1);
        Part P(I, cur);
        std::vector<int> rm = di == 0 ? destroyRandom(P, q, rng) : di == 1 ? destroyWorstOnce(P, q, rng, 3.0) : di == 2 ? destroyShawFull(P, q, rng, ctx, maxD, maxQ)
                              : di == 3 ? destroyCluster(P, q, rng, ctx) : destroyHistory(P, q);
        P.reindex();
        repairRegretFull(P, rm, ri == 0 ? 1 : ri == 1 ? 2 : ri == 2 ? 3 : (int)P.routes.size(), ni == 1 ? eta * maxD : 0.0, rng);
        P.reindex();
        Solution cand = P.toSolution();
        const double oc = objective(I, cand, o.lamFleet);
        record(cand, oc);
        double score = 0;
        if (oc < curO - temp * std::log(1.0 - rng.random())) {
            const bool isNew = visited.emplace(hashOf(cand), 1).second;
            if (oc < T.bestO - EPS) { T.offer(cand, oc); score = s1; }
            else if (isNew) score = oc < curO - EPS ? s2 : s3;
            cur = std::move(cand);
            curO = oc;
            if (visited.size() > 2000000) visited.clear();
        }
        scd[di] += score; scr[ri] += score; scn[ni] += score; usd[di]++; usr[ri]++; usn[ni]++; totD[di]++; totR[ri]++; totN[ni]++;
        if (it % seg == 0) { adapt(wd, scd, usd); adapt(wr, scr, usr); adapt(wn, scn, usn); }
        temp *= cool;
        if (T.timeUp(o.timeLimit)) break;
    }
    Result res = T.finish(std::min(it, o.iters));
    const char* dn[] = {"removal:random", "removal:worst", "removal:shaw", "removal:cluster", "removal:history"};
    const char* rn[] = {"insertion:greedy", "insertion:regret2", "insertion:regret3", "insertion:regretm"};
    const char* nn[] = {"noise:off", "noise:on"};
    for (int i = 0; i < ND; i++) res.extra.push_back({dn[i], totD[i], wd[i]});
    for (int i = 0; i < NR; i++) res.extra.push_back({rn[i], totR[i], wr[i]});
    for (int i = 0; i < NN; i++) res.extra.push_back({nn[i], totN[i], wn[i]});
    return res;
}

// ===========================================================================================
// ga-f: hybrid genetic search with adaptive diversity control (Vidal, Crainic, Gendreau, Lahrichi & Rei 2012; Vidal 2022).
//   Added to the compact memetic GA:
//   * two sub-populations, feasible and infeasible: the education (local search) prices excess load and excess duration with
//     penalties, so offspring may be infeasible; an infeasible offspring is repaired with probability 0.5 (education with ten
//     times the penalties).  The penalties are adapted every 100 offspring to keep about 20 % of them feasible (x 1.2 / x 0.85).
//   * biased fitness = rank by penalised cost + (1 - nElite / size) x rank by diversity contribution (mean broken-pairs distance
//     to the nClose = 5 closest individuals of the sub-population; depot genes that differ count as well).
//   * parents by binary tournament on the biased fitness over both sub-populations; survivor selection when a sub-population
//     reaches mu + lambda = 25 + 40: clones first, then the worst biased fitness, down to mu.
//   * diversification: after 4000 offspring without a new best solution the best third of each sub-population is kept and
//     4 mu new individuals are generated.
//   Kept from the compact GA: giant tour + depot genes, order crossover (depot genes inherited uniformly), optimal split.
//   One iteration = one educated individual (initial individuals included).
// ===========================================================================================
struct HInd {
    Solution sol;
    std::vector<int> tour;
    std::vector<char> dep;
    double dist = 0, exQ = 0, exD = 0, pc = 0, bf = 0;
    int fleetEx = 0;
    bool feasible() const { return exQ < 1e-9 && exD < 1e-9; }
};
struct HPop {
    std::vector<HInd> v;
    std::vector<std::vector<int>> D;                                 // pairwise distances
    bool stale = true;
    static int distance(const HInd& a, const HInd& b) {
        int d = brokenPairs(a.tour, b.tour);
        for (size_t c = 0; c < a.dep.size(); c++) d += a.dep[c] != b.dep[c];
        return d;
    }
    void add(HInd x) {
        std::vector<int> row;
        for (const HInd& y : v) row.push_back(distance(x, y));
        for (size_t i = 0; i < D.size(); i++) D[i].push_back(row[i]);
        row.push_back(0);
        D.push_back(std::move(row));
        v.push_back(std::move(x));
        stale = true;
    }
    void removeAt(int i) {
        v.erase(v.begin() + i);
        D.erase(D.begin() + i);
        for (auto& r : D) r.erase(r.begin() + i);
        stale = true;
    }
    double closeness(int i, int nClose) const {
        std::vector<int> d;
        for (int j = 0; j < (int)v.size(); j++) if (j != i) d.push_back(D[i][j]);
        if (d.empty()) return 0;
        const int k = std::min(nClose, (int)d.size());
        std::partial_sort(d.begin(), d.begin() + k, d.end());
        double s = 0;
        for (int j = 0; j < k; j++) s += d[j];
        return s / k;
    }
    void rank(int nElite, int nClose) {
        if (!stale) return;
        stale = false;
        const int P = (int)v.size();
        if (P == 0) return;
        if (P == 1) { v[0].bf = 0; return; }
        std::vector<int> byCost(P), byDiv(P);
        std::iota(byCost.begin(), byCost.end(), 0);
        byDiv = byCost;
        std::vector<double> cl(P);
        for (int i = 0; i < P; i++) cl[i] = closeness(i, nClose);
        std::stable_sort(byCost.begin(), byCost.end(), [&](int a, int b) { return v[a].pc < v[b].pc; });
        std::stable_sort(byDiv.begin(), byDiv.end(), [&](int a, int b) { return cl[a] > cl[b]; });
        const double wDiv = P > nElite ? 1.0 - (double)nElite / P : 0.0;
        for (int r = 0; r < P; r++) v[byCost[r]].bf = (double)r / (P - 1);
        for (int r = 0; r < P; r++) v[byDiv[r]].bf += wDiv * r / (P - 1);
    }
    void survivors(int mu, int nElite, int nClose) {
        while ((int)v.size() > mu) {
            rank(nElite, nClose);
            int worst = -1;
            bool worstClone = false;
            for (int i = 0; i < (int)v.size(); i++) {
                bool clone = false;
                for (int j = 0; j < (int)v.size() && !clone; j++) clone = j != i && D[i][j] == 0;
                if (worst < 0 || (clone && !worstClone) || (clone == worstClone && v[i].bf > v[worst].bf)) { worst = i; worstClone = clone; }
            }
            removeAt(worst);
        }
    }
    void keepBest(int k) {
        while ((int)v.size() > k) {
            int w = 0;
            for (int i = 1; i < (int)v.size(); i++) if (v[i].pc > v[w].pc) w = i;
            removeAt(w);
        }
    }
};

Result runGaFull(const Instance& I, const MethodOptions& o, const LSContext& ctx, const Solution& init) {
    PyRandom rng(o.seed);
    Tracker T(I, init, o.lamFleet, o.target);
    const DepotOrder dord = buildDepotOrder(I);
    // ---- parameters (Vidal et al. 2012; Vidal 2022)
    const int mu = std::max(4, ipar(o, "mu", 25)), lambda = std::max(1, ipar(o, "lambda", 40)), nElite = ipar(o, "nelite", 4), nClose = ipar(o, "nclose", 5),
              divAfter = ipar(o, "divafter", 4000), window = 100;
    const double targetFeasible = par(o, "target", 0.2), penInc = 1.2, penDec = 0.85, repairProb = par(o, "repair", 0.5), repairMult = 10;
    const double arc = std::max(init.cost(I), 1e-9) / (I.n + std::max(1, init.nRoutes())), penQ0 = arc / std::max(I.totalDemand() / I.n, 1e-9);
    double penQ = penQ0, penD = 1.0;
    const Preset ps = preset(presetOf(I, o));

    HPop F, X;                                                       // feasible / infeasible
    int nonImp = 0, okQ = 0, okD = 0, seen = 0;
    auto pcOf = [&](const HInd& x) { return x.dist + penQ * x.exQ + penD * x.exD + o.lamFleet * x.fleetEx; };
    auto educate = [&](const Solution& s, double mult) {
        LSState ls(ctx, s, o.lamFleet);
        ls.setSoft(penQ * mult, penD * mult);
        RunOptions ro;
        ro.ops = ps.ops & ~o.excludeOps;
        ro.relocLen = ps.relocLen;
        ro.swapLen = ps.swapLen;
        if (ro.ops) ls.run(ro, rng);
        account(ls);
        HInd x;
        x.sol = ls.toSolution();
        x.dist = x.sol.cost(I);
        x.exQ = ls.excessLoad();
        x.exD = ls.excessDuration();
        x.fleetEx = ls.excess();
        x.tour = giantTour(I, x.sol);
        x.dep = depGenes(I, dord, x.sol);
        return x;
    };
    auto insert = [&](HInd x) {
        x.pc = pcOf(x);
        HPop& pop = x.feasible() ? F : X;
        if (x.feasible() && T.offer(x.sol, x.dist + o.lamFleet * x.fleetEx)) nonImp = 0;
        pop.add(std::move(x));
        if ((int)pop.v.size() >= mu + lambda) pop.survivors(mu, nElite, nClose);
    };
    auto process = [&](const Solution& raw) {
        HInd x = educate(raw, 1);
        okQ += x.exQ < 1e-9;
        okD += x.exD < 1e-9;
        seen++;
        const bool feas = x.feasible();
        Solution keep;
        if (!feas) keep = x.sol;
        insert(std::move(x));
        if (!feas && rng.random() < repairProb) {
            HInd y = educate(keep, repairMult);
            if (y.feasible()) insert(std::move(y));
        }
    };
    auto randomIndividual = [&]() {
        std::vector<int> tour(I.n);
        std::iota(tour.begin(), tour.end(), 0);
        rng.shuffle(tour);
        std::vector<char> dep(I.n, 0);
        if (I.t > 1) for (int c = 0; c < I.n; c++) dep[c] = rng.random() >= 0.75 ? 1 : 0;
        return decodeChrom(I, dord, tour, dep);
    };
    auto pick = [&]() -> const HInd& {                               // binary tournament over both sub-populations
        const int total = (int)F.v.size() + (int)X.v.size();
        auto at = [&](int k) -> const HInd& { return k < (int)F.v.size() ? F.v[k] : X.v[k - F.v.size()]; };
        const HInd& a = at(rng.randrange(total));
        const HInd& b = at(rng.randrange(total));
        return a.bf <= b.bf ? a : b;
    };

    int initLeft = 4 * mu, it;
    for (it = 1; it <= o.iters; it++) {
        T.it = it;
        nonImp++;
        if (it == 1) process(init.nonEmpty());
        else if (initLeft > 0 || F.v.size() + X.v.size() < 2) { process(randomIndividual()); if (initLeft > 0) initLeft--; }
        else {
            F.rank(nElite, nClose);
            X.rank(nElite, nClose);
            const HInd* p1 = &pick();
            const HInd* p2 = &pick();
            for (int t = 0; t < 5 && p2 == p1; t++) p2 = &pick();
            std::vector<int> tour = orderCrossover(p1->tour, p2->tour, rng);
            std::vector<char> dep(I.n);
            for (int c = 0; c < I.n; c++) dep[c] = rng.random() < 0.5 ? p1->dep[c] : p2->dep[c];
            process(decodeChrom(I, dord, tour, dep));
        }
        if (seen >= window) {                                        // penalty adaptation towards the target share of feasible offspring
            const double fq = (double)okQ / seen, fd = (double)okD / seen;
            if (fq < targetFeasible - 0.05) penQ = std::min(penQ * penInc, penQ0 * 1e4); else if (fq > targetFeasible + 0.05) penQ = std::max(penQ * penDec, penQ0 * 1e-2);
            if (fd < targetFeasible - 0.05) penD = std::min(penD * penInc, 1e4); else if (fd > targetFeasible + 0.05) penD = std::max(penD * penDec, 1e-2);
            okQ = okD = seen = 0;
            for (HInd& x : X.v) x.pc = pcOf(x);
            X.stale = true;
        }
        if (nonImp >= divAfter) {                                    // diversification
            F.keepBest(mu / 3);
            X.keepBest(mu / 3);
            initLeft = 4 * mu;
            nonImp = 0;
        }
        if (T.timeUp(o.timeLimit)) break;
    }
    return T.finish(std::min(it, o.iters));
}

// ===========================================================================================
// Shared by grasp-f, grasppr-f and aco-f: savings pairs over the granular neighbourhood and merge-based construction.
// ===========================================================================================
struct SavData {
    struct Pair { int i, j; double s; };
    std::vector<Pair> pairs;                                         // customers of the same depot group, sorted by decreasing saving
    std::vector<int> grp;                                            // depot group of every customer (its nearest depot)
    std::unordered_map<long long, int> index;                        // undirected pair -> position in `pairs`
    long long key(int a, int b, int n) const { return a < b ? (long long)a * n + b : (long long)b * n + a; }
};
SavData buildSavings(const Instance& I, LSContext& ctx, const DepotOrder& dord, int kn = 25) {
    ctx.ensureAdjacency();
    SavData S;
    S.grp.resize(I.n);
    for (int c = 0; c < I.n; c++) S.grp[c] = dord[c][0];
    for (int i = 0; i < I.n; i++) {
        const int dn = I.n + S.grp[i];
        int taken = 0;
        for (int j : ctx.adj[i]) {
            if (j == i) continue;
            if (taken++ >= kn) break;
            if (S.grp[j] != S.grp[i] || S.index.count(S.key(i, j, I.n))) continue;
            S.index[S.key(i, j, I.n)] = 0;
            S.pairs.push_back({i, j, I.d(dn, i) + I.d(dn, j) - I.d(i, j)});
        }
    }
    std::stable_sort(S.pairs.begin(), S.pairs.end(), [](const SavData::Pair& a, const SavData::Pair& b) { return a.s > b.s; });
    for (int k = 0; k < (int)S.pairs.size(); k++) S.index[S.key(S.pairs[k].i, S.pairs[k].j, I.n)] = k;
    return S;
}
struct SavState {                                                    // one route per customer, merged at their end points
    const Instance& I;
    const SavData& S;
    std::vector<int> rid;
    std::vector<std::vector<int>> routes;
    std::vector<double> load, dur;
    SavState(const Instance& inst, const SavData& s) : I(inst), S(s), rid(inst.n), routes(inst.n), load(inst.n), dur(inst.n) {
        for (int c = 0; c < I.n; c++) {
            rid[c] = c;
            routes[c] = {c};
            load[c] = I.demand[c];
            dur[c] = 2 * I.d(I.n + S.grp[c], c) + I.service[c];
        }
    }
    bool valid(const SavData::Pair& p) const {                       // once false, false for the rest of the construction
        const int ri = rid[p.i], rj = rid[p.j], d = S.grp[p.i];
        if (ri == rj) return false;
        if ((routes[ri].front() != p.i && routes[ri].back() != p.i) || (routes[rj].front() != p.j && routes[rj].back() != p.j)) return false;
        if (load[ri] + load[rj] > I.capacity[d] + EPS) return false;
        return dur[ri] + dur[rj] - p.s <= I.maxDuration[d] + EPS;
    }
    void merge(const SavData::Pair& p) {
        const int ri = rid[p.i], rj = rid[p.j];
        std::vector<int>& a = routes[ri];
        std::vector<int>& b = routes[rj];
        if (a.back() != p.i) std::reverse(a.begin(), a.end());
        if (b.front() != p.j) std::reverse(b.begin(), b.end());
        for (int c : b) rid[c] = ri;
        a.insert(a.end(), b.begin(), b.end());
        load[ri] += load[rj];
        dur[ri] += dur[rj] - p.s;
        b.clear();
    }
    Solution solution() const {
        Solution s;
        for (int c = 0; c < I.n; c++) if (!routes[c].empty()) { s.routes.push_back(routes[c]); s.depots.push_back(S.grp[routes[c][0]]); }
        return repairFleetLimit(I, s);
    }
};
// biased-randomised savings (Juan et al. 2010 in spirit): the list is scanned in order of decreasing saving and a feasible merge
// is skipped with probability alpha; further passes without skipping complete the solution
Solution savingsRandomised(const Instance& I, const SavData& S, double alpha, PyRandom& rng) {
    SavState st(I, S);
    std::vector<char> dead(S.pairs.size(), 0);
    for (int pass = 0;; pass++) {
        bool merged = false;
        for (size_t k = 0; k < S.pairs.size(); k++) {
            if (dead[k]) continue;
            if (!st.valid(S.pairs[k])) { dead[k] = 1; continue; }
            if (pass == 0 && rng.random() < alpha) continue;
            st.merge(S.pairs[k]);
            dead[k] = 1;
            merged = true;
        }
        if (pass > 0 && !merged) break;
    }
    return st.solution();
}
// reactive choice of alpha (Prais & Ribeiro 2000): probabilities proportional to (best / mean cost obtained with that alpha)^10
struct ReactiveAlpha {
    std::vector<double> alphas = {0.05, 0.1, 0.2, 0.3, 0.4}, prob, sum, cnt;
    double best = INF;
    int since = 0;
    ReactiveAlpha() : prob(alphas.size(), 1.0), sum(alphas.size(), 0.0), cnt(alphas.size(), 0.0) {}
    int pick(PyRandom& rng) const {
        double tot = 0;
        for (double p : prob) tot += p;
        double x = rng.random() * tot;
        for (int i = 0; i < (int)prob.size(); i++) { x -= prob[i]; if (x <= 0) return i; }
        return (int)prob.size() - 1;
    }
    void report(int i, double cost) {
        sum[i] += cost;
        cnt[i]++;
        best = std::min(best, cost);
        if (++since >= 50) {
            since = 0;
            for (size_t k = 0; k < prob.size(); k++) if (cnt[k] > 0) prob[k] = std::pow(best / (sum[k] / cnt[k]), 10.0);
        }
    }
};

// ===========================================================================================
// ils-f: multi-start iterated local search in the scheme of ILS-RVND (Subramanian et al. 2010; Penna, Subramanian & Ochi 2013).
//   Each start: randomised construction + descent; then perturbation by 2-4 random swap(1,1) / shift(1,0) moves + descent of
//   the touched routes; only improvements of the incumbent of the start are accepted; after n + 5 x (number of routes)
//   consecutive failures a new start is made.  The local search is a randomised variable neighbourhood descent (rvnd below).
//   Not included: the set-partitioning phase of Subramanian et al. (2013), which needs a mixed-integer solver.
//   One iteration = one perturbation + descent, or one new start.
// ===========================================================================================
const std::vector<Kind> KINDS_PERTURB = {K_SWAP, K_RELOCATE};
// Randomised variable neighbourhood descent (Subramanian et al. 2010): the inter-route neighbourhoods shift(1,0), shift(2,0),
// swap(1,1), swap(2,*), 2-opt* (and the depot moves when there are several depots) are tried in random order; a neighbourhood
// that does not improve leaves the list; after an improvement the intra-route search (2-opt, Or-opt) is run and all
// neighbourhoods return to the list.  Each neighbourhood is searched to its own local optimum by the granular engine.
Solution rvnd(const LSContext& ctx, const Solution& s, const MethodOptions& o, PyRandom& rng, const std::vector<int>* dirty = nullptr) {
    struct NB { OpSet ops; int relocLen, swapLen; };
    std::vector<NB> all = {{opBit(OP_RELOCATE), 1, 1}, {opBit(OP_RELOCATE), 2, 1}, {opBit(OP_SWAP), 1, 1}, {opBit(OP_SWAP), 1, 2}, {opBit(OP_TWO_OPT_STAR), 1, 1}};
    if (ctx.inst.t > 1) all.push_back({opBit(OP_ROUTE_DEPOT) | opBit(OP_DEPOT_SWAP), 1, 1});
    LSState ls(ctx, s, o.lamFleet);
    auto run = [&](OpSet ops, int rl, int sl) {
        ops &= ~o.excludeOps;
        if (ops == 0) return false;
        RunOptions ro;
        ro.ops = ops;
        ro.relocLen = rl;
        ro.swapLen = sl;
        ro.dirty = dirty;
        const double before = ls.objective();
        ls.run(ro, rng);
        return ls.objective() < before - EPS;
    };
    run(opBit(OP_TWO_OPT) | opBit(OP_OR_OPT), 1, 1);
    std::vector<NB> list = all;
    while (!list.empty()) {
        const int k = rng.randrange((int)list.size());
        if (run(list[k].ops, list[k].relocLen, list[k].swapLen)) { run(opBit(OP_TWO_OPT) | opBit(OP_OR_OPT), 1, 1); list = all; }
        else list.erase(list.begin() + k);
    }
    account(ls);
    return ls.toSolution();
}
Result runIlsFull(const Instance& I, const MethodOptions& o, const LSContext& ctx, const Solution& init) {
    PyRandom rng(o.seed);
    Tracker T(I, init, o.lamFleet, o.target);
    const DepotOrder dord = buildDepotOrder(I);
    const int maxFails = I.n + (int)(par(o, "failfactor", 5) * std::max(1, init.nRoutes()));
    const int pMin = std::max(1, ipar(o, "pmin", 2)), pRange = std::max(1, ipar(o, "prange", 3));
    const double restartAlpha = par(o, "alpha", 0.3);
    Solution cur = init.nonEmpty();
    double curO = objective(I, cur, o.lamFleet);
    int fails = 0, it;
    for (it = 1; it <= o.iters; it++) {
        T.it = it;
        if (it == 1 || fails >= maxFails) {                           // (re)start
            cur = rvnd(ctx, it == 1 ? cur : graspConstruct(I, dord, rng, restartAlpha), o, rng);
            curO = objective(I, cur, o.lamFleet);
            T.offer(cur, curO);
            fails = 0;
        } else {
            Walk w(ctx, cur);
            std::vector<int> touched = shake(w, pMin + rng.randrange(pRange), rng, KINDS_PERTURB);
            if (touched.empty()) { fails++; continue; }
            Solution cand = rvnd(ctx, w.toSolution(), o, rng, &touched);
            const double oc = objective(I, cand, o.lamFleet);
            if (oc < curO - EPS) { cur = std::move(cand); curO = oc; fails = 0; T.offer(cur, curO); } else fails++;
        }
        if (T.timeUp(o.timeLimit)) break;
    }
    return T.finish(std::min(it, o.iters));
}

// ===========================================================================================
// sa-f: simulated annealing with the cooling scheme of Osman (1993).  Start and final temperature are the largest and the
//   smallest positive cost change among 5000 sampled moves (the latter not below 0.2 % of an average arc); the temperature falls after every move by T <- T / (1 + beta T)
//   (beta chosen so that the final temperature is reached at the last iteration); when 20 n consecutive moves are rejected the
//   temperature is reset to max(T_reset / 2, T at which the best solution was found) (non-monotonic cooling).
//   Moves: the random relocate / Or-opt / swap / 2-opt / 2-opt* / new-route moves of the compact method (Osman scans the
//   1-interchange neighbourhood systematically).  One iteration = one sampled move.
// ===========================================================================================
Result runSaFull(const Instance& I, const MethodOptions& o, const LSContext& ctx, const Solution& init) {
    PyRandom rng(o.seed);
    Walk w(ctx, init);
    Tracker T(I, init, o.lamFleet, o.target);
    Cand m;
    double dmax = 0, dmin = INF;
    for (int k = 0; k < 5000; k++)
        if (w.randomMove(rng, KINDS_NEW, m) && m.delta > 1e-9) { dmax = std::max(dmax, m.delta); dmin = std::min(dmin, m.delta); }
    if (dmax <= 0) { dmax = 1; dmin = 0.01; }
    dmin = std::max(dmin, std::min(dmax, par(o, "tfloor", 0.002) * w.cost / (I.n + init.nRoutes())));   // not below 0.2 % of an average arc
    const double Ts = dmax, Tf = dmin, beta = (Ts - Tf) / (std::max(1, o.iters) * Ts * Tf);
    const int resetAfter = std::max(1, (int)(par(o, "resetfactor", 20) * I.n));
    double temp = Ts, Treset = Ts, Tbest = Tf, cur = w.cost;          // Tbest: temperature at which the best solution was found
    int rejected = 0, it;
    for (it = 1; it <= o.iters; it++) {
        T.it = it;
        bool taken = false;
        if (w.randomMove(rng, KINDS_NEW, m) && (m.delta <= 0 || rng.random() < std::exp(-m.delta / temp))) {
            w.apply(m);
            cur += m.delta;
            taken = true;
            if (cur + o.lamFleet * w.excess() < T.bestO - 1e-7) {
                T.offer(w.toSolution());
                w.cost = w.recomputeCost();
                cur = w.cost;
                Tbest = temp;
            }
        }
        temp = temp / (1 + beta * temp);
        if (taken) rejected = 0;
        else if (++rejected >= resetAfter) { Treset = std::max(Treset / 2, Tbest); temp = Treset; rejected = 0; }
        if (T.timeUp(o.timeLimit)) break;
    }
    return T.finish(std::min(it, o.iters));
}

// ===========================================================================================
// ta-f: backtracking adaptive threshold accepting (Tarantilis, Kiranoudis & Vassiliadis 2002), in spirit.  The threshold is
//   lowered by a constant factor after every window of n moves in which at least one move was accepted, and raised again
//   ("backtracking") after a window without any accepted move, so the search cannot freeze; it never exceeds its start value.
//   The lowering factor is set so that an uninterrupted descent reaches 1 % of the start threshold at the last iteration.
//   Moves as in the compact method.  One iteration = one sampled move.
// ===========================================================================================
Result runTaFull(const Instance& I, const MethodOptions& o, const LSContext& ctx, const Solution& init) {
    PyRandom rng(o.seed);
    Walk w(ctx, init);
    Tracker T(I, init, o.lamFleet, o.target);
    const double th0 = 0.05 * w.cost / (I.n + init.nRoutes());
    const int window = std::max(50, I.n);
    const double lower = std::pow(0.01, std::min(1.0, (double)window / std::max(1, o.iters))), raise = std::pow(lower, -20.0);
    double th = th0, cur = w.cost;
    int accepted = 0, it;
    Cand m;
    for (it = 1; it <= o.iters; it++) {
        T.it = it;
        if (w.randomMove(rng, KINDS_NEW, m) && m.delta <= th) {
            w.apply(m);
            cur += m.delta;
            if (std::fabs(m.delta) > 1e-12) accepted++;
            if (cur + o.lamFleet * w.excess() < T.bestO - 1e-7) {
                T.offer(w.toSolution());
                w.cost = w.recomputeCost();
                cur = w.cost;
            }
        }
        if (it % window == 0) {
            th = accepted > 0 ? th * lower : std::min(th0, th * raise);
            accepted = 0;
        }
        if (T.timeUp(o.timeLimit)) break;
    }
    return T.finish(std::min(it, o.iters));
}

// ===========================================================================================
// rtr-f: record-to-record travel in the scheme of Li, Golden & Wasil (2005).  Uphill phase: the customers are visited in
//   random order; for each one the best move with its granular neighbours (one-point move, two-point exchange, 2-opt within
//   and between routes) is executed if the resulting cost is at most (1 + 0.01) x record.  After 30 sweeps a downhill phase
//   (descent with the same moves) follows; if the record was not improved in 5 consecutive rounds, the best solution is
//   perturbed by max(3, n / 50) random moves.  One iteration = one customer visit.
// ===========================================================================================
Result runRtrFull(const Instance& I, const MethodOptions& o, const LSContext& ctx, const Solution& init) {
    PyRandom rng(o.seed);
    Tracker T(I, init, o.lamFleet, o.target);
    const double dev = 0.01;
    const int sweeps = 30, maxFails = 5, n = I.n;
    const OpSet ops = (opBit(OP_RELOCATE) | opBit(OP_SWAP) | opBit(OP_TWO_OPT) | opBit(OP_TWO_OPT_STAR) | opBit(OP_OR_OPT)) & ~o.excludeOps;
    RunOptions ro;
    ro.ops = ops;
    ro.relocLen = 1;
    ro.swapLen = 1;
    auto st = std::make_unique<LSState>(ctx, init, o.lamFleet);
    double cur = st->objective(), record = cur;
    std::vector<int> order(n);
    std::iota(order.begin(), order.end(), 0);
    int fails = 0, step = 0, it;
    bool improvedRound = false;
    for (it = 1; it <= o.iters; it++) {
        T.it = it;
        if (step % n == 0) rng.shuffle(order);
        const int u = order[step % n];
        Move mv;
        mv.delta = INF;
        for (int v : ctx.nbrs[u]) st->scanPair(u, v, ops, 1, 1, mv);
        if (mv.valid && std::fabs(mv.delta) > 1e-9 && cur + mv.delta <= record * (1 + dev)) {
            st->applyMove(mv);
            cur += mv.delta;
            if (cur < record - EPS) { record = cur = st->objective(); T.offer(st->toSolution(), cur); improvedRound = true; }
        }
        if (++step >= sweeps * n) {                                   // downhill phase, then bookkeeping of the round
            step = 0;
            if (ops) st->run(ro, rng);
            cur = st->objective();
            if (cur < record - EPS) { record = cur; T.offer(st->toSolution(), cur); improvedRound = true; }
            fails = improvedRound ? 0 : fails + 1;
            improvedRound = false;
            if (fails >= maxFails) {
                Walk w(ctx, T.best);
                shake(w, std::max(3, n / 50), rng, KINDS);
                account(*st);
                st = std::make_unique<LSState>(ctx, w.toSolution(), o.lamFleet);
                cur = st->objective();
                fails = 0;
            }
        }
        if (T.timeUp(o.timeLimit)) break;
    }
    account(*st);
    return T.finish(std::min(it, o.iters));
}

// ===========================================================================================
// vns-f: variable neighbourhood search with CROSS-exchange shaking and annealing-type acceptance (Hemmelmayr, Doerner & Hartl
//   2009; Pirkwieser & Raidl 2008).  Shaking k = 1..6: a segment of up to k customers of one route is exchanged with a segment
//   of 0..k customers of a nearby route (0: the segment is moved).  Then descent of the two routes.  A worse solution is
//   accepted with probability exp(-increase / T), T falling linearly from 0.5 % of the average arc cost x n ... to zero.
//   k returns to 1 after an improvement of the incumbent, otherwise k <- k mod 6 + 1.  One iteration = shaking + descent.
// ===========================================================================================
bool shakeCross(const Instance& I, const LSContext& ctx, Solution& s, int k, PyRandom& rng, std::vector<int>& touched) {
    std::vector<int> where(I.n, -1), pos(I.n, 0);
    for (size_t r = 0; r < s.routes.size(); r++) for (size_t p = 0; p < s.routes[r].size(); p++) { where[s.routes[r][p]] = (int)r; pos[s.routes[r][p]] = (int)p; }
    auto ok = [&](int d, const std::vector<int>& nodes) {
        double q = 0, sv = 0;
        for (int c : nodes) { q += I.demand[c]; sv += I.service[c]; }
        return q <= I.capacity[d] + EPS && (nodes.empty() ? 0.0 : tourCost(I, I.n + d, nodes)) + sv <= I.maxDuration[d] + EPS;
    };
    for (int attempt = 0; attempt < 30; attempt++) {
        const int u = rng.randrange(I.n);
        if (ctx.nbrs[u].empty()) return false;
        const int v = ctx.nbrs[u][rng.randrange((int)ctx.nbrs[u].size())], ra = where[u], rb = where[v];
        if (ra == rb) continue;
        const std::vector<int>& a = s.routes[ra];
        const std::vector<int>& b = s.routes[rb];
        const int i = pos[u], j = pos[v];
        const int la = 1 + rng.randrange(std::min(k, (int)a.size() - i)), lb = rng.randrange(std::min(k, (int)b.size() - j) + 1);
        std::vector<int> na(a.begin(), a.begin() + i), nb(b.begin(), b.begin() + j);
        na.insert(na.end(), b.begin() + j, b.begin() + j + lb);
        na.insert(na.end(), a.begin() + i + la, a.end());
        nb.insert(nb.end(), a.begin() + i, a.begin() + i + la);
        nb.insert(nb.end(), b.begin() + j + lb, b.end());
        if (!ok(s.depots[ra], na) || !ok(s.depots[rb], nb)) continue;
        touched = na;
        touched.insert(touched.end(), nb.begin(), nb.end());
        s.routes[ra] = std::move(na);
        s.routes[rb] = std::move(nb);
        s = s.nonEmpty();
        return true;
    }
    return false;
}
Result runVnsFull(const Instance& I, const MethodOptions& o, const LSContext& ctx, const Solution& init) {
    PyRandom rng(o.seed);
    Tracker T(I, init, o.lamFleet, o.target);
    const int kmax = std::max(1, ipar(o, "kmax", 6));
    T.it = 1;                                                        // the first descent belongs to iteration 1
    Solution cur = descend(ctx, init.nonEmpty(), o, rng);
    double curO = objective(I, cur, o.lamFleet);
    T.offer(cur, curO);
    const double T0 = par(o, "t0", 0.5) * cur.cost(I) / (I.n + cur.nRoutes());      // default: half an average arc
    int k = 1, it;
    std::vector<int> touched;
    for (it = 1; it <= o.iters; it++) {
        T.it = it;
        Solution cand = cur;
        if (shakeCross(I, ctx, cand, k, rng, touched)) {
            cand = descend(ctx, cand, o, rng, &touched);
            const double oc = objective(I, cand, o.lamFleet), temp = T0 * (1.0 - (double)(it - 1) / o.iters);
            const bool better = oc < curO - EPS;
            if (better || (temp > 0 && rng.random() < std::exp(-(oc - curO) / temp))) { cur = std::move(cand); curO = oc; T.offer(cur, curO); }
            k = better ? 1 : k % kmax + 1;
        } else k = k % kmax + 1;
        if (T.timeUp(o.timeLimit)) break;
    }
    return T.finish(std::min(it, o.iters));
}

// ===========================================================================================
// grasp-f: reactive GRASP (Prais & Ribeiro 2000) with a randomised Clarke-Wright construction (savings list scanned in order,
//   a feasible merge skipped with probability alpha) instead of the randomised nearest-neighbour construction, then descent.
//   alpha is drawn from {0.05, 0.1, 0.2, 0.3, 0.4} with probabilities updated every 50 iterations.
//   One iteration = one construction + descent.
// ===========================================================================================
Result runGraspFull(const Instance& I, const MethodOptions& o, LSContext& ctx, const Solution& init) {
    PyRandom rng(o.seed);
    Tracker T(I, init, o.lamFleet, o.target);
    const DepotOrder dord = buildDepotOrder(I);
    const SavData S = buildSavings(I, ctx, dord);
    ReactiveAlpha ra;
    int it;
    for (it = 1; it <= o.iters; it++) {
        T.it = it;
        const int a = ra.pick(rng);
        Solution s = descend(ctx, savingsRandomised(I, S, ra.alphas[a], rng), o, rng);
        const double os = objective(I, s, o.lamFleet);
        ra.report(a, os);
        T.offer(s, os);
        if (T.timeUp(o.timeLimit)) break;
    }
    return T.finish(std::min(it, o.iters));
}

// ===========================================================================================
// lns-f: large neighbourhood search of Shaw (1998).  Removal by relatedness (distance, demand, same route; determinism 6);
//   the number of removed customers starts at 2 and grows by one after 100 consecutive failures, up to min(30, n / 2), then
//   starts again; only improving solutions are accepted.  Re-insertion: branch and bound with limited discrepancy (2 discrepancies
//   up to 6 removed customers, 1 up to 12, 0 = farthest insertion above), bounded by the cost of the current solution.
//   One iteration = one removal + re-insertion.
// ===========================================================================================
// Re-insertion of Shaw (1998): branch and bound with limited discrepancy.  The customer whose cheapest insertion is most
// expensive is inserted next (farthest insertion); its insertion points (the best position of every route) are tried in order
// of increasing cost, and taking the i-th best costs i discrepancies; a branch is cut when its cost reaches the best known.
struct Lds {
    const Instance& I;
    double lam, best;
    Solution bestSol;
    bool found = false;
    void rec(Part& P, std::vector<int>& left, int disc, double partial) {
        if (left.empty()) {
            const double obj = partial + lam * excessOf(I, P.cnt);
            if (obj < best - EPS) { best = obj; bestSol = P.toSolution(); found = true; }
            return;
        }
        if (partial >= best - EPS) return;
        int pk = -1;
        double far = -INF;
        std::vector<std::tuple<double, int, int>> opts, cand;
        for (int k = 0; k < (int)left.size(); k++) {
            cand.clear();
            for (int r = 0; r < (int)P.routes.size(); r++) {
                int p;
                double dl;
                if (!P.routes[r].nodes.empty() && bestPos(I, P.routes[r], left[k], p, dl)) cand.emplace_back(dl, r, p);
            }
            std::sort(cand.begin(), cand.end());
            const double cheapest = cand.empty() ? INF : std::get<0>(cand[0]);
            if (cheapest > far) { far = cheapest; pk = k; opts = cand; }
            if (cand.empty()) break;                                 // a customer that fits nowhere is handled first
        }
        const int c = left[pk];
        left.erase(left.begin() + pk);
        if (opts.empty()) {                                          // open a route
            openRoute(P, c);
            rec(P, left, disc, partial + P.routes.back().cost);
            P.cnt[P.routes.back().d]--;
            P.routes.pop_back();
            P.where[c] = -1;
        } else {
            for (int x = 0; x <= disc && x < (int)opts.size(); x++) {
                const int r = std::get<1>(opts[x]);
                insertAt(P, r, std::get<2>(opts[x]), c);
                rec(P, left, disc - x, partial + std::get<0>(opts[x]));
                P.remove(c);
            }
        }
        left.insert(left.begin() + pk, c);
    }
};

Result runLnsFull(const Instance& I, const MethodOptions& o, LSContext& ctx, const Solution& init) {
    ctx.ensureAdjacency();
    PyRandom rng(o.seed);
    Tracker T(I, init, o.lamFleet, o.target);
    double maxD = 0, maxQ = 1e-9;
    for (double d : I.dist) maxD = std::max(maxD, d);
    for (double q : I.demand) maxQ = std::max(maxQ, q);
    const int qMin = std::min(2, I.n - 1), qMax = std::max(qMin, std::min(30, I.n / 2)), growAfter = 100;
    Solution cur = init.nonEmpty();
    double curO = objective(I, cur, o.lamFleet);
    int q = qMin, fails = 0, it;
    for (it = 1; it <= o.iters; it++) {
        T.it = it;
        Part P(I, cur);
        std::vector<int> rm = destroyShawFull(P, q, rng, ctx, maxD, maxQ);
        P.reindex();
        double partial = 0;
        for (const PR& R : P.routes) partial += R.cost;
        Lds lds{I, o.lamFleet, curO, Solution()};
        lds.rec(P, rm, q <= 6 ? 2 : q <= 12 ? 1 : 0, partial);      // only solutions better than the current one are kept
        if (lds.found) { cur = lds.bestSol.nonEmpty(); curO = lds.best; fails = 0; T.offer(cur, curO); }
        else if (++fails >= growAfter) { fails = 0; q = q >= qMax ? qMin : q + 1; }
        if (T.timeUp(o.timeLimit)) break;
    }
    return T.finish(std::min(it, o.iters));
}

// ===========================================================================================
// ig-f: iterated greedy with local search and constant-temperature acceptance (Ruiz & Stuetzle 2007).  Destruction: 10 % of
//   the customers (2..60) at random; construction: cheapest insertion in random order; then descent of the re-inserted
//   customers (the compact method has no local search); a worse solution is accepted with probability
//   exp(-increase / T), T = 0.05 x average arc cost (constant).  One iteration = destruction + construction + descent.
// ===========================================================================================
Result runIgFull(const Instance& I, const MethodOptions& o, LSContext& ctx, const Solution& init) {
    ctx.ensureAdjacency();
    PyRandom rng(o.seed);
    Tracker T(I, init, o.lamFleet, o.target);
    const int q = std::max(std::min(2, I.n - 1), std::min(60, (int)std::floor(par(o, "frac", 0.1) * I.n + 0.5)));
    T.it = 1;                                                        // the first descent belongs to iteration 1
    Solution cur = descend(ctx, init.nonEmpty(), o, rng);
    double curO = objective(I, cur, o.lamFleet);
    T.offer(cur, curO);
    const double temp = par(o, "temp", 0.05) * cur.cost(I) / (I.n + cur.nRoutes());
    int it;
    for (it = 1; it <= o.iters; it++) {
        T.it = it;
        Part P(I, cur);
        std::vector<int> rm = destroyRandom(P, q, rng);
        P.reindex();
        repairGreedy(P, rm, rng);
        P.reindex();
        Solution cand = descend(ctx, P.toSolution(), o, rng, &rm);
        const double oc = objective(I, cand, o.lamFleet);
        if (oc < curO - EPS || rng.random() < std::exp(-(oc - curO) / temp)) { cur = std::move(cand); curO = oc; T.offer(cur, curO); }
        if (T.timeUp(o.timeLimit)) break;
    }
    return T.finish(std::min(it, o.iters));
}

// ===========================================================================================
// aco-f: savings-based ant system (Reimann, Doerner & Hartl 2004, D-Ants without the decomposition) with rank-based pheromone
//   update (Bullnheimer, Hartl & Strauss 1999).  An ant builds its solution by merging routes: among the k = n / 4 best
//   feasible merges, ranked by tau^alpha x saving^beta (alpha = beta = 5), one is drawn with probability proportional to that
//   value; the solution is improved by descent.  After 10 ants: tau <- 0.95 tau; the 5 best ants of the colony deposit
//   (6 - rank) / cost on their arcs and the best solution found so far deposits 6 / cost.
//   One iteration = one ant (construction + descent).
// ===========================================================================================
Result runAcoFull(const Instance& I, const MethodOptions& o, LSContext& ctx, const Solution& init) {
    PyRandom rng(o.seed);
    Tracker T(I, init, o.lamFleet, o.target);
    const DepotOrder dord = buildDepotOrder(I);
    const SavData S = buildSavings(I, ctx, dord);
    const double alpha = 5, beta = 5, rho = 0.95;
    const int ants = 10, sigma = 6, k = std::max(2, I.n / 4), M = (int)S.pairs.size();
    Solution globS = init.nonEmpty();
    double globO = objective(I, globS, o.lamFleet);
    std::vector<double> tau(M, 0.5 * sigma * (sigma + 1) / ((1 - rho) * globO)), weight(M);
    std::vector<int> order(M);
    auto refresh = [&]() {
        for (int p = 0; p < M; p++) weight[p] = std::pow(tau[p], alpha) * std::pow(std::max(S.pairs[p].s, 1e-9), beta);
        std::iota(order.begin(), order.end(), 0);
        std::stable_sort(order.begin(), order.end(), [&](int a, int b) { return weight[a] > weight[b]; });
    };
    refresh();
    auto deposit = [&](const Solution& s, double amount) {
        for (const auto& r : s.routes)
            for (size_t p = 0; p + 1 < r.size(); p++) {
                auto f = S.index.find(S.key(r[p], r[p + 1], I.n));
                if (f != S.index.end()) tau[f->second] += amount;
            }
    };
    auto construct = [&]() {
        SavState st(I, S);
        std::vector<int> live = order, top;
        while (true) {
            top.clear();
            size_t keep = 0;
            double tot = 0;
            for (size_t x = 0; x < live.size(); x++) {
                const int p = live[x];
                if ((int)top.size() < k) {
                    if (!st.valid(S.pairs[p])) continue;                 // dropped for good
                    top.push_back(p);
                    tot += weight[p];
                }
                live[keep++] = p;
            }
            live.resize(keep);
            if (top.empty()) break;
            double x = rng.random() * tot;
            int pick = top.back();
            for (int p : top) { x -= weight[p]; if (x <= 0) { pick = p; break; } }
            st.merge(S.pairs[pick]);
        }
        return st.solution();
    };
    std::vector<std::pair<double, Solution>> colony;
    int it;
    for (it = 1; it <= o.iters; it++) {
        T.it = it;
        Solution s = descend(ctx, construct(), o, rng);
        const double os = objective(I, s, o.lamFleet);
        T.offer(s, os);
        colony.push_back({os, std::move(s)});
        if ((int)colony.size() >= ants) {
            std::stable_sort(colony.begin(), colony.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
            if (colony[0].first < globO) { globO = colony[0].first; globS = colony[0].second; }
            for (double& x : tau) x *= rho;
            for (int r = 0; r < sigma - 1 && r < (int)colony.size(); r++) deposit(colony[r].second, (sigma - 1 - r) / colony[r].first);
            deposit(globS, sigma / globO);
            colony.clear();
            refresh();
        }
        if (T.timeUp(o.timeLimit)) break;
    }
    return T.finish(std::min(it, o.iters));
}

// ===========================================================================================
// grasppr-f: GRASP with evolutionary path relinking (Resende & Ribeiro 2005; Resende & Werneck 2004).  Construction: randomised
//   Clarke-Wright (alpha = 0.2) + descent.  The local optimum is relinked with an elite solution drawn with probability
//   proportional to its distance, in both directions (back-and-forward), at 20 / 40 / 60 / 80 % of each path; the best
//   intermediate solution is improved by descent.  Elite pool of 10: a solution enters if it is better than the worst and
//   differs from every member; it replaces the most similar member that is not better.  One iteration in ten relinks two elite
//   solutions instead of building a new one (evolutionary path relinking).  One iteration = one construction or elite pair.
// ===========================================================================================
Result runGraspPrFull(const Instance& I, const MethodOptions& o, LSContext& ctx, const Solution& init) {
    PyRandom rng(o.seed);
    Tracker T(I, init, o.lamFleet, o.target);
    const DepotOrder dord = buildDepotOrder(I);
    const SavData S = buildSavings(I, ctx, dord);
    const int eliteSize = std::max(2, ipar(o, "elite", 10));
    const double alphaSav = par(o, "alpha", 0.2), pEvo = par(o, "pevo", 0.1);
    const std::vector<double> fracs = {0.2, 0.4, 0.6, 0.8};
    std::vector<Ind> elite = {makeInd(I, dord, init.nonEmpty(), o.lamFleet)};
    auto add = [&](Ind x) {
        int nearest = -1, nd = std::numeric_limits<int>::max(), worst = 0;
        for (int i = 0; i < (int)elite.size(); i++) {
            const int d = brokenPairs(x.tour, elite[i].tour);
            if (d == 0 || std::fabs(elite[i].o - x.o) < 1e-7) return;
            if (elite[i].o >= x.o && d < nd) { nd = d; nearest = i; }
            if (elite[i].o > elite[worst].o) worst = i;
        }
        if ((int)elite.size() < eliteSize) { elite.push_back(std::move(x)); return; }
        if (x.o < elite[worst].o - EPS) elite[nearest >= 0 ? nearest : worst] = std::move(x);
    };
    auto relinkBoth = [&](const std::vector<int>& a, const std::vector<int>& b) {
        Solution bs;
        double bo = INF;
        for (int dir = 0; dir < 2; dir++)
            for (const auto& t : relink(dir ? b : a, dir ? a : b, fracs)) {
                Solution s = decodeTour(I, dord, t);
                const double os = objective(I, s, o.lamFleet);
                if (os < bo) { bo = os; bs = std::move(s); }
            }
        if (bo < INF) {
            Solution s = descend(ctx, bs, o, rng);
            T.offer(s);
            add(makeInd(I, dord, s, o.lamFleet));
        }
    };
    int it;
    for (it = 1; it <= o.iters; it++) {
        T.it = it;
        if (elite.size() >= 3 && rng.random() < pEvo) {
            const int a = rng.randrange((int)elite.size());
            int b = rng.randrange((int)elite.size() - 1);
            if (b >= a) b++;
            const std::vector<int> ta = elite[a].tour, tb = elite[b].tour;
            relinkBoth(ta, tb);
        } else {
            Solution x = descend(ctx, savingsRandomised(I, S, alphaSav, rng), o, rng);
            T.offer(x);
            Ind X = makeInd(I, dord, x, o.lamFleet);
            std::vector<double> wgt;
            double tot = 0;
            for (const Ind& e : elite) { wgt.push_back(1.0 + brokenPairs(X.tour, e.tour)); tot += wgt.back(); }
            double r = rng.random() * tot;
            int y = (int)elite.size() - 1;
            for (int i = 0; i < (int)elite.size(); i++) { r -= wgt[i]; if (r <= 0) { y = i; break; } }
            const std::vector<int> tx = X.tour, ty = elite[y].tour;
            add(std::move(X));
            relinkBoth(tx, ty);
        }
        if (T.timeUp(o.timeLimit)) break;
    }
    return T.finish(std::min(it, o.iters));
}

// ===========================================================================================
// scatter-f: scatter search with a common-arc combination and dynamic reference-set update (Russell & Chiang 2006; Marti,
//   Laguna & Glover 2006).  Combination of A and B: the customers of A that have no arc in common with B are removed and
//   re-inserted by regret-2 insertion (the common structure is kept); then descent.  The child enters the reference set
//   (10 solutions) at once if it is better than the worst member and not a duplicate, and is then combined with all members.
//   When no new pair is left, the 4 worst members are replaced by the most distant of 8 new random-start local optima.
//   One iteration = one combination or one new start.
// ===========================================================================================
Result runScatterFull(const Instance& I, const MethodOptions& o, LSContext& ctx, const Solution& init) {
    ctx.ensureAdjacency();
    PyRandom rng(o.seed);
    Tracker T(I, init, o.lamFleet, o.target);
    const DepotOrder dord = buildDepotOrder(I);
    const int b = std::max(3, ipar(o, "b", 10)), bDiv = std::max(1, std::min(b - 2, ipar(o, "bdiv", 4))), N = I.n + I.t;
    long nextId = 0;
    struct Member { Ind ind; long id; };
    std::vector<Member> ref, fresh;
    std::vector<std::pair<long, long>> queue;
    auto member = [&](const Solution& s) { return Member{makeInd(I, dord, s, o.lamFleet), nextId++}; };
    auto pairsWith = [&](long id) { for (const Member& m : ref) if (m.id != id) queue.push_back({id, m.id}); };
    auto find = [&](long id) -> const Member* { for (const Member& m : ref) if (m.id == id) return &m; return nullptr; };
    auto combine = [&](const Solution& A, const Solution& B) {
        std::unordered_map<long long, char> arcsB;
        auto key = [&](int x, int y) { return x < y ? (long long)x * N + y : (long long)y * N + x; };
        for (size_t r = 0; r < B.routes.size(); r++) {
            int a = I.n + B.depots[r];
            for (int c : B.routes[r]) { arcsB[key(a, c)] = 1; a = c; }
            arcsB[key(a, I.n + B.depots[r])] = 1;
        }
        Part P(I, A);
        std::vector<int> rm;
        for (const PR& R : P.routes) {
            const int dn = I.n + R.d, L = (int)R.nodes.size();
            for (int i = 0; i < L; i++) {
                const int c = R.nodes[i], pv = i ? R.nodes[i - 1] : dn, nx = i + 1 < L ? R.nodes[i + 1] : dn;
                if (!arcsB.count(key(pv, c)) && !arcsB.count(key(c, nx))) rm.push_back(c);
            }
        }
        if (rm.empty()) { std::vector<int> all(I.n); std::iota(all.begin(), all.end(), 0); rng.shuffle(all); rm.assign(all.begin(), all.begin() + std::max(1, I.n / 10)); }
        for (int c : rm) P.remove(c);
        P.reindex();
        repairRegretFull(P, rm, 2, 0.0, rng);
        P.reindex();
        return descend(ctx, P.toSolution(), o, rng, &rm);
    };
    ref.push_back(member(init.nonEmpty()));
    int toStart = b - 1, it;
    for (it = 1; it <= o.iters; it++) {
        T.it = it;
        if (toStart > 0) {                                            // initial reference set / regeneration of the diverse part
            Solution s = randomStart(I, ctx, dord, o, rng);
            T.offer(s);
            fresh.push_back(member(s));
            if (--toStart == 0) {
                const int want = std::min(b - (int)ref.size(), (int)fresh.size());
                for (int k = 0; k < want; k++) {                      // most distant from the current set first
                    int bi = 0, bd = -1;
                    for (int i = 0; i < (int)fresh.size(); i++) {
                        int md = std::numeric_limits<int>::max();
                        for (const Member& m : ref) md = std::min(md, brokenPairs(fresh[i].ind.tour, m.ind.tour));
                        if (md > bd) { bd = md; bi = i; }
                    }
                    ref.push_back(std::move(fresh[bi]));
                    fresh.erase(fresh.begin() + bi);
                    pairsWith(ref.back().id);
                }
                fresh.clear();
                for (int i = (int)queue.size() - 1; i > 0; i--) std::swap(queue[i], queue[rng.randrange(i + 1)]);
            }
        } else if (queue.empty()) {
            std::stable_sort(ref.begin(), ref.end(), [](const Member& x, const Member& y) { return x.ind.o < y.ind.o; });
            if ((int)ref.size() > b - bDiv) ref.resize(b - bDiv);
            toStart = 2 * bDiv;
        } else {
            const std::pair<long, long> pr = queue.back();
            queue.pop_back();
            const Member* A = find(pr.first);
            const Member* B = find(pr.second);
            if (!A || !B) { it--; continue; }                         // a member of the pair has left the set
            if (rng.random() < 0.5) std::swap(A, B);
            Solution s = combine(A->ind.sol, B->ind.sol);
            T.offer(s);
            Member child = member(s);
            int worst = 0;
            bool dup = false;
            for (int i = 0; i < (int)ref.size(); i++) {
                if (std::fabs(ref[i].ind.o - child.ind.o) < 1e-7) dup = true;
                if (ref[i].ind.o > ref[worst].ind.o) worst = i;
            }
            if (!dup && child.ind.o < ref[worst].ind.o - EPS) {
                const long id = child.id;
                ref[worst] = std::move(child);
                pairsWith(id);
            }
        }
        if (T.timeUp(o.timeLimit)) break;
    }
    return T.finish(std::min(it, o.iters));
}

// ===========================================================================================
// Random-key methods with a local-search decoder.  In the compact methods the fitness of a key vector is the cost of the
// decoded solution; here every decoded solution is improved by descent and the fitness is the cost after the descent.
//   brkga-f: multi-parent BRKGA (Andrade et al. 2021): 3 parents (2 elite, 1 non-elite), a gene comes from the parent of rank r
//            with probability proportional to 1 / r; the improved solution is written back to the keys (chromosome adjustment);
//            population 30, elite 20 %, mutants 15 %; restart (best individual kept) after 30 generations without improvement.
//   pso-f:   GLNPSO (Ai & Kachitvichyanukul 2009): the velocity is drawn towards the personal best, the global best, the local
//            best (ring of 5 particles) and the near-neighbour best (fitness-distance ratio per dimension), with
//            c_p = c_g = 0.5, c_l = c_n = 1.5, inertia 0.9 -> 0.4; swarm 20.  Their solution representation SR-2 is not used:
//            the keys are decoded as in the compact method.
//   de-f:    differential evolution with local search (Teoh, Ponnambalam & Kanagaraj 2015 in spirit): DE/rand/1/bin, F drawn
//            from [0.4, 0.9] for every trial vector, CR = 0.9; the improved solution is written back to the trial vector.
//   One iteration = one decoded and improved key vector.
// ===========================================================================================
Result runKeysFull(const Instance& I, const MethodOptions& o, const LSContext& ctx, const Solution& init, const std::string& id) {
    PyRandom rng(o.seed);
    Tracker T(I, init, o.lamFleet, o.target);
    const DepotOrder dord = buildDepotOrder(I);
    const int dim = keyLen(I);
    auto eval = [&](Keys& k, bool writeBack) {
        Solution s = descend(ctx, decodeKeys(I, dord, k), o, rng);
        const double os = objective(I, s, o.lamFleet);
        T.offer(s, os);
        if (writeBack) k = seedKeys(I, dord, s);
        return os;
    };
    struct KI { Keys keys; double o; };
    auto byO = [](const KI& a, const KI& b) { return a.o < b.o; };
    int it = 1;
    if (id == "brkga-f") {
        const int pop = std::max(6, ipar(o, "pop", 30)), ne = std::max(2, std::min(pop - 2, (int)std::floor(par(o, "pe", 0.2) * pop + 0.5))),
                  nm = std::max(0, std::min(pop - ne - 1, (int)std::floor(par(o, "pm", 0.15) * pop + 0.5))), restartAfter = std::max(1, ipar(o, "restart", 30));
        std::vector<KI> cur, nxt;
        bool evolving = false;
        double genBest = INF;
        int stale = 0;
        auto individual = [&](Keys k) { const double f = eval(k, true); return KI{std::move(k), f}; };
        for (; it <= o.iters; it++) {
            T.it = it;
            if (!evolving) {
                cur.push_back(individual(cur.empty() ? seedKeys(I, dord, init) : randomKeys(I, rng)));
                if ((int)cur.size() >= pop) { std::stable_sort(cur.begin(), cur.end(), byO); nxt.assign(cur.begin(), cur.begin() + ne); evolving = true; genBest = cur[0].o; }
            } else {
                const int idx = (int)nxt.size() - ne, nOff = pop - ne;
                Keys k;
                if (idx >= nOff - nm) k = randomKeys(I, rng);
                else {
                    int e1 = rng.randrange(ne), e2 = rng.randrange(ne - 1);
                    if (e2 >= e1) e2++;
                    int par[3] = {e1, e2, ne + rng.randrange(pop - ne)};
                    std::sort(par, par + 3);                          // the population is sorted: a lower index is a better rank
                    k.resize(dim);
                    for (int i = 0; i < dim; i++) { const double u = rng.random() * (11.0 / 6.0); k[i] = cur[par[u < 1.0 ? 0 : u < 1.5 ? 1 : 2]].keys[i]; }
                }
                nxt.push_back(individual(std::move(k)));
                if ((int)nxt.size() >= pop) {
                    cur = std::move(nxt);
                    std::stable_sort(cur.begin(), cur.end(), byO);
                    if (cur[0].o < genBest - EPS) { genBest = cur[0].o; stale = 0; } else stale++;
                    if (stale >= restartAfter) { cur.resize(1); evolving = false; stale = 0; nxt.clear(); }
                    else nxt.assign(cur.begin(), cur.begin() + ne);
                }
            }
            if (T.timeUp(o.timeLimit)) break;
        }
    } else if (id == "pso-f") {
        const int swarm = 20, ring = 2;
        const double w0 = 0.9, w1 = 0.4, cp = 0.5, cg = 0.5, cl = 1.5, cn = 1.5, vmax = 0.5;
        struct Pt { Keys x, v, pb; double po; };
        std::vector<Pt> sw;
        int gb = 0;
        for (; it <= o.iters; it++) {
            T.it = it;
            if ((int)sw.size() < swarm) {
                Pt p;
                p.x = sw.empty() ? seedKeys(I, dord, init) : randomKeys(I, rng);
                p.v.resize(dim);
                for (int i = 0; i < dim; i++) p.v[i] = (rng.random() - 0.5) * vmax;
                p.po = eval(p.x, false);
                p.pb = p.x;
                sw.push_back(std::move(p));
                if (sw.back().po < sw[gb].po) gb = (int)sw.size() - 1;
            } else {
                const int pi = (it - 1) % swarm;
                Pt& p = sw[pi];
                int lb = pi;                                          // local best: ring neighbourhood
                for (int d = -ring; d <= ring; d++) { const int j = ((pi + d) % swarm + swarm) % swarm; if (sw[j].po < sw[lb].po) lb = j; }
                const double w = w0 + (w1 - w0) * std::min(1.0, (double)it / o.iters);
                for (int i = 0; i < dim; i++) {
                    int nb = pi;                                      // near-neighbour best: largest fitness-distance ratio in this dimension
                    double best = 0;
                    for (int j = 0; j < swarm; j++) {
                        if (j == pi) continue;
                        const double fdr = (p.po - sw[j].po) / (std::fabs(sw[j].pb[i] - p.x[i]) + 1e-9);
                        if (fdr > best) { best = fdr; nb = j; }
                    }
                    double v = w * p.v[i] + cp * rng.random() * (p.pb[i] - p.x[i]) + cg * rng.random() * (sw[gb].pb[i] - p.x[i]) +
                               cl * rng.random() * (sw[lb].pb[i] - p.x[i]) + cn * rng.random() * (sw[nb].pb[i] - p.x[i]);
                    v = std::max(-vmax, std::min(vmax, v));
                    p.v[i] = v;
                    const double x = p.x[i] + v;
                    p.x[i] = x - std::floor(x);
                }
                const double f = eval(p.x, false);
                if (f < p.po) { p.po = f; p.pb = p.x; }
                if (f < sw[gb].po) gb = pi;
            }
            if (T.timeUp(o.timeLimit)) break;
        }
    } else {                                                          // de-f
        const int np = std::max(4, ipar(o, "np", 20));
        const double CR = par(o, "cr", 0.9), fMin = par(o, "fmin", 0.4), fSpan = par(o, "fspan", 0.5);
        std::vector<KI> pop;
        int tgt = 0;
        for (; it <= o.iters; it++) {
            T.it = it;
            if ((int)pop.size() < np) {
                Keys x = pop.empty() ? seedKeys(I, dord, init) : randomKeys(I, rng);
                const double f = eval(x, true);
                pop.push_back({std::move(x), f});
            } else {
                const int i = tgt;
                tgt = (tgt + 1) % np;
                int a, b, c;
                do { a = rng.randrange(np); } while (a == i);
                do { b = rng.randrange(np); } while (b == i || b == a);
                do { c = rng.randrange(np); } while (c == i || c == a || c == b);
                const int jr = rng.randrange(dim);
                const double F = fMin + fSpan * rng.random();
                Keys t(dim);
                for (int j = 0; j < dim; j++) {
                    if (rng.random() < CR || j == jr) { const double v = pop[a].keys[j] + F * (pop[b].keys[j] - pop[c].keys[j]); t[j] = v - std::floor(v); }
                    else t[j] = pop[i].keys[j];
                }
                const double f = eval(t, true);
                if (f <= pop[i].o) pop[i] = {std::move(t), f};
            }
            if (T.timeUp(o.timeLimit)) break;
        }
    }
    return T.finish(std::min(it, o.iters));
}
