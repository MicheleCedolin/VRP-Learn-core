#include "vrplearn/exact.hpp"

#include <algorithm>
#include <stdexcept>

namespace vrp {

namespace {
constexpr double BIG = 1e300;
constexpr double TOL = 1e-9;

struct DepotTables {
    std::vector<double> f;                       // best tour cost per subset (BIG = infeasible)
    std::vector<int> lastOf;                     // best last customer per subset
    std::vector<int> par;                        // Held-Karp parent: par[S*n + last]
    std::vector<double> h;                       // cover with at most m routes
    std::vector<std::vector<int>> choice;        // choice[j][S]: subset T used by the j-th route
    std::vector<std::vector<double>> dp;         // dp[j][S] = exactly j routes
    std::vector<int> bestJ;                      // j achieving h[S]
};
}  // namespace

ExactResult solveExact(const Instance& I) {
    const int n = I.n, t = I.t;
    if (n > 14) throw std::invalid_argument("solveExact: n too large for the subset DP");
    const unsigned FULL = (1u << n) - 1, NS = 1u << n;
    std::vector<double> loadS(NS, 0.0), svcS(NS, 0.0);
    for (unsigned S = 1; S < NS; S++) {
        int low = __builtin_ctz(S);
        loadS[S] = loadS[S & (S - 1)] + I.demand[low];
        svcS[S] = svcS[S & (S - 1)] + I.service[low];
    }

    std::vector<DepotTables> T(t);
    for (int d = 0; d < t; d++) {
        DepotTables& Td = T[d];
        const int dn = n + d;
        std::vector<double> path(static_cast<size_t>(NS) * n, BIG);
        Td.par.assign(static_cast<size_t>(NS) * n, -1);
        for (int v = 0; v < n; v++) path[(size_t(1) << v) * n + v] = I.d(dn, v);
        for (unsigned S = 1; S < NS; S++)
            for (int last = 0; last < n; last++) {
                double cur = path[size_t(S) * n + last];
                if (cur >= BIG || !(S >> last & 1)) continue;
                for (int nx = 0; nx < n; nx++) {
                    if (S >> nx & 1) continue;
                    unsigned S2 = S | (1u << nx);
                    double c = cur + I.d(last, nx);
                    if (c < path[size_t(S2) * n + nx]) { path[size_t(S2) * n + nx] = c; Td.par[size_t(S2) * n + nx] = last; }
                }
            }
        Td.f.assign(NS, BIG);
        Td.lastOf.assign(NS, -1);
        Td.f[0] = 0;
        for (unsigned S = 1; S < NS; S++) {
            if (loadS[S] > I.capacity[d] + TOL) continue;
            double best = BIG;
            int bl = -1;
            for (int last = 0; last < n; last++) {
                if (!(S >> last & 1) || path[size_t(S) * n + last] >= BIG) continue;
                double c = path[size_t(S) * n + last] + I.d(last, dn);
                if (c < best) { best = c; bl = last; }
            }
            if (bl >= 0 && best + svcS[S] <= I.maxDuration[d] + TOL) { Td.f[S] = best; Td.lastOf[S] = bl; }
        }
        const int m = std::min(I.fleet[d], n);
        Td.dp.assign(m + 1, std::vector<double>(NS, BIG));
        Td.choice.assign(m + 1, std::vector<int>(NS, 0));
        Td.dp[0][0] = 0;
        for (int j = 1; j <= m; j++)
            for (unsigned S = 1; S < NS; S++) {
                unsigned lowbit = S & (~S + 1);
                // T ranges over subsets of S that contain the lowest customer of S (canonical order of routes)
                unsigned rest = S ^ lowbit;
                for (unsigned sub = rest;; sub = (sub - 1) & rest) {
                    unsigned Tset = sub | lowbit;
                    double ft = Td.f[Tset];
                    if (ft < BIG) {
                        double prev = Td.dp[j - 1][S ^ Tset];
                        if (prev < BIG && ft + prev < Td.dp[j][S]) { Td.dp[j][S] = ft + prev; Td.choice[j][S] = (int)Tset; }
                    }
                    if (sub == 0) break;
                }
            }
        Td.h.assign(NS, BIG);
        Td.bestJ.assign(NS, 0);
        Td.h[0] = 0;
        for (unsigned S = 1; S < NS; S++)
            for (int j = 1; j <= m; j++)
                if (Td.dp[j][S] < Td.h[S]) { Td.h[S] = Td.dp[j][S]; Td.bestJ[S] = j; }
    }

    // combine depots
    std::vector<std::vector<double>> G(t, std::vector<double>(NS, BIG));
    std::vector<std::vector<unsigned>> pick(t, std::vector<unsigned>(NS, 0));
    for (unsigned S = 0; S < NS; S++) { G[0][S] = T[0].h[S]; pick[0][S] = S; }
    for (int d = 1; d < t; d++)
        for (unsigned S = 0; S < NS; S++)
            for (unsigned sub = S;; sub = (sub - 1) & S) {                 // sub = customers served by depot d
                double hd = T[d].h[sub], gp = G[d - 1][S ^ sub];
                if (hd < BIG && gp < BIG && hd + gp < G[d][S]) { G[d][S] = hd + gp; pick[d][S] = sub; }
                if (sub == 0) break;
            }
    ExactResult res;
    if (G[t - 1][FULL] >= BIG) return res;
    res.feasible = true;
    res.cost = G[t - 1][FULL];

    // reconstruct
    unsigned S = FULL;
    for (int d = t - 1; d >= 0; d--) {
        unsigned mine = d == 0 ? S : pick[d][S];
        unsigned rest = S ^ mine;
        const DepotTables& Td = T[d];
        unsigned U = mine;
        int j = Td.bestJ[U];
        while (U && j > 0) {
            unsigned Tset = (unsigned)Td.choice[j][U];
            std::vector<int> route;
            int cur = Td.lastOf[Tset];
            unsigned Sub = Tset;
            while (cur >= 0) {
                route.push_back(cur);
                int p = Td.par[size_t(Sub) * n + cur];
                Sub &= ~(1u << cur);
                cur = p;
            }
            std::reverse(route.begin(), route.end());
            res.sol.routes.push_back(route);
            res.sol.depots.push_back(d);
            U ^= Tset;
            j--;
        }
        S = rest;
    }
    return res;
}

}  // namespace vrp
