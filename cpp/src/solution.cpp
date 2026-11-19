#include "vrplearn/solution.hpp"

#include <cstdio>
#include <sstream>

#include "vrplearn/pyrandom.hpp"

namespace vrp {

double tourCost(const Instance& inst, int dn, const std::vector<int>& nodes) {
    if (nodes.empty()) return 0.0;
    double c = inst.d(dn, nodes.front()) + inst.d(nodes.back(), dn);
    for (size_t k = 0; k + 1 < nodes.size(); k++) c += inst.d(nodes[k], nodes[k + 1]);
    return c;
}

Solution Solution::nonEmpty() const {
    Solution out;
    for (size_t r = 0; r < routes.size(); r++)
        if (!routes[r].empty()) { out.routes.push_back(routes[r]); out.depots.push_back(depots[r]); }
    return out;
}

double Solution::routeCost(const Instance& inst, size_t r) const {
    return tourCost(inst, inst.depotNode(depots[r]), routes[r]);
}

double Solution::cost(const Instance& inst) const {
    std::vector<double> c(routes.size());
    for (size_t r = 0; r < routes.size(); r++) c[r] = routeCost(inst, r);
    return pysum(c);
}

int Solution::nRoutes() const {
    int k = 0;
    for (const auto& r : routes) k += !r.empty();
    return k;
}

std::vector<int> Solution::routesPerDepot(const Instance& inst) const {
    std::vector<int> cnt(inst.t, 0);
    for (size_t r = 0; r < routes.size(); r++)
        if (!routes[r].empty() && depots[r] >= 0 && depots[r] < inst.t) cnt[depots[r]]++;   // bad indices are reported by violations()
    return cnt;
}

std::vector<std::string> Solution::violations(const Instance& inst) const {
    std::vector<std::string> out;
    char buf[256];
    std::vector<int> seen(inst.n, 0);
    for (size_t r = 0; r < routes.size(); r++) {
        int d = depots[r];
        if (d < 0 || d >= inst.t) {
            std::snprintf(buf, sizeof buf, "route %zu: bad depot %d", r, d);
            out.push_back(buf);
            continue;
        }
        double load = 0, svc = 0;
        bool badNode = false;
        for (int v : routes[r]) {
            if (v < 0 || v >= inst.n) {
                badNode = true;
                std::snprintf(buf, sizeof buf, "route %zu: node %d is not a customer", r, v);
                out.push_back(buf);
            } else {
                seen[v]++;
                load += inst.demand[v];
                svc += inst.service[v];
            }
        }
        if (load > inst.capacity[d] + CHECK_EPS) {
            std::snprintf(buf, sizeof buf, "route %zu (depot %d): load %g > Q %g", r, d, load, inst.capacity[d]);
            out.push_back(buf);
        }
        if (badNode) continue;                  // distances are undefined for non-customer nodes
        double dur = routeCost(inst, r) + svc;
        if (dur > inst.maxDuration[d] + CHECK_EPS) {
            std::snprintf(buf, sizeof buf, "route %zu (depot %d): duration %g > D %g", r, d, dur, inst.maxDuration[d]);
            out.push_back(buf);
        }
    }
    for (int v = 0; v < inst.n; v++)
        if (seen[v] != 1) {
            std::snprintf(buf, sizeof buf, "customer %d visited %d times", v, seen[v]);
            out.push_back(buf);
        }
    std::vector<int> cnt = routesPerDepot(inst);
    for (int d = 0; d < inst.t; d++)
        if (cnt[d] > inst.fleet[d]) {
            std::snprintf(buf, sizeof buf, "depot %d: %d routes > fleet %d", d, cnt[d], inst.fleet[d]);
            out.push_back(buf);
        }
    return out;
}

std::string Solution::toText(const Instance& inst) const {
    std::ostringstream os;
    int k = 0;
    for (size_t r = 0; r < routes.size(); r++) {
        ++k;
        if (routes[r].empty()) continue;
        os << "Route #" << k << " (depot " << depots[r] + 1 << "):";
        for (int v : routes[r]) os << ' ' << v + 1;
        os << '\n';
    }
    char buf[64];
    std::snprintf(buf, sizeof buf, "Cost %.4f", cost(inst));
    os << buf;
    return os.str();
}

Solution Solution::fromText(const std::string& text) {
    Solution s;
    std::istringstream in(text);
    std::string line;
    while (std::getline(in, line)) {
        if (line.rfind("Route #", 0) != 0) continue;
        size_t p = line.find("(depot ");
        size_t q = line.find("):");
        int d = 0;
        size_t from;
        if (p != std::string::npos && q != std::string::npos) {
            d = std::stoi(line.substr(p + 7, q - p - 7)) - 1;
            from = q + 2;
        } else {                                                   // CVRPLIB / HGS-CVRP style "Route #k: ids" (single depot)
            const size_t c = line.find(':');
            if (c == std::string::npos) continue;
            from = c + 1;
        }
        std::istringstream ns(line.substr(from));
        std::vector<int> nodes;
        int v;
        while (ns >> v) nodes.push_back(v - 1);
        s.routes.push_back(nodes);
        s.depots.push_back(d);
    }
    return s;
}

}  // namespace vrp
