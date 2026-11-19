// Solution container, objective evaluation and an independent feasibility checker.
#pragma once
#include <string>
#include <vector>

#include "vrplearn/instance.hpp"

namespace vrp {

constexpr double CHECK_EPS = 1e-7;

struct Solution {
    std::vector<std::vector<int>> routes;
    std::vector<int> depots;            // depots[r] = depot index (0..t-1) route r starts/ends at

    Solution nonEmpty() const;
    double routeCost(const Instance& inst, size_t r) const;
    double cost(const Instance& inst) const;
    int nRoutes() const;
    std::vector<int> routesPerDepot(const Instance& inst) const;

    // Independent check (shares no code with the local search): partition, capacity, duration
    // (travel + service) and per-depot fleet.  Empty result = feasible.
    std::vector<std::string> violations(const Instance& inst) const;
    bool feasible(const Instance& inst) const { return violations(inst).empty(); }

    std::string toText(const Instance& inst) const;             // "Route #k (depot d): c1 c2 ..." + "Cost ..."
    static Solution fromText(const std::string& text);          // inverse of toText (1-based ids)
};

// cost of the closed tour depot -> nodes -> depot (0 for an empty route)
double tourCost(const Instance& inst, int depotNode, const std::vector<int>& nodes);

}  // namespace vrp
