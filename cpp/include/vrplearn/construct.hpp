// Constructive heuristics (single and multi depot).
//   savings        Clarke & Wright (1964), per depot after nearest-depot assignment
//   sweep          Gillett & Miller (1974); with several depots: assignment + sweep per depot
//   nearestNeighbor sequential nearest neighbour
//   cheapestInsertion parallel cheapest insertion across all routes of all depots
// All of them finish with a fleet repair (dissolve routes of over-used depots).
#pragma once
#include <optional>
#include <vector>

#include "vrplearn/instance.hpp"
#include "vrplearn/solution.hpp"

namespace vrp {

std::vector<int> assignNearestDepot(const Instance& inst);                  // customer -> depot (ties: lowest index)
Solution savings(const Instance& inst, double lambda = 1.0);                // = savingsMd for t == 1
Solution savingsMd(const Instance& inst, double lambda = 1.0, const std::vector<int>* assignment = nullptr);
Solution sweep(const Instance& inst, int nStarts = 12);
Solution nearestNeighbor(const Instance& inst);
Solution cheapestInsertion(const Instance& inst, std::optional<unsigned long long> shuffleSeed = std::nullopt);
// Dissolves routes of depots that use more vehicles than their fleet (cheapest feasible re-insertion of the customers); no-op if the fleet is respected.
Solution repairFleetLimit(const Instance& inst, Solution sol);

}  // namespace vrp
