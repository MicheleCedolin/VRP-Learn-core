// Exact solver for TINY instances (n <= ~13), used to certify the heuristics in the test-suite.
//
// It is a set-partitioning DP, not a competitive exact method (branch-cut-and-price handles n~200,
// see Pecin et al. 2017):
//   f_d[S]      cheapest feasible tour of customer set S from depot d  (Held-Karp over subsets,
//               with the capacity Q_d and duration D_d limits of that depot)
//   h_d[S]      cheapest cover of S by at most m_d routes of depot d
//   G_d[S]      cheapest cover of S by the depots 0..d     G_d[S] = min_T G_{d-1}[S \ T] + h_d[T]
// For t = 1 this is the classical CVRP; for t > 1 it is the MDVRP with per-depot fleets, capacities
// and duration limits.  Complexity O(t * (2^n n^2 + n 3^n)).
#pragma once
#include "vrplearn/instance.hpp"
#include "vrplearn/solution.hpp"

namespace vrp {

struct ExactResult {
    bool feasible = false;
    double cost = 0;
    Solution sol;
};

ExactResult solveExact(const Instance& inst);

}  // namespace vrp
