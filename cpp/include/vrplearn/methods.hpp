// C++ versions of eighteen methods of the browser Solver (docs/assets/methods.js); SISR and ILS are in metaheuristics.hpp.
//   trajectory        sa (Kirkpatrick et al. 1983; Osman 1993), ta (Dueck & Scheuer 1990), rtr (Dueck 1993),
//                     tabu (Glover 1989; Gendreau, Hertz & Laporte 1994), vns (Mladenovic & Hansen 1997)
//   multi-start       grasp (Feo & Resende 1995), msls (random decodes + local search)
//   destroy & repair  alns (Ropke & Pisinger 2006), lns (Shaw 1998), ig (Ruiz & Stuetzle 2007)
//   population        ga (memetic, Prins 2004), aco (MAX-MIN ant system, Stuetzle & Hoos 2000), grasppr (GRASP + path relinking,
//                     Resende & Ribeiro 2005), scatter (Glover 1998; Marti, Laguna & Glover 2006)
//   random keys       brkga (Goncalves & Resende 2011), pso (Kennedy & Eberhart 1995), de (Storn & Price 1997),
//                     rkrs (random search: the control baseline of the family, not a published method)
// They follow the same schemes and default parameters as the JavaScript versions and work for the CVRP (t = 1) and the MDVRP
// (any customer may be served from any depot).  Random numbers come from PyRandom, the browser engine uses mulberry32, so single
// runs differ between the two; the tests compare feasibility, monotone history and solution quality statistically.
// Teaching-grade: not tuned to reproduce published results.
#pragma once
#include <map>
#include <string>
#include <vector>

#include "vrplearn/metaheuristics.hpp"

namespace vrp {

struct MethodOptions {
    int iters = -1;                              // -1: the default of the method (defaultIters)
    unsigned long long seed = 0;
    double target = 0;                           // > 0: stop as soon as a feasible solution costs <= target (proven optima only)
    double timeLimit = 0;                        // seconds, 0 = none
    int k = 12;                                  // granular neighbourhood of random moves and local search (as in the JS engine)
    std::string preset;                          // operators of the descents: basic | hgs | full; empty = hgs for t = 1, full for t > 1
    double lamFleet = 1e6;                       // penalty per surplus route over a depot's fleet
    const Solution* init = nullptr;              // default: savings
    LSContext* ctx = nullptr;                    // default: built internally with k
    OpSet excludeOps = 0;                        // local-search operators switched off in every descent (ablation studies)
    std::map<std::string, double> params;        // tunable parameters by name (see methodParams); missing names keep their default
};

// Tunable parameters of a method (for automatic configuration, e.g. irace): name, type ('r' real, 'i' integer), default, range,
// whether the range should be sampled on a log scale, and a one-line description.  runMethod rejects unknown names.
// Covers the twenty methods of the study (scripts/hpc/methods_main.txt); sisr is included (its values go to SisrOptions).
struct ParamSpec { std::string name; char type; double def, lo, hi; bool logScale; std::string doc; };
const std::vector<ParamSpec>& methodParams(const std::string& id);
// runMethod fills Result::opStats / opGain with the moves applied and the cost decrease obtained per local-search operator,
// summed over all descents of the run (random moves and removal / insertion operators are not counted there).

const std::vector<std::string>& methodIds();     // the eighteen ids above, in catalogue order
// Upgraded variants that follow the published routing versions more closely (C++ only; the compact methods above stay as they are
// and remain the twins of the browser engine).  Run with runMethod like the others.
//   tabu-f   granular tabu search (Toth & Vigo 2003): steepest move over the granular neighbourhood, arc-based tabu list,
//            capacity / duration violations allowed with self-adjusting penalties (Gendreau, Hertz & Laporte 1994)
//   alns-f   ALNS with the removal / insertion operators, noise, scores and annealing schedule of Pisinger & Ropke (2007)
//   ga-f     hybrid genetic search: feasible / infeasible sub-populations, biased fitness with diversity contribution,
//            adaptive penalties (Vidal et al. 2012; Vidal 2022)
//   ils-f (ILS-RVND scheme), sa-f (Osman cooling), ta-f (backtracking threshold), rtr-f (Li, Golden & Wasil), vns-f (CROSS shaking),
//   grasp-f (reactive, randomised savings), lns-f (Shaw), ig-f (with local search), aco-f (savings-based ants, rank update),
//   grasppr-f (evolutionary path relinking), scatter-f (common-arc combination), brkga-f / pso-f / de-f (local-search decoder):
//   see the comments in src/methods_full.inc for the schemes, sources and parameters.
const std::vector<std::string>& fullMethodIds();
int defaultIters(const std::string& id);
Result runMethod(const std::string& id, const Instance& inst, const MethodOptions& opt = MethodOptions());

// Exposed for the unit tests.
// Optimal split of a giant tour (Prins 2004) for one depot: best partition of `tour` into consecutive routes that respect
// capacity and the duration limit; the order of the tour is kept.
std::vector<std::vector<int>> splitTour(const Instance& inst, int depot, const std::vector<int>& tour);

}  // namespace vrp
