// Local-search engine and neighbourhood operators for the CVRP and the MDVRP.
//
// Design
//  * Every route stores prefix sums of demand, service time and edge cost, so a move is evaluated
//    in O(1): a candidate route is a concatenation of "pieces" (route, i, j, reversed) =
//    route.nodes[i..j) (optionally reversed).  The same machinery evaluates relocate, swap /
//    CROSS-exchange, 2-opt, 2-opt*, Or-opt.
//  * The depot of a route is part of the route.  Every edge cost is evaluated with the route's own
//    depot, so one code path handles the single-depot case and the multi-depot case.  The explicit
//    structural differences are:
//      - 2-opt* with reversed heads is only generated for two routes of the SAME depot;
//      - capacity / duration limits are looked up per depot;
//      - fleet limits are per depot and soft (penalty lamFleet per surplus route);
//      - MDVRP-only moves: route re-homing (route_depot) and depot_swap.
//  * Neighbourhoods are granular: only pairs (u, v) with v among the k nearest customers of u are
//    examined (Toth & Vigo 2003).  A queue of "dirty" customers keeps re-optimisation cheap.
//  * First improvement over customers, best improvement over the move variants of a pair.
#pragma once
#include <array>
#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "vrplearn/instance.hpp"
#include "vrplearn/pyrandom.hpp"
#include "vrplearn/solution.hpp"

namespace vrp {

enum Op : int {
    OP_RELOCATE = 0, OP_SWAP, OP_TWO_OPT, OP_TWO_OPT_STAR, OP_OR_OPT, OP_SWAP_STAR,
    OP_ROUTE_DEPOT, OP_DEPOT_SWAP, OP_NEW_ROUTE,                     // OP_NEW_ROUTE: reported, not selectable
    OP_COUNT
};
const char* opName(int op);
std::optional<Op> opFromName(const std::string& name);

using OpSet = unsigned;
constexpr OpSet opBit(Op o) { return 1u << o; }

struct Preset { OpSet ops; int relocLen; int swapLen; };
// "basic": classical operators; "hgs": + SWAP* and segments of length 2; "full": + depot moves
Preset preset(const std::string& name);
OpSet opsFromNames(const std::vector<std::string>& names);

// Instance-level data shared by all local-search runs.
struct LSContext {
    const Instance& inst;
    int k;                                       // neighbourhood size actually used (<= n-1)
    std::vector<std::vector<int>> nbrs;          // k nearest customers (stable order by distance, then index)
    std::vector<std::vector<int>> adj;           // all customers by distance (incl. itself first); lazily filled

    explicit LSContext(const Instance& inst, int k = 20);
    void ensureAdjacency();                      // needed by the string-removal ruin
};

struct RouteData {
    static long nextUid;                         // deterministic identities (never use addresses)
    int d, dn;
    std::vector<int> nodes;
    std::vector<double> cq, cs, ce;              // prefix sums of demand, service, edge cost
    double load = 0, svc = 0, cost = 0;
    long stamp = 0;
    bool dead = false;
    long uid;

    RouteData(int d, int dn, std::vector<int> nodes);
    void refresh(const Instance& inst, long clock);
};

struct Piece { RouteData* R; int i, j; bool rev; };
struct Pieces {
    Piece p[4];
    int k = 0;
    void add(RouteData* R, int i, int j, bool rev) { p[k++] = Piece{R, i, j, rev}; }
};

struct Move {
    Op op = OP_COUNT;
    RouteData *A = nullptr, *B = nullptr;
    Pieces pa, pb;
    double delta = 0;
    bool valid = false;
    int ss[4] = {0, 0, 0, 0};                    // SWAP*: i, j, insertion position in A, in B
};

struct RunOptions {
    OpSet ops = 0;                               // 0 -> hgs preset
    int relocLen = 2;
    int swapLen = 2;
    const std::vector<int>* dirty = nullptr;     // customers to start from (default: all, shuffled)
    bool allowNewRoute = true;
    long maxPasses = 10000;
};

class LSState {
public:
    LSState(const LSContext& ctx, const Solution& sol, double lamFleet = 1e6);

    void run(const RunOptions& opt, PyRandom& rng);
    Solution toSolution() const;
    double cost() const;                         // total distance of the live routes
    int excess() const;                          // total fleet surplus over all depots
    double objective() const { return cost() + lam_ * excess(); }

    std::array<long, OP_COUNT> stats{};          // moves applied per operator
    std::array<double, OP_COUNT> gain{};         // total cost decrease per operator
    // Test hook: if set, every applied move is re-costed from scratch and compared with the
    // predicted delta (fleet penalty excluded); the largest absolute discrepancy is kept here.
    bool verify = false;
    double maxVerifyError = 0;

    // ---- soft-constraint mode and steepest-move interface (used by the granular tabu search tabu-f; off by default) ----
    // In soft mode capacity and duration may be violated: a move delta then includes penQ per unit of excess load and penD per
    // unit of excess duration.  SWAP* and the depot-level passes keep hard constraints in soft mode (they are only applied when
    // the routes they produce are feasible, and their gain is then a lower bound of the penalised gain).
    void setSoft(double penQ, double penD) { soft_ = true; penQ_ = penQ; penD_ = penD; }
    double excessLoad() const;                   // total load above capacity over the live routes
    double excessDuration() const;               // total duration above the limit over the live routes
    std::function<bool(const Move&)> admissible; // optional filter applied to every candidate that would become the best (tabu status)
    // best move between u and v / best new-route move of u; a candidate replaces `best` if its delta is smaller than best.delta
    void scanPair(int u, int v, OpSet ops, int relocLen, int swapLen, Move& best) { bestPair(u, v, ops, relocLen, swapLen, best); }
    void scanNewRoute(int u, int relocLen, Move& best) { newRouteMove(u, relocLen, best); }   // `best` is valid until the next call
    std::vector<RouteData*> applyMove(const Move& m) { return apply(m); }
    const RouteData* routeOf(int c) const { return where_[c]; }
    int posOf(int c) const { return pos_[c]; }

    // exposed for the unit tests of the depot-level moves
    bool routeDepotPass(OpSet ops);

private:
    const LSContext& ctx_;
    const Instance& inst_;
    double lam_;
    bool soft_ = false;
    double penQ_ = 0, penD_ = 0;
    double softPen(int d, double q, double dur) const;
    long clock_ = 0;
    std::vector<std::unique_ptr<RouteData>> routes_;
    std::vector<std::unique_ptr<RouteData>> tmp_;          // candidate fresh routes of new_route moves
    std::vector<RouteData*> where_;
    std::vector<int> pos_, count_;

    void addRoute(std::unique_ptr<RouteData> R);
    void reindex(RouteData* R);
    bool eval(int dn, const Pieces& ps, double& cost, double& q, double& s, int& cnt) const;
    double fleetPen(const RouteData* A, int na, const RouteData* B, int nb) const;
    bool delta2(RouteData* A, const Pieces& pa, RouteData* B, const Pieces& pb, double& out) const;
    bool delta1(RouteData* A, const Pieces& pa, double& out) const;
    static std::vector<int> build(const Pieces& ps);
    std::vector<RouteData*> apply(const Move& m);
    void bestPair(int u, int v, OpSet ops, int relocLen, int swapLen, Move& best);
    void newRouteMove(int u, int relocLen, Move& best);
    void swapStar(RouteData* A, RouteData* B, Move& best) const;
    std::vector<RouteData*> applySwapStar(const Move& m);
};

// Convenience wrapper: run a descent from `sol`.  `ops` (if non-empty) overrides the preset's operators.
struct LSResult { Solution sol; std::array<long, OP_COUNT> stats; std::array<double, OP_COUNT> gain; };
LSResult localSearch(const Instance& inst, const Solution& sol, const std::string& presetName = "hgs",
                     const std::vector<std::string>* ops = nullptr, int k = 20, unsigned long long seed = 0,
                     const LSContext* ctx = nullptr);

}  // namespace vrp
