// Instance model, file readers / writer and a seeded random generator.
//
// Node indexing used everywhere in the library
//   customers : 0 .. n-1
//   depots    : n .. n+t-1     (depot d  <->  node n+d)
// A single-depot CVRP is an instance with t == 1; the MDVRP is t > 1.  Everything downstream
// (construction, local search, metaheuristics, exact DP) is written for the general case.
#pragma once
#include <limits>
#include <optional>
#include <string>
#include <vector>

namespace vrp {

constexpr double INF = std::numeric_limits<double>::infinity();
constexpr int BIG_FLEET = 1000000;

struct Instance {
    std::string name;
    int n = 0;                          // customers
    int t = 1;                          // depots
    std::vector<double> x, y;           // coordinates, size n + t
    std::vector<double> demand;         // size n
    std::vector<double> service;        // size n
    std::vector<double> capacity;       // size t   vehicle capacity per depot
    std::vector<double> maxDuration;    // size t   travel + service limit per route (INF = none)
    std::vector<int> fleet;             // size t   vehicles per depot (BIG_FLEET = unlimited)
    std::vector<double> dist;           // (n+t)^2, row-major
    bool rounded = false;               // TSPLIB nint convention
    std::string source;

    int size() const { return n + t; }
    double d(int i, int j) const { return dist[static_cast<size_t>(i) * (n + t) + j]; }
    const double* row(int i) const { return dist.data() + static_cast<size_t>(i) * (n + t); }
    int depotNode(int dep) const { return n + dep; }
    const char* kind() const { return t == 1 ? "CVRP" : "MDVRP"; }
    double totalDemand() const;
    int minVehicles() const;            // trivial capacity lower bound on the number of routes
    std::string describe() const;
    void buildDistances();              // fills `dist` from x, y and `rounded`
};

Instance readVrplib(const std::string& path, std::optional<bool> rounded = std::nullopt);
Instance readCordeau(const std::string& path);
Instance readInstance(const std::string& path);
void writeCordeau(const Instance& inst, const std::string& path);

struct GenParams {
    int n = 50, t = 1;
    unsigned long long seed = 0;
    double capacity = -1;               // < 0: derived from avgRouteSize
    bool clustered = false;
    int grid = 1000;
    double avgRouteSize = 8.0;
    double fleetSlack = 1.5;
    double maxDuration = INF;
    bool rounded = false;
};
// Random instance in the spirit of Uchoa et al. (2017); reproducible for a given seed.
Instance generate(const GenParams& p);

}  // namespace vrp
