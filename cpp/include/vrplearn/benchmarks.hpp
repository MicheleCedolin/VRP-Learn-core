// Benchmark catalogue (data/catalogue.tsv) and loader.
//
// Instances are NOT redistributed in this repository (the MDVRP files come from a GPL-2.0
// repository; CVRPLIB sets have their own terms).  scripts/fetch_instances.sh downloads them into
// data/external/.  `bks` holds best-known / optimal values that were checked against a published
// source; entries whose status is not "verified*" must not be quoted as BKS.
#pragma once
#include <optional>
#include <string>
#include <vector>

#include "vrplearn/instance.hpp"

namespace vrp {

struct BenchEntry {
    std::string name, kind, file, url;
    std::optional<double> bks;
    std::string bksStatus, bksSource;
};

// Directory holding data/catalogue.tsv and data/external/ (env VRPLEARN_DATA, else searched upwards from the cwd).
std::string dataDir();
const std::vector<BenchEntry>& catalogue();
const BenchEntry* findEntry(const std::string& name);
std::string pathOf(const std::string& name);                     // <dataDir>/external/<file>
bool available(const std::string& name);
Instance loadBenchmark(const std::string& name);                  // throws if the file has not been fetched
std::optional<double> gapPct(double cost, const std::string& name);

}  // namespace vrp
