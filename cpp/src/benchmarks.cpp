#include "vrplearn/benchmarks.hpp"

#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <sys/stat.h>

namespace vrp {

namespace {
bool exists(const std::string& p) {
    struct stat st;
    return stat(p.c_str(), &st) == 0;
}

std::vector<std::string> splitTabs(const std::string& line) {
    std::vector<std::string> out;
    std::string cur;
    for (char c : line) {
        if (c == '\t') { out.push_back(cur); cur.clear(); }
        else cur += c;
    }
    out.push_back(cur);
    return out;
}
}  // namespace

std::string dataDir() {
    if (const char* env = std::getenv("VRPLEARN_DATA")) return env;
    std::string prefix = "";
    for (int up = 0; up < 6; up++) {
        if (exists(prefix + "data/catalogue.tsv")) return prefix + "data";
        prefix += "../";
    }
    throw std::runtime_error("data/catalogue.tsv not found (run from inside the repository or set VRPLEARN_DATA)");
}

const std::vector<BenchEntry>& catalogue() {
    static std::vector<BenchEntry> entries;
    if (entries.empty()) {
        std::ifstream f(dataDir() + "/catalogue.tsv");
        if (!f) throw std::runtime_error("cannot read catalogue.tsv");
        std::string line;
        std::getline(f, line);                                       // header
        while (std::getline(f, line)) {
            if (line.empty()) continue;
            auto c = splitTabs(line);
            if (c.size() < 7) continue;
            BenchEntry e;
            e.name = c[0]; e.kind = c[1]; e.file = c[2]; e.url = c[3];
            if (c[4] != "-") e.bks = std::stod(c[4]);
            e.bksStatus = c[5]; e.bksSource = c[6];
            entries.push_back(e);
        }
    }
    return entries;
}

const BenchEntry* findEntry(const std::string& name) {
    for (const auto& e : catalogue()) if (e.name == name) return &e;
    return nullptr;
}

std::string pathOf(const std::string& name) {
    const BenchEntry* e = findEntry(name);
    if (!e) throw std::invalid_argument("unknown benchmark instance " + name);
    return dataDir() + "/external/" + e->file;
}

bool available(const std::string& name) { return exists(pathOf(name)); }

Instance loadBenchmark(const std::string& name) {
    std::string p = pathOf(name);
    if (!exists(p)) throw std::runtime_error(p + " missing - run scripts/fetch_instances.sh " + name);
    return readInstance(p);
}

std::optional<double> gapPct(double cost, const std::string& name) {
    const BenchEntry* e = findEntry(name);
    if (!e || !e->bks) return std::nullopt;
    return 100.0 * (cost / *e->bks - 1.0);
}

}  // namespace vrp
