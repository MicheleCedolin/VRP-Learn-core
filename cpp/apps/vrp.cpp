// Command line interface:  vrp --help
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <sys/stat.h>

#include "vrplearn/benchmarks.hpp"
#include "vrplearn/construct.hpp"
#include "vrplearn/exact.hpp"
#include "vrplearn/instance.hpp"
#include "vrplearn/localsearch.hpp"
#include "vrplearn/methods.hpp"
#include "vrplearn/metaheuristics.hpp"

using namespace vrp;

namespace {

void usage() {
    std::puts(
        "usage: vrp <instance> [options]            solve an instance (file, or catalogue name such as p01)\n"
        "       vrp --generate N DEPOTS SEED [opts] solve a random instance\n"
        "       vrp check <instance> <solution>     verify a solution file with the independent checker\n"
        "       vrp info <instance>...              one tab-separated line of statistics per instance\n"
        "options:\n"
        "  --construct savings|sweep|nearest|insertion   starting solution (default savings)\n"
        "  --method ID   none: construction only; ls: descent with --preset; ils; sisr (default);\n"
        "                sa ta rtr tabu vns grasp msls alns lns ig ga aco grasppr scatter brkga pso de rkrs\n"
        "  --preset basic|hgs|full                       operator set of ls / ils (default hgs)\n"
        "  --iters N     iterations (default: sisr 20000, ils 500, sa 2e6, tabu 5e4, vns 300, alns 1e4, ga 400)\n"
        "  --seed S      random seed (default 0)         --k K    neighbourhood size (default 20; 12 for sa..ga)\n"
        "  --time SEC    wall-clock limit of the search  --out F  write the solution to F\n"
        "  --write-instance F   write the (generated) instance in Cordeau format\n"
        "  --exact       solve exactly with the subset DP (n <= 14)");
}

bool fileExists(const std::string& p) {
    struct stat st;
    return stat(p.c_str(), &st) == 0;
}

Instance loadAny(const std::string& arg) {
    if (!fileExists(arg) && findEntry(arg)) return loadBenchmark(arg);
    return readInstance(arg);
}

double since(std::chrono::steady_clock::time_point t0) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

std::string slurp(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("cannot open " + path);
    std::ostringstream ss;
    ss << f.rdbuf();
    return ss.str();
}

int check(const std::string& instArg, const std::string& solPath) {
    Instance inst = loadAny(instArg);
    Solution sol = Solution::fromText(slurp(solPath));
    auto v = sol.violations(inst);
    std::printf("%s: cost %.6f, routes %d, %s\n", inst.name.c_str(), sol.cost(inst), sol.nRoutes(),
                v.empty() ? "FEASIBLE" : "INFEASIBLE");
    for (size_t k = 0; k < v.size() && k < 5; k++) std::printf("  %s\n", v[k].c_str());
    return v.empty() ? 0 : 1;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        std::string instArg, construct = "savings", method = "sisr", presetName = "hgs", outPath, writeInstance;
        int iters = -1, k = 20;
        bool kSet = false;
        unsigned long long seed = 0;
        double timeLimit = 0;
        bool generate = false, exact = false;
        GenParams gp;
        std::vector<std::string> pos;
        for (int a = 1; a < argc; a++) {
            std::string s = argv[a];
            auto next = [&]() -> std::string {
                if (a + 1 >= argc) throw std::runtime_error("missing value after " + s);
                return argv[++a];
            };
            if (s == "-h" || s == "--help") { usage(); return 0; }
            else if (s == "--generate") {
                generate = true;
                gp.n = std::stoi(next()); gp.t = std::stoi(next()); gp.seed = std::stoull(next());
            }
            else if (s == "--construct") construct = next();
            else if (s == "--method") method = next();
            else if (s == "--preset") presetName = next();
            else if (s == "--iters") iters = std::stoi(next());
            else if (s == "--seed") seed = std::stoull(next());
            else if (s == "--k") { k = std::stoi(next()); kSet = true; }
            else if (s == "--time") timeLimit = std::stod(next());
            else if (s == "--out") outPath = next();
            else if (s == "--write-instance") writeInstance = next();
            else if (s == "--exact") exact = true;
            else if (!s.empty() && s[0] == '-') throw std::runtime_error("unknown option " + s);
            else pos.push_back(s);
        }
        if (!pos.empty() && pos[0] == "check") {
            if (pos.size() != 3) { usage(); return 2; }
            return check(pos[1], pos[2]);
        }
        if (!pos.empty() && pos[0] == "info") {
            for (size_t q = 1; q < pos.size(); q++) {
                Instance I = loadAny(pos[q]);
                std::printf("%s\t%s\t%d\t%d\t%.12g\t%d\t%.12g\t%d\n", pos[q].c_str(), I.kind(), I.n, I.t, I.capacity[0], I.fleet[0],
                            I.maxDuration[0], I.minVehicles());
            }
            return 0;
        }
        if (!generate && pos.size() != 1) { usage(); return 2; }
        Instance inst = generate ? vrp::generate(gp) : loadAny(pos[0]);
        if (!writeInstance.empty()) writeCordeau(inst, writeInstance);
        std::puts(inst.describe().c_str());

        auto t0 = std::chrono::steady_clock::now();
        Solution sol;
        if (exact) {
            ExactResult er = solveExact(inst);
            if (!er.feasible) { std::puts("exact: infeasible"); return 1; }
            sol = er.sol;
            method = "exact";
        } else {
            if (construct == "savings") sol = savingsMd(inst);
            else if (construct == "sweep") sol = sweep(inst);
            else if (construct == "nearest") sol = nearestNeighbor(inst);
            else if (construct == "insertion") sol = cheapestInsertion(inst);
            else throw std::runtime_error("unknown construction " + construct);
            std::printf("construction (%s): cost %.2f, routes %d, %.2fs\n", construct.c_str(), sol.cost(inst), sol.nRoutes(), since(t0));
            if (method == "ls") sol = localSearch(inst, sol, presetName, nullptr, k, seed).sol;
            else if (method == "ils") {
                IlsOptions o;
                o.iters = iters < 0 ? 500 : iters; o.seed = seed; o.preset = presetName; o.k = k; o.timeLimit = timeLimit; o.init = &sol;
                sol = ils(inst, o).best;
            } else if (method == "sisr") {
                SisrOptions o;
                o.iters = iters < 0 ? 20000 : iters; o.seed = seed; o.timeLimit = timeLimit; o.init = &sol;
                sol = sisr(inst, o).best;
            } else if (std::find(methodIds().begin(), methodIds().end(), method) != methodIds().end() ||
                       std::find(fullMethodIds().begin(), fullMethodIds().end(), method) != fullMethodIds().end()) {
                MethodOptions o;
                o.iters = iters; o.seed = seed; o.timeLimit = timeLimit; o.init = &sol;
                if (kSet) o.k = k;
                if (presetName != "hgs") o.preset = presetName;
                sol = runMethod(method, inst, o).best;
            } else if (method != "none") throw std::runtime_error("unknown method " + method);
        }
        auto v = sol.violations(inst);
        std::printf("%s: cost %.2f, routes %d, total %.2fs, %s", method.c_str(), sol.cost(inst), sol.nRoutes(), since(t0),
                    v.empty() ? "FEASIBLE" : "VIOLATIONS:");
        for (size_t q = 0; q < v.size() && q < 3; q++) std::printf(" %s;", v[q].c_str());
        std::puts("");
        if (!outPath.empty()) {
            std::ofstream f(outPath);
            f << sol.toText(inst) << '\n';
        }
        return v.empty() ? 0 : 1;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 2;
    }
}
