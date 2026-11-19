#include "vrplearn/instance.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>

#include "vrplearn/pyrandom.hpp"

namespace vrp {

double Instance::totalDemand() const {
    double s = 0;
    for (double q : demand) s += q;
    return s;
}

int Instance::minVehicles() const {
    double q = *std::max_element(capacity.begin(), capacity.end());
    return static_cast<int>(std::ceil(totalDemand() / q));
}

std::string Instance::describe() const {
    char buf[512];
    int fmax = *std::max_element(fleet.begin(), fleet.end());
    double dmax = *std::max_element(maxDuration.begin(), maxDuration.end());
    std::string fl = fmax >= BIG_FLEET ? "unlimited" : std::to_string(fmax);
    char dur[64];
    if (std::isinf(dmax)) std::snprintf(dur, sizeof dur, "none"); else std::snprintf(dur, sizeof dur, "%g", dmax);
    std::snprintf(buf, sizeof buf, "%s: %s, n=%d, depots=%d, Q=%g, fleet/depot=%s, D=%s, LB#routes=%d, rounded=%s",
                  name.c_str(), kind(), n, t, *std::max_element(capacity.begin(), capacity.end()), fl.c_str(), dur,
                  minVehicles(), rounded ? "true" : "false");
    return buf;
}

void Instance::buildDistances() {
    const int N = n + t;
    dist.assign(static_cast<size_t>(N) * N, 0.0);
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++) {
            double dx = x[i] - x[j], dy = y[i] - y[j];
            double v = std::sqrt(dx * dx + dy * dy);
            if (rounded) v = std::floor(v + 0.5);          // TSPLIB / CVRPLIB "nint"
            dist[static_cast<size_t>(i) * N + j] = v;
        }
}

namespace {

std::string slurp(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("cannot open " + path);
    std::ostringstream ss;
    ss << f.rdbuf();
    return ss.str();
}

std::vector<std::string> splitLines(const std::string& s) {
    std::vector<std::string> out;
    std::string cur;
    for (char c : s) {
        if (c == '\n') { out.push_back(cur); cur.clear(); }
        else if (c != '\r') cur += c;
    }
    if (!cur.empty()) out.push_back(cur);
    return out;
}

std::vector<std::string> tokens(const std::string& s) {
    std::vector<std::string> out;
    std::istringstream is(s);
    std::string w;
    while (is >> w) out.push_back(w);
    return out;
}

std::string strip(const std::string& s) {
    size_t a = 0, b = s.size();
    while (a < b && std::isspace(static_cast<unsigned char>(s[a]))) a++;
    while (b > a && std::isspace(static_cast<unsigned char>(s[b - 1]))) b--;
    return s.substr(a, b - a);
}

bool endsWith(const std::string& s, const std::string& suf) {
    return s.size() >= suf.size() && s.compare(s.size() - suf.size(), suf.size(), suf) == 0;
}

// Python str.isupper(): at least one cased character and all cased characters upper-case
bool isUpper(const std::string& s) {
    bool cased = false;
    for (unsigned char c : s) {
        if (std::islower(c)) return false;
        if (std::isupper(c)) cased = true;
    }
    return cased;
}

std::string stem(const std::string& path) {
    size_t slash = path.find_last_of("/\\");
    std::string base = slash == std::string::npos ? path : path.substr(slash + 1);
    size_t dot = base.find_last_of('.');
    return dot == std::string::npos || dot == 0 ? base : base.substr(0, dot);
}

std::string basename(const std::string& path) {
    size_t slash = path.find_last_of("/\\");
    return slash == std::string::npos ? path : path.substr(slash + 1);
}

}  // namespace

Instance readVrplib(const std::string& path, std::optional<bool> rounded) {
    std::map<std::string, std::string> header;
    std::map<int, std::pair<double, double>> coords;
    std::map<int, double> demand;
    std::vector<int> depots;
    std::string section;
    for (const std::string& raw : splitLines(slurp(path))) {
        std::string line = strip(raw);
        if (line.empty() || line == "EOF") continue;
        size_t colon = line.find(':');
        if (colon != std::string::npos && !std::isdigit(static_cast<unsigned char>(line[0])) &&
            isUpper(strip(line.substr(0, colon))) && !endsWith(line, "SECTION")) {
            header[strip(line.substr(0, colon))] = strip(line.substr(colon + 1));
            section.clear();
            continue;
        }
        if (endsWith(line, "SECTION")) { section = line; continue; }
        std::vector<std::string> p = tokens(line);
        if (section == "NODE_COORD_SECTION" && p.size() >= 3)
            coords[std::stoi(p[0])] = {std::stod(p[1]), std::stod(p[2])};
        else if (section == "DEMAND_SECTION" && p.size() >= 2)
            demand[std::stoi(p[0])] = std::stod(p[1]);
        else if (section == "DEPOT_SECTION") {
            int v = std::stoi(p[0]);
            if (v != -1) depots.push_back(v);
        }
    }
    if (depots.empty()) depots.push_back(1);
    int depotId = depots[0];
    Instance I;
    std::vector<int> ids;
    for (auto& kv : coords) if (kv.first != depotId) ids.push_back(kv.first);   // std::map: sorted by id
    I.n = static_cast<int>(ids.size());
    I.t = 1;
    for (int id : ids) { I.x.push_back(coords[id].first); I.y.push_back(coords[id].second); I.demand.push_back(demand[id]); }
    I.x.push_back(coords[depotId].first);
    I.y.push_back(coords[depotId].second);
    // optional route-duration limit and uniform service time (CMT6-10, 13, 14 and several Golden instances):
    // travel time + service times of a route must not exceed DISTANCE
    I.service.assign(I.n, header.count("SERVICE_TIME") ? std::stod(header["SERVICE_TIME"]) : 0.0);
    std::string ewt = header.count("EDGE_WEIGHT_TYPE") ? header["EDGE_WEIGHT_TYPE"] : "EUC_2D";
    I.name = header.count("NAME") ? header["NAME"] : stem(path);
    // Rounding convention of CVRPLIB when the caller does not decide: EUC_2D distances are rounded to the nearest integer
    // (sets A, B, E, F, M, P, X, XL ...), except in the sets that are always used with real distances: CMT, Golden, Li, tai.
    auto startsWith = [&](const char* p) { return I.name.rfind(p, 0) == 0; };
    const bool realSet = startsWith("CMT") || startsWith("Golden") || startsWith("Li_") || startsWith("tai");
    I.rounded = rounded.has_value() ? *rounded : (ewt == "EUC_2D" && !realSet);
    if (!header.count("CAPACITY")) throw std::runtime_error(path + ": no CAPACITY");
    I.capacity = {std::stod(header["CAPACITY"])};
    I.maxDuration = {header.count("DISTANCE") ? std::stod(header["DISTANCE"]) : INF};
    I.fleet = {BIG_FLEET};
    I.source = basename(path);
    I.buildDistances();
    return I;
}

Instance readCordeau(const std::string& path) {
    std::vector<std::vector<std::string>> rows;
    for (const std::string& ln : splitLines(slurp(path))) {
        std::vector<std::string> tk = tokens(ln);
        if (!tk.empty()) rows.push_back(tk);
    }
    if (rows.empty() || rows[0].size() < 4) throw std::runtime_error(path + ": bad Cordeau header");
    int typ = std::stoi(rows[0][0]), m = std::stoi(rows[0][1]);
    Instance I;
    I.n = std::stoi(rows[0][2]);
    I.t = std::stoi(rows[0][3]);
    if (typ != 0 && typ != 2)
        throw std::runtime_error(basename(path) + ": Cordeau type " + std::to_string(typ) + " not supported (only 0=VRP, 2=MDVRP)");
    if (static_cast<int>(rows.size()) < 1 + I.t + I.n + I.t) throw std::runtime_error(path + ": truncated file");
    for (int d = 0; d < I.t; d++) {
        double D = std::stod(rows[1 + d][0]);
        I.maxDuration.push_back(D > 0 ? D : INF);
        I.capacity.push_back(std::stod(rows[1 + d][1]));
        I.fleet.push_back(m);
    }
    const int base = 1 + I.t;
    for (int i = 0; i < I.n; i++) {
        const auto& r = rows[base + i];
        I.x.push_back(std::stod(r[1]));
        I.y.push_back(std::stod(r[2]));
        I.service.push_back(std::stod(r[3]));
        I.demand.push_back(std::stod(r[4]));
    }
    for (int d = 0; d < I.t; d++) {
        const auto& r = rows[base + I.n + d];
        I.x.push_back(std::stod(r[1]));
        I.y.push_back(std::stod(r[2]));
    }
    I.name = stem(path);
    I.source = basename(path);
    I.buildDistances();
    return I;
}

void writeCordeau(const Instance& I, const std::string& path) {
    int fmax = *std::max_element(I.fleet.begin(), I.fleet.end());
    int m = fmax < BIG_FLEET ? fmax : std::max(I.minVehicles() + 2, 2);
    std::ofstream f(path);
    if (!f) throw std::runtime_error("cannot write " + path);
    char buf[256];
    f << (I.t == 1 ? 0 : 2) << ' ' << m << ' ' << I.n << ' ' << I.t << '\n';
    for (int d = 0; d < I.t; d++) {
        std::snprintf(buf, sizeof buf, "%.12g %.12g", std::isinf(I.maxDuration[d]) ? 0.0 : I.maxDuration[d], I.capacity[d]);
        f << buf << '\n';
    }
    for (int i = 0; i < I.n; i++) {
        std::snprintf(buf, sizeof buf, "%d %.12g %.12g %.12g %.12g 1 1 1", i + 1, I.x[i], I.y[i], I.service[i], I.demand[i]);
        f << buf << '\n';
    }
    for (int d = 0; d < I.t; d++) {
        std::snprintf(buf, sizeof buf, "%d %.12g %.12g 0 0 0 0", I.n + d + 1, I.x[I.n + d], I.y[I.n + d]);
        f << buf << '\n';
    }
}

Instance readInstance(const std::string& path) {
    std::string text = slurp(path);
    if (text.find("NODE_COORD_SECTION") != std::string::npos || text.substr(0, 200).find("DIMENSION") != std::string::npos)
        return readVrplib(path);
    return readCordeau(path);
}

Instance generate(const GenParams& p) {
    PyRandom rng(p.seed);
    Instance I;
    I.n = p.n;
    I.t = p.t;
    const double grid = p.grid;
    if (p.clustered) {
        int k = std::max(2, p.n / 12);
        std::vector<std::pair<double, double>> centers;
        for (int c = 0; c < k; c++) centers.push_back({double(rng.randrange(p.grid)), double(rng.randrange(p.grid))});
        for (int i = 0; i < p.n; i++) {
            auto c = centers[rng.randrange(k)];
            I.x.push_back(std::min(grid, std::max(0.0, c.first + rng.normal() * grid * 0.06)));
            I.y.push_back(std::min(grid, std::max(0.0, c.second + rng.normal() * grid * 0.06)));
        }
    } else {
        for (int i = 0; i < p.n; i++) {
            I.x.push_back(rng.randrange(p.grid + 1));
            I.y.push_back(rng.randrange(p.grid + 1));
        }
    }
    if (p.t == 1) {
        I.x.push_back(grid / 2);
        I.y.push_back(grid / 2);
    } else {
        // depots spread on a jittered ring so that each has a catchment area
        double jitter = rng.uniform(0, 0.4);
        for (int d = 0; d < p.t; d++) {
            double ang = 2 * M_PI * (double(d) / p.t) + jitter;
            I.x.push_back(grid / 2 + 0.28 * grid * std::cos(ang));
            I.y.push_back(grid / 2 + 0.28 * grid * std::sin(ang));
        }
    }
    double sum = 0;
    for (int i = 0; i < p.n; i++) { double q = 1 + rng.randrange(10); I.demand.push_back(q); sum += q; }
    double cap = p.capacity > 0 ? p.capacity : std::ceil(sum / p.n * p.avgRouteSize);
    int kmin = static_cast<int>(std::ceil(sum / cap));
    int perDepot = std::max(2, static_cast<int>(std::ceil(p.fleetSlack * kmin / p.t)) + 1);
    I.service.assign(p.n, 0.0);
    I.capacity.assign(p.t, cap);
    I.maxDuration.assign(p.t, p.maxDuration);
    I.fleet.assign(p.t, perDepot);
    I.rounded = p.rounded;
    I.name = "gen-n" + std::to_string(p.n) + "-t" + std::to_string(p.t) + "-s" + std::to_string(p.seed);
    I.source = "vrplearn::generate";
    I.buildDistances();
    return I;
}

}  // namespace vrp
