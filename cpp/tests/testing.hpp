// Minimal self-contained test harness (no external dependencies).
#pragma once
#include <cmath>
#include <cstdio>
#include <functional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace testing {

struct Failure : std::runtime_error { using std::runtime_error::runtime_error; };
struct Skip : std::runtime_error { using std::runtime_error::runtime_error; };

struct Case { std::string name; std::function<void()> fn; };
inline std::vector<Case>& registry() { static std::vector<Case> r; return r; }
struct Registrar { Registrar(const char* n, std::function<void()> f) { registry().push_back({n, std::move(f)}); } };

inline std::string where(const char* file, int line) { return std::string(file) + ":" + std::to_string(line); }

}  // namespace testing

#define TEST(name)                                                       \
    static void test_##name();                                           \
    static testing::Registrar reg_##name(#name, test_##name);            \
    static void test_##name()

#define CHECK(cond)                                                                                         \
    do { if (!(cond)) throw testing::Failure(testing::where(__FILE__, __LINE__) + ": CHECK(" #cond ")"); } while (0)

#define CHECK_MSG(cond, msg)                                                                                \
    do { if (!(cond)) { std::ostringstream os_; os_ << msg;                                                 \
        throw testing::Failure(testing::where(__FILE__, __LINE__) + ": CHECK(" #cond ") " + os_.str()); } } while (0)

#define CHECK_NEAR(a, b, tol)                                                                               \
    do { double a_ = (a), b_ = (b);                                                                         \
        if (!(std::fabs(a_ - b_) <= (tol))) { std::ostringstream os_; os_ << a_ << " vs " << b_;            \
            throw testing::Failure(testing::where(__FILE__, __LINE__) + ": CHECK_NEAR(" #a ", " #b ") " + os_.str()); } } while (0)

#define SKIP(msg) throw testing::Skip(msg)
