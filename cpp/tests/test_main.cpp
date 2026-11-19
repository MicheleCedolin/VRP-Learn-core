#include <cstdio>
#include <cstring>
#include <string>

#include "testing.hpp"

int main(int argc, char** argv) {
    std::string filter = argc > 1 ? argv[1] : "";
    int passed = 0, failed = 0, skipped = 0;
    for (const auto& c : testing::registry()) {
        if (!filter.empty() && c.name.find(filter) == std::string::npos) continue;
        try {
            c.fn();
            ++passed;
            std::printf("[ ok ] %s\n", c.name.c_str());
        } catch (const testing::Skip& s) {
            ++skipped;
            std::printf("[skip] %s (%s)\n", c.name.c_str(), s.what());
        } catch (const std::exception& e) {
            ++failed;
            std::printf("[FAIL] %s: %s\n", c.name.c_str(), e.what());
        }
        std::fflush(stdout);
    }
    std::printf("\n%d passed, %d failed, %d skipped\n", passed, failed, skipped);
    return failed ? 1 : 0;
}
