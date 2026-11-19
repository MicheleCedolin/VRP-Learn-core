// Small utilities that make the C++ code bit-for-bit reproducible against the frozen
// reference tables in cpp/tests/reference/:
//   * PyRandom : MT19937 with CPython's seeding (init_by_array) and the exact algorithms of
//                random.random / randrange / randint / uniform / shuffle / choices.
//   * PySetInt : the open-addressing hash table of CPython's `set` for non-negative ints, so that
//                the *iteration order* of small integer sets is the same (the ILS feeds the order of
//                a set into a shuffle).
//   * pysum    : Neumaier-compensated summation, as used by Python >= 3.12 `sum()` for floats.
//   * pymod    : Python's float modulo (result has the sign of the divisor).
// None of this is needed for correctness of the algorithms; it makes results identical to the
// reference tables, which is how the code is verified (cpp/tests/compare_reference.sh).
#pragma once
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace vrp {

class PyRandom {
public:
    explicit PyRandom(uint64_t seed = 0) {
        uint32_t key[2];
        int len = 1;
        key[0] = static_cast<uint32_t>(seed & 0xffffffffULL);
        if (seed >> 32) {
            key[1] = static_cast<uint32_t>(seed >> 32);
            len = 2;
        }
        initByArray(key, len);
    }

    uint32_t genrand() {
        static const uint32_t mag01[2] = {0x0U, 0x9908b0dfU};
        uint32_t y;
        if (mti_ >= N) {
            int kk;
            for (kk = 0; kk < N - M; kk++) {
                y = (mt_[kk] & UPPER) | (mt_[kk + 1] & LOWER);
                mt_[kk] = mt_[kk + M] ^ (y >> 1) ^ mag01[y & 1U];
            }
            for (; kk < N - 1; kk++) {
                y = (mt_[kk] & UPPER) | (mt_[kk + 1] & LOWER);
                mt_[kk] = mt_[kk + (M - N)] ^ (y >> 1) ^ mag01[y & 1U];
            }
            y = (mt_[N - 1] & UPPER) | (mt_[0] & LOWER);
            mt_[N - 1] = mt_[M - 1] ^ (y >> 1) ^ mag01[y & 1U];
            mti_ = 0;
        }
        y = mt_[mti_++];
        y ^= (y >> 11);
        y ^= (y << 7) & 0x9d2c5680U;
        y ^= (y << 15) & 0xefc60000U;
        y ^= (y >> 18);
        return y;
    }

    // random.random(): 53-bit float in [0, 1)
    double random() {
        uint32_t a = genrand() >> 5, b = genrand() >> 6;
        return (a * 67108864.0 + b) * (1.0 / 9007199254740992.0);
    }

    uint32_t getrandbits(int k) { return k == 0 ? 0U : (genrand() >> (32 - k)); }

    // random.Random._randbelow (n > 0)
    uint32_t randbelow(uint32_t n) {
        int k = 0;
        for (uint32_t v = n; v; v >>= 1) ++k;                 // n.bit_length()
        uint32_t r = getrandbits(k);
        while (r >= n) r = getrandbits(k);
        return r;
    }

    int randrange(int n) { return static_cast<int>(randbelow(static_cast<uint32_t>(n))); }
    int randint(int a, int b) { return a + static_cast<int>(randbelow(static_cast<uint32_t>(b - a + 1))); }
    double uniform(double a, double b) { return a + (b - a) * random(); }

    template <class T>
    void shuffle(std::vector<T>& x) {
        for (size_t i = x.size(); i-- > 1;) {
            size_t j = randbelow(static_cast<uint32_t>(i + 1));
            std::swap(x[i], x[j]);
        }
    }

    // random.choices(range(len(w)), weights=w)[0] for integer weights
    int choice(const std::vector<int>& w) {
        std::vector<double> cum;
        double s = 0;
        for (int v : w) { s += v; cum.push_back(s); }
        double total = cum.back() + 0.0;
        double x = random() * total;
        int hi = static_cast<int>(w.size()) - 1, lo = 0;      // bisect_right(cum, x, 0, hi)
        while (lo < hi) {
            int mid = (lo + hi) / 2;
            if (x < cum[mid]) hi = mid; else lo = mid + 1;
        }
        return lo;
    }

    // standard normal (Box-Muller); only used by the instance generator
    double normal() {
        double u1 = 1.0 - random(), u2 = random();
        return std::sqrt(-2.0 * std::log(u1)) * std::cos(6.283185307179586 * u2);
    }

private:
    static constexpr int N = 624, M = 397;
    static constexpr uint32_t UPPER = 0x80000000U, LOWER = 0x7fffffffU;
    uint32_t mt_[N];
    int mti_ = N + 1;

    void initGenrand(uint32_t s) {
        mt_[0] = s;
        for (mti_ = 1; mti_ < N; mti_++)
            mt_[mti_] = 1812433253U * (mt_[mti_ - 1] ^ (mt_[mti_ - 1] >> 30)) + static_cast<uint32_t>(mti_);
    }
    void initByArray(const uint32_t* key, int len) {
        initGenrand(19650218U);
        int i = 1, j = 0;
        for (int k = (N > len ? N : len); k; k--) {
            mt_[i] = (mt_[i] ^ ((mt_[i - 1] ^ (mt_[i - 1] >> 30)) * 1664525U)) + key[j] + static_cast<uint32_t>(j);
            i++; j++;
            if (i >= N) { mt_[0] = mt_[N - 1]; i = 1; }
            if (j >= len) j = 0;
        }
        for (int k = N - 1; k; k--) {
            mt_[i] = (mt_[i] ^ ((mt_[i - 1] ^ (mt_[i - 1] >> 30)) * 1566083941U)) - static_cast<uint32_t>(i);
            i++;
            if (i >= N) { mt_[0] = mt_[N - 1]; i = 1; }
        }
        mt_[0] = 0x80000000U;
        mti_ = N;
    }
};

// ---------------------------------------------------------------------------------------------
class PySetInt {
public:
    PySetInt() : table_(8), mask_(7) {}

    void add(long long key) {
        std::size_t i = static_cast<std::size_t>(key) & mask_, perturb = static_cast<std::size_t>(key);
        Slot* freeslot = nullptr;
        while (true) {
            std::size_t probes = (i + LINEAR_PROBES <= mask_) ? LINEAR_PROBES : 0;
            std::size_t k = i;
            do {
                Slot& e = table_[k];
                if (e.state == 0) {
                    if (freeslot) { place(*freeslot, key); ++used_; return; }
                    place(e, key);
                    ++fill_; ++used_;
                    if (fill_ * 5 >= mask_ * 3) resize(used_ > 50000 ? used_ * 2 : used_ * 4);
                    return;
                }
                if (e.state == 1 && e.key == key) return;
                if (e.state == 2 && !freeslot) freeslot = &e;
                ++k;
            } while (probes--);
            perturb >>= 5;
            i = (i * 5 + 1 + perturb) & mask_;
        }
    }

    void discard(long long key) {
        std::size_t i = static_cast<std::size_t>(key) & mask_, perturb = static_cast<std::size_t>(key);
        while (true) {
            std::size_t probes = (i + LINEAR_PROBES <= mask_) ? LINEAR_PROBES : 0;
            std::size_t k = i;
            do {
                Slot& e = table_[k];
                if (e.state == 0) return;
                if (e.state == 1 && e.key == key) { e.state = 2; --used_; return; }
                ++k;
            } while (probes--);
            perturb >>= 5;
            i = (i * 5 + 1 + perturb) & mask_;
        }
    }

    // set.pop(): first active slot at or after the finger
    long long pop() {
        std::size_t k = finger_ & mask_;
        while (table_[k].state != 1) { if (++k > mask_) k = 0; }
        long long key = table_[k].key;
        table_[k].state = 2; --used_; finger_ = k + 1;
        return key;
    }

    std::size_t size() const { return used_; }
    bool empty() const { return used_ == 0; }

    std::vector<long long> items() const {
        std::vector<long long> out;
        for (const Slot& s : table_) if (s.state == 1) out.push_back(s.key);
        return out;
    }

private:
    struct Slot { long long key = 0; int state = 0; };        // 0 empty, 1 active, 2 dummy
    static constexpr std::size_t LINEAR_PROBES = 9;
    std::vector<Slot> table_;
    std::size_t mask_, fill_ = 0, used_ = 0, finger_ = 0;

    static void place(Slot& s, long long key) { s.key = key; s.state = 1; }

    static void insertClean(std::vector<Slot>& table, std::size_t mask, long long key) {
        std::size_t i = static_cast<std::size_t>(key) & mask, perturb = static_cast<std::size_t>(key);
        while (true) {
            if (table[i].state == 0) { place(table[i], key); return; }
            if (i + LINEAR_PROBES <= mask) {
                for (std::size_t j = 1; j <= LINEAR_PROBES; j++)
                    if (table[i + j].state == 0) { place(table[i + j], key); return; }
            }
            perturb >>= 5;
            i = (i * 5 + 1 + perturb) & mask;
        }
    }

    void resize(std::size_t minused) {
        std::size_t newsize = 8;
        while (newsize <= minused) newsize <<= 1;
        std::vector<Slot> nt(newsize);
        for (const Slot& s : table_) if (s.state == 1) insertClean(nt, newsize - 1, s.key);
        table_.swap(nt);
        mask_ = newsize - 1;
        fill_ = used_;
    }
};

// ---------------------------------------------------------------------------------------------
inline double pysum(const double* a, std::size_t n) {
    double f = 0.0, c = 0.0;
    for (std::size_t k = 0; k < n; k++) {
        double x = a[k], t = f + x;
        if (std::fabs(f) >= std::fabs(x)) c += (f - t) + x; else c += (x - t) + f;
        f = t;
    }
    if (c != 0.0 && std::isfinite(c)) f += c;
    return f;
}
inline double pysum(const std::vector<double>& v) { return pysum(v.data(), v.size()); }

inline double pymod(double x, double y) {
    double m = std::fmod(x, y);
    if (m != 0.0) { if ((y < 0) != (m < 0)) m += y; }
    else m = std::copysign(0.0, y);
    return m;
}

}  // namespace vrp
