#ifndef SPYRALPP_NUMPY_PCG64_H
#define SPYRALPP_NUMPY_PCG64_H
// Integer seed -> NumPy 2.1.3 SeedSequence(pool_size=4) -> PCG64 XSL RR.
// Adapted from numpy/random/bit_generator.pyx and src/pcg64/pcg64.h.
// See notices/numpy-random.txt for MIT attribution and pinned source links.
#include <cstdint>
#include <array>
namespace rk4 {
class NumpyPcg64 {
    using U128 = __uint128_t;
    U128 state_ = 0, increment_ = 0;
    static constexpr U128 Multiplier() {
        return (U128(2549297995355413924ULL) << 64) | 4865540595714422341ULL;
    }
    static uint32_t Hash(uint32_t value, uint32_t& constant) {
        value ^= constant;
        constant *= 0x931e8875u;
        value *= constant;
        return value ^ (value >> 16);
    }
    static uint32_t Mix(uint32_t x, uint32_t y) {
        uint32_t value = 0xca01f9ddu * x - 0x4973f715u * y;
        return value ^ (value >> 16);
    }
public:
    explicit NumpyPcg64(uint64_t seed = 20261026) { Seed(seed); }
    void Seed(uint64_t seed) {
        const uint32_t entropy[2] = {uint32_t(seed), uint32_t(seed >> 32)};
        const int words = entropy[1] ? 2 : 1;
        uint32_t pool[4], constant = 0x43b0d7e5u;
        for (int i=0; i<4; ++i) pool[i] = Hash(i<words ? entropy[i] : 0, constant);
        for (int i=0; i<4; ++i)
            for (int j=0; j<4; ++j)
                if (i != j) pool[j] = Mix(pool[j], Hash(pool[i], constant));
        uint32_t generated[8]; constant = 0x8b51f9ddu;
        for (int i=0; i<8; ++i) {
            uint32_t value = pool[i%4] ^ constant;
            constant *= 0x58f38dedu;
            value *= constant;
            generated[i] = value ^ (value >> 16);
        }
        uint64_t words64[4];
        for (int i=0; i<4; ++i)
            words64[i] = uint64_t(generated[2*i]) | (uint64_t(generated[2*i+1]) << 32);
        const U128 initial = (U128(words64[0]) << 64) | words64[1];
        const U128 sequence = (U128(words64[2]) << 64) | words64[3];
        increment_ = (sequence << 1) | 1;
        state_ = increment_;
        state_ += initial;
        state_ = state_ * Multiplier() + increment_;
    }
    uint64_t Raw() {
        state_ = state_ * Multiplier() + increment_;
        const uint64_t value = uint64_t(state_ >> 64) ^ uint64_t(state_);
        const unsigned rotation = unsigned(state_ >> 122);
        return (value >> rotation) | (value << ((-rotation) & 63));
    }
    double Random() { return (Raw() >> 11) * (1.0 / 9007199254740992.0); }
};
}
#endif
