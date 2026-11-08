#include "reco/NumpyPcg64.h"
#include <iostream>
#include <iomanip>
int main() {
    uint64_t seed; int n;
    std::cout << std::setprecision(17);
    while (std::cin >> seed >> n) {
        rk4::NumpyPcg64 raw(seed), uniform(seed);
        for (int i=0;i<n;++i) std::cout << raw.Raw() << ' ' << uniform.Random() << '\n';
    }
}
