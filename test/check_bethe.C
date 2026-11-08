R__ADD_INCLUDE_PATH(src)
R__LOAD_LIBRARY(build/librk4.so)
#include "reco/ParticleId.h"
#include <cstdio>
using namespace rk4;
void check_bethe() {
    std::printf("%8s %12s %12s %12s %12s\n", "Brho", "p", "d", "t", "alpha");
    for (double b : {0.1, 0.2, 0.4, 0.8}) {
        std::printf("%8.2f %12.3f %12.3f %12.3f %12.3f\n", b,
                    ParticleId::ExpectedDedx(b, 1, 1), ParticleId::ExpectedDedx(b, 2, 1),
                    ParticleId::ExpectedDedx(b, 3, 1), ParticleId::ExpectedDedx(b, 4, 2));
    }
}
