R__ADD_INCLUDE_PATH(src)
R__LOAD_LIBRARY(build/librk4.so)
#include "reco/Estimator.h"
#include "TSystem.h"
#include <cmath>
#include <limits>
#include <cstdio>
void testReinschDirections() {
    rk4::Cluster forward,backward;
    forward.direction=rk4::Direction::Forward;
    backward.direction=rk4::Direction::Backward;
    for (int i=0;i<150;++i) {
        double t=0.009*i;
        rk4::CloudPoint p;
        p.x=80*std::sin(t);p.y=80*(1-std::cos(t));p.z=100+90*t;p.integral=1000;
        forward.points.push_back(p);p.z=800-90*t;backward.points.push_back(p);
    }
    rk4::Estimator e; e.SetUseReinsch(true);rk4::EstimateResult f,b;
    bool ok=e.Estimate(forward,2.85,f) && e.Estimate(backward,2.85,b);
    ok=ok && f.ok && b.ok && f.polar<M_PI/2 && b.polar>M_PI/2;
    ok=ok && std::fabs(f.brho/b.brho-1)<1e-7 && std::fabs(f.dEdx/b.dEdx-1)<1e-7;
    ok=ok && std::fabs(f.radius-80)/80<0.03;
    rk4::EstimateResult invalid;
    ok=ok && !e.Estimate(forward,0,invalid) && !invalid.ok;
    forward.points[10].z=std::numeric_limits<double>::quiet_NaN();
    ok=ok && !e.Estimate(forward,2.85,invalid) && !invalid.ok;
    std::printf("Reinsch forward/backward symmetry and invalid input: %s\n",ok?"PASS":"FAIL");
    gSystem->Exit(ok?0:1);
}
