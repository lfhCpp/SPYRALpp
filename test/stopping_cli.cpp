#include "sim/Stopping.h"
#include <iostream>
int main() {
    Element h; h.a=1; h.z=1; h.s=2;
    GasTarget gas({h},600);
    StoppingPower power;
    if(!power.Init(12,4,gas))return 1;
    for(double z:{0.,.1,.5,1.})std::cout<<z<<' '<<power.EnergyAfter(153.91616,z)<<std::endl;
}
