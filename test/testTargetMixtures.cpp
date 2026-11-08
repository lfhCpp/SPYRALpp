#include "sim/TargetDensity.h"
#include "sim/Stopping.h"
#include "sim/ParticleDynamics.h"
#include "core/MassData.h"
#include <cmath>
#include <limits>
#include <stdexcept>
#include <iostream>
static void Check(bool good,const char* text){if(!good)throw std::runtime_error(text);}
int main(){
    const std::vector<Element> h2{{1,1,2}},he{{4,2,1}};
    GasMixTarget hydrogen({h2},{1.},600.);
    const auto before=hydrogen.compound();const double before_density=hydrogen.density();
    GasMixTarget helium({he},{1.},760.);
    Check(hydrogen.comps().size()==1 && hydrogen.compound().size()==1,"Gas instances share composition");
    Check(hydrogen.density()==before_density && hydrogen.compound()[0].s==before[0].s,"Second target changed first target");
    GasMixTarget mix({h2,he},{.999,.001},600.);
    const auto atoms=mix.compound();
    Check(atoms.size()==2,"Trace gas disappeared");
    Check(atoms[0].a==1 && atoms[0].z==1 && std::abs(atoms[0].s-1.998)<1e-14,"Hydrogen mole count incorrect");
    Check(atoms[1].a==4 && atoms[1].z==2 && std::abs(atoms[1].s-.001)<1e-14,"Helium isotope/count incorrect");
    const double mass=.999*2*MassData::GetNuclearMass(1,1)+.001*MassData::GetNuclearMass(4,2);
    Check(std::abs(mix.molarMass()-mass)<1e-13,"Molecular mass is not mole-fraction weighted");
    Check(std::abs(mix.density()-600.*mass/(Constants::gasConstant*Constants::roomTemp))<1e-15,"Ideal gas density incorrect");
    Check(std::abs(MassData::CompoundAvgMass(atoms)-mass/1.999)<1e-13,"Fractional counts truncated in average atomic mass");
    ParticleDynamics dynamics;dynamics.ConfigureTarget(atoms,mix.density());
    const auto material=dynamics.GetTargetMaterial();
    Check(material.ncomponents()==2,"Dynamics ignored a component");
    const double expected_he=4*.001/(1.998+4*.001);
    Check(std::abs(material.weight_fraction(1)-expected_he)<1e-13,"Catima interpreted mole fraction as mass fraction");
    StoppingPower stopping;Check(stopping.Init(12,4,mix),"Mixed target stopping initialization");
    const double exit=stopping.EnergyAfter(150.,1.);
    Check(exit>0 && exit<150,"Mixed target produced invalid energy loss");
    bool threw=false;
    try{GasMixTarget invalid({h2,he},{std::numeric_limits<double>::quiet_NaN(),.1},600.);}catch(const TargetError&){threw=true;}
    Check(threw,"NaN gas fraction accepted");
    threw=false;
    try{GasTarget invalid(h2,std::numeric_limits<double>::infinity());}catch(const TargetError&){threw=true;}
    Check(threw,"Infinite gas pressure accepted");
    std::cout<<"Gas composition, independence, trace fraction and Catima conversion passed\n";
}
