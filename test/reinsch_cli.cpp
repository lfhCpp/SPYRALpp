#include "reco/Reinsch.h"
#include <iostream>
#include <iomanip>
int main(){size_t n; double lam; while(std::cin>>n>>lam){std::vector<double>x(n),y(n),f;for(size_t i=0;i<n;++i) std::cin>>x[i]>>y[i]; if(!rk4::ReinschSmooth(x,y,lam,f)) return 2;for(double v:f)std::cout<<std::setprecision(17)<<v<<' '; std::cout<<'\n';}}
