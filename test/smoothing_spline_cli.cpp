#include "reco/NaturalSmoothingSpline.h"
#include <iostream>
#include <iomanip>
int main(){int n,m;double lambda;
while(std::cin>>n>>m>>lambda){
    std::vector<double> z(n),y(n);for(int i=0;i<n;++i)std::cin>>z[i]>>y[i];
    rk4::NaturalSmoothingSpline spline;if(!spline.Fit(z,y,lambda))return 2;
    for(int i=0;i<m;++i){double x;std::cin>>x;std::cout<<std::setprecision(17)<<spline.Evaluate(x)<<' ';}
    std::cout<<'\n';
}return 0;}
