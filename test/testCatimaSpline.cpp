#include "catima/spline.h"
#include "catima/storage.h"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <stdexcept>
int main() {
    catima::tridiagonal_matrix<5> matrix;
    std::array<double,5> rhs{};
    const std::array<double,5> truth{1,-2,3,-4,5};
    for(int i=0;i<5;++i){
        matrix(i,i)=5+i;
        if(i)matrix(i,i-1)=.2+i;
        if(i<4)matrix(i,i+1)=.7+i;
    }
    const auto& read=matrix;
    for(int i=0;i<5;++i)for(int j=0;j<5;++j)rhs[i]+=read(i,j)*truth[j];
    const auto solution=matrix.trig_solve(rhs);
    for(int i=0;i<5;++i)if(std::abs(solution[i]-truth[i])>1e-12)
        throw std::runtime_error("Catima nonsymmetric tridiagonal solve failed");
    catima::EnergyTable<12> x(-2,2);
    std::vector<double> y;
    for(int i=0;i<12;++i)y.push_back(std::log(x[i])+std::sin(.1*x[i]));
    const catima::cspline_special<decltype(x)> spline(x,y);
    std::cout<<std::setprecision(17);
    for(int i=0;i<=250;++i){
        const double t=std::pow(10.,-2.+4.*i/250.);
        std::cout<<t<<' '<<spline(t)<<' '<<spline.deriv(t)<<'\n';
    }
}
