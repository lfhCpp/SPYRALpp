#ifndef RK4_RECO_NATURAL_SMOOTHING_SPLINE_H
#define RK4_RECO_NATURAL_SMOOTHING_SPLINE_H
#include "Reinsch.h"
#include <algorithm>
#include <stdexcept>
#include <vector>
namespace rk4 {
// min sum_i (f(z_i)-y_i)^2 + lambda * integral(f''(z)^2 dz).
// If z_new=a*z, use lambda_new=a^3*lambda to preserve the same curve.
class NaturalSmoothingSpline {
public:
    bool Fit(const std::vector<double>& z,const std::vector<double>& y,double lambda){
        valid_=false;
        if(!ReinschSmooth(z,y,lambda,value_,&second_))return false;
        z_=z;valid_=true;return true;
    }
    const std::vector<double>& Values() const{return value_;}
    double Evaluate(double t) const{
        if(!valid_)throw std::logic_error("Smoothing spline is not fitted");
        const auto end=std::upper_bound(z_.begin(),z_.end(),t);
        const size_t i=end==z_.begin()?0:std::min<size_t>(end-z_.begin()-1,z_.size()-2);
        const double h=z_[i+1]-z_[i],a=(z_[i+1]-t)/h,b=(t-z_[i])/h;
        return a*value_[i]+b*value_[i+1]+((a*a*a-a)*second_[i]+(b*b*b-b)*second_[i+1])*h*h/6;
    }
private:
    bool valid_=false;
    std::vector<double> z_,value_,second_;
};
}
#endif
