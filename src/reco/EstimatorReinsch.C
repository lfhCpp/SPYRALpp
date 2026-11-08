#include "Estimator.h"
#include "NaturalSmoothingSpline.h"
#include "Geometry.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace rk4 {
namespace {
constexpr double Pi=3.1415926535897932384626433832795;
struct Point { double x,y,z,q; };

using Spline=NaturalSmoothingSpline;

// Spyral geometry.circle.least_squares_circle: centered 2x2 normal equations,
// then the mean distance to the fitted center (not the RMS distance).
bool Circle(const std::vector<Point>& points,size_t count,double& cx,double& cy,double& radius) {
    std::vector<double> x,y;x.reserve(count);y.reserve(count);
    for(size_t i=0;i<count;++i){x.push_back(points[i].x);y.push_back(points[i].y);}
    return FitCircleSpyral(x,y,cx,cy,radius);
}

// Python selects the nearest member of linspace(0,2*pi,100000).
// Evaluate the neighboring grid angles directly instead of sorting 100000
// distances. Preserve the sampled vertex: a continuous projection changes
// vertex cuts and polar angles for short arcs.
void Vertex(double cx,double cy,double radius,double& vx,double& vy) {
    const double step=2*Pi/99999;
    double theta=std::atan2(-cy,-cx);if(theta<0)theta+=2*Pi;
    const long center=static_cast<long>(std::floor(theta/step));
    double best=std::numeric_limits<double>::infinity();
    long best_index=100000;
    for(long candidate=center-1;candidate<=center+2;++candidate){
        if(candidate<0||candidate>99999)continue;
        const double angle=candidate==99999?2*Pi:candidate*step;
        const double x=cx+std::cos(angle)*radius,y=cy+std::sin(angle)*radius;
        const double distance=x*x+y*y;
        if(distance<best||(distance==best&&candidate<best_index)){
            best=distance;best_index=candidate;vx=x;vy=y;
        }
    }
}

enum class Pass { Accepted,Rejected,WrongDirection };
Pass EstimatePass(const std::vector<Point>& smooth,bool backward,const Spline& sx,const Spline& sy,
                  double field,double beam_radius,double pad_radius,EstimateResult& result) {
    std::vector<Point> d=smooth;
    if(backward)std::reverse(d.begin(),d.end());
    const size_t n=d.size();
    double vx=d[0].x,vy=d[0].y;
    double farthest=-1;size_t initial=0;
    for(size_t i=1;i<n;++i){
        const double rho=std::hypot(d[i].x-vx,d[i].y-vy);
        if(rho>farthest){farthest=rho;initial=i;}
    }
    // Upstream argmax is over data[1:], then [:maximum+1]. Preserve that slice.
    double cx=0,cy=0,radius=0;
    if(!Circle(d,initial,cx,cy,radius))return Pass::Rejected;
    Vertex(cx,cy,radius,vx,vy);
    std::vector<double> rho(n);
    size_t maximum=0;farthest=-1;
    for(size_t i=0;i<n;++i){
        rho[i]=std::hypot(d[i].x-vx,d[i].y-vy);
        if(rho[i]>farthest){farthest=rho[i];maximum=i;}
    }
    const size_t regression=std::min(n,std::max<size_t>(10,maximum/2));
    double mz=0,mr=0;
    for(size_t i=0;i<regression;++i){mz+=d[i].z;mr+=rho[i];}
    mz/=regression;mr/=regression;
    double zz=0,zr=0;
    for(size_t i=0;i<regression;++i){
        const double z=d[i].z-mz;zz+=z*z;zr+=z*(rho[i]-mr);
    }
    if(zz==0)return Pass::Rejected;
    const double slope=zr/zz;
    if(slope==0||!std::isfinite(slope))return Pass::Rejected;
    if(std::hypot(vx,vy)>beam_radius)return Pass::Rejected;
    double polar=std::atan(slope);
    if((polar>0&&backward)||(polar<0&&!backward))return Pass::WrongDirection;
    if(backward)polar+=Pi;
    double charge=d[0].q;size_t cutoff=maximum;
    for(size_t i=0;i<maximum;++i){
        if(std::hypot(d[i+1].x,d[i+1].y)>pad_radius){cutoff=i+1;break;}
        charge+=d[i+1].q;
    }
    if(charge==d[0].q)return Pass::Rejected;
    // np.linspace(...,1000) has 999 line segments.
    const double begin=d[0].z,end=d[cutoff].z;
    double px=sx.Evaluate(begin),py=sy.Evaluate(begin),pz=begin,length=0;
    for(int i=1;i<1000;++i){
        const double z=i==999?end:begin+(end-begin)*i/999;
        const double x=sx.Evaluate(z),y=sy.Evaluate(z);
        length+=std::sqrt((x-px)*(x-px)+(y-py)*(y-py)+(z-pz)*(z-pz));
        px=x;py=y;pz=z;
    }
    if(!(length>0)||!std::isfinite(length))return Pass::Rejected;
    EstimateResult out;
    out.vertex_x=vx;out.vertex_y=vy;out.vertex_z=-(mr-slope*mz)/slope;
    out.center_x=cx;out.center_y=cy;out.radius=radius;
    out.polar=polar;
    out.azimuthal=std::atan2(vy-cy,vx-cx);
    if(out.azimuthal<0)out.azimuthal+=2*Pi;
    out.azimuthal+=0.5*Pi;if(out.azimuthal>2*Pi)out.azimuthal-=2*Pi;
    out.brho=field*radius*0.001/std::abs(std::sin(polar));
    out.charge=charge;out.arclength=length;out.dEdx=charge/length;
    out.direction=static_cast<int>(backward?Direction::Backward:Direction::Forward);
    out.ok=std::isfinite(out.brho)&&out.brho>0&&std::isfinite(out.dEdx)&&out.dEdx>0;
    if(!out.ok)return Pass::Rejected;
    result=out;return Pass::Accepted;
}
} // namespace

bool Estimator::EstimateReinsch(const Cluster& cluster,double field,EstimateResult& out) const {
    out=EstimateResult{};
    if(cluster.points.size()<static_cast<size_t>(std::max(5,min_trajectory_points_))||!std::isfinite(field)||field==0)return false;
    std::vector<Point> points;
    for(const auto& p:cluster.points){
        if(!std::isfinite(p.x)||!std::isfinite(p.y)||!std::isfinite(p.z)||!std::isfinite(p.integral))return false;
        points.push_back({p.x,p.y,p.z,p.integral});
    }
    std::stable_sort(points.begin(),points.end(),[](const Point&a,const Point&b){return a.z<b.z;});
    std::vector<double> z,x,y,q;
    for(const auto& p:points){z.push_back(p.z);x.push_back(p.x);y.push_back(p.y);q.push_back(p.q);}
    Spline sx,sy,sq;
    if(!sx.Fit(z,x,reinsch_lambda_)||!sy.Fit(z,y,reinsch_lambda_)||!sq.Fit(z,q,reinsch_lambda_))return false;
    for(size_t i=0;i<points.size();++i){points[i].x=sx.Values()[i];points[i].y=sy.Values()[i];points[i].q=sq.Values()[i];}
    // Spyral 1.0 chooses direction from the smoothed endpoints and retries the
    // opposite direction only when the regression sign is inconsistent.
    const bool backward=!(std::hypot(points.front().x,points.front().y)<std::hypot(points.back().x,points.back().y));
    // The signed field controls trajectory curvature, but magnetic rigidity is
    // a magnitude. Keep the signed input at the API boundary while using |B|
    // for Bρ so Spyral's negative effective-field convention remains valid.
    const double field_magnitude=std::abs(field);
    const Pass status=EstimatePass(points,backward,sx,sy,field_magnitude,beam_region_radius_mm_,small_pad_radius_mm_,out);
    if(status==Pass::Accepted)return true;
    return status==Pass::WrongDirection &&
           EstimatePass(points,!backward,sx,sy,field_magnitude,beam_region_radius_mm_,small_pad_radius_mm_,out)==Pass::Accepted;
}
} // namespace rk4
