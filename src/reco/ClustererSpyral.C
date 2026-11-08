// Port of attpc_spyral 1.0 core/clusterize.py and core/cluster.py (GPL-3.0).
// LOF formula follows scikit-learn neighbors/_lof.py (BSD-3-Clause);
// see notices/scikit-learn-hdbscan.txt for copyright conditions/disclaimer.
#include "Clusterer.h"
#include "Geometry.h"
#include "NumpyArgsort.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <map>
#include <numeric>
#include <utility>
#include <vector>
namespace rk4 {
namespace {
constexpr double Pi=3.1415926535897932384626433832795;
struct Circle { double x=0,y=0,r=0;bool valid=false; };
Circle Fit(const Cluster& cluster){
    std::vector<double>x,y;for(const auto&p:cluster.points){x.push_back(p.x);y.push_back(p.y);}
    Circle circle;circle.valid=FitCircleSpyral(x,y,circle.x,circle.y,circle.r);return circle;
}
double Overlap(const Circle&a,const Circle&b){
    const double d=std::hypot(a.x-b.x,a.y-b.y);
    if(d>a.r+b.r)return 0;
    if(d<=std::abs(a.r-b.r))return Pi*std::min(a.r*a.r,b.r*b.r);
    const double t1=std::clamp((d*d+a.r*a.r-b.r*b.r)/(2*d*a.r),-1.0,1.0);
    const double t2=std::clamp((d*d+b.r*b.r-a.r*a.r)/(2*d*b.r),-1.0,1.0);
    const double product=(-d+a.r+b.r)*(d+a.r-b.r)*(d-a.r+b.r)*(d+a.r+b.r);
    return a.r*a.r*std::acos(t1)+b.r*b.r*std::acos(t2)-.5*std::sqrt(std::max(0.0,product));
}
void OverlapStep(std::vector<Cluster>& clusters,int min_size,double ratio,
                 std::vector<std::vector<int>>& source_ids,int iteration,
                 std::vector<OriginalMergeDecision>* decisions,
                 const ExperimentalOriginalMergeVeto& veto){
    const bool track_sources=decisions||static_cast<bool>(veto);
    const size_t n=clusters.size();
    std::vector<Circle> circles;for(const auto& c:clusters)circles.push_back(Fit(c));
    std::map<int,std::vector<size_t>> groups;
    for(size_t i=0;i<n;++i)groups[clusters[i].label]={i};
    std::map<int,std::vector<int>> group_sources;
    if(track_sources)for(size_t i=0;i<n;++i)group_sources[clusters[i].label]=source_ids[i];
    for(size_t i=0;i<n;++i){
        if(!circles[i].valid||circles[i].r<10||static_cast<int>(clusters[i].points.size())<min_size)continue;
        for(size_t j=0;j<n;++j){
            if(i==j||!circles[j].valid||static_cast<int>(clusters[j].points.size())<min_size)continue;
            const double intersection=Overlap(circles[i],circles[j]);
            const double denominator=Pi*std::min(circles[i].r*circles[i].r,circles[j].r*circles[j].r);
            const double fraction=denominator>0?intersection/denominator:0.0;
            if(intersection<=ratio*denominator)continue;
            auto& current=groups.at(clusters[i].label);
            const bool already_joined=std::find(current.begin(),current.end(),j)!=current.end();
            OriginalMergeDecision decision;
            if(decisions || veto){
                auto& d=decision;
                d.iteration=iteration;d.left_index=i;d.right_index=j;
                d.left_raw_cluster_ids=group_sources.at(clusters[i].label);
                d.right_raw_cluster_ids=group_sources.at(clusters[j].label);
                d.left_points=static_cast<int>(std::accumulate(current.begin(),current.end(),size_t{0},[&](size_t npoints,size_t k){return npoints+clusters[k].points.size();}));
                const auto& right_group=groups.at(clusters[j].label);
                d.right_points=static_cast<int>(std::accumulate(right_group.begin(),right_group.end(),size_t{0},[&](size_t npoints,size_t k){return npoints+clusters[k].points.size();}));
                d.left_circle_x=circles[i].x;d.left_circle_y=circles[i].y;d.left_circle_r=circles[i].r;
                d.right_circle_x=circles[j].x;d.right_circle_y=circles[j].y;d.right_circle_r=circles[j].r;
                d.overlap_fraction=fraction;d.threshold=ratio;d.selected_for_merge=!already_joined;
                d.reason=already_joined?"already_in_same_merge_group":"circle_overlap_threshold_pass";
            }
            if(already_joined)continue;
            if(veto && veto(decision)){
                decision.selected_for_merge=false;
                decision.vetoed=true;
                decision.reason="experimental_original_merge_veto";
                if(decisions)decisions->push_back(std::move(decision));
                continue;
            }
            if(decisions)decisions->push_back(std::move(decision));
            const int from=clusters[j].label,to=clusters[i].label;
            auto members=groups.at(from);groups.erase(from);
            if(track_sources){
                auto merged_sources=group_sources.at(to);
                const auto& from_sources=group_sources.at(from);
                merged_sources.insert(merged_sources.end(),from_sources.begin(),from_sources.end());
                std::sort(merged_sources.begin(),merged_sources.end());
                merged_sources.erase(std::unique(merged_sources.begin(),merged_sources.end()),merged_sources.end());
                group_sources[to]=std::move(merged_sources);group_sources.erase(from);
            }
            for(size_t k:members)clusters[k].label=to;
            current.insert(current.end(),members.begin(),members.end());
        }
    }
    std::vector<Cluster> joined;
    std::vector<std::vector<int>> joined_sources;
    for(const auto& entry:groups){
        Cluster c;c.label=entry.first;
        for(size_t index:entry.second)c.points.insert(c.points.end(),clusters[index].points.begin(),clusters[index].points.end());
        joined.push_back(std::move(c));
        if(track_sources)joined_sources.push_back(group_sources.at(entry.first));
    }
    clusters.swap(joined);
    if(track_sources)source_ids.swap(joined_sources);
}
// sklearn LOF with contamination='auto': k excludes self, offset=-1.5.
void RemoveOutliers(Cluster& cluster,double factor){
    const size_t n=cluster.points.size();if(n<2)return;
    const size_t k=std::min(n-1,std::max<size_t>(2,static_cast<size_t>(factor*n)));
    using Neighbor=std::pair<double,size_t>;
    std::vector<std::vector<Neighbor>> neighbors(n);
    std::vector<double> distance_k(n),density(n);
    for(size_t i=0;i<n;++i){
        std::vector<Neighbor> row;row.reserve(n-1);
        for(size_t j=0;j<n;++j){
            if(i==j)continue;
            const auto&a=cluster.points[i];const auto&b=cluster.points[j];
            const double x=a.x-b.x,y=a.y-b.y,z=a.z-b.z;
            row.push_back({std::sqrt(x*x+y*y+z*z),j});
        }
        std::partial_sort(row.begin(),row.begin()+k,row.end());row.resize(k);
        distance_k[i]=row.back().first;neighbors[i]=std::move(row);
    }
    for(size_t i=0;i<n;++i){
        double reach=0;for(const auto& neighbor:neighbors[i])reach+=std::max(neighbor.first,distance_k[neighbor.second]);
        density[i]=1/(reach/k+1e-10);
    }
    std::vector<CloudPoint> kept;kept.reserve(n);
    for(size_t i=0;i<n;++i){
        double lof=0;for(const auto& neighbor:neighbors[i])lof+=density[neighbor.second]/density[i];
        if(lof/k<=1.5)kept.push_back(cluster.points[i]);
    }
    cluster.points.swap(kept);
}
} // namespace
int Clusterer::FinalizeSpyral(std::vector<Cluster>& built,std::vector<Cluster>& clusters,
                              ClusterSnapshotTiming* timing,
                              const ExperimentalClusterMerge& experimental_merge,
                              std::vector<OriginalMergeDecision>* original_decisions,
                              bool apply_lof,
                              const ExperimentalOriginalMergeVeto& original_merge_veto){
    using Clock = std::chrono::steady_clock;
    const size_t before=built.size();
    const auto merge_started = Clock::now();
    if(params_.enable_join){
        std::vector<std::vector<int>> source_ids;
        if(original_decisions||original_merge_veto){source_ids.reserve(built.size());for(const auto& cluster:built)source_ids.push_back({cluster.label});}
        size_t previous=0;
        int iteration=0;
        while(built.size()>1&&built.size()!=previous){
            previous=built.size();++iteration;
            OverlapStep(built,params_.overlap_min_cluster_size,params_.overlap_ratio,
                        source_ids,iteration,original_decisions,original_merge_veto);
        }
    }
    if(timing)timing->original_merge_seconds=params_.enable_join?
        std::chrono::duration<double>(Clock::now()-merge_started).count():0.0;
    if(experimental_merge){
        const auto experimental_started=Clock::now();
        experimental_merge(built);
        if(timing)timing->experimental_merge_seconds=
            std::chrono::duration<double>(Clock::now()-experimental_started).count();
    }
    const auto postprocess_started = Clock::now();
    if(timing)timing->lof_seconds=0.0;
    last_join_count_=static_cast<int>(before-built.size());
    for(auto& cluster:built){
        if(cluster.points.size()<2)continue;
        std::vector<double> z;for(const auto&p:cluster.points)z.push_back(p.z);
        const auto order=NumpyArgsort(z);const auto original=cluster.points;
        for(size_t i=0;i<order.size();++i)cluster.points[i]=original[order[i]];
        if(apply_lof){
            const auto lof_started=Clock::now();
            RemoveOutliers(cluster,params_.outlier_scale_factor);
            if(timing)timing->lof_seconds+=std::chrono::duration<double>(Clock::now()-lof_started).count();
        }
        cluster.direction=Direction::None; // Python 1.0 estimator infers this after smoothing.
        cluster.original_label=cluster.label;
        cluster.label=static_cast<int>(clusters.size()); // Canonical cluster_index in Python artifacts.
        clusters.push_back(std::move(cluster));
    }
    if(timing)timing->postprocess_seconds=
        std::chrono::duration<double>(Clock::now()-postprocess_started).count();
    return static_cast<int>(clusters.size());
}
} // namespace rk4
