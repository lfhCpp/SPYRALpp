#ifndef SPYRALPP_APP_RECO_CONFIG_H
#define SPYRALPP_APP_RECO_CONFIG_H
#include "io/Config.h"
#include "reco/Cluster.h"
#include "reco/Estimator.h"
#include <cmath>
#include <stdexcept>
namespace rk4 {
inline ClusterParams ReadClusterParams(const Config& cfg){
    ClusterParams p;
    p.eps_mm=cfg.GetDouble("cluster.eps_mm",p.eps_mm);
    p.min_points=cfg.GetInt("cluster.min_points",p.min_points);
    p.min_cluster_size=cfg.GetInt("cluster.min_cluster_size",p.min_cluster_size);
    p.use_hdbscan=cfg.GetBool("cluster.use_hdbscan",p.use_hdbscan);
    p.min_cloud_size=cfg.GetInt("cluster.min_cloud_size",p.min_cloud_size);
    p.hdbscan_min_samples=cfg.GetInt("cluster.hdbscan_min_samples",p.min_points);
    p.min_size_scale_factor=cfg.GetDouble("cluster.min_size_scale_factor",p.min_size_scale_factor);
    p.hdbscan_z_scale=cfg.GetDouble("cluster.hdbscan_z_scale",p.hdbscan_z_scale);
    p.hdbscan_selection_epsilon=cfg.GetDouble("cluster.cluster_selection_epsilon",p.hdbscan_selection_epsilon);
    p.hdbscan_max_points=cfg.GetInt("cluster.hdbscan_max_points",p.hdbscan_max_points);
    p.hdbscan_allow_single=cfg.GetBool("cluster.hdbscan_allow_single",p.hdbscan_allow_single);
    p.overlap_min_cluster_size=cfg.GetInt("cluster.overlap_min_cluster_size",p.overlap_min_cluster_size);
    p.overlap_ratio=cfg.GetDouble("cluster.overlap_ratio",p.overlap_ratio);
    p.outlier_scale_factor=cfg.GetDouble("cluster.outlier_scale_factor",p.outlier_scale_factor);
    p.adaptive_eps=cfg.GetBool("cluster.adaptive_eps",p.adaptive_eps);
    p.eps_scale=cfg.GetDouble("cluster.eps_scale",p.eps_scale);
    p.eps_max_mm=cfg.GetDouble("cluster.eps_max_mm",p.eps_max_mm);
    p.enable_join=cfg.GetBool("cluster.enable_join",p.enable_join);
    p.join_max_gap_mm=cfg.GetDouble("cluster.join_gap_mm",p.join_max_gap_mm);
    p.join_max_angle_deg=cfg.GetDouble("cluster.join_angle_deg",p.join_max_angle_deg);
    p.join_edge_points=cfg.GetInt("cluster.join_edge_points",p.join_edge_points);
    if(p.min_points<1||p.min_cluster_size<2||p.min_cloud_size<1||p.hdbscan_min_samples<1||p.hdbscan_max_points<0||
       !(p.eps_mm>0)||!std::isfinite(p.eps_mm)||p.min_size_scale_factor<0||!std::isfinite(p.min_size_scale_factor)||
       !(p.hdbscan_z_scale>0)||!std::isfinite(p.hdbscan_z_scale)||p.hdbscan_selection_epsilon<0||!std::isfinite(p.hdbscan_selection_epsilon)||
       p.outlier_scale_factor<0||!std::isfinite(p.outlier_scale_factor)||p.overlap_ratio<0||p.overlap_ratio>1||!std::isfinite(p.overlap_ratio))
        throw std::invalid_argument("Invalid clustering configuration");
    return p;
}
inline void ConfigureEstimator(Estimator& estimator,const Config& cfg){
    const auto method=cfg.GetString("estimate.method","reinsch");
    const double lambda=cfg.GetDouble("estimate.reinsch_lambda",100);
    const int points=cfg.GetInt("estimate.min_total_trajectory_points",30);
    if(method!="reinsch"&&method!="legacy")throw std::invalid_argument("Unknown estimate.method: "+method);
    if(!std::isfinite(lambda)||lambda<0||points<5)throw std::invalid_argument("Invalid estimation parameters");
    estimator.SetUseReinsch(method=="reinsch",lambda);
    estimator.SetMinTrajectoryPoints(points);
    estimator.SetSmallPadRadius(cfg.GetDouble("detector.small_pad_radius",152));
    estimator.SetBeamRegionRadius(cfg.GetDouble("detector.beam_region_radius",25));
    estimator.SetSmoothing(cfg.GetBool("cluster.smooth",true),cfg.GetDouble("cluster.smooth_lambda",1));
}
}
#endif
