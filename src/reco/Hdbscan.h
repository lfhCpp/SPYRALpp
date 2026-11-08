#ifndef RK4_RECO_HDBSCAN_H
#define RK4_RECO_HDBSCAN_H
#include "Cluster.h"
#include <vector>
namespace rk4 {
// Exact mutual-reachability MST, condensed tree, EOM stability and epsilon
// selection, following sklearn.cluster.HDBSCAN. min_samples includes self.
class HdbscanClusterer {
public:
    int Run(const PointCloud& cloud,int min_samples,int min_cluster_size,
            bool allow_single_cluster=false,double selection_epsilon=0.0,double z_scale=1.0);
    std::vector<Cluster> clusters_;
    std::vector<int> point_labels_;
};
}
#endif
