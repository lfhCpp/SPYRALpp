#include "reco/Clusterer.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

namespace {
void Require(bool condition, const char* message) {
    if (!condition) {
        std::cerr << "clusterer_snapshot_regression: " << message << '\n';
        std::exit(1);
    }
}

bool SamePoint(const rk4::CloudPoint& a, const rk4::CloudPoint& b) {
    return a.pad == b.pad && a.tb == b.tb && a.x == b.x && a.y == b.y &&
           a.z == b.z && a.integral == b.integral && a.scale == b.scale;
}

bool SameClusters(const std::vector<rk4::Cluster>& a,
                  const std::vector<rk4::Cluster>& b) {
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i].label != b[i].label || a[i].original_label != b[i].original_label ||
            a[i].direction != b[i].direction || a[i].points.size() != b[i].points.size())
            return false;
        for (size_t j = 0; j < a[i].points.size(); ++j)
            if (!SamePoint(a[i].points[j], b[i].points[j])) return false;
    }
    return true;
}
}

int main() {
    rk4::ClusterParams params;
    params.min_points = 3;
    params.min_cluster_size = 8;
    params.min_cloud_size = 8;
    params.hdbscan_min_samples = 3;
    params.min_size_scale_factor = 0.0;
    params.hdbscan_selection_epsilon = 0.0;
    params.hdbscan_z_scale = 1.0;
    params.outlier_scale_factor = 0.0;

    // Two complementary arcs of one circle exercise the overlap merger and
    // produce a stable HDBSCAN partition with source-order labels.
    rk4::PointCloud cloud;
    constexpr int count = 80;
    for (int i = 0; i < count; ++i) {
        const double angle = 2.0 * std::acos(-1.0) * i / count;
        rk4::CloudPoint p;
        p.x = 45.0 * std::cos(angle);
        p.y = 45.0 * std::sin(angle);
        p.z = 0.7 * i;
        p.integral = 1.0;
        p.pad = i;
        p.tb = i * 2.0;
        cloud.points.push_back(p);
    }

    rk4::Clusterer clusterer;
    clusterer.SetParams(params);
    std::vector<rk4::Cluster> production, raw, replay_a, replay_d;
    std::vector<int> labels;
    clusterer.Run(cloud, production);
    clusterer.RunHdbscanSnapshot(cloud, raw, labels);
    clusterer.FinalizeHdbscanSnapshot(raw, true, replay_a);
    clusterer.FinalizeHdbscanSnapshot(raw, false, replay_d);

    Require(labels.size() == cloud.points.size(), "point label map must preserve input row count");
    for (size_t i = 0; i < labels.size(); ++i) {
        if (labels[i] >= 0) {
            Require(static_cast<size_t>(labels[i]) < raw.size(), "point label must index raw HDBSCAN clusters");
            const auto& cluster = raw[static_cast<size_t>(labels[i])];
            bool found = false;
            for (const auto& point : cluster.points)
                if (SamePoint(point, cloud.points[i])) { found = true; break; }
            Require(found, "point label must map back to the original pointcloud row");
        }
    }

    Require(SameClusters(production, replay_a),
            "snapshot replay with original merge must exactly match production Run");

    std::vector<rk4::Cluster> fragmented(2), merged_fragments, unmerged_fragments, experimental_fragments;
    fragmented[0].label = 0;
    fragmented[1].label = 1;
    constexpr int arc_points = 48;
    for (int i = 0; i < arc_points; ++i) {
        const double angle = 2.0 * std::acos(-1.0) * i / arc_points;
        for (int k = 0; k < 2; ++k) {
            rk4::CloudPoint p;
            p.x = 45.0 * std::cos(angle);
            p.y = 45.0 * std::sin(angle);
            p.z = 0.2 * i;
            p.integral = 1.0;
            p.pad = i;
            p.tb = i;
            fragmented[static_cast<size_t>(k)].points.push_back(p);
        }
    }
    clusterer.FinalizeHdbscanSnapshot(fragmented, true, merged_fragments);
    clusterer.FinalizeHdbscanSnapshot(fragmented, false, unmerged_fragments);
    Require(merged_fragments.size() == 1, "overlapping circle fragments must merge in mode A");
    Require(unmerged_fragments.size() == 2, "mode D must retain the two input fragments");
    size_t retained_a = 0, retained_d = 0;
    for (const auto& cluster : merged_fragments) retained_a += cluster.points.size();
    for (const auto& cluster : unmerged_fragments) retained_d += cluster.points.size();
    Require(retained_a == retained_d, "merge toggle must not change LOF point retention");

    int experimental_calls = 0;
    rk4::ClusterSnapshotTiming experimental_timing;
    const rk4::ExperimentalClusterMerge test_merge = [&](std::vector<rk4::Cluster>& clusters) {
        ++experimental_calls;
        if (clusters.size() != 2) return;
        clusters[0].points.insert(clusters[0].points.end(), clusters[1].points.begin(), clusters[1].points.end());
        clusters.erase(clusters.begin() + 1);
    };
    clusterer.FinalizeHdbscanSnapshot(fragmented, false, experimental_fragments,
                                      &experimental_timing, test_merge);
    Require(experimental_calls == 1, "experimental merge hook must execute once");
    Require(experimental_fragments.size() == 1,
            "experimental merge hook must run before common postprocessing");
    Require(experimental_timing.original_merge_seconds == 0.0,
            "Mode B experimental merge must leave original merger disabled");
    return 0;
}
