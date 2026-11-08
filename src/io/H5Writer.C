#include "H5Writer.h"

#include <highfive/H5File.hpp>

#include <algorithm>
#include <cstdio>
#include <string>

namespace rk4 {

bool Hdf5Available() { return true; }

bool WriteKinematicsHdf5(const std::string& path, const std::vector<Event>& events,
                         int chunk_size) {
    if (chunk_size <= 0) chunk_size = 1000000;

    try {
        HighFive::File file(path, HighFive::File::Truncate);

        auto data = file.createGroup("data");
        const long long n_events = static_cast<long long>(events.size());
        data.createAttribute("n_events", n_events);
        data.createAttribute("chunk_size", static_cast<long long>(chunk_size));

        // 粒子编号在同一生成器内固定，取第一条事件即可
        std::vector<long long> proton_numbers;
        std::vector<long long> mass_numbers;
        if (!events.empty()) {
            for (const auto& p : events.front().particles) {
                proton_numbers.push_back(p.Z);
                mass_numbers.push_back(p.A);
            }
        }
        data.createAttribute("proton_numbers", proton_numbers);
        data.createAttribute("mass_numbers", mass_numbers);

        const long long n_chunks = (n_events + chunk_size - 1) / chunk_size;
        data.createAttribute("n_chunks", n_chunks);

        for (long long c = 0; c < n_chunks; ++c) {
            auto chunk = data.createGroup("chunk_" + std::to_string(c));
            const long long begin = c * chunk_size;
            const long long end = std::min(begin + static_cast<long long>(chunk_size), n_events);

            chunk.createAttribute("min_event", begin);
            chunk.createAttribute("max_event", end - 1);

            for (long long i = begin; i < end; ++i) {
                const Event& ev = events[static_cast<size_t>(i)];
                const size_t n_particles = ev.particles.size();

                std::vector<std::vector<double>> matrix(
                    n_particles, std::vector<double>(4, 0.0));
                for (size_t j = 0; j < n_particles; ++j) {
                    const Particle& p = ev.particles[j];
                    matrix[j][0] = p.px;
                    matrix[j][1] = p.py;
                    matrix[j][2] = p.pz;
                    matrix[j][3] = p.E;
                }

                auto dataset = chunk.createDataSet("event_" + std::to_string(i), matrix);
                dataset.createAttribute("vertex_x", ev.vx);
                dataset.createAttribute("vertex_y", ev.vy);
                dataset.createAttribute("vertex_z", ev.vz);
            }
        }
    } catch (const std::exception& e) {
        std::fprintf(stderr, "[H5Writer] 写出失败：%s\n", e.what());
        return false;
    }

    return true;
}


bool WritePointCloudsHdf5(const std::string& path, const std::vector<PointCloud>& clouds,
                          const std::vector<std::vector<long long>>& labels) {
    try {
        HighFive::File file(path, HighFive::File::Truncate);
        auto cloud_group = file.createGroup("cloud");

        // 事件范围：取实际写出的最小/最大事件号
        long long min_event = -1;
        long long max_event = -1;
        if (!clouds.empty()) {
            min_event = clouds.front().event;
            max_event = clouds.front().event;
            for (const auto& c : clouds) {
                if (c.event < min_event) min_event = c.event;
                if (c.event > max_event) max_event = c.event;
            }
        }
        cloud_group.createAttribute("min_event", min_event);
        cloud_group.createAttribute("max_event", max_event);

        for (size_t i = 0; i < clouds.size(); ++i) {
            const PointCloud& pc = clouds[i];
            const size_t n = pc.points.size();

            // (M, 8) = [x, y, z, amplitude, integral, pad, tb, scale]
            std::vector<std::vector<double>> matrix(n, std::vector<double>(8, 0.0));
            for (size_t k = 0; k < n; ++k) {
                const CloudPoint& p = pc.points[k];
                matrix[k][0] = p.x;
                matrix[k][1] = p.y;
                matrix[k][2] = p.z;
                matrix[k][3] = p.amplitude;
                matrix[k][4] = p.integral;
                matrix[k][5] = static_cast<double>(p.pad);
                matrix[k][6] = p.tb;
                matrix[k][7] = p.scale;
            }

            auto dset = cloud_group.createDataSet("cloud_" + std::to_string(pc.event), matrix);
            dset.createAttribute("orig_run", static_cast<long long>(0));
            dset.createAttribute("orig_event", static_cast<long long>(pc.event));
            dset.createAttribute("ic_amplitude", -1.0);
            dset.createAttribute("ic_integral", -1.0);
            dset.createAttribute("ic_centroid", -1.0);
            dset.createAttribute("ic_multiplicity", -1.0);

            // labels：没有时写全 -1
            std::vector<long long> lab(n, -1);
            if (i < labels.size() && labels[i].size() == n) {
                lab = labels[i];
            }
            cloud_group.createDataSet("labels_" + std::to_string(pc.event), lab);
        }
    } catch (const std::exception& e) {
        std::fprintf(stderr, "[H5Writer] 点云写出失败：%s\n", e.what());
        return false;
    }

    return true;
}

} // namespace rk4
