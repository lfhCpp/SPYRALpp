// Reproducible estimation from shared clusters, independent of clustering choices.
#include "reco/Estimator.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>

int main(int argc, char** argv) {
    if (argc < 3 || argc > 5) {
        std::cerr << "Usage: rk4_estimate_csv clusters.csv estimates.csv [bfield_T=2.85] [lambda=100]\n";
        return 1;
    }
    try {
        const auto start=std::chrono::steady_clock::now();
        const double field=argc>3?std::stod(argv[3]):2.85;
        const double lambda=argc>4?std::stod(argv[4]):100.0;
        if (!std::isfinite(field) || field<=0 || !std::isfinite(lambda) || lambda<0)
            throw std::runtime_error("Field must be positive; lambda must be nonnegative; both must be finite");
        std::error_code path_error;
        if (std::filesystem::equivalent(argv[1],argv[2],path_error) ||
            std::filesystem::absolute(argv[1]).lexically_normal() ==
                std::filesystem::absolute(argv[2]).lexically_normal())
            throw std::runtime_error("Input and output must be different files");
        std::ifstream input(argv[1]);
        if (!input) throw std::runtime_error("Cannot open input CSV");
        std::string line;
        if (!std::getline(input,line)) throw std::runtime_error("Missing CSV header");
        if (!line.empty() && line.back()=='\r') line.pop_back();
        if (line!="event,cluster,x,y,z,integral") throw std::runtime_error("Expected header: event,cluster,x,y,z,integral");
        std::map<std::pair<long long,int>,rk4::Cluster> groups;
        size_t line_number=1;
        while (std::getline(input,line)) {
            ++line_number;
            if (line.empty()) continue;
            std::replace(line.begin(),line.end(),',',' ');
            std::istringstream row(line);
            long long event; int label; rk4::CloudPoint point;
            std::string extra;
            if (!(row>>event>>label>>point.x>>point.y>>point.z>>point.integral) || row>>extra)
                throw std::runtime_error("Malformed cluster CSV at line "+std::to_string(line_number));
            groups[{event,label}].points.push_back(point);
        }
        // Open the output after the input has been validated and fully read.
        std::ofstream output(argv[2]);
        if (!output) throw std::runtime_error("Cannot open output CSV");
        output<<"event,cluster,brho,dedx,polar_deg,arclength,npts,vertex_x,vertex_y,vertex_z,direction\n"<<std::setprecision(17);
        rk4::Estimator estimator; estimator.SetUseReinsch(true,lambda);
        size_t accepted=0;
        for (const auto& entry:groups) {
            rk4::EstimateResult result;
            if (!estimator.Estimate(entry.second,field,result)) continue;
            output<<entry.first.first<<','<<entry.first.second<<','<<result.brho<<','<<result.dEdx<<','
                  <<result.polar*180/std::acos(-1.0)<<','<<result.arclength<<','<<entry.second.points.size()<<','
                  <<result.vertex_x<<','<<result.vertex_y<<','<<result.vertex_z<<','<<result.direction<<'\n';
            ++accepted;
        }
        output.close();
        if (!output) throw std::runtime_error("Failed to write estimates CSV");
        const double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        std::cout<<"Reinsch estimation: accepted="<<accepted<<" / "<<groups.size()<<", elapsed="<<elapsed<<" s\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr<<error.what()<<'\n';return 1;
    }
}
