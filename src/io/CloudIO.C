#include "CloudIO.h"
#include "EventMetaIO.h"

#include "TFile.h"
#include "TTree.h"

#include <cstdio>
#include <cmath>
#include <vector>

namespace rk4 {
namespace {

struct WriteBuffer {
    EventMetadata metadata;
    Long64_t event = 0;
    Int_t npoints = 0;
    std::vector<double> x, y, z;
    std::vector<double> amplitude, integral;
    std::vector<int> pad;
    std::vector<double> tb, scale;
};

// 读取 vector 分支时，ROOT 要求传入「指向指针的指针」，
// 由 ROOT 自己管理 vector 的生命周期。
struct ReadBuffer {
    Long64_t event = 0;
    Int_t npoints = 0;
    std::vector<double>* x = nullptr;
    std::vector<double>* y = nullptr;
    std::vector<double>* z = nullptr;
    std::vector<double>* amplitude = nullptr;
    std::vector<double>* integral = nullptr;
    std::vector<int>* pad = nullptr;
    std::vector<double>* tb = nullptr;
    std::vector<double>* scale = nullptr;
};

void BindWrite(TTree& tree, WriteBuffer& b) {
    BindEventMetadata(tree,b.metadata);
    tree.Branch("event", &b.event, "event/L");
    tree.Branch("npoints", &b.npoints, "npoints/I");
    tree.Branch("x", &b.x);
    tree.Branch("y", &b.y);
    tree.Branch("z", &b.z);
    tree.Branch("amplitude", &b.amplitude);
    tree.Branch("integral", &b.integral);
    tree.Branch("pad", &b.pad);
    tree.Branch("tb", &b.tb);
    tree.Branch("scale", &b.scale);
}

bool BindRead(TTree& tree, ReadBuffer& b) {
    const auto bind = [&](const char* name, auto* value) {
        if (!tree.GetBranch(name) || tree.SetBranchAddress(name,value)<0) {
            std::fprintf(stderr,"[CloudIO] Missing or incompatible branch: %s\n",name);
            return false;
        }
        return true;
    };
    return bind("event",&b.event) && bind("npoints",&b.npoints) &&
           bind("x",&b.x) && bind("y",&b.y) && bind("z",&b.z) &&
           bind("amplitude",&b.amplitude) && bind("integral",&b.integral) &&
           bind("pad",&b.pad) && bind("tb",&b.tb) && bind("scale",&b.scale);
}

} // namespace

bool WritePointClouds(const std::string& path, const std::vector<PointCloud>& clouds,
                      const std::string& tree_name) {
    for (const auto& cloud : clouds) if (cloud.pending_dither) {
        std::fprintf(stderr,"[CloudIO] Cannot serialize unresolved time dithering\n");
        return false;
    }
    TFile file(path.c_str(), "RECREATE");
    if (file.IsZombie()) {
        std::fprintf(stderr, "[CloudIO] 无法创建 %s\n", path.c_str());
        return false;
    }

    TTree tree(tree_name.c_str(), "AT-TPC point clouds");
    WriteBuffer buf;
    BindWrite(tree, buf);

    for (const auto& cloud : clouds) {
        buf.event = cloud.event;
        buf.metadata=cloud.metadata;
        if(buf.metadata.orig_event<0)buf.metadata.orig_event=cloud.event;
        buf.npoints = static_cast<int>(cloud.points.size());

        buf.x.clear();
        buf.y.clear();
        buf.z.clear();
        buf.amplitude.clear();
        buf.integral.clear();
        buf.pad.clear();
        buf.tb.clear();
        buf.scale.clear();
        buf.x.reserve(cloud.points.size());
        buf.y.reserve(cloud.points.size());
        buf.z.reserve(cloud.points.size());
        buf.amplitude.reserve(cloud.points.size());
        buf.integral.reserve(cloud.points.size());
        buf.pad.reserve(cloud.points.size());
        buf.tb.reserve(cloud.points.size());
        buf.scale.reserve(cloud.points.size());

        for (const auto& p : cloud.points) {
            buf.x.push_back(p.x);
            buf.y.push_back(p.y);
            buf.z.push_back(p.z);
            buf.amplitude.push_back(p.amplitude);
            buf.integral.push_back(p.integral);
            buf.pad.push_back(p.pad);
            buf.tb.push_back(p.tb);
            buf.scale.push_back(p.scale);
        }
        tree.Fill();
    }

    tree.Write();
    file.Close();
    return true;
}

bool ReadPointClouds(const std::string& path, std::vector<PointCloud>& clouds,
                     const std::string& tree_name) {
    clouds.clear();
    TFile file(path.c_str(), "READ");
    if (file.IsZombie()) {
        std::fprintf(stderr, "[CloudIO] 无法打开 %s\n", path.c_str());
        return false;
    }

    TTree* tree = dynamic_cast<TTree*>(file.Get(tree_name.c_str()));
    if (tree == nullptr) {
        std::fprintf(stderr, "[CloudIO] 没有树 %s\n", tree_name.c_str());
        return false;
    }

    ReadBuffer buf;
    if (!BindRead(*tree, buf)) return false;

    std::vector<PointCloud> parsed;
    const Long64_t n = tree->GetEntries();
    parsed.reserve(static_cast<size_t>(n));

    for (Long64_t i = 0; i < n; ++i) {
        if (tree->GetEntry(i)<=0) {
            std::fprintf(stderr,"[CloudIO] Cannot read entry %lld\n",i);
            return false;
        }
        PointCloud cloud;
        cloud.event = buf.event;
        cloud.metadata=ReadEventMetadata(*tree,cloud.event);

        const size_t m = buf.x != nullptr ? buf.x->size() : 0;
        const auto size_ok = [&](const auto* column) { return column && column->size()==m; };
        if (buf.npoints<0 || static_cast<size_t>(buf.npoints)!=m ||
            !size_ok(buf.x)||!size_ok(buf.y)||!size_ok(buf.z)||!size_ok(buf.amplitude)||
            !size_ok(buf.integral)||!size_ok(buf.pad)||!size_ok(buf.tb)||!size_ok(buf.scale)) {
            std::fprintf(stderr,"[CloudIO] Inconsistent point columns in entry %lld\n",i);
            return false;
        }
        cloud.points.reserve(m);
        for (size_t j = 0; j < m; ++j) {
            CloudPoint p;
            p.x = (*buf.x)[j];
            p.y = (*buf.y)[j];
            p.z = (*buf.z)[j];
            p.amplitude = (*buf.amplitude)[j];
            p.integral = (*buf.integral)[j];
            p.pad = (*buf.pad)[j];
            p.tb = (*buf.tb)[j];
            p.scale = (*buf.scale)[j];
            if (!std::isfinite(p.x)||!std::isfinite(p.y)||!std::isfinite(p.z)||
                !std::isfinite(p.tb)||!std::isfinite(p.amplitude)||!std::isfinite(p.integral)||
                !std::isfinite(p.scale)) {
                std::fprintf(stderr,"[CloudIO] Nonfinite point in entry %lld\n",i);
                return false;
            }
            cloud.points.push_back(p);
        }
        parsed.push_back(std::move(cloud));
    }

    file.Close();
    clouds.swap(parsed);
    return true;
}

} // namespace rk4
