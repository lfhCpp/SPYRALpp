//----------------------------------------------------------------------
// 点云 HDF5 输出测试：验证布局与 attpc_engine 的 SpyralWriter 一致
//   root -l -b -q test/testH5Cloud.C
//----------------------------------------------------------------------
R__ADD_INCLUDE_PATH(src)
R__ADD_INCLUDE_PATH(extern/hdf5/include)
R__LOAD_LIBRARY(build/librk4.so)

#include "io/H5Writer.h"
#include "sim/Simulator.h"

#include "TSystem.h"

#include <hdf5.h>

#include <cstdio>
#include <vector>

using namespace rk4;

static int g_fail = 0;

static void Check(bool ok, const char* msg) {
    if (!ok) {
        std::printf("  [FAIL] %s\n", msg);
        ++g_fail;
    }
}

void testH5Cloud() {
    std::printf("========== 点云 HDF5 输出测试 ==========\n");

    if (!Hdf5Available()) {
        std::printf("========== SKIP（未启用 HDF5） ==========\n");
        gSystem->Exit(0);
    }

    // 造两个点云
    std::vector<PointCloud> clouds;
    for (int e = 0; e < 2; ++e) {
        PointCloud pc;
        pc.event = 5 + e;
        for (int k = 0; k < 10; ++k) {
            CloudPoint p;
            p.pad = 100 + k;
            p.tb = k * 2.0;
            p.x = 1.0 * k;
            p.y = 2.0 * k;
            p.z = 3.0 * k;
            p.amplitude = 400.0 + k;
            p.integral = 900.0 + k;
            p.scale = 1.0;
            pc.points.push_back(p);
        }
        clouds.push_back(pc);
    }

    const char* path = "/tmp/rk4_h5cloud_test.h5";
    gSystem->Unlink(path);
    Check(WritePointCloudsHdf5(path, clouds), "写出点云 HDF5");

    hid_t f = H5Fopen(path, H5F_ACC_RDONLY, H5P_DEFAULT);
    Check(f >= 0, "打开文件");

    hid_t g = H5Gopen2(f, "/cloud", H5P_DEFAULT);
    Check(g >= 0, "/cloud 组存在");

    if (g >= 0) {
        long long mn = -99, mx = -99;
        hid_t a = H5Aopen(g, "min_event", H5P_DEFAULT);
        if (a >= 0) { H5Aread(a, H5T_NATIVE_LLONG, &mn); H5Aclose(a); }
        hid_t b = H5Aopen(g, "max_event", H5P_DEFAULT);
        if (b >= 0) { H5Aread(b, H5T_NATIVE_LLONG, &mx); H5Aclose(b); }
        std::printf("    min_event=%lld max_event=%lld (期望 5/6)\n", mn, mx);
        Check(mn == 5 && mx == 6, "min/max_event 正确");
        H5Gclose(g);
    }

    for (int e = 0; e < 2; ++e) {
        char name[64];
        std::snprintf(name, sizeof(name), "/cloud/cloud_%d", 5 + e);
        hid_t d = H5Dopen2(f, name, H5P_DEFAULT);
        Check(d >= 0, "cloud_N 存在");
        if (d >= 0) {
            hid_t sp = H5Dget_space(d);
            hsize_t dims[2] = {0, 0};
            H5Sget_simple_extent_dims(sp, dims, nullptr);
            std::printf("    cloud_%d 形状 = (%llu, %llu)\n", 5 + e,
                        (unsigned long long)dims[0], (unsigned long long)dims[1]);
            Check(dims[0] == 10 && dims[1] == 8, "cloud_N 是 (10,8)");

            char an[32];
            std::snprintf(an, sizeof(an), "orig_event");
            hid_t oa = H5Aopen(d, an, H5P_DEFAULT);
            Check(oa >= 0, "有 orig_event 属性");
            if (oa >= 0) {
                long long v = -1;
                H5Aread(oa, H5T_NATIVE_LLONG, &v);
                Check(v == 5 + e, "orig_event 值正确");
                H5Aclose(oa);
            }
            H5Sclose(sp);
            H5Dclose(d);
        }

        std::snprintf(name, sizeof(name), "/cloud/labels_%d", 5 + e);
        hid_t l = H5Dopen2(f, name, H5P_DEFAULT);
        Check(l >= 0, "labels_N 存在");
        if (l >= 0) H5Dclose(l);
    }

    H5Fclose(f);
    gSystem->Unlink(path);

    if (g_fail == 0) {
        std::printf("========== 全部通过 ==========\n");
        gSystem->Exit(0);
    } else {
        std::printf("========== 失败 %d 项 ==========\n", g_fail);
        gSystem->Exit(1);
    }
}
