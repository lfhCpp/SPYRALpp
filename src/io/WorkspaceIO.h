#ifndef SPYRALPP_WORKSPACE_IO_H
#define SPYRALPP_WORKSPACE_IO_H
#include "../reco/Cluster.h"
#include <string>
#include <vector>
#include <variant>
#include <cstdint>
namespace rk4 {
struct StageFileInfo {
    long long run=-1,first_event=0,last_event=-1,raw_events=0,excluded_ic=0,good_downscale=0;
    uint64_t seed=20261026;
    std::string waveform_method="baseline";
};
struct PointcloudData {StageFileInfo info;std::vector<PointCloud> events;};
struct ClusterEvent {long long event=-1;EventMetadata metadata;std::vector<Cluster> clusters;};
struct ClusterData {StageFileInfo info;std::vector<ClusterEvent> events;};
enum class ColumnType {Integer,Real,Text};
struct Column {std::string name;ColumnType type;};
using Cell=std::variant<long long,double,std::string>;
struct ScalarTable {
    StageFileInfo info;
    std::vector<Column> columns;
    std::vector<std::vector<Cell>> rows;
    size_t Index(const std::string& name)const;
    double Number(size_t row,const std::string& name)const;
    long long Integer(size_t row,const std::string& name)const;
};
const std::vector<Column>& EstimationColumns();
const std::vector<Column>& SolverColumns();
PointcloudData ReadPointcloudStage(const std::string& path);
void WritePointcloudStage(const std::string& path,const PointcloudData& data);
ClusterData ReadClusterStage(const std::string& path);
void WriteClusterStage(const std::string& path,const ClusterData& data);
ScalarTable ReadScalarStage(const std::string& path,const std::string& stage);
void WriteScalarStage(const std::string& path,const std::string& stage,const ScalarTable& table);
// Native implementations, also used by the explicit format-conversion interface.
PointcloudData ReadPointcloudH5(const std::string& path);
void WritePointcloudH5(const std::string& path,const PointcloudData& data);
ClusterData ReadClusterH5(const std::string& path);
void WriteClusterH5(const std::string& path,const ClusterData& data);
}
#endif
